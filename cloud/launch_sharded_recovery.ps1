param(
    [string]$ProjectId = "project-7f7ff467-5d61-4072-8f4",
    [string]$Account = "shreshth.rach1@gmail.com",
    [string]$LinuxUser = "shreshth_rach1",
    [Parameter(Mandatory = $true)]
    [string]$Bundle,
    [switch]$ConfirmCredit
)

$ErrorActionPreference = "Stop"
function Assert-NativeSuccess {
    param([string]$Action)
    if ($LASTEXITCODE -ne 0) {
        throw "$Action failed with exit code $LASTEXITCODE."
    }
}
if (-not $ConfirmCredit) {
    throw "Re-run with -ConfirmCredit only after verifying active welcome credit in Cloud Billing."
}
if (-not (Test-Path -LiteralPath $Bundle)) {
    throw "Bundle not found: $Bundle"
}

$workers = @(
    @{ Id = "00"; Name = "alderaan-rec-00"; Zone = "us-central1-a" },
    @{ Id = "01"; Name = "alderaan-rec-01"; Zone = "us-east1-b" }
)

$active = (gcloud.cmd auth list --filter="status:ACTIVE" --format="value(account)").Trim()
Assert-NativeSuccess "gcloud account lookup"
if ($active -ne $Account) {
    throw "Active gcloud account is '$active', expected '$Account'."
}
gcloud.cmd config set project $ProjectId | Out-Host
Assert-NativeSuccess "gcloud project selection"
$billing = (gcloud.cmd billing projects describe $ProjectId --format="value(billingEnabled)").Trim()
Assert-NativeSuccess "gcloud billing lookup"
if ($billing -ne "True") {
    throw "Billing is not enabled for $ProjectId."
}

Write-Host "Bundle SHA256:"
Get-FileHash -Algorithm SHA256 -LiteralPath $Bundle | Format-List | Out-Host
Write-Host "Creating two standard C3D workers. Max runtime is 24h and termination action is STOP."

$created = @()
try {
    foreach ($worker in $workers) {
        $existing = gcloud.cmd compute instances list `
            --project=$ProjectId `
            --filter="name=$($worker.Name)" `
            --format="value(name)"
        Assert-NativeSuccess "worker existence check for $($worker.Name)"
        if ($existing) {
            throw "Instance already exists: $($worker.Name). Refusing to reuse it implicitly."
        }
        gcloud.cmd compute instances create $worker.Name `
            --project=$ProjectId `
            --zone=$($worker.Zone) `
            --machine-type=c3d-highcpu-16 `
            --provisioning-model=STANDARD `
            --image-family=ubuntu-2204-lts `
            --image-project=ubuntu-os-cloud `
            --boot-disk-size=60GB `
            --boot-disk-type=pd-balanced `
            --boot-disk-auto-delete `
            --max-run-duration=24h `
            --instance-termination-action=STOP `
            --no-restart-on-failure `
            --maintenance-policy=MIGRATE `
            --metadata=enable-guest-attributes=TRUE | Out-Host
        Assert-NativeSuccess "create $($worker.Name)"
        $created += $worker
    }

    Start-Sleep -Seconds 75
    foreach ($worker in $workers) {
        $bundleName = Split-Path -Leaf $Bundle
        $remoteBundle = "/home/$LinuxUser/$bundleName"
        gcloud.cmd compute scp `
            $Bundle `
            "$LinuxUser@$($worker.Name):$remoteBundle" `
            --project=$ProjectId `
            --zone=$($worker.Zone) `
            --quiet | Out-Host
        Assert-NativeSuccess "upload bundle to $($worker.Name)"

        $remoteDir = "/home/$LinuxUser/cloud_published_inventory_missing_batch"
        $unit = "alderaan-shard-$($worker.Id)"
        $remote = "cd /home/$LinuxUser && rm -rf cloud_published_inventory_missing_batch && python3 -m zipfile -e '$remoteBundle' . && test -f '$remoteDir/target_shards/targets_shard_$($worker.Id).csv' && chmod +x '$remoteDir/'*.sh && sudo systemd-run --unit=$unit --property=User=$LinuxUser --property=WorkingDirectory=$remoteDir --setenv=HOME=/home/$LinuxUser /bin/bash '$remoteDir/bootstrap_shard.sh' $($worker.Id) 14"
        gcloud.cmd compute ssh "$LinuxUser@$($worker.Name)" `
            --project=$ProjectId `
            --zone=$($worker.Zone) `
            --command=$remote | Out-Host
        Assert-NativeSuccess "launch shard $($worker.Id)"
        Write-Host "Launched shard $($worker.Id) on $($worker.Name)."
    }
}
catch {
    Write-Warning $_
    Write-Warning "Stopping every worker created by this invocation."
    foreach ($worker in $created) {
        gcloud.cmd compute instances stop $worker.Name `
            --project=$ProjectId `
            --zone=$($worker.Zone) `
            --quiet | Out-Host
    }
    throw
}

Write-Host "Both shards launched. Do not delete either VM or disk before local archive verification."
Write-Host "Use check_sharded_recovery.ps1 to monitor guest status."
