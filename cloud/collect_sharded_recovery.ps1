param(
    [string]$ProjectId = "project-7f7ff467-5d61-4072-8f4",
    [string]$LinuxUser = "shreshth_rach1",
    [string]$DestinationRoot = (Join-Path $PSScriptRoot "recovery_results")
)

$ErrorActionPreference = "Stop"
function Assert-NativeSuccess {
    param([string]$Action)
    if ($LASTEXITCODE -ne 0) {
        throw "$Action failed with exit code $LASTEXITCODE."
    }
}
$workers = @(
    @{ Id = "00"; Name = "alderaan-rec-00"; Zone = "us-central1-a" },
    @{ Id = "01"; Name = "alderaan-rec-01"; Zone = "us-east1-b" }
)

$stamp = Get-Date -Format "yyyyMMddTHHmmss"
$destination = Join-Path $DestinationRoot "recovery_shards_$stamp"
New-Item -ItemType Directory -Path $destination | Out-Null

foreach ($worker in $workers) {
    $state = (gcloud.cmd compute instances describe $worker.Name `
        --project=$ProjectId `
        --zone=$($worker.Zone) `
        --format="value(status)").Trim()
    Assert-NativeSuccess "describe $($worker.Name)"
    if ($state -ne "TERMINATED") {
        throw (
            "Refusing to collect shard $($worker.Id) while $($worker.Name) is $state. " +
            "Wait for the runner to finish and stop before collecting its final archive."
        )
    }
    gcloud.cmd compute instances start $worker.Name `
        --project=$ProjectId `
        --zone=$($worker.Zone) | Out-Host
    Assert-NativeSuccess "start $($worker.Name)"
    Start-Sleep -Seconds 60

    $archiveName = "alderaan_shard_$($worker.Id)_results.tar.gz"
    gcloud.cmd compute scp `
        "$LinuxUser@$($worker.Name):~/$archiveName" `
        (Join-Path $destination $archiveName) `
        --project=$ProjectId `
        --zone=$($worker.Zone) | Out-Host
    Assert-NativeSuccess "download archive for shard $($worker.Id)"
    gcloud.cmd compute scp `
        "$LinuxUser@$($worker.Name):~/$archiveName.sha256" `
        (Join-Path $destination "$archiveName.sha256") `
        --project=$ProjectId `
        --zone=$($worker.Zone) | Out-Host
    Assert-NativeSuccess "download hash for shard $($worker.Id)"

    $expected = (Get-Content (Join-Path $destination "$archiveName.sha256")).Split()[0].ToLower()
    $actual = (Get-FileHash -Algorithm SHA256 (Join-Path $destination $archiveName)).Hash.ToLower()
    if ($actual -ne $expected) {
        throw "Hash mismatch for shard $($worker.Id): expected $expected, got $actual"
    }
    gcloud.cmd compute instances stop $worker.Name `
        --project=$ProjectId `
        --zone=$($worker.Zone) `
        --quiet | Out-Host
    Assert-NativeSuccess "stop $($worker.Name)"
}

$archives = Get-ChildItem $destination -Filter "alderaan_shard_*_results.tar.gz" |
    Sort-Object Name |
    ForEach-Object { $_.FullName }
if ($archives.Count -ne 2) {
    throw "Expected two verified archives, found $($archives.Count)."
}

$projectRoot = Split-Path -Parent $PSScriptRoot
$merger = Join-Path $projectRoot "scripts\merge_sharded_recovery_archives.py"
$targets = Join-Path $PSScriptRoot "recovery_142\targets_missing_launchable.csv"
$merged = Join-Path $destination "alderaan_results_published_inventory_missing_merged.tar.gz"
$staging = Join-Path $destination "merged_staging"
& python $merger `
    @archives `
    --targets $targets `
    --staging $staging `
    --output $merged
if ($LASTEXITCODE -ne 0) {
    throw "Shard merge failed. VMs and disks have been retained."
}

Write-Host "Verified shard archives and merged recovery archive:"
Write-Host $merged
Get-FileHash -Algorithm SHA256 $merged | Format-List | Out-Host
Write-Host "Do not delete the two VMs until postprocessing validates the merged archive."
