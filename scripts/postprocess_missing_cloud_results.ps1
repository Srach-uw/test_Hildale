param(
    [Parameter(Mandatory=$true)]
    [string]$TarPath,
    [string]$RunId = "sagear_missing",
    [string]$Python = "C:\Users\shres\anaconda3\python.exe"
)

$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectDir = Join-Path $ScriptDir "alderaan_project"
$OutputsDir = Join-Path $ScriptDir "outputs"
$CodexOutputs = "C:\Users\shres\Documents\Codex\2026-06-29\8\outputs"

New-Item -ItemType Directory -Force -Path $ProjectDir | Out-Null
New-Item -ItemType Directory -Force -Path $OutputsDir | Out-Null
New-Item -ItemType Directory -Force -Path $CodexOutputs | Out-Null

Write-Host "Extracting $TarPath into $ProjectDir"
tar -xzf $TarPath -C $ProjectDir

$summary = Join-Path $OutputsDir "eccentricity_posterior_summary_$RunId.csv"
$coverage = Join-Path $OutputsDir "eccentricity_posterior_coverage_$RunId.csv"
$merged = Join-Path $OutputsDir "eccentricity_posterior_summary_merged_$RunId.csv"
$mergedCoverage = Join-Path $OutputsDir "eccentricity_posterior_coverage_merged_$RunId.csv"
$posteriorSubdir = "eccentricity_posteriors_$RunId"

Push-Location $ScriptDir
try {
    & $Python .\extract_eccentricity_posteriors.py `
        --sample .\outputs\canonical_sample_old_astropy_rawcc.csv `
        --run-id $RunId `
        --posterior-subdir $posteriorSubdir `
        --summary-out $summary `
        --coverage-out $coverage

    & $Python .\merge_posterior_summaries.py `
        --new $summary `
        --out $merged `
        --coverage-out $mergedCoverage

    & $Python .\hierarchical_rayleigh.py --summary $merged

    foreach ($path in @($summary, $coverage, $merged, $mergedCoverage, ".\outputs\rayleigh_population_fit.csv", ".\outputs\rayleigh_population_fit_transit_selection.csv")) {
        if (Test-Path $path) {
            Copy-Item $path -Destination (Join-Path $CodexOutputs (Split-Path -Leaf $path)) -Force
        }
    }
}
finally {
    Pop-Location
}

Write-Host "Done. Key outputs:"
Write-Host "  $summary"
Write-Host "  $coverage"
Write-Host "  $merged"
Write-Host "  $mergedCoverage"
