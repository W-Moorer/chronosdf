param(
    [string]$BuildConfigDir = "..\..\build\chrono_sdf_contact_vs\Release",
    [string]$OutputDir = "build\chrono_sdf_contact_vs\Release\results_p4_w1_frozen",
    [int]$NumSeeds = 20,
    [UInt64]$BaseSeed = 100,
    [bool]$AllowP2ValidationFailure = $false,
    [bool]$RunM3 = $false,
    [bool]$GeneratePlots = $true,
    [bool]$PackageArtifacts = $true
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$runner = Join-Path $scriptDir "run_experiments.ps1"

if (-not (Test-Path $runner)) {
    throw "Missing runner script: $runner"
}

Write-Host "Running frozen P4-W1 protocol..."
Write-Host "BuildConfigDir=$BuildConfigDir"
Write-Host "OutputDir=$OutputDir"
Write-Host "NumSeeds=$NumSeeds BaseSeed=$BaseSeed"
Write-Host "AllowP2ValidationFailure=$AllowP2ValidationFailure"

& $runner `
    -BuildConfigDir $BuildConfigDir `
    -OutputDir $OutputDir `
    -NumSeeds $NumSeeds `
    -BaseSeed $BaseSeed `
    -RunP2 $true `
    -P2SceneDts "0.001,0.002,0.003" `
    -P2ScaleCounts "8,16,32,64" `
    -RunM3 $RunM3 `
    -ValidateP1 $true `
    -AllowP2ValidationFailure $AllowP2ValidationFailure `
    -GeneratePlots $GeneratePlots `
    -GenerateSummary $true `
    -WriteMetadata $true `
    -PackageArtifacts $PackageArtifacts

if ($LASTEXITCODE -ne 0) {
    throw "run_experiments.ps1 failed with code $LASTEXITCODE"
}

Write-Host "Frozen P4-W1 run completed."
