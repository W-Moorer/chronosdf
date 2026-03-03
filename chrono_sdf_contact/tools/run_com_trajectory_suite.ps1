param(
    [string]$BuildConfigDir = "..\..\build\chrono_sdf_contact_vs\Release",
    [string]$OutputDir = "build\chrono_sdf_contact_vs\Release\results_com_trajectory_validation",
    [UInt64]$Seed = 100,
    [int]$MeshSteps = 1200,
    [double]$MeshDt = 0.001,
    [double]$M2Duration = 1.5,
    [string]$M2Dts = "0.001,0.002,0.004",
    [double]$P2SceneDuration = 1.2,
    [string]$P2SceneDts = "0.001,0.002,0.003",
    [double]$P2ScaleDuration = 0.8,
    [double]$P2ScaleDt = 0.0015,
    [string]$P2ScaleCounts = "8,16,32,64",
    [string]$PythonExe = "python"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..\..")).Path
$binDir = (Resolve-Path (Join-Path $scriptDir $BuildConfigDir)).Path
$outDir = if ([System.IO.Path]::IsPathRooted($OutputDir)) { $OutputDir } else { Join-Path $repoRoot $OutputDir }
New-Item -ItemType Directory -Path $outDir -Force | Out-Null

$meshExe = Join-Path $binDir "demo_MBS_sdf_mesh_baseline_compare.exe"
$m2Exe = Join-Path $binDir "demo_MBS_sdf_sdf_stability_compare.exe"
$p2SceneExe = Join-Path $binDir "demo_MBS_sdf_p2_scene_matrix.exe"
$p2ScaleExe = Join-Path $binDir "demo_MBS_sdf_p2_scaleup.exe"
$analysisScript = Join-Path $scriptDir "analyze_com_trajectories.py"

foreach ($exe in @($meshExe, $m2Exe, $p2SceneExe, $p2ScaleExe)) {
    if (-not (Test-Path $exe)) {
        throw "Executable not found: $exe"
    }
}
if (-not (Test-Path $analysisScript)) {
    throw "Missing analysis script: $analysisScript"
}

$meshCsv = Join-Path $outDir "sdf_mesh_baseline_compare.csv"
$meshTrajCsv = Join-Path $outDir "sdf_mesh_baseline_compare_trajectory.csv"
$m2Csv = Join-Path $outDir "sdf_sdf_stability_compare.csv"
$m2TrajCsv = Join-Path $outDir "sdf_sdf_stability_compare_trajectory.csv"
$p2SceneCsv = Join-Path $outDir "sdf_p2_scene_matrix.csv"
$p2SceneTrajCsv = Join-Path $outDir "sdf_p2_scene_matrix_trajectory.csv"
$p2ScaleCsv = Join-Path $outDir "sdf_p2_scaleup.csv"
$p2ScaleTrajCsv = Join-Path $outDir "sdf_p2_scaleup_trajectory.csv"

Write-Host "Running M1 mesh baseline-vs-sdf trajectory capture..."
& $meshExe --csv $meshCsv --trajectory-csv $meshTrajCsv --steps "$MeshSteps" --dt "$MeshDt" --seed "$Seed"
if ($LASTEXITCODE -ne 0) {
    throw "demo_MBS_sdf_mesh_baseline_compare failed with code $LASTEXITCODE"
}

Write-Host "Running M2 SDF-SDF trajectory capture..."
& $m2Exe --csv $m2Csv --trajectory-csv $m2TrajCsv --duration "$M2Duration" --dts "$M2Dts" --seed "$Seed"
if ($LASTEXITCODE -ne 0) {
    throw "demo_MBS_sdf_sdf_stability_compare failed with code $LASTEXITCODE"
}

Write-Host "Running P2 scene-matrix trajectory capture..."
& $p2SceneExe --csv $p2SceneCsv --trajectory-csv $p2SceneTrajCsv --duration "$P2SceneDuration" --dts "$P2SceneDts" --seed "$Seed"
if ($LASTEXITCODE -ne 0) {
    throw "demo_MBS_sdf_p2_scene_matrix failed with code $LASTEXITCODE"
}

Write-Host "Running P2 scale-up trajectory capture..."
& $p2ScaleExe --csv $p2ScaleCsv --trajectory-csv $p2ScaleTrajCsv --duration "$P2ScaleDuration" --dt "$P2ScaleDt" --counts "$P2ScaleCounts" --seed "$Seed"
if ($LASTEXITCODE -ne 0) {
    throw "demo_MBS_sdf_p2_scaleup failed with code $LASTEXITCODE"
}

Write-Host "Analyzing COM trajectories and generating figures..."
& $PythonExe $analysisScript `
    --mesh-traj-csv $meshTrajCsv `
    --scene-traj-csv $p2SceneTrajCsv `
    --scale-traj-csv $p2ScaleTrajCsv `
    --m2-traj-csv $m2TrajCsv `
    --out-dir $outDir `
    --out-summary "com_trajectory_diff_summary.md"
if ($LASTEXITCODE -ne 0) {
    throw "analyze_com_trajectories.py failed with code $LASTEXITCODE"
}

Write-Host "COM trajectory suite completed."
Write-Host "Output directory: $outDir"
