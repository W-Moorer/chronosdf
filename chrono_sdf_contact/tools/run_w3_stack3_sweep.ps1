param(
    [string]$BuildConfigDir = "..\..\build\chrono_sdf_contact_vs\Release",
    [string]$OutputDir = "build\chrono_sdf_contact_vs\Release\results_w3_stack3_sweep",
    [int]$NumSeeds = 20,
    [UInt64]$BaseSeed = 100,
    [double]$Duration = 1.2,
    [double]$Dt = 0.003,
    [string[]]$KmaxList = @("32", "64", "96", "128"),
    [string[]]$CellSizeList = @("0.04", "0.05", "0.06"),
    [string[]]$ClusteringModes = @("on", "off"),
    [string[]]$NarrowphaseModes = @("default")
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Split-TokenList {
    param([object]$Value)
    $chunks = New-Object System.Collections.Generic.List[string]
    if ($null -eq $Value) {
        return @()
    }
    if ($Value -is [System.Array]) {
        foreach ($item in $Value) {
            if ($null -ne $item) {
                $chunks.Add([string]$item)
            }
        }
    } else {
        $chunks.Add([string]$Value)
    }

    $tokens = New-Object System.Collections.Generic.List[string]
    foreach ($chunk in $chunks) {
        foreach ($token in ($chunk -split "[,;\s]+")) {
            $t = $token.Trim()
            if ([string]::IsNullOrWhiteSpace($t)) {
                continue
            }
            $tokens.Add($t)
        }
    }
    return $tokens.ToArray()
}

function Parse-IntList {
    param([object]$Text)
    $out = New-Object System.Collections.Generic.List[int]
    foreach ($t in (Split-TokenList $Text)) {
        $v = 0
        if (-not [int]::TryParse($t, [ref]$v) -or $v -le 0) {
            throw "Invalid int token in list: '$t'"
        }
        $out.Add($v)
    }
    if ($out.Count -eq 0) {
        throw "Empty int list: $Text"
    }
    return $out.ToArray()
}

function Parse-DoubleList {
    param([object]$Text)
    $out = New-Object System.Collections.Generic.List[double]
    foreach ($t in (Split-TokenList $Text)) {
        $v = 0.0
        if (-not [double]::TryParse($t, [System.Globalization.NumberStyles]::Float, [System.Globalization.CultureInfo]::InvariantCulture, [ref]$v) -or $v -le 0.0) {
            throw "Invalid double token in list: '$t'"
        }
        $out.Add($v)
    }
    if ($out.Count -eq 0) {
        throw "Empty double list: $Text"
    }
    return $out.ToArray()
}

function Parse-ClusteringModes {
    param([object]$Text)
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($token in (Split-TokenList $Text)) {
        $t = $token.ToLowerInvariant()
        if ($t -ne "on" -and $t -ne "off") {
            throw "Invalid clustering mode token: '$t' (expected on/off)"
        }
        if (-not $out.Contains($t)) {
            $out.Add($t)
        }
    }
    if ($out.Count -eq 0) {
        throw "Empty clustering mode list: $Text"
    }
    return $out.ToArray()
}

function Parse-NarrowphaseModes {
    param([object]$Text)
    $valid = @("default", "no_cache", "no_cull", "no_cache_no_cull")
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($token in (Split-TokenList $Text)) {
        $t = $token.ToLowerInvariant()
        if (-not ($valid -contains $t)) {
            throw "Invalid narrowphase mode token: '$t' (expected one of: $($valid -join ', '))"
        }
        if (-not $out.Contains($t)) {
            $out.Add($t)
        }
    }
    if ($out.Count -eq 0) {
        throw "Empty narrowphase mode list: $Text"
    }
    return $out.ToArray()
}

function Mean-Value {
    param([double[]]$Values)
    if ($null -eq $Values -or $Values.Count -eq 0) {
        return $null
    }
    return ($Values | Measure-Object -Average).Average
}

if ($NumSeeds -lt 1) {
    throw "NumSeeds must be >= 1"
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..\..")).Path
$binDir = (Resolve-Path (Join-Path $scriptDir $BuildConfigDir)).Path
$outDir = if ([System.IO.Path]::IsPathRooted($OutputDir)) { $OutputDir } else { Join-Path $repoRoot $OutputDir }
New-Item -ItemType Directory -Path $outDir -Force | Out-Null

$exe = Join-Path $binDir "demo_MBS_sdf_p2_scene_matrix.exe"
if (-not (Test-Path $exe)) {
    throw "Executable not found: $exe"
}

$kmaxValues = Parse-IntList $KmaxList
$cellValues = Parse-DoubleList $CellSizeList
$clusterModes = Parse-ClusteringModes $ClusteringModes
$narrowphaseModes = Parse-NarrowphaseModes $NarrowphaseModes

$detailRows = New-Object System.Collections.Generic.List[object]
$summaryRows = New-Object System.Collections.Generic.List[object]

Write-Host "Running W3 stack3 sweep..."
Write-Host "OutputDir=$outDir"
Write-Host "NumSeeds=$NumSeeds BaseSeed=$BaseSeed Dt=$Dt"
Write-Host ("Kmax=[{0}] Cell=[{1}] Cluster=[{2}] NP=[{3}]" -f `
    ($kmaxValues -join ","), ($cellValues -join ","), ($clusterModes -join ","), ($narrowphaseModes -join ","))
Write-Host "List args hint: use comma-delimited strings, e.g. -NarrowphaseModes ""default,no_cache,no_cull"""

foreach ($clusterMode in $clusterModes) {
    foreach ($narrowphaseMode in $narrowphaseModes) {
        foreach ($kmax in $kmaxValues) {
            foreach ($cell in $cellValues) {
            $sdfExploded = New-Object System.Collections.Generic.List[double]
            $baseExploded = New-Object System.Collections.Generic.List[double]
            $sdfPen = New-Object System.Collections.Generic.List[double]
            $basePen = New-Object System.Collections.Generic.List[double]
            $sdfWall = New-Object System.Collections.Generic.List[double]
            $baseWall = New-Object System.Collections.Generic.List[double]
            $disableCluster = ($clusterMode -eq "off")
            $disableCache = ($narrowphaseMode -eq "no_cache" -or $narrowphaseMode -eq "no_cache_no_cull")
            $disableCull = ($narrowphaseMode -eq "no_cull" -or $narrowphaseMode -eq "no_cache_no_cull")

            for ($i = 0; $i -lt $NumSeeds; $i++) {
                $seed = $BaseSeed + [UInt64]$i
                $csvPath = Join-Path $outDir ("stack3_seed_{0:D3}_{1}_mode_{2}_{3}_k{4}_c{5}.csv" -f `
                    $i, $seed, $clusterMode, $narrowphaseMode, $kmax, ($cell.ToString("0.###", [System.Globalization.CultureInfo]::InvariantCulture)))

                $cmdArgs = New-Object System.Collections.Generic.List[string]
                $cmdArgs.Add("--csv")
                $cmdArgs.Add($csvPath)
                $cmdArgs.Add("--duration")
                $cmdArgs.Add("$Duration")
                $cmdArgs.Add("--dts")
                $cmdArgs.Add("$Dt")
                $cmdArgs.Add("--kmax")
                $cmdArgs.Add("$kmax")
                $cmdArgs.Add("--cell-size")
                $cmdArgs.Add("$cell")
                $cmdArgs.Add("--seed")
                $cmdArgs.Add("$seed")
                if ($disableCluster) {
                    $cmdArgs.Add("--no-cluster")
                }
                if ($disableCache) {
                    $cmdArgs.Add("--no-cache")
                }
                if ($disableCull) {
                    $cmdArgs.Add("--no-bound-cull")
                }

                & $exe $cmdArgs.ToArray()
                if ($LASTEXITCODE -ne 0 -and $LASTEXITCODE -ne 3) {
                    throw "demo_MBS_sdf_p2_scene_matrix failed with code $LASTEXITCODE (seed=$seed, mode=$clusterMode, np=$narrowphaseMode, kmax=$kmax, cell=$cell)"
                }
                if (-not (Test-Path $csvPath)) {
                    throw "Expected CSV not found: $csvPath"
                }

                $rows = Import-Csv -Path $csvPath | Where-Object {
                    $_.scene -eq "stack3_boxes" -and [Math]::Abs([double]$_.step_size - $Dt) -lt 1e-12
                }
                if ($rows.Count -lt 2) {
                    throw "Missing stack3 rows in CSV: $csvPath"
                }

                $base = $rows | Where-Object { $_.mode -eq "chrono_baseline" } | Select-Object -First 1
                $sdf = $rows | Where-Object { $_.mode -eq "sdf_contact" } | Select-Object -First 1
                if ($null -eq $base -or $null -eq $sdf) {
                    throw "Missing baseline or sdf row in CSV: $csvPath"
                }

                $baseExploded.Add([double]$base.exploded)
                $sdfExploded.Add([double]$sdf.exploded)
                $basePen.Add([double]$base.max_penetration)
                $sdfPen.Add([double]$sdf.max_penetration)
                $baseWall.Add([double]$base.wall_time_s)
                $sdfWall.Add([double]$sdf.wall_time_s)

                $detailRows.Add([pscustomobject]@{
                    seed = "$seed"
                    clustering_mode = $clusterMode
                    clustering_enabled = if ($disableCluster) { 0 } else { 1 }
                    narrowphase_mode = $narrowphaseMode
                    cache_enabled = if ($disableCache) { 0 } else { 1 }
                    bound_cull_enabled = if ($disableCull) { 0 } else { 1 }
                    kmax = $kmax
                    cell_size = $cell
                    mode = "chrono_baseline"
                    exploded = [double]$base.exploded
                    max_penetration = [double]$base.max_penetration
                    wall_time_s = [double]$base.wall_time_s
                })
                $detailRows.Add([pscustomobject]@{
                    seed = "$seed"
                    clustering_mode = $clusterMode
                    clustering_enabled = if ($disableCluster) { 0 } else { 1 }
                    narrowphase_mode = $narrowphaseMode
                    cache_enabled = if ($disableCache) { 0 } else { 1 }
                    bound_cull_enabled = if ($disableCull) { 0 } else { 1 }
                    kmax = $kmax
                    cell_size = $cell
                    mode = "sdf_contact"
                    exploded = [double]$sdf.exploded
                    max_penetration = [double]$sdf.max_penetration
                    wall_time_s = [double]$sdf.wall_time_s
                })
            }

            $baseExplodedMean = Mean-Value $baseExploded.ToArray()
            $sdfExplodedMean = Mean-Value $sdfExploded.ToArray()
            $basePenMean = Mean-Value $basePen.ToArray()
            $sdfPenMean = Mean-Value $sdfPen.ToArray()
            $baseWallMean = Mean-Value $baseWall.ToArray()
            $sdfWallMean = Mean-Value $sdfWall.ToArray()
            $overheadRatio = $null
            if ($baseWallMean -ne $null -and $baseWallMean -gt 1e-12) {
                $overheadRatio = $sdfWallMean / $baseWallMean
            }

            $summaryRows.Add([pscustomobject]@{
                clustering_mode = $clusterMode
                clustering_enabled = if ($disableCluster) { 0 } else { 1 }
                narrowphase_mode = $narrowphaseMode
                cache_enabled = if ($disableCache) { 0 } else { 1 }
                bound_cull_enabled = if ($disableCull) { 0 } else { 1 }
                kmax = $kmax
                cell_size = $cell
                dt = $Dt
                num_seeds = $NumSeeds
                baseline_exploded_rate = $baseExplodedMean
                sdf_exploded_rate = $sdfExplodedMean
                exploded_delta_sdf_minus_baseline = if ($baseExplodedMean -ne $null -and $sdfExplodedMean -ne $null) { $sdfExplodedMean - $baseExplodedMean } else { $null }
                baseline_max_penetration_mean = $basePenMean
                sdf_max_penetration_mean = $sdfPenMean
                baseline_wall_time_mean = $baseWallMean
                sdf_wall_time_mean = $sdfWallMean
                overhead_ratio_sdf_over_baseline = $overheadRatio
            })

            Write-Host ("mode={0}, np={1}, kmax={2}, cell={3}, baseline_exploded={4:P1}, sdf_exploded={5:P1}, overhead={6:N3}" -f `
                $clusterMode, $narrowphaseMode, $kmax, $cell, $baseExplodedMean, $sdfExplodedMean, $overheadRatio)
            }
        }
    }
}

$detailCsv = Join-Path $outDir "w3_stack3_dt003_detail.csv"
$summaryCsv = Join-Path $outDir "w3_stack3_dt003_summary.csv"

$detailRows | Export-Csv -Path $detailCsv -NoTypeInformation -Encoding UTF8
$summaryRows | Export-Csv -Path $summaryCsv -NoTypeInformation -Encoding UTF8

Write-Host "W3 sweep completed."
Write-Host "Detail CSV: $detailCsv"
Write-Host "Summary CSV: $summaryCsv"
