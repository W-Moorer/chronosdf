param(
    [string]$BuildConfigDir = "..\..\build\chrono_sdf_contact_vs\Release",
    [string]$OutputDir = ".\results",
    [int]$NumSeeds = 1,
    [UInt64]$BaseSeed = 42,
    [int]$MeshSteps = 1200,
    [double]$MeshDt = 0.001,
    [double]$SdfSdfDuration = 1.5,
    [string]$SdfSdfDts = "0.001,0.002,0.004",
    [bool]$RunP2 = $true,
    [double]$P2SceneDuration = 1.2,
    [string]$P2SceneDts = "0.001,0.002,0.003",
    [double]$P2ScaleDuration = 0.8,
    [double]$P2ScaleDt = 0.0015,
    [string]$P2ScaleCounts = "8,16,32,64",
    [bool]$RunM3 = $true,
    [string]$M3EpsList = "4000,2000,1000,500,250,125",
    [int]$M3Queries = 30000,
    [int]$M3MaxDepth = 8,
    [double]$M3Kn = 60000,
    [double]$M3CError = 0.5,
    [double]$M3Band = 0.25,
    [double]$M3HalfExtent = 1.0,
    [double]$M3Radius = 0.45,
    [bool]$ValidateP1 = $true,
    [double]$P1MinFitR2 = 0.95,
    [double]$P1MinSlope = 1.5,
    [double]$P1MaxSlope = 2.5,
    [double]$P1MaxBoundViolation = 1e-8,
    [double]$P1MinMapR2 = 0.4,
    [double]$P1MaxMapRelMae = 0.6,
    [double]$P1MaxRmseRatio = 1.0,
    [bool]$AllowP2ValidationFailure = $false,
    [string]$PythonExe = "python",
    [bool]$GeneratePlots = $true,
    [bool]$GenerateSummary = $true,
    [bool]$WriteMetadata = $true,
    [bool]$PackageArtifacts = $false
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path $scriptDir "..\..")).Path
$binDir = (Resolve-Path (Join-Path $scriptDir $BuildConfigDir)).Path
$outDir = if ([System.IO.Path]::IsPathRooted($OutputDir)) { $OutputDir } else { Join-Path $repoRoot $OutputDir }
New-Item -ItemType Directory -Path $outDir -Force | Out-Null

$meshExe = Join-Path $binDir "demo_MBS_sdf_mesh_baseline_compare.exe"
$sdfSdfExe = Join-Path $binDir "demo_MBS_sdf_sdf_stability_compare.exe"
$p2SceneExe = Join-Path $binDir "demo_MBS_sdf_p2_scene_matrix.exe"
$p2ScaleExe = Join-Path $binDir "demo_MBS_sdf_p2_scaleup.exe"
$m3Exe = Join-Path $binDir "demo_MBS_sdf_octree_pareto.exe"

if (-not (Test-Path $meshExe)) {
    throw "Executable not found: $meshExe"
}
if (-not (Test-Path $sdfSdfExe)) {
    throw "Executable not found: $sdfSdfExe"
}
if ($RunP2 -and -not (Test-Path $p2SceneExe)) {
    throw "Executable not found: $p2SceneExe"
}
if ($RunP2 -and -not (Test-Path $p2ScaleExe)) {
    throw "Executable not found: $p2ScaleExe"
}
if ($RunM3 -and -not (Test-Path $m3Exe)) {
    throw "Executable not found: $m3Exe"
}

$meshCsv = Join-Path $outDir "sdf_mesh_baseline_compare.csv"
$sdfSdfCsv = Join-Path $outDir "sdf_sdf_stability_compare.csv"
$p2SceneCsv = Join-Path $outDir "sdf_p2_scene_matrix.csv"
$p2ScaleCsv = Join-Path $outDir "sdf_p2_scaleup.csv"
$m3Csv = Join-Path $outDir "sdf_octree_pareto.csv"
$m3ParetoCsv = Join-Path $outDir "sdf_octree_pareto_front.csv"
$summaryMd = Join-Path $outDir "results_summary.md"
$plotScript = Join-Path $scriptDir "plot_results.py"
$plotP2Script = Join-Path $scriptDir "plot_p2_results.py"
$plotM3Script = Join-Path $scriptDir "plot_m3_pareto.py"
$summaryScript = Join-Path $scriptDir "summarize_results.py"
$validateP1Script = Join-Path $scriptDir "validate_p1.py"
$p1ValidationMd = Join-Path $outDir "p1_validation.md"

if ($NumSeeds -lt 1) {
    throw "NumSeeds must be >= 1"
}

function Invoke-DemoChecked {
    param(
        [string]$ExePath,
        [string]$Name,
        [string[]]$CommandArgs,
        [int[]]$AllowedExitCodes = @(0),
        [string[]]$RequiredOutputs = @()
    )
    & $ExePath @CommandArgs
    if ($AllowedExitCodes -contains $LASTEXITCODE) {
        foreach ($outPath in $RequiredOutputs) {
            if (-not [string]::IsNullOrWhiteSpace($outPath) -and -not (Test-Path $outPath)) {
                throw "$Name did not produce expected output: $outPath"
            }
        }
        if ($LASTEXITCODE -ne 0) {
            Write-Warning "$Name exited with accepted code $LASTEXITCODE"
        }
        return
    }
    if ($LASTEXITCODE -ne 0) {
        throw "$Name failed with code $LASTEXITCODE"
    }
}

function Get-CMakeCacheValue {
    param(
        [string]$CachePath,
        [string]$Key
    )
    if (-not (Test-Path $CachePath)) {
        return $null
    }
    $line = Select-String -Path $CachePath -Pattern ("^{0}:[^=]+=.*" -f [regex]::Escape($Key)) | Select-Object -First 1
    if (-not $line) {
        return $null
    }
    $text = $line.Line
    $idx = $text.IndexOf("=")
    if ($idx -lt 0) {
        return $null
    }
    return $text.Substring($idx + 1)
}

function Get-CommandOutputOrNull {
    param(
        [string]$Command,
        [string[]]$CommandArgs
    )
    try {
        $result = & $Command @CommandArgs 2>$null
        $lastExitVar = Get-Variable -Name LASTEXITCODE -ErrorAction SilentlyContinue
        if ($lastExitVar -and [int]$lastExitVar.Value -ne 0) {
            return $null
        }
        if ($null -eq $result) {
            return ""
        }
        return (($result -join "`n").Trim())
    } catch {
        return $null
    }
}

function Add-ArtifactIfExists {
    param(
        [System.Collections.Generic.List[object]]$List,
        [string]$Path
    )
    if (-not (Test-Path $Path)) {
        return
    }
    $item = Get-Item $Path
    $hash = (Get-FileHash -Path $Path -Algorithm SHA256).Hash
    $List.Add([ordered]@{
        name = $item.Name
        path = $item.FullName
        size_bytes = $item.Length
        sha256 = $hash
    })
}

function Merge-CsvWithSeed {
    param(
        [string[]]$CsvPaths,
        [UInt64[]]$Seeds,
        [string]$OutPath
    )
    if ($CsvPaths.Count -ne $Seeds.Count) {
        throw "CSV path count does not match seed count for $OutPath"
    }

    $merged = @()
    for ($idx = 0; $idx -lt $CsvPaths.Count; $idx++) {
        $rows = Import-Csv -Path $CsvPaths[$idx]
        foreach ($r in $rows) {
            $ordered = [ordered]@{ seed = "$($Seeds[$idx])" }
            foreach ($p in $r.PSObject.Properties) {
                $ordered[$p.Name] = $p.Value
            }
            $merged += [pscustomobject]$ordered
        }
    }
    if ($merged.Count -eq 0) {
        throw "No rows merged for $OutPath"
    }
    $merged | Export-Csv -Path $OutPath -NoTypeInformation -Encoding UTF8
}

Write-Host "Running M1/M2 demos with $NumSeeds seed(s), base seed=$BaseSeed..."
$meshSeedCsvs = New-Object System.Collections.Generic.List[string]
$stabilitySeedCsvs = New-Object System.Collections.Generic.List[string]
$p2SceneSeedCsvs = New-Object System.Collections.Generic.List[string]
$p2ScaleSeedCsvs = New-Object System.Collections.Generic.List[string]
$m3SeedCsvs = New-Object System.Collections.Generic.List[string]
$m3ParetoSeedCsvs = New-Object System.Collections.Generic.List[string]
$seedValues = New-Object System.Collections.Generic.List[UInt64]

for ($i = 0; $i -lt $NumSeeds; $i++) {
    $seed = $BaseSeed + [UInt64]$i
    $seedValues.Add($seed)

    $seedDir = if ($NumSeeds -eq 1) { $outDir } else { Join-Path $outDir ("seed_{0:D3}_{1}" -f $i, $seed) }
    New-Item -ItemType Directory -Path $seedDir -Force | Out-Null

    $meshSeedCsv = Join-Path $seedDir "sdf_mesh_baseline_compare.csv"
    $stabilitySeedCsv = Join-Path $seedDir "sdf_sdf_stability_compare.csv"
    $meshSeedCsvs.Add($meshSeedCsv)
    $stabilitySeedCsvs.Add($stabilitySeedCsv)

    Write-Host "Running mesh baseline compare (seed=$seed)..."
    Invoke-DemoChecked -ExePath $meshExe -Name "demo_MBS_sdf_mesh_baseline_compare" -CommandArgs @(
        "--csv", $meshSeedCsv, "--steps", "$MeshSteps", "--dt", "$MeshDt", "--seed", "$seed"
    ) -RequiredOutputs @($meshSeedCsv)

    Write-Host "Running SDF-SDF stability compare (seed=$seed)..."
    Invoke-DemoChecked -ExePath $sdfSdfExe -Name "demo_MBS_sdf_sdf_stability_compare" -CommandArgs @(
        "--csv", $stabilitySeedCsv, "--duration", "$SdfSdfDuration", "--dts", "$SdfSdfDts", "--seed", "$seed"
    ) -RequiredOutputs @($stabilitySeedCsv)

    if ($RunP2) {
        $p2SceneSeedCsv = Join-Path $seedDir "sdf_p2_scene_matrix.csv"
        $p2ScaleSeedCsv = Join-Path $seedDir "sdf_p2_scaleup.csv"
        $p2SceneSeedCsvs.Add($p2SceneSeedCsv)
        $p2ScaleSeedCsvs.Add($p2ScaleSeedCsv)
        $p2AllowedCodes = if ($AllowP2ValidationFailure) { @(0, 3) } else { @(0) }

        Write-Host "Running P2 scene matrix benchmark (seed=$seed)..."
        Invoke-DemoChecked -ExePath $p2SceneExe -Name "demo_MBS_sdf_p2_scene_matrix" -CommandArgs @(
            "--csv", $p2SceneSeedCsv, "--duration", "$P2SceneDuration", "--dts", "$P2SceneDts", "--seed", "$seed"
        ) -AllowedExitCodes $p2AllowedCodes -RequiredOutputs @($p2SceneSeedCsv)

        Write-Host "Running P2 scale-up benchmark (seed=$seed)..."
        Invoke-DemoChecked -ExePath $p2ScaleExe -Name "demo_MBS_sdf_p2_scaleup" -CommandArgs @(
            "--csv", $p2ScaleSeedCsv, "--duration", "$P2ScaleDuration", "--dt", "$P2ScaleDt", "--counts", "$P2ScaleCounts", "--seed", "$seed"
        ) -AllowedExitCodes $p2AllowedCodes -RequiredOutputs @($p2ScaleSeedCsv)
    }

    if ($RunM3) {
        $m3SeedCsv = Join-Path $seedDir "sdf_octree_pareto.csv"
        $m3SeedParetoCsv = Join-Path $seedDir "sdf_octree_pareto_front.csv"
        $m3SeedCsvs.Add($m3SeedCsv)
        $m3ParetoSeedCsvs.Add($m3SeedParetoCsv)

        Write-Host "Running M3 octree Pareto demo (seed=$seed)..."
        Invoke-DemoChecked -ExePath $m3Exe -Name "demo_MBS_sdf_octree_pareto" -CommandArgs @(
            "--csv", $m3SeedCsv, "--pareto-csv", $m3SeedParetoCsv, "--eps-list", $M3EpsList, "--queries", "$M3Queries",
            "--max-depth", "$M3MaxDepth", "--kn", "$M3Kn", "--c-error", "$M3CError", "--band", "$M3Band",
            "--half-extent", "$M3HalfExtent", "--radius", "$M3Radius", "--seed", "$seed"
        ) -RequiredOutputs @($m3SeedCsv, $m3SeedParetoCsv)
    }
}

Merge-CsvWithSeed -CsvPaths $meshSeedCsvs.ToArray() -Seeds $seedValues.ToArray() -OutPath $meshCsv
Merge-CsvWithSeed -CsvPaths $stabilitySeedCsvs.ToArray() -Seeds $seedValues.ToArray() -OutPath $sdfSdfCsv
if ($RunP2) {
    Merge-CsvWithSeed -CsvPaths $p2SceneSeedCsvs.ToArray() -Seeds $seedValues.ToArray() -OutPath $p2SceneCsv
    Merge-CsvWithSeed -CsvPaths $p2ScaleSeedCsvs.ToArray() -Seeds $seedValues.ToArray() -OutPath $p2ScaleCsv
}

if ($RunM3) {
    Merge-CsvWithSeed -CsvPaths $m3SeedCsvs.ToArray() -Seeds $seedValues.ToArray() -OutPath $m3Csv
    Merge-CsvWithSeed -CsvPaths $m3ParetoSeedCsvs.ToArray() -Seeds $seedValues.ToArray() -OutPath $m3ParetoCsv
}

if ($GeneratePlots) {
    Write-Host "Generating plots..."
    & $PythonExe $plotScript --mesh-csv $meshCsv --stability-csv $sdfSdfCsv --out-dir $outDir
    if ($LASTEXITCODE -ne 0) {
        throw "plot_results.py failed with code $LASTEXITCODE"
    }

    if ($RunM3) {
        & $PythonExe $plotM3Script --csv $m3Csv --pareto-csv $m3ParetoCsv --out (Join-Path $outDir "m3_octree_pareto.png")
        if ($LASTEXITCODE -ne 0) {
            throw "plot_m3_pareto.py failed with code $LASTEXITCODE"
        }
    }

    if ($RunP2) {
        & $PythonExe $plotP2Script --scene-csv $p2SceneCsv --scale-csv $p2ScaleCsv --out-dir $outDir
        if ($LASTEXITCODE -ne 0) {
            throw "plot_p2_results.py failed with code $LASTEXITCODE"
        }
    }
}

if ($GenerateSummary) {
    Write-Host "Generating markdown summary..."
    if ($RunM3 -and $RunP2) {
        & $PythonExe $summaryScript --mesh-csv $meshCsv --stability-csv $sdfSdfCsv --p2-scene-csv $p2SceneCsv --p2-scale-csv $p2ScaleCsv --m3-csv $m3Csv --out $summaryMd
    } elseif ($RunM3) {
        & $PythonExe $summaryScript --mesh-csv $meshCsv --stability-csv $sdfSdfCsv --m3-csv $m3Csv --out $summaryMd
    } elseif ($RunP2) {
        & $PythonExe $summaryScript --mesh-csv $meshCsv --stability-csv $sdfSdfCsv --p2-scene-csv $p2SceneCsv --p2-scale-csv $p2ScaleCsv --out $summaryMd
    } else {
        & $PythonExe $summaryScript --mesh-csv $meshCsv --stability-csv $sdfSdfCsv --out $summaryMd
    }
    if ($LASTEXITCODE -ne 0) {
        throw "summarize_results.py failed with code $LASTEXITCODE"
    }
}

if ($RunM3 -and $ValidateP1) {
    Write-Host "Running P1 validation checks..."
    & $PythonExe $validateP1Script --csv $m3Csv --out $p1ValidationMd `
        --min-fit-r2 $P1MinFitR2 --min-slope $P1MinSlope --max-slope $P1MaxSlope `
        --max-bound-violation $P1MaxBoundViolation --min-map-r2 $P1MinMapR2 `
        --max-map-rel-mae $P1MaxMapRelMae --max-rmse-ratio $P1MaxRmseRatio
    if ($LASTEXITCODE -ne 0) {
        throw "validate_p1.py failed with code $LASTEXITCODE"
    }
}

if ($WriteMetadata) {
    $manifestJson = Join-Path $outDir "artifact_manifest.json"
    $manifestMd = Join-Path $outDir "artifact_manifest.md"
    $cachePath = Join-Path (Split-Path -Parent $binDir) "CMakeCache.txt"

    $cpuName = $null
    $logicalCores = $null
    $memoryBytes = $null
    $osCaption = $null
    $osVersion = $null
    try {
        $cpu = Get-CimInstance Win32_Processor | Select-Object -First 1
        if ($cpu) {
            $cpuName = $cpu.Name
            $logicalCores = $cpu.NumberOfLogicalProcessors
        }
        $cs = Get-CimInstance Win32_ComputerSystem
        if ($cs) {
            $memoryBytes = [UInt64]$cs.TotalPhysicalMemory
        }
        $os = Get-CimInstance Win32_OperatingSystem
        if ($os) {
            $osCaption = $os.Caption
            $osVersion = $os.Version
        }
    } catch {
        # Ignore metadata probe failures; keep run script robust.
    }

    $gitCommit = Get-CommandOutputOrNull -Command "git" -CommandArgs @("-C", $repoRoot, "rev-parse", "HEAD")
    $gitBranch = Get-CommandOutputOrNull -Command "git" -CommandArgs @("-C", $repoRoot, "rev-parse", "--abbrev-ref", "HEAD")
    $gitStatus = Get-CommandOutputOrNull -Command "git" -CommandArgs @("-C", $repoRoot, "status", "--porcelain")
    $gitDescribe = Get-CommandOutputOrNull -Command "git" -CommandArgs @("-C", $repoRoot, "describe", "--always", "--dirty")

    $buildType = Get-CMakeCacheValue -CachePath $cachePath -Key "CMAKE_BUILD_TYPE"
    if ([string]::IsNullOrWhiteSpace($buildType)) {
        $buildType = "MultiConfig"
    }
    $cxxCompiler = Get-CMakeCacheValue -CachePath $cachePath -Key "CMAKE_CXX_COMPILER"
    $linkerPath = Get-CMakeCacheValue -CachePath $cachePath -Key "CMAKE_LINKER"
    if ([string]::IsNullOrWhiteSpace($cxxCompiler)) {
        $cxxCompiler = $linkerPath
    }
    $cxxCompilerVersion = Get-CMakeCacheValue -CachePath $cachePath -Key "CMAKE_CXX_COMPILER_VERSION"
    if ([string]::IsNullOrWhiteSpace($cxxCompilerVersion) -and -not [string]::IsNullOrWhiteSpace($cxxCompiler)) {
        $m = [regex]::Match($cxxCompiler, "MSVC[\\/](\d+\.\d+\.\d+)")
        if ($m.Success) {
            $cxxCompilerVersion = $m.Groups[1].Value
        }
    }

    $artifacts = New-Object 'System.Collections.Generic.List[object]'
    Add-ArtifactIfExists -List $artifacts -Path $meshCsv
    Add-ArtifactIfExists -List $artifacts -Path $sdfSdfCsv
    if ($RunP2) {
        Add-ArtifactIfExists -List $artifacts -Path $p2SceneCsv
        Add-ArtifactIfExists -List $artifacts -Path $p2ScaleCsv
    }
    if ($RunM3) {
        Add-ArtifactIfExists -List $artifacts -Path $m3Csv
        Add-ArtifactIfExists -List $artifacts -Path $m3ParetoCsv
    }
    if ($GenerateSummary) {
        Add-ArtifactIfExists -List $artifacts -Path $summaryMd
    }
    if ($RunM3 -and $ValidateP1) {
        Add-ArtifactIfExists -List $artifacts -Path $p1ValidationMd
    }

    $plotCandidates = @(
        (Join-Path $outDir "mesh_compare_summary.png"),
        (Join-Path $outDir "sdf_sdf_stability_summary.png"),
        (Join-Path $outDir "m3_octree_pareto.png"),
        (Join-Path $outDir "p2_scene_penetration.png"),
        (Join-Path $outDir "p2_scene_walltime.png"),
        (Join-Path $outDir "p2_scaleup_summary.png")
    )
    foreach ($plot in $plotCandidates) {
        Add-ArtifactIfExists -List $artifacts -Path $plot
    }

    $manifest = [ordered]@{
        generated_utc = (Get-Date).ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ssZ")
        repo_root = $repoRoot
        output_dir = $outDir
        bin_dir = $binDir
        git = [ordered]@{
            commit = $gitCommit
            describe = $gitDescribe
            branch = $gitBranch
            dirty = [bool](-not [string]::IsNullOrWhiteSpace($gitStatus))
        }
        build = [ordered]@{
            cmake_cache = $cachePath
            generator = Get-CMakeCacheValue -CachePath $cachePath -Key "CMAKE_GENERATOR"
            build_type = $buildType
            cxx_compiler = $cxxCompiler
            cxx_compiler_version = $cxxCompilerVersion
            cxx_flags_release = Get-CMakeCacheValue -CachePath $cachePath -Key "CMAKE_CXX_FLAGS_RELEASE"
            linker = $linkerPath
            chrono_dir = Get-CMakeCacheValue -CachePath $cachePath -Key "Chrono_DIR"
        }
        system = [ordered]@{
            machine_name = $env:COMPUTERNAME
            os_caption = $osCaption
            os_version = $osVersion
            cpu_name = $cpuName
            logical_cores = $logicalCores
            total_memory_bytes = $memoryBytes
        }
        run_config = [ordered]@{
            NumSeeds = $NumSeeds
            BaseSeed = $BaseSeed
            MeshSteps = $MeshSteps
            MeshDt = $MeshDt
            SdfSdfDuration = $SdfSdfDuration
            SdfSdfDts = $SdfSdfDts
            RunP2 = $RunP2
            P2SceneDuration = $P2SceneDuration
            P2SceneDts = $P2SceneDts
            P2ScaleDuration = $P2ScaleDuration
            P2ScaleDt = $P2ScaleDt
            P2ScaleCounts = $P2ScaleCounts
            RunM3 = $RunM3
            M3EpsList = $M3EpsList
            M3Queries = $M3Queries
            M3MaxDepth = $M3MaxDepth
            M3Kn = $M3Kn
            M3CError = $M3CError
            M3Band = $M3Band
            M3HalfExtent = $M3HalfExtent
            M3Radius = $M3Radius
            ValidateP1 = $ValidateP1
            AllowP2ValidationFailure = $AllowP2ValidationFailure
            GeneratePlots = $GeneratePlots
            GenerateSummary = $GenerateSummary
        }
        artifacts = $artifacts
    }

    $manifest | ConvertTo-Json -Depth 8 | Set-Content -Path $manifestJson -Encoding UTF8

    $mdLines = @()
    $mdLines += "# Artifact Manifest"
    $mdLines += ""
    $mdLines += "- Generated (UTC): $($manifest.generated_utc)"
    $mdLines += "- Repo: $($manifest.repo_root)"
    $mdLines += "- Output: $($manifest.output_dir)"
    $mdLines += "- Git commit: $($manifest.git.commit)"
    $mdLines += "- Git describe: $($manifest.git.describe)"
    $mdLines += "- Git branch: $($manifest.git.branch)"
    $mdLines += "- Dirty worktree: $($manifest.git.dirty)"
    $mdLines += "- Compiler: $($manifest.build.cxx_compiler)"
    $mdLines += "- Compiler version: $($manifest.build.cxx_compiler_version)"
    $mdLines += "- Build type: $($manifest.build.build_type)"
    $mdLines += "- CPU: $($manifest.system.cpu_name)"
    $mdLines += "- Logical cores: $($manifest.system.logical_cores)"
    $mdLines += "- Total memory (bytes): $($manifest.system.total_memory_bytes)"
    $mdLines += ""
    $mdLines += "## Artifacts"
    $mdLines += ""
    $mdLines += "| file | size_bytes | sha256 |"
    $mdLines += "| --- | ---: | --- |"
    foreach ($a in $artifacts) {
        $mdLines += "| $($a.name) | $($a.size_bytes) | $($a.sha256) |"
    }
    Set-Content -Path $manifestMd -Value ($mdLines -join "`n") -Encoding UTF8
}

if ($PackageArtifacts) {
    $bundlePath = Join-Path $outDir "artifact_bundle.zip"
    if (Test-Path $bundlePath) {
        Remove-Item -Path $bundlePath -Force
    }
    Compress-Archive -Path (Join-Path $outDir "*") -DestinationPath $bundlePath -CompressionLevel Optimal
}

Write-Host "Finished."
Write-Host "Seeds: $NumSeeds (base=$BaseSeed)"
Write-Host "Mesh CSV: $meshCsv"
Write-Host "SDF-SDF CSV: $sdfSdfCsv"
if ($RunP2) {
    Write-Host "P2 Scene CSV: $p2SceneCsv"
    Write-Host "P2 Scale CSV: $p2ScaleCsv"
}
if ($RunM3) {
    Write-Host "M3 CSV: $m3Csv"
    Write-Host "M3 Pareto CSV: $m3ParetoCsv"
}
if ($GenerateSummary) {
    Write-Host "Summary MD: $summaryMd"
}
if ($RunM3 -and $ValidateP1) {
    Write-Host "P1 Validation MD: $p1ValidationMd"
}
if ($WriteMetadata) {
    Write-Host "Manifest JSON: $(Join-Path $outDir 'artifact_manifest.json')"
    Write-Host "Manifest MD: $(Join-Path $outDir 'artifact_manifest.md')"
}
if ($PackageArtifacts) {
    Write-Host "Artifact Bundle: $(Join-Path $outDir 'artifact_bundle.zip')"
}
