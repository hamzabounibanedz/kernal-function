Param(
    [ValidateSet('Release', 'Debug')]
    [string]$Configuration = 'Release',
    [switch]$Console
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root

function Write-Step($msg) { Write-Host "==> $msg" -ForegroundColor Cyan }
function Write-Info($msg) { Write-Host "    $msg" -ForegroundColor DarkGray }

function Sanitize-Env {
    Write-Step "Sanitizing environment (removing MSYS2/MinGW from build)"
    $patterns = @('\msys64\', '/msys64/', '\mingw64\', '/mingw64/')
    foreach ($name in 'PATH', 'INCLUDE', 'LIB', 'LIBPATH', 'C_INCLUDE_PATH', 'CPLUS_INCLUDE_PATH', 'CPATH') {
        $val = [Environment]::GetEnvironmentVariable($name, 'Process')
        if (-not $val) { continue }
        $parts = $val -split ';' | Where-Object { $_ -and $_.Trim() -ne '' }
        $filtered = @()
        foreach ($p in $parts) {
            $bad = $false
            foreach ($pat in $patterns) { if ($p -like "*${pat}*") { $bad = $true; break } }
            if (-not $bad) { $filtered += $p }
        }
        [Environment]::SetEnvironmentVariable($name, ($filtered -join ';'), 'Process')
    }
    foreach ($n in 'MSYSTEM', 'MINGW_PREFIX', 'MSYS2_PATH_TYPE', 'CHERE_INVOKING', 'PKG_CONFIG', 'PKG_CONFIG_PATH') {
        [Environment]::SetEnvironmentVariable($n, $null, 'Process')
    }
}

Sanitize-Env

function Find-CMake {
    $cmd = Get-Command cmake.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    $candidates = @(
        'C:\\Program Files\\CMake\\bin\\cmake.exe',
        'C:\\Program Files\\Kitware\\CMake\\bin\\cmake.exe',
        'C:\\Program Files (x86)\\CMake\\bin\\cmake.exe'
    )
    foreach ($c in $candidates) { if (Test-Path $c) { return $c } }
    $vsRoots = @('C:\\Program Files\\Microsoft Visual Studio\\2022', 'C:\\Program Files (x86)\\Microsoft Visual Studio\\2022', 'C:\\VS2022-BT')
    foreach ($rootVS in $vsRoots) {
        if (Test-Path $rootVS) {
            $found = Get-ChildItem -Path $rootVS -Recurse -Filter cmake.exe -ErrorAction SilentlyContinue |
            Where-Object { $_.FullName -match "\\CMake\\bin\\cmake\\.exe$" } |
            Select-Object -First 1
            if ($found) { return $found.FullName }
        }
    }
    return $null
}

$CMAKE = Find-CMake
if (-not $CMAKE) {
    Write-Host "ERROR: CMake not found in PATH or common locations." -ForegroundColor Red
    Write-Host "Install it once, then rerun this script:" -ForegroundColor Yellow
    Write-Host "    winget install --accept-package-agreements --accept-source-agreements Kitware.CMake" -ForegroundColor DarkGray
    exit 1
}

function Find-Windeployqt {
    $cmd = Get-Command windeployqt.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    if ($env:QTDIR) {
        $cand = Join-Path $env:QTDIR 'bin/windeployqt.exe'
        if (Test-Path $cand) { return $cand }
    }
    if ($env:QT_DIR) {
        $cand = Join-Path $env:QT_DIR 'bin/windeployqt.exe'
        if (Test-Path $cand) { return $cand }
    }
    $localQt = Join-Path $root '.qt-runtime'
    if (Test-Path $localQt) {
        $found = Get-ChildItem -Path $localQt -Filter windeployqt.exe -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($found) { return $found.FullName }
    }
    $qtRoot = 'C:\Qt'
    if (Test-Path $qtRoot) {
        $found = Get-ChildItem -Path $qtRoot -Filter windeployqt.exe -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($found) { return $found.FullName }
    }
    return $null
}

function Find-QtPlatformDll {
    $candidates = @()
    if ($env:QT_PLUGIN_PATH) {
        $candidates += Join-Path $env:QT_PLUGIN_PATH 'platforms/qwindows.dll'
    }
    if ($env:QTDIR) {
        $candidates += Join-Path $env:QTDIR 'plugins/platforms/qwindows.dll'
    }
    if ($env:QT_DIR) {
        $candidates += Join-Path $env:QT_DIR 'plugins/platforms/qwindows.dll'
    }
    foreach ($p in $candidates) { if (Test-Path $p) { return $p } }

    $localQt = Join-Path $root '.qt-runtime'
    if (Test-Path $localQt) {
        $hit = Get-ChildItem -Path $localQt -Recurse -Filter qwindows.dll -ErrorAction SilentlyContinue |
        Select-Object -First 1
        if ($hit) { return $hit.FullName }
    }

    $roots = @('C:\Qt', 'C:\Program Files\Qt', 'C:\Program Files (x86)\Qt')
    foreach ($r in $roots) {
        if (Test-Path $r) {
            $hit = Get-ChildItem -Path $r -Recurse -Filter qwindows.dll -ErrorAction SilentlyContinue |
            Select-Object -First 1
            if ($hit) { return $hit.FullName }
        }
    }
    return $null
}

function Find-QtPluginsRoot {
    $candidates = @()
    if ($env:QT_PLUGIN_PATH) { $candidates += $env:QT_PLUGIN_PATH }
    if ($env:QTDIR) { $candidates += (Join-Path $env:QTDIR 'plugins') }
    if ($env:QT_DIR) { $candidates += (Join-Path $env:QT_DIR 'plugins') }
    foreach ($p in $candidates) {
        if (Test-Path (Join-Path $p 'platforms/qwindows.dll')) { return $p }
    }
    $localQt = Join-Path $root '.qt-runtime'
    if (Test-Path $localQt) {
        $hit = Get-ChildItem -Path $localQt -Recurse -Directory -Filter platforms -ErrorAction SilentlyContinue |
        Where-Object { Test-Path (Join-Path $_.FullName 'qwindows.dll') } |
        Select-Object -First 1
        if ($hit) { return Split-Path $hit.FullName -Parent }
    }
    $roots = @('C:\Qt', 'C:\Program Files\Qt', 'C:\Program Files (x86)\Qt')
    foreach ($r in $roots) {
        if (Test-Path $r) {
            $hit = Get-ChildItem -Path $r -Recurse -Directory -Filter platforms -ErrorAction SilentlyContinue |
            Where-Object { Test-Path (Join-Path $_.FullName 'qwindows.dll') } |
            Select-Object -First 1
            if ($hit) { return Split-Path $hit.FullName -Parent }
        }
    }
    return $null
}

function Sync-QtRuntime {
    param([Parameter(Mandatory = $true)][string]$PluginsRoot,
        [Parameter(Mandatory = $true)][string]$TargetBin)
    Write-Info "Synchronizing Qt runtime into $TargetBin"
    $qtBin = Join-Path (Split-Path $PluginsRoot -Parent) 'bin'
    if (Test-Path $qtBin) {
        # Remove possibly mismatched Qt6*.dll in app bin
        Get-ChildItem -Path $TargetBin -Filter 'Qt6*.dll' -ErrorAction SilentlyContinue | ForEach-Object {
            try { Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue } catch {}
        }
        foreach ($d in 'libEGL.dll', 'libGLESv2.dll', 'd3dcompiler_47.dll', 'opengl32sw.dll') {
            $dst = Join-Path $TargetBin $d
            if (Test-Path $dst) { try { Remove-Item $dst -Force -ErrorAction SilentlyContinue } catch {} }
        }
        # Copy matching Qt6*.dll into app bin
        Get-ChildItem -Path $qtBin -Filter 'Qt6*.dll' -ErrorAction SilentlyContinue | Copy-Item -Destination $TargetBin -Force
        foreach ($d in 'libEGL.dll', 'libGLESv2.dll', 'd3dcompiler_47.dll', 'opengl32sw.dll') {
            $src = Join-Path $qtBin $d; if (Test-Path $src) { Copy-Item $src $TargetBin -Force }
        }
    }
    # Replace plugin subfolders with those from PluginsRoot
    foreach ($sub in 'platforms', 'styles', 'imageformats', 'iconengines') {
        $src = Join-Path $PluginsRoot $sub
        $dst = Join-Path $TargetBin $sub
        if (Test-Path $dst) { try { Remove-Item $dst -Recurse -Force -ErrorAction SilentlyContinue } catch {} }
        if (Test-Path $src) {
            New-Item -ItemType Directory -Force -Path $dst | Out-Null
            Copy-Item (Join-Path $src '*') $dst -Recurse -Force -ErrorAction SilentlyContinue
        }
    }
}

function Find-QtDll {
    param([Parameter(Mandatory = $true)][string]$Name)
    $candidates = @()
    if ($env:QTDIR) { $candidates += Join-Path $env:QTDIR     ("bin/" + $Name) }
    if ($env:QT_DIR) { $candidates += Join-Path $env:QT_DIR    ("bin/" + $Name) }
    foreach ($p in $candidates) { if (Test-Path $p) { return $p } }

    $localQt = Join-Path $root '.qt-runtime'
    if (Test-Path $localQt) {
        $hit = Get-ChildItem -Path $localQt -Recurse -Filter $Name -ErrorAction SilentlyContinue |
        Where-Object { $_.FullName -match "\\bin\\$([regex]::Escape($Name))$" } |
        Select-Object -First 1
        if ($hit) { return $hit.FullName }
    }

    $roots = @('C:\\Qt', 'C:\\Program Files\\Qt', 'C:\\Program Files (x86)\\Qt')
    foreach ($r in $roots) {
        if (Test-Path $r) {
            $hit = Get-ChildItem -Path $r -Recurse -Filter $Name -ErrorAction SilentlyContinue |
            Where-Object { $_.FullName -match "\\bin\\$([regex]::Escape($Name))$" } |
            Select-Object -First 1
            if ($hit) { return $hit.FullName }
        }
    }
    return $null
}

# Configure once (idempotent)
$qtPrefix = Join-Path $root ".qt-runtime/6.7.3/msvc2019_64"
if (Test-Path $qtPrefix) {
    $qtCMake = Join-Path $qtPrefix 'lib/cmake/Qt6'
    if (Test-Path "$root/build_run") {
        Write-Step "Cleaning existing build_run to avoid mixed Qt versions"
        try { Remove-Item -Recurse -Force "$root/build_run" } catch {}
    }
    Write-Step "Configuring with local Qt: $qtPrefix"
    & $CMAKE -S "$root" -B "$root/build_run" -G "Visual Studio 17 2022" -A x64 -D CMAKE_PREFIX_PATH="$qtPrefix" -D Qt6_DIR="$qtCMake"
}
else {
    if (!(Test-Path "$root/build_run/CMakeCache.txt")) {
        Write-Step "Configuring (preset: windows-msvc)"
        & $CMAKE --preset windows-msvc
    }
    else {
        Write-Step "Using existing CMake cache in build_run"
    }
}

# Build
if ($Configuration -eq 'Debug') {
    Write-Step "Building Debug"
    & $CMAKE --build --preset windows-debug
    $bin = Join-Path $root 'build_run/Debug'
}
else {
    Write-Step "Building Release"
    & $CMAKE --build --preset windows-release
    $bin = Join-Path $root 'build_run/Release'
}

# Paths
$exe = if ($Console) { 'kernel_console.exe' } else { 'kernel_ui.exe' }
$path = Join-Path $bin $exe
if (!(Test-Path $path)) { throw "Executable not found: $path" }

# Ensure Qt platform plugins (GUI only)
if (-not $Console) {
    $platformDll = Join-Path $bin 'platforms/qwindows.dll'
    if (!(Test-Path $platformDll)) {
        Write-Step "Qt platform plugin missing"
        $wdq = Find-Windeployqt
        if ($wdq) {
            Write-Info "Using windeployqt at $wdq"
            $mode = if ($Configuration -eq 'Debug') { '--debug' } else { '--release' }
            try { & $wdq $mode --verbose 0 "$path" } catch { Write-Info "windeployqt failed, attempting manual copy" }
        }
        else {
            Write-Info "windeployqt not found, attempting manual copy"
        }
        if (!(Test-Path $platformDll)) {
            $plat = Find-QtPlatformDll
            if ($plat) {
                Write-Info "Found platform plugin: $plat"
                $dir = Split-Path $platformDll -Parent
                New-Item -ItemType Directory -Force -Path $dir | Out-Null
                Copy-Item $plat $dir -Force
            }
        }
        if (!(Test-Path $platformDll)) {
            # Copy entire platforms folder and common plugin folders from plugins root
            $pluginsRoot = Find-QtPluginsRoot
            if ($pluginsRoot) {
                Write-Info "Copying plugins from $pluginsRoot"
                foreach ($sub in 'platforms', 'styles', 'imageformats', 'iconengines') {
                    $src = Join-Path $pluginsRoot $sub
                    if (Test-Path $src) {
                        $dst = Join-Path $bin $sub
                        New-Item -ItemType Directory -Force -Path $dst | Out-Null
                        Copy-Item (Join-Path $src '*') $dst -Recurse -Force -ErrorAction SilentlyContinue
                    }
                }
                # Also ensure core Qt DLLs match the plugin's Qt by copying from the sibling bin
                $qtBin = Join-Path (Split-Path $pluginsRoot -Parent) 'bin'
                if (Test-Path $qtBin) {
                    $dlls = @('Qt6Core.dll', 'Qt6Gui.dll', 'Qt6Widgets.dll', 'Qt6Concurrent.dll', 'Qt6OpenGL.dll', 'Qt6OpenGLWidgets.dll', 'd3dcompiler_47.dll', 'opengl32sw.dll', 'libEGL.dll', 'libGLESv2.dll')
                    foreach ($d in $dlls) {
                        $srcDll = Join-Path $qtBin $d
                        if (Test-Path $srcDll) { Copy-Item $srcDll $bin -Force }
                    }
                }
            }
        }
        if (!(Test-Path $platformDll)) { throw "Qt platform plugin still missing: $platformDll" }
    }

    # Ensure QtCharts DLL if the binary links it
    $chartsName = if ($Configuration -eq 'Debug') { 'Qt6Chartsd.dll' } else { 'Qt6Charts.dll' }
    $chartsBin = Join-Path $bin $chartsName
    if (!(Test-Path $chartsBin)) {
        $chartsSrc = Find-QtDll -Name $chartsName
        if ($chartsSrc) {
            Write-Info "Copying $chartsName from $chartsSrc"
            Copy-Item $chartsSrc $bin -Force
        }
        else {
            Write-Info "$chartsName not found (ok if built without Charts). If needed, install via aqt:"
            Write-Info "python -m aqt install-qt windows desktop 6.7.3 win64_msvc2019_64 -O .qt-runtime --modules qtcharts"
        }
    }
}

# Launch
if ($Console) {
    Write-Step "Running console app"
    & $path
}
else {
    # Ensure Qt loads plugins from the deployed directory or local Qt
    $pluginsRoot = Find-QtPluginsRoot
    if ($pluginsRoot) {
        Sync-QtRuntime -PluginsRoot $pluginsRoot -TargetBin $bin
        $qtBin = Join-Path (Split-Path $pluginsRoot -Parent) 'bin'
        if (Test-Path $qtBin) { $env:PATH = "$qtBin;$env:PATH" }
        # Prefer plugins next to the exe after sync
        $env:QT_PLUGIN_PATH = $bin
        $env:QT_QPA_PLATFORM_PLUGIN_PATH = (Join-Path $bin 'platforms')
    }
    else {
        $env:QT_QPA_PLATFORM_PLUGIN_PATH = (Join-Path $bin 'platforms')
    }
    $env:QT_DEBUG_PLUGINS = '1'
    Write-Step "Running GUI (attached)"
    & $path
}


