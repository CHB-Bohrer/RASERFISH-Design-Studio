#!/bin/bash
# RASER-FISH Design Studio 1.6: fresh all-in-one install into the folder containing this script.
# No previous Codex installation is needed. No scientific code or inputs are edited.
# Run from the directory containing studio.py: bash install.sh
set -eEuo pipefail
OPEN_AFTER=1
DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --no-open) OPEN_AFTER=0 ;;
        --dry-run) DRY_RUN=1 ;;
        --help|-h)
            echo 'Run from the new complete Studio folder: bash install.sh'
            echo 'Installs local software, downloads/indexes local genome + RNA, tests, then opens Studio.'
            echo '  --no-open   Install and verify without starting the browser server.'
            echo '  --dry-run   Print destinations only. No files are changed or downloaded.'
            exit 0 ;;
        *) echo "Unknown option: $arg" >&2; exit 2 ;;
    esac
done
ROOT="$(cd "$(dirname "$0")" && pwd -P)"
cd "$ROOT"
for file in studio.py design_engine.py sequence_tools.py engine/qc_core.py start_studio.sh setup/studio_install_helper.py; do
    if [ ! -f "$ROOT/$file" ]; then
        echo "STOP: $file was not found in the current folder." >&2
        echo 'Keep the complete extracted folder together and run its install.sh.' >&2
        exit 2
    fi
done
case "$ROOT" in *[[:space:]]*) echo 'STOP: The installation path must not contain spaces or other whitespace.' >&2; exit 2;; esac
for folder in .qc_env .tools references; do
    if [ -L "$ROOT/$folder" ]; then
        echo "STOP: $folder is a symbolic link. This installer expects actual folders here." >&2
        echo 'No linked files were changed. Review the existing shortcut before proceeding.' >&2
        exit 2
    fi
done
if [ "$DRY_RUN" -eq 1 ]; then
    echo "Studio stays at: $ROOT"
    echo "Software environment: $ROOT/.qc_env"
    echo "Private package manager/cache: $ROOT/.tools"
    echo "Human DNA + RNA: $ROOT/references/GCF_000001405.40"
    echo "New design results: $ROOT/design_runs"
    echo 'No link to an old installation. No old files are deleted. No probe design is run.'
    exit 0
fi
case "$(uname -s)" in Darwin) OS=MacOSX;; Linux) OS=Linux;; *) echo 'STOP: This installer supports macOS and Linux.' >&2; exit 2;; esac
ARCH="$(uname -m)"
case "$OS:$ARCH" in MacOSX:arm64|MacOSX:x86_64|Linux:x86_64|Linux:aarch64) ;; *) echo "Unsupported platform: $OS $ARCH" >&2; exit 2;; esac
if [ "$OS" = MacOSX ] && [ "$(sysctl -n sysctl.proc_translated 2>/dev/null || true)" = 1 ]; then
    echo 'STOP: Terminal is running under Rosetta. Use native Apple Silicon Terminal for this new installation.' >&2
    exit 2
fi
if ! command -v curl >/dev/null 2>&1; then echo 'STOP: curl is needed to download the installer.' >&2; exit 2; fi
FREE_KB="$(df -Pk "$ROOT" | awk 'END {print $4}')"
case "$FREE_KB" in ''|*[!0-9]*) echo 'STOP: Could not determine available disk space.' >&2; exit 2;; esac
if [ ! -f "$ROOT/references/GCF_000001405.40/reference_manifest.json" ] && [ "$FREE_KB" -lt 15728640 ]; then
    echo 'STOP: Less than 15 GiB free. Plan approximately 30 GB of free space for software, references and working files.' >&2
    exit 2
fi
mkdir -p "$ROOT/.tools/studio_setup" "$ROOT/setup_logs"
LOCK="$ROOT/.tools/studio_setup/install.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
    echo "STOP: Another or interrupted installation has a lock: $LOCK" >&2
    echo 'Check for an active installer before removing that empty lock directory.' >&2
    exit 2
fi
trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
STAMP="$(date -u +%Y%m%dT%H%M%SZ)_$$"
LOGDIR="$ROOT/setup_logs/install_here_$STAMP"
mkdir -p "$LOGDIR"
exec > >(tee -a "$LOGDIR/INSTALL.log") 2>&1
trap 'echo "INSTALLATION STOPPED at line $LINENO. Keep $LOGDIR/INSTALL.log. Do not delete or move your old installation." >&2' ERR
ENV="$ROOT/.qc_env"
MF="$ROOT/.tools/miniforge"
export CONDA_PKGS_DIRS="$ROOT/.tools/conda_pkgs"
export CONDA_ENVS_PATH="$ROOT/.tools/conda_envs"
export PIP_CACHE_DIR="$ROOT/.tools/pip_cache"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export PIP_DISABLE_PIP_VERSION_CHECK=1
export PIP_CONFIG_FILE=/dev/null
unset PYTHONPATH PYTHONHOME || true
mkdir -p "$CONDA_PKGS_DIRS" "$CONDA_ENVS_PATH" "$PIP_CACHE_DIR"
LOCATION_FILE="$ROOT/.tools/studio_setup/environment_location.txt"
if [ -f "$LOCATION_FILE" ] && [ "$(cat "$LOCATION_FILE")" != "$ROOT" ]; then
    echo 'STOP: This installation folder has moved since setup. Reinstall at its final location rather than using broken absolute paths.' >&2
    exit 2
fi
printf '%s\n' "$ROOT" > "$LOCATION_FILE"
echo 'RASER-FISH: one-folder installation'
echo "Destination: $ROOT"
echo 'This downloads software and public human DNA/RNA references into this folder.'
echo 'It does NOT run an experimental-target design, change QC thresholds, or replace existing probes.'
echo 'Allow about 30 GB of free space as a planning estimate; actual size varies.'
echo ''
echo '[1/5] Install the private package manager.'
if [ ! -x "$MF/bin/conda" ]; then
    if [ -e "$MF" ]; then echo "STOP: An incomplete Miniforge folder exists: $MF. Keep the log and review it before retrying." >&2; exit 2; fi
    VERSION='26.7.2-0'
    FILE="Miniforge3-${VERSION}-${OS}-${ARCH}.sh"
    URL="https://github.com/conda-forge/miniforge/releases/download/${VERSION}"
    curl --fail --location --retry 3 --connect-timeout 20 --output "$ROOT/.tools/$FILE" "$URL/$FILE"
    curl --fail --location --retry 3 --connect-timeout 20 --output "$ROOT/.tools/$FILE.sha256" "$URL/$FILE.sha256"
    EXPECTED="$(awk 'NR==1 {print $1}' "$ROOT/.tools/$FILE.sha256")"
    if command -v shasum >/dev/null 2>&1; then
        ACTUAL="$(shasum -a 256 "$ROOT/.tools/$FILE" | awk '{print $1}')"
    else
        ACTUAL="$(sha256sum "$ROOT/.tools/$FILE" | awk '{print $1}')"
    fi
    case "$EXPECTED" in ''|*[!0-9a-fA-F]*) echo 'STOP: Invalid installer checksum record.' >&2; exit 2;; esac
    if [ "${#EXPECTED}" -ne 64 ] || [ "$EXPECTED" != "$ACTUAL" ]; then
        echo 'STOP: Installer checksum mismatch. Nothing from this download will be executed.' >&2
        exit 2
    fi
    cp "$ROOT/.tools/$FILE.sha256" "$LOGDIR/miniforge_installer.sha256"
    bash "$ROOT/.tools/$FILE" -b -p "$MF"
fi

echo '[2/5] Install the local Python, BLAST, NCBI Datasets and DNA-folding software.'
if [ ! -f "$ENV/conda-meta/history" ]; then
    if [ -e "$ENV" ]; then
        echo "STOP: An incomplete/non-Conda .qc_env already exists at $ENV. Nothing was overwritten." >&2
        exit 2
    fi
    "$MF/bin/conda" create -y -p "$ENV" --override-channels \
        -c conda-forge -c bioconda --strict-channel-priority \
        'python=3.11' 'biopython=1.86' 'numpy>=1.26,<3' pip blast ncbi-datasets-cli
else
    echo 'A local Conda environment is already present. Ensuring the required packages are installed.'
    "$MF/bin/conda" install -y -p "$ENV" --override-channels \
        -c conda-forge -c bioconda --strict-channel-priority \
        'python=3.11' 'biopython=1.86' 'numpy>=1.26,<3' pip blast ncbi-datasets-cli
fi
export PATH="$ENV/bin:$PATH"
"$ENV/bin/python" -m pip install --only-binary=:all: 'ViennaRNA==2.7.2'
"$MF/bin/conda" list -p "$ENV" --explicit > "$LOGDIR/conda_explicit.txt"
"$ENV/bin/python" -m pip freeze > "$LOGDIR/pip_freeze.txt"

# Keep the inspectable Python helper local; no scientific source files are edited.

HELPER="$ROOT/setup/studio_install_helper.py"
"$ENV/bin/python" "$HELPER" snapshot --root "$ROOT" --run "$LOGDIR"
echo '[3/5] Test the actual Studio folding and BLAST engines.'
"$ENV/bin/python" "$HELPER" engines --root "$ROOT" --run "$LOGDIR"
echo 'Running all bundled Python regression tests before the large reference download.'
"$ENV/bin/python" -m unittest discover -s "$ROOT/tests" -v > "$LOGDIR/REGRESSION_TESTS.txt" 2>&1
tail -n 6 "$LOGDIR/REGRESSION_TESTS.txt"
echo '[4/5] Download and index human DNA and annotated RNA inside this Studio folder.'
"$ENV/bin/python" "$HELPER" references --root "$ROOT" --run "$LOGDIR"
echo '[5/5] Verify Studio finds its local software/references and all protected inputs are unchanged.'
"$ENV/bin/python" "$HELPER" finish --root "$ROOT" --run "$LOGDIR"

printf '\nInstallation finished. Software and both search databases are HERE:\n%s\n\n' "$ROOT"
echo "Setup record: $LOGDIR/INSTALLATION_READY.json"
echo 'Normal future start (no reinstall or reference download):'
echo '  bash start_studio.sh'
echo 'Or double-click Open_Studio.command in Finder.'
echo 'Do not rename or relocate this installed folder. Keep Dropbox files available offline.'
# Release the installation lock before running the long-lived browser server.
rmdir "$LOCK"
trap - EXIT ERR
if [ "$OPEN_AFTER" -eq 1 ]; then
    echo 'Opening the live designer. Keep this Terminal open; Ctrl-C stops the server.'
    exec bash "$ROOT/start_studio.sh"
fi
