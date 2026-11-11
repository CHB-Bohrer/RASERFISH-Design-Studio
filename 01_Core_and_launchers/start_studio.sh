#!/bin/bash
# Reuse this installation's runtime and references; never install or download.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
cd "$HERE"
QC_ROOT="$HERE"
PORT=8769
NO_OPEN=0
while [ "$#" -gt 0 ]; do
    case "$1" in
        --qc-root)
            [ "$#" -ge 2 ] || { echo 'Missing --qc-root directory.' >&2; exit 2; }
            QC_ROOT="$(cd "$2" && pwd -P)"; shift 2 ;;
        --port)
            [ "$#" -ge 2 ] || { echo 'Missing --port number.' >&2; exit 2; }
            PORT="$2"; shift 2 ;;
        --no-open) NO_OPEN=1; shift ;;
        *) echo 'Usage: bash start_studio.sh [--qc-root /path/to/installed/Studio] [--port 8769] [--no-open]' >&2; exit 2 ;;
    esac
done
if [ ! -x "$QC_ROOT/.qc_env/bin/python" ]; then
    echo 'The existing .qc_env runtime was not found in the selected Studio folder.' >&2
    echo 'Use: bash start_studio.sh --qc-root /path/to/installed/Studio' >&2
    echo 'First-time setup: bash install.sh (only for a new installation).' >&2
    exit 2
fi
LOCATION="$QC_ROOT/.tools/studio_setup/environment_location.txt"
if [ -f "$LOCATION" ] && [ "$(cat "$LOCATION")" != "$QC_ROOT" ]; then
    echo 'The runtime location does not match this folder. Select its existing original location with --qc-root.' >&2
    exit 2
fi
if [ ! -f "$QC_ROOT/references/GCF_000001405.40/reference_manifest.json" ]; then
    echo 'The existing human DNA/RNA reference manifest was not found in the selected folder.' >&2
    exit 2
fi
export PATH="$QC_ROOT/.qc_env/bin:$PATH"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
unset PYTHONPATH PYTHONHOME || true
set -- "$QC_ROOT/.qc_env/bin/python" -u "$HERE/studio_launcher.py" --qc-root "$QC_ROOT" --port "$PORT"
if [ "$NO_OPEN" = 1 ]; then set -- "$@" --no-open; fi
exec "$@"
