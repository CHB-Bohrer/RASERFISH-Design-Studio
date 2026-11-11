#!/bin/bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
cd "$HERE"
if bash "$HERE/install.sh"; then
    exit 0
else
    status=$?
    echo "Stopped with exit $status. Keep the error text above."
    read -r -p "Press Return to close this window. " answer || true
    exit "$status"
fi
