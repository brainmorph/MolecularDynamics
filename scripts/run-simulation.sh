#!/usr/bin/env bash
# Run the core LJ simulation from a settings file. Writes the trajectory and log next to the
# paths named in the settings (default: out/ in the repo root).
#
# Usage: scripts/run-simulation.sh [settings.json]
#   settings.json  run settings (default: examples/argon-8.json)
#
# Environment:
#   CONFIGURATION  build configuration, Release (default) or Debug
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
settings="${1:-$repo/examples/argon-8.json}"
configuration="${CONFIGURATION:-Release}"

if [[ ! -f "$settings" ]]; then
    echo "error: no settings file '$settings'." >&2
    exit 1
fi

dotnet build "$repo/src/MD.Cli" -c "$configuration" -nologo -v quiet >/dev/null || {
    echo "error: build failed. Run 'dotnet build src/MD.Cli' to see the errors." >&2
    exit 1
}

exec dotnet run -c "$configuration" --no-build --project "$repo/src/MD.Cli" -- run "$settings"
