#!/usr/bin/env bash
# Open the most recently written trajectory in the Godot viewer.
#
# Usage: scripts/view-latest.sh [search-dir]
#   search-dir  where to look for *.traj files (default: out/ in the repo root)
#
# Environment:
#   GODOT       path to the Godot 4.7 Mono binary (default: the local install on the Desktop)
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
search="${1:-$repo/out}"
godot="${GODOT:-/home/purple-eden/Desktop/Godot_v4.7-stable_mono_linux_x86_64/Godot_v4.7-stable_mono_linux.x86_64}"

if [[ ! -d "$search" ]]; then
    echo "error: no directory '$search'. Run a simulation first, e.g.:" >&2
    echo "  dotnet run -c Release --project src/MD.Cli -- run examples/argon-8.json" >&2
    exit 1
fi
if [[ ! -x "$godot" ]]; then
    echo "error: Godot not found at '$godot'. Set GODOT to the Godot 4.7 Mono binary." >&2
    exit 1
fi

# Newest *.traj by modification time.
latest="$(find "$search" -type f -name '*.traj' -printf '%T@ %p\n' | sort -n | tail -n 1 | cut -d' ' -f2-)"
if [[ -z "$latest" ]]; then
    echo "error: no .traj files under '$search'." >&2
    exit 1
fi

echo "Trajectory: $latest ($(date -r "$latest" '+%Y-%m-%d %H:%M:%S'))"
echo "Building viewer..."
dotnet build "$repo/viewer" -nologo -v quiet >/dev/null || {
    echo "error: viewer build failed. Run 'dotnet build viewer' to see the errors." >&2
    exit 1
}

exec "$godot" --path "$repo/viewer" -- "$latest"
