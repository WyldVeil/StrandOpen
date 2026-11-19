#!/bin/bash
# Before publishing: nothing in the tracked files (or a release folder) may
# point at the author, their data or this machine. The patterns live outside
# the repository (one per line in ~/.config/strandopen-private-patterns), so
# the check itself reveals nothing.
#   tools/check_public.sh [FOLDER]     (default: the tracked files)
set -uo pipefail
cd "$(dirname "$0")/.."
PAT="${STRANDOPEN_PRIVATE_PATTERNS:-$HOME/.config/strandopen-private-patterns}"
[ -s "$PAT" ] || { echo "no pattern file at $PAT"; exit 2; }
if [ $# -gt 0 ]; then
    hits=$(grep -rIil -f "$PAT" "$1" 2>/dev/null)
else
    hits=$(git grep -Iil -f "$PAT" -- . 2>/dev/null)
fi
if [ -n "$hits" ]; then
    echo "private details found in:"; echo "$hits"
    exit 1
fi
echo "clean: no private patterns found"
