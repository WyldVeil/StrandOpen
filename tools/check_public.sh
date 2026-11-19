#!/bin/bash
# Before publishing: nothing in the tracked files (or a release folder) may
# point at the author, their data or this machine. The patterns live outside
# the repository (one per line in ~/.config/strandopen-private-patterns), so
# the check itself reveals nothing. Binary files are read too (a compiled
# program or a Python cache can carry a path), with the patterns matched as
# whole words there (a short pattern turns up by chance inside machine code
# and symbol names), and every file is searched for the home folder's path.
#   tools/check_public.sh [FOLDER]     (default: the tracked files)
set -uo pipefail
cd "$(dirname "$0")/.."
PAT="${STRANDOPEN_PRIVATE_PATTERNS:-$HOME/.config/strandopen-private-patterns}"
[ -s "$PAT" ] || { echo "no pattern file at $PAT"; exit 2; }
if [ $# -gt 0 ]; then
    hits=$( { grep -rIil -f "$PAT" "$1"; grep -railw -f "$PAT" "$1"; grep -rlaF -- "$HOME" "$1"; } 2>/dev/null | sort -u)
else
    hits=$( { git grep -Iil -f "$PAT" -- .; git grep -ailw -f "$PAT" -- .; git grep -laF -- "$HOME" -- .; } 2>/dev/null | sort -u)
fi
if [ -n "$hits" ]; then
    echo "private details found in:"; echo "$hits"
    exit 1
fi
echo "clean: no private patterns found"
