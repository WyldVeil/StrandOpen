#!/bin/bash
# Builds the Windows release, on Linux: dist/StrandOpen-<version>-windows-x86_64/
# and its .zip. Everything needed to run is inside (the program, the import,
# the bundled tools), so nobody installs anything: unzip, run StrandOpen.exe.
#   tools/release_windows.sh
set -euo pipefail
cd "$(dirname "$0")/.."
VER=$(grep -o 'ST_VERSION "[^"]*"' src/app.h | cut -d'"' -f2)
NAME="StrandOpen-$VER-windows-x86_64"
D="dist/$NAME"
./build.sh windows
[ -x build/win/runtime/env/python.exe ] || tools/bundle_runtime_windows.sh
rm -rf "$D" "dist/$NAME.zip"
mkdir -p "$D/share"
cp build/windows/StrandOpen.exe "$D/"
rsync -a --exclude '__pycache__/' --exclude '*.pyc' --exclude '*.log' pipeline/ "$D/pipeline/"
rsync -a data/ "$D/data/"
cp -a build/win/runtime "$D/runtime"
find "$D" -type d -name __pycache__ -prune -exec rm -rf {} +
find "$D" -type f -name '*.pyc' -delete
cp README.md LICENSE "$D/"
[ -d docs ] && rsync -a docs/ "$D/docs/"
cp assets/fonts/Inter-LICENSE.txt assets/fonts/JetBrainsMono-OFL.txt "$D/share/"
cat > "$D/THIRD_PARTY.txt" <<'EOF'
Strand Open includes, unchanged or compiled in:
  raylib 6.2 (zlib licence) https://github.com/raysan5/raylib
  SQLite 3.54.1 (public domain) https://sqlite.org
  Inter (SIL Open Font Licence 1.1) https://rsms.me/inter/
  JetBrains Mono (SIL Open Font Licence 1.1) https://www.jetbrains.com/lp/mono/
  Natural Earth 1:50m land outlines (public domain) https://www.naturalearthdata.com/
runtime\ holds the bundled programs (Python, bcftools, samtools, plink2, Java,
busybox and more), each under its own licence: see runtime\LICENSES\SUMMARY.txt.
The reference data downloaded at setup comes from its publishers and stays under
their terms (listed in README.md).
EOF
cat > "$D/START HERE.txt" <<'EOF'
Strand Open for Windows

Run StrandOpen.exe. Everything it needs is in this folder; nothing to install.
Keep the folder together (StrandOpen.exe looks for runtime\ beside it). You can
move the whole folder anywhere, for example to C:\Programs\StrandOpen, and pin
StrandOpen.exe to the Start menu or taskbar.

Windows may say it protected your PC the first time, because the program is not
signed with a paid certificate: choose "More info", then "Run anyway".

Your data and the reference downloads go to %LOCALAPPDATA%\StrandOpen unless
you choose another folder on the People page. Needs Windows 10 (version 1903 or
later) or Windows 11, 64-bit.
EOF
# unix line ends would show as one line in Notepad of older Windows: write CRLF
for f in "START HERE.txt" THIRD_PARTY.txt; do sed -i 's/$/\r/' "$D/$f"; done

# nothing from the build machine, in any file, binaries included (the result
# captured, not piped into grep -q: see tools/release.sh)
leaks=$(grep -rlaF -- "$HOME" "$D" || true)
if [ -n "$leaks" ]; then
    echo "the release mentions the build machine's home folder:"; echo "$leaks" | head; exit 1
fi
tools/check_public.sh "$D" || { echo "release stopped: private details in $D"; exit 1; }
# Explorer's own unzip fails on paths over 260 characters: keep ours well under,
# leaving room for where the person unzips it
longest=$(cd dist && find "$NAME" | awk '{ print length($0) }' | sort -n | tail -1)
echo "longest path inside the zip: $longest characters"
[ "$longest" -le 170 ] || { echo "paths too long for Windows' unzip"; cd dist && find "$NAME" | awk 'length($0) > 170' | head; exit 1; }
(cd dist && zip -q -r -9 "$NAME.zip" "$NAME")
du -sh "$D" "dist/$NAME.zip"
echo "release ready: dist/$NAME.zip"
