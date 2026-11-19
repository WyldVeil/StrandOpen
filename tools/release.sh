#!/bin/bash
# Builds the release: dist/StrandOpen-<version>-linux-x86_64/ and its .tar.gz.
# Everything needed to run is inside (the program, the import, the bundled
# tools), so nobody installs anything: unpack, run ./strandopen.
#   tools/release.sh
set -euo pipefail
cd "$(dirname "$0")/.."
VER=$(grep -o 'ST_VERSION "[^"]*"' src/app.h | cut -d'"' -f2)
NAME="StrandOpen-$VER-linux-x86_64"
D="dist/$NAME"
./build.sh portable
[ -x build/runtime/bin/python3 ] || tools/bundle_runtime.sh
rm -rf "$D" "dist/$NAME.tar.gz"
mkdir -p "$D"
cp build/portable/strandopen "$D/strandopen"
rsync -a --exclude '__pycache__/' --exclude '*.pyc' --exclude '*.log' pipeline/ "$D/pipeline/"
rsync -a data/ "$D/data/"
cp -a build/runtime "$D/runtime"
cp README.md LICENSE "$D/"
[ -d docs ] && rsync -a docs/ "$D/docs/"
mkdir -p "$D/share"
cp assets/icon/strand-256.png "$D/share/strandopen.png"
cat > "$D/THIRD_PARTY.txt" <<'EOF'
Strand Open includes, unchanged or compiled in:
  raylib 6.2 (zlib licence) https://github.com/raysan5/raylib
  SQLite 3.54.1 (public domain) https://sqlite.org
  Inter (SIL Open Font Licence 1.1) https://rsms.me/inter/
  JetBrains Mono (SIL Open Font Licence 1.1) https://www.jetbrains.com/lp/mono/
  Natural Earth 1:50m land outlines (public domain) https://www.naturalearthdata.com/
runtime/ holds the bundled programs (bcftools, samtools, plink2, Python, Java and
more), each under its own licence: see runtime/LICENSES/SUMMARY.txt.
The reference data downloaded at setup comes from its publishers and stays under
their terms (listed in README.md).
EOF
cp assets/fonts/Inter-LICENSE.txt "$D/share/Inter-LICENSE.txt"
cp assets/fonts/JetBrainsMono-OFL.txt "$D/share/JetBrainsMono-OFL.txt"
cat > "$D/add-to-menu.sh" <<'EOF'
#!/bin/sh
# Adds Strand Open to your application menu (for this user only).
D=$(cd "$(dirname "$0")" && pwd)
A="${XDG_DATA_HOME:-$HOME/.local/share}/applications"
mkdir -p "$A"
cat > "$A/strandopen.desktop" <<DESK
[Desktop Entry]
Type=Application
Name=Strand Open
Comment=Explore your DNA from a chip or genome file
Exec="$D/strandopen"
Icon=$D/share/strandopen.png
Categories=Science;Biology;Education;
Terminal=false
DESK
echo "Added: $A/strandopen.desktop"
EOF
chmod +x "$D/add-to-menu.sh" "$D/strandopen"
# nothing from the build machine
if grep -rIl "$HOME" "$D" | grep -q . || strings "$D/strandopen" | grep -q "$HOME"; then
    echo "the release mentions the build machine's home folder:"; grep -rIl "$HOME" "$D" | head; strings "$D/strandopen" | grep "$HOME" | head -3; exit 1
fi
tar -C dist -czf "dist/$NAME.tar.gz" "$NAME"
du -sh "$D" "dist/$NAME.tar.gz"
echo "release ready: dist/$NAME.tar.gz"
