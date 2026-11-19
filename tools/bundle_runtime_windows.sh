#!/bin/bash
# Builds build/win/runtime/: every program the import needs on Windows, made
# on Linux, so nobody installs anything.
#   env/    a conda-forge win-64 environment: Python 3.12 with numpy, scipy,
#           numba and matplotlib (for hmmix), OpenBLAS; trimmed. Same versions
#           as the Linux bundle. sitecustomize.py (tools/windows) lets the
#           pipeline's Unix-style commands run with the tools below.
#   jre/    a Java runtime cut down (jlink, from the Windows jmods) to the
#           modules PharmCAT uses
#   tools/  bcftools, samtools, bgzip, tabix, htsfile and the DLLs they load
#           (MSYS2 UCRT64 builds), bcftools' plugins, plink2 (its own Windows
#           build) and busybox-w32 (sh, awk, sort... for the shell steps)
#   LICENSES/  each component's licence, a summary, and busybox's source
#   tools/bundle_runtime_windows.sh [OUT]       (default build/win/runtime)
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=${1:-build/win/runtime}
WB=build/win
ENV=/tmp/strandopen-winenv                 # neutral path: conda writes its location into a few files
CACHE=$PWD/build/mamba-root-win
LINUX_JLINK=${LINUX_JLINK:-/tmp/strandopen-runtime/bin/jlink}   # same JDK version as the Windows jmods
MSYS=https://repo.msys2.org/mingw/ucrt64
PLINK=https://s3.amazonaws.com/plink2-assets/beta1/plink2_win64_20261107.zip
BUSYBOX=https://frippery.org/files/busybox/busybox64u.exe
BUSYBOX_COMMIT=5f09fb446                   # the source of that build (BusyBox v1.39.0-FRP-6354-g5f09fb446)
mkdir -p "$WB/dl" "$WB/msys2/pkgs"

# ---- the conda-forge environment
if [ ! -x "$ENV/python.exe" ]; then
    MAMBA_ROOT_PREFIX="$CACHE" micromamba create -y -q --platform win-64 -p "$ENV" -c conda-forge --strict-channel-priority \
        "python=3.12.16" "numpy=2.5.5" "scipy=1.18.1" "numba=0.68.1" "matplotlib-base=3.11.2" "openjdk=21.0.13.1" "libblas=*=*openblas" pip
fi
# micromamba, making a Windows environment on Linux, puts noarch Python packages
# in the Linux layout: move them where Windows Python looks
if [ -d "$ENV/lib/python3.12/site-packages" ]; then
    cp -a "$ENV/lib/python3.12/site-packages/." "$ENV/Lib/site-packages/"
    rm -rf "$ENV/lib/python3.12"
fi
if [ ! -f "$ENV/Lib/site-packages/main.py" ]; then
    T=$(mktemp -d)
    "${PYLINUX:-/tmp/strandopen-runtime/bin/python}" -m pip install -q --no-deps --no-compile --target "$T" hmmix==0.8.0
    rm -rf "$T/bin"
    cp -a "$T/." "$ENV/Lib/site-packages/"
    rm -rf "$T"
fi

rm -rf "$OUT"
mkdir -p "$OUT/tools/plugins" "$OUT/LICENSES"

# ---- Java: the modules PharmCAT uses (as tools/bundle_runtime.sh)
"$LINUX_JLINK" --module-path "$ENV/Library/lib/jvm/jmods" \
    --add-modules java.base,java.desktop,java.naming,java.scripting,java.security.jgss,java.sql,java.xml.crypto,jdk.dynalink,jdk.unsupported,java.logging,java.management,jdk.crypto.ec,jdk.zipfs,jdk.charsets,jdk.localedata \
    --strip-debug --no-header-files --no-man-pages --compress=zip-6 --output "$OUT/jre"

# ---- the environment, trimmed
rsync -a --delete \
    --exclude 'Library/lib/jvm/' --exclude 'Library/bin/java*' --exclude 'Library/bin/jlink*' --exclude 'Library/bin/jdeps*' \
    --exclude 'Library/bin/jar.exe' --exclude 'Library/bin/jcmd*' --exclude 'Library/bin/jshell*' --exclude 'Library/bin/keytool*' \
    --exclude 'include/' --exclude 'Library/include/' --exclude 'Library/share/doc/' --exclude 'Library/share/man/' \
    --exclude 'Library/lib/*.lib' --exclude 'Library/lib/cmake/' --exclude 'Library/lib/pkgconfig/' --exclude 'libs/' \
    --exclude 'conda-meta/' --exclude 'Lib/test/' --exclude 'Lib/idlelib/' --exclude 'Lib/tkinter/' --exclude 'Lib/turtledemo/' \
    --exclude 'Lib/ensurepip/' --exclude 'tcl/' --exclude 'Library/lib/tcl*' --exclude 'Library/lib/tk*' \
    --exclude 'site-packages/*/tests/' --exclude 'site-packages/*/*/tests/' --exclude 'site-packages/*/*/*/tests/' \
    --exclude 'site-packages/numba/tests/' --exclude '__pycache__/' --exclude '*.pyc' \
    "$ENV/" "$OUT/env/"
cp tools/windows/sitecustomize.py "$OUT/env/Lib/site-packages/sitecustomize.py"

# ---- bcftools, samtools, htslib from MSYS2, with exactly the DLLs they load
cd "$WB/msys2"
[ -s ucrt64.db ] || curl -sf -o ucrt64.db "$MSYS/ucrt64.db"
python3 - <<'PY'
import tarfile, re, json
db = tarfile.open("ucrt64.db")
pk, prov = {}, {}
for m in db.getmembers():
    if not m.name.endswith("/desc"):
        continue
    f, cur = {}, None
    for line in db.extractfile(m).read().decode().splitlines():
        if line.startswith("%") and line.endswith("%"):
            cur = line.strip("%"); f[cur] = []
        elif line and cur:
            f[cur].append(line)
    pk[f["NAME"][0]] = f
    for p in f.get("PROVIDES", []):
        prov[re.split("[<>=]", p)[0]] = f["NAME"][0]
want, seen = ["mingw-w64-ucrt-x86_64-%s" % n for n in ("bcftools", "samtools", "htslib")], []
while want:
    n = want.pop()
    n = n if n in pk else prov.get(n, n)
    if n in seen or n not in pk:
        continue
    seen.append(n)
    want += [re.split("[<>=]", d)[0] for d in pk[n].get("DEPENDS", [])]
json.dump([(n, pk[n]["FILENAME"][0], pk[n]["VERSION"][0], pk[n].get("LICENSE", ["?"])) for n in seen], open("closure.json", "w"))
PY
python3 -c "import json; [print(f) for n, f, v, l in json.load(open('closure.json'))]" | while read -r f; do
    [ -s "pkgs/$f" ] || curl -sf -o "pkgs/$f" "$MSYS/$f"
done
rm -rf root && mkdir root
for f in pkgs/*.zst; do tar --zstd -xf "$f" -C root; done
cd - >/dev/null
T="$OUT/tools"
cp "$WB/msys2/root/ucrt64/bin/"{bcftools,samtools,bgzip,tabix,htsfile}.exe "$T/"
cp "$WB/msys2/root/ucrt64/libexec/bcftools/"*.dll "$T/plugins/"
python3 - "$WB/msys2/root/ucrt64/bin" "$T" <<'PY'
import os, shutil, subprocess, sys
B, T = sys.argv[1], sys.argv[2]
have = {f.lower(): f for f in os.listdir(B)}
todo = [os.path.join(T, f) for f in os.listdir(T) if f.endswith(".exe")]
todo += [os.path.join(T, "plugins", f) for f in os.listdir(os.path.join(T, "plugins"))]
need = set()
while todo:
    out = subprocess.run(["x86_64-w64-mingw32-objdump", "-p", todo.pop()], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "DLL Name:" in line:
            d = line.split("DLL Name:")[1].strip().lower()
            if d in have and d not in need:
                need.add(d)
                todo.append(os.path.join(B, have[d]))
for d in sorted(need):
    shutil.copy(os.path.join(B, have[d]), T)
print("MSYS2 DLLs bundled:", len(need))
PY

# ---- plink2 and busybox-w32
[ -s "$WB/dl/plink2_win64.zip" ] || curl -sfL -o "$WB/dl/plink2_win64.zip" "$PLINK"
unzip -o -q "$WB/dl/plink2_win64.zip" plink2.exe -d "$T"
[ -s "$WB/dl/busybox64u.exe" ] || curl -sfL -o "$WB/dl/busybox64u.exe" "$BUSYBOX"
cp "$WB/dl/busybox64u.exe" "$T/busybox.exe"
SRC="$WB/dl/busybox-w32-$BUSYBOX_COMMIT.tar.gz"
[ -s "$SRC" ] || curl -sfL -o "$SRC" "https://github.com/rmyorston/busybox-w32/archive/$BUSYBOX_COMMIT.tar.gz"

# ---- licences
L="$OUT/LICENSES"
mkdir -p "$L/conda-forge" "$L/msys2" "$L/source"
{
    echo "Programs and libraries bundled in runtime\\, each under its own licence. None is modified."
    echo
    echo "env\\ and jre\\: conda-forge (https://conda-forge.org) win-64 packages; the source for each is"
    echo "linked from its package page there. hmmix comes from PyPI (https://pypi.org/project/hmmix/)."
    echo "tools\\: MSYS2 UCRT64 packages (https://www.msys2.org); the source package for each version"
    echo "below is at https://repo.msys2.org/mingw/sources/ (named mingw-w64-<name>-<version>.src.tar.zst)."
    echo "tools\\plink2.exe: PLINK 2.0 (GPL-3.0), https://www.cog-genomics.org/plink/2.0/, source at"
    echo "https://github.com/chrchang/plink-ng."
    echo "tools\\busybox.exe: busybox-w32 (GPL-2.0), https://frippery.org/busybox/; its source is included"
    echo "in source\\ (commit $BUSYBOX_COMMIT of https://github.com/rmyorston/busybox-w32)."
    echo
    printf "%-34s %-24s %s\n" package version licence
} > "$L/SUMMARY.txt"
for j in "$ENV"/conda-meta/*.json; do
    python3 - "$j" "$CACHE/pkgs" "$L/conda-forge" >> "$L/SUMMARY.txt" <<'PY'
import json, os, shutil, sys
m = json.load(open(sys.argv[1]))
name, ver = m.get("name", "?"), m.get("version", "?")
print("%-34s %-24s %s" % (name, ver, m.get("license", "unknown")))
src = os.path.join(m.get("extracted_package_dir") or os.path.join(sys.argv[2], m.get("fn", "").replace(".conda", "").replace(".tar.bz2", "")), "info", "licenses")
if os.path.isdir(src):
    shutil.copytree(src, os.path.join(sys.argv[3], name), dirs_exist_ok=True)
PY
done
python3 - "$WB/msys2/closure.json" "$WB/msys2/root/ucrt64/share/licenses" "$L/msys2" >> "$L/SUMMARY.txt" <<'PY'
import json, os, shutil, sys
for n, f, v, lic in json.load(open(sys.argv[1])):
    short = n.replace("mingw-w64-ucrt-x86_64-", "")
    print("%-34s %-24s %s" % ("msys2:" + short, v, ", ".join(lic)))
    src = os.path.join(sys.argv[2], short)
    if os.path.isdir(src):
        shutil.copytree(src, os.path.join(sys.argv[3], short), dirs_exist_ok=True)
PY
printf "%-34s %-24s %s\n" "plink2" "2.0.0-b.1.4 (20261107)" "GPL-3.0-or-later" "busybox-w32" "1.39.0-FRP-6354" "GPL-2.0-only" >> "$L/SUMMARY.txt"
cp /usr/share/licenses/spdx/GPL-3.0-or-later.txt "$L/GPL-3.0.txt"
cp /usr/share/licenses/spdx/GPL-2.0-only.txt "$L/GPL-2.0.txt"
cp "$SRC" "$L/source/"
du -sh "$OUT" "$OUT"/*
