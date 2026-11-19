#!/bin/bash
# Builds runtime/: every program the import needs, so nobody installs anything.
#   env/   a conda-forge + bioconda environment: Python 3.12 with numpy,
#          scipy, numba and matplotlib (for hmmix), bcftools, samtools,
#          htslib, plink2; trimmed of tests, headers, docs and Perl
#   jre/   a Java runtime cut down (jlink) to the modules PharmCAT uses
#   bin/   small launchers that find their files relative to themselves, so
#          the folder can live anywhere
#   LICENSES/  each package's licence files and a summary
# conda-forge builds against glibc 2.17, so the result runs on any current
# Linux. The environment is built at a neutral path (/tmp/strandopen-runtime)
# because conda writes its own location into some files.
#   tools/bundle_runtime.sh [OUT]       (default build/runtime)
set -euo pipefail
cd "$(dirname "$0")/.."
ENV=/tmp/strandopen-runtime
OUT=${1:-build/runtime}
CACHE=${MAMBA_ROOT_PREFIX:-$PWD/build/mamba-root}
PKGS="python=3.12 numpy scipy numba matplotlib-base pip bcftools samtools htslib plink2 openjdk=21"
if [ ! -x "$ENV/bin/bcftools" ]; then
    MAMBA_ROOT_PREFIX="$CACHE" micromamba create -y -p "$ENV" -c conda-forge -c bioconda --strict-channel-priority $PKGS
    "$ENV/bin/python" -m pip install --no-deps --no-cache-dir hmmix==0.8.0
fi
rm -rf "$OUT"
mkdir -p "$OUT/bin" "$OUT/LICENSES"

# ---- Java: only the modules PharmCAT uses (jdeps on pharmcat-3.4.0-all.jar,
# plus the ones it loads as services: the zip file system, charsets, locales)
"$ENV/bin/jlink" --add-modules java.base,java.desktop,java.naming,java.scripting,java.security.jgss,java.sql,java.xml.crypto,jdk.dynalink,jdk.unsupported,java.logging,java.management,jdk.crypto.ec,jdk.zipfs,jdk.charsets,jdk.localedata \
    --strip-debug --no-header-files --no-man-pages --compress=zip-6 --output "$OUT/jre"

# ---- the environment, trimmed
rsync -a --delete \
    --exclude 'lib/jvm/' --exclude 'include/' --exclude 'share/man/' --exclude 'share/doc/' --exclude 'share/info/' \
    --exclude 'share/locale/' --exclude 'share/terminfo/' --exclude 'share/gtk-doc/' --exclude 'lib/perl5/' --exclude 'bin/perl*' \
    --exclude 'conda-meta/' --exclude '*.a' --exclude 'lib/pkgconfig/' --exclude 'lib/cmake/' --exclude 'fonts/' \
    --exclude 'lib/python3.12/test/' --exclude 'lib/python3.12/idlelib/' --exclude 'lib/python3.12/tkinter/' \
    --exclude 'lib/python3.12/turtledemo/' --exclude 'lib/python3.12/ensurepip/' --exclude 'lib/python3.12/lib2to3/' \
    --exclude 'site-packages/*/tests/' --exclude 'site-packages/*/*/tests/' --exclude 'site-packages/*/*/*/tests/' \
    --exclude 'site-packages/numba/tests/' --exclude '__pycache__/' \
    "$ENV/" "$OUT/env/"
rm -f "$OUT/env/bin/java" "$OUT/env/bin/javac" "$OUT/env/bin/jlink" "$OUT/env/bin/jdeps"

# ---- launchers
launcher() {   # launcher <name> <command line, with $R for the runtime folder>
    cat > "$OUT/bin/$1" <<EOF
#!/bin/sh
# Strand Open runtime: $1, run from the bundled files wherever this folder is
R=\$(cd "\$(dirname "\$0")/.." && pwd)
export BCFTOOLS_PLUGINS="\$R/env/libexec/bcftools"
export SSL_CERT_FILE="\$R/env/ssl/cacert.pem"
export PYTHONNOUSERSITE=1
export NUMBA_CACHE_DIR="\${XDG_CACHE_HOME:-\$HOME/.cache}/StrandOpen/numba"
exec $2 "\$@"
EOF
    chmod +x "$OUT/bin/$1"
}
launcher python3 '"$R/env/bin/python3.12"'
launcher python '"$R/env/bin/python3.12"'
for t in bcftools samtools bgzip tabix htsfile plink2; do launcher "$t" "\"\$R/env/bin/$t\""; done
launcher java '"$R/jre/bin/java"'
launcher hmmix '"$R/env/bin/python3.12" -c "import sys; from main import main; sys.argv[0] = \"hmmix\"; sys.exit(main())"'

# ---- licences: every package's own files, and a summary
{
    echo "Programs and libraries bundled in runtime/, from conda-forge (https://conda-forge.org)"
    echo "and bioconda (https://bioconda.github.io), each under its own licence; source for each"
    echo "is linked from its package page there. hmmix comes from PyPI (https://pypi.org/project/hmmix/)."
    echo
    printf "%-28s %-22s %s\n" package version licence
} > "$OUT/LICENSES/SUMMARY.txt"
for j in "$ENV"/conda-meta/*.json; do
    "$ENV/bin/python" - "$j" "$CACHE/pkgs" "$OUT/LICENSES" >> "$OUT/LICENSES/SUMMARY.txt" <<'PY'
import json, os, shutil, sys
m = json.load(open(sys.argv[1]))
name, ver = m.get("name", "?"), m.get("version", "?")
print("%-28s %-22s %s" % (name, ver, m.get("license", "unknown")))
# micromamba keeps each extracted package under a per-channel folder
src = os.path.join(m.get("extracted_package_dir") or os.path.join(sys.argv[2], m.get("fn", "").replace(".conda", "").replace(".tar.bz2", "")),
                   "info", "licenses")
if os.path.isdir(src):
    shutil.copytree(src, os.path.join(sys.argv[3], name), dirs_exist_ok=True)
PY
done
printf "%-28s %-22s %s\n" hmmix 0.8.0 MIT >> "$OUT/LICENSES/SUMMARY.txt"
d=$(ls -d "$ENV"/lib/python3.12/site-packages/hmmix-0.8.0.dist-info 2>/dev/null || true)
[ -n "$d" ] && mkdir -p "$OUT/LICENSES/hmmix" && cp "$d"/licenses/* "$d"/LICENSE* "$OUT/LICENSES/hmmix/" 2>/dev/null || true

# ---- nothing from the machine that built it
if grep -rIl "$HOME" "$OUT" 2>/dev/null | head -1 | grep -q .; then
    echo "the runtime mentions the build machine's home folder:"; grep -rIl "$HOME" "$OUT" | head; exit 1
fi
du -sh "$OUT" "$OUT"/*
echo "runtime ready: $OUT"
