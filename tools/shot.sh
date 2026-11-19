#!/bin/bash
# Headless screenshots of Strand Open, never on the real desktop:
#   tools/shot.sh <out.png> <page[:arg]> [frames] [w h]
# A private Xvfb, software GL, silent audio stubs and a throwaway HOME for
# every run (the engine's own test rules: CONVENTIONS.md, Build + test).
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="$(realpath -m "$1")"; PAGE="${2:-overview}"; FRAMES="${3:-40}"; W="${4:-1600}"; H="${5:-1000}"
T=$(mktemp -d -p "${SHOT_TMP:-${TMPDIR:-/tmp}}" shot.XXXXXX)
mkdir -p "$T/home" "$T/stub"
for n in libasound.so libasound.so.2 libpulse.so libpulse.so.0 libjack.so libjack.so.0; do : > "$T/stub/$n"; done
exec 3>"$T/display"
Xvfb -displayfd 3 -noreset -nolisten tcp -screen 0 "${W}x${H}x24" >/dev/null 2>&1 &
XPID=$!
for _ in $(seq 50); do [[ -s "$T/display" ]] && break; sleep 0.1; done
DISP=":$(cat "$T/display")"
set +e
env -u WAYLAND_DISPLAY -u XDG_SESSION_TYPE DISPLAY="$DISP" HOME="$T/home" XDG_CONFIG_HOME="$T/home/.config" \
    LIBGL_ALWAYS_SOFTWARE=1 STRAND_NO_OPEN=1 LD_LIBRARY_PATH="$T/stub" PULSE_SERVER=unix:/nonexistent \
    STRANDOPEN_ROOT="${STRANDOPEN_ROOT:?set STRANDOPEN_ROOT to a data root}" \
    timeout 300 ${BIN:-./build/strandopen} --size "$W" "$H" --enter "$PAGE" --shot "$FRAMES" "$OUT" ${PROFILE:+--profile "$PROFILE"} > "$T/out.txt" 2> "$T/err.txt"
RC=$?
set -e
kill "$XPID" 2>/dev/null || true
wait "$XPID" 2>/dev/null || true
if [[ ! -s "$OUT" ]]; then echo "no shot (rc $RC)"; tail -20 "$T/err.txt"; tail -5 "$T/out.txt"; fi
grep -i -E 'error|fail|assert|abort' "$T/err.txt" | head -5 || true
[[ -n "${SHOT_KEEP:-}" ]] && cp -r "$T/home" "$SHOT_KEEP"
rm -rf "$T"
[[ -s "$OUT" ]] && echo "shot: $OUT"
