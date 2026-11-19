#!/bin/bash
# Strand Open build. Everything it compiles is in this folder: raylib 6.2
# (vendor/raylib, zlib licence) and SQLite (vendor/sqlite, public domain) are
# compiled once into build/obj; the app itself is one translation unit,
# game.c, which includes every module.
#
#   ./build.sh            -> build/strandopen
#   ./build.sh sanitize   -> build/strandopen-san (AddressSanitizer + UBSan)
#   ./build.sh clean      remove build/
#
# Needs: a C compiler with C23 and #embed (gcc 15 or clang 19), and the X11
# and OpenGL development headers (libx11, mesa).
set -euo pipefail
cd "$(dirname "$0")"
ROOT="$PWD"
RAYLIB="$ROOT/vendor/raylib/src"
B="${BUILD_DIR:-$ROOT/build}"
CC="${CC:-cc}"
mkdir -p "$B/obj"
STD=(-std=gnu23 -ffp-contract=off)
RLMODS="rcore rshapes rtextures rtext rmodels rglfw raudio"
RLFLAGS=(-std=gnu99 -O2 -DPLATFORM_DESKTOP -DGRAPHICS_API_OPENGL_33 -D_GLFW_X11
         -DSUPPORT_FILEFORMAT_JPG=1 -DSUPPORT_SCREEN_CAPTURE=0 -DSUPPORT_GPU_SKINNING=0
         -I"$RAYLIB" -I"$RAYLIB/external/glfw/include")
SQLFLAGS=(-O2 -DSQLITE_THREADSAFE=2 -DSQLITE_OMIT_LOAD_EXTENSION -DSQLITE_DEFAULT_MEMSTATUS=0
          -DSQLITE_ENABLE_FTS5 -DSQLITE_DQS=0 -DSQLITE_DEFAULT_MMAP_SIZE=1073741824 -DSQLITE_MAX_MMAP_SIZE=17179869184)

objs() {   # raylib and sqlite, compiled once and kept until their flags change
    local key="$CC ${RLFLAGS[*]} | ${SQLFLAGS[*]}" stamp="$B/obj/.flags" m pids=()
    if [[ ! -f "$stamp" || "$(cat "$stamp")" != "$key" ]]; then rm -f "$B"/obj/*.o; fi
    for m in $RLMODS; do
        [[ -s "$B/obj/$m.o" && "$B/obj/$m.o" -nt "$RAYLIB/$m.c" ]] && continue
        echo "-- raylib $m"
        "$CC" "${RLFLAGS[@]}" -c "$RAYLIB/$m.c" -o "$B/obj/$m.o" & pids+=($!)
    done
    if [[ ! -s "$B/obj/sqlite3.o" || "$B/obj/sqlite3.o" -ot "$ROOT/vendor/sqlite/sqlite3.c" ]]; then
        echo "-- sqlite"
        "$CC" "${SQLFLAGS[@]}" -c "$ROOT/vendor/sqlite/sqlite3.c" -o "$B/obj/sqlite3.o" & pids+=($!)
    fi
    local p rc=0
    for p in ${pids[@]+"${pids[@]}"}; do wait "$p" || rc=1; done
    [[ $rc -eq 0 ]] || { rm -f "$B/obj/.flags"; echo "object compile failed"; exit 1; }
    echo "$key" > "$stamp"
    OBJS=()
    for m in $RLMODS; do OBJS+=("$B/obj/$m.o"); done
    OBJS+=("$B/obj/sqlite3.o")
}

WARN=(-Wall -Wextra -Wno-unused-function -Wno-unused-parameter -Wno-missing-field-initializers -Wno-sign-compare -Wno-format-truncation)
INC=(-I"$RAYLIB" -I"$ROOT/vendor/sqlite" -I"$ROOT")
LIBS=(-lX11 -lm -lpthread -ldl -lrt)
case "${1:-native}" in
    native)
        objs
        "$CC" "${STD[@]}" -O2 -g "${WARN[@]}" "${INC[@]}" -o "$B/strandopen" "$ROOT/game.c" "${OBJS[@]}" "${LIBS[@]}"
        ls -la "$B/strandopen" ;;
    portable)
        # the release build: conda-forge GCC against a glibc 2.17 sysroot, so the
        # binary runs on any current Linux (made by tools/release.sh; BUILDENV
        # is the compiler environment it creates)
        BE="${BUILDENV:-/tmp/strandopen-buildenv}"
        CC="$BE/bin/x86_64-conda-linux-gnu-gcc"
        # source paths (__FILE__ in vendored asserts) relative, so the binary names no folder of the build machine
        PFX=(-ffile-prefix-map="$ROOT"=.)
        RLFLAGS+=(-I"$BE/include" "${PFX[@]}")
        SQLFLAGS+=("${PFX[@]}")
        B="$B/portable"
        mkdir -p "$B/obj"
        objs
        "$CC" "${STD[@]}" -O2 "${PFX[@]}" "${WARN[@]}" "${INC[@]}" -I"$BE/include" -o "$B/strandopen" "$ROOT/game.c" "${OBJS[@]}" \
            -L"$BE/lib" -static-libgcc "${LIBS[@]}"
        echo "newest glibc symbol needed: $(objdump -T "$B/strandopen" | grep -o 'GLIBC_[0-9.]*' | sort -V | tail -1)"
        ls -la "$B/strandopen" ;;
    sanitize)
        objs
        "$CC" "${STD[@]}" -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer "${WARN[@]}" "${INC[@]}" \
            -o "$B/strandopen-san" "$ROOT/game.c" "${OBJS[@]}" "${LIBS[@]}"
        ls -la "$B/strandopen-san" ;;
    clean) rm -rf "$B" ;;
    *) echo "usage: $0 [native|sanitize|clean]"; exit 1 ;;
esac
