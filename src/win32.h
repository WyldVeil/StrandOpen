/* ========================================================================
   src/win32.h - the Windows calls Strand Open needs, behind plain C types.

   src/win32.c is compiled on its own (only for Windows): <windows.h> and
   raylib.h cannot share a file (both define Rectangle, CloseWindow,
   DrawText and more). Strings are UTF-8 (the program's manifest sets the
   UTF-8 code page) and paths may use / or \.
   ======================================================================== */
#ifndef STRAND_WIN32_H
#define STRAND_WIN32_H
#include <stdbool.h>
#include <stddef.h>

enum { W32_CONFIG, W32_DATA, W32_DOCUMENTS };

bool w32ExePath(char *out, size_t n);                       /* this program, with / separators */
bool w32Folder(int which, char *out, size_t n);             /* %APPDATA%, %LOCALAPPDATA%, Documents */
/* the bundled runtime on PATH and its variables set, for every child (runtime: the folder) */
void w32RuntimeEnv(const char *runtime, const char *cacheDir);
/* a child program with no console window. outFile: its output and errors go
   there and it is left to run on its own (it outlives this program);
   otherwise its output comes back on *fdOut (read() it). */
bool w32Spawn(char *const argv[], const char *outFile, int *fdOut, unsigned *pidOut, void **hOut);
int w32Wait(void *h, bool stop);                            /* exit code, -1 if unknown; stop ends it first */
bool w32Alive(unsigned pid);
double w32FreeBytes(const char *path);
bool w32Replace(const char *from, const char *to);         /* rename over an existing file */
bool w32PickFile(void *owner, char *out, size_t n);
bool w32PickFolder(void *owner, const char *title, char *out, size_t n);
void w32Open(const char *pathOrUrl);                        /* the default program or browser */
bool w32Relaunch(char *const argv[]);                       /* start a new copy; the caller then exits */
long long w32Pread(int fd, void *buf, size_t n, long long off); /* a read at an offset, -1 on failure */
/* the bundled command-line tools (bcftools, samtools, plink2) read paths in the
   system's legacy code page: false when this path has characters it lacks */
bool w32ToolsCanRead(const char *path);

#endif
