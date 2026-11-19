/* ========================================================================
   src/platform.c - what Strand Open needs from the system and the window.

   The application loop (window, frames, the command-line test harness),
   a small settings file, child processes with a pipe, threads, atomic file
   writes, logging and checked allocation. Everything sits on raylib and
   POSIX; nothing else is needed at run time.

   Command line:
     --size W H            initial window size
     --enter PAGE[:ARG]    open a page at start (see stEnter in app.c)
     --shot N FILE         after N frames save a PNG of the window and quit
     --profile SLUG        open this person (see profile.c)
   ======================================================================== */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"

#define NELEMS(a) (sizeof(a) / sizeof((a)[0]))
static inline float clampf(float v, float lo, float hi) { return v < lo ? lo : v > hi ? hi : v; }
static inline int clampi(int v, int lo, int hi) { return v < lo ? lo : v > hi ? hi : v; }
static inline float clamp01(float v) { return clampf(v, 0.0f, 1.0f); }

/* ---- allocation: out of memory ends the program with a message -------- */
static void *xalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) { fputs("Strand Open: out of memory\n", stderr); abort(); }
    return p;
}
static void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n ? n : 1, sz ? sz : 1);
    if (!p) { fputs("Strand Open: out of memory\n", stderr); abort(); }
    return p;
}
static void *xrealloc(void *old, size_t n) {
    void *p = realloc(old, n ? n : 1);
    if (!p) { fputs("Strand Open: out of memory\n", stderr); abort(); }
    return p;
}
static char *xstrdup(const char *s) {
    size_t n = strlen(s ? s : "") + 1;
    char *p = xalloc(n);
    memcpy(p, s ? s : "", n);
    return p;
}

/* ---- folders ------------------------------------------------------------ */
static char g_appDir[PATH_MAX];      /* where the binary lives (bundled tools: runtime/bin beside build/) */
static char g_cfgDir[PATH_MAX];      /* $XDG_CONFIG_HOME/StrandOpen */
static char g_argProfile[64];        /* --profile SLUG */
static bool g_quitNow;               /* the close button or Ctrl+Q: leave the frame loop, shut down as usual */

static bool mkdirs(const char *path, mode_t mode) {
    char tmp[PATH_MAX];
    snprintf(tmp, sizeof tmp, "%s", path);
    for (char *p = tmp + 1; *p; p++) {
        if (*p != '/') continue;
        *p = 0;
        if (mkdir(tmp, mode) != 0 && errno != EEXIST) return false;
        *p = '/';
    }
    return mkdir(tmp, mode) == 0 || errno == EEXIST;
}
static void platformDirs(const char *argv0) {
    char exe[PATH_MAX];
    ssize_t n = readlink("/proc/self/exe", exe, sizeof exe - 1);
    if (n <= 0) { snprintf(exe, sizeof exe, "%s", argv0); n = (ssize_t)strlen(exe); }
    exe[n] = 0;
    char *slash = strrchr(exe, '/');
    if (slash) *slash = 0;
    snprintf(g_appDir, sizeof g_appDir, "%s", exe);
    const char *x = getenv("XDG_CONFIG_HOME"), *h = getenv("HOME");
    if (x && *x) snprintf(g_cfgDir, sizeof g_cfgDir, "%s/StrandOpen", x);
    else snprintf(g_cfgDir, sizeof g_cfgDir, "%s/.config/StrandOpen", h && *h ? h : "/tmp");
    mkdirs(g_cfgDir, 0700);
}

/* ---- logging: standard error and a log file beside the settings -------- */
static FILE *g_logFile;
__attribute__((format(printf, 1, 2))) static void logBoth(const char *fmt, ...) {
    char line[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(line, sizeof line, fmt, ap);
    va_end(ap);
    fprintf(stderr, "%s\n", line);
    if (!g_logFile && *g_cfgDir) {
        char p[PATH_MAX + 16];
        snprintf(p, sizeof p, "%s/log.txt", g_cfgDir);
        g_logFile = fopen(p, "w");
    }
    if (g_logFile) { fprintf(g_logFile, "%s\n", line); fflush(g_logFile); }
}

/* ---- atomic file writes: write to tmp, flush, sync, rename over final -- */
[[nodiscard]] static bool commitFile(FILE *f, const char *tmp, const char *final) {
    bool ok = fflush(f) == 0 && fsync(fileno(f)) == 0;
    ok = (fclose(f) == 0) && ok;
    if (ok) ok = rename(tmp, final) == 0;
    if (!ok) unlink(tmp);
    return ok;
}

/* ---- settings: key=value lines, read at start, written on change ------ */
typedef struct { char k[64]; char v[512]; } CfgItem;
static CfgItem g_cfg[256];
static int g_cfgN;
static char g_cfgPath[PATH_MAX + 32];

static void cfgLoad(void) {
    snprintf(g_cfgPath, sizeof g_cfgPath, "%s/settings.kv", g_cfgDir);
    FILE *f = fopen(g_cfgPath, "r");
    if (!f) return;
    char line[700];
    while (g_cfgN < (int)NELEMS(g_cfg) && fgets(line, sizeof line, f)) {
        char *eq = strchr(line, '=');
        if (!eq || line[0] == '#') continue;
        *eq = 0;
        line[strcspn(line, "\r\n")] = 0;
        char *v = eq + 1;
        v[strcspn(v, "\r\n")] = 0;
        snprintf(g_cfg[g_cfgN].k, sizeof g_cfg[0].k, "%s", line);
        snprintf(g_cfg[g_cfgN].v, sizeof g_cfg[0].v, "%s", v);
        g_cfgN++;
    }
    fclose(f);
}
static const char *cfgGet(const char *k, const char *def) {
    for (int i = 0; i < g_cfgN; i++) if (!strcmp(g_cfg[i].k, k)) return g_cfg[i].v;
    return def;
}
static void cfgSet(const char *k, const char *v) {
    for (int i = 0; i < g_cfgN; i++) if (!strcmp(g_cfg[i].k, k)) { snprintf(g_cfg[i].v, sizeof g_cfg[i].v, "%s", v); return; }
    if (g_cfgN >= (int)NELEMS(g_cfg)) return;
    snprintf(g_cfg[g_cfgN].k, sizeof g_cfg[0].k, "%s", k);
    snprintf(g_cfg[g_cfgN].v, sizeof g_cfg[0].v, "%s", v);
    g_cfgN++;
}
static bool cfgSave(void) {
    char tmp[PATH_MAX + 48];
    snprintf(tmp, sizeof tmp, "%s.tmp%d", g_cfgPath, (int)getpid());
    FILE *f = fopen(tmp, "w");
    if (!f) return false;
    for (int i = 0; i < g_cfgN; i++) fprintf(f, "%s=%s\n", g_cfg[i].k, g_cfg[i].v);
    return commitFile(f, tmp, g_cfgPath);
}

/* ---- child processes ---------------------------------------------------
   spawnReader starts argv[0] (a full path or a name found by findTool) with
   stdin and stderr on /dev/null and its stdout on a pipe; reap closes the
   pipe, optionally stops the child, and waits for it. */
typedef struct { pid_t pid; int fd; } Proc;

static bool findTool(const char *name, char *out, size_t n) {
    char p[PATH_MAX + 64];
    /* bundled first: <app>/../runtime/bin and <app>/runtime/bin */
    const char *pre[2] = { "%s/../runtime/bin/%s", "%s/runtime/bin/%s" };
    for (int i = 0; i < 2; i++) {
        snprintf(p, sizeof p, pre[i], g_appDir, name);
        if (access(p, X_OK) == 0) { snprintf(out, n, "%s", p); return true; }
    }
    const char *path = getenv("PATH");
    if (!path) return false;
    char buf[8192];
    snprintf(buf, sizeof buf, "%s", path);
    for (char *save = NULL, *d = strtok_r(buf, ":", &save); d; d = strtok_r(NULL, ":", &save)) {
        snprintf(p, sizeof p, "%s/%s", d, name);
        if (access(p, X_OK) == 0) { snprintf(out, n, "%s", p); return true; }
    }
    return false;
}

static Proc spawnReader(char *const argv[]) {
    Proc pr = { 0, -1 };
    int fds[2];
    if (pipe(fds) != 0) return pr;
    char exe[PATH_MAX];
    if (strchr(argv[0], '/')) snprintf(exe, sizeof exe, "%s", argv[0]);
    else if (!findTool(argv[0], exe, sizeof exe)) { close(fds[0]); close(fds[1]); return pr; }
    pid_t pid = fork();
    if (pid < 0) { close(fds[0]); close(fds[1]); return pr; }
    if (pid == 0) {
        int nul = open("/dev/null", O_RDWR);
        if (nul >= 0) { dup2(nul, 0); dup2(nul, 2); }
        dup2(fds[1], 1);
        close(fds[0]);
        close(fds[1]);
        setpgid(0, 0);
        execv(exe, argv);
        _exit(127);
    }
    close(fds[1]);
    fcntl(fds[0], F_SETFD, FD_CLOEXEC);
    pr.pid = pid;
    pr.fd = fds[0];
    return pr;
}
static int reap(Proc *p, bool stop) {
    if (p->fd >= 0) { close(p->fd); p->fd = -1; }
    if (p->pid <= 0) return -1;
    if (stop) kill(-p->pid, SIGTERM);
    int st = 0;
    while (waitpid(p->pid, &st, 0) < 0 && errno == EINTR) {}
    p->pid = 0;
    return WIFEXITED(st) ? WEXITSTATUS(st) : -1;
}

/* ---- threads ------------------------------------------------------------ */
static bool thrStart(pthread_t *t, void *(*fn)(void *), void *arg, const char *what) {
    int e = pthread_create(t, NULL, fn, arg);
    if (e) logBoth("could not start the %s thread (%s)", what, strerror(e));
    return e == 0;
}

/* ---- the application loop ----------------------------------------------- */
typedef struct {
    const char *title;
    const unsigned char *iconPng;
    int iconPngLen;
    int w, h;
    void (*init)(void);
    void (*frame)(float dt, double now);
    void (*shutdown)(void);
    void (*enter)(const char *what);
} AppDef;

static int appMain(const AppDef *d, int argc, char **argv) {
    platformDirs(argc > 0 ? argv[0] : "strandopen");
    cfgLoad();
    int w = d->w, h = d->h, shotFrames = 0;
    const char *enter = NULL, *shotPath = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--size") && i + 2 < argc) { w = atoi(argv[i + 1]); h = atoi(argv[i + 2]); i += 2; }
        else if (!strcmp(argv[i], "--enter") && i + 1 < argc) enter = argv[++i];
        else if (!strcmp(argv[i], "--shot") && i + 2 < argc) { shotFrames = atoi(argv[i + 1]); shotPath = argv[i + 2]; i += 2; }
        else if (!strcmp(argv[i], "--profile") && i + 1 < argc) snprintf(g_argProfile, sizeof g_argProfile, "%s", argv[++i]);
    }
    SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT | FLAG_MSAA_4X_HINT);
    InitWindow(w > 0 ? w : 1600, h > 0 ? h : 1000, d->title);
    if (d->iconPng && d->iconPngLen > 0) {
        Image icon = LoadImageFromMemory(".png", d->iconPng, d->iconPngLen);
        if (icon.data) { SetWindowIcon(icon); UnloadImage(icon); }
    }
    d->init();
    if (enter && d->enter) d->enter(enter);
    int frames = 0;
    while (!WindowShouldClose() && !g_quitNow) {
        BeginDrawing();
        d->frame(GetFrameTime(), GetTime());
        bool shoot = shotPath && ++frames >= shotFrames;
        if (shoot) {
            rlDrawRenderBatchActive();          /* everything drawn so far, before the buffers swap */
            Image img = LoadImageFromScreen();
            ExportImage(img, shotPath);
            UnloadImage(img);
        }
        EndDrawing();
        if (shoot) break;
    }
    d->shutdown();
    CloseWindow();
    if (g_logFile) fclose(g_logFile);
    return 0;
}
