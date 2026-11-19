/* ========================================================================
   src/db.c - the genome database (SQLite, read only).

   The main thread has its own connection for small queries (a gene's row,
   a variant's detail: a millisecond or two). Heavy ones (the variant table
   under a filter and a sort, the browser's window at a wide zoom) go to two
   worker threads, each with a connection of its own, through SLOTS: one
   request per slot, a newer request replaces the older one, and a query
   still running for an older request is interrupted by SQLite's progress
   handler as soon as the slot moves on. Results are flat Rows blocks (one
   allocation) that the main thread takes and frees.

   threads: g_dbSlots and g_dbMeta* are guarded by g_dbLock (main thread
   writes requests and takes results; workers take requests and write
   results). g_dbSlotGen[] is atomic so a worker's progress handler can see
   a request superseded without the lock.
   ======================================================================== */

#include <pthread.h>

static sqlite3 *g_db;
static char g_dbPath[1024];

typedef struct { char *key, *val; } DbMeta;
static DbMeta *g_dbMeta;
static int g_dbMetaN;

/* ---- rows ------------------------------------------------------------ */
typedef struct { char **cells; size_t n, cap; size_t bytes; } RowBuild;
static void rbAdd(RowBuild *b, const char *s) {
    if (b->n == b->cap) {
        b->cap = b->cap ? b->cap * 2 : 256;
        b->cells = xrealloc(b->cells, b->cap * sizeof *b->cells);
    }
    b->cells[b->n++] = s ? xstrdup(s) : NULL;
    if (s) b->bytes += strlen(s) + 1;
}
/* run one statement into a Rows block; NULL params allowed */
static Rows *dbRun(sqlite3 *db, const char *sql, int n, const char *const *params) {
    double t0 = (GetTime)();
    sqlite3_stmt *st = NULL;
    Rows *out = NULL;
    RowBuild b = { 0 };
    char err[200] = "";
    int ncol = 0, nrow = 0;
    char *names[64] = { 0 };
    if (sqlite3_prepare_v2(db, sql, -1, &st, NULL) != SQLITE_OK) {
        snprintf(err, sizeof err, "%s", sqlite3_errmsg(db));
        goto done;
    }
    for (int i = 0; i < n; i++) sqlite3_bind_text(st, i + 1, params[i], -1, SQLITE_TRANSIENT);
    ncol = sqlite3_column_count(st);
    if (ncol > 64) ncol = 64;
    for (int c = 0; c < ncol; c++) names[c] = xstrdup(sqlite3_column_name(st, c) ? sqlite3_column_name(st, c) : "");
    int rcode;
    while ((rcode = sqlite3_step(st)) == SQLITE_ROW) {
        for (int c = 0; c < ncol; c++) rbAdd(&b, (const char *)sqlite3_column_text(st, c));
        nrow++;
    }
    if (rcode != SQLITE_DONE) snprintf(err, sizeof err, "%s", sqlite3_errmsg(db));
done:
    if (st) sqlite3_finalize(st);
    {
        size_t namesBytes = 0;
        for (int c = 0; c < ncol; c++) namesBytes += strlen(names[c]) + 1;
        size_t ptrs = (size_t)ncol * (size_t)nrow + (size_t)ncol;
        size_t total = sizeof(Rows) + ptrs * sizeof(char *) + b.bytes + namesBytes;
        out = xalloc(total);
        memset(out, 0, sizeof *out);
        out->ncol = ncol;
        out->nrow = nrow;
        out->cell = (char **)(out + 1);
        out->colName = out->cell + (size_t)ncol * (size_t)nrow;
        char *w = (char *)(out->colName + ncol);
        for (size_t i = 0; i < b.n; i++) {
            if (b.cells[i]) {
                size_t l = strlen(b.cells[i]) + 1;
                memcpy(w, b.cells[i], l);
                out->cell[i] = w;
                w += l;
                free(b.cells[i]);
            } else out->cell[i] = NULL;
        }
        for (int c = 0; c < ncol; c++) {
            size_t l = strlen(names[c]) + 1;
            memcpy(w, names[c], l);
            out->colName[c] = w;
            w += l;
            free(names[c]);
        }
        free(b.cells);
        snprintf(out->err, sizeof out->err, "%s", err);
        out->ms = ((GetTime)() - t0) * 1000.0;
    }
    return out;
}
static const char *rc(const Rows *r, int row, int col) {
    if (!r || row < 0 || row >= r->nrow || col < 0 || col >= r->ncol) return "";
    const char *s = r->cell[(size_t)row * (size_t)r->ncol + (size_t)col];
    return s ? s : "";
}
static bool rnull(const Rows *r, int row, int col) {
    if (!r || row < 0 || row >= r->nrow || col < 0 || col >= r->ncol) return true;
    return r->cell[(size_t)row * (size_t)r->ncol + (size_t)col] == NULL;
}
static double rf(const Rows *r, int row, int col) { return rnull(r, row, col) ? NAN : atof(rc(r, row, col)); }
static long long ri(const Rows *r, int row, int col) { return rnull(r, row, col) ? 0 : atoll(rc(r, row, col)); }
static void rowsFree(Rows **r) { if (*r) { free(*r); *r = NULL; } }

/* ---- the main thread's connection ------------------------------------ */
static sqlite3 *dbOpenConn(void) {
    sqlite3 *db = NULL;
    if (sqlite3_open_v2(g_dbPath, &db, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, NULL) != SQLITE_OK) {
        if (db) sqlite3_close(db);
        return NULL;
    }
    sqlite3_exec(db, "PRAGMA mmap_size=8000000000; PRAGMA cache_size=-200000; PRAGMA temp_store=MEMORY; PRAGMA query_only=1;", NULL, NULL, NULL);
    return db;
}
static Rows *dbqv(const char *sql, int n, const char *const *params) {
    if (!g_db) { Rows *r = xcalloc(1, sizeof *r); snprintf(r->err, sizeof r->err, "no database"); return r; }
    Rows *r = dbRun(g_db, sql, n, params);
    if (r->err[0]) logBoth("query failed: %s", r->err);
    return r;
}
static Rows *dbq(const char *sql, ...) {
    const char *p[16];
    int n = 0;
    va_list ap;
    va_start(ap, sql);
    const char *s;
    while (n < 16 && (s = va_arg(ap, const char *)) != NULL) p[n++] = s;
    va_end(ap);
    return dbqv(sql, n, p);
}
static const char *dbMeta(const char *key) {
    for (int i = 0; i < g_dbMetaN; i++) if (!strcmp(g_dbMeta[i].key, key)) return g_dbMeta[i].val;
    return "";
}
/* the schema does not change while the app runs: answers are kept */
static bool dbHasTable(const char *name) {
    static struct { char name[32]; bool ok; } memo[64];
    static int memoN;
    for (int i = 0; i < memoN; i++) if (!strcmp(memo[i].name, name)) return memo[i].ok;
    Rows *r = dbq("SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?", name, NULL);
    bool ok = r->nrow > 0;
    rowsFree(&r);
    if (memoN < 64 && strlen(name) < 32) { snprintf(memo[memoN].name, sizeof memo[memoN].name, "%s", name); memo[memoN++].ok = ok; }
    return ok;
}

/* a column of a table (for databases built before the column existed) */
static bool dbColumnExists(const char *table, const char *col) {
    if (!dbHasTable(table)) return false;
    Rows *r = dbq("SELECT 1 FROM pragma_table_info(?) WHERE name=?", table, col, NULL);
    bool ok = r->nrow > 0;
    rowsFree(&r);
    return ok;
}

/* ---- slots and workers ------------------------------------------------ */
typedef struct {
    char *sql;
    char *params[12];
    int n;
    unsigned gen;          /* bumped by each request */
    bool pending, running;
    Rows *result;
} DbSlot;
static DbSlot g_dbSlots[SLOT__N];
static _Atomic unsigned g_dbSlotGen[SLOT__N];
static pthread_mutex_t g_dbLock = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t g_dbCond = PTHREAD_COND_INITIALIZER;
static atomic_bool g_dbQuit;
static pthread_t g_dbThreads[2];
static int g_dbThreadN;

typedef struct { int slot; unsigned gen; } DbWorkCtx;
static int dbProgress(void *u) {
    DbWorkCtx *c = u;
    return atomic_load(&g_dbSlotGen[c->slot]) != c->gen || atomic_load(&g_dbQuit);
}
static void *dbWorker(void *arg) {
    (void)arg;
    sqlite3 *db = dbOpenConn();
    if (!db) return NULL;
    DbWorkCtx ctx = { 0 };
    sqlite3_progress_handler(db, 2000, dbProgress, &ctx);
    for (;;) {
        pthread_mutex_lock(&g_dbLock);
        int pick = -1;
        while (!atomic_load(&g_dbQuit)) {
            for (int i = 0; i < SLOT__N; i++) if (g_dbSlots[i].pending && !g_dbSlots[i].running) { pick = i; break; }
            if (pick >= 0) break;
            pthread_cond_wait(&g_dbCond, &g_dbLock);
        }
        if (atomic_load(&g_dbQuit)) { pthread_mutex_unlock(&g_dbLock); break; }
        DbSlot *s = &g_dbSlots[pick];
        char *sql = s->sql;
        char *params[12];
        int n = s->n;
        memcpy(params, s->params, sizeof params);
        s->sql = NULL;
        memset(s->params, 0, sizeof s->params);
        s->pending = false;
        s->running = true;
        unsigned gen = s->gen;
        pthread_mutex_unlock(&g_dbLock);

        ctx.slot = pick;
        ctx.gen = gen;
        Rows *r = dbRun(db, sql, n, (const char *const *)params);
        free(sql);
        for (int i = 0; i < n; i++) free(params[i]);

        pthread_mutex_lock(&g_dbLock);
        s->running = false;
        if (s->gen == gen && !r->err[0]) { rowsFree(&s->result); s->result = r; }
        else if (s->gen == gen) { rowsFree(&s->result); s->result = r; }
        else free(r);
        pthread_cond_broadcast(&g_dbCond);
        pthread_mutex_unlock(&g_dbLock);
    }
    sqlite3_close(db);
    return NULL;
}
static void dbAsync(StSlot slot, const char *sql, int n, const char *const *params) {
    pthread_mutex_lock(&g_dbLock);
    DbSlot *s = &g_dbSlots[slot];
    free(s->sql);
    for (int i = 0; i < s->n; i++) free(s->params[i]);
    s->sql = xstrdup(sql);
    s->n = n > 12 ? 12 : n;
    for (int i = 0; i < s->n; i++) s->params[i] = xstrdup(params[i] ? params[i] : "");
    s->gen++;
    atomic_store(&g_dbSlotGen[slot], s->gen);
    s->pending = true;
    rowsFree(&s->result);
    pthread_cond_broadcast(&g_dbCond);
    pthread_mutex_unlock(&g_dbLock);
    if (g_dbThreadN == 0) {   /* no workers: run it here */
        pthread_mutex_lock(&g_dbLock);
        char *q = s->sql; s->sql = NULL; s->pending = false;
        char *ps[12]; int pn = s->n; memcpy(ps, s->params, sizeof ps); memset(s->params, 0, sizeof s->params);
        pthread_mutex_unlock(&g_dbLock);
        Rows *r = dbRun(g_db, q, pn, (const char *const *)ps);
        free(q);
        for (int i = 0; i < pn; i++) free(ps[i]);
        pthread_mutex_lock(&g_dbLock);
        s->result = r;
        pthread_mutex_unlock(&g_dbLock);
    }
}
static Rows *dbTake(StSlot slot) {
    pthread_mutex_lock(&g_dbLock);
    Rows *r = NULL;
    if (!g_dbSlots[slot].pending && !g_dbSlots[slot].running) { r = g_dbSlots[slot].result; g_dbSlots[slot].result = NULL; }
    pthread_mutex_unlock(&g_dbLock);
    if (r && r->err[0] && strcmp(r->err, "interrupted")) logBoth("async query failed: %s", r->err);
    return r;
}
static bool dbBusy(StSlot slot) {
    pthread_mutex_lock(&g_dbLock);
    bool b = g_dbSlots[slot].pending || g_dbSlots[slot].running;
    pthread_mutex_unlock(&g_dbLock);
    return b;
}
static void dbPoll(void) {}

static bool dbOpen(const char *path) {
    snprintf(g_dbPath, sizeof g_dbPath, "%s", path);
    struct stat st;
    if (stat(path, &st) != 0) return false;
    g_db = dbOpenConn();
    if (!g_db) return false;
    Rows *m = dbq("SELECT key, value FROM meta", NULL);
    g_dbMetaN = m->nrow;
    g_dbMeta = xcalloc((size_t)(m->nrow + 1), sizeof *g_dbMeta);
    for (int i = 0; i < m->nrow; i++) { g_dbMeta[i].key = xstrdup(rc(m, i, 0)); g_dbMeta[i].val = xstrdup(rc(m, i, 1)); }
    rowsFree(&m);
    for (int i = 0; i < 2; i++)
        if (thrStart(&g_dbThreads[g_dbThreadN], dbWorker, NULL, "db worker")) g_dbThreadN++;
    return true;
}
static void dbShutdown(void) {
    atomic_store(&g_dbQuit, true);
    pthread_mutex_lock(&g_dbLock);
    pthread_cond_broadcast(&g_dbCond);
    pthread_mutex_unlock(&g_dbLock);
    for (int i = 0; i < g_dbThreadN; i++) pthread_join(g_dbThreads[i], NULL);
    g_dbThreadN = 0;
    if (g_db) sqlite3_close(g_db);
    g_db = NULL;
}
