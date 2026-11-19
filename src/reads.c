/* ========================================================================
   src/reads.c - aligned reads for the genome browser.

   The sample's reads live in a CRAM (the import pipeline aligns the raw
   FASTQ to GRCh38). A window's reads come from `samtools view` run as a
   child process on a worker thread (the engine's media.c pattern: a
   helper's stdout through a pipe, never a library linked in), parsed from
   SAM into compact reads: start, end, strand, mapping quality, flags, and
   the places a read differs from the reference (mismatched bases,
   insertions, deletions). Reads are packed into rows for the pileup and
   counted per base for the allele-aware coverage track.

   threads: g_reads is guarded by g_readsLock. The worker reads only its
   request snapshot (chrom, start, end, generation, paths) and writes the
   result under the lock; the main thread takes it.
   ======================================================================== */

typedef struct { int pos; uint16_t len; char kind; char base; } ReadEv;   /* kind: X mismatch, I insertion, D deletion, S soft clip */
typedef struct {
    int start, end;      /* 1-based inclusive reference span */
    int ev0, evN;        /* events in the result's ev[] */
    int row;
    uint8_t mapq;
    uint8_t rev, dup, paired, mateOther;
} ReadRec;
typedef struct {
    int chrom, start, end;      /* the window asked for */
    ReadRec *r; int n;
    ReadEv *ev; int evN;
    int rows;
    int *depth;                 /* per base, start..end */
    int (*acgt)[4];             /* per base allele counts */
    char *ref;                  /* the reference over start..end */
    unsigned gen;
    char err[160];
} ReadsResult;

static pthread_mutex_t g_readsLock = PTHREAD_MUTEX_INITIALIZER;
static ReadsResult *g_readsNew, *g_readsCur;
static unsigned g_readsGen;
static _Atomic int g_readsBusy;
static char g_cramPath[1200];
static int g_readsAvail = -1;

static bool readsAvailable(void) {
    if (g_readsAvail < 0) {
        snprintf(g_cramPath, sizeof g_cramPath, "%s/align/reads.cram", g_dataDir);
        char crai[1300];
        snprintf(crai, sizeof crai, "%s.crai", g_cramPath);
        struct stat st;
        char tool[1024];
        g_readsAvail = stat(g_cramPath, &st) == 0 && stat(crai, &st) == 0 && findTool("samtools", tool, sizeof tool);
    }
    return g_readsAvail == 1;
}

typedef struct { int chrom, start, end; unsigned gen; char cram[1200], ref[1100]; } ReadsReq;
static void readsFree(ReadsResult *r) {
    if (!r) return;
    free(r->r); free(r->ev); free(r->depth); free(r->acgt); free(r->ref); free(r);
}
static int baseIdx(char b) { return b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : b == 'T' ? 3 : -1; }

static void *readsWorker(void *arg) {
    ReadsReq q = *(ReadsReq *)arg;
    free(arg);
    ReadsResult *res = xcalloc(1, sizeof *res);
    res->chrom = q.chrom; res->start = q.start; res->end = q.end; res->gen = q.gen;
    int span = q.end - q.start + 1;
    res->depth = xcalloc((size_t)span, sizeof(int));
    res->acgt = xcalloc((size_t)span, sizeof(int[4]));
    /* reference with a margin (reads overhang the window) */
    int m0 = q.start - 600 < 1 ? 1 : q.start - 600, m1 = q.end + 600;
    int refCap = m1 - m0 + 1;
    char *ref = xalloc((size_t)refCap + 1);
    int refN = refSeq(q.chrom, m0, m1, ref, refCap);
    ref[refN] = 0;
    res->ref = xalloc((size_t)span + 1);
    for (int i = 0; i < span; i++) { int k = q.start + i - m0; res->ref[i] = k >= 0 && k < refN ? ref[k] : 'N'; }
    res->ref[span] = 0;

    char region[64];
    snprintf(region, sizeof region, "%s:%d-%d", chromName(q.chrom), q.start, q.end);
    char *argv[] = { "samtools", "view", "-F", "0x904", "-T", q.ref, q.cram, region, NULL };
    Proc p = spawnReader(argv);
    if (!p.pid) { snprintf(res->err, sizeof res->err, "could not start samtools"); goto out; }
    int cap = 4096, evCap = 16384;
    res->r = xalloc(sizeof(ReadRec) * (size_t)cap);
    res->ev = xalloc(sizeof(ReadEv) * (size_t)evCap);
    size_t bufCap = 1 << 20, have = 0;
    char *buf = xalloc(bufCap);
    bool eof = false;
    while (!eof || have) {
        if (!eof) {
            if (have == bufCap) { bufCap *= 2; buf = xrealloc(buf, bufCap); }
            ssize_t got = read(p.fd, buf + have, bufCap - have);
            if (got < 0 && errno == EINTR) continue;
            if (got <= 0) eof = true; else have += (size_t)got;
        }
        char *nl;
        size_t used = 0;
        while ((nl = memchr(buf + used, '\n', have - used))) {
            *nl = 0;
            char *line = buf + used;
            used = (size_t)(nl - buf) + 1;
            /* QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ */
            char *f[11];
            int nf = 0;
            for (char *s = line; nf < 11; ) {
                f[nf++] = s;
                char *t = strchr(s, '\t');
                if (!t) break;
                *t = 0;
                s = t + 1;
            }
            if (nf < 10) continue;
            int flag = atoi(f[1]);
            int pos = atoi(f[3]);
            if (res->n >= 60000) continue;
            if (res->n == cap) { cap *= 2; res->r = xrealloc(res->r, sizeof(ReadRec) * (size_t)cap); }
            ReadRec *rd = &res->r[res->n];
            memset(rd, 0, sizeof *rd);
            rd->mapq = (uint8_t)clampi(atoi(f[4]), 0, 255);
            rd->rev = (flag & 16) != 0;
            rd->dup = (flag & 1024) != 0;
            rd->paired = (flag & 1) != 0;
            rd->mateOther = (flag & 1) && !(flag & 8) && strcmp(f[6], "=") != 0;
            rd->ev0 = res->evN;
            const char *cig = f[5], *seq = f[9];
            int rp = pos, qp = 0;
            while (*cig && *cig != '*') {
                int len = (int)strtol(cig, (char **)&cig, 10);
                char op = *cig++;
                if (res->evN + len + 4 >= evCap) { evCap = evCap * 2 + len; res->ev = xrealloc(res->ev, sizeof(ReadEv) * (size_t)evCap); }
                switch (op) {
                case 'M': case '=': case 'X':
                    for (int i = 0; i < len; i++) {
                        int at = rp + i;
                        char b = seq[qp + i];
                        int k = at - m0;
                        char rb = k >= 0 && k < refN ? ref[k] : 'N';
                        if (at >= q.start && at <= q.end && !rd->dup) {
                            res->depth[at - q.start]++;
                            int bi = baseIdx(b);
                            if (bi >= 0) res->acgt[at - q.start][bi]++;
                        }
                        if (b != rb && rb != 'N' && b != 'N')
                            res->ev[res->evN++] = (ReadEv){ at, 1, 'X', b };
                    }
                    rp += len; qp += len; break;
                case 'I': res->ev[res->evN++] = (ReadEv){ rp, (uint16_t)len, 'I', 0 }; qp += len; break;
                case 'D': case 'N':
                    res->ev[res->evN++] = (ReadEv){ rp, (uint16_t)len, op == 'D' ? 'D' : 'N', 0 };
                    rp += len; break;
                case 'S': res->ev[res->evN++] = (ReadEv){ qp == 0 ? pos : rp, (uint16_t)len, 'S', 0 }; qp += len; break;
                default: break;
                }
            }
            rd->start = pos;
            rd->end = rp - 1;
            rd->evN = res->evN - rd->ev0;
            res->n++;
        }
        memmove(buf, buf + used, have - used);
        have -= used;
        if (eof && have) break;
    }
    free(buf);
    reap(&p, false);
    /* pack into rows: greedy, by start (samtools gives them sorted) */
    {
        int rowsCap = 512;
        int *rowEnd = xcalloc((size_t)rowsCap, sizeof(int));
        int rows = 0;
        for (int i = 0; i < res->n; i++) {
            ReadRec *rd = &res->r[i];
            int k = 0;
            while (k < rows && rowEnd[k] + 3 >= rd->start) k++;
            if (k == rows) {
                if (rows == rowsCap) { rd->row = -1; continue; }
                rows++;
            }
            rowEnd[k] = rd->end;
            rd->row = k;
        }
        res->rows = rows;
        free(rowEnd);
    }
out:
    free(ref);
    pthread_mutex_lock(&g_readsLock);
    if (res->gen == g_readsGen) { readsFree(g_readsNew); g_readsNew = res; }
    else readsFree(res);
    pthread_mutex_unlock(&g_readsLock);
    atomic_fetch_sub(&g_readsBusy, 1);
    return NULL;
}

static void readsRequest(int chrom, long long start, long long end) {
    if (!readsAvailable()) return;
    ReadsReq *q = xcalloc(1, sizeof *q);
    q->chrom = chrom;
    q->start = (int)(start < 1 ? 1 : start);
    q->end = (int)end;
    pthread_mutex_lock(&g_readsLock);
    q->gen = ++g_readsGen;
    pthread_mutex_unlock(&g_readsLock);
    snprintf(q->cram, sizeof q->cram, "%s", g_cramPath);
    snprintf(q->ref, sizeof q->ref, "%s/GRCh38.fa", g_refDir);
    pthread_t t;
    atomic_fetch_add(&g_readsBusy, 1);
    if (thrStart(&t, readsWorker, q, "reads")) pthread_detach(t);
    else { atomic_fetch_sub(&g_readsBusy, 1); free(q); }
}
/* the newest finished window, or NULL; the caller keeps showing the old one */
static ReadsResult *readsCurrent(void) {
    pthread_mutex_lock(&g_readsLock);
    if (g_readsNew) { readsFree(g_readsCur); g_readsCur = g_readsNew; g_readsNew = NULL; }
    pthread_mutex_unlock(&g_readsLock);
    return g_readsCur;
}
static int readsState(void) {
    if (!readsAvailable()) return 0;
    if (atomic_load(&g_readsBusy) > 0) return 1;
    return g_readsCur ? (g_readsCur->err[0] ? 3 : 2) : 0;
}
