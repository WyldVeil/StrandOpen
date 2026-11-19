/* ========================================================================
   src/genome.c - chromosomes, cytobands, the reference sequence, and the
   formatting every page shares.

   The reference is GRCh38 (Ensembl primary assembly, the same sequence the
   sample's calls were made against), read straight from the FASTA through
   its .fai index with pread: a window of a few kilobases costs a page or
   two of disk.
   ======================================================================== */

#include <fcntl.h>

static char g_dataDir[1024];     /* the open person's folder */
static char g_root[1024], g_refDir[1024];   /* data root; shared reference data */
typedef struct { long long len, off; int lineBases, lineBytes; bool ok; } StFai;
static StFai g_fai[ST_NCHROM + 1];
static int g_faFd = -1;

static const char *const g_chromNames[ST_NCHROM + 1] = { "", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT" };

static int chromIndex(const char *name) {
    if (!name) return 0;
    while (*name == ' ') name++;
    if (!strncasecmp(name, "chr", 3)) name += 3;
    char b[4] = { 0 };
    int n = 0;
    while (name[n] && n < 3 && (isalnum((unsigned char)name[n]))) { b[n] = (char)toupper((unsigned char)name[n]); n++; }
    if (!strcmp(b, "M")) return 25;
    for (int i = 1; i <= ST_NCHROM; i++) if (!strcmp(b, g_chromNames[i])) return i;
    return 0;
}
static const char *chromName(int idx) { return idx >= 1 && idx <= ST_NCHROM ? g_chromNames[idx] : "?"; }

static void genomeInit(void) {
    for (int i = 1; i <= ST_NCHROM; i++) snprintf(g_chrom[i].name, sizeof g_chrom[i].name, "%s", g_chromNames[i]);
    Rows *r = dbq("SELECT idx, length FROM chrom", NULL);
    for (int i = 0; i < r->nrow; i++) {
        int c = (int)ri(r, i, 0);
        if (c >= 1 && c <= ST_NCHROM) g_chrom[c].len = ri(r, i, 1);
    }
    rowsFree(&r);
    r = dbq("SELECT chrom, start, end, name, stain FROM cytoband ORDER BY chrom, start", NULL);
    g_bands = xcalloc((size_t)(r->nrow + 1), sizeof *g_bands);
    for (int i = 0; i < r->nrow; i++) {
        StBand *b = &g_bands[g_bandN++];
        b->chrom = (int)ri(r, i, 0);
        b->start = ri(r, i, 1);
        b->end = ri(r, i, 2);
        snprintf(b->name, sizeof b->name, "%s", rc(r, i, 3));
        snprintf(b->stain, sizeof b->stain, "%s", rc(r, i, 4));
    }
    rowsFree(&r);
    /* the reference FASTA and its index */
    char fa[1100], fai[1110];
    snprintf(fa, sizeof fa, "%s/GRCh38.fa", g_refDir);
    snprintf(fai, sizeof fai, "%s.fai", fa);
    FILE *f = fopen(fai, "r");
    if (f) {
        char line[512];
        while (fgets(line, sizeof line, f)) {
            char name[64];
            long long len, off;
            int lb, lw;
            if (sscanf(line, "%63s %lld %lld %d %d", name, &len, &off, &lb, &lw) != 5) continue;
            int c = chromIndex(name);
            if (c && !strcmp(name, g_chromNames[c])) g_fai[c] = (StFai){ len, off, lb, lw, true };
        }
        fclose(f);
        g_faFd = open(fa, O_RDONLY | O_CLOEXEC);
    }
}
static bool refAvailable(void) { return g_faFd >= 0; }

/* reference bases start..end (1-based, inclusive) into out (upper case);
   returns the count written */
static int refSeq(int chrom, long long start, long long end, char *out, int cap) {
    if (g_faFd < 0 || chrom < 1 || chrom > ST_NCHROM || !g_fai[chrom].ok) return 0;
    StFai *x = &g_fai[chrom];
    if (start < 1) start = 1;
    if (end > x->len) end = x->len;
    if (end < start) return 0;
    long long n = end - start + 1;
    if (n > cap) n = cap;
    long long s0 = start - 1, e0 = s0 + n - 1;
    long long o0 = x->off + s0 / x->lineBases * x->lineBytes + s0 % x->lineBases;
    long long o1 = x->off + e0 / x->lineBases * x->lineBytes + e0 % x->lineBases;
    size_t span = (size_t)(o1 - o0 + 1);
    char *buf = xalloc(span);
    ssize_t got = pread(g_faFd, buf, span, (off_t)o0);
    int k = 0;
    for (ssize_t i = 0; i < got && k < n; i++) {
        char c = buf[i];
        if (c == '\n' || c == '\r') continue;
        out[k++] = (char)toupper((unsigned char)c);
    }
    free(buf);
    return k;
}

/* ---- formatting --------------------------------------------------------
   Rotating scratch buffers: a returned string lives until 32 more calls. */
static char *scratch(void) {
    static char bufs[32][512];
    static int k;
    k = (k + 1) & 31;
    return bufs[k];
}
static const char *sstr(const char *fmt, ...) {
    char *b = scratch();
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(b, 512, fmt, ap);
    va_end(ap);
    return b;
}
static const char *fmtInt(long long v) {
    char tmp[32];
    snprintf(tmp, sizeof tmp, "%lld", v < 0 ? -v : v);
    char *b = scratch();
    int n = (int)strlen(tmp), w = 0;
    if (v < 0) b[w++] = '-';
    for (int i = 0; i < n; i++) {
        b[w++] = tmp[i];
        if ((n - i - 1) % 3 == 0 && i < n - 1) b[w++] = ',';
    }
    b[w] = 0;
    return b;
}
static const char *fmtBig(double v) {
    double a = fabs(v);
    if (a >= 1e9) return sstr("%.2fG", v / 1e9);
    if (a >= 1e6) return sstr("%.2fM", v / 1e6);
    if (a >= 1e4) return sstr("%.1fk", v / 1e3);
    return sstr("%.0f", v);
}
static const char *fmtBp(double v) {
    double a = fabs(v);
    if (a >= 1e6) return sstr("%.2f Mb", v / 1e6);
    if (a >= 1e3) return sstr("%.1f kb", v / 1e3);
    return sstr("%.0f bp", v);
}
static const char *fmtPos(int chrom, long long pos) { return sstr("chr%s:%s", chromName(chrom), fmtInt(pos)); }
static const char *fmtAf(double af) {
    if (isnan(af)) return "not seen";
    if (af <= 0) return "0";
    if (af >= 0.01) return sstr("%.1f%%", af * 100);
    return sstr("1 in %s", fmtInt((long long)llround(1.0 / af)));
}
static const char *zygName(int zyg) {
    switch (zyg) {
    case 1: return "Heterozygous";
    case 2: return "Homozygous";
    case 3: return "Hemizygous";
    case 4: return "Het (two alts)";
    default: return "Reference";
    }
}
static const char *zygShort(int zyg) {
    switch (zyg) {
    case 1: return "het";
    case 2: return "hom";
    case 3: return "hemi";
    case 4: return "het*";
    default: return "ref";
    }
}
static const char *clinName(int code) {
    /* 100 + class: a rare chip call set aside as unconfirmed (pipeline step 5c) */
    if (code > 100 && code <= 111) return sstr("Unconfirmed chip call (%s)", clinName(code - 100));
    switch (code) {
    case 1: return "Benign";
    case 2: return "Likely benign";
    case 3: return "Uncertain";
    case 4: return "Conflicting";
    case 5: return "Likely pathogenic";
    case 6: return "Pathogenic";
    case 7: return "Risk factor";
    case 8: return "Drug response";
    case 9: return "Protective";
    case 10: return "Association";
    case 11: return "Other";
    default: return "";
    }
}
static Color clinColor(int code) {
    if (code > 100) return TH->faint;
    switch (code) {
    case 1: case 2: return TH->ok;
    case 3: return TH->dim;
    case 4: return TH->warn;
    case 5: return (Color){ 222, 132, 88, 255 };
    case 6: return TH->bad;
    case 7: return TH->warn;
    case 8: return TH->info;
    case 9: return TH->ok;
    case 10: return TH->purple;
    default: return TH->faint;
    }
}
static Color impactColor(int impact) {
    switch (impact) {
    case 3: return TH->bad;
    case 2: return TH->warn;
    case 1: return TH->info;
    default: return TH->faint;
    }
}
static const char *csqLabel(const char *csq) {
    static const struct { const char *k, *v; } m[] = {
        { "stop_gained", "Stop gained" }, { "frameshift", "Frameshift" }, { "splice_acceptor", "Splice acceptor" },
        { "splice_donor", "Splice donor" }, { "start_lost", "Start lost" }, { "stop_lost", "Stop lost" },
        { "inframe_deletion", "In-frame deletion" }, { "inframe_insertion", "In-frame insertion" },
        { "inframe_altering", "In-frame change" }, { "missense", "Missense" }, { "splice_region", "Splice region" },
        { "coding_sequence", "Coding" }, { "feature_elongation", "Elongation" }, { "start_retained", "Start retained" },
        { "stop_retained", "Stop retained" }, { "synonymous", "Synonymous" }, { "5_prime_utr", "5' UTR" },
        { "3_prime_utr", "3' UTR" }, { "non_coding", "Non-coding RNA" }, { "intron", "Intron" },
    };
    if (!csq || !*csq) return "Intergenic";
    for (size_t i = 0; i < NELEMS(m); i++) if (!strcmp(m[i].k, csq)) return m[i].v;
    return csq;
}
static const char *kindName(int kind) {
    switch (kind) {
    case 0: return "SNV";
    case 1: return "Insertion";
    case 2: return "Deletion";
    default: return "Complex";
    }
}
/* parse "chr7:117,559,590", "7:117559590-117560000", "7 117559590" */
static bool parseLocus(const char *s, int *chrom, long long *a, long long *b) {
    char buf[128];
    int n = 0;
    for (const char *p = s; *p && n < 126; p++) if (*p != ',') buf[n++] = *p;
    buf[n] = 0;
    char *p = buf;
    while (*p == ' ') p++;
    char cname[16];
    int k = 0;
    if (!strncasecmp(p, "chr", 3)) p += 3;
    while (*p && *p != ':' && *p != ' ' && *p != '-' && k < 15) cname[k++] = *p++;
    cname[k] = 0;
    int c = chromIndex(cname);
    if (!c || (*p != ':' && *p != ' ')) return false;
    p++;
    char *e;
    long long x = strtoll(p, &e, 10);
    if (e == p) return false;
    long long y = x;
    while (*e == ' ') e++;
    if (*e == '-') { char *e2; y = strtoll(e + 1, &e2, 10); if (e2 == e + 1) y = x; }
    *chrom = c;
    *a = x < y ? x : y;
    *b = x < y ? y : x;
    return x > 0;
}

/* the found / not-found filter of the trait, well-known and remarkable lists */
static const char *const g_foundOpts[] = { "All", "Found in you", "Not found" };
/* a trait or well-known variant: does the person carry the effect allele? */
static bool traitFound(Rows *r, int i, int copiesCol, int resultCol) {
    if (!rnull(r, i, copiesCol)) return ri(r, i, copiesCol) > 0;
    const char *res = rc(r, i, resultCol);
    return strncasecmp(res, "not ", 4) && strncasecmp(res, "no ", 3) && strcmp(res, "PiMM");
}

/* ---- what kind of data this person has ---- */
static bool dkChip(void) { return !strcmp(dbMeta("call_model"), "chip"); }
static bool dkGvcf(void) { return !strcmp(dbMeta("call_model"), "gvcf"); }
static bool dkReads(void) {   /* aligned reads to show: a CRAM next to the database (none in files this edition imports) */
    static int has = -1;
    if (has < 0) has = access(sstr("%s/align/reads.cram", g_dataDir), R_OK) == 0;
    return has == 1;
}
static const char *dkName(void) { return *dbMeta("seq_kind") ? dbMeta("seq_kind") : "DNA data"; }
