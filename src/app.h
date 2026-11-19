/* ========================================================================
   src/app.h - Strand's shared types and every cross-module prototype.

   game.c includes the engine, then this header, then each src/ module, so
   include order among the modules does not matter (a module may call a
   function defined in one included later). Strand's own symbols carry the
   st / g_st prefix; widgets are gx*, drawing helpers dr*.
   ======================================================================== */

#include <sqlite3.h>
#include <math.h>
#include <stdatomic.h>

#define ST_VERSION "1.2.0"

/* ---- pages ---------------------------------------------------------- */
typedef enum {
    PG_OVERVIEW, PG_BROWSER, PG_VARIANTS, PG_GENE, PG_HEALTH, PG_CARRIER, PG_PGX, PG_TRAITS,
    PG_PRS, PG_ANCESTRY, PG_LINEAGE, PG_RAW, PG_ABOUT, PG_RARE, PG_RARITY, PG_PEOPLE, PG_BRIEF, PG_DIET, PG_REPORT, PG__N
} StPage;

/* ---- fonts ---------------------------------------------------------- */
typedef enum { F_REG, F_SEMI, F_MONO, F__N } StFace;
static void fontsInit(void);
static void fontsReload(float scale);
static Font *fontAt(StFace f, int px);
static int  txtW(StFace f, int px, const char *s);
static int  txtWn(StFace f, int px, const char *s, int n);
static void txt(StFace f, int px, float x, float y, Color c, const char *s);
static void txtn(StFace f, int px, float x, float y, Color c, const char *s, int n);
static void txtf(StFace f, int px, float x, float y, Color c, const char *fmt, ...) __attribute__((format(printf, 6, 7)));
static void txtIn(StFace f, int px, Rectangle r, Color c, int align, const char *s);
static int  txtWrap(StFace f, int px, float w, const char *s, int (*line)(const char *s, int n, void *u), void *u);
static float txtPara(StFace f, int px, float x, float y, float w, float lh, Color c, const char *s, int maxLines);
static float txtParaH(StFace f, int px, float w, float lh, const char *s, int maxLines);

/* ---- theme ---------------------------------------------------------- */
typedef struct {
    Color bg, side, panel, panel2, card, line, lineHi, text, dim, faint, accent, accent2, accentDim,
          ok, warn, bad, info, purple, hover, press, field, sel, selDim, shadow, onAccent, zebra;
    int fs, fsSmall, fsTiny, fsH1, fsH2, fsH3, rowH, pad;
} StTheme;
static StTheme g_t;
#define TH (&g_t)

/* ---- query results (one malloc'd block: cell pointers then strings) --- */
typedef struct {
    int ncol, nrow;
    char **cell;      /* nrow * ncol, row-major; NULL for SQL NULL */
    char **colName;   /* ncol */
    double ms;        /* how long the query took */
    char err[200];
} Rows;
static const char *rc(const Rows *r, int row, int col);      /* "" for NULL */
static double rf(const Rows *r, int row, int col);           /* NAN for NULL */
static long long ri(const Rows *r, int row, int col);        /* 0 for NULL */
static bool rnull(const Rows *r, int row, int col);
static void rowsFree(Rows **r);

/* ---- database ------------------------------------------------------- */
typedef enum {
    SLOT_VARTABLE, SLOT_VARCOUNT, SLOT_BROWSE_VAR, SLOT_BROWSE_GENE, SLOT_BROWSE_COV, SLOT_BROWSE_BLOCK,
    SLOT_SEARCH, SLOT_GENE, SLOT_GENEVARS, SLOT_DETAIL, SLOT_DETAIL2, SLOT_HEALTH, SLOT_CARRIER, SLOT_PAGE,
    SLOT_PAGE2, SLOT_PAGE3, SLOT_GWAS, SLOT__N
} StSlot;
static bool dbOpen(const char *path);
static Rows *dbq(const char *sql, ...);          /* synchronous, main thread; NULL-terminated text params */
static Rows *dbqv(const char *sql, int n, const char *const *params);
static const char *dbMeta(const char *key);      /* cached meta table; "" if absent */
static bool dbHasTable(const char *name);
/* the open person's kind of data (meta call_model): a genotyping chip, a VCF
   of differences only, or a gVCF with reference blocks; and whether reads exist */
static bool dkChip(void);
static bool dkGvcf(void);
static bool dkReads(void);
static const char *dkName(void);   /* "23andMe chip", "whole genome (gVCF)", ... */
static bool dbColumnExists(const char *table, const char *col);
static const char *pctLabel(double pct);   /* "37th", "top 1%" */
static void dbAsync(StSlot slot, const char *sql, int n, const char *const *params);
static Rows *dbTake(StSlot slot);               /* finished result for this slot (caller frees), else NULL */
static bool dbBusy(StSlot slot);
static void dbPoll(void);
static void dbShutdown(void);

/* ---- genome model --------------------------------------------------- */
#define ST_NCHROM 25
typedef struct { char name[4]; long long len; } StChrom;
typedef struct { int chrom; long long start, end; char name[16]; char stain[10]; } StBand;
static StChrom g_chrom[ST_NCHROM + 1];   /* 1-based */
static StBand *g_bands;
static int g_bandN;
static int chromIndex(const char *name);  /* "chr7" "7" "X" -> 1..25, 0 if unknown */
static const char *chromName(int idx);
static void genomeInit(void);
static int  refSeq(int chrom, long long start, long long end, char *out, int cap);   /* 1-based inclusive */
static bool refAvailable(void);

/* ---- navigation ----------------------------------------------------- */
typedef struct {
    StPage page;
    int geneId;               /* PG_GENE */
    long long varId;          /* selected variant */
    int chrom; long long start, end;  /* PG_BROWSER locus */
} StLoc;
static void navGo(StLoc l);
static void navPage(StPage p);
static void navGene(int geneId);
static void navGeneSym(const char *sym);
static void stEnter(const char *name);      /* go to a page by name, "page[:arg]" (as --enter) */
static void rareFocus(const char *key);       /* open the remarkable entry with this key */
static void navLocus(int chrom, long long start, long long end);
static void navVariant(long long varId);
static StLoc *navCur(void);
static void toast(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void openUrl(const char *url);
static void copyText(const char *s);
static void exportReport(void);
static bool reportDocsDir(char *out, size_t n);
static void variantsShow(int scope, const char *text, long long sel);

/* ---- shared formatting ---------------------------------------------- */
static const char *fmtInt(long long v);                  /* 1,234,567 (rotating buffers) */
static const char *fmtBig(double v);                     /* 1.2M, 3.4G */
static const char *fmtPos(int chrom, long long pos);     /* chr7:117,559,590 */
static const char *fmtAf(double af);                     /* 0.12% / 1 in 4,200 / not seen */
static const char *zygName(int zyg);
static const char *clinName(int code);
static Color clinColor(int code);
static Color impactColor(int impact);
static const char *csqLabel(const char *csq);
static Color phenoColor(const char *ph);
static const char *sstr(const char *fmt, ...) __attribute__((format(printf, 1, 2)));  /* rotating scratch string */

/* ---- pages ---------------------------------------------------------- */
static void pageOverview(Rectangle r);
static void pageBrowser(Rectangle r);
static void pageVariants(Rectangle r);
static void pageGene(Rectangle r);
static void pageHealth(Rectangle r);
static void pageCarrier(Rectangle r);
static void pagePgx(Rectangle r);
static void pageTraits(Rectangle r);
static void pagePrs(Rectangle r);
static void pageAncestry(Rectangle r);
static void pageLineage(Rectangle r);
static void pageRaw(Rectangle r);
static void pageAbout(Rectangle r);
static void variantPanel(Rectangle r, long long varId, bool showBrowse);   /* the variant detail card */
static float variantPanelH(void);
static void pageEnter(StPage p);    /* a page became current */
static void pageMissing(Rectangle r, const char *title, const char *why);
static float pageHeader(Rectangle r, const char *title, const char *sub);

/* ---- reads (CRAM via samtools) --------------------------------------- */
static bool readsAvailable(void);
static void readsRequest(int chrom, long long start, long long end);
static int  readsState(void);    /* 0 none, 1 loading, 2 ready, 3 failed */
