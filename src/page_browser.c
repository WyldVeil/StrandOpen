/* ========================================================================
   src/page_browser.c - the genome browser.

   One chromosome window at a time, from a whole chromosome down to single
   bases. Tracks, top to bottom: the chromosome's ideogram with the window
   boxed, a ruler, coverage, the sample's variants, genes (canonical
   transcripts: exons, coding parts, strand), the reference sequence when a
   base is wide enough to read, and the aligned reads when the window is
   small enough to fetch them. What is drawn depends on the window's width:

     wider than 3 Mb    variant density and coverage from 10 kb bins, genes
                        as plain spans (protein-coding only)
     3 Mb to 200 kb     each variant, gene structure, binned coverage
     under 200 kb       coverage from the gVCF's reference blocks
     under 25 kb        reads (when the CRAM exists)
     about 150 bp       the bases themselves, in reads and reference

   The view is a continuous [a, b) over 1-based positions; base p spans
   [p, p + 1). Data is fetched for three window widths around the view,
   after the view settles for a moment, so panning shows what is already
   loaded and the edges fill in.
   ======================================================================== */

typedef struct {
    int chrom;
    double a, b;
    StLoc applied;              /* the nav location last turned into a view */
    /* what has been fetched */
    int fChrom; double fa, fb; int fMode;
    bool fetchWanted; double changedAt;
    Rows *vars, *genes, *exons, *cov, *blocks;
    long long selVar;
    int selGene;
    bool dragging, dragMoved;
    float dragX0;
    double dragA0, dragB0;
    char locus[160];
    GxId locusId;
    bool showReads, showSeq, showCov, onlyCoding;
    float readsScroll;
    double readsAskA, readsAskB; int readsAskChrom;
    double readsWantAt;
} Browser;
static Browser g_br = { .selVar = -1, .selGene = -1, .showReads = true, .showSeq = true, .showCov = true };

static const Color g_baseCol[4] = { { 108, 191, 108, 255 }, { 96, 150, 220, 255 }, { 230, 162, 60, 255 }, { 224, 93, 93, 255 } };
static Color baseColor(char b) {
    int i = baseIdx(b);
    return i >= 0 ? g_baseCol[i] : TH->faint;
}

static double brW(void) { return g_br.b - g_br.a; }
static void brClamp(void) {
    double len = (double)g_chrom[g_br.chrom].len + 1;
    double w = g_br.b - g_br.a;
    if (w < 40) { double c = (g_br.a + g_br.b) / 2; g_br.a = c - 20; g_br.b = c + 20; w = 40; }
    if (w > len) { g_br.a = 1; g_br.b = len; return; }
    if (g_br.a < 1) { g_br.b += 1 - g_br.a; g_br.a = 1; }
    if (g_br.b > len) { g_br.a -= g_br.b - len; g_br.b = len; }
}
static void brSet(int chrom, double a, double b) {
    if (chrom < 1 || chrom > ST_NCHROM) return;
    g_br.chrom = chrom;
    g_br.a = a;
    g_br.b = b;
    brClamp();
    g_br.changedAt = (GetTime)();
    g_br.fetchWanted = true;
    /* keep the history entry current (panning is not a new step) */
    StLoc *l = navCur();
    l->chrom = chrom;
    l->start = (long long)g_br.a;
    l->end = (long long)g_br.b;
    g_br.applied = *l;
    if (g_gx.editId != g_br.locusId)
        snprintf(g_br.locus, sizeof g_br.locus, "chr%s:%s-%s", chromName(chrom), fmtInt((long long)g_br.a), fmtInt((long long)g_br.b - 1));
}
static void brZoom(double factor, double at) {
    double a = at - (at - g_br.a) * factor, b = at + (g_br.b - at) * factor;
    brSet(g_br.chrom, a, b);
}

static int brMode(void) {
    double w = brW();
    return w > 3e6 ? 0 : w > 2e5 ? 1 : 2;
}
static void brFetch(void) {
    double w = brW();
    double fa = fmax(1, g_br.a - w), fb = g_br.b + w;
    int mode = brMode();
    g_br.fChrom = g_br.chrom; g_br.fa = fa; g_br.fb = fb; g_br.fMode = mode;
    char c[8], sa[24], sb[24], ba[24], bb[24];
    snprintf(c, sizeof c, "%d", g_br.chrom);
    snprintf(sa, sizeof sa, "%lld", (long long)fa);
    snprintf(sb, sizeof sb, "%lld", (long long)fb);
    snprintf(ba, sizeof ba, "%lld", (long long)fa / 10000);
    snprintf(bb, sizeof bb, "%lld", (long long)fb / 10000);
    const char *p3[] = { c, sa, sb };
    const char *pb[] = { c, ba, bb };
    if (mode >= 1)
        dbAsync(SLOT_BROWSE_VAR, "SELECT id, pos, end, ref, alt, zyg, impact, clin, gene, hgvs_p, csq, pass, dp, rsid, af "
                                 "FROM variant WHERE chrom=?1 AND pos BETWEEN ?2 AND ?3 ORDER BY pos LIMIT 60000", 3, p3);
    else
        rowsFree(&g_br.vars);
    if (mode == 0)
        dbAsync(SLOT_BROWSE_GENE, "SELECT id, symbol, start, end, strand, biotype, cds_start, cds_end FROM gene "
                                  "WHERE chrom=?1 AND start<=?3 AND end>=?2 AND biotype='protein_coding' ORDER BY start LIMIT 4000", 3, p3);
    else
        dbAsync(SLOT_BROWSE_GENE, "SELECT id, symbol, start, end, strand, biotype, cds_start, cds_end FROM gene "
                                  "WHERE chrom=?1 AND start<=?3 AND end>=?2 ORDER BY biotype!='protein_coding', start LIMIT 4000", 3, p3);
    if (mode >= 1)
        dbAsync(SLOT_DETAIL2, "SELECT e.gene, e.start, e.end, e.kind FROM exon e JOIN gene g ON g.id=e.gene "
                              "WHERE g.chrom=?1 AND g.start<=?3 AND g.end>=?2 ORDER BY e.gene, e.start", 3, p3);
    dbAsync(SLOT_BROWSE_COV, "SELECT bin, depth, callable, nhet+nhom+nhemi FROM covbin WHERE chrom=?1 AND bin BETWEEN ?2 AND ?3 ORDER BY bin", 3, pb);
    if (mode == 2)
        dbAsync(SLOT_BROWSE_BLOCK, "SELECT start, end, dp, pass FROM refblock WHERE chrom=?1 AND start BETWEEN ?2 - 300000 AND ?3 "
                                   "AND end >= ?2 ORDER BY start", 3, p3);
    else
        rowsFree(&g_br.blocks);
}
static void brTakeResults(void) {
    Rows *t;
    if ((t = dbTake(SLOT_BROWSE_VAR))) { rowsFree(&g_br.vars); g_br.vars = t; }
    if ((t = dbTake(SLOT_BROWSE_GENE))) { rowsFree(&g_br.genes); g_br.genes = t; }
    if ((t = dbTake(SLOT_DETAIL2))) { rowsFree(&g_br.exons); g_br.exons = t; }
    if ((t = dbTake(SLOT_BROWSE_COV))) { rowsFree(&g_br.cov); g_br.cov = t; }
    if ((t = dbTake(SLOT_BROWSE_BLOCK))) { rowsFree(&g_br.blocks); g_br.blocks = t; }
}

/* ---- geometry ---- */
static Rectangle g_brArea;    /* the tracks' horizontal extent */
static float brX(double pos) { return g_brArea.x + (float)((pos - g_br.a) / brW() * g_brArea.width); }
static double brPosAt(float x) { return g_br.a + (double)(x - g_brArea.x) / g_brArea.width * brW(); }

static double niceStep(double span, int ticks) {
    double raw = span / ticks, p = pow(10, floor(log10(raw)));
    double m = raw / p;
    return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}
static const char *tickLabel(double v, double step) {
    if (step >= 1e6) return sstr("%.0f Mb", v / 1e6);
    if (step >= 1e5) return sstr("%.1f Mb", v / 1e6);
    if (step >= 1e3) return sstr("%s kb", fmtInt((long long)(v / 1e3)));
    return fmtInt((long long)v);
}

/* ---- tracks ---- */
static void trackLabel(float y, const char *name) {
    txt(F_SEMI, TH->fsTiny, g_brArea.x - 104, y, TH->faint, name);
}
static void drawIdeoNav(Rectangle r) {
    Rectangle ir = { r.x, r.y + 4, r.width, 14 };
    ideogramH(ir, g_br.chrom);
    float x0 = ideoX(ir, g_br.chrom, (long long)g_br.a), x1 = ideoX(ir, g_br.chrom, (long long)g_br.b);
    if (x1 - x0 < 3) { float m = (x0 + x1) / 2; x0 = m - 1.5f; x1 = m + 1.5f; }
    drRect((Rectangle){ x0, ir.y - 3, x1 - x0, ir.height + 6 }, drFade(TH->bad, 0.25f));
    DrawRectangleLinesEx((Rectangle){ x0, ir.y - 3, x1 - x0, ir.height + 6 }, 1.5f, TH->bad);
    bool h;
    GxId id = gxId("ideo nav");
    gxBehave(id, (Rectangle){ r.x, r.y, r.width, r.height }, &h, NULL);
    if (h) {
        g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
        long long p = (long long)((g_gx.mouse.x - ir.x) / ir.width * (double)g_chrom[g_br.chrom].len);
        const char *band = "";
        for (int j = 0; j < g_bandN; j++) if (g_bands[j].chrom == g_br.chrom && p >= g_bands[j].start && p <= g_bands[j].end) band = g_bands[j].name;
        gxTip(gxId("ideo tip"), r, sstr("%s%s  %s\nClick to jump here", chromName(g_br.chrom), band, fmtBp((double)p)));
    }
    if (g_gx.active == id && g_gx.down[0]) {
        double p = (g_gx.mouse.x - ir.x) / ir.width * (double)g_chrom[g_br.chrom].len;
        double w = brW();
        brSet(g_br.chrom, p - w / 2, p + w / 2);
    }
    /* band names under the bar at wide zoom */
    txt(F_REG, TH->fsTiny, r.x, r.y + 22, TH->faint, sstr("chr%s", chromName(g_br.chrom)));
}
static void drawRuler(Rectangle r) {
    double step = niceStep(brW(), (int)(r.width / 110));
    double first = ceil(g_br.a / step) * step;
    drHLine(r.x, r.y + r.height - 1, r.width, TH->line);
    for (double v = first; v < g_br.b; v += step) {
        float x = brX(v);
        drVLine(x, r.y + r.height - 7, 7, TH->lineHi);
        const char *s = tickLabel(v, step);
        txt(F_REG, TH->fsTiny, x - (float)txtW(F_REG, TH->fsTiny, s) / 2, r.y + 2, TH->dim, s);
    }
    const char *w = fmtBp(brW());
    txt(F_SEMI, TH->fsTiny, r.x + r.width - (float)txtW(F_SEMI, TH->fsTiny, w), r.y - 14, TH->accent, w);
}

static void drawCoverage(Rectangle r, ReadsResult *rd) {
    trackLabel(r.y + 2, "COVERAGE");
    drRect(r, (Color){ 24, 27, 34, 255 });
    float maxDp = 60;
    double w = brW();
    Color fill = drFade(TH->accent2, 0.55f);
    if (rd && rd->chrom == g_br.chrom && rd->start <= g_br.a + 1 && rd->end >= g_br.b - 1 && w <= 25000) {
        int mx = 1;
        for (int i = 0; i < rd->end - rd->start + 1; i++) mx = rd->depth[i] > mx ? rd->depth[i] : mx;
        maxDp = fmaxf(30, (float)mx * 1.1f);
        for (long long p = (long long)g_br.a; p <= (long long)g_br.b; p++) {
            int k = (int)(p - rd->start);
            if (k < 0 || k > rd->end - rd->start) continue;
            float x0 = floorf(brX((double)p)), x1 = floorf(brX((double)p + 1));
            if (x1 <= x0) x1 = x0 + 1;
            int d = rd->depth[k];
            float h = r.height * fminf(1, (float)d / maxDp);
            /* allele-aware: alt share coloured when it is above a fifth */
            char rb = rd->ref[k];
            int altN = 0, tot = 0;
            for (int b = 0; b < 4; b++) { tot += rd->acgt[k][b]; if ("ACGT"[b] != rb) altN += rd->acgt[k][b]; }
            if (tot > 0 && altN * 5 >= tot) {
                float y = r.y + r.height;
                for (int b = 0; b < 4; b++) {
                    float hh = h * (float)rd->acgt[k][b] / (float)tot;
                    drRect((Rectangle){ x0, y - hh, x1 - x0, hh }, "ACGT"[b] == rb ? fill : g_baseCol[b]);
                    y -= hh;
                }
            } else drRect((Rectangle){ x0, r.y + r.height - h, x1 - x0, h }, fill);
        }
        txt(F_REG, TH->fsTiny, r.x + 4, r.y + 2, TH->faint, sstr("reads, 0-%.0f", maxDp));
    } else if (g_br.blocks && w <= 2e5) {
        maxDp = 40;
        for (int i = 0; i < g_br.blocks->nrow; i++) {
            double s = (double)ri(g_br.blocks, i, 0), e = (double)ri(g_br.blocks, i, 1) + 1;
            if (e < g_br.a || s > g_br.b) continue;
            float x0 = fmaxf(r.x, brX(s)), x1 = fminf(r.x + r.width, brX(e));
            float d = (float)ri(g_br.blocks, i, 2);
            float h = r.height * fminf(1, d / maxDp);
            drRect((Rectangle){ x0, r.y + r.height - h, fmaxf(1, x1 - x0), h }, ri(g_br.blocks, i, 3) ? fill : drFade(TH->warn, 0.5f));
        }
        if (g_br.vars) for (int i = 0; i < g_br.vars->nrow; i++) {
            double p = (double)ri(g_br.vars, i, 1);
            if (p < g_br.a || p > g_br.b) continue;
            float d = (float)ri(g_br.vars, i, 12);
            float h = r.height * fminf(1, d / maxDp);
            float x0 = brX(p), x1 = brX(p + 1);
            drRect((Rectangle){ x0, r.y + r.height - h, fmaxf(1.5f, x1 - x0), h }, drFade(TH->accent, 0.8f));
        }
        txt(F_REG, TH->fsTiny, r.x + 4, r.y + 2, TH->faint, dkGvcf() ? "gVCF minimum depth, 0-40 (amber: low quality)"
                                                              : dkChip() ? "chip: tested positions only, no depth" : "no read depth in this file");
    } else if (g_br.cov) {
        maxDp = 40;
        for (int i = 0; i < g_br.cov->nrow; i++) {
            double s = (double)ri(g_br.cov, i, 0) * 10000 + 1, e = s + 10000;
            if (e < g_br.a || s > g_br.b) continue;
            float x0 = fmaxf(r.x, brX(s)), x1 = fminf(r.x + r.width, brX(e));
            float d = (float)rf(g_br.cov, i, 1), cf = (float)rf(g_br.cov, i, 2);
            float h = r.height * fminf(1, d * cf / maxDp);
            drRect((Rectangle){ x0, r.y + r.height - h, fmaxf(1, x1 - x0 + 0.5f), h }, cf < 0.5f ? drFade(TH->warn, 0.6f) : fill);
        }
        txt(F_REG, TH->fsTiny, r.x + 4, r.y + 2, TH->faint, "mean depth per 10 kb, 0-40");
    }
}

static long long g_brHoverVar = -1;
static void drawVariants(Rectangle r) {
    trackLabel(r.y + 2, "VARIANTS");
    g_brHoverVar = -1;
    double w = brW();
    if (brMode() == 0 || !g_br.vars) {
        /* density from the bins */
        if (!g_br.cov) return;
        float mx = 1;
        for (int i = 0; i < g_br.cov->nrow; i++) mx = fmaxf(mx, (float)ri(g_br.cov, i, 3));
        for (int i = 0; i < g_br.cov->nrow; i++) {
            double s = (double)ri(g_br.cov, i, 0) * 10000 + 1, e = s + 10000;
            if (e < g_br.a || s > g_br.b) continue;
            float x0 = fmaxf(r.x, brX(s)), x1 = fminf(r.x + r.width, brX(e));
            float v = (float)ri(g_br.cov, i, 3) / mx;
            float h = (r.height - 6) * v;
            drRect((Rectangle){ x0, r.y + r.height - 3 - h, fmaxf(1, x1 - x0 + 0.5f), h }, heatColor(v * 1.2f));
        }
        txt(F_REG, TH->fsTiny, r.x + 4, r.y + 2, TH->faint, "variants per 10 kb");
        return;
    }
    float pxPerBase = (float)(r.width / w);
    float mid = r.y + r.height / 2;
    drHLine(r.x, mid, r.width, drFade(TH->line, 0.6f));
    float bestD = 6;
    for (int i = 0; i < g_br.vars->nrow; i++) {
        double p = (double)ri(g_br.vars, i, 1);
        if (p + 1 < g_br.a || p > g_br.b) continue;
        long long id = ri(g_br.vars, i, 0);
        int zyg = (int)ri(g_br.vars, i, 5), imp = (int)ri(g_br.vars, i, 6), clin = (int)ri(g_br.vars, i, 7);
        float x0 = brX(p), x1 = brX((double)ri(g_br.vars, i, 2) + 1);
        float wpx = fmaxf(2, x1 - x0 - (pxPerBase > 4 ? 1 : 0));
        Color c = imp >= 3 ? TH->bad : imp == 2 ? TH->warn : imp == 1 ? TH->info : (Color){ 120, 150, 190, 255 };
        if (!ri(g_br.vars, i, 11)) c = drFade(c, 0.45f);
        float top = r.y + 6, h = r.height - 12;
        if (zyg == 1 || zyg == 4) {   /* het: the lower half filled */
            drRect((Rectangle){ x0, top, wpx, h }, drFade(c, 0.25f));
            drRect((Rectangle){ x0, top + h / 2, wpx, h / 2 }, c);
        } else drRect((Rectangle){ x0, top, wpx, h }, c);
        if (clin == 5 || clin == 6) DrawRectangleLinesEx((Rectangle){ x0 - 2, top - 3, wpx + 4, h + 6 }, 1.5f, TH->bad);
        if (id == g_br.selVar) DrawRectangleLinesEx((Rectangle){ x0 - 3, r.y + 1, wpx + 6, r.height - 2 }, 2, TH->text);
        float d = fabsf(g_gx.mouse.x - (x0 + wpx / 2));
        if (CheckCollisionPointRec(g_gx.mouse, r) && d < bestD) { bestD = d; g_brHoverVar = i; }
    }
    if (g_brHoverVar >= 0) {
        int i = (int)g_brHoverVar;
        const char *gene = rc(g_br.vars, i, 8);
        char tip[400];
        snprintf(tip, sizeof tip, "%s %s>%s  %s\n%s%s%s  %s\n%s%s",
                 fmtPos(g_br.chrom, ri(g_br.vars, i, 1)), rc(g_br.vars, i, 3), rc(g_br.vars, i, 4), zygName((int)ri(g_br.vars, i, 5)),
                 *gene ? gene : "intergenic", *gene ? " " : "", rc(g_br.vars, i, 9), csqLabel(rc(g_br.vars, i, 10)),
                 ri(g_br.vars, i, 13) ? sstr("rs%lld  ", ri(g_br.vars, i, 13)) : "",
                 ri(g_br.vars, i, 7) ? clinName((int)ri(g_br.vars, i, 7)) : "");
        gxTip(gxId("var tip") ^ (GxId)ri(g_br.vars, i, 0), r, tip);
        g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
    }
}

/* genes: packed into rows by their drawn extent (name included) */
typedef struct { int idx; float x0, x1; int row; } GeneLay;
static GeneLay g_gl[4000];
static int g_glN, g_glRows;
static int g_brHoverGene = -1;
static void layoutGenes(Rectangle r) {
    g_glN = 0;
    g_glRows = 0;
    if (!g_br.genes) return;
    float rowEnd[12];
    for (int i = 0; i < g_br.genes->nrow && g_glN < 4000; i++) {
        double s = (double)ri(g_br.genes, i, 2), e = (double)ri(g_br.genes, i, 3) + 1;
        if (e < g_br.a || s > g_br.b) continue;
        if (brMode() >= 1 && strcmp(rc(g_br.genes, i, 5), "protein_coding") && brW() > 1.5e6) continue;
        float x0 = brX(s), x1 = brX(e);
        float lw = (float)txtW(F_SEMI, TH->fsTiny, rc(g_br.genes, i, 1)) + 10;
        float cx0 = fminf(x0, fmaxf(r.x, x0)), cx1 = fmaxf(x1, fmaxf(x0, r.x) + lw);
        int k = 0;
        while (k < g_glRows && rowEnd[k] > cx0 - 6) k++;
        if (k >= 12) continue;
        if (k == g_glRows) g_glRows++;
        rowEnd[k] = cx1;
        g_gl[g_glN++] = (GeneLay){ i, x0, x1, k };
    }
}
static float genesH(void) { return 16 + (float)(g_glRows > 0 ? g_glRows : 1) * 32; }
static void drawGenes(Rectangle r) {
    trackLabel(r.y + 2, "GENES");
    g_brHoverGene = -1;
    if (!g_br.genes) return;
    int e0 = 0;
    for (int k = 0; k < g_glN; k++) {
        GeneLay *L = &g_gl[k];
        int i = L->idx;
        int gid = (int)ri(g_br.genes, i, 0);
        bool coding = !strcmp(rc(g_br.genes, i, 5), "protein_coding");
        Color c = coding ? (Color){ 110, 178, 206, 255 } : (Color){ 150, 140, 190, 255 };
        if (gid == g_br.selGene) c = TH->warn;
        float y = r.y + 10 + (float)L->row * 32;
        float x0 = fmaxf(r.x, L->x0), x1 = fminf(r.x + r.width, L->x1);
        float mid = y + 8;
        drRect((Rectangle){ x0, mid - 0.5f, fmaxf(1, x1 - x0), 1.5f }, drFade(c, 0.8f));
        /* exons */
        bool drew = false;
        if (g_br.exons) {
            while (e0 < g_br.exons->nrow && ri(g_br.exons, e0, 0) < gid) e0++;
            /* the exon rows are ordered by gene id; genes here are by start, so search */
            int lo = 0, hi = g_br.exons->nrow;
            while (lo < hi) { int m = (lo + hi) / 2; if (ri(g_br.exons, m, 0) < gid) lo = m + 1; else hi = m; }
            for (int j = lo; j < g_br.exons->nrow && ri(g_br.exons, j, 0) == gid; j++) {
                double s = (double)ri(g_br.exons, j, 1), e = (double)ri(g_br.exons, j, 2) + 1;
                if (e < g_br.a || s > g_br.b) continue;
                float ex0 = fmaxf(r.x, brX(s)), ex1 = fminf(r.x + r.width, brX(e));
                if (ex1 <= ex0) continue;
                bool cds = ri(g_br.exons, j, 3) == 1;
                float hh = cds ? 14 : 8;
                drRect((Rectangle){ ex0, mid - hh / 2, fmaxf(1, ex1 - ex0), hh }, cds ? c : drFade(c, 0.7f));
                drew = true;
            }
        }
        if (!drew && brMode() == 0) drRect((Rectangle){ x0, mid - 4, fmaxf(1, x1 - x0), 8 }, drFade(c, 0.85f));
        /* strand chevrons along the intron line */
        int strand = (int)ri(g_br.genes, i, 4);
        if (x1 - x0 > 30) {
            for (float x = x0 + 12; x < x1 - 8; x += 40) {
                float d = strand > 0 ? 3 : -3;
                drLine(x - d, mid - 3, x + d, mid, 1.2f, drFade(c, 0.9f));
                drLine(x + d, mid, x - d, mid + 3, 1.2f, drFade(c, 0.9f));
            }
        }
        const char *name = rc(g_br.genes, i, 1);
        float lx = fmaxf(r.x + 2, L->x0);
        if (lx + 40 < r.x + r.width) txt(F_SEMI, TH->fsTiny, lx, mid + 9, gid == g_br.selGene ? TH->warn : TH->dim, name);
        Rectangle hit = { x0, y, fmaxf(6, x1 - x0), 30 };
        if (CheckCollisionPointRec(g_gx.mouse, hit) && gxInClip(g_gx.mouse) && g_brHoverVar < 0) {
            g_brHoverGene = i;
            gxTip(gxId("gene tip") ^ (GxId)gid, hit, sstr("%s  (%s strand)\n%s\nchr%s:%s-%s\nClick to select, double-click for the gene page",
                  name, strand > 0 ? "+" : "-", rc(g_br.genes, i, 5), chromName(g_br.chrom), fmtInt(ri(g_br.genes, i, 2)), fmtInt(ri(g_br.genes, i, 3))));
        }
    }
}

static void drawSequence(Rectangle r) {
    trackLabel(r.y + 4, "SEQUENCE");
    long long a = (long long)floor(g_br.a), b = (long long)ceil(g_br.b);
    if (b - a > 2000) return;
    char seq[2100];
    int n = refSeq(g_br.chrom, a, b, seq, (int)sizeof seq - 1);
    float pxb = (float)(g_brArea.width / brW());
    for (int i = 0; i < n; i++) {
        float x0 = brX((double)(a + i)), x1 = brX((double)(a + i + 1));
        Color c = baseColor(seq[i]);
        if (pxb >= 7) {
            drRect((Rectangle){ x0 + 0.5f, r.y + 2, x1 - x0 - 1, r.height - 4 }, drFade(c, 0.2f));
            char s[2] = { seq[i], 0 };
            txtIn(F_MONO, pxb >= 14 ? 13 : pxb >= 9 ? 11 : 9, (Rectangle){ x0 - 2, r.y, x1 - x0 + 4, r.height }, c, 1, s);
        } else drRect((Rectangle){ x0, r.y + 4, fmaxf(1, x1 - x0), r.height - 8 }, c);
    }
}

static void drawReads(Rectangle r, ReadsResult *rd) {
    trackLabel(r.y + 2, "READS");
    double w = brW();
    if (!readsAvailable()) {
        txt(F_REG, TH->fsSmall, r.x + 4, r.y + 4, TH->faint, dkReads() ? "Reads appear here once the import has aligned them." : "No sequencing reads for this person: their data has calls only.");
        return;
    }
    if (w > 25000) { txt(F_REG, TH->fsSmall, r.x + 4, r.y + 4, TH->faint, "Zoom in below 25 kb to see the reads."); return; }
    if (!rd || rd->chrom != g_br.chrom) {
        txt(F_REG, TH->fsSmall, r.x + 4, r.y + 4, TH->faint, readsState() == 1 ? "Loading reads..." : "");
        return;
    }
    if (rd->err[0]) { txt(F_REG, TH->fsSmall, r.x + 4, r.y + 4, TH->bad, rd->err); return; }
    float pxb = (float)(r.width / w);
    /* rows that hold a read inside the view, renumbered without gaps */
    static int remap[2048];
    int nvis = 0;
    for (int k = 0; k < rd->rows && k < 2048; k++) remap[k] = -1;
    for (int i = 0; i < rd->n; i++) {
        ReadRec *q = &rd->r[i];
        if (q->row < 0 || q->row >= 2048 || q->end + 1 < g_br.a || q->start > g_br.b) continue;
        remap[q->row] = 0;
    }
    for (int k = 0; k < rd->rows && k < 2048; k++) if (remap[k] == 0) remap[k] = ++nvis;
    int rowsShown = nvis;
    float rowH = rowsShown * 11 <= r.height ? 10 : rowsShown * 6 <= r.height ? 5 : 3;
    float gapH = rowH >= 10 ? 1 : 0.5f;
    int visRows = (int)(r.height / (rowH + gapH));
    float maxScroll = fmaxf(0, (float)(rowsShown - visRows) * (rowH + gapH));
    if (gxHover(r) && g_gx.wheel != 0 && (IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT))) g_br.readsScroll -= g_gx.wheel * 60;
    g_br.readsScroll = clampf(g_br.readsScroll, 0, maxScroll);
    gxClipPush(r);
    Color fwd = { 112, 120, 136, 255 }, rev = { 96, 104, 124, 255 };
    for (int i = 0; i < rd->n; i++) {
        ReadRec *q = &rd->r[i];
        if (q->row < 0 || q->row >= 2048) continue;
        if (q->end + 1 < g_br.a || q->start > g_br.b) continue;
        float y = r.y + (float)(remap[q->row] - 1) * (rowH + gapH) - g_br.readsScroll;
        if (y + rowH < r.y || y > r.y + r.height) continue;
        float x0 = brX(q->start), x1 = brX((double)q->end + 1);
        Color c = q->rev ? rev : fwd;
        if (q->mapq == 0) c = (Color){ 60, 64, 74, 255 };
        if (q->dup) c = drFade(c, 0.35f);
        drRect((Rectangle){ x0, y, x1 - x0, rowH }, c);
        if (pxb > 2 && rowH >= 10) {   /* strand point */
            if (q->rev) DrawTriangle((Vector2){ x0, y + rowH / 2 }, (Vector2){ x0 + 4, y + rowH }, (Vector2){ x0 + 4, y }, c);
            else DrawTriangle((Vector2){ x1 + 4, y + rowH / 2 }, (Vector2){ x1, y }, (Vector2){ x1, y + rowH }, c);
        }
        for (int e = 0; e < q->evN; e++) {
            ReadEv *ev = &rd->ev[q->ev0 + e];
            if (ev->kind == 'X') {
                float bx0 = brX(ev->pos), bx1 = brX((double)ev->pos + 1);
                Color bc = baseColor(ev->base);
                if (q->mapq < 10) bc = drFade(bc, 0.4f);
                drRect((Rectangle){ bx0, y, fmaxf(1, bx1 - bx0), rowH }, bc);
                if (pxb >= 7 && rowH >= 10) {
                    char s[2] = { ev->base, 0 };
                    txtIn(F_MONO, pxb >= 9 ? 10 : 9, (Rectangle){ bx0 - 2, y, bx1 - bx0 + 4, rowH }, (Color){ 20, 20, 24, 255 }, 1, s);
                }
            } else if (ev->kind == 'D') {
                float dx0 = brX(ev->pos), dx1 = brX((double)ev->pos + ev->len);
                drRect((Rectangle){ dx0, y, dx1 - dx0, rowH }, TH->bg);
                drRect((Rectangle){ dx0, y + rowH / 2 - 0.5f, dx1 - dx0, 1 }, TH->text);
            } else if (ev->kind == 'I') {
                float ix = brX(ev->pos);
                drRect((Rectangle){ ix - 1, y - 1, 2, rowH + 2 }, TH->purple);
            }
        }
    }
    gxClipPop();
    if (maxScroll > 0) txt(F_REG, TH->fsTiny, r.x + r.width - 220, r.y + r.height - 14, TH->faint, "Shift+wheel scrolls the reads");
}

/* ---- the page ---- */
static void browserApplyNav(void) {
    StLoc *l = navCur();
    if (!memcmp(l, &g_br.applied, sizeof *l) && g_br.chrom) return;
    int c = l->chrom ? l->chrom : 7;
    double a = l->start ? (double)l->start : 117480000, b = l->end ? (double)l->end : 117680000;
    if (l->varId) g_br.selVar = l->varId;
    brSet(c, a, b + 1);
}
static void sidePanel(Rectangle r) {
    drRect(r, TH->panel);
    drVLine(r.x, r.y, r.height, TH->line);
    if (g_br.selVar > 0) {
        if (gxIconButton(gxId("close var"), (Rectangle){ r.x + r.width - 38, r.y + 10, 28, 28 }, icClose, false, "Close")) g_br.selVar = -1;
        variantPanel((Rectangle){ r.x, r.y, r.width, r.height }, g_br.selVar, false);
        return;
    }
    txt(F_SEMI, TH->fsH3, r.x + 20, r.y + 20, TH->text, "Nothing selected");
    txtPara(F_REG, TH->fsSmall, r.x + 20, r.y + 50, r.width - 40, 19, TH->dim,
            "Click a variant to see what it is, how well it was called, how common it is and what ClinVar says about it. "
            "Click a gene to select it; double-click it for its page.\n\n"
            "Drag to pan. Scroll to zoom around the pointer. Double-click to zoom in. Arrow keys pan, plus and minus zoom.\n\n"
            "Variant marks: a half-filled bar is heterozygous (one copy), a full bar homozygous (both copies). "
            "Red: high impact; amber: changes the protein; blue: other coding; slate: non-coding. "
            "A red outline means ClinVar classes it pathogenic.", 30);
}
static void pageBrowser(Rectangle r) {
    if (!g_br.locusId) g_br.locusId = gxId("browser locus");
    browserApplyNav();
    brTakeResults();
    ReadsResult *rd = readsCurrent();
    double now = (GetTime)();
    if (g_br.fetchWanted && now - g_br.changedAt > 0.12) {
        double w = brW();
        bool outside = g_br.fChrom != g_br.chrom || g_br.a < g_br.fa || g_br.b > g_br.fb || brMode() != g_br.fMode || w * 3 < (g_br.fb - g_br.fa) / 4;
        if (outside) brFetch();
        g_br.fetchWanted = false;
        if (w <= 25000 && readsAvailable()) {
            bool covered = rd && rd->chrom == g_br.chrom && rd->start <= g_br.a && rd->end >= g_br.b && (rd->end - rd->start) < w * 6;
            if (!covered) readsRequest(g_br.chrom, (long long)(g_br.a - w), (long long)(g_br.b + w));
        }
    }
    float sideW = r.width > 1300 ? 440 : 380;
    Rectangle main = { r.x, r.y, r.width - sideW, r.height };
    /* toolbar */
    float tx = main.x + 20, ty = main.y + 14;
    static int chromSel = 0;
    static const char *chromItems[25];
    if (!chromItems[0]) for (int i = 0; i < 25; i++) chromItems[i] = sstr("chr%s", g_chromNames[i + 1]), chromItems[i] = strdup(chromItems[i]);
    chromSel = g_br.chrom - 1;
    if (gxCombo(gxId("chrom combo"), (Rectangle){ tx, ty, 92, 32 }, chromItems, 25, &chromSel)) brSet(chromSel + 1, 1, (double)g_chrom[chromSel + 1].len + 1);
    tx += 102;
    Rectangle lr = { tx, ty, clampf(main.width - 700, 180, 330), 32 };
    int res = gxTextField(g_br.locusId, lr, g_br.locus, sizeof g_br.locus, "chr7:117,559,000-117,560,000 or a gene", false);
    if (res == 1) {
        int c; long long a, b;
        if (parseLocus(g_br.locus, &c, &a, &b)) { if (b - a < 40) { a -= 30; b += 30; } brSet(c, (double)a, (double)b + 1); }
        else {
            Rows *g = dbq("SELECT chrom, start, end, id FROM gene WHERE symbol=? COLLATE NOCASE ORDER BY biotype!='protein_coding' LIMIT 1", g_br.locus, NULL);
            if (g->nrow) { double pad = ((double)ri(g, 0, 2) - (double)ri(g, 0, 1)) * 0.08 + 500; brSet((int)ri(g, 0, 0), (double)ri(g, 0, 1) - pad, (double)ri(g, 0, 2) + pad); g_br.selGene = (int)ri(g, 0, 3); }
            else toast("Not a position or gene: %s", g_br.locus);
            rowsFree(&g);
        }
        g_gx.editId = 0;
    }
    tx += lr.width + 10;
    double mid = (g_br.a + g_br.b) / 2;
    if (gxIconButton(gxId("zoom out"), (Rectangle){ tx, ty, 32, 32 }, icZoomOut, false, "Zoom out (minus)")) brZoom(3, mid);
    tx += 36;
    if (gxIconButton(gxId("zoom in"), (Rectangle){ tx, ty, 32, 32 }, icZoomIn, false, "Zoom in (plus)")) brZoom(1.0 / 3, mid);
    tx += 44;
    if (gxButtonEx(gxId("whole chrom"), (Rectangle){ tx, ty, 92, 32 }, "Whole chr", BTN_GHOST, true, "Show the whole chromosome"))
        brSet(g_br.chrom, 1, (double)g_chrom[g_br.chrom].len + 1);
    tx += 98;
    /* step through your variants */
    static int stepKind = 1;
    static const char *const stepItems[] = { "Any variant", "Coding", "In ClinVar", "Rare (< 1%)" };
    gxCombo(gxId("step kind"), (Rectangle){ tx, ty, 116, 32 }, stepItems, 4, &stepKind);
    tx += 122;
    int step = 0;
    if (gxIconButton(gxId("step prev"), (Rectangle){ tx, ty, 32, 32 }, icBack, false, "Previous variant (Page Up)")) step = -1;
    if (gxIconButton(gxId("step next"), (Rectangle){ tx + 34, ty, 32, 32 }, icFwd, false, "Next variant (Page Down)")) step = 1;
    if (!g_gx.keysTaken && IsKeyPressed(KEY_PAGE_DOWN)) step = 1;
    if (!g_gx.keysTaken && IsKeyPressed(KEY_PAGE_UP)) step = -1;
    if (step) {
        static const char *const cond[] = { "1", "impact>=1", "clin>0", "(af IS NULL OR af<0.01)" };
        double cpos = (g_br.a + g_br.b) / 2;
        if (g_br.selVar > 0) {
            Rows *sv = dbq("SELECT chrom, pos FROM variant WHERE id=?", sstr("%lld", g_br.selVar), NULL);
            if (sv->nrow && ri(sv, 0, 0) == g_br.chrom) cpos = (double)ri(sv, 0, 1);
            rowsFree(&sv);
        }
        char sql[300];
        if (step > 0) snprintf(sql, sizeof sql, "SELECT id, chrom, pos FROM variant WHERE ((chrom=?1 AND pos>?2) OR chrom>?1) AND %s ORDER BY chrom, pos LIMIT 1", cond[stepKind]);
        else snprintf(sql, sizeof sql, "SELECT id, chrom, pos FROM variant WHERE ((chrom=?1 AND pos<?2) OR chrom<?1) AND %s ORDER BY chrom DESC, pos DESC LIMIT 1", cond[stepKind]);
        Rows *nx = dbq(sql, sstr("%d", g_br.chrom), sstr("%lld", (long long)cpos), NULL);
        if (nx->nrow) {
            double w = fmax(120, fmin(brW(), 2000));
            double p = (double)ri(nx, 0, 2);
            g_br.selVar = ri(nx, 0, 0);
            brSet((int)ri(nx, 0, 1), p - w / 2, p + w / 2);
        } else toast("No more variants that way");
        rowsFree(&nx);
    }
    tx += 76;
    if (g_br.selGene > 0 && gxButtonEx(gxId("open gene"), (Rectangle){ tx, ty, 96, 32 }, "Gene page", BTN_NORMAL, true, NULL)) navGene(g_br.selGene);

    /* tracks */
    float labelW = 110;
    g_brArea = (Rectangle){ main.x + labelW, main.y, main.width - labelW - 20, main.height };
    float y = main.y + 60;
    drawIdeoNav((Rectangle){ g_brArea.x, y, g_brArea.width, 36 });
    y += 54;
    drawRuler((Rectangle){ g_brArea.x, y, g_brArea.width, 24 });
    y += 30;
    float tracksTop = y;
    drawCoverage((Rectangle){ g_brArea.x, y, g_brArea.width, 64 }, rd);
    y += 72;
    drawVariants((Rectangle){ g_brArea.x, y, g_brArea.width, 40 });
    y += 48;
    Rectangle gr = { g_brArea.x, y, g_brArea.width, 0 };
    layoutGenes(gr);
    gr.height = fminf(genesH(), 12 * 32 + 16);
    gxClipPush((Rectangle){ g_brArea.x - 110, gr.y, g_brArea.width + 110, gr.height });
    drawGenes(gr);
    gxClipPop();
    y += gr.height + 6;
    double w = brW();
    if (w <= 1500 && g_br.showSeq) {
        drawSequence((Rectangle){ g_brArea.x, y, g_brArea.width, 22 });
        y += 30;
    }
    Rectangle readsR = { g_brArea.x, y, g_brArea.width, main.y + main.height - y - 12 };
    if (readsR.height > 40) drawReads(readsR, rd);
    float tracksBottom = main.y + main.height;

    /* hover guide */
    Rectangle tracks = { g_brArea.x, tracksTop, g_brArea.width, tracksBottom - tracksTop };
    bool overTracks = gxHover(tracks);
    if (overTracks) {
        double p = brPosAt(g_gx.mouse.x);
        drVLine(g_gx.mouse.x, tracksTop, tracksBottom - tracksTop, drFade(TH->text, 0.18f));
        const char *s = sstr("%s:%s", chromName(g_br.chrom), fmtInt((long long)p));
        txt(F_SEMI, TH->fsTiny, g_gx.mouse.x + 6, tracksTop - 16, TH->dim, s);
    }
    /* interaction: drag to pan, wheel to zoom, double-click to zoom in, click to select */
    GxId tid = gxId("tracks");
    bool h, held;
    bool clicked = gxBehave(tid, tracks, &h, &held);
    if (h && g_gx.pressed[0]) { g_br.dragging = true; g_br.dragMoved = false; g_br.dragX0 = g_gx.mouse.x; g_br.dragA0 = g_br.a; g_br.dragB0 = g_br.b; }
    if (g_br.dragging) {
        if (!g_gx.down[0]) g_br.dragging = false;
        else {
            float dx = g_gx.mouse.x - g_br.dragX0;
            if (fabsf(dx) > 3) g_br.dragMoved = true;
            if (g_br.dragMoved) {
                double d = (double)dx / g_brArea.width * (g_br.dragB0 - g_br.dragA0);
                brSet(g_br.chrom, g_br.dragA0 - d, g_br.dragB0 - d);
                g_gx.cursor = MOUSE_CURSOR_RESIZE_EW;
            }
        }
    }
    if (h && g_gx.wheel != 0 && !(IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT)))
        brZoom(g_gx.wheel > 0 ? 0.8 : 1.25, brPosAt(g_gx.mouse.x));
    if (h && g_gx.dbl) {
        if (g_brHoverGene >= 0 && g_br.genes) navGene((int)ri(g_br.genes, g_brHoverGene, 0));
        else brZoom(1.0 / 3, brPosAt(g_gx.mouse.x));
    }
    if (clicked && !g_br.dragMoved) {
        if (g_brHoverVar >= 0 && g_br.vars) g_br.selVar = ri(g_br.vars, (int)g_brHoverVar, 0);
        else if (g_brHoverGene >= 0 && g_br.genes) g_br.selGene = (int)ri(g_br.genes, g_brHoverGene, 0);
    }
    if (!g_gx.keysTaken) {
        if (IsKeyPressed(KEY_LEFT) || IsKeyPressedRepeat(KEY_LEFT)) brSet(g_br.chrom, g_br.a - w * 0.2, g_br.b - w * 0.2);
        if (IsKeyPressed(KEY_RIGHT) || IsKeyPressedRepeat(KEY_RIGHT)) brSet(g_br.chrom, g_br.a + w * 0.2, g_br.b + w * 0.2);
        if (IsKeyPressed(KEY_EQUAL) || IsKeyPressed(KEY_KP_ADD)) brZoom(1.0 / 2, mid);
        if (IsKeyPressed(KEY_MINUS) || IsKeyPressed(KEY_KP_SUBTRACT)) brZoom(2, mid);
    }
    sidePanel((Rectangle){ r.x + r.width - sideW, r.y, sideW, r.height });
}
