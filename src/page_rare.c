/* ========================================================================
   src/page_rare.c - rare variants, health or not.

   Three views: the remarkable-variant catalog (short sleepers, naturally
   low cholesterol, infection resistance...) checked against the genome;
   every rare protein-changing variant in plain words (what the gene does,
   what the change does to it, how much a broken copy tends to matter); and
   genes hit twice, with read-backed phase where the reads reach both.
   Tables from pipeline/16_rare.py.
   ======================================================================== */

static int g_rTab;
static GxScroll g_rScroll;
static Rows *g_rRem, *g_rVar, *g_rGene, *g_rMat, *g_rSfs, *g_rUa;
static int g_rUaOpen = -1;
static bool g_rAsked;
static int g_rVerdict;           /* 0 likely disruptive, 1 possibly, 2 all */
static bool g_rShowSuspect;
static int g_rOpen = -1;         /* expanded remarkable row */
static int g_rFound;             /* remarkable: all, found, not found */
static char g_rFocus[64];        /* a remarkable entry to scroll to (from search) */
static float g_rTop;
static void rareFocus(const char *key) {
    snprintf(g_rFocus, sizeof g_rFocus, "%s", key);
    g_rTab = 0;
    g_rFound = 0;
    navPage(PG_RARE);
}
static long long g_rStat[6];
/* 16_rare.py sets meta rare_has_af to '0' when the variant table has no
   1000 Genomes frequencies: its rarity tables are then made but left empty,
   and remarkable leaves its gene entries unchecked (status -1) */
static const char g_rNoAfKg[] = "Rarity cannot be judged without the 1000 Genomes reference. Download it from Setup on the People page, then "
                                "re-import this person.";
static const char g_rNoAfChip[] = "This needs sequencing data. A genotyping chip cannot find rare variants reliably: in a study of nearly 50,000 people, "
                                  "most very rare variants reported by consumer chips were wrong when checked by sequencing (Weedon et al. 2021, BMJ). "
                                  "The remarkable variants above are checked on the positions the chip did read.";
static bool rNoAf(void) { return !strcmp(dbMeta("rare_has_af"), "0"); }
#define g_rNoAf (dkChip() ? g_rNoAfChip : g_rNoAfKg)

static void rAsk(void) {
    g_rAsked = true;
    if (!dbHasTable("rare_var")) return;
    /* carried cards by status and evidence; the rest by category, so each
       heading of the not-found list comes once. rsids came after the table. */
    g_rRem = dbHasTable("remarkable")
        ? dbq(dbColumnExists("remarkable", "rsids")
              ? "SELECT key, cat, gene, title, status, copies, found, vid, effect, absent, rarity, source, pmid, conf, need, af, kind, rsids FROM remarkable "
                "ORDER BY status<=0, CASE WHEN status>0 THEN status END DESC, CASE WHEN status<=0 THEN cat END, "
                "CASE conf WHEN 'strong' THEN 0 WHEN 'moderate' THEN 1 ELSE 2 END, cat, title"
              : "SELECT key, cat, gene, title, status, copies, found, vid, effect, absent, rarity, source, pmid, conf, need, af, kind, '' FROM remarkable "
                "ORDER BY status<=0, CASE WHEN status>0 THEN status END DESC, CASE WHEN status<=0 THEN cat END, "
                "CASE conf WHEN 'strong' THEN 0 WHEN 'moderate' THEN 1 ELSE 2 END, cat, title", NULL) : NULL;
    g_rVar = dbq("SELECT vid, gene_id, symbol, csq, hgvs, zyg, af, am, loeuf, clin, verdict, plain, change, means, disease, inheritance, "
                 "ad_ref, ad_alt, suspect FROM rare_var ORDER BY verdict DESC, suspect, CASE WHEN af IS NULL THEN 0 ELSE af END, coalesce(am, 0) DESC", NULL);
    /* copies are counted from the gene's rare_var rows (15: changes on both
       copies, homozygous or a different change on each as 16_rare.py counts
       them, 16: on the only copy, 17: of those, possibly or likely
       disruptive), so a database from before rare_gene.nhemi, whose nhom
       also counted hemizygous calls, reads the same */
    g_rGene = dbq("SELECT gene_id, symbol, n, nhom, vids, changes, phase, phase_reads, plain, disease, inheritance, loeuf, worst, suspect, ndmg, "
                  "(SELECT count(*) FROM rare_var v WHERE v.gene_id=g.gene_id AND instr(','||g.vids||',', ','||v.vid||',') AND v.zyg IN (2,4)), "
                  "(SELECT count(*) FROM rare_var v WHERE v.gene_id=g.gene_id AND instr(','||g.vids||',', ','||v.vid||',') AND v.zyg=3), "
                  "(SELECT count(*) FROM rare_var v WHERE v.gene_id=g.gene_id AND instr(','||g.vids||',', ','||v.vid||',') AND v.zyg IN (2,3,4) AND v.verdict>=2) AS nbare "
                  "FROM rare_gene g ORDER BY suspect, (inheritance LIKE '%AR%' OR inheritance LIKE '%XLR%') AND (ndmg>=2 OR nbare>0) DESC, "
                  "(inheritance LIKE '%AR%' OR inheritance LIKE '%XLR%') DESC, ndmg DESC, worst DESC, n", NULL);
    g_rMat = dbq("SELECT grp, label, ord, fclass, n FROM rare_matrix ORDER BY ord, fclass", NULL);
    if (dbHasTable("unanswerable")) g_rUa = dbq("SELECT title, answer, explanation, sources, pmids FROM unanswerable ORDER BY ord", NULL);
    g_rSfs = dbq("SELECT lo, hi, n FROM rare_sfs ORDER BY lo", NULL);
    Rows *s = dbq("SELECT (SELECT count(*) FROM rare_var WHERE suspect=0), (SELECT count(*) FROM rare_var WHERE verdict=3 AND suspect=0), "
                  "(SELECT sum(n) FROM rare_matrix WHERE fclass=0), (SELECT sum(n) FROM rare_matrix WHERE fclass=1), "
                  "(SELECT sum(n) FROM rare_matrix WHERE fclass=2), (SELECT count(*) FROM rare_gene WHERE suspect=0)", NULL);
    for (int i = 0; i < 6 && s->nrow; i++) g_rStat[i] = ri(s, 0, i);
    rowsFree(&s);
}

static Color confColor(const char *c) {
    return !strcmp(c, "strong") ? TH->ok : !strcmp(c, "moderate") ? TH->warn : !strcmp(c, "disputed") ? TH->bad : TH->dim;
}

/* ---- remarkable variants ---- */
static void remLinks(Rows *r, int i) {
    lkBegin(rc(r, i, 3));
    if (!rnull(r, i, 7)) lkVariant(ri(r, i, 7));
    lkGene(rc(r, i, 2));
    const char *rs = rc(r, i, 17);
    for (const char *p = rs; *p && (p = strstr(p, "rs")); p += 2) {
        long long n = atoll(p + 2);
        if (n > 0 && rnull(r, i, 7)) lkUrl(sstr("rs%lld on SNPedia", n), "snpedia.com", sstr("https://www.snpedia.com/index.php/Rs%lld", n));
    }
    lkPubmed(rc(r, i, 12), "The key study");
    lkOpen();
}
static float remCard(Rectangle c, Rows *r, int i, bool draw) {
    float w = c.width - 40;
    const char *found = rc(r, i, 6);
    float h = 18 + 30 + 24 + txtParaH(F_REG, TH->fs, w, 21, rc(r, i, 8), 9) + 10 + (*found ? 22 : 0) + 22 + 24;
    if (!draw) return h;
    c.height = h;
    gxCard(c);
    int st = (int)ri(r, i, 4);
    Color tone = st == 2 ? TH->accent : TH->warn;
    drRound((Rectangle){ c.x, c.y + 14, 4, h - 28 }, 2, tone);
    float x = c.x + 20, y = c.y + 18;
    txt(F_REG, TH->fsSmall, x, y, TH->dim, sstr("%s  /  %s", rc(r, i, 1), rc(r, i, 2)));
    const char *cf = !strcmp(rc(r, i, 13), "disputed") ? "disputed: not replicated" : sstr("%s evidence", rc(r, i, 13));
    gxPill(c.x + c.width - gxPillW(cf, TH->fsTiny) - 18, y - 2, cf, confColor(rc(r, i, 13)), TH->fsTiny);
    y += 24;
    txtIn(F_SEMI, TH->fsH3, (Rectangle){ x, y, w, 26 }, TH->text, 0, rc(r, i, 3));
    y += 30;
    const char *verdict = st == 2 ? "You carry it" : "You carry one copy; the effect needs two";
    txt(F_SEMI, TH->fs, x, y, tone, verdict);
    y += 24;
    y += txtPara(F_REG, TH->fs, x, y, w, 21, TH->text, rc(r, i, 8), 9) + 10;
    if (*found) { txtIn(F_MONO, 12, (Rectangle){ x, y, w, 18 }, TH->dim, 0, found); y += 22; }
    txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ x, y + 1, w - 230, 18 }, TH->faint, 0, sstr("%s.  %s", rc(r, i, 10), rc(r, i, 11)));
    if (!rnull(r, i, 7) && gxLink(gxIdN(41000 + i), c.x + c.width - 120, y, TH->fsSmall, F_SEMI, dkReads() ? "See the reads" : "In the browser", TH->accent, NULL)) navVariant(ri(r, i, 7));
    if (*rc(r, i, 12) && gxLink(gxIdN(42000 + i), c.x + c.width - 200, y, TH->fsSmall, F_SEMI, "Paper", TH->accent, "Open the study on PubMed"))
        openUrl(sstr("https://pubmed.ncbi.nlm.nih.gov/%s/", rc(r, i, 12)));
    if (gxBoxClick(gxIdN(41500 + i), c)) remLinks(r, i);
    return h;
}

static void tabRemarkable(float x, float *y, float W, Rectangle view) {
    Rows *r = g_rRem;
    if (!r || !r->nrow) {
        *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->faint, "The remarkable-variant catalogue has not been checked yet. "
                      "The rare-variant step did not run for this person. Re-import them from the People page.", 2) + 20;
        return;
    }
    int carried = 0, half = 0;
    for (int i = 0; i < r->nrow; i++) { carried += ri(r, i, 4) == 2; half += ri(r, i, 4) == 1; }
    *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->dim,
                  "A few single variants have striking, well-documented effects: people who thrive on six hours of sleep, have naturally low cholesterol "
                  "or shrug off certain infections. Most are very rare. Each one below is checked against your genome; the evidence tag says how well the effect is established.", 4) + 14;
    txt(F_SEMI, TH->fsH3, x, *y + 6, TH->text, sstr("You carry %d of the %d%s", carried, r->nrow, half ? sstr(", and one copy of %d more", half) : ""));
    gxSegmented(gxId("rem found"), (Rectangle){ x + W - 360, *y, 360, 32 }, g_foundOpts, 3, &g_rFound);
    *y += 48;
    /* carried: cards in two columns */
    float colW = (W - 16) / 2, colY[2] = { *y, *y };
    for (int i = 0; i < r->nrow && g_rFound != 2; i++) {
        if (ri(r, i, 4) <= 0) continue;
        int col = colY[1] < colY[0] ? 1 : 0;
        Rectangle c = { x + (colW + 16) * (float)col, colY[col], colW, 0 };
        if (g_rFocus[0] && !strcmp(g_rFocus, rc(r, i, 0))) { gxScrollTo(&g_rScroll, colY[col] - g_rTop - 20); g_rFocus[0] = 0; }
        float h = remCard(c, r, i, false);
        if (colY[col] + h > view.y && colY[col] < view.y + view.height) remCard(c, r, i, true);
        colY[col] += h + 14;
    }
    *y = fmaxf(colY[0], colY[1]) + 14;
    if (g_rFound == 1) return;
    /* the rest, compact, by category */
    txt(F_SEMI, TH->fsH3, x, *y, TH->text, "Not found in your genome");
    *y += 30;
    *y += txtPara(F_REG, TH->fsSmall, x, *y, W, 19, TH->dim, "Click one to read what it does.", 1) + 8;
    const char *lastCat = "";
    for (int i = 0; i < r->nrow; i++) {
        int st = (int)ri(r, i, 4);
        if (st > 0) continue;
        if (strcmp(lastCat, rc(r, i, 1))) {
            lastCat = rc(r, i, 1);
            *y += 8;
            txt(F_SEMI, TH->fsTiny, x, *y, TH->faint, sstr("%s", lastCat));
            *y += 20;
        }
        if (g_rFocus[0] && !strcmp(g_rFocus, rc(r, i, 0))) { g_rOpen = i; gxScrollTo(&g_rScroll, *y - g_rTop - 60); g_rFocus[0] = 0; }
        bool open = g_rOpen == i;
        /* the absent text says you do not carry it: not true where nothing was
           called, nor for a gene entry, which reads rare_var, left unchecked
           for want of 1000 Genomes frequencies */
        bool unchecked = st < 0 && rNoAf() && strcmp(rc(r, i, 16), "alleles");
        const char *absent = unchecked ? (dkChip() ? "This entry could not be checked: it looks for rare changes anywhere in the gene, which needs "
                                                     "sequencing data rather than a chip."
                                                   : "This entry could not be checked: it looks for rare changes in the gene, and without 1000 Genomes "
                                                     "frequencies rarity cannot be judged.")
                           : st < 0 ? "Your genome has no confident call at the positions checked, so whether you carry this variant is unknown." : rc(r, i, 9);
        float h = 40;
        if (open) h += txtParaH(F_REG, TH->fsSmall, W - 48, 19, rc(r, i, 8), 9) + txtParaH(F_REG, TH->fsSmall, W - 48, 19, absent, 3) + 46;
        Rectangle row = { x, *y, W, h };
        bool hov;
        if (gxBehave(gxIdN(43000 + i), (Rectangle){ x, *y, W, 40 }, &hov, NULL)) g_rOpen = open ? -1 : i;
        if (hov || open) { drRound(row, 6, open ? TH->panel2 : TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
        DrawCircleV((Vector2){ x + 14, *y + 20 }, 3.5f, st < 0 ? TH->faint : TH->line);
        txtIn(F_SEMI, TH->fs, (Rectangle){ x + 28, *y + 10, W * 0.45f, 20 }, TH->text, 0, rc(r, i, 3));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ x + 28 + W * 0.45f, *y + 11, W * 0.3f, 20 }, TH->dim, 0, rc(r, i, 10));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ x + W - 200, *y + 11, 186, 20 }, st < 0 ? TH->warn : TH->faint, 2,
              unchecked ? "not checked: no frequencies" : st < 0 ? "not called: no data here" : "not carried");
        if (open) {
            float ty = *y + 44;
            ty += txtPara(F_REG, TH->fsSmall, x + 28, ty, W - 48, 19, TH->text, rc(r, i, 8), 9) + 6;
            ty += txtPara(F_REG, TH->fsSmall, x + 28, ty, W - 48, 19, TH->dim, absent, 3) + 6;
            txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ x + 28, ty, W - 260, 18 }, TH->faint, 0, sstr("%s gene.  %s", rc(r, i, 2), rc(r, i, 11)));
            if (*rc(r, i, 12) && gxLink(gxIdN(44000 + i), x + W - 80, ty, TH->fsSmall, F_SEMI, "Paper", TH->accent, "Open the study on PubMed"))
                openUrl(sstr("https://pubmed.ncbi.nlm.nih.gov/%s/", rc(r, i, 12)));
            if (gxLink(gxIdN(44500 + i), x + W - 170, ty, TH->fsSmall, F_SEMI, "Links", TH->accent, "Gene and variant pages on the web")) remLinks(r, i);
        }
        *y += h + 2;
    }
}

/* questions genes cannot answer (yet): popular ones with no usable test */
static void rUnanswerable(float x, float *y, float W) {
    Rows *u = g_rUa;
    if (!u || !u->nrow) return;
    *y += 26;
    txt(F_SEMI, TH->fsH2, x, *y, TH->text, "What genes cannot tell you (yet)");
    *y += 34;
    *y += txtPara(F_REG, TH->fsSmall, x, *y, W, 19, TH->dim,
                  "Popular questions with no genetic test behind them: either no variant has been found, or the effects are so spread out and small "
                  "that they say nothing about one person. Click one for what research has found.", 3) + 10;
    for (int i = 0; i < u->nrow; i++) {
        bool open = g_rUaOpen == i;
        float h = 46;
        if (open) h += txtParaH(F_REG, TH->fsSmall, W - 48, 19, rc(u, i, 2), 8) + 40;
        Rectangle row = { x, *y, W, h };
        bool hov;
        if (gxBehave(gxIdN(48000 + i), (Rectangle){ x, *y, W, 46 }, &hov, NULL)) g_rUaOpen = open ? -1 : i;
        if (hov || open) { drRound(row, 6, open ? TH->panel2 : TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
        txtIn(F_SEMI, TH->fs, (Rectangle){ x + 16, *y + 6, W * 0.42f, 20 }, TH->text, 0, rc(u, i, 0));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ x + 16, *y + 26, W - 32, 18 }, TH->dim, 0, rc(u, i, 1));
        if (open) {
            float ty = *y + 52;
            ty += txtPara(F_REG, TH->fsSmall, x + 16, ty, W - 48, 19, TH->text, rc(u, i, 2), 8) + 8;
            txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ x + 16, ty, W - 200, 18 }, TH->faint, 0, rc(u, i, 3));
            char first[24] = "";
            sscanf(rc(u, i, 4), "%23[0-9]", first);
            if (first[0] && gxLink(gxIdN(49000 + i), x + W - 170, ty, TH->fsSmall, F_SEMI, "Main paper", TH->accent, "Open it on PubMed"))
                openUrl(sstr("https://pubmed.ncbi.nlm.nih.gov/%s/", first));
        }
        *y += h + 2;
    }
}

/* ---- your rare changes ---- */
static void rStat(Rectangle c, const char *label, const char *value, const char *sub, Color tone) {
    gxCard(c);
    drRound((Rectangle){ c.x + 16, c.y + 18, 3, 16 }, 1.5f, tone);
    txt(F_REG, TH->fsSmall, c.x + 26, c.y + 18, TH->dim, label);
    txt(F_SEMI, 26, c.x + 18, c.y + 44, TH->text, value);
    txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 18, c.y + 82, c.width - 30, 16 }, TH->faint, 0, sub);
}

static void rSpectrum(Rectangle c) {
    gxCard(c);
    gxCardTitle(c, "How common your variants are", "Every variant you carry, by how often it turns up among 2,590 people in 1000 Genomes");
    if (!g_rSfs || !g_rSfs->nrow) return;
    int n = g_rSfs->nrow;
    double novel = 0, mx = 1;
    for (int i = 0; i < n; i++) mx = fmax(mx, rf(g_rSfs, i, 2));
    for (int i = 0; g_rMat && i < g_rMat->nrow; i++) if (ri(g_rMat, i, 3) == 0) novel += rf(g_rMat, i, 4);
    mx = fmax(mx, novel);
    Rectangle a = { c.x + 22, c.y + 76, c.width - 44, c.height - 130 };
    float slot = a.width / (float)(n + 2);
    float lmx = (float)log10(mx + 1);
    for (int i = -1; i < n; i++) {
        double v = i < 0 ? novel : rf(g_rSfs, i, 2);
        float hh = a.height * (float)log10(v + 1) / lmx;
        float bx = i < 0 ? a.x : a.x + slot * (float)(i + 2);
        Color col = i < 0 ? TH->purple : rf(g_rSfs, i, 1) <= 0.0011 ? TH->bad : rf(g_rSfs, i, 1) <= 0.0101 ? TH->warn : rf(g_rSfs, i, 1) <= 0.051 ? TH->accent2 : TH->accent;
        Rectangle b = { bx + 2, a.y + a.height - hh, slot - 4, hh };
        drRect(b, drFade(col, 0.85f));
        if (gxHover((Rectangle){ bx, a.y, slot, a.height }))
            gxTip(gxIdN(45000 + i), b, i < 0 ? sstr("%s confident heterozygous variants not in the panel", fmtInt((long long)v))
                                              : sstr("%s variants seen in %s to %s of chromosomes", fmtInt((long long)v), fmtAf(rf(g_rSfs, i, 0)), fmtAf(rf(g_rSfs, i, 1))));
    }
    drHLine(a.x, a.y + a.height, a.width, TH->line);
    float ly = a.y + a.height + 6;
    txtIn(F_REG, TH->fsTiny, (Rectangle){ a.x - 6, ly, slot + 12, 16 }, TH->dim, 1, "not seen");
    static const struct { double f; const char *l; } ticks[] = { { 0.001, "1 in 1,000" }, { 0.01, "1 in 100" }, { 0.1, "1 in 10" }, { 1.0, "everyone" } };
    for (size_t k = 0; k < NELEMS(ticks); k++) {
        float t = (float)((log10(ticks[k].f) + 3.5) / 3.5);
        float tx = a.x + slot * 2 + (a.width - slot * 2) * clamp01(t);
        drVLine(tx, a.y + a.height, 4, TH->line);
        txtIn(F_REG, TH->fsTiny, (Rectangle){ tx - 50, ly, 100, 16 }, TH->faint, 1, ticks[k].l);
    }
    txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 22, c.y + c.height - 26, c.width - 44, 16 }, TH->faint, 0,
          "Heights on a log scale. Red: under 1 in 1,000 chromosomes; amber: under 1 in 100; the rest are shared with many people.");
}

/* the matrix's five value columns, for a card this wide */
static float rMatrixColW(float width) { return (width - 40 - 190) / 5; }
/* do its header and widest count fit those columns without being cut? */
static bool rMatrixFits(float width) {
    int need = txtW(F_SEMI, TH->fsTiny + 1, "Common");
    for (int i = 0; g_rMat && i < g_rMat->nrow; i++) {
        int n = txtW(F_REG, TH->fsSmall, fmtInt(ri(g_rMat, i, 4)));
        if (n > need) need = n;
    }
    return rMatrixColW(width) - 8 >= (float)need;
}

static void rMatrix(Rectangle c) {
    gxCard(c);
    gxCardTitle(c, "By effect and rarity", "Every call the variant caller passed, before the artefact checks used below");
    static const char *const cols[] = { "Unseen", "< 0.1%", "0.1 - 1%", "1 - 5%", "Common" };
    float x0 = c.x + 20, lw = 190, cw = rMatrixColW(c.width), y = c.y + 70;
    txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ x0, y, lw - 12, 18 }, TH->faint, 2, "Share of chromosomes:");
    for (int k = 0; k < 5; k++) txtIn(F_SEMI, TH->fsTiny + 1, (Rectangle){ x0 + lw + cw * (float)k, y, cw - 8, 18 }, TH->dim, 2, cols[k]);
    y += 24;
    drHLine(x0, y - 4, c.width - 40, TH->line);
    if (!g_rMat) return;
    int lastOrd = -1;
    for (int i = 0; i < g_rMat->nrow; i++) {
        int ord = (int)ri(g_rMat, i, 2);
        if (ord == lastOrd) continue;
        lastOrd = ord;
        Color lc = ord == 0 ? TH->bad : ord <= 2 ? TH->warn : TH->text;
        txtIn(F_REG, TH->fsSmall, (Rectangle){ x0, y, lw, 20 }, lc, 0, rc(g_rMat, i, 1));
        for (int j = i; j < g_rMat->nrow && ri(g_rMat, j, 2) == ord; j++) {
            int fc = (int)ri(g_rMat, j, 3);
            txtIn(F_REG, TH->fsSmall, (Rectangle){ x0 + lw + cw * (float)fc, y, cw - 8, 20 }, fc <= 1 && ord <= 1 ? TH->text : TH->dim, 2, fmtInt(ri(g_rMat, j, 4)));
        }
        y += 23;
    }
    const char *hom = dbMeta("rare_unmatched_hom");
    txtPara(F_REG, TH->fsTiny + 1, x0, c.y + c.height - 62, c.width - 40, 16, TH->faint,
            sstr("\"Unseen\": not in 1000 Genomes; only clean heterozygous calls count (depth 15 to 70, a quarter to three quarters of reads showing the variant). "
                 "Another %s homozygous calls are missing from the panel too, but those are common variants at sites the panel filtered out.", fmtInt(atoll(hom ? hom : "0"))), 3);
}

static const char *const g_rVerdicts[] = { "Likely disruptive", "Possibly disruptive", "All rare changes" };

static float rVarRow(Rectangle row, Rows *r, int i, bool draw) {
    float w = row.width - 44;
    const char *plainTxt = rc(r, i, 11), *chg = rc(r, i, 12), *means = rc(r, i, 13), *dis = rc(r, i, 14);
    float h = 12 + 24 + (*plainTxt ? 20 : 0) + txtParaH(F_REG, TH->fsSmall, w, 19, chg, 3) + (*means ? txtParaH(F_REG, TH->fsSmall - 1, w, 18, means, 3) : 0)
              + (*dis ? 18 : 0) + 14;
    if (!draw) return h;
    row.height = h;
    bool hov;
    if (gxBehave(gxIdN(46000 + i), row, &hov, NULL)) lkVariantMenu(ri(r, i, 0));
    drRound(row, 8, hov ? TH->hover : TH->card);
    if (hov) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
    int verdict = (int)ri(r, i, 10);
    Color tone = ri(r, i, 18) ? TH->faint : verdict == 3 ? TH->bad : verdict == 2 ? TH->warn : TH->dim;
    drRound((Rectangle){ row.x, row.y + 10, 3, h - 20 }, 1.5f, tone);
    float x = row.x + 20, y = row.y + 12;
    const char *sym = rc(r, i, 2);
    txt(F_SEMI, TH->fs + 1, x, y, TH->text, sym);
    float gx = x + (float)txtW(F_SEMI, TH->fs + 1, sym) + 10;
    txtIn(F_REG, TH->fs, (Rectangle){ gx, y + 1, 260, 20 }, TH->dim, 0, *rc(r, i, 4) ? rc(r, i, 4) : csqLabel(rc(r, i, 3)));
    float px = row.x + row.width - 16;
    /* af is a share of chromosomes (alleles), not of people who carry it */
    double af = rf(r, i, 6);
    const char *fr = rnull(r, i, 6) ? "not in 1000G" : af <= 0 ? "very rare in 1000G"
                   : af < 0.01 ? sstr("%s chromosomes", fmtAf(af)) : sstr("%s of chromosomes", fmtAf(af));
    px -= gxPillW(fr, TH->fsTiny); gxPill(px, y, fr, rnull(r, i, 6) ? TH->purple : TH->accent2, TH->fsTiny); px -= 6;
    const char *zs = zygName((int)ri(r, i, 5));
    px -= gxPillW(zs, TH->fsTiny); gxPill(px, y, zs, TH->purple, TH->fsTiny); px -= 6;
    const char *cs = csqLabel(rc(r, i, 3));
    px -= gxPillW(cs, TH->fsTiny); gxPill(px, y, cs, tone, TH->fsTiny); px -= 6;
    int clin = (int)ri(r, i, 9);
    if (clin) { const char *cn = sstr("ClinVar: %s", clinName(clin)); px -= gxPillW(cn, TH->fsTiny); gxPill(px, y, cn, clinColor(clin), TH->fsTiny); }
    y += 26;
    if (*plainTxt) { txtIn(F_REG, TH->fsSmall, (Rectangle){ x, y, w, 18 }, TH->text, 0, plainTxt); y += 20; }
    y += txtPara(F_REG, TH->fsSmall, x, y, w, 19, TH->dim, chg, 3);
    if (*means) y += txtPara(F_REG, TH->fsSmall - 1, x, y, w, 18, ri(r, i, 18) ? TH->warn : TH->faint, means, 3);
    if (*dis) txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ x, y, w, 18 }, TH->faint, 0, sstr("Conditions linked to %s: %s%s", sym, dis, *rc(r, i, 15) ? sstr(" (%s)", rc(r, i, 15)) : ""));
    return h;
}

/* row heights of the rare-changes list: wrapping every row's text each
   frame stalls the page with all rare changes shown, so each row is
   measured once for a given width and font size (g_rVar is loaded once) */
static float *g_rVarH;
static float g_rVarHW, g_rVarHScale;
static int g_rVarHPx;
static float rVarRowH(Rectangle row, Rows *r, int i) {
    if (!g_rVarH || g_rVarHW != row.width || g_rVarHScale != g_fontScale || g_rVarHPx != TH->fsSmall) {
        free(g_rVarH);
        g_rVarH = xalloc(sizeof(float) * (size_t)r->nrow);
        for (int k = 0; k < r->nrow; k++) g_rVarH[k] = -1;
        g_rVarHW = row.width, g_rVarHScale = g_fontScale, g_rVarHPx = TH->fsSmall;
    }
    if (g_rVarH[i] < 0) g_rVarH[i] = rVarRow(row, r, i, false);
    return g_rVarH[i];
}

static void tabRareChanges(float x, float *y, float W, Rectangle view) {
    if (rNoAf()) { *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->faint, g_rNoAf, 3) + 20; return; }
    float cw = (W - 3 * 14) / 4;
    rStat((Rectangle){ x, *y, cw, 106 }, "Rare protein changes", fmtInt(g_rStat[0]), "under 1% or not in 1000G, quality-checked", TH->accent);
    rStat((Rectangle){ x + (cw + 14), *y, cw, 106 }, "Likely disruptive", fmtInt(g_rStat[1]), "truncating in a sensitive gene, or AlphaMissense high", TH->bad);
    rStat((Rectangle){ x + 2 * (cw + 14), *y, cw, 106 }, "Not in 1000 Genomes", fmtInt(g_rStat[2]), "confident calls not in 1000 Genomes, nearly all non-coding", TH->purple);
    rStat((Rectangle){ x + 3 * (cw + 14), *y, cw, 106 }, "Under 1 in 1,000 chromosomes", fmtInt(g_rStat[3]), sstr("and %s more under 1 in 100", fmtInt(g_rStat[4])), TH->warn);
    *y += 106 + 16;
    float lw = W * 0.56f;
    bool wide = W > 1000 && rMatrixFits(W - lw - 8);
    if (wide) {
        rSpectrum((Rectangle){ x, *y, lw - 8, 350 });
        rMatrix((Rectangle){ x + lw + 8, *y, W - lw - 8, 350 });
        *y += 350 + 22;
    } else {
        rSpectrum((Rectangle){ x, *y, W, 300 }); *y += 314;
        rMatrix((Rectangle){ x, *y, W, 330 }); *y += 344;
    }
    txt(F_SEMI, TH->fsH2, x, *y, TH->text, "Your rare protein changes, in plain words");
    *y += 36;
    *y += txtPara(F_REG, TH->fsSmall, x, *y, W, 19, TH->dim,
                  "Each rare change to a protein you carry: what the gene does, what the change does to it, and how much a broken copy tends to matter. "
                  "Everyone carries a hundred or so broken genes and hundreds of rare changes; nearly all are harmless. Click one for the full record and the reads.", 3) + 12;
    gxSegmented(gxId("rare verdict"), (Rectangle){ x, *y, 480, 32 }, g_rVerdicts, 3, &g_rVerdict);
    gxCheckbox(gxId("rare suspect"), (Rectangle){ x + 500, *y + 6, 360, 22 }, &g_rShowSuspect, "Include calls in hard-to-map genes");
    *y += 46;
    Rows *r = g_rVar;
    int shown = 0;
    for (int i = 0; r && i < r->nrow; i++) {
        int v = (int)ri(r, i, 10);
        if (!g_rShowSuspect && ri(r, i, 18)) continue;
        if ((g_rVerdict == 0 && v < 3) || (g_rVerdict == 1 && v < 2)) continue;
        Rectangle row = { x, *y, W, 0 };
        float h = rVarRowH(row, r, i);
        if (*y + h > view.y && *y < view.y + view.height) rVarRow(row, r, i, true);
        *y += h + 8;
        shown++;
    }
    if (!shown) { txt(F_REG, TH->fs, x, *y, TH->faint, "Nothing in this group."); *y += 30; }
}

/* ---- two changes in one gene ---- */
/* rare_gene.phase in plain words. Older databases left it '' for every
   unphased pair, so there the read count tells "disagree" from "too few". */
static const char *rPhaseText(const char *ph, long long reads, long long nboth, long long nonly) {
    if (nboth || !strcmp(ph, "hom")) return "Both copies changed";
    if (nonly || !strcmp(ph, "hemi")) return "Only copy changed (hemizygous)";
    if (!strcmp(ph, "cis")) return sstr("Same copy (cis), %lld reads", reads);
    if (!strcmp(ph, "trans")) return sstr("Different copies (trans), %lld reads", reads);
    if (!strcmp(ph, "far")) return "Phase unknown: too far apart for one read";
    if (!strcmp(ph, "notsnv")) return "Phase not checked: not a single-letter change";
    if (!strcmp(ph, "noreads")) return "Phase not checked: no reads on file";
    if (!strcmp(ph, "multi")) return "Phase not checked: more than two changes";
    if (!strcmp(ph, "mixed") || (!*ph && reads >= 4)) return sstr("Phase unclear: %lld reads disagree", reads);
    if (!strcmp(ph, "few") || (!*ph && reads > 0))
        return reads == 0 ? "Phase unknown: no read spans both" : reads == 1 ? "Phase unknown: only 1 read spans both"
                                                                             : sstr("Phase unknown: only %lld reads span both", reads);
    return "Phase unknown";
}

static void tabTwoHits(float x, float *y, float W, Rectangle view) {
    if (rNoAf()) { *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->faint, g_rNoAf, 3) + 20; return; }
    *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->dim,
                  "For a recessive condition, what matters is whether both copies of a gene are changed. Two different rare changes in one gene could sit on the same copy "
                  "(cis: the other copy still works) or on different copies (trans: neither copy is untouched). Where single reads span both changes, the reads settle it.", 4) + 14;
    Rows *r = g_rGene;
    if (!r || !r->nrow) { txt(F_REG, TH->fs, x, *y, TH->faint, "No gene has two rare protein changes."); *y += 30; return; }
    gxCheckbox(gxId("two suspect"), (Rectangle){ x, *y, 360, 22 }, &g_rShowSuspect, "Include hard-to-map genes");
    *y += 36;
    for (int i = 0; i < r->nrow; i++) {
        bool sus = ri(r, i, 13);
        if (sus && !g_rShowSuspect) continue;
        const char *ph = rc(r, i, 6), *inh = rc(r, i, 10);
        /* recessive: both copies must be hit (X-linked recessive too, for two X copies);
           dominant-like: one changed copy can be enough */
        bool ar = strstr(inh, "AR") || strstr(inh, "XLR");
        bool dom = !ar && (strstr(inh, "AD") || strstr(inh, "XL") || strstr(inh, "YL"));
        bool hasDis = *rc(r, i, 9);
        long long nh = ri(r, i, 15), nonly = ri(r, i, 16), bare = ri(r, i, 17);   /* both copies, only copy, disruptive of those */
        const char *phase = rPhaseText(ph, ri(r, i, 7), nh, nonly);
        Color pc = nh || nonly || !strcmp(ph, "trans") ? TH->warn : !strcmp(ph, "cis") ? TH->ok : TH->dim;
        long long dmg = ri(r, i, 14);
        bool quiet = !nonly && !ar && !dom && !hasDis;
        const char *meaning = nonly ? (bare ? "You have one copy of this gene (it lies on the X or Y chromosome) and a change that may harm the protein is on it, "
                                              "so no untouched copy is left: worth a closer look (and the reads)."
                                            : "You have one copy of this gene (it lies on the X or Y chromosome), so a change has no second copy to cover for it, "
                                              "but no change here looks damaging.") :
                              dom ? "Its recorded condition can show with one changed copy, so the two-copy reasoning above does not apply: "
                                    "see each change under Your rare changes." :
                              !ar && hasDis ? "A condition is recorded for this gene, but not as recessive, so the pairing says little: see each change under Your rare changes." :
                              quiet ? "No condition is recorded for this gene, so these changes are unlikely to matter much." :
                              bare ? "The gene has a recessive condition and changes that may harm the protein leave neither copy untouched: "
                                     "worth a closer look (and the reads)." :
                              (!strcmp(ph, "cis") && !nh) ? "Both changes are on one copy, so the other copy is intact: like being a carrier." :
                              dmg >= 2 ? (!strcmp(ph, "trans") ? "The gene has a recessive condition, both changes may harm the protein and the reads put them on different copies: "
                                                                 "worth a closer look (and the reads)."
                                                               : "The gene has a recessive condition, both changes may harm the protein and they could be on different copies: "
                                                                 "worth a closer look (and the reads).") :
                              dmg == 1 ? "The gene has a recessive condition, but only one of the changes looks damaging, so at least one copy probably still works." :
                              "The gene has a recessive condition, but no change here looks damaging.";
        if (!sus && ((ar && (bare || (dmg >= 2 && strcmp(ph, "cis")))) || (nonly && bare && hasDis))) pc = TH->bad;
        float h = 16 + 26 + (*rc(r, i, 8) ? 20 : 0) + 20 + txtParaH(F_REG, TH->fsSmall, W - 40, 19, meaning, 2) + (*rc(r, i, 9) ? 20 : 0) + 6;
        Rectangle c = { x, *y, W, h };
        if (*y + h > view.y && *y < view.y + view.height) {
            gxCard(c);
            if (!sus && ((ar && dmg >= 1) || bare)) drRound((Rectangle){ c.x, c.y + 14, 4, h - 28 }, 2, pc);
            float tx = c.x + 20, ty = c.y + 16;
            if (gxLink(gxIdN(47000 + i), tx, ty, TH->fsH3, F_SEMI, rc(r, i, 1), TH->accent, "Open the gene page")) navGene((int)ri(r, i, 0));
            gxPill(c.x + c.width - gxPillW(phase, TH->fsTiny) - 18, ty, phase, pc, TH->fsTiny);
            if (sus) txt(F_REG, TH->fsSmall, tx + (float)txtW(F_SEMI, TH->fsH3, rc(r, i, 1)) + 12, ty + 3, TH->warn, "hard-to-map gene: likely artefacts");
            ty += 26;
            if (*rc(r, i, 8)) { txtIn(F_REG, TH->fsSmall, (Rectangle){ tx, ty, W - 40, 18 }, TH->text, 0, rc(r, i, 8)); ty += 20; }
            txtIn(F_REG, TH->fsSmall, (Rectangle){ tx, ty, W - 40, 18 }, TH->dim, 0, rc(r, i, 5));
            ty += 20;
            ty += txtPara(F_REG, TH->fsSmall, tx, ty, W - 40, 19, quiet ? TH->dim : TH->text, meaning, 2);
            if (*rc(r, i, 9)) txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ tx, ty, W - 40, 18 }, TH->faint, 0, sstr("%s (%s)", rc(r, i, 9), inh));
            if (gxBoxClick(gxIdN(47500 + i), c)) lkGeneMenu(rc(r, i, 1));
        }
        *y += h + 10;
    }
}

static const char *const g_rTabs[] = { "Remarkable variants", "Your rare changes", "Two changes in one gene" };
static void pageRare(Rectangle r) {
    if (!g_rAsked) rAsk();
    if (!g_rVar) {
        pageMissing(r, "Rare and remarkable", rNoAf() ? g_rNoAf
                    : "The rare-variant step did not run for this person. Re-import them from the People page.");
        return;
    }
    float x = r.x + 32, W = fminf(1240, r.width - 64);
    txt(F_SEMI, TH->fsH1, x, r.y + 22, TH->text, "Rare and remarkable");
    if (gxTabs(gxId("rare tabs"), x - 14, r.y + 58, 40, g_rTabs, 3, &g_rTab)) g_rScroll = (GxScroll){ 0 };
    drHLine(r.x, r.y + 98, r.width, TH->line);
    Rectangle body = { r.x, r.y + 99, r.width, r.height - 99 };
    float top = gxScrollBegin(&g_rScroll, body);
    g_rTop = top;
    float y = top + 22;
    switch (g_rTab) {
    case 0: tabRemarkable(x, &y, W, body); rUnanswerable(x, &y, W); break;
    case 1: tabRareChanges(x, &y, W, body); break;
    case 2: tabTwoHits(x, &y, W, body); break;
    }
    y += 20;
    txtPara(F_REG, TH->fsTiny + 1, x, y, W, 17, TH->faint,
            sstr("Rarity is judged against 1000 Genomes (2,590 people); gnomAD, with far more people, can be checked from any variant's panel. "
                 "Predictions such as AlphaMissense are guides, not verdicts, and single %s calls can be wrong.", dkChip() ? "chip" : "sequencing"), 3);
    y += 60;
    gxScrollEnd(&g_rScroll, body, y - top);
}
