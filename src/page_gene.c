/* ========================================================================
   src/page_gene.c - genes: an index of the ones worth a look, and a page
   for each gene (what it does, how intolerant it is of damage, the
   conditions linked to it, its transcript with the sample's variants on it,
   and the variants themselves).
   ======================================================================== */

typedef struct {
    int id;
    Rows *g, *dis, *vars, *exons, *acmg, *pgx;
    GxTable vt;
    long long sel;
    GxScroll scroll;
    int vfilter;
} GenePage;
static GenePage g_gp = { .id = -1, .vt = { .sel = -1, .sortCol = -1 }, .sel = -1 };
static Rows *g_giRows;
static GxTable g_giT = { .sel = -1, .sortCol = -1 };
static char g_giText[80];
static int g_giScope;
static bool g_giDirty = true;

static void gpLoad(int id) {
    if (g_gp.id == id) return;
    rowsFree(&g_gp.g); rowsFree(&g_gp.dis); rowsFree(&g_gp.vars); rowsFree(&g_gp.exons); rowsFree(&g_gp.acmg); rowsFree(&g_gp.pgx);
    g_gp.id = id;
    g_gp.sel = -1;
    g_gp.vt = (GxTable){ .sel = -1, .sortCol = -1 };
    g_gp.scroll = (GxScroll){ 0 };
    const char *ids = sstr("%d", id);
    g_gp.g = dbq("SELECT id, ensg, symbol, name, biotype, chrom, start, end, strand, tx, mane, cds_start, cds_end, nexon, hgnc, "
                 "protein, function, uniprot_disease, uniprot, pli, loeuf, mis_z, nvar, ncoding FROM gene WHERE id=?", ids, NULL);
    g_gp.dis = dbq("SELECT disease, inheritance, disease_id FROM gene_disease WHERE gene=? ORDER BY disease", ids, NULL);
    g_gp.exons = dbq("SELECT start, end, kind FROM exon WHERE gene=? ORDER BY start", ids, NULL);
    if (g_gp.g->nrow) {
        g_gp.acmg = dbHasTable("acmg") ? dbq("SELECT category, disease, inheritance, report FROM acmg WHERE gene=?", rc(g_gp.g, 0, 2), NULL) : NULL;
        g_gp.pgx = dbHasTable("pgx_gene") ? dbq("SELECT diplotype, phenotype FROM pgx_gene WHERE gene=?", rc(g_gp.g, 0, 2), NULL) : NULL;
    }
    g_gp.vars = dbq("SELECT id, pos, ref, alt, hgvs_p, csq, impact, zyg, af, clin, am, dp, rsid, pass FROM variant WHERE gene_id=? "
                    "ORDER BY impact DESC, clin IN (5,6) DESC, pos", ids, NULL);
}

/* the canonical transcript, exons to scale against the gene's span, with
   the sample's variants as lollipops */
static void geneDiagram(Rectangle r) {
    Rows *g = g_gp.g;
    long long gs = ri(g, 0, 6), ge = ri(g, 0, 7);
    if (ge <= gs) return;
    float mid = r.y + r.height * 0.62f;
    double span = (double)(ge - gs + 1);
    #define GX(p) (r.x + (float)((double)((p) - gs) / span * r.width))
    drRect((Rectangle){ r.x, mid - 0.75f, r.width, 1.5f }, TH->lineHi);
    if (g_gp.exons) for (int i = 0; i < g_gp.exons->nrow; i++) {
        float x0 = GX(ri(g_gp.exons, i, 0)), x1 = GX(ri(g_gp.exons, i, 1) + 1);
        bool cds = ri(g_gp.exons, i, 2) == 1;
        float h = cds ? 16 : 9;
        drRect((Rectangle){ x0, mid - h / 2, fmaxf(1.5f, x1 - x0), h }, cds ? TH->accent : drFade(TH->accent, 0.55f));
    }
    int strand = (int)ri(g, 0, 8);
    txt(F_REG, TH->fsTiny, r.x, mid + 14, TH->faint, strand > 0 ? "5'" : "3'");
    txt(F_REG, TH->fsTiny, r.x + r.width - 12, mid + 14, TH->faint, strand > 0 ? "3'" : "5'");
    /* lollipops: coding and splice variants, plus anything in ClinVar */
    Rows *v = g_gp.vars;
    int hover = -1;
    for (int i = 0; v && i < v->nrow; i++) {
        int imp = (int)ri(v, i, 6), clin = (int)ri(v, i, 9);
        if (imp < 1 && !(clin >= 5 && clin <= 7)) continue;
        float x = GX(ri(v, i, 1));
        float h = imp >= 3 ? 56 : imp == 2 ? 44 : 30;
        Color c = (clin == 5 || clin == 6) ? TH->bad : impactColor(imp);
        drRect((Rectangle){ x - 0.5f, mid - h, 1, h - 8 }, drFade(c, 0.7f));
        DrawCircleV((Vector2){ x, mid - h }, ri(v, i, 7) == 2 ? 5.5f : 4.0f, c);
        if (ri(v, i, 7) != 2) DrawCircleV((Vector2){ x, mid - h }, 2.0f, TH->card);
        if (Vector2Distance(g_gx.mouse, (Vector2){ x, mid - h }) < 7) hover = i;
    }
    if (hover >= 0) {
        gxTip(gxId("gene lolli") ^ (GxId)ri(v, hover, 0), (Rectangle){ g_gx.mouse.x - 8, g_gx.mouse.y - 8, 16, 16 },
              sstr("%s  %s\n%s, %s\nClick for details", *rc(v, hover, 4) ? rc(v, hover, 4) : sstr("%s>%s", rc(v, hover, 2), rc(v, hover, 3)),
                   csqLabel(rc(v, hover, 5)), zygName((int)ri(v, hover, 7)), fmtAf(rf(v, hover, 8))));
        g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
        if (g_gx.pressed[0]) g_gp.sel = ri(v, hover, 0);
    }
    #undef GX
}

static void gvCell(int row, int col, Rectangle c, bool sel, void *u) {
    (void)sel; (void)u;
    Rows *v = g_gp.vars;
    int fs = TH->fsSmall;
    switch (col) {
    case 0: txtIn(F_MONO, 12, c, TH->dim, 0, fmtInt(ri(v, row, 1))); break;
    case 1: txtIn(F_REG, fs, c, TH->text, 0, *rc(v, row, 4) ? rc(v, row, 4) : sstr("%.20s>%.20s", rc(v, row, 2), rc(v, row, 3))); break;
    case 2: DrawCircleV((Vector2){ c.x + 4, c.y + c.height / 2 }, 3.5f, impactColor((int)ri(v, row, 6)));
            txtIn(F_REG, fs, (Rectangle){ c.x + 14, c.y, c.width - 14, c.height }, TH->text, 0, csqLabel(rc(v, row, 5))); break;
    case 3: txtIn(F_REG, fs, c, TH->dim, 0, zygShort((int)ri(v, row, 7))); break;
    case 4: txtIn(F_REG, fs, c, rnull(v, row, 8) ? TH->warn : TH->text, 2, fmtAf(rf(v, row, 8))); break;
    case 5: { int cl = (int)ri(v, row, 9); if (cl) txtIn(F_REG, fs, c, clinColor(cl), 0, clinName(cl)); break; }
    }
}
static const GxCol g_gvCols[] = { { "POSITION", 110, 0, false, NULL }, { "CHANGE", -2, 0, false, NULL }, { "EFFECT", 140, 0, false, NULL },
                                  { "COPIES", 60, 0, false, NULL }, { "FREQUENCY", 90, 2, false, NULL }, { "CLINVAR", 120, 0, false, NULL } };

static void geneDetail(Rectangle r) {
    Rows *g = g_gp.g;
    if (!g || !g->nrow) { pageMissing(r, "Gene not found", "That gene is not in the database."); return; }
    float sideW = g_gp.sel > 0 ? 440 : 0;
    Rectangle main = { r.x, r.y, r.width - sideW, r.height };
    float top = gxScrollBegin(&g_gp.scroll, main);
    float x = main.x + 32, W = main.width - 64, y = top + 24;
    const char *sym = rc(g, 0, 2);
    txt(F_SEMI, 30, x, y, TH->text, sym);
    float px = x + (float)txtW(F_SEMI, 30, sym) + 16;
    px += gxPill(px, y + 8, !strcmp(rc(g, 0, 4), "protein_coding") ? "protein coding" : rc(g, 0, 4), TH->accent2, TH->fsTiny) + 6;
    if (g_gp.acmg && g_gp.acmg->nrow) px += gxPill(px, y + 8, "ACMG secondary-findings gene", TH->warn, TH->fsTiny) + 6;
    if (g_gp.pgx && g_gp.pgx->nrow) px += gxPill(px, y + 8, "pharmacogene", TH->purple, TH->fsTiny) + 6;
    y += 42;
    txt(F_REG, TH->fsH3, x, y, TH->dim, *rc(g, 0, 3) ? rc(g, 0, 3) : rc(g, 0, 15));
    y += 28;
    int chrom = (int)ri(g, 0, 5);
    const char *loc = sstr("chr%s:%s-%s  (%s strand, %s, %s exons)", chromName(chrom), fmtInt(ri(g, 0, 6)), fmtInt(ri(g, 0, 7)),
                           ri(g, 0, 8) > 0 ? "+" : "-", fmtBp((double)(ri(g, 0, 7) - ri(g, 0, 6))), fmtInt(ri(g, 0, 13)));
    txt(F_MONO, 12, x, y, TH->faint, loc);
    y += 26;
    float bx = x;
    if (gxButtonEx(gxId("gp browse"), (Rectangle){ bx, y, 150, 30 }, "Open in browser", BTN_PRIMARY, true, NULL)) {
        double pad = (double)(ri(g, 0, 7) - ri(g, 0, 6)) * 0.06 + 300;
        g_br.selGene = g_gp.id;
        navLocus(chrom, (long long)((double)ri(g, 0, 6) - pad), (long long)((double)ri(g, 0, 7) + pad));
    }
    bx += 158;
    if (gxButtonEx(gxId("gp ens"), (Rectangle){ bx, y, 90, 30 }, "Ensembl", BTN_GHOST, true, "Open the gene in Ensembl (your web browser)"))
        openUrl(sstr("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s", rc(g, 0, 1)));
    bx += 96;
    if (*rc(g, 0, 18) && gxButtonEx(gxId("gp uni"), (Rectangle){ bx, y, 90, 30 }, "UniProt", BTN_GHOST, true, NULL))
        openUrl(sstr("https://www.uniprot.org/uniprotkb/%s/entry", rc(g, 0, 18)));
    bx += 96;
    if (gxButtonEx(gxId("gp gnomad"), (Rectangle){ bx, y, 90, 30 }, "gnomAD", BTN_GHOST, true, NULL))
        openUrl(sstr("https://gnomad.broadinstitute.org/gene/%s?dataset=gnomad_r4", rc(g, 0, 1)));
    bx += 96;
    if (gxButtonEx(gxId("gp more"), (Rectangle){ bx, y, 120, 30 }, "More links...", BTN_GHOST, true, "GeneCards, OMIM, NCBI, GTEx, MedlinePlus and more"))
        lkGeneMenu(rc(g, 0, 2));
    y += 48;

    float colW = (W - 20) / 2;
    /* what it does */
    Rectangle c1 = { x, y, colW, 0 };
    float fh = txtParaH(F_REG, TH->fsSmall, colW - 36, 19, *rc(g, 0, 16) ? rc(g, 0, 16) : "No function summary in UniProt for this gene.", 14);
    c1.height = 64 + fh + 16;
    /* constraint and conditions */
    Rectangle c2 = { x + colW + 20, y, colW, c1.height };
    gxCard(c1);
    gxCardTitle(c1, "What it does", *rc(g, 0, 15) ? rc(g, 0, 15) : NULL);
    txtPara(F_REG, TH->fsSmall, c1.x + 18, c1.y + 64, colW - 36, 19, TH->dim, *rc(g, 0, 16) ? rc(g, 0, 16) : "No function summary in UniProt for this gene.", 14);
    gxCard(c2);
    gxCardTitle(c2, "How much it tolerates damage", "gnomAD v4.1 constraint");
    float cy = c2.y + 66;
    if (!rnull(g, 0, 20)) {
        double lo = rf(g, 0, 20), pli = rf(g, 0, 19), mz = rf(g, 0, 21);
        txt(F_REG, TH->fsSmall, c2.x + 18, cy, TH->dim, "Loss of function (LOEUF)");
        Rectangle b = { c2.x + 200, cy + 3, c2.width - 290, 10 };
        gxMeter(b, (float)(1 - fmin(lo, 2) / 2), lo < 0.35 ? TH->bad : lo < 0.7 ? TH->warn : TH->ok);
        txtIn(F_SEMI, TH->fsSmall, (Rectangle){ b.x + b.width + 8, cy - 2, 70, 20 }, TH->text, 0, sstr("%.2f", lo));
        cy += 26;
        txtPara(F_REG, TH->fsSmall, c2.x + 18, cy, c2.width - 36, 19, TH->text,
                lo < 0.35 ? "Highly intolerant: healthy people almost never carry a broken copy, so a damaging variant here could matter."
                : lo < 0.7 ? "Somewhat intolerant of losing a copy." : "Tolerant: many healthy people carry a broken copy.", 3);
        cy += 44;
        cy = vpRow(c2.x + 18, cy, c2.width - 36, "pLI", sstr("%.2f (1 = intolerant)", pli), TH->text);
        if (!isnan(mz)) cy = vpRow(c2.x + 18, cy, c2.width - 36, "Missense Z", sstr("%.2f (above 3.09: few missense variants tolerated)", mz), TH->text);
    } else txt(F_REG, TH->fsSmall, c2.x + 18, cy, TH->faint, "No constraint data for this gene.");
    y += c1.height + 20;

    /* conditions */
    int nd = g_gp.dis ? g_gp.dis->nrow : 0;
    if (nd || (g_gp.acmg && g_gp.acmg->nrow) || *rc(g, 0, 17)) {
        float h = 64 + (float)(nd > 12 ? 12 : nd) * 24;
        if (*rc(g, 0, 17)) h += txtParaH(F_REG, TH->fsSmall, W - 36, 19, rc(g, 0, 17), 8) + 16;
        if (g_gp.acmg && g_gp.acmg->nrow) h += 50;
        Rectangle c = { x, y, W, h };
        gxCard(c);
        gxCardTitle(c, "Conditions linked to this gene", "HPO / OMIM / Orphanet; inheritance: AD dominant, AR recessive, XL X-linked");
        float ry = c.y + 64;
        if (g_gp.acmg && g_gp.acmg->nrow) {
            txtPara(F_REG, TH->fsSmall, c.x + 18, ry, W - 36, 19, TH->warn,
                    sstr("On the ACMG list of genes a clinical genome reports secondary findings in: %s (%s). Reported: %s.",
                         rc(g_gp.acmg, 0, 1), rc(g_gp.acmg, 0, 2), rc(g_gp.acmg, 0, 3)), 2);
            ry += 50;
        }
        for (int i = 0; i < nd && i < 12; i++) {
            const char *inh = rc(g_gp.dis, i, 1), *did = rc(g_gp.dis, i, 2);
            Rectangle row = { c.x + 10, ry - 1, W - 20, 24 };
            bool hv;
            if (gxBehave(gxIdN(31000 + i), row, &hv, NULL)) {
                if (!strncmp(did, "OMIM:", 5)) openUrl(sstr("https://omim.org/entry/%s", did + 5));
                else if (!strncmp(did, "ORPHA:", 6)) openUrl(sstr("https://www.orpha.net/en/disease/detail/%s", did + 6));
            }
            if (hv) { drRound(row, 5, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; gxTip(gxIdN(31500 + i), row, "Open this condition in OMIM or Orphanet (your web browser)"); }
            gxPill(c.x + 18, ry + 1, *inh ? inh : "?", TH->accent2, TH->fsTiny);
            txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 80, ry, W - 200, 22 }, TH->text, 0, rc(g_gp.dis, i, 0));
            txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + W - 130, ry, 110, 22 }, hv ? TH->accent : TH->faint, 2, did);
            ry += 24;
        }
        if (nd > 12) txt(F_REG, TH->fsTiny, c.x + 18, ry, TH->faint, sstr("and %d more", nd - 12));
        if (*rc(g, 0, 17)) txtPara(F_REG, TH->fsSmall, c.x + 18, ry + 8, W - 36, 19, TH->dim, rc(g, 0, 17), 8);
        y += h + 20;
    }
    /* transcript diagram */
    Rectangle dc = { x, y, W, 160 };
    gxCard(dc);
    gxCardTitle(dc, "Your variants on the main transcript", sstr("%s%s; coding and ClinVar-listed variants as stems (filled: both copies)", rc(g, 0, 9), ri(g, 0, 10) ? " (MANE Select)" : ""));
    geneDiagram((Rectangle){ dc.x + 24, dc.y + 62, dc.width - 48, 80 });
    y += 180;
    /* variants table */
    int nv = g_gp.vars ? g_gp.vars->nrow : 0;
    txt(F_SEMI, TH->fsH3, x, y, TH->text, sstr("Your variants in %s", sym));
    txt(F_REG, TH->fsSmall, x + 220, y + 2, TH->dim, sstr("%s in total, %s affecting the protein or its splicing", fmtInt(nv), fmtInt(ri(g, 0, 23))));
    y += 30;
    float th = fminf(600, 34 + (float)nv * 30 + 4);
    if (nv) {
        int clicked = gxTable(gxId("gene vars"), &g_gp.vt, (Rectangle){ x, y, W, th }, g_gvCols, 6, nv, 30, gvCell, NULL, NULL);
        if (clicked >= 0) g_gp.sel = ri(g_gp.vars, clicked, 0);
    } else txt(F_REG, TH->fsSmall, x, y, TH->faint, "You match the reference genome across this gene.");
    y += th + 30;
    gxScrollEnd(&g_gp.scroll, main, y - top);
    if (sideW > 0) {
        Rectangle sp = { r.x + r.width - sideW, r.y, sideW, r.height };
        drRect(sp, TH->panel);
        drVLine(sp.x, sp.y, sp.height, TH->line);
        if (gxIconButton(gxId("gp close"), (Rectangle){ sp.x + sp.width - 38, sp.y + 10, 28, 28 }, icClose, false, "Close")) g_gp.sel = -1;
        else variantPanel(sp, g_gp.sel, true);
    }
}

/* ---- the index ---- */
static const char *const g_giScopes[] = { "Worth a look", "ACMG genes", "Pharmacogenes", "Most variants" };
static void giQuery(void) {
    static const char *const sql[] = {
        "SELECT g.id, g.symbol, g.name, g.chrom, g.nvar, g.ncoding, g.loeuf, "
        " (SELECT count(*) FROM variant v WHERE v.gene_id=g.id AND (v.impact=3 OR v.clin IN (5,6))) AS hi "
        "FROM gene g WHERE g.id IN (SELECT gene_id FROM variant WHERE impact=3 OR clin IN (5,6,7) OR (impact=2 AND (af IS NULL OR af<0.01))) "
        "ORDER BY hi DESC, g.loeuf LIMIT 1500",
        "SELECT g.id, g.symbol, g.name, g.chrom, g.nvar, g.ncoding, g.loeuf, 0 FROM gene g WHERE g.symbol IN (SELECT gene FROM acmg) "
        "AND g.biotype='protein_coding' ORDER BY g.symbol",
        "SELECT g.id, g.symbol, g.name, g.chrom, g.nvar, g.ncoding, g.loeuf, 0 FROM gene g WHERE g.symbol IN (SELECT gene FROM pgx_gene) "
        "AND g.biotype='protein_coding' ORDER BY g.symbol",
        "SELECT g.id, g.symbol, g.name, g.chrom, g.nvar, g.ncoding, g.loeuf, 0 FROM gene g WHERE g.biotype='protein_coding' ORDER BY g.ncoding DESC LIMIT 1500",
    };
    if (g_giText[0]) {
        char like[100];
        snprintf(like, sizeof like, "%s%%", g_giText);
        const char *p[] = { like };
        dbAsync(SLOT_GENE, "SELECT id, symbol, name, chrom, nvar, ncoding, loeuf, 0 FROM gene WHERE symbol LIKE ?1 ORDER BY length(symbol), symbol LIMIT 1500", 1, p);
    } else dbAsync(SLOT_GENE, sql[g_giScope], 0, NULL);
}
static void giCell(int row, int col, Rectangle c, bool sel, void *u) {
    (void)sel; (void)u;
    Rows *r = g_giRows;
    switch (col) {
    case 0: txtIn(F_SEMI, TH->fs, c, TH->text, 0, rc(r, row, 1)); break;
    case 1: txtIn(F_REG, TH->fsSmall, c, TH->dim, 0, rc(r, row, 2)); break;
    case 2: txtIn(F_REG, TH->fsSmall, c, TH->dim, 0, sstr("chr%s", chromName((int)ri(r, row, 3)))); break;
    case 3: txtIn(F_REG, TH->fsSmall, c, TH->text, 2, fmtInt(ri(r, row, 4))); break;
    case 4: txtIn(F_REG, TH->fsSmall, c, TH->text, 2, fmtInt(ri(r, row, 5))); break;
    case 5: if (!rnull(r, row, 6)) { double lo = rf(r, row, 6); txtIn(F_REG, TH->fsSmall, c, lo < 0.35 ? TH->bad : lo < 0.7 ? TH->warn : TH->dim, 2, sstr("%.2f", lo)); } break;
    case 6: if (ri(r, row, 7)) txtIn(F_SEMI, TH->fsSmall, c, TH->bad, 2, fmtInt(ri(r, row, 7))); break;
    }
}
static const GxCol g_giCols[] = { { "GENE", 120, 0, false, NULL }, { "NAME", -3, 0, false, NULL }, { "CHROM", 70, 0, false, NULL },
                                  { "VARIANTS", 90, 2, false, NULL }, { "CODING", 80, 2, false, NULL },
                                  { "LOEUF", 70, 2, false, "Lower = less tolerant of losing a copy" }, { "NOTABLE", 80, 2, false, "High-impact or pathogenic variants" } };
static void pageGene(Rectangle r) {
    StLoc *l = navCur();
    if (l->geneId > 0) { gpLoad(l->geneId); geneDetail(r); return; }
    Rows *t;
    if ((t = dbTake(SLOT_GENE))) { rowsFree(&g_giRows); g_giRows = t; g_giT.scroll = (GxScroll){ 0 }; }
    if (g_giDirty) { g_giDirty = false; giQuery(); }
    float x = r.x + 24, y = r.y + 18;
    txt(F_SEMI, TH->fsH2, x, y, TH->text, "Genes");
    txt(F_REG, TH->fsSmall, x + 84, y + 6, TH->dim, "Click a gene for its page. 'Worth a look' lists genes where you carry a high-impact, pathogenic or rare protein-changing variant.");
    y += 40;
    if (gxTextField(gxId("gi text"), (Rectangle){ x, y, 240, 32 }, g_giText, sizeof g_giText, "Gene symbol", true)) g_giDirty = true;
    if (gxSegmented(gxId("gi scope"), (Rectangle){ x + 252, y, 520, 32 }, g_giScopes, 4, &g_giScope)) { g_giText[0] = 0; g_giDirty = true; }
    y += 46;
    int n = g_giRows ? g_giRows->nrow : 0;
    bool dbl;
    int clicked = gxTable(gxId("gi table"), &g_giT, (Rectangle){ r.x + 16, y, r.width - 32, r.y + r.height - y - 12 }, g_giCols, 7, n, 32, giCell, NULL, &dbl);
    if (clicked >= 0) navGene((int)ri(g_giRows, clicked, 0));
}
