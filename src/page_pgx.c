/* ========================================================================
   src/page_pgx.c - medications (pharmacogenomics, from PharmCAT), the
   curated traits, and the maternal and paternal lineages.
   ======================================================================== */

/* ---- medications ---- */
static Rows *g_pxGenes, *g_pxDrugs, *g_pxPgkb;
static int g_pxSrc, g_pxSel = -1, g_pxTab;
static char g_pxQ[80];
static GxScroll g_pxScroll;
static Color phenoColor(const char *ph) {
    if (!*ph || strstr(ph, "No Result") || strstr(ph, "Indeterminate")) return TH->faint;
    if (strstr(ph, "Poor") || strstr(ph, "Decreased") || strstr(ph, "Deficient") || strstr(ph, "Possible") || strstr(ph, "non-responsive")) return TH->bad;
    if (strstr(ph, "Intermediate") || strstr(ph, "Rapid") || strstr(ph, "Ultrarapid") || strstr(ph, "Increased") || strstr(ph, "Uncertain")) return TH->warn;
    return TH->ok;
}
static const char *const g_pxSrcs[] = { "CPIC", "DPWG", "FDA label", "All" };
static void pxSelect(long long i) { g_pxSel = g_pxSel == (int)i ? -1 : (int)i; }
/* a recommendation's source ("CPIC Guideline Annotation", "FDA Label Annotation"...) in plain words */
static const char *pxSrcName(const char *s) {
    return strstr(s, "CPIC") ? "CPIC guideline" : strstr(s, "DPWG") ? "DPWG guideline" : strstr(s, "Label") ? "FDA drug label"
         : strstr(s, "FDA") ? "FDA PGx table" : s;
}
static const char *const g_pxTabs[] = { "Genes", "Drugs", "PharmGKB annotations" };
static void pagePgx(Rectangle r) {
    if (!dbHasTable("pgx_gene")) { pageMissing(r, "Medications", "The medicines step did not run for this person. Re-import them from the People page."); return; }
    if (!g_pxGenes) {
        g_pxGenes = dbq("SELECT gene, diplotype, phenotype, activity, called, functions, notes, uncallable FROM pgx_gene "
                        "ORDER BY called DESC, gene", NULL);
        g_pxDrugs = dbq("SELECT drug, gene, source, phenotype, recommendation, implication, strength, guideline, url FROM pgx_drug "
                        "WHERE recommendation!='' ORDER BY drug, source", NULL);
        g_pxPgkb = dbHasTable("pgkb") ? dbq("SELECT rsid, gene, level, drugs, genotype, text, category FROM pgkb WHERE level IN ('1A','1B','2A','2B') "
                                            "ORDER BY level, gene", NULL) : NULL;
    }
    float top = gxScrollBegin(&g_pxScroll, r);
    float x = r.x + 32, W = fminf(1300, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Medications",
                         sstr("How your genes may change the way some medicines work, from %s. Never change a medicine on your own: "
                              "show these to a doctor or pharmacist.", dbMeta("pgx_tool")));
    /* what stands out: genes whose result changes guidance */
    {
        char items[1600] = "";
        int n = 0;
        for (int i = 0; i < g_pxGenes->nrow; i++) {
            const char *ph = rc(g_pxGenes, i, 2);
            if (!ri(g_pxGenes, i, 4) || !*ph) continue;
            Color c = phenoColor(ph);
            bool notable = c.r == TH->bad.r && c.g == TH->bad.g;
            notable = notable || (c.r == TH->warn.r && c.g == TH->warn.g) || strstr(ph, "positive");
            if (!notable || strstr(ph, "Uncertain")) continue;
            if (n) strncat(items, "   |   ", sizeof items - strlen(items) - 1);
            strncat(items, sstr("%s %s: %s", rc(g_pxGenes, i, 0), rc(g_pxGenes, i, 1), ph), sizeof items - strlen(items) - 1);
            n++;
        }
        if (n) {
            float h = 56 + txtParaH(F_REG, TH->fs, W - 40, 22, items, 6);
            Rectangle c = { x, y, W, h };
            gxCard(c);
            drRound((Rectangle){ c.x, c.y + 12, 4, h - 24 }, 2, TH->warn);
            txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, "What stands out for you");
            txtPara(F_REG, TH->fs, c.x + 20, c.y + 44, W - 40, 22, TH->dim, items, 6);
            y += h + 18;
        }
    }
    if (gxButtonEx(gxId("pgx export"), (Rectangle){ x + W - 260, y - 10, 260, 32 }, "Summary for a doctor (HTML)", BTN_PRIMARY, true,
                   "Writes a printable summary of your medication genes, HLA type and clinical findings into your Documents folder and opens it"))
        exportReport();
    if (gxTabs(gxId("pgx tabs"), x - 14, y - 6, 40, g_pxTabs, 3, &g_pxTab)) g_pxSel = -1;
    y += 46;
    if (g_pxTab == 0) {
        int cols = W > 1100 ? 4 : 3;
        float cw = (W - 16 * (float)(cols - 1)) / (float)cols, ch = 132;
        for (int i = 0; i < g_pxGenes->nrow; i++) {
            Rectangle c = { x + (cw + 16) * (float)(i % cols), y + (ch + 16) * (float)(i / cols), cw, ch };
            bool h;
            if (gxBehave(gxIdN(30000 + i), c, &h, NULL)) {
                const char *gn = rc(g_pxGenes, i, 0);
                lkBegin(sstr("%s %s", gn, rc(g_pxGenes, i, 1)));
                lkHead("In Strand Open");
                lkFn(g_pxSel == i ? "Hide the PharmCAT notes" : "Show the PharmCAT notes", pxSelect, i);
                lkPgxGene(gn);
                lkGene(gn);
                lkOpen();
            }
            gxCard(c);
            if (h || g_pxSel == i) drRoundLine(c, 10, TH->accent);
            if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
            const char *ph = rc(g_pxGenes, i, 2);
            bool called = ri(g_pxGenes, i, 4);
            txt(F_SEMI, TH->fsH3, c.x + 16, c.y + 14, TH->text, rc(g_pxGenes, i, 0));
            txtIn(F_SEMI, TH->fsH3, (Rectangle){ c.x + 16, c.y + 40, cw - 32, 24 }, called ? TH->text : TH->faint, 0, called ? rc(g_pxGenes, i, 1) : "Not called");
            const char *phs = *ph ? ph : called ? "No phenotype assigned" : "";
            if (*phs) txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 16, c.y + 70, cw - 32, 20 }, phenoColor(ph), 0, phs);
            if (*rc(g_pxGenes, i, 3)) txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 16, c.y + 96, cw - 32, 18 }, TH->faint, 0, sstr("activity score %s", rc(g_pxGenes, i, 3)));
            if (!called) {
                const char *why = !strcmp(rc(g_pxGenes, i, 0), "CYP2D6") ? "needs copy-number calling from reads" :
                                  !strncmp(rc(g_pxGenes, i, 0), "HLA", 3) ? "needs HLA typing" : "not determinable from these calls";
                txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 16, c.y + 96, cw - 32, 18 }, TH->faint, 0, why);
            }
        }
        y += (ch + 16) * (float)((g_pxGenes->nrow + cols - 1) / cols) + 10;
        if (g_pxSel >= 0) {
            float h = 40 + txtParaH(F_REG, TH->fsSmall, W - 40, 19, rc(g_pxGenes, g_pxSel, 6), 20) + 60;
            Rectangle c = { x, y, W, h };
            gxCard(c);
            txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, sstr("%s  %s", rc(g_pxGenes, g_pxSel, 0), rc(g_pxGenes, g_pxSel, 1)));
            float yy = c.y + 46;
            if (*rc(g_pxGenes, g_pxSel, 5)) { txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 20, yy, W - 40, 20 }, TH->dim, 0, sstr("Allele function: %s", rc(g_pxGenes, g_pxSel, 5))); yy += 26; }
            txtPara(F_REG, TH->fsSmall, c.x + 20, yy, W - 40, 19, TH->dim, *rc(g_pxGenes, g_pxSel, 6) ? rc(g_pxGenes, g_pxSel, 6) : "No notes from PharmCAT.", 20);
            y += h + 16;
        }
    } else if (g_pxTab == 1) {
        if (gxSegmented(gxId("pgx src"), (Rectangle){ x, y, 380, 32 }, g_pxSrcs, 4, &g_pxSrc)) {}
        gxTextField(gxId("pgx q"), (Rectangle){ x + 396, y, 260, 32 }, g_pxQ, sizeof g_pxQ, "Find a medicine", true);
        y += 48;
        static const char *const srcKey[] = { "CPIC", "DPWG", "FDA", "" };
        for (int i = 0; i < g_pxDrugs->nrow; i++) {
            if (*srcKey[g_pxSrc] && !strstr(rc(g_pxDrugs, i, 2), srcKey[g_pxSrc])) continue;
            if (g_pxQ[0] && !strcasestr(rc(g_pxDrugs, i, 0), g_pxQ)) continue;
            const char *rec = rc(g_pxDrugs, i, 4);
            bool std = strcasestr(rec, "standard") || strcasestr(rec, "label-recommended") || strcasestr(rec, "no reason to avoid") ||
                       strcasestr(rec, "no adjustment") || strcasestr(rec, "per standard") || strcasestr(rec, "no action");
            float h = 34 + txtParaH(F_REG, TH->fsSmall, W - 290, 19, rec, 8) + 8 + txtParaH(F_REG, TH->fsTiny + 1, W - 290, 17, rc(g_pxDrugs, i, 5), 4) + 12;
            Rectangle c = { x, y, W, h };
            gxCard(c);
            drRound((Rectangle){ c.x, c.y + 12, 4, h - 24 }, 2, std ? TH->ok : TH->warn);
            txt(F_SEMI, TH->fsH3, c.x + 18, c.y + 12, TH->text, rc(g_pxDrugs, i, 0));
            txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 18, c.y + 36, 220, 16 }, TH->faint, 0, rc(g_pxDrugs, i, 2));
            txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 18, c.y + 54, 220, 16 }, TH->dim, 0, rc(g_pxDrugs, i, 3));
            if (*rc(g_pxDrugs, i, 6)) txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 18, c.y + 72, 220, 16 }, TH->faint, 0, sstr("strength: %s", rc(g_pxDrugs, i, 6)));
            float ty = c.y + 12;
            ty += txtPara(F_REG, TH->fsSmall, c.x + 260, ty, W - 290, 19, TH->text, rec, 8) + 8;
            txtPara(F_REG, TH->fsTiny + 1, c.x + 260, ty, W - 290, 17, TH->dim, rc(g_pxDrugs, i, 5), 4);
            if (gxBoxClick(gxIdN(30500 + i), c)) {
                char g1[40] = "";
                sscanf(rc(g_pxDrugs, i, 1), "%39[A-Za-z0-9-]", g1);
                lkBegin(rc(g_pxDrugs, i, 0));
                if (*rc(g_pxDrugs, i, 8)) lkUrl(sstr("This advice on ClinPGx (%s)", pxSrcName(rc(g_pxDrugs, i, 2))), "clinpgx.org", rc(g_pxDrugs, i, 8));
                lkDrug(rc(g_pxDrugs, i, 0));
                if (g1[0]) lkPgxGene(g1);
                lkOpen();
            }
            y += h + 10;
        }
    } else {
        Rows *p = g_pxPgkb;
        y += txtPara(F_REG, TH->fsSmall, x, y, W, 19, TH->dim, "PharmGKB's clinical annotations (evidence levels 1A to 2B) for single variants, with the text for your own genotype. "
                     "Level 1 means a guideline or strong replicated evidence; level 2, moderate.", 3) + 12;
        for (int i = 0; p && i < p->nrow; i++) {
            float h = 30 + txtParaH(F_REG, TH->fsSmall, W - 300, 19, rc(p, i, 5), 6) + 14;
            Rectangle c = { x, y, W, h };
            gxCard(c);
            txt(F_SEMI, TH->fs, c.x + 16, c.y + 12, TH->text, rc(p, i, 3));
            txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 16, c.y + 34, 260, 16 }, TH->faint, 0, sstr("%s  %s  your genotype %s", rc(p, i, 1), rc(p, i, 0), rc(p, i, 4)));
            gxPill(c.x + 16, c.y + 56, sstr("level %s", rc(p, i, 2)), rc(p, i, 2)[0] == '1' ? TH->accent : TH->accent2, TH->fsTiny);
            txtPara(F_REG, TH->fsSmall, c.x + 280, c.y + 12, W - 300, 19, TH->text, rc(p, i, 5), 6);
            if (gxBoxClick(gxIdN(30900 + i), c)) {
                lkBegin(sstr("%s %s", rc(p, i, 1), rc(p, i, 0)));
                lkRsid(atoll(rc(p, i, 0) + (rc(p, i, 0)[0] == 'r' ? 2 : 0)));
                lkUrl(sstr("%s on ClinPGx (PharmGKB)", rc(p, i, 0)), "clinpgx.org", sstr("https://www.clinpgx.org/search?query=%s", rc(p, i, 0)));
                lkOpen();
            }
            y += h + 10;
        }
    }
    y += 30;
    gxScrollEnd(&g_pxScroll, r, y - top);
}

/* ---- traits ---- */
static Rows *g_trRows;
static int g_trCat;
static GxScroll g_trScroll;
static const char *const g_trCats[] = { "All", "Appearance", "Food and drink", "Body" };
static int g_trFound;
static void pageTraits(Rectangle r) {
    if (!dbHasTable("trait")) { pageMissing(r, "Traits", "The traits step did not run for this person. Re-import them from the People page."); return; }
    if (!g_trRows) g_trRows = dbq("SELECT id, cat, title, result, detail, conf, genotype, about, source, status, copies FROM trait WHERE cat!='Health' ORDER BY ord", NULL);
    float top = gxScrollBegin(&g_trScroll, r);
    float x = r.x + 32, W = fminf(1300, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Traits",
                         "What a few well-studied variants say about how you look, taste and respond. Most traits depend on many genes and on life, "
                         "so these are tendencies; the confidence tag says how much one variant tells.");
    gxSegmented(gxId("trait cat"), (Rectangle){ x, y - 4, 520, 32 }, g_trCats, 4, &g_trCat);
    bool wrap = 536 + 360 > W;   /* the found filter on its own row when both do not fit */
    gxSegmented(gxId("trait found"), (Rectangle){ wrap ? x : x + 536, wrap ? y + 36 : y - 4, 360, 32 }, g_foundOpts, 3, &g_trFound);
    y += wrap ? 84 : 44;
    int cols = W > 1100 ? 3 : 2;
    float cw = (W - 16 * (float)(cols - 1)) / (float)cols;
    float colY[3] = { y, y, y };
    for (int i = 0; i < g_trRows->nrow; i++) {
        if (g_trCat && strcmp(rc(g_trRows, i, 1), g_trCats[g_trCat])) continue;
        if (g_trFound && traitFound(g_trRows, i, 10, 3) != (g_trFound == 1)) continue;
        int col = 0;
        for (int k = 1; k < cols; k++) if (colY[k] < colY[col]) col = k;
        float h = 82 + txtParaH(F_REG, TH->fsSmall, cw - 36, 19, rc(g_trRows, i, 4), 8) + 10 + txtParaH(F_REG, TH->fsTiny + 1, cw - 36, 17, rc(g_trRows, i, 7), 6) + 46;
        Rectangle c = { x + (cw + 16) * (float)col, colY[col], cw, h };
        gxCard(c);
        const char *conf = rc(g_trRows, i, 5);
        Color cc = !strcmp(conf, "high") ? TH->ok : !strcmp(conf, "moderate") ? TH->warn : TH->faint;
        txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 18, c.y + 14, cw - 140, 20 }, TH->dim, 0, rc(g_trRows, i, 2));
        gxPill(c.x + cw - gxPillW(sstr("%s confidence", conf), TH->fsTiny) - 16, c.y + 12, sstr("%s confidence", conf), cc, TH->fsTiny);
        txtIn(F_SEMI, TH->fsH2, (Rectangle){ c.x + 18, c.y + 40, cw - 36, 28 }, strcmp(rc(g_trRows, i, 9), "nocall") ? TH->text : TH->faint, 0, rc(g_trRows, i, 3));
        float ty = c.y + 80;
        ty += txtPara(F_REG, TH->fsSmall, c.x + 18, ty, cw - 36, 19, TH->text, rc(g_trRows, i, 4), 8) + 10;
        ty += txtPara(F_REG, TH->fsTiny + 1, c.x + 18, ty, cw - 36, 17, TH->faint, rc(g_trRows, i, 7), 6) + 8;
        txtIn(F_MONO, 11, (Rectangle){ c.x + 18, ty, cw - 36, 16 }, TH->dim, 0, rc(g_trRows, i, 6));
        txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 18, ty + 18, cw - 36, 16 }, TH->faint, 0, rc(g_trRows, i, 8));
        if (gxBoxClick(gxIdN(31500 + i), c)) { lkBegin(rc(g_trRows, i, 2)); lkTrait(rc(g_trRows, i, 2), rc(g_trRows, i, 6)); lkOpen(); }
        colY[col] += h + 16;
    }
    y = fmaxf(colY[0], fmaxf(colY[1], colY[2])) + 20;
    if (!g_trCat || g_trCat == 3) {
        static Rows *hla, *del;
        if (!hla && dbHasTable("hla")) hla = dbq("SELECT gene, a1, a2 FROM hla WHERE gene IN ('HLA-A','HLA-B','HLA-C','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DPB1') "
                                                 "ORDER BY instr('HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQA1 HLA-DQB1 HLA-DPB1', gene)", NULL);
        if (!del && dbHasTable("genedel")) del = dbq("SELECT gene, call, note FROM genedel", NULL);
        if (hla && hla->nrow) {
            Rectangle c = { x, y, W, 0 };
            const char *note = "Your HLA (tissue) type, read from the sequencing reads of the MHC region. HLA genes present pieces of proteins to the immune system; "
                               "they matter for transplants, some drug reactions (see Medications) and many autoimmune conditions. "
                               "DQ2.5 is DQA1*05 with DQB1*02; DQ8 is DQA1*03 with DQB1*03:02: the coeliac types.";
            c.height = 70 + 26 * 2 + txtParaH(F_REG, TH->fsSmall, W - 40, 19, note, 5) + 10;
            gxCard(c);
            gxCardTitle(c, "HLA type", "two-field alleles, typed with T1K against IPD-IMGT/HLA");
            float hx = c.x + 20, hy = c.y + 64;
            for (int i = 0; i < hla->nrow; i++) {
                float cw2 = (W - 40) / (float)hla->nrow;
                txt(F_SEMI, TH->fsTiny, hx + cw2 * (float)i, hy, TH->faint, rc(hla, i, 0));
                txt(F_MONO, 13, hx + cw2 * (float)i, hy + 18, TH->text, rc(hla, i, 1));
                txt(F_MONO, 13, hx + cw2 * (float)i, hy + 36, TH->text, rc(hla, i, 2));
            }
            txtPara(F_REG, TH->fsSmall, c.x + 20, hy + 62, W - 40, 19, TH->dim, note, 5);
            y += c.height + 16;
        }
        if (del && del->nrow) {
            Rectangle c = { x, y, W, 64 + 30 * (float)del->nrow };
            gxCard(c);
            gxCardTitle(c, "Whole-gene deletions", "common deletions read from coverage");
            for (int i = 0; i < del->nrow; i++) {
                float yy = c.y + 60 + 30 * (float)i;
                txt(F_SEMI, TH->fs, c.x + 20, yy, TH->text, rc(del, i, 0));
                txt(F_SEMI, TH->fsSmall, c.x + 120, yy + 1, strstr(rc(del, i, 1), "both") ? TH->warn : TH->ok, rc(del, i, 1));
                txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 320, yy - 2, W - 340, 20 }, TH->dim, 0, rc(del, i, 2));
                if (gxBoxClick(gxIdN(31900 + i), (Rectangle){ c.x + 8, yy - 6, W - 16, 30 })) lkGeneMenu(rc(del, i, 0));
            }
            y += c.height + 16;
        }
    }
    gxScrollEnd(&g_trScroll, r, y - top);
}

/* ---- lineage ---- */
static Rows *g_lnMt, *g_lnMtC, *g_lnY;
static float g_lnNeed[2];   /* content height each column used last frame (0 maternal, 1 paternal) */
static GxScroll g_lnScroll;
static void lineageColumn(Rectangle c, bool mt) {
    gxCard(c);
    const char *hg = dbMeta(mt ? "mt_haplogroup" : "y_haplogroup");
    txt(F_REG, TH->fsSmall, c.x + 20, c.y + 18, TH->dim, mt ? "Maternal line: mitochondrial DNA" : "Paternal line: Y chromosome");
    txt(F_SEMI, 34, c.x + 20, c.y + 40, TH->text, *hg ? hg : "-");
    float y = c.y + 92;
    if (!*hg) {
        /* nothing to trace: no calls for this line in the data */
        bool female = !strcmp(dbMeta("sex"), "female");
        txtPara(F_REG, TH->fsSmall, c.x + 20, y, c.width - 40, 19, TH->dim,
                mt ? "This person's file has no mitochondrial calls, so the maternal line cannot be traced. Some sequencing VCFs leave "
                     "the mitochondrial genome out; a chip export or the provider's full file usually includes it."
                   : female ? "Only men carry a Y chromosome, so a woman's file has no paternal line to trace. Her father's, brother's or "
                              "paternal uncle's DNA would show it."
                            : "This person's file has no Y chromosome calls, so the paternal line cannot be traced.", 5);
        return;
    }
    if (mt) {
        y += txtPara(F_REG, TH->fsSmall, c.x + 20, y, c.width - 40, 19, TH->dim,
                     sstr("Passed down from mother to child, unchanged except for rare mutations, so it traces your mother's mother's mother's line back tens of thousands "
                          "of years. Match quality %.0f%% against PhyloTree 17.", atof(dbMeta("mt_quality")) * 100), 4) + 14;
    } else {
        const char *f = dbMeta("y_formed"), *t = dbMeta("y_tmrca");
        y += txtPara(F_REG, TH->fsSmall, c.x + 20, y, c.width - 40, 19, TH->dim,
                     sstr("Passed from father to son. Your branch formed about %s years ago%s, on the YFull tree; the SNPs that define it are marked below.",
                          fmtInt(atoll(f)), *t && strcmp(t, f) ? sstr(" and its known members share an ancestor about %s years ago", fmtInt(atoll(t))) : ""), 4) + 14;
    }
    /* the path, root at the top */
    Rows *p = mt ? g_lnMt : g_lnY;
    int n = p ? p->nrow : 0;
    float lineX = c.x + 30;
    float y0 = y;
    for (int i = 0; i < n; i++) {
        const char *name = rc(p, i, 1);
        bool last = i == n - 1;
        bool hasDer = mt ? *rc(p, i, 3) != 0 : ri(p, i, 2) > 0;
        const char *ev = mt ? (*rc(p, i, 3) ? rc(p, i, 3) : rc(p, i, 2)) : rc(p, i, 4);
        float h = 22 + (*ev ? fminf(54, txtParaH(F_MONO, 11, c.width - 90, 15, ev, 3)) : 0) + 8;
        if (!mt && !hasDer && !last) h = 26;
        DrawCircleV((Vector2){ lineX, y + 9 }, last ? 7 : 4.5f, last ? TH->accent : hasDer ? TH->accent2 : TH->lineHi);
        if (i + 1 < n) drRect((Rectangle){ lineX - 1, y + 14, 2, h - 6 }, TH->lineHi);
        txt(last ? F_SEMI : F_REG, last ? TH->fs + 1 : TH->fs, lineX + 18, y, last ? TH->text : hasDer ? TH->text : TH->dim, name);
        float nx = lineX + 18 + (float)txtW(last ? F_SEMI : F_REG, last ? TH->fs + 1 : TH->fs, name) + 10;
        if (!mt) {
            const char *meta = ri(p, i, 2) || ri(p, i, 3) ? sstr("%lld derived, %lld ancestral", ri(p, i, 2), ri(p, i, 3)) : "no callable SNPs";
            const char *f = rc(p, i, 6);
            txt(F_REG, TH->fsTiny, nx, y + 3, TH->faint, sstr("%s%s", meta, *f && strcmp(f, "-") ? sstr("   formed %s ybp", fmtInt(atoll(f))) : ""));
        } else if (*rc(p, i, 4)) txt(F_REG, TH->fsTiny, nx, y + 3, TH->warn, sstr("not seen: %s", rc(p, i, 4)));
        if (*ev && (mt || hasDer)) txtPara(F_MONO, 11, lineX + 18, y + 22, c.width - 90, 15, mt ? TH->accent2 : TH->dim, ev, 3);
        y += h;
    }
    (void)y0;
    y += 10;
    if (mt) {
        const char *priv = dbMeta("mt_private");
        y += txtPara(F_REG, TH->fsSmall, c.x + 20, y, c.width - 40, 19, TH->dim,
                     sstr("Your own differences beyond %s (private or hotspot mutations): %s", hg, *priv ? priv : "none"), 4) + 10;
        if (g_lnMtC && g_lnMtC->nrow > 1) {
            txt(F_SEMI, TH->fsTiny, c.x + 20, y, TH->faint, "CLOSEST ALTERNATIVES");
            y += 20;
            for (int i = 1; i < g_lnMtC->nrow && i < 6; i++) {
                txt(F_REG, TH->fsSmall, c.x + 20, y, TH->dim, rc(g_lnMtC, i, 1));
                txt(F_REG, TH->fsSmall, c.x + 140, y, TH->faint, sstr("%.1f%%", rf(g_lnMtC, i, 2) * 100));
                y += 20;
            }
        }
    }
    g_lnNeed[mt ? 0 : 1] = y - c.y + 16;
}
/* ---- where your lines come from (step 19): stories, and ancient people on a map ---- */
static Rows *g_lnStory, *g_lnAnc[2], *g_lnNode[2];
static int g_lnLine;           /* 0 paternal (Y), 1 maternal (mtDNA) */
static bool g_lnAsked;
static const char *const g_lnLines[2] = { "Y", "MT" };
static int g_lnSpan[2];        /* branches from your major clade down to your own */

/* how close an ancient person is, relative to the length of your path below the major clade:
   2 your branch or nearly, 1 the nearer half, 0 the rest */
static int lnBand(int gap) {
    int span = g_lnSpan[g_lnLine] > 0 ? g_lnSpan[g_lnLine] : 1;
    return gap <= span / 5 ? 2 : gap <= span / 2 ? 1 : 0;
}

static void lnAsk(void) {
    g_lnAsked = true;
    if (!dbHasTable("lin_story")) return;
    g_lnStory = dbq("SELECT line, ord, text FROM lin_story ORDER BY line, ord", NULL);
    for (int k = 0; k < 2; k++) {
        g_lnAnc[k] = dbq("SELECT culture, place, country, lat, lon, years, hg, shared, shared_depth, gap, doi, grp FROM lin_ancient "
                         "WHERE line=? ORDER BY shared_depth DESC, years DESC", g_lnLines[k], NULL);
        g_lnNode[k] = dbq("SELECT depth, name, formed, tmrca, n_ancient FROM lin_node WHERE line=? ORDER BY depth DESC", g_lnLines[k], NULL);
        if (g_lnNode[k] && g_lnNode[k]->nrow) g_lnSpan[k] = (int)(ri(g_lnNode[k], 0, 0) - ri(g_lnNode[k], g_lnNode[k]->nrow - 1, 0));
    }
    if (g_lnAnc[0] && !g_lnAnc[0]->nrow && g_lnAnc[1] && g_lnAnc[1]->nrow) g_lnLine = 1;
}

static float lnStoryCard(float x, float y, float W) {
    if (!g_lnStory || !g_lnStory->nrow) return 0;
    float cw = (W - 20) / 2, H = 0;
    for (int pass = 0; pass < 2; pass++)
        for (int k = 0; k < 2; k++) {
            const char *line = k == 0 ? "MT" : "Y";        /* maternal on the left, as in the columns below */
            const char *hg = dbMeta(k == 0 ? "mt_haplogroup" : "y_haplogroup");
            float hh = 58;
            Rectangle c = { x + (cw + 20) * (float)k, y, cw, H };
            if (pass == 1) {
                gxCard(c);
                txt(F_REG, TH->fsSmall, c.x + 20, c.y + 18, TH->dim, k == 0 ? "Where your maternal line comes from" : "Where your paternal line comes from");
                txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 36, TH->text, *hg ? hg : (k == 0 ? "Not in this data" : "No Y chromosome"));
            }
            if (!*hg) {
                /* no story for a line the data does not have: say why instead of an empty card */
                bool female = !strcmp(dbMeta("sex"), "female");
                const char *why = k == 0 ? "This file has no mitochondrial calls, so the maternal line cannot be traced."
                                : female ? "Only men carry a Y chromosome, so a woman has no paternal line of her own to trace. Her father's, a "
                                           "brother's or a paternal uncle's DNA would show it: the same Y passes unchanged from father to son."
                                         : "This file has no Y chromosome calls, so the paternal line cannot be traced.";
                if (pass == 1) txtPara(F_REG, TH->fsSmall, c.x + 20, c.y + hh, cw - 40, 20, TH->dim, why, 6);
                hh += txtParaH(F_REG, TH->fsSmall, cw - 40, 20, why, 6) + 10;
            }
            for (int i = 0; i < g_lnStory->nrow; i++) {
                if (strcmp(rc(g_lnStory, i, 0), line)) continue;
                const char *t = rc(g_lnStory, i, 2);
                if (pass == 1) txtPara(F_REG, TH->fsSmall, c.x + 20, c.y + hh, cw - 40, 20, i == g_lnStory->nrow - 1 || !strncmp(t, "The ancient", 11) ? TH->text : TH->dim, t, 12);
                hh += txtParaH(F_REG, TH->fsSmall, cw - 40, 20, t, 12) + 10;
            }
            if (pass == 0) H = fmaxf(H, hh + 10);
        }
    return H;
}

/* the view for a set of points: the middle 96% of them, padded */
static GeoView lnView(Rows *r, Rectangle mr) {
    int n = r ? r->nrow : 0;
    if (!n) return geoFit(mr, -12, 40, 35, 65, 0);
    float *la = malloc(sizeof(float) * (size_t)n * 2), *lo = la + n;
    for (int i = 0; i < n; i++) { la[i] = (float)rf(r, i, 3); lo[i] = (float)rf(r, i, 4); }
    qsort(la, (size_t)n, sizeof(float), geoCmpF);
    qsort(lo, (size_t)n, sizeof(float), geoCmpF);
    int a = n / 50, b = n - 1 - n / 50;
    GeoView v = geoFit(mr, lo[a], lo[b], la[a], la[b], 4);
    free(la);
    return v;
}

static float lnAncientCard(float x, float y, float W) {
    Rows *rr[2] = { g_lnAnc[0], g_lnAnc[1] };
    if ((!rr[0] || !rr[0]->nrow) && (!rr[1] || !rr[1]->nrow)) return 0;
    float H = 640;
    Rectangle c = { x, y, W, H };
    gxCard(c);
    static const char *const segs[2] = { "Paternal line", "Maternal line" };
    gxSegmented(gxId("ln line"), (Rectangle){ c.x + W - 300, c.y + 16, 280, 30 }, segs, 2, &g_lnLine);
    Rows *r = rr[g_lnLine];
    int n = r ? r->nrow : 0;
    const char *major = dbMeta(g_lnLine == 0 ? "lin_y_major" : "lin_mt_major");
    txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 18, TH->text, g_lnLine == 0 ? "Ancient men on your paternal line" : "Ancient people on your maternal line");
    txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 20, c.y + 46, W - 360, 20 }, TH->dim, 0,
          n ? sstr("%s ancient people share a branch with you, from %s down to your own. Brighter dots share a more recent branch.", fmtInt(n), major)
            : "No ancient person in the database shares a branch with this line yet.");
    float mw = W * 0.56f;
    Rectangle mr = { c.x + 20, c.y + 80, mw, H - 120 };
    drRound(mr, 6, drMix(TH->card, TH->bg, 0.5f));
    static GeoView v[2];
    static float vw[2], vh[2];
    if (vw[g_lnLine] != mr.width || vh[g_lnLine] != mr.height) { v[g_lnLine] = lnView(r, mr); vw[g_lnLine] = mr.width; vh[g_lnLine] = mr.height; }
    GeoView gv = v[g_lnLine];
    gv.r = mr;
    gxClipPush(mr);
    geoLand(&gv, drMix(TH->panel2, TH->lineHi, 0.35f));
    /* dots: the farther relatives first, so the closest end on top */
    int hot = -1;
    float hotD = 64;
    for (int pass = 0; pass < 3; pass++)
        for (int i = n - 1; i >= 0; i--) {
            int band = lnBand((int)ri(r, i, 9));
            if (band != pass) continue;
            Vector2 p = geoXY(&gv, (float)rf(r, i, 4), (float)rf(r, i, 3));
            Color col = band == 2 ? TH->warn : band == 1 ? TH->accent : drFade(TH->dim, 0.55f);
            DrawCircleV(p, band == 2 ? 5.5f : band == 1 ? 4.0f : 2.6f, col);
            if (band == 2) DrawRing(p, 5.5f, 7.0f, 0, 360, 16, drFade(TH->text, 0.6f));
            float dx = p.x - g_gx.mouse.x, dy = p.y - g_gx.mouse.y, d = dx * dx + dy * dy;
            if (gxHover(mr) && d < hotD) { hotD = d; hot = i; }
        }
    if (hot >= 0) {
        Vector2 p = geoXY(&gv, (float)rf(r, hot, 4), (float)rf(r, hot, 3));
        DrawRing(p, 7, 9, 0, 360, 20, TH->text);
        const char *t1 = sstr("%s, %s", rc(r, hot, 1), rc(r, hot, 2));
        const char *t2 = sstr("%s%sabout %s years ago, %s (shares %s)", rc(r, hot, 0), *rc(r, hot, 0) ? ", " : "", fmtInt(ri(r, hot, 5)), rc(r, hot, 6), rc(r, hot, 7));
        float tw = fmaxf((float)txtW(F_SEMI, TH->fsSmall, t1), (float)txtW(F_REG, TH->fsTiny + 1, t2)) + 20;
        Rectangle tb = { fminf(p.x + 12, mr.x + mr.width - tw - 6), fmaxf(mr.y + 6, p.y - 50), tw, 44 };
        drRound(tb, 6, TH->panel2);
        drRoundLine(tb, 6, TH->lineHi);
        txtIn(F_SEMI, TH->fsSmall, (Rectangle){ tb.x + 10, tb.y + 4, tw - 20, 18 }, TH->text, 0, t1);
        txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ tb.x + 10, tb.y + 23, tw - 20, 16 }, TH->dim, 0, t2);
    }
    gxClipPop();
    txtIn(F_REG, TH->fsTiny, (Rectangle){ mr.x, mr.y + mr.height + 6, mr.width, 14 }, TH->faint, 0,
          "Map: Natural Earth. Places are burial sites; a dot is one ancient genome.");
    /* the list: closest branches first */
    float lx = mr.x + mr.width + 24, lw = c.x + W - 20 - lx, ly = c.y + 80;
    const char *last = "";
    int groups = 0, inGroup = 0, more = 0;
    float bottom = c.y + H - 44;          /* room for the source line */
    for (int i = 0; i < n; i++) {
        const char *sh = rc(r, i, 7);
        if (strcmp(sh, last)) {
            if (more) { txt(F_REG, TH->fsTiny + 1, lx + 12, ly, TH->faint, sstr("and %d more", more)); ly += 20; more = 0; }
            if (++groups > 7 || ly + 46 + 22 > bottom) break;
            last = sh; inGroup = 0; more = 0;
            int cnt = 0;
            for (int k = i; k < n && !strcmp(rc(r, k, 7), sh); k++) cnt++;
            const char *formed = "";
            for (int k = 0; g_lnNode[g_lnLine] && k < g_lnNode[g_lnLine]->nrow; k++)
                if (!strcmp(rc(g_lnNode[g_lnLine], k, 1), sh) && atoi(rc(g_lnNode[g_lnLine], k, 2)) > 0)
                    formed = sstr(", formed about %s years ago", fmtInt(atoll(rc(g_lnNode[g_lnLine], k, 2))));
            ly += 6;
            txtIn(F_SEMI, TH->fsSmall, (Rectangle){ lx, ly, lw, 18 }, lnBand((int)ri(r, i, 9)) == 2 ? TH->warn : TH->text, 0,
                  sstr("Shares %s%s", sh, formed));
            txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ lx, ly + 18, lw, 16 }, TH->faint, 0, sstr("%d ancient %s", cnt, cnt == 1 ? "person" : "people"));
            ly += 40;
        }
        if (inGroup >= 4 || ly + 22 > bottom) { more++; continue; }
        inGroup++;
        Rectangle row = { lx, ly - 2, lw, 22 };
        bool h;
        if (gxBehave(gxIdN(36000 + i), row, &h, NULL)) {
            lkBegin(sstr("%s, %s", rc(r, i, 1), rc(r, i, 2)));
            if (*rc(r, i, 10)) lkUrl("The study", "doi.org", rc(r, i, 10));
            lkUrl("The place on OpenStreetMap", "openstreetmap.org",
                  sstr("https://www.openstreetmap.org/?mlat=%.4f&mlon=%.4f#map=9/%.4f/%.4f", rf(r, i, 3), rf(r, i, 4), rf(r, i, 3), rf(r, i, 4)));
            lkOpen();
        }
        if (h) { drRound(row, 4, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
                 Vector2 p = geoXY(&gv, (float)rf(r, i, 4), (float)rf(r, i, 3));
                 if (CheckCollisionPointRec(p, mr)) DrawRing(p, 7, 9, 0, 360, 20, TH->text); }
        txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ lx + 12, ly, lw - 12, 18 }, TH->text, 0,
              sstr("%s, %s - about %s years ago%s", *rc(r, i, 0) ? rc(r, i, 0) : rc(r, i, 11), rc(r, i, 2), fmtInt(ri(r, i, 5)),
                   strcmp(rc(r, i, 6), sh) ? sstr(" (%s)", rc(r, i, 6)) : ""));
        ly += 22;
    }
    if (more && ly + 18 <= bottom + 20) txt(F_REG, TH->fsTiny + 1, lx + 12, ly, TH->faint, sstr("and %d more", more));
    txtIn(F_REG, TH->fsTiny, (Rectangle){ lx, c.y + H - 30, lw, 14 }, TH->faint, 0, dbMeta("lineage_source"));
    return H;
}

/* ---- what each parent gave you: inheritance arithmetic for one genome ----
   Exact from the chromosome lengths (GRCh38), the gene table and the person's
   sex and calls; nothing here needs a relative's DNA, and the one thing that
   would (the real share from each grandparent) is said to be unknowable. */
static Rows *g_lnGenes, *g_lnSexVars;
static float lnParentsCard(float x, float y, float W, bool draw) {
    if (!g_lnGenes) {
        g_lnGenes = dbq("SELECT chrom, count(*) FROM gene WHERE biotype='protein_coding' GROUP BY chrom", NULL);
        g_lnSexVars = dbq("SELECT chrom, count(*) FROM variant WHERE chrom>=23 AND pass=1 AND zyg>0 GROUP BY chrom", NULL);
    }
    double A = 0, X = (double)g_chrom[23].len, Y = (double)g_chrom[24].len, M = (double)g_chrom[25].len;
    for (int c = 1; c <= 22; c++) A += (double)g_chrom[c].len;
    long long gA = 0, gX = 0, gY = 0, gM = 0, vX = 0, vY = 0, vM = 0;
    for (int i = 0; g_lnGenes && i < g_lnGenes->nrow; i++) {
        int c = (int)ri(g_lnGenes, i, 0);
        long long n = ri(g_lnGenes, i, 1);
        if (c <= 22) gA += n; else if (c == 23) gX = n; else if (c == 24) gY = n; else if (c == 25) gM = n;
    }
    for (int i = 0; g_lnSexVars && i < g_lnSexVars->nrow; i++) {
        int c = (int)ri(g_lnSexVars, i, 0);
        if (c == 23) vX = ri(g_lnSexVars, i, 1); else if (c == 24) vY = ri(g_lnSexVars, i, 1); else if (c == 25) vM = ri(g_lnSexVars, i, 1);
    }
    bool male = !strcmp(dbMeta("sex"), "male"), female = !strcmp(dbMeta("sex"), "female");
    double mom = male ? 100.0 * (A + X) / (2 * A + X + Y) : 50.0;
    static char p1[700], p2[600], p3[700], p4[600];
    if (male) {
        snprintf(p1, sizeof p1, "Your mother gave you one copy of every autosome and your X chromosome (%.0f million letters); your father gave you "
                 "one copy of every autosome and your Y (%.0f million). The X is nearly three times the size of the Y, so a man's DNA is a little "
                 "more his mother's. Counting protein-coding genes: %s copies from your mother (%s on the X) and %s from your father (%s on the Y).",
                 X / 1e6, Y / 1e6, fmtInt(gA + gX), fmtInt(gX), fmtInt(gA + gY), fmtInt(gY));
        snprintf(p2, sizeof p2, "Your X is a mix of your mother's two X chromosomes, one from each of her parents, so it carries DNA from both your "
                 "maternal grandparents and none from your father's side. Your Y came from your father, his father and so on back, unchanged "
                 "except for new mutations.");
    } else if (female) {
        snprintf(p1, sizeof p1, "Each parent gave you exactly half of your chromosomes: one copy of every autosome and one X each, so %s "
                 "protein-coding gene copies from each.", fmtInt(gA + gX));
        snprintf(p2, sizeof p2, "The X from your father is the one he got from his mother, passed on whole, so it all comes from your paternal "
                 "grandmother (a mix of her own two X chromosomes): your father's father gave you no X DNA. The X from your mother is a mix "
                 "from both her parents.");
    } else {
        snprintf(p1, sizeof p1, "A woman gets an X from each parent, so her chromosomes are half her mother's and half her father's. A man gets his "
                 "mother's X (%.0f million letters) and his father's much smaller Y (%.0f million), so about %.1f%% of his DNA is his mother's.",
                 X / 1e6, Y / 1e6, 100.0 * (A + X) / (2 * A + X + Y));
        p2[0] = 0;
    }
    snprintf(p3, sizeof p3, "All of your mitochondrial DNA came from your mother: only %s letters and %lld genes, but copied several hundred "
             "to thousands of times in each cell (most in heart, muscle and brain) and over 100,000 times in an egg. %s",
             fmtInt((long long)M), gM ? gM : 13,
             (vX || vY || vM) ? sstr("In your data: %s%s%s.",
                                     vX ? sstr("%s variant%s on your X%s", fmtInt(vX), vX == 1 ? "" : "s", male ? " (all from your mother)" : " (from both parents)") : "",
                                     vY ? sstr("%s%s on your Y (from your father)", vX ? ", " : "", fmtInt(vY)) : "",
                                     vM ? sstr("%s%s in your mitochondria (from your mother)", (vX || vY) ? ", " : "", fmtInt(vM)) : "") : "");
    snprintf(p4, sizeof p4, "On average a quarter of your DNA comes from each grandparent, but the real share varies, more on your father's "
             "side because men pass chromosomes on in bigger pieces: in simulations with real crossover data, about 1 person in 200 has "
             "under 10 percent from one of their father's parents. Your own split cannot be read from one genome; it would take a "
             "grandparent's or a parent's DNA.");
    float tw = W - 300;
    float h = 66 + txtParaH(F_REG, TH->fsSmall, tw, 19, p1, 8) + 10 + (p2[0] ? txtParaH(F_REG, TH->fsSmall, tw, 19, p2, 6) + 10 : 0)
              + txtParaH(F_REG, TH->fsSmall, tw, 19, p3, 6) + 10 + txtParaH(F_REG, TH->fsSmall, tw, 19, p4, 6) + 40;
    h = fmaxf(h, 230);
    if (!draw) return h;
    Rectangle c = { x, y, W, h };
    gxCard(c);
    gxCardTitle(c, "What each parent gave you", "From your own genome: no relative's DNA needed");
    /* the split, as two big numbers and a bar */
    float bx = c.x + 20, by = c.y + 70;
    txt(F_REG, TH->fsSmall, bx, by, TH->dim, "From your mother");
    txt(F_SEMI, TH->fsH1, bx, by + 18, TH->text, sstr("%.1f%%", mom));
    txt(F_REG, TH->fsSmall, bx + 130, by, TH->dim, "From your father");
    txt(F_SEMI, TH->fsH1, bx + 130, by + 18, TH->text, sstr("%.1f%%", 100 - mom));
    Rectangle bar = { bx, by + 66, 230, 10 };
    drRound(bar, 5, TH->accent2);
    drRound((Rectangle){ bar.x, bar.y, bar.width * (float)(mom / 100), bar.height }, 5, TH->accent);
    txtPara(F_REG, TH->fsTiny + 1, bx, by + 86, 240, 16, TH->faint,
            male ? "of your nuclear DNA, by length; plus all your mitochondria from your mother" :
            female ? "of your nuclear DNA; plus all your mitochondria from your mother" : "for a woman; a man's split depends on his X and Y", 4);
    float tx = c.x + 280, ty = c.y + 66;
    ty += txtPara(F_REG, TH->fsSmall, tx, ty, tw, 19, TH->text, p1, 8) + 10;
    if (p2[0]) ty += txtPara(F_REG, TH->fsSmall, tx, ty, tw, 19, TH->text, p2, 6) + 10;
    ty += txtPara(F_REG, TH->fsSmall, tx, ty, tw, 19, TH->text, p3, 6) + 10;
    txtPara(F_REG, TH->fsSmall, tx, ty, tw, 19, TH->dim, p4, 6);
    txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 20, c.y + c.height - 24, W - 40, 14 }, TH->faint, 0,
          "Chromosome lengths from GRCh38 and genes from Ensembl; mitochondrial copies: Filograna et al. 2020 (FEBS Letters); grandparent "
          "shares: Coop lab simulations from observed crossovers (2013).");
    return h;
}

static void pageLineage(Rectangle r) {
    if (!dbHasTable("hap_mt")) { pageMissing(r, "Lineage", "The lineage step did not run for this person. Re-import them from the People page."); return; }
    if (!g_lnMt) {
        g_lnMt = dbq("SELECT depth, name, polys, found, missing FROM hap_mt ORDER BY depth", NULL);
        g_lnMtC = dbq("SELECT rank, name, score FROM hap_mt_cand ORDER BY rank", NULL);
        g_lnY = dbq("SELECT depth, name, nder, nanc, der, anc, formed, tmrca, isogg FROM hap_y WHERE depth>0 ORDER BY depth", NULL);
    }
    float top = gxScrollBegin(&g_lnScroll, r);
    float x = r.x + 32, W = r.width - 64;
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Lineage",
                         "Two lines of descent that pass through the generations almost unchanged: your mother's line in the mitochondria "
                         "and your father's line on the Y chromosome. Each is one ancestor in each generation, out of thousands.");
    if (!g_lnAsked) lnAsk();
    float sh = lnStoryCard(x, y, W);
    if (sh > 0) y += sh + 20;
    float ah = lnAncientCard(x, y, W);
    if (ah > 0) y += ah + 24;
    float cw = (W - 20) / 2;
    /* a line missing from the data (no MT calls, or no Y) only needs room for the explanation */
    float hm = *dbMeta("mt_haplogroup") ? 160 + 60 * (float)(g_lnMt ? g_lnMt->nrow : 0) + 260 : 200;
    float hy = *dbMeta("y_haplogroup") ? 160 + 50 * (float)(g_lnY ? g_lnY->nrow : 0) + 40 : 200;
    float h = fmaxf(fmaxf(hm, hy), fmaxf(g_lnNeed[0], g_lnNeed[1]));
    lineageColumn((Rectangle){ x, y, cw, h }, true);
    lineageColumn((Rectangle){ x + cw + 20, y, cw, h }, false);
    for (int k = 0; k < 2; k++) {
        Rectangle c = { x + (cw + 20) * (float)k, y, cw, h };
        const char *hg = dbMeta(k == 0 ? "mt_haplogroup" : "y_haplogroup");
        if (*hg && gxBoxClick(gxIdN(32500 + k), c)) {
            lkBegin(sstr("%s haplogroup %s", k == 0 ? "Maternal" : "Paternal", hg));
            if (k == 0) {
                lkUrl(sstr("%s on the YFull mtDNA tree", hg), "yfull.com", sstr("https://www.yfull.com/mtree/%s/", hg));
                lkUrl(sstr("%s on the FamilyTreeDNA tree", hg), "familytreedna.com", sstr("https://discover.familytreedna.com/mtdna/%s/tree", hg));
                lkUrl("PhyloTree (the reference tree)", "phylotree.org", "https://www.phylotree.org/tree/index.htm");
            } else {
                lkUrl(sstr("%s on the YFull tree", hg), "yfull.com", sstr("https://www.yfull.com/tree/%s/", hg));
                lkUrl(sstr("%s on the FamilyTreeDNA tree", hg), "familytreedna.com", sstr("https://discover.familytreedna.com/y-dna/%s/tree", hg));
                lkUrl("Where it is found (FamilyTreeDNA globe)", "familytreedna.com", sstr("https://discover.familytreedna.com/y-dna/%s/globe", hg));
            }
            lkUrl(sstr("Haplogroup %.1s (%s) on Wikipedia", hg, k == 0 ? "mtDNA" : "Y-DNA"), "wikipedia.org",
                  sstr("https://en.wikipedia.org/w/index.php?search=Haplogroup+%.1s+%s", hg, k == 0 ? "(mtDNA)" : "(Y-DNA)"));
            lkOpen();
        }
    }
    y += h + 20;
    /* what each parent gave you: last, so the stories and the map stay at the top */
    float ph = lnParentsCard(x, y, W, false);
    lnParentsCard(x, y, W, true);
    y += ph + 30;
    gxScrollEnd(&g_lnScroll, r, y - top);
}
