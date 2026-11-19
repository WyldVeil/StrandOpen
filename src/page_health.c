/* ========================================================================
   src/page_health.c - clinical findings and carrier status.

   A ClinVar classification is only half the story: what a variant means for
   the person depends on how the condition is inherited and how many copies
   they carry. Each finding here is read against the gene's inheritance (HPO)
   so a single copy in a recessive gene says "carrier" and a dominant one says
   "may matter for you". Tabs: pathogenic and likely pathogenic, risk factors,
   the ACMG secondary-findings genes, well-known variants, structural changes,
   and the GWAS research associations.
   ======================================================================== */

static int g_hTab;
static int g_hFound;
static Rows *g_hChipRare;   /* chip calls of rare ClinVar-pathogenic variants, set aside (step 5c) */
static Rows *g_hPath, *g_hRisk, *g_hAcmg, *g_hNotable, *g_hSv, *g_hDel, *g_hConf, *g_hRare;
static bool g_hAsked;
static bool g_hRareRv;   /* rare_var (step 16) was there when g_hRare was read */
static GxScroll g_hScroll;
static GxTable g_hGwasT = { .sel = -1, .sortCol = -1 };
static Rows *g_hGwas;
static char g_hGwasQ[120];
static bool g_hGwasDirty = true;
static int g_hGwasMine = 1;

static const char *inhMeaning(const char *inh, int zyg, bool *relevant) {
    bool ad = strstr(inh, "AD") != NULL, ar = strstr(inh, "AR") != NULL, xl = strstr(inh, "XL") != NULL || strstr(inh, "XLR") != NULL;
    *relevant = false;
    if (zyg == 2 || zyg == 3) {
        *relevant = ar || xl || ad;
        if (zyg == 3 && xl) return "One copy on your single X: for an X-linked condition this can be enough to cause it.";
        if (zyg == 3) return "One copy, on a chromosome you have only one of (X or Y): there is no second copy to make up for it. Worth discussing with a doctor.";
        if (ar) return "Two copies: for a recessive condition this can be enough to cause it. Worth discussing with a doctor.";
        return "Two copies of a disease-linked variant.";
    }
    if (ad && !ar) { *relevant = true; return "Dominant condition: one copy can be enough. Penetrance varies; worth discussing with a doctor or genetic counsellor."; }
    if (ar && !ad) return "Recessive condition: with one copy you are a carrier and usually unaffected. A partner who also carries one could pass both on.";
    if (ad && ar) { *relevant = true; return "Linked to both dominant and recessive forms: one copy may matter depending on the condition."; }
    return "Inheritance not recorded for this gene.";
}

static void hAsk(void) {
    g_hAsked = true;
    rowsFree(&g_hPath); rowsFree(&g_hRisk); rowsFree(&g_hAcmg); rowsFree(&g_hNotable); rowsFree(&g_hSv); rowsFree(&g_hDel); rowsFree(&g_hConf); rowsFree(&g_hChipRare);
    g_hPath = dbq("SELECT v.id, v.gene, v.hgvs_p, v.csq, v.zyg, v.clin, v.clin_sig, v.clin_dn, v.clin_stars, v.af, v.am, v.chrom, v.pos, "
                  "(SELECT group_concat(DISTINCT d.inheritance) FROM gene_disease d WHERE d.gene=v.gene_id), v.pass, v.dp, v.ad_alt, v.ad_ref, v.rsid "
                  "FROM variant v WHERE v.clin IN (5,6) ORDER BY v.clin DESC, v.clin_stars DESC", NULL);
    g_hChipRare = dbq("SELECT v.id, v.gene, v.hgvs_p, v.csq, v.zyg, v.clin_dn FROM variant v WHERE v.clin IN (105,106) ORDER BY v.gene", NULL);
    g_hConf = dbq("SELECT v.id, v.gene, v.hgvs_p, v.csq, v.zyg, v.clin_sig, v.clin_dn, v.af FROM variant v WHERE v.clin=4 AND "
                  "v.clin_sig LIKE '%athogenic%' AND v.impact>=2 AND (v.af IS NULL OR v.af < 0.02) ORDER BY v.af LIMIT 80", NULL);
    g_hRisk = dbq("SELECT v.id, v.gene, v.hgvs_p, v.csq, v.zyg, v.clin, v.clin_sig, v.clin_dn, v.clin_stars, v.af, v.rsid FROM variant v "
                  "WHERE v.clin IN (7,9,10) OR ((v.clin_flags & 1) AND v.clin NOT IN (1,2)) ORDER BY v.clin_stars DESC, v.clin LIMIT 400", NULL);
    if (dbHasTable("acmg"))
        g_hAcmg = dbq("SELECT a.gene, a.category, a.disease, a.inheritance, a.report, g.id, "
                      "(SELECT count(*) FROM variant v WHERE v.gene_id=g.id AND v.clin IN (5,6)), "
                      "(SELECT count(*) FROM variant v WHERE v.gene_id=g.id AND v.impact>=2 AND (v.af IS NULL OR v.af<0.01)), "
                      "(SELECT count(*) FROM variant v WHERE v.gene_id=g.id AND v.impact>=1) "
                      "FROM acmg a LEFT JOIN gene g ON g.symbol=a.gene AND g.biotype='protein_coding' ORDER BY a.category, a.gene", NULL);
    if (dbHasTable("trait")) g_hNotable = dbq("SELECT title, result, detail, genotype, about, source, conf, status, copies FROM trait WHERE cat='Health' ORDER BY ord", NULL);
    if (dbHasTable("sv"))
        g_hSv = dbq("SELECT chrom, pos, end, type, len, gt, genes, exons, support, depth_ratio, prec, filter FROM sv WHERE pass=1 AND exons>0 "
                    "AND type IN ('DEL','DUP','INV') ORDER BY support DESC, exons DESC LIMIT 300", NULL);
    if (dbHasTable("genedel")) g_hDel = dbq("SELECT gene, call, callable, depth, expected, note FROM genedel", NULL);
    rowsFree(&g_hRare);
    /* rare_var (step 16) adds the plain reading and the suspect flag; without it the candidates
       are still listed, with empty columns 14-16 */
    bool hasRv = g_hRareRv = dbHasTable("rare_var");
    char sql[2048];
    snprintf(sql, sizeof sql,
             "SELECT v.id, v.gene, v.hgvs_p, v.csq, v.zyg, v.af, v.am, g.loeuf, "
             "(SELECT group_concat(DISTINCT d.inheritance) FROM gene_disease d WHERE d.gene=v.gene_id), "
             "(SELECT d.disease FROM gene_disease d WHERE d.gene=v.gene_id LIMIT 1), v.pass, v.dp, v.ad_alt, v.ad_ref, %s "
             "FROM variant v JOIN gene g ON g.id=v.gene_id %s WHERE (v.clin=0 OR v.clin=3) AND (v.af IS NULL OR v.af < 0.001) AND v.pass=1 "
             "AND v.ad_alt * 5 >= (v.ad_ref + v.ad_alt) AND ("
             " (v.impact=3 AND g.loeuf < 0.5) OR (v.csq='missense' AND v.am > 0.564 AND v.gene_id IN (SELECT gene FROM gene_disease)))"
             " ORDER BY %s v.impact DESC, g.loeuf LIMIT 200",
             hasRv ? "coalesce(rv.plain, ''), coalesce(rv.change, ''), coalesce(rv.suspect, 0)" : "'', '', 0",
             hasRv ? "LEFT JOIN rare_var rv ON rv.vid=v.id" : "",
             hasRv ? "coalesce(rv.suspect, 0)," : "");
    g_hRare = dbq(sql, NULL);
}

/* one finding as a card; returns its height */
static float findingCard(Rectangle c, Rows *r, int i, bool draw) {
    const char *gene = rc(r, i, 1);
    const char *inh = rc(r, i, 13);
    bool rel;
    const char *meaning = inhMeaning(inh, (int)ri(r, i, 4), &rel);
    char cond[600];
    snprintf(cond, sizeof cond, "%s", rc(r, i, 7));
    for (char *p = cond; *p; p++) if (*p == '|') { p[0] = ';'; }
    float h = 18 + 28 + txtParaH(F_REG, TH->fsSmall, c.width - 40, 19, cond, 3) + 8 + txtParaH(F_REG, TH->fsSmall, c.width - 40, 19, meaning, 3) + 52;
    if (!draw) return h;
    c.height = h;
    gxCard(c);
    if (rel) drRound((Rectangle){ c.x, c.y + 14, 4, c.height - 28 }, 2, TH->bad);
    float x = c.x + 20, y = c.y + 18;
    int clin = (int)ri(r, i, 5);
    if (gxLink(gxIdN(20000 + i), x, y, TH->fsH3, F_SEMI, *gene ? gene : "intergenic", TH->accent, "Open the gene page"))
        { Rows *g = dbq("SELECT gene_id FROM variant WHERE id=?", rc(r, i, 0), NULL); if (g->nrow) navGene((int)ri(g, 0, 0)); rowsFree(&g); }
    float gx = x + (float)txtW(F_SEMI, TH->fsH3, *gene ? gene : "intergenic") + 10;
    txtIn(F_REG, TH->fs, (Rectangle){ gx, y - 1, c.width - (gx - c.x) - 330, 20 }, TH->text, 0,
          *rc(r, i, 2) ? rc(r, i, 2) : csqLabel(rc(r, i, 3)));
    float px = c.x + c.width - 20;
    const char *zs = zygName((int)ri(r, i, 4));
    px -= gxPillW(clinName(clin), TH->fsTiny); gxPill(px, y - 2, clinName(clin), clinColor(clin), TH->fsTiny); px -= 6;
    px -= gxPillW(zs, TH->fsTiny); gxPill(px, y - 2, zs, TH->purple, TH->fsTiny); px -= 6;
    if (*inh) { px -= gxPillW(inh, TH->fsTiny); gxPill(px, y - 2, inh, TH->accent2, TH->fsTiny); }
    y += 28;
    y += txtPara(F_REG, TH->fsSmall, x, y, c.width - 40, 19, TH->text, cond, 3) + 8;
    y += txtPara(F_REG, TH->fsSmall, x, y, c.width - 40, 19, rel ? TH->warn : TH->dim, meaning, 3) + 8;
    int st = (int)ri(r, i, 8);
    txt(F_REG, TH->fsTiny, x, y + 6, TH->faint, sstr("%s   ClinVar review %d/4   frequency %s   %s", fmtPos((int)ri(r, i, 11), ri(r, i, 12)), st,
                                                     fmtAf(rf(r, i, 9)), ri(r, i, 14) ? "passed filters" : "flagged call: check the reads"));
    if (gxButtonEx(gxIdN(21000 + i), (Rectangle){ c.x + c.width - 230, y, 100, 26 }, "Details", BTN_NORMAL, true, NULL)) {
        variantsShow(0, gene, ri(r, i, 0));
    }
    if (gxButtonEx(gxIdN(22000 + i), (Rectangle){ c.x + c.width - 122, y, 100, 26 }, "See reads", BTN_GHOST, true, "Open the position in the genome browser"))
        navVariant(ri(r, i, 0));
    if (gxBoxClick(gxIdN(20500 + i), c)) { lkBegin(sstr("%s %s", gene, rc(r, i, 2))); lkVariant(ri(r, i, 0)); lkGene(gene); lkOpen(); }
    return h;
}

static void tabPathogenic(float x, float *y, float W) {
    Rows *r = g_hPath;
    int nrel = 0;
    for (int i = 0; r && i < r->nrow; i++) { bool rel; inhMeaning(rc(r, i, 13), (int)ri(r, i, 4), &rel); nrel += rel; }
    char summary[400];
    snprintf(summary, sizeof summary, "%d variant%s you carry %s classed pathogenic or likely pathogenic in ClinVar. %s",
             r ? r->nrow : 0, r && r->nrow == 1 ? "" : "s", r && r->nrow == 1 ? "is" : "are",
             nrel ? sstr("%d may be relevant to your own health (marked in red); the others make you a carrier.", nrel)
                  : "None is expected to affect your own health on its own: each is a single copy in a recessive gene.");
    *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->text, summary, 3) + 16;
    if (dkChip() && r && r->nrow)
        *y += txtPara(F_SEMI, TH->fsSmall, x, *y, W, 19, TH->warn,
                      "These are chip calls. Chips can report rare disease variants that are not there (in Weedon et al. 2021, BMJ, most very rare "
                      "chip calls were not confirmed by sequencing), so a clinical test may make sense before acting on any of these.", 4) + 14;
    for (int i = 0; r && i < r->nrow; i++) {
        float h = findingCard((Rectangle){ x, *y, W, 0 }, r, i, false);
        findingCard((Rectangle){ x, *y, W, h }, r, i, true);
        *y += h + 12;
    }
    /* chip data: rare calls ClinVar lists as disease-causing, set aside by the import */
    if (dkChip() && g_hChipRare && g_hChipRare->nrow) {
        *y += 12;
        txt(F_SEMI, TH->fsH3, x, *y, TH->text, sstr("Unconfirmed chip calls (%d)", g_hChipRare->nrow));
        *y += 26;
        *y += txtPara(F_REG, TH->fsSmall, x, *y, W, 19, TH->dim,
                      "The chip reported these rare variants, which ClinVar lists as disease-causing. Chip calls this rare are usually wrong: in a "
                      "study of nearly 50,000 people most were not there when the same people were sequenced (Weedon et al. 2021, BMJ). So they "
                      "are not counted above and the rest of Strand Open treats these positions as untested; a clinical test may make sense before "
                      "reading anything into them.", 5) + 10;
        for (int i = 0; i < g_hChipRare->nrow; i++) {
            Rectangle row = { x, *y, W, 30 };
            bool h;
            if (gxBehave(gxIdN(23500 + i), row, &h, NULL)) lkVariantMenu(ri(g_hChipRare, i, 0));
            if (h) { drRound(row, 6, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
            txtIn(F_SEMI, TH->fsSmall, (Rectangle){ row.x + 10, row.y, 90, 30 }, TH->text, 0, rc(g_hChipRare, i, 1));
            txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 100, row.y, 200, 30 }, TH->dim, 0, *rc(g_hChipRare, i, 2) ? rc(g_hChipRare, i, 2) : csqLabel(rc(g_hChipRare, i, 3)));
            txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 310, row.y, 60, 30 }, TH->dim, 0, zygShort((int)ri(g_hChipRare, i, 4)));
            char cond[300];
            snprintf(cond, sizeof cond, "%s", rc(g_hChipRare, i, 5));
            for (char *p = cond; *p; p++) if (*p == '|') *p = ';';
            txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 380, row.y, W - 390, 30 }, TH->faint, 0, cond);
            *y += 30;
        }
    }
    if (g_hConf && g_hConf->nrow) {
        *y += 12;
        txt(F_SEMI, TH->fsH3, x, *y, TH->text, "Conflicting interpretations that include pathogenic");
        *y += 26;
        *y += txtPara(F_REG, TH->fsSmall, x, *y, W, 19, TH->dim, "Rare protein-changing variants that some ClinVar submitters call pathogenic and others do not. Usually these end up benign or uncertain.", 3) + 10;
        for (int i = 0; i < g_hConf->nrow; i++) {
            Rectangle row = { x, *y, W, 30 };
            bool h;
            if (gxBehave(gxIdN(23000 + i), row, &h, NULL)) lkVariantMenu(ri(g_hConf, i, 0));
            if (h) { drRound(row, 6, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
            txtIn(F_SEMI, TH->fsSmall, (Rectangle){ row.x + 10, row.y, 90, 30 }, TH->text, 0, rc(g_hConf, i, 1));
            txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 100, row.y, 200, 30 }, TH->dim, 0, *rc(g_hConf, i, 2) ? rc(g_hConf, i, 2) : csqLabel(rc(g_hConf, i, 3)));
            txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 310, row.y, 60, 30 }, TH->dim, 0, zygShort((int)ri(g_hConf, i, 4)));
            char cond[300];
            snprintf(cond, sizeof cond, "%s", rc(g_hConf, i, 6));
            for (char *p = cond; *p; p++) if (*p == '|') *p = ';';
            txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 380, row.y, W - 480, 30 }, TH->faint, 0, cond);
            txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + W - 100, row.y, 90, 30 }, TH->dim, 2, fmtAf(rf(g_hConf, i, 7)));
            *y += 30;
        }
    }
}
static void tabRisk(float x, float *y, float W) {
    Rows *r = g_hRisk;
    *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->dim, "Variants ClinVar lists as risk factors, protective or associated with a condition. Most raise or lower risk modestly; "
                  "they are common and not diagnostic.", 3) + 14;
    for (int i = 0; r && i < r->nrow; i++) {
        Rectangle row = { x, *y, W, 46 };
        bool h;
        if (gxBehave(gxIdN(24000 + i), row, &h, NULL)) lkVariantMenu(ri(r, i, 0));
        if (h) { drRound(row, 6, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
        int clin = (int)ri(r, i, 5);
        txtIn(F_SEMI, TH->fs, (Rectangle){ row.x + 12, row.y + 4, 110, 20 }, TH->text, 0, *rc(r, i, 1) ? rc(r, i, 1) : "-");
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 124, row.y + 4, 220, 20 }, TH->dim, 0, *rc(r, i, 2) ? rc(r, i, 2) : ri(r, i, 10) ? sstr("rs%lld", ri(r, i, 10)) : csqLabel(rc(r, i, 3)));
        gxPill(row.x + 350, row.y + 3, zygShort((int)ri(r, i, 4)), TH->purple, TH->fsTiny);
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + W - 200, row.y + 4, 190, 20 }, clinColor(clin), 2, rc(r, i, 6));
        char cond[400];
        snprintf(cond, sizeof cond, "%s", rc(r, i, 7));
        for (char *p = cond; *p; p++) if (*p == '|') *p = ';';
        txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ row.x + 12, row.y + 24, W - 120, 18 }, TH->faint, 0, cond);
        txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ row.x + W - 110, row.y + 24, 100, 18 }, TH->faint, 2, fmtAf(rf(r, i, 9)));
        *y += 48;
    }
}
static void tabAcmg(float x, float *y, float W) {
    Rows *r = g_hAcmg;
    if (!r) { txt(F_REG, TH->fs, x, *y, TH->faint, "The ACMG step did not run for this person: re-import them from the People page."); return; }
    int found = 0;
    for (int i = 0; i < r->nrow; i++) found += ri(r, i, 6) > 0;
    *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->text,
                  sstr("ACMG SF v3.3 names %d genes where a clinical genome reports a pathogenic variant even when nobody asked, because acting on it "
                       "can prevent serious disease (cancers, heart conditions, a few metabolic ones). %s", r->nrow,
                       found ? sstr("You have ClinVar pathogenic or likely pathogenic variants in %d of them (below, in red).", found)
                             : "You have no ClinVar pathogenic or likely pathogenic variant in any of them."), 4) + 8;
    *y += txtPara(F_REG, TH->fsSmall, x, *y, W, 19, TH->dim, "Rare protein-changing variants in these genes are counted too: most are harmless; "
                  "a few are unrecognised risks. 'Coding' counts every coding variant, common ones included.", 3) + 14;
    const char *cat = "";
    float colW = (W - 24) / 3;
    int col = 0;
    float rowY = *y;
    for (int i = 0; i < r->nrow; i++) {
        if (strcmp(cat, rc(r, i, 1))) {
            if (col) { rowY += 64; col = 0; }
            cat = rc(r, i, 1);
            rowY += 8;
            gxSectionTitle(x, rowY, cat);
            rowY += 22;
        }
        Rectangle c = { x + (colW + 12) * (float)col, rowY, colW, 56 };
        bool h;
        if (gxBehave(gxIdN(25000 + i), c, &h, NULL)) lkGeneMenu(rc(r, i, 0));
        drRound(c, 8, h ? TH->hover : TH->card);
        drRoundLine(c, 8, ri(r, i, 6) ? TH->bad : TH->line);
        if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
        txt(F_SEMI, TH->fs, c.x + 12, c.y + 8, ri(r, i, 6) ? TH->bad : TH->text, rc(r, i, 0));
        txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 12, c.y + 30, c.width - 24, 18 }, TH->faint, 0, rc(r, i, 2));
        const char *st = ri(r, i, 6) ? sstr("%lld pathogenic", ri(r, i, 6)) : ri(r, i, 7) ? sstr("%lld rare coding", ri(r, i, 7)) : sstr("%lld coding", ri(r, i, 8));
        txtIn(F_SEMI, TH->fsTiny, (Rectangle){ c.x + c.width - 130, c.y + 8, 118, 18 }, ri(r, i, 6) ? TH->bad : ri(r, i, 7) ? TH->warn : TH->faint, 2, st);
        if (++col == 3) { col = 0; rowY += 64; }
    }
    *y = rowY + 70;
}
static void tabNotable(float x, float *y, float W) {
    Rows *r = g_hNotable;
    if (!r) { txt(F_REG, TH->fs, x, *y, TH->faint, "The traits step did not run for this person: re-import them from the People page."); return; }
    *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->dim, "Well-studied variants that most genetic services report. Each is read from your own calls at the exact position.", 2) + 14;
    gxSegmented(gxId("notable found"), (Rectangle){ x, *y, 360, 32 }, g_foundOpts, 3, &g_hFound);
    *y += 46;
    float colW = (W - 16) / 2;
    float colY[2] = { *y, *y };
    for (int i = 0; i < r->nrow; i++) {
        if (g_hFound && traitFound(r, i, 8, 1) != (g_hFound == 1)) continue;
        int col = colY[0] <= colY[1] ? 0 : 1;
        float cx = x + (colW + 16) * (float)col;
        float h = 60 + txtParaH(F_REG, TH->fsSmall, colW - 36, 19, rc(r, i, 2), 6) + 10 + txtParaH(F_REG, TH->fsTiny + 1, colW - 36, 17, rc(r, i, 4), 6) + 34;
        Rectangle c = { cx, colY[col], colW, h };
        gxCard(c);
        txt(F_SEMI, TH->fs, c.x + 18, c.y + 16, TH->dim, rc(r, i, 0));
        txt(F_SEMI, TH->fsH2, c.x + 18, c.y + 36, TH->text, rc(r, i, 1));
        float ty = c.y + 66;
        ty += txtPara(F_REG, TH->fsSmall, c.x + 18, ty, colW - 36, 19, TH->text, rc(r, i, 2), 6) + 10;
        ty += txtPara(F_REG, TH->fsTiny + 1, c.x + 18, ty, colW - 36, 17, TH->faint, rc(r, i, 4), 6) + 6;
        txtIn(F_MONO, 11, (Rectangle){ c.x + 18, ty, colW - 36, 18 }, TH->faint, 0, rc(r, i, 3));
        if (gxBoxClick(gxIdN(29500 + i), c)) { lkBegin(rc(r, i, 0)); lkTrait(rc(r, i, 0), rc(r, i, 3)); lkOpen(); }
        colY[col] += h + 14;
    }
    *y = fmaxf(colY[0], colY[1]);
}
static void tabRare(float x, float *y, float W) {
    if (dkChip()) {
        *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->dim,
                      "Rare variants that could disrupt a gene need sequencing data to find. A genotyping chip reads common, well-known positions, and "
                      "the rare calls it does make are mostly wrong when checked by sequencing (Weedon et al. 2021, BMJ), so none are listed here. "
                      "Known disease variants the chip did read appear under ClinVar.", 5) + 14;
        return;
    }
    *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->dim,
                  "Rare variants nobody has classified that predictions suggest could disrupt a gene (calls with under a fifth of reads showing the variant are left out): protein-truncating variants in genes that "
                  "healthy people rarely carry broken, and missense changes AlphaMissense rates likely pathogenic in genes linked to a condition. "
                  "Most such variants turn out harmless; these are leads to check (look at the reads first), not findings.", 4) + 14;
    if (!g_hRareRv)
        *y += txtPara(F_REG, TH->fsSmall, x, *y, W, 19, TH->warn,
                      "The check for doubtful calls has not run, so none are flagged here, and the plain readings are missing. "
                      "Re-import this person from the People page to add them.", 4) + 12;
    Rows *r = g_hRare;
    if (!r || !r->nrow) { txt(F_REG, TH->fs, x, *y, TH->faint, "None found."); *y += 30; return; }
    for (int i = 0; i < r->nrow; i++) {
        const char *plainTxt = rc(r, i, 14), *chg = rc(r, i, 15);
        bool sus = ri(r, i, 16) != 0;
        float h = 48 + (*plainTxt ? 19 : 0) + (*chg ? txtParaH(F_REG, TH->fsSmall - 1, W - 40, 18, chg, 2) : 0) + (sus ? 18 : 0) + 6;
        Rectangle row = { x, *y, W, h };
        bool hov;
        if (gxBehave(gxIdN(28000 + i), row, &hov, NULL)) lkVariantMenu(ri(r, i, 0));
        if (hov) { drRound(row, 6, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
        const char *csq = rc(r, i, 3);
        bool lof = strcmp(csq, "missense") != 0;
        DrawCircleV((Vector2){ row.x + 12, row.y + 15 }, 4, sus ? TH->faint : lof ? TH->bad : TH->warn);
        txtIn(F_SEMI, TH->fs, (Rectangle){ row.x + 24, row.y + 4, 100, 20 }, TH->text, 0, rc(r, i, 1));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 128, row.y + 4, 260, 20 }, TH->dim, 0, *rc(r, i, 2) ? rc(r, i, 2) : csqLabel(csq));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 392, row.y + 4, 140, 20 }, TH->dim, 0, csqLabel(csq));
        gxPill(row.x + 540, row.y + 3, zygShort((int)ri(r, i, 4)), TH->purple, TH->fsTiny);
        const char *score = lof ? sstr("LOEUF %.2f", rf(r, i, 7)) : sstr("AlphaMissense %.2f", rf(r, i, 6));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + W - 330, row.y + 4, 160, 20 }, TH->text, 2, score);
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + W - 160, row.y + 4, 150, 20 }, rnull(r, i, 5) ? TH->warn : TH->dim, 2, fmtAf(rf(r, i, 5)));
        float ty = row.y + 25;
        if (*plainTxt) { txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 24, ty, W - 40, 18 }, TH->text, 0, plainTxt); ty += 19; }
        if (*chg) ty += txtPara(F_REG, TH->fsSmall - 1, row.x + 24, ty, W - 40, 18, TH->dim, chg, 2);
        if (sus) { txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ row.x + 24, ty, W - 40, 18 }, TH->warn, 0,
                         "Suspect call: hard-to-map DNA, a repeat, a cluster of odd calls or unusual read counts. Check the reads before trusting it."); ty += 18; }
        int adr = (int)ri(r, i, 13), ada = (int)ri(r, i, 12);
        txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ row.x + 24, ty, W - 40, 18 }, TH->faint, 0,
              sstr("%s%s   reads %d ref / %d alt", *rc(r, i, 9) ? rc(r, i, 9) : "no condition recorded for this gene",
                   *rc(r, i, 8) ? sstr(" (%s)", rc(r, i, 8)) : "", adr, ada));
        *y += h + 4;
    }
    *y += 10;
    if (gxButtonEx(gxId("rare page"), (Rectangle){ x, *y, 260, 32 }, "Rare and remarkable variants", BTN_NORMAL, true, NULL)) navPage(PG_RARE);
    *y += 40;
}
static void tabSv(float x, float *y, float W) {
    if (!dkGvcf() && !dbHasTable("sv")) {
        *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->dim, dkChip()
                      ? "Larger changes (deletions and duplications of whole stretches of DNA) need sequencing data: a chip reads single positions."
                      : "Larger changes (deletions and duplications of whole stretches of DNA) come from a separate structural-variant file, which this "
                        "person's data does not include.", 3) + 14;
        return;
    }
    *y += txtPara(F_REG, TH->fs, x, *y, W, 21, TH->dim, "Larger changes: deletions and duplications of whole stretches of DNA that touch protein-coding exons, "
                  "and common whole-gene deletions read from coverage. Short-read callers report many large structural calls that are not real, "
                  "so each is checked against read depth: 'supported' means the depth agrees.", 4) + 14;
    if (g_hDel && g_hDel->nrow) {
        txt(F_SEMI, TH->fsH3, x, *y, TH->text, "Common gene deletions");
        *y += 30;
        for (int i = 0; i < g_hDel->nrow; i++) {
            Rectangle c = { x, *y, W, 54 };
            gxCard(c);
            txt(F_SEMI, TH->fs, c.x + 16, c.y + 10, TH->text, rc(g_hDel, i, 0));
            const char *call = rc(g_hDel, i, 1);
            txt(F_SEMI, TH->fs, c.x + 110, c.y + 10, strstr(call, "both") ? TH->warn : TH->ok, call);
            txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ c.x + 16, c.y + 30, W - 32, 18 }, TH->dim, 0, rc(g_hDel, i, 5));
            txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + W - 300, c.y + 10, 284, 18 }, TH->faint, 2,
                  sstr("%.0f%% of coding bases called, depth %.1f", rf(g_hDel, i, 2) * 100, rf(g_hDel, i, 3)));
            *y += 62;
        }
        *y += 10;
    }
    txt(F_SEMI, TH->fsH3, x, *y, TH->text, "Structural variants touching coding exons");
    *y += 30;
    Rows *r = g_hSv;
    for (int i = 0; r && i < r->nrow && i < 120; i++) {
        Rectangle row = { x, *y, W, 44 };
        bool h;
        int chrom = (int)ri(r, i, 0);
        if (gxBehave(gxIdN(26000 + i), row, &h, NULL)) {
            lkBegin(sstr("%s %s, %s", rc(r, i, 3), fmtPos(chrom, ri(r, i, 1)), fmtBp((double)ri(r, i, 4))));
            lkRegion(chrom, ri(r, i, 1), ri(r, i, 2));
            char g1[40] = "";
            sscanf(rc(r, i, 6), "%39[A-Za-z0-9-]", g1);
            if (g1[0]) lkGene(g1);
            lkOpen();
        }
        if (h) { drRound(row, 6, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
        bool sup = ri(r, i, 8) == 1, unk = rnull(r, i, 8);
        txtIn(F_SEMI, TH->fsSmall, (Rectangle){ row.x + 10, row.y + 4, 60, 18 }, TH->text, 0, rc(r, i, 3));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 70, row.y + 4, 220, 18 }, TH->dim, 0, sstr("%s  %s", fmtPos(chrom, ri(r, i, 1)), fmtBp((double)ri(r, i, 4))));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 300, row.y + 4, 60, 18 }, TH->dim, 0, rc(r, i, 5));
        txtIn(F_SEMI, TH->fsSmall, (Rectangle){ row.x + W - 230, row.y + 4, 220, 18 }, unk ? TH->dim : sup ? TH->warn : TH->faint, 2,
              unk ? (ri(r, i, 10) ? "precise breakpoints" : "imprecise") : sup ? sstr("supported (depth %.2fx)", rf(r, i, 9)) : sstr("not supported (depth %.2fx)", rf(r, i, 9)));
        txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ row.x + 10, row.y + 24, W - 20, 18 }, TH->faint, 0, sstr("%lld coding exons: %s", ri(r, i, 7), rc(r, i, 6)));
        *y += 46;
    }
}
static void gwasCell(int row, int col, Rectangle c, bool sel, void *u) {
    (void)sel; (void)u;
    Rows *r = g_hGwas;
    switch (col) {
    case 0: txtIn(F_REG, TH->fsSmall, c, TH->text, 0, rc(r, row, 6)); break;
    case 1: txtIn(F_MONO, 12, c, TH->dim, 0, sstr("%s-%s", rc(r, row, 0), rc(r, row, 3))); break;
    case 2: { int cp = (int)ri(r, row, 4); txtIn(F_SEMI, TH->fsSmall, c, cp == 0 ? TH->faint : cp == 1 ? TH->warn : TH->bad, 1, sstr("%d of %lld", cp, ri(r, row, 5))); break; }
    case 3: txtIn(F_REG, TH->fsSmall, c, TH->dim, 2, rc(r, row, 10)); break;
    case 4: txtIn(F_REG, TH->fsSmall, c, TH->dim, 2, sstr("1e-%.0f", rf(r, row, 9))); break;
    case 5: txtIn(F_REG, TH->fsSmall, c, TH->dim, 0, rc(r, row, 16)); break;
    case 6: txtIn(F_REG, TH->fsSmall, c, TH->faint, 0, sstr("%s %s", rc(r, row, 14), rc(r, row, 15))); break;
    }
}
static const GxCol g_gwasCols[] = { { "TRAIT", -3, 0, false, NULL }, { "VARIANT", 150, 0, false, NULL }, { "YOUR COPIES", 90, 1, false, "Copies you carry of the risk (or effect) allele" },
                                    { "OR / BETA", 80, 2, false, NULL }, { "P", 64, 2, false, NULL }, { "GENE", -1, 0, false, NULL }, { "STUDY", 150, 0, false, NULL } };
static const char *const g_gwasMine[] = { "All", "You carry it", "Two copies" };
static void tabGwas(Rectangle r) {
    Rows *t;
    if ((t = dbTake(SLOT_GWAS))) { rowsFree(&g_hGwas); g_hGwas = t; g_hGwasT.scroll = (GxScroll){ 0 }; }
    if (g_hGwasDirty) {
        g_hGwasDirty = false;
        char where[200];
        snprintf(where, sizeof where, "%s", g_hGwasMine == 1 ? "g.copies>0" : g_hGwasMine == 2 ? "g.copies>=2" : "1");
        if (g_hGwasQ[0]) {
            char fts[200];
            snprintf(fts, sizeof fts, "\"%s\"*", g_hGwasQ);
            for (char *c = fts + 1; *c; c++) if (*c == '"' && c[1] != '*') *c = ' ';
            const char *p[] = { fts };
            dbAsync(SLOT_GWAS, sstr("SELECT g.* FROM gwas g WHERE g.rowid IN (SELECT rowid FROM gwas_fts WHERE gwas_fts MATCH ?1) AND %s ORDER BY g.mlog DESC LIMIT 3000", where), 1, p);
        } else dbAsync(SLOT_GWAS, sstr("SELECT g.* FROM gwas g WHERE %s ORDER BY g.mlog DESC LIMIT 3000", where), 0, NULL);
    }
    float x = r.x, y = r.y;
    txtPara(F_REG, TH->fsSmall, x, y, r.width, 19, TH->dim, "Genome-wide significant findings from the GWAS Catalog at positions where your genotype is known. "
            "Each is a small average effect seen across thousands of people; on its own it predicts little for one person. Use the polygenic scores page for combined estimates.", 3);
    y += 48;
    if (gxTextField(gxId("gwas q"), (Rectangle){ x, y, 320, 32 }, g_hGwasQ, sizeof g_hGwasQ, "Search traits or genes (e.g. migraine, height)", true)) g_hGwasDirty = true;
    if (gxSegmented(gxId("gwas mine"), (Rectangle){ x + 332, y, 330, 32 }, g_gwasMine, 3, &g_hGwasMine)) g_hGwasDirty = true;
    txt(F_REG, TH->fsSmall, x + 680, y + 8, TH->faint, dbBusy(SLOT_GWAS) ? "searching..." : sstr("%s associations", fmtInt(g_hGwas ? g_hGwas->nrow : 0)));
    y += 44;
    int n = g_hGwas ? g_hGwas->nrow : 0;
    int clicked = gxTable(gxId("gwas t"), &g_hGwasT, (Rectangle){ x - 8, y, r.width + 16, r.y + r.height - y }, g_gwasCols, 7, n, 30, gwasCell, NULL, NULL);
    if (clicked >= 0 && *rc(g_hGwas, clicked, 13)) openUrl(sstr("https://pubmed.ncbi.nlm.nih.gov/%s/", rc(g_hGwas, clicked, 13)));
}

static const char *const g_hTabs[] = { "Pathogenic", "Risk factors", "ACMG genes", "Well-known variants", "Rare, possibly damaging", "Structural", "Research (GWAS)" };
static void pageHealth(Rectangle r) {
    if (!g_hAsked) hAsk();
    float x = r.x + 32, W = fminf(1240, r.width - 64);
    txt(F_SEMI, TH->fsH1, x, r.y + 22, TH->text, "Clinical findings");
    if (gxTabs(gxId("health tabs"), x - 14, r.y + 58, 40, g_hTabs, 7, &g_hTab)) g_hScroll = (GxScroll){ 0 };
    drHLine(r.x, r.y + 98, r.width, TH->line);
    Rectangle body = { r.x, r.y + 99, r.width, r.height - 99 };
    if (g_hTab == 6) { tabGwas((Rectangle){ x, body.y + 18, W, body.height - 30 }); return; }
    float top = gxScrollBegin(&g_hScroll, body);
    float y = top + 22;
    switch (g_hTab) {
    case 0: tabPathogenic(x, &y, W); break;
    case 1: tabRisk(x, &y, W); break;
    case 2: tabAcmg(x, &y, W); break;
    case 3: tabNotable(x, &y, W); break;
    case 4: tabRare(x, &y, W); break;
    case 5: tabSv(x, &y, W); break;
    }
    y += 20;
    txtPara(F_REG, TH->fsTiny + 1, x, y, W, 17, TH->faint, "Strand Open is for personal exploration, not diagnosis. Sequencing calls can be wrong, ClinVar classifications change, "
            "and most conditions depend on much more than one variant. Confirm anything important with a clinical-grade test and a doctor or genetic counsellor.", 3);
    y += 60;
    gxScrollEnd(&g_hScroll, body, y - top);
}

/* ---- carrier status ---- */
static Rows *g_cRows, *g_cScreen;
static GxScroll g_cScroll;
static void pageCarrier(Rectangle r) {
    if (!g_cRows) {
        g_cRows = dbq("SELECT v.id, v.gene, v.hgvs_p, v.csq, v.zyg, v.clin, v.clin_dn, v.clin_stars, v.af, "
                      "(SELECT group_concat(DISTINCT d.inheritance) FROM gene_disease d WHERE d.gene=v.gene_id), v.chrom "
                      "FROM variant v WHERE v.clin IN (5,6) ORDER BY v.gene", NULL);
        g_cScreen = dbq("SELECT count(DISTINCT gene) FROM gene_disease WHERE inheritance LIKE '%AR%'", NULL);
    }
    float top = gxScrollBegin(&g_cScroll, r);
    float x = r.x + 32, W = fminf(1100, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Carrier status",
                         "Recessive conditions need a broken copy of a gene from each parent. Carriers have one broken copy and one working one, "
                         "and are usually healthy, but two carriers of the same condition have a one in four chance of an affected child each pregnancy.");
    int n = 0;
    for (int i = 0; i < g_cRows->nrow; i++) {
        const char *inh = rc(g_cRows, i, 9);
        if (!strstr(inh, "AR") && !strstr(inh, "XLR") && !strstr(inh, "XL")) continue;
        int zyg = (int)ri(g_cRows, i, 4);
        float h = 104;
        Rectangle c = { x, y, W, h };
        gxCard(c);
        drRound((Rectangle){ c.x, c.y + 14, 4, h - 28 }, 2, zyg == 1 ? TH->warn : TH->bad);
        txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, rc(g_cRows, i, 1));
        txt(F_REG, TH->fs, c.x + 20 + (float)txtW(F_SEMI, TH->fsH3, rc(g_cRows, i, 1)) + 10, c.y + 17, TH->dim, *rc(g_cRows, i, 2) ? rc(g_cRows, i, 2) : csqLabel(rc(g_cRows, i, 3)));
        const char *status = zyg == 1 ? "Carrier (one copy)" : zyg == 3 ? "Hemizygous: may be affected" : "Two copies: may be affected";
        gxPill(c.x + c.width - gxPillW(status, TH->fsTiny) - 18, c.y + 14, status, zyg == 1 ? TH->warn : TH->bad, TH->fsTiny);
        char cond[400];
        snprintf(cond, sizeof cond, "%s", rc(g_cRows, i, 6));
        for (char *p = cond; *p; p++) if (*p == '|') *p = ';';
        txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 20, c.y + 44, W - 40, 20 }, TH->text, 0, cond);
        txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ c.x + 20, c.y + 68, W - 260, 20 }, TH->faint, 0,
              sstr("Inheritance %s. Seen in %s of chromosomes. ClinVar review %lld/4.", inh, fmtAf(rf(g_cRows, i, 8)), ri(g_cRows, i, 7)));
        if (gxButtonEx(gxIdN(27000 + i), (Rectangle){ c.x + c.width - 120, c.y + 62, 100, 28 }, "Details", BTN_NORMAL, true, NULL)) navVariant(ri(g_cRows, i, 0));
        if (gxBoxClick(gxIdN(27500 + i), c)) lkVariantMenu(ri(g_cRows, i, 0));
        y += h + 12;
        n++;
    }
    if (!n) {
        Rectangle c = { x, y, W, 90 };
        gxCard(c);
        txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 18, TH->ok, "No carrier findings");
        txtPara(F_REG, TH->fsSmall, c.x + 20, c.y + 46, W - 40, 19, TH->dim, "No ClinVar pathogenic or likely pathogenic variant in a recessive gene.", 2);
        y += 102;
    }
    y += 10;
    long long screened = g_cScreen && g_cScreen->nrow ? ri(g_cScreen, 0, 0) : 0;
    y += txtPara(F_REG, TH->fsSmall, x, y, W, 19, TH->dim,
                 sstr("Screened: every gene with a recessive condition in HPO (%s genes), for variants ClinVar classes pathogenic or likely pathogenic. "
                      "Not covered: variants ClinVar has not classified, copy-number changes such as SMA (SMN1) and the alpha-thalassaemia deletions, "
                      "and repeat expansions such as fragile X. A clinical carrier screen covers those.", fmtInt(screened)), 5);
    y += 40;
    gxScrollEnd(&g_cScroll, r, y - top);
}
