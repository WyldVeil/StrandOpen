/* ========================================================================
   src/page_variants.c - the variant explorer and the variant detail panel.

   The explorer filters and sorts all 4.6 million calls. One query (on a
   worker) returns the ids of the matching variants in the chosen order, up
   to 250,000 of them; the table then fetches the visible rows by id in
   pages of 100 (a primary-key lookup, well under a millisecond), so
   scrolling never re-runs the sort.

   The detail panel is shared with the genome browser and the gene page.
   ======================================================================== */

/* ---- the detail panel ---- */
typedef struct { long long id; Rows *v; Rows *g, *gw, *pk; GxScroll scroll; } VarPanel;
static VarPanel g_vp = { .id = -1 };
static const char *const VCOLS = "v.id, v.chrom, v.pos, v.ref, v.alt, v.rsid, v.kind, v.zyg, v.filter, v.pass, v.qual, v.dp, "
    "v.ad_ref, v.ad_alt, v.gq, v.gene, v.gene_id, v.csq, v.impact, v.tx, v.aa, v.hgvs_p, v.dna, v.clin, v.clin_sig, v.clin_dn, "
    "v.clin_stars, v.clin_id, v.af, v.af_eur, v.af_afr, v.af_eas, v.af_sas, v.af_amr, v.am, v.am_class";
enum { V_ID, V_CHROM, V_POS, V_REF, V_ALT, V_RSID, V_KIND, V_ZYG, V_FILTER, V_PASS, V_QUAL, V_DP, V_ADR, V_ADA, V_GQ,
       V_GENE, V_GENEID, V_CSQ, V_IMPACT, V_TX, V_AA, V_HGVSP, V_DNA, V_CLIN, V_CLINSIG, V_CLINDN, V_CLINSTARS, V_CLINID,
       V_AF, V_AFEUR, V_AFAFR, V_AFEAS, V_AFSAS, V_AFAMR, V_AM, V_AMCLASS };

static void vpLoad(long long id) {
    if (g_vp.id == id && g_vp.v) return;
    rowsFree(&g_vp.v);
    rowsFree(&g_vp.g);
    rowsFree(&g_vp.gw);
    rowsFree(&g_vp.pk);
    g_vp.id = id;
    g_vp.scroll = (GxScroll){ 0 };
    char sql[1200];
    snprintf(sql, sizeof sql, "SELECT %s FROM variant v WHERE v.id=?", VCOLS);
    g_vp.v = dbq(sql, sstr("%lld", id), NULL);
    if (g_vp.v->nrow && !rnull(g_vp.v, 0, V_GENEID))
        g_vp.g = dbq("SELECT id, symbol, name, pli, loeuf, mis_z, function, strand FROM gene WHERE id=?", rc(g_vp.v, 0, V_GENEID), NULL);
    if (g_vp.v->nrow && ri(g_vp.v, 0, V_RSID)) {
        const char *rs = sstr("rs%lld", ri(g_vp.v, 0, V_RSID));
        if (dbHasTable("gwas")) g_vp.gw = dbq("SELECT trait, risk, copies, effect, mlog, author, year, pmid FROM gwas WHERE rsid=? ORDER BY mlog DESC LIMIT 12", rs, NULL);
        if (dbHasTable("pgkb")) g_vp.pk = dbq("SELECT drugs, level, genotype, text FROM pgkb WHERE rsid=? ORDER BY level LIMIT 8", rs, NULL);
    }
}
static float vpSection(float x, float y, float w, const char *title) {
    drHLine(x, y, w, TH->line);
    txt(F_SEMI, TH->fsTiny, x, y + 12, TH->faint, title);
    return y + 34;
}
static float vpRow(float x, float y, float w, const char *k, const char *v, Color vc) {
    txt(F_REG, TH->fsSmall, x, y, TH->dim, k);
    float vx = x + 128;
    float h = txtPara(F_REG, TH->fsSmall, vx, y, w - 128, 18, vc, v, 6);
    return y + fmaxf(22, h + 4);
}
static void afBar(float x, float y, float w, const char *label, double af, Color c) {
    txt(F_REG, TH->fsSmall, x, y, TH->dim, label);
    Rectangle b = { x + 128, y + 4, w - 128 - 80, 9 };
    /* log scale from 1 in 10,000 to 1 */
    float v = isnan(af) || af <= 0 ? 0 : clamp01((float)((log10(af) + 4) / 4));
    gxMeter(b, v, c);
    txtIn(F_SEMI, TH->fsSmall, (Rectangle){ b.x + b.width + 6, y - 2, 74, 20 }, TH->text, 2, fmtAf(af));
}
static float g_vpH;
static float variantPanelH(void) { return g_vpH; }
static void variantPanel(Rectangle r, long long varId, bool showBrowse) {
    vpLoad(varId);
    Rows *v = g_vp.v;
    if (!v || !v->nrow) { txt(F_REG, TH->fs, r.x + 20, r.y + 20, TH->faint, "Variant not found."); return; }
    float top = gxScrollBegin(&g_vp.scroll, r);
    float x = r.x + 20, w = r.width - 40, y = top + 18;
    int chrom = (int)ri(v, 0, V_CHROM);
    long long pos = ri(v, 0, V_POS);
    const char *gene = rc(v, 0, V_GENE);
    /* heading */
    if (*gene) {
        if (gxLink(gxId("vp gene"), x, y, TH->fsH2, F_SEMI, gene, TH->accent, "Open the gene page"))
            navGene((int)ri(v, 0, V_GENEID));
        float gx = x + (float)txtW(F_SEMI, TH->fsH2, gene) + 10;
        txtIn(F_SEMI, TH->fsH3, (Rectangle){ gx, y, w - (gx - x) - 34, (float)TH->fsH2 }, TH->text, 0, rc(v, 0, V_HGVSP));
    } else {
        txt(F_SEMI, TH->fsH2, x, y, TH->text, "Intergenic variant");
        /* the nearest protein-coding gene, and which side of it */
        static long long nearFor = -1;
        static char nearTxt[160];
        static int nearId;
        if (nearFor != varId) {
            nearFor = varId;
            nearTxt[0] = 0;
            Rows *ng = dbq("SELECT id, symbol, start, end, strand FROM gene WHERE chrom=?1 AND biotype='protein_coding' "
                           "ORDER BY min(abs(start-?2), abs(end-?2)) LIMIT 1", sstr("%d", chrom), sstr("%lld", pos), NULL);
            if (ng->nrow) {
                long long gs = ri(ng, 0, 2), ge = ri(ng, 0, 3);
                int str = (int)ri(ng, 0, 4);
                long long d = pos < gs ? gs - pos : pos - ge;
                bool before = pos < gs;
                const char *side = (before == (str > 0)) ? "upstream of" : "downstream of";
                snprintf(nearTxt, sizeof nearTxt, "%s %s", fmtBp((double)d), side);
                nearId = (int)ri(ng, 0, 0);
                snprintf(nearTxt + strlen(nearTxt), sizeof nearTxt - strlen(nearTxt), "|%s", rc(ng, 0, 1));
            }
            rowsFree(&ng);
        }
        if (nearTxt[0]) {
            char *bar = strchr(nearTxt, '|');
            if (bar) {
                *bar = 0;
                float tx = x + (float)txtW(F_SEMI, TH->fsH2, "Intergenic variant") + 12;
                txt(F_REG, TH->fsSmall, tx, y + 6, TH->dim, nearTxt);
                float lx = tx + (float)txtW(F_REG, TH->fsSmall, nearTxt) + 5;
                if (gxLink(gxId("vp near"), lx, y + 6, TH->fsSmall, F_SEMI, bar + 1, TH->accent, "Open the gene page")) navGene(nearId);
                *bar = '|';
            }
        }
    }
    y += (float)TH->fsH2 + 10;
    char allele[200];
    snprintf(allele, sizeof allele, "%s  %.40s > %.40s", fmtPos(chrom, pos), rc(v, 0, V_REF), rc(v, 0, V_ALT));
    txt(F_MONO, 12, x, y, TH->dim, allele);
    y += 22;
    if (ri(v, 0, V_RSID)) {
        const char *rs = sstr("rs%lld", ri(v, 0, V_RSID));
        if (gxLink(gxId("vp rs"), x, y, TH->fsSmall, F_SEMI, rs, TH->accent2, "Open this rsID in dbSNP (your web browser)"))
            openUrl(sstr("https://www.ncbi.nlm.nih.gov/snp/%s", rs));
        y += 22;
    }
    float px = x;
    int zyg = (int)ri(v, 0, V_ZYG), imp = (int)ri(v, 0, V_IMPACT), clin = (int)ri(v, 0, V_CLIN);
    px += gxPill(px, y, zygName(zyg), TH->purple, TH->fsTiny) + 6;
    px += gxPill(px, y, csqLabel(rc(v, 0, V_CSQ)), impactColor(imp), TH->fsTiny) + 6;
    if (clin) px += gxPill(px, y, clinName(clin), clinColor(clin), TH->fsTiny) + 6;
    if (!ri(v, 0, V_PASS)) gxPill(px, y, "Filtered", TH->warn, TH->fsTiny);
    y += 34;

    /* what it means, in a sentence */
    {
        char say[600];
        const char *copies = zyg == 2 ? "both copies" : zyg == 3 ? "your single copy" : "one of your two copies";
        double af = rf(v, 0, V_AF);
        const char *rar = isnan(af) ? "not seen among the 2,590 people of 1000 Genomes" :
                          af < 0.001 ? "very rare" : af < 0.01 ? "rare" : af < 0.05 ? "uncommon" : af < 0.5 ? "common" : "the more common allele";
        snprintf(say, sizeof say, "You carry this %s in %s. It is %s%s.",
                 imp >= 2 ? "protein-changing variant" : imp == 1 ? "coding variant" : "variant", copies, rar,
                 isnan(af) ? "" : sstr(" (%s of chromosomes)", fmtAf(af)));
        y += txtPara(F_REG, TH->fsSmall, x, y, w, 19, TH->text, say, 4) + 10;
    }

    y = vpSection(x, y, w, "CALL QUALITY");
    {
        int adr = (int)ri(v, 0, V_ADR), ada = (int)ri(v, 0, V_ADA);
        y = vpRow(x, y, w, "Read depth", sstr("%lld reads (%d reference, %d alternate)", ri(v, 0, V_DP), adr, ada), TH->text);
        if (adr + ada > 0) {
            float ab = (float)ada / (float)(adr + ada);
            txt(F_REG, TH->fsSmall, x, y, TH->dim, "Allele balance");
            Rectangle b = { x + 128, y + 3, w - 128 - 60, 10 };
            float fr[2] = { ab, 1 - ab };
            Color cc[2] = { TH->accent, TH->lineHi };
            stackBar(b, fr, cc, 2);
            txtIn(F_SEMI, TH->fsSmall, (Rectangle){ b.x + b.width + 6, y - 2, 54, 20 }, TH->text, 2, sstr("%.0f%%", ab * 100));
            y += 22;
            const char *note = zyg == 1 && (ab < 0.25f || ab > 0.75f) ? "Unusual for one copy (expected near 50%): check the reads." :
                               zyg == 2 && ab < 0.85f ? "Lower than expected for two copies." : "";
            if (*note) y = vpRow(x, y, w, "", note, TH->warn);
        }
        y = vpRow(x, y, w, "Call quality", sstr("QUAL %.0f, genotype quality %lld", rf(v, 0, V_QUAL), ri(v, 0, V_GQ)), TH->text);
        const char *flt = rc(v, 0, V_FILTER);
        y = vpRow(x, y, w, "Filter", !strcmp(flt, "PASS") ? "Passed all filters" :
                  sstr("%s: %s", flt, strstr(flt, "LowDP") ? "low read depth" : strstr(flt, "LowQUAL") ? "low call quality" : strstr(flt, "LowGQ") ? "low genotype quality" : "flagged"),
                  !strcmp(flt, "PASS") ? TH->ok : TH->warn);
    }
    y = vpSection(x, y + 4, w, "HOW COMMON (1000 GENOMES)");
    if (rnull(v, 0, V_AF)) {
        y = vpRow(x, y, w, "", "Not seen in 1000 Genomes: rare, private to you or your family, or a calling artefact. Check the reads.", TH->warn);
    } else {
        afBar(x, y, w, "Everyone", rf(v, 0, V_AF), TH->accent); y += 22;
        afBar(x, y, w, "European", rf(v, 0, V_AFEUR), (Color){ 110, 178, 206, 255 }); y += 22;
        afBar(x, y, w, "African", rf(v, 0, V_AFAFR), (Color){ 232, 196, 110, 255 }); y += 22;
        afBar(x, y, w, "East Asian", rf(v, 0, V_AFEAS), (Color){ 143, 192, 120, 255 }); y += 22;
        afBar(x, y, w, "South Asian", rf(v, 0, V_AFSAS), (Color){ 180, 142, 200, 255 }); y += 22;
        afBar(x, y, w, "American", rf(v, 0, V_AFAMR), (Color){ 220, 104, 112, 255 }); y += 22;
        txt(F_REG, TH->fsTiny, x + 128, y, TH->faint, "share of chromosomes carrying it (log scale)");
        y += 18;
    }
    if (clin) {
        y = vpSection(x, y + 4, w, "CLINVAR");
        int st = (int)ri(v, 0, V_CLINSTARS);
        y = vpRow(x, y, w, "Classification", sstr("%s", rc(v, 0, V_CLINSIG)), clinColor(clin));
        char stars[64];
        snprintf(stars, sizeof stars, "%s (%d of 4)", st >= 3 ? "expert reviewed" : st == 2 ? "several submitters agree" : st == 1 ? "one submitter, criteria given" : "no review criteria", st);
        y = vpRow(x, y, w, "Review", stars, TH->text);
        char cond[700];
        snprintf(cond, sizeof cond, "%s", rc(v, 0, V_CLINDN));
        for (char *c = cond; *c; c++) if (*c == '|') *c = '\n';
        y = vpRow(x, y, w, "Conditions", cond, TH->text);
        if (gxButtonEx(gxId("vp clinvar"), (Rectangle){ x + 128, y, 140, 28 }, "Open in ClinVar", BTN_NORMAL, true, NULL))
            openUrl(sstr("https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/", rc(v, 0, V_CLINID)));
        y += 38;
    }
    if (*rc(v, 0, V_CSQ) && strcmp(rc(v, 0, V_CSQ), "intron")) {
        y = vpSection(x, y + 4, w, "EFFECT ON THE GENE");
        y = vpRow(x, y, w, "Consequence", csqLabel(rc(v, 0, V_CSQ)), impactColor(imp));
        if (*rc(v, 0, V_HGVSP)) y = vpRow(x, y, w, "Protein", rc(v, 0, V_HGVSP), TH->text);
        if (*rc(v, 0, V_DNA)) y = vpRow(x, y, w, "DNA", rc(v, 0, V_DNA), TH->text);
        y = vpRow(x, y, w, "Transcript", *rc(v, 0, V_TX) ? rc(v, 0, V_TX) : "-", TH->dim);
        if (!rnull(v, 0, V_AM)) {
            double am = rf(v, 0, V_AM);
            txt(F_REG, TH->fsSmall, x, y, TH->dim, "AlphaMissense");
            Rectangle b = { x + 128, y + 4, w - 128 - 60, 10 };
            drRound(b, 5, TH->field);
            drRect((Rectangle){ b.x, b.y, b.width * 0.34f, b.height }, drFade(TH->ok, 0.35f));
            drRect((Rectangle){ b.x + b.width * 0.34f, b.y, b.width * 0.224f, b.height }, drFade(TH->warn, 0.30f));
            drRect((Rectangle){ b.x + b.width * 0.564f, b.y, b.width * 0.436f, b.height }, drFade(TH->bad, 0.35f));
            float mx = b.x + b.width * (float)am;
            drRect((Rectangle){ mx - 1.5f, b.y - 4, 3, b.height + 8 }, TH->text);
            txtIn(F_SEMI, TH->fsSmall, (Rectangle){ b.x + b.width + 6, y - 2, 54, 20 }, TH->text, 2, sstr("%.2f", am));
            y += 22;
            const char *cls = rc(v, 0, V_AMCLASS);
            y = vpRow(x, y, w, "", sstr("Predicted %s. A computer prediction of how much the amino acid change disrupts the protein, not a diagnosis.",
                      !strcmp(cls, "likely_pathogenic") ? "likely damaging" : !strcmp(cls, "likely_benign") ? "likely tolerated" : "ambiguous"), TH->dim);
        }
    }
    if (g_vp.g && g_vp.g->nrow) {
        Rows *g = g_vp.g;
        y = vpSection(x, y + 4, w, "THE GENE");
        y = vpRow(x, y, w, "Name", rc(g, 0, 2), TH->text);
        if (!rnull(g, 0, 4)) {
            double lo = rf(g, 0, 4);
            y = vpRow(x, y, w, "Constraint", sstr("LOEUF %.2f: %s", lo, lo < 0.35 ? "people rarely carry broken copies (intolerant)" :
                      lo < 0.7 ? "somewhat intolerant of loss" : "tolerates loss of one copy"), lo < 0.35 ? TH->warn : TH->text);
        }
        if (*rc(g, 0, 6)) {
            y += 2;
            y += txtPara(F_REG, TH->fsSmall, x, y, w, 18, TH->dim, rc(g, 0, 6), 7) + 8;
        }
    }
    if (g_vp.pk && g_vp.pk->nrow) {
        y = vpSection(x, y + 4, w, "MEDICINES (PHARMGKB)");
        for (int i = 0; i < g_vp.pk->nrow; i++) {
            txt(F_SEMI, TH->fsSmall, x, y, TH->text, sstr("%s  (level %s, you: %s)", rc(g_vp.pk, i, 0), rc(g_vp.pk, i, 1), rc(g_vp.pk, i, 2)));
            y += 20;
            y += txtPara(F_REG, TH->fsSmall, x, y, w, 18, TH->dim, rc(g_vp.pk, i, 3), 6) + 8;
        }
    }
    if (g_vp.gw && g_vp.gw->nrow) {
        y = vpSection(x, y + 4, w, "RESEARCH ASSOCIATIONS (GWAS CATALOG)");
        for (int i = 0; i < g_vp.gw->nrow; i++) {
            int cp = (int)ri(g_vp.gw, i, 2);
            float h = txtParaH(F_REG, TH->fsSmall, w - 110, 18, rc(g_vp.gw, i, 0), 2);
            txtPara(F_REG, TH->fsSmall, x, y, w - 110, 18, TH->text, rc(g_vp.gw, i, 0), 2);
            txtIn(F_SEMI, TH->fsTiny + 1, (Rectangle){ x + w - 106, y, 106, 18 }, cp ? TH->warn : TH->faint, 2,
                  sstr("%s: %d cop%s", rc(g_vp.gw, i, 1), cp, cp == 1 ? "y" : "ies"));
            y += h + 2;
            txt(F_REG, TH->fsTiny, x, y, TH->faint, sstr("%s %s, p 1e-%.0f%s%s", rc(g_vp.gw, i, 5), rc(g_vp.gw, i, 6), rf(g_vp.gw, i, 4),
                                                       *rc(g_vp.gw, i, 3) ? ", effect " : "", rc(g_vp.gw, i, 3)));
            y += 22;
        }
    }
    /* actions */
    y += 8;
    float bx = x;
    if (showBrowse) { if (gxButtonEx(gxId("vp browse"), (Rectangle){ bx, y, 150, 30 }, "Show in browser", BTN_PRIMARY, true, NULL)) navVariant(varId); bx += 158; }
    if (gxButtonEx(gxId("vp copy"), (Rectangle){ bx, y, 70, 30 }, "Copy", BTN_NORMAL, true, "Copy a one-line summary")) {
        copyText(sstr("%s %s>%s %s %s %s %s", fmtPos(chrom, pos), rc(v, 0, V_REF), rc(v, 0, V_ALT), gene, rc(v, 0, V_HGVSP),
                      zygName(zyg), ri(v, 0, V_RSID) ? sstr("rs%lld", ri(v, 0, V_RSID)) : ""));
    }
    bx += 78;
    if (gxButtonEx(gxId("vp gnomad"), (Rectangle){ bx, y, 86, 30 }, "gnomAD", BTN_GHOST, true, "Look it up in gnomAD (your web browser)"))
        openUrl(sstr("https://gnomad.broadinstitute.org/variant/%s-%lld-%s-%s?dataset=gnomad_r4", chromName(chrom), pos, rc(v, 0, V_REF), rc(v, 0, V_ALT)));
    bx += 92;
    if (bx + 120 > x + w) { bx = x; y += 38; }   /* a narrow panel: the last button on a second row */
    if (gxButtonEx(gxId("vp more"), (Rectangle){ bx, y, 120, 30 }, "More links...", BTN_GHOST, true, "SNPedia, ClinVar, PubMed, Ensembl, UCSC and more"))
        lkVariantMenu(varId);
    y += 46;
    g_vpH = y - top;
    gxScrollEnd(&g_vp.scroll, r, y - top);
}

/* ---- the explorer ---- */
typedef struct {
    char text[120];
    int scope, zyg, chrom, freq;
    bool passOnly;
    GxTable t;
    Rows *ids;              /* matching ids in order */
    long long total;        /* all matches (may exceed the ids kept) */
    bool dirty;
    double changedAt;
    /* page cache: 100 rows per page */
    struct { int page; Rows *rows; unsigned age; } cache[24];
    unsigned age;
    long long sel;
} VarExp;
static VarExp g_vx = { .dirty = true, .sel = -1, .t = { .sel = -1, .sortCol = -1 } };
static const char *const g_vxScopes[] = { "All", "Coding", "Protein-altering", "High impact", "In ClinVar", "Pathogenic" };
static const char *const g_vxZyg[] = { "Any zygosity", "Heterozygous", "Homozygous", "Hemizygous" };
static const char *const g_vxFreq[] = { "Any frequency", "Rare (< 1%)", "Very rare (< 0.1%)", "Not in 1000G", "Common (> 5%)" };
static const char *g_vxChroms[26];

static const GxCol g_vxCols[] = {
    { "POSITION", 150, 0, true, NULL }, { "GENE", 96, 0, true, NULL }, { "CHANGE", -2.2f, 0, false, NULL },
    { "EFFECT", 130, 0, true, NULL }, { "COPIES", 66, 0, true, "het: one copy, hom: both, hemi: single-copy chromosome" },
    { "FREQUENCY", 92, 2, true, "How common in 1000 Genomes (all populations)" }, { "CLINVAR", 120, 0, true, NULL },
    { "ALPHAMISS.", 86, 2, true, "AlphaMissense pathogenicity (missense only)" }, { "DEPTH", 60, 2, true, NULL },
};
static const char *const g_vxSortSql[] = { "id", "gene", "", "impact", "zyg", "af", "clin", "am", "dp" };

static char g_vxLastWhere[1200] = "1", g_vxLastOrder[80] = "id", g_vxLastP[4][140];
static int g_vxLastN;
static void vxExport(void) {
    char dir[1024];
    if (!reportDocsDir(dir, sizeof dir)) { toast("No Documents folder to write into"); return; }
    time_t now = time(NULL);
    struct tm tmv;
    localtime_r(&now, &tmv);
    char stamp[32], path[1200], tmp[1300];
    strftime(stamp, sizeof stamp, "%Y-%m-%d %H%M", &tmv);
    snprintf(path, sizeof path, "%s/Strand Open variants %s.csv", dir, stamp);
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    char sql[1700];
    snprintf(sql, sizeof sql, "SELECT CASE chrom WHEN 23 THEN 'X' WHEN 24 THEN 'Y' WHEN 25 THEN 'MT' ELSE chrom END, pos, ref, alt, "
             "CASE WHEN rsid>0 THEN 'rs'||rsid ELSE '' END, gene, csq, hgvs_p, CASE zyg WHEN 1 THEN 'het' WHEN 2 THEN 'hom' WHEN 3 THEN 'hemi' ELSE 'het2' END, "
             "filter, dp, ad_ref, ad_alt, af, af_eur, clin_sig, clin_dn, clin_stars, am FROM variant WHERE %s ORDER BY %s LIMIT 250000", g_vxLastWhere, g_vxLastOrder);
    const char *ps[4] = { g_vxLastP[0], g_vxLastP[1], g_vxLastP[2], g_vxLastP[3] };
    Rows *r = dbqv(sql, g_vxLastN, ps);
    FILE *f = fopen(tmp, "w");
    if (!f) { rowsFree(&r); toast("Could not write the file"); return; }
    fputs("chrom,pos,ref,alt,rsid,gene,consequence,protein,zygosity,filter,depth,ref_reads,alt_reads,af_1000g,af_eur,clinvar,conditions,clinvar_stars,alphamissense\n", f);
    for (int i = 0; i < r->nrow; i++) {
        for (int c = 0; c < r->ncol; c++) {
            const char *v = rc(r, i, c);
            bool q = strpbrk(v, ",\"\n|") != NULL;
            if (q) fputc('"', f);
            for (const char *k = v; *k; k++) { if (*k == '"') fputc('"', f); fputc(*k == '|' ? ';' : *k, f); }
            if (q) fputc('"', f);
            fputc(c + 1 < r->ncol ? ',' : '\n', f);
        }
    }
    int n = r->nrow;
    rowsFree(&r);
    if (!commitFile(f, tmp, path)) { toast("Could not save the file"); return; }
    toast("Saved %s variants to Documents", fmtInt(n));
}
static void vxQuery(void) {
    char where[1200] = "1";
    const char *p[4];
    int n = 0;
    char like[140], rs[40];
    const char *t = g_vx.text;
    while (*t == ' ') t++;
    if (*t) {
        if ((t[0] == 'r' || t[0] == 'R') && (t[1] == 's' || t[1] == 'S') && isdigit((unsigned char)t[2])) {
            snprintf(rs, sizeof rs, "%s", t + 2);
            strcat(where, " AND rsid=?");
            p[n++] = rs;
        } else {
            snprintf(like, sizeof like, "%s", t);
            strcat(where, " AND gene_id IN (SELECT id FROM gene WHERE symbol=? COLLATE NOCASE)");
            p[n++] = like;
        }
    }
    static const char *const scopeSql[] = { "", " AND impact>=1", " AND impact>=2", " AND impact=3", " AND clin>0", " AND clin IN (5,6)" };
    strcat(where, scopeSql[g_vx.scope]);
    static const char *const zygSql[] = { "", " AND zyg IN (1,4)", " AND zyg=2", " AND zyg=3" };
    strcat(where, zygSql[g_vx.zyg]);
    static const char *const freqSql[] = { "", " AND (af IS NULL OR af<0.01)", " AND (af IS NULL OR af<0.001)", " AND af IS NULL", " AND af>0.05" };
    strcat(where, freqSql[g_vx.freq]);
    if (g_vx.chrom > 0) strcat(where, sstr(" AND chrom=%d", g_vx.chrom));
    if (g_vx.passOnly) strcat(where, " AND pass=1");
    char order[80] = "id";
    if (g_vx.t.sortCol >= 0 && *g_vxSortSql[g_vx.t.sortCol]) {
        const char *c = g_vxSortSql[g_vx.t.sortCol];
        bool desc = g_vx.t.sortDesc;
        if (!strcmp(c, "af")) snprintf(order, sizeof order, "af IS NULL %s, af %s, id", desc ? "ASC" : "DESC", desc ? "DESC" : "ASC");
        else snprintf(order, sizeof order, "%s %s, id", c, desc ? "DESC" : "ASC");
    }
    char sql[1600];
    /* kept for the CSV export */
    snprintf(g_vxLastWhere, sizeof g_vxLastWhere, "%s", where);
    snprintf(g_vxLastOrder, sizeof g_vxLastOrder, "%s", order);
    g_vxLastN = n;
    for (int i = 0; i < n; i++) snprintf(g_vxLastP[i], sizeof g_vxLastP[i], "%s", p[i]);
    snprintf(sql, sizeof sql, "SELECT id FROM variant WHERE %s ORDER BY %s LIMIT 250000", where, order);
    dbAsync(SLOT_VARTABLE, sql, n, p);
    snprintf(sql, sizeof sql, "SELECT count(*) FROM variant WHERE %s", where);
    if (!strcmp(where, "1")) dbAsync(SLOT_VARCOUNT, "SELECT value FROM meta WHERE key='n_variants'", 0, NULL);
    else dbAsync(SLOT_VARCOUNT, sql, n, p);
    for (int i = 0; i < 24; i++) rowsFree(&g_vx.cache[i].rows), g_vx.cache[i].page = -1;
}
static Rows *vxPage(int page) {
    for (int i = 0; i < 24; i++) if (g_vx.cache[i].rows && g_vx.cache[i].page == page) { g_vx.cache[i].age = ++g_vx.age; return g_vx.cache[i].rows; }
    if (!g_vx.ids) return NULL;
    int a = page * 100, b = a + 100 > g_vx.ids->nrow ? g_vx.ids->nrow : a + 100;
    if (a >= b) return NULL;
    char *sql = xalloc(2400 + 24 * 100);
    int w = snprintf(sql, 2400, "SELECT id, chrom, pos, ref, alt, gene, hgvs_p, csq, impact, zyg, af, clin, am, dp, rsid, pass FROM variant WHERE id IN (");
    for (int i = a; i < b; i++) w += snprintf(sql + w, 24, "%s%s", i > a ? "," : "", rc(g_vx.ids, i, 0));
    snprintf(sql + w, 8, ")");
    Rows *r = dbq(sql, NULL);
    free(sql);
    /* back into the ids' order */
    Rows *o = r;
    if (r->nrow > 1) {
        int *ord = xalloc(sizeof(int) * (size_t)(b - a));
        for (int i = a; i < b; i++) {
            ord[i - a] = -1;
            for (int j = 0; j < r->nrow; j++) if (!strcmp(rc(r, j, 0), rc(g_vx.ids, i, 0))) { ord[i - a] = j; break; }
        }
        size_t cells = (size_t)r->ncol * (size_t)r->nrow;
        char **tmp = xalloc(sizeof(char *) * cells);
        memcpy(tmp, r->cell, sizeof(char *) * cells);
        for (int i = 0; i < b - a && i < r->nrow; i++)
            if (ord[i] >= 0) memcpy(r->cell + (size_t)i * (size_t)r->ncol, tmp + (size_t)ord[i] * (size_t)r->ncol, sizeof(char *) * (size_t)r->ncol);
        free(tmp);
        free(ord);
    }
    int slot = 0;
    for (int i = 1; i < 24; i++) if (g_vx.cache[i].age < g_vx.cache[slot].age) slot = i;
    rowsFree(&g_vx.cache[slot].rows);
    g_vx.cache[slot].rows = o;
    g_vx.cache[slot].page = page;
    g_vx.cache[slot].age = ++g_vx.age;
    return o;
}
static void vxCell(int row, int col, Rectangle c, bool sel, void *u) {
    (void)u; (void)sel;
    Rows *pg = vxPage(row / 100);
    int i = row % 100;
    if (!pg || i >= pg->nrow) { if (col == 0) txtIn(F_REG, TH->fsSmall, c, TH->faint, 0, "..."); return; }
    int fs = TH->fsSmall;
    switch (col) {
    case 0: txtIn(F_MONO, 12, c, TH->dim, 0, sstr("%s:%s", chromName((int)ri(pg, i, 1)), fmtInt(ri(pg, i, 2)))); break;
    case 1: txtIn(F_SEMI, fs, c, *rc(pg, i, 5) ? TH->text : TH->faint, 0, *rc(pg, i, 5) ? rc(pg, i, 5) : "-"); break;
    case 2: {
        const char *hp = rc(pg, i, 6);
        if (*hp) txtIn(F_REG, fs, c, TH->text, 0, hp);
        else txtIn(F_MONO, 12, c, TH->dim, 0, sstr("%.24s>%.24s", rc(pg, i, 3), rc(pg, i, 4)));
        break;
    }
    case 3: {
        int imp = (int)ri(pg, i, 8);
        DrawCircleV((Vector2){ c.x + 4, c.y + c.height / 2 }, 3.5f, impactColor(imp));
        txtIn(F_REG, fs, (Rectangle){ c.x + 14, c.y, c.width - 14, c.height }, TH->text, 0, csqLabel(rc(pg, i, 7)));
        break;
    }
    case 4: txtIn(F_REG, fs, c, TH->dim, 0, zygShort((int)ri(pg, i, 9))); break;
    case 5: txtIn(F_REG, fs, c, rnull(pg, i, 10) ? TH->warn : TH->text, 2, fmtAf(rf(pg, i, 10))); break;
    case 6: { int cl = (int)ri(pg, i, 11); if (cl) txtIn(F_REG, fs, c, clinColor(cl), 0, clinName(cl)); break; }
    case 7: if (!rnull(pg, i, 12)) { double am = rf(pg, i, 12); txtIn(F_REG, fs, c, am > 0.564 ? TH->bad : am > 0.34 ? TH->warn : TH->dim, 2, sstr("%.2f", am)); } break;
    case 8: txtIn(F_REG, fs, c, ri(pg, i, 15) ? TH->dim : TH->warn, 2, sstr("%lld", ri(pg, i, 13))); break;
    }
}
static void variantsShow(int scope, const char *text, long long sel) {
    g_vx.scope = scope;
    snprintf(g_vx.text, sizeof g_vx.text, "%s", text ? text : "");
    g_vx.dirty = true;
    g_vx.changedAt = 0;
    g_vx.sel = sel;
    navPage(PG_VARIANTS);
}
static bool n_matches_ok(void) { return g_vx.ids && g_vx.ids->nrow > 0; }
static void pageVariants(Rectangle r) {
    if (!g_vxChroms[0]) { g_vxChroms[0] = "All chromosomes"; for (int i = 1; i <= 25; i++) g_vxChroms[i] = strdup(sstr("chr%s", g_chromNames[i])); }
    Rows *t;
    if ((t = dbTake(SLOT_VARTABLE))) { rowsFree(&g_vx.ids); g_vx.ids = t; g_vx.t.scroll = (GxScroll){ 0 }; g_vx.t.sel = -1; }
    if ((t = dbTake(SLOT_VARCOUNT))) { g_vx.total = t->nrow ? ri(t, 0, 0) : 0; rowsFree(&t); }
    if (g_vx.dirty && (GetTime)() - g_vx.changedAt > 0.25) { g_vx.dirty = false; vxQuery(); }

    float sideW = g_vx.sel > 0 ? (r.width > 1400 ? 460 : 400) : 0;
    Rectangle main = { r.x, r.y, r.width - sideW, r.height };
    float x = main.x + 24, y = main.y + 18;
    txt(F_SEMI, TH->fsH2, x, y, TH->text, "Variants");
    const char *cnt = dbBusy(SLOT_VARTABLE) ? "searching..." :
                      g_vx.total > (g_vx.ids ? g_vx.ids->nrow : 0) ? sstr("%s matches (first %s shown in this order)", fmtInt(g_vx.total), fmtInt(g_vx.ids ? g_vx.ids->nrow : 0))
                                                                     : sstr("%s matches", fmtInt(g_vx.total));
    txt(F_REG, TH->fsSmall, x + 110, y + 6, TH->dim, cnt);
    y += 40;
    bool ch = false;
    if (gxTextField(gxId("vx text"), (Rectangle){ x, y, 220, 32 }, g_vx.text, sizeof g_vx.text, "Gene or rsID", true)) ch = true;
    float sx = x + 232;
    float segW = fminf(560, main.width - 640);
    if (gxSegmented(gxId("vx scope"), (Rectangle){ sx, y, segW, 32 }, g_vxScopes, 6, &g_vx.scope)) ch = true;
    sx += segW + 10;
    if (gxCombo(gxId("vx zyg"), (Rectangle){ sx, y, 136, 32 }, g_vxZyg, 4, &g_vx.zyg)) ch = true;
    sx += 144;
    if (gxCombo(gxId("vx freq"), (Rectangle){ sx, y, 150, 32 }, g_vxFreq, 5, &g_vx.freq)) ch = true;
    y += 40;
    if (gxCombo(gxId("vx chrom"), (Rectangle){ x, y, 150, 30 }, g_vxChroms, 26, &g_vx.chrom)) ch = true;
    if (gxCheckbox(gxId("vx pass"), (Rectangle){ x + 162, y + 3, 200, 24 }, &g_vx.passOnly, "Only calls that passed filters")) ch = true;
    if (gxButtonEx(gxId("vx export"), (Rectangle){ x + 380, y, 150, 30 }, "Export CSV", BTN_NORMAL, n_matches_ok(), "Save the matching variants (up to 250,000, in this order) as a CSV file in your Documents folder"))
        vxExport();
    if (ch) { g_vx.dirty = true; g_vx.changedAt = (GetTime)(); }
    y += 40;
    Rectangle tr = { main.x + 16, y, main.width - 32, main.y + main.height - y - 12 };
    int n = g_vx.ids ? g_vx.ids->nrow : 0;
    bool dbl;
    int clicked = gxTable(gxId("vx table"), &g_vx.t, tr, g_vxCols, (int)NELEMS(g_vxCols), n, 30, vxCell, NULL, &dbl);
    if (g_vx.t.sortChanged) { g_vx.dirty = true; g_vx.changedAt = 0; }
    if (clicked >= 0) {
        Rows *pg = vxPage(clicked / 100);
        if (pg && clicked % 100 < pg->nrow) {
            g_vx.sel = ri(pg, clicked % 100, 0);
            if (dbl) navVariant(g_vx.sel);
        }
    }
    /* arrow keys move the selection */
    if (!g_gx.keysTaken && n && (IsKeyPressed(KEY_DOWN) || IsKeyPressedRepeat(KEY_DOWN) || IsKeyPressed(KEY_UP) || IsKeyPressedRepeat(KEY_UP))) {
        int cur = g_vx.t.sel < 0 ? -1 : g_vx.t.sel;
        cur += (IsKeyPressed(KEY_DOWN) || IsKeyPressedRepeat(KEY_DOWN)) ? 1 : -1;
        cur = clampi(cur, 0, n - 1);
        g_vx.t.sel = cur;
        Rows *pg = vxPage(cur / 100);
        if (pg && cur % 100 < pg->nrow) g_vx.sel = ri(pg, cur % 100, 0);
        float rowTop = (float)cur * 30, view = g_vx.t.scroll.viewH;
        if (rowTop < g_vx.t.scroll.target) g_vx.t.scroll.target = rowTop;
        if (rowTop + 30 > g_vx.t.scroll.target + view) g_vx.t.scroll.target = rowTop + 30 - view;
    }
    if (!n && !dbBusy(SLOT_VARTABLE)) txtIn(F_REG, TH->fs, (Rectangle){ tr.x, tr.y + 60, tr.width, 30 }, TH->faint, 1, "No variants match these filters.");
    if (sideW > 0) {
        Rectangle sp = { r.x + r.width - sideW, r.y, sideW, r.height };
        drRect(sp, TH->panel);
        drVLine(sp.x, sp.y, sp.height, TH->line);
        if (gxIconButton(gxId("vx close"), (Rectangle){ sp.x + sp.width - 38, sp.y + 10, 28, 28 }, icClose, false, "Close")) g_vx.sel = -1;
        else variantPanel(sp, g_vx.sel, true);
    }
}
