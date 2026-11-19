/* ========================================================================
   src/page_report.c - the report generator: tick what to include, get a
   PDF in the Documents folder. Each section reads the tables the pages
   use; a section whose data is missing is greyed out.
   ======================================================================== */

typedef struct {
    const char *label, *about, *table;   /* table: the data the section needs */
    bool doctor;                          /* in the "for a doctor" preset */
} RpSection;

enum { RP_BRIEF, RP_MEDS, RP_CLIN, RP_WELLKNOWN, RP_DIET, RP_PRS, RP_MIND, RP_TRAITS, RP_ANCESTRY, RP_LINEAGE, RP_REMARK, RP_RAREST, RP_DATA, RP__N };
static const RpSection g_rpSec[RP__N] = {
    { "Summary (In brief)", "the headline findings, each with what it means", "brief", true },
    { "Medicines", "pharmacogene results and the guidance that differs from standard use", "pgx_gene", true },
    { "Clinical findings", "variants ClinVar classes as pathogenic, with what one or two copies mean", "variant", true },
    { "Well-known variants", "haemochromatosis, clotting factors, APOE, coeliac HLA types and more", "trait", true },
    { "Diet and supplements", "what genes can and cannot say, graded by evidence", "diet", false },
    { "Polygenic scores", "health scores and where you sit among people of your ancestry", "prs", false },
    { "Mind and behaviour scores", "the scores shown last in the app, each with what it can and cannot say", "prs", false },
    { "Traits", "looks, taste and body", "trait", false },
    { "Ancestry", "composition, parents' relatedness and Neanderthal DNA", "anc_admix", false },
    { "Lineage", "maternal and paternal lines and where they come from", "lin_story", false },
    { "Remarkable variants", "notable variants you carry, in plain words", "remarkable", false },
    { "Your rarest variants", "the rarest protein changes and two-copy genotypes", "rarity", false },
    { "About the data", "what the file was, how it was read, and the sources", "meta", true },
};
static bool g_rpOn[RP__N] = { true, true, true, true, true, true, false, true, true, true, true, false, true };
static bool g_rpGenotypes = true;
static bool g_rpDark;                 /* dark pages: light text on near-black */
static GxScroll g_rpScroll;
static bool g_rpLoaded;

/* the choices are remembered: a string of 0s and 1s, one per section, then the two options */
static void rpLoadChoices(void) {
    if (g_rpLoaded) return;
    g_rpLoaded = true;
    const char *v = cfgGet("report_sections", "");
    if (strlen(v) == RP__N) for (int i = 0; i < RP__N; i++) g_rpOn[i] = v[i] == '1';
    g_rpGenotypes = strcmp(cfgGet("report_genotypes", "1"), "0") != 0;
    g_rpDark = !strcmp(cfgGet("report_dark", "0"), "1");
}
static void rpSaveChoices(void) {
    char v[RP__N + 1];
    for (int i = 0; i < RP__N; i++) v[i] = g_rpOn[i] ? '1' : '0';
    v[RP__N] = 0;
    cfgSet("report_sections", v);
    cfgSet("report_genotypes", g_rpGenotypes ? "1" : "0");
    cfgSet("report_dark", g_rpDark ? "1" : "0");
    cfgSave();
}
static char g_rpLast[1200];

/* ---- building the PDF ---- */
/* a cell's text, tidied: list separators ('|' or ';') become "; " */
static const char *rpCell(char *out, size_t cap, const char *in) {
    size_t n = 0;
    for (const char *q = in; *q && n + 3 < cap; q++) {
        if (*q == '|' || *q == ';') {
            while (n && out[n - 1] == ' ') n--;
            if (n) { out[n++] = ';'; out[n++] = ' '; }
            while (q[1] == ' ') q++;
            continue;
        }
        out[n++] = *q;
    }
    while (n && (out[n - 1] == ' ' || out[n - 1] == ';')) n--;
    out[n] = 0;
    return out;
}
/* a chip that misses positions can leave dozens of possible star-allele results: the first three, then a count */
static void rpShortenOr(char *s, size_t cap) {
    int n = 0;
    char *cut = NULL;
    for (char *q = strstr(s, " or "); q; q = strstr(q + 4, " or ")) if (++n == 3) cut = q;
    if (n < 5 || !cut) return;
    snprintf(cut, cap - (size_t)(cut - s), " or %d others (not enough positions tested to tell them apart)", n - 2);
}
static void rpTable(Pdf *p, const char *sql, const char *const *head, const float *w, int ncol, int monoCol) {
    Rows *r = dbq(sql, NULL);
    if (!r || !r->nrow) { pdfPara(p, "None.", 9, PDF_DIM, PDF_REG); rowsFree(&r); return; }
    static char clean[8][2400];   /* long guidance texts must not be cut */
    const char *cells[8];
    for (int i = 0; i < r->nrow; i++) {
        for (int c = 0; c < ncol && c < 8; c++) { cells[c] = rpCell(clean[c], sizeof clean[c], rc(r, i, c)); rpShortenOr(clean[c], sizeof clean[c]); }
        if (i == 0) pdfTableBegin(p, head, w, ncol, cells, monoCol);
        pdfRow(p, cells, w, ncol, monoCol);
    }
    pdfTableEnd(p);
    rowsFree(&r);
}

static void rpBuild(Pdf *p) {
    char day[32];
    time_t now = time(NULL);
    struct tm tmv;
    localtime_r(&now, &tmv);
    strftime(day, sizeof day, "%d %B %Y", &tmv);
    const char *who = *dbMeta("sample_id") ? dbMeta("sample_id") : "This person";
    snprintf(p->footer, sizeof p->footer, "Genome report: %s, %s. Not a diagnosis.", who, day);
    pdfNewPage(p);
    pdfPara(p, "Genome report", 22, PDF_HEAD, PDF_BOLD);
    pdfPara(p, sstr("%s  |  %s  |  %s", who, *dbMeta("seq_kind") ? dbMeta("seq_kind") : "DNA data", day), 10, PDF_DIM, PDF_REG);
    pdfPara(p, "Made on this computer from the person's own DNA file. This report is for curiosity and for preparing questions, not a "
               "medical test: anything that might matter for health deserves a clinical-grade confirmation and a conversation with a "
               "doctor, pharmacist or genetic counsellor.", 9, PDF_DIM, PDF_REG);
    bool chip = !strcmp(dbMeta("call_model"), "chip");

    if (g_rpOn[RP_BRIEF] && dbHasTable("brief")) {
        pdfHeading(p, "Summary", 1);
        Rows *r = dbq("SELECT section, title, text, meaning FROM brief ORDER BY ord", NULL);
        const char *sec = "";
        for (int i = 0; r && i < r->nrow; i++) {
            if (strcmp(sec, rc(r, i, 0))) { sec = rc(r, i, 0); pdfHeading(p, sec, 2); }
            pdfPara(p, rc(r, i, 1), 9.5f, PDF_INK, PDF_BOLD);
            if (*rc(r, i, 2)) pdfPara(p, rc(r, i, 2), 9, PDF_INK, PDF_REG);
            if (*rc(r, i, 3)) pdfPara(p, sstr("What it means: %s", rc(r, i, 3)), 8.5f, PDF_DIM, PDF_REG);
        }
        rowsFree(&r);
    }
    if (g_rpOn[RP_MEDS] && dbHasTable("pgx_gene")) {
        pdfHeading(p, "Medicines", 1);
        pdfHeading(p, "Pharmacogene results", 2);
        static const char *const h1[] = { "Gene", "Result", "Phenotype" };
        static const float w1[] = { 90, 200, 209 };
        rpTable(p, "SELECT gene, diplotype, phenotype FROM pgx_gene WHERE called=1 ORDER BY gene", h1, w1, 3, -1);
        if (dbHasTable("pgx_drug")) {
            pdfHeading(p, "Guidance that differs from standard use", 2);
            static const char *const h2[] = { "Medicine", "Source", "Your result", "Guidance" };
            static const float w2[] = { 95, 45, 110, 249 };
            rpTable(p, "SELECT drug, CASE WHEN source LIKE 'CPIC%' THEN 'CPIC' WHEN source LIKE 'DPWG%' THEN 'DPWG' ELSE source END, phenotype, recommendation "
                       "FROM pgx_drug WHERE recommendation!='' AND source NOT LIKE 'FDA%' "
                       "AND (phenotype LIKE '%Poor%' OR phenotype LIKE '%Intermediate%' OR phenotype LIKE '%Rapid%' OR phenotype LIKE '%Decreased%' "
                       "OR phenotype LIKE '%Increased%' OR phenotype LIKE '%positive%' OR phenotype LIKE '%Deficien%') "
                       "AND recommendation NOT LIKE '%does not provide a recommendation%' AND recommendation NOT LIKE 'Initiate therapy with recommended starting dose%' "
                       "AND recommendation NOT LIKE 'Prescribe desired starting dose%' AND recommendation NOT LIKE '%per standard dosing%' "
                       "AND recommendation NOT LIKE 'No recommendation%' AND recommendation NOT LIKE '%no reason to avoid%' "
                       "AND recommendation NOT LIKE 'Initiate standard%' AND recommendation NOT LIKE '%no action is needed%' AND recommendation NOT LIKE '%no action is required%' "
                       "AND recommendation NOT LIKE '%use at standard dose%' AND recommendation NOT LIKE 'No action recommended%' AND recommendation NOT LIKE 'Initiate therapy with standard%' "
                       "GROUP BY drug, source ORDER BY drug", h2, w2, 4, -1);
        }
        Rows *hn = dbHasTable("hla") ? dbq("SELECT count(*) FROM hla", NULL) : NULL;
        bool typed = hn && hn->nrow && ri(hn, 0, 0) > 0;
        rowsFree(&hn);
        if (typed) {
            pdfHeading(p, "HLA type", 2);
            static const char *const h3[] = { "Gene", "Allele 1", "Allele 2" };
            static const float w3[] = { 120, 190, 189 };
            rpTable(p, "SELECT gene, a1, a2 FROM hla WHERE gene IN ('HLA-A','HLA-B','HLA-C','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DPB1')", h3, w3, 3, -1);
        }
    }
    if (g_rpOn[RP_CLIN] && dbHasTable("variant")) {
        pdfHeading(p, "Clinical findings", 1);
        pdfPara(p, "Variants ClinVar classes as pathogenic or likely pathogenic. AR: recessive (one copy means carrier status); AD: dominant "
                   "(one copy can be enough).", 8.5f, PDF_DIM, PDF_REG);
        static const char *const h[] = { "Gene", "Variant", "Copies", "Conditions", "Inheritance" };
        static const float w[] = { 60, 95, 65, 219, 60 };
        rpTable(p, g_rpGenotypes
                   ? "SELECT v.gene, coalesce(v.hgvs_p, 'chr'||(CASE v.chrom WHEN 23 THEN 'X' WHEN 24 THEN 'Y' WHEN 25 THEN 'MT' ELSE v.chrom END)||':'||v.pos||' '||v.ref||'>'||v.alt), "
                     "CASE v.zyg WHEN 1 THEN 'one copy' WHEN 2 THEN 'two copies' WHEN 3 THEN 'single copy' ELSE 'one of each' END, v.clin_dn, "
                     "replace((SELECT group_concat(DISTINCT nullif(d.inheritance, '')) FROM gene_disease d WHERE d.gene=v.gene_id), ',', ', ') FROM variant v WHERE v.clin IN (5,6) ORDER BY v.gene"
                   : "SELECT v.gene, coalesce(v.hgvs_p, ''), CASE v.zyg WHEN 1 THEN 'one copy' WHEN 2 THEN 'two copies' WHEN 3 THEN 'single copy' ELSE 'one of each' END, v.clin_dn, "
                     "replace((SELECT group_concat(DISTINCT nullif(d.inheritance, '')) FROM gene_disease d WHERE d.gene=v.gene_id), ',', ', ') FROM variant v WHERE v.clin IN (5,6) ORDER BY v.gene",
                h, w, 5, -1);
        if (chip) {
            Rows *n = dbq("SELECT count(*) FROM variant WHERE clin IN (105,106)", NULL);
            long long k = n && n->nrow ? ri(n, 0, 0) : 0;
            rowsFree(&n);
            if (k) pdfPara(p, sstr("%lld rare chip call%s that ClinVar lists as disease-causing %s set aside as unconfirmed: chip calls this rare "
                                   "are usually wrong (Weedon et al. 2021), so a clinical test may make sense before reading anything into %s.",
                                   k, k == 1 ? "" : "s", k == 1 ? "was" : "were", k == 1 ? "it" : "them"), 8.5f, PDF_DIM, PDF_REG);
        }
    }
    if (g_rpOn[RP_WELLKNOWN] && dbHasTable("trait")) {
        pdfHeading(p, "Well-known variants", 1);
        static const char *const h[] = { "Variant", "Result", "Genotype" };
        static const float w[] = { 190, 150, 159 };
        static const char *const h2[] = { "Variant", "Result", "What it means" };
        static const float w2[] = { 170, 120, 209 };
        if (g_rpGenotypes) rpTable(p, "SELECT title, result, genotype FROM trait WHERE cat='Health' ORDER BY ord", h, w, 3, 2);
        else rpTable(p, "SELECT title, result, detail FROM trait WHERE cat='Health' ORDER BY ord", h2, w2, 3, -1);
    }
    if (g_rpOn[RP_DIET] && dbHasTable("diet")) {
        pdfHeading(p, "Diet and supplements", 1);
        Rows *r = dbq("SELECT ord, topic, grade, result, genotype, says, means, source FROM diet ORDER BY ord", NULL);
        for (int i = 0; r && i < r->nrow; i++) {
            bool intro = ri(r, i, 0) < 0;
            pdfPara(p, intro ? rc(r, i, 1) : sstr("%s: %s  (evidence: %s)", rc(r, i, 1), rc(r, i, 3), rc(r, i, 2)), 9.5f, PDF_INK, PDF_BOLD);
            if (g_rpGenotypes && *rc(r, i, 4)) pdfPara(p, rc(r, i, 4), 8, PDF_DIM, PDF_MONO);
            pdfPara(p, rc(r, i, 5), 9, PDF_INK, PDF_REG);
            pdfPara(p, sstr("%s%s", intro ? "" : "What it means: ", rc(r, i, 6)), 8.5f, PDF_DIM, PDF_REG);
            if (*rc(r, i, 7)) pdfPara(p, sstr("Sources: %s", rc(r, i, 7)), 7.5f, PDF_DIM, PDF_REG);
        }
        rowsFree(&r);
    }
    for (int pass = 0; pass < 2; pass++) {
        int sec = pass == 0 ? RP_PRS : RP_MIND;
        if (!g_rpOn[sec] || !dbHasTable("prs")) continue;
        pdfHeading(p, pass == 0 ? "Polygenic scores" : "Mind and behaviour scores", 1);
        if (pass == 1)
            pdfPara(p, sstr("Each score adds up the small effects of many common variants and is placed among %s. These scores explain only a "
                            "few percent of the differences between people; upbringing, circumstances and chance matter far more, so they say "
                            "little about any one person.", *dbMeta("prs_reference") ? dbMeta("prs_reference") : "the Europeans of 1000 Genomes"),
                    8.5f, PDF_DIM, PDF_REG);
        else pdfPara(p, sstr("Each score adds up the small effects of many common variants and is placed among %s. A score shifts the odds; it "
                        "does not decide them, and most people with a high score never develop the condition.%s%s",
                        *dbMeta("prs_reference") ? dbMeta("prs_reference") : "the Europeans of 1000 Genomes",
                        *dbMeta("prs_portability") ? " " : "", dbMeta("prs_portability")), 8.5f, PDF_DIM, PDF_REG);
        static const char *const h[] = { "Score", "Area", "Percentile" };
        static const float w[] = { 260, 159, 80 };
        static const char *const hm[] = { "Score", "Percentile", "What it can and cannot say" };
        static const float wm[] = { 120, 60, 319 };
        if (pass == 0)
            rpTable(p, "SELECT trait, category, printf('%.0f', pct_eur) FROM prs WHERE category != 'Mind and behaviour' ORDER BY category, trait", h, w, 3, -1);
        else
            rpTable(p, "SELECT trait, printf('%.0f', pct_eur), note FROM prs WHERE category = 'Mind and behaviour' ORDER BY trait", hm, wm, 3, -1);
    }
    if (g_rpOn[RP_TRAITS] && dbHasTable("trait")) {
        pdfHeading(p, "Traits", 1);
        static const char *const h[] = { "Trait", "Result", "Genotype" };
        static const float w[] = { 200, 150, 149 };
        static const char *const h2[] = { "Trait", "Result", "What it means" };
        static const float w2[] = { 170, 130, 199 };
        if (g_rpGenotypes) rpTable(p, "SELECT title, result, genotype FROM trait WHERE cat != 'Health' ORDER BY ord", h, w, 3, 2);
        else rpTable(p, "SELECT title, result, detail FROM trait WHERE cat != 'Health' ORDER BY ord", h2, w2, 3, -1);
    }
    if (g_rpOn[RP_ANCESTRY] && dbHasTable("anc_admix")) {
        pdfHeading(p, "Ancestry", 1);
        if (*dbMeta("ancestry_summary")) pdfPara(p, dbMeta("ancestry_summary"), 9.5f, PDF_INK, PDF_BOLD);
        static const char *const h[] = { "Population", "Share", "95% interval" };
        static const float w[] = { 260, 100, 139 };
        rpTable(p, "SELECT label, printf('%.0f%%', frac*100), printf('%.0f to %.0f%%', lo*100, hi*100) FROM anc_admix WHERE level='pop' AND frac >= 0.01 ORDER BY frac DESC", h, w, 3, -1);
        if (*dbMeta("roh_summary")) { pdfHeading(p, "Did the parents share recent ancestors?", 2); pdfPara(p, dbMeta("roh_summary"), 9, PDF_INK, PDF_REG); }
        if (*dbMeta("archaic_frac")) {
            pdfHeading(p, "Neanderthal and other archaic DNA", 2);
            pdfPara(p, sstr("About %.1f%% of the genome sits in %s stretches inherited from archaic humans. %s", atof(dbMeta("archaic_frac")) * 100,
                            dbMeta("archaic_segments"), dbMeta("archaic_ref")), 9, PDF_INK, PDF_REG);
        } else if (*dbMeta("archaic_unmeasurable")) {
            pdfHeading(p, "Neanderthal and other archaic DNA", 2);
            pdfPara(p, dbMeta("archaic_unmeasurable"), 9, PDF_INK, PDF_REG);
        }
    }
    if (g_rpOn[RP_LINEAGE] && dbHasTable("lin_story")) {
        pdfHeading(p, "Lineage", 1);
        for (int k = 0; k < 2; k++) {
            const char *hg = dbMeta(k == 0 ? "mt_haplogroup" : "y_haplogroup");
            pdfHeading(p, sstr("%s: %s", k == 0 ? "Maternal line (mitochondrial DNA)" : "Paternal line (Y chromosome)",
                               *hg ? hg : (k == 0 ? "not in this data" : "none in this data")), 2);
            if (!*hg) continue;
            Rows *r = dbq("SELECT text FROM lin_story WHERE line=? ORDER BY ord", k == 0 ? "MT" : "Y", NULL);
            for (int i = 0; r && i < r->nrow; i++) pdfPara(p, rc(r, i, 0), 9, PDF_INK, PDF_REG);
            rowsFree(&r);
        }
    }
    if (g_rpOn[RP_REMARK] && dbHasTable("remarkable")) {
        pdfHeading(p, "Remarkable variants", 1);
        Rows *r = dbq("SELECT title, found, effect, rarity FROM remarkable WHERE status > 0 ORDER BY status DESC, cat", NULL);
        if (!r || !r->nrow) pdfPara(p, "None of the catalogue's variants were found.", 9, PDF_DIM, PDF_REG);
        for (int i = 0; r && i < r->nrow; i++) {
            pdfPara(p, rc(r, i, 0), 9.5f, PDF_INK, PDF_BOLD);
            if (g_rpGenotypes && *rc(r, i, 1)) pdfPara(p, sstr("Found: %s", rc(r, i, 1)), 8, PDF_DIM, PDF_REG);
            pdfPara(p, rc(r, i, 2), 9, PDF_INK, PDF_REG);
            if (*rc(r, i, 3)) pdfPara(p, sstr("How common: %s", rc(r, i, 3)), 8, PDF_DIM, PDF_REG);
        }
        rowsFree(&r);
    }
    if (g_rpOn[RP_RAREST] && dbHasTable("rarity")) {
        pdfHeading(p, "Your rarest variants", 1);
        pdfHeading(p, "Rarest protein changes", 2);
        static const char *const hp[] = { "Gene", "Change", "Carried by", "What the change does", "The gene" };
        static const float wp[] = { 58, 78, 92, 160, 111 };
        Rows *r = dbq("SELECT gene, change, world, note, plain FROM rarity WHERE kind='protein' ORDER BY rank LIMIT 12", NULL);
        if (!r || !r->nrow) pdfPara(p, "None.", 9, PDF_DIM, PDF_REG);
        for (int i = 0; r && i < r->nrow; i++) {
            char who[64];
            snprintf(who, sizeof who, "%s people worldwide", oneIn(rf(r, i, 2)));
            const char *cells[5] = { rc(r, i, 0), rc(r, i, 1), who, rc(r, i, 3), rc(r, i, 4) };
            if (i == 0) pdfTableBegin(p, hp, wp, 5, cells, -1);
            pdfRow(p, cells, wp, 5, -1);
        }
        pdfTableEnd(p);
        rowsFree(&r);
        pdfHeading(p, "Rarest two-copy genotypes", 2);
        static const char *const hh[] = { "Gene", "Change", "How rare", "The gene" };
        static const float wh[] = { 70, 100, 150, 179 };
        rpTable(p, "SELECT gene, change, note, plain FROM rarity WHERE kind='homo' ORDER BY rank LIMIT 12", hh, wh, 4, -1);
    }
    if (g_rpOn[RP_DATA]) {
        pdfHeading(p, "About the data", 1);
        pdfPara(p, sstr("%s, on GRCh38. %s variant calls.", *dbMeta("seq_kind") ? dbMeta("seq_kind") : "DNA data", fmtInt(atoll(dbMeta("n_variants")))),
                9, PDF_INK, PDF_REG);
        pdfPara(p, "Sources: GRCh38 and gene models (Ensembl 115), ClinVar (NCBI), AlphaMissense, gnomAD constraint, UniProt, the Human Phenotype "
                   "Ontology, PharmCAT with CPIC and DPWG guidelines, PharmGKB, the GWAS Catalog, the PGS Catalog, the 1000 Genomes 30x panel, "
                   "PhyloTree, the YFull tree, and the Allen Ancient DNA Resource.", 8.5f, PDF_DIM, PDF_REG);
        pdfPara(p, chip ? "Chip calls can be wrong, above all for rare variants. ClinVar classifications are revised over time, and risk estimates "
                          "from common variants are averages over populations."
                        : "Calls from short-read sequencing can be wrong, especially in repetitive DNA and for insertions and deletions. ClinVar "
                          "classifications are revised over time, and risk estimates from common variants are averages over populations.",
                8.5f, PDF_DIM, PDF_REG);
    }
}

static void rpCreate(void) {
    char dir[1024], path[1300];
    if (!reportDocsDir(dir, sizeof dir)) { toast("No Documents folder to write into"); return; }
    time_t now = time(NULL);
    struct tm tmv;
    localtime_r(&now, &tmv);
    char day[32];
    strftime(day, sizeof day, "%Y-%m-%d", &tmv);
    char who[120];
    snprintf(who, sizeof who, "%s", *dbMeta("sample_id") ? dbMeta("sample_id") : "genome");
    for (char *c = who; *c; c++) if (*c == '/' || *c == '\\' || *c == ':') *c = '-';
    snprintf(path, sizeof path, "%s/Genome report %s %s.pdf", dir, who, day);
    rpSaveChoices();
    Pdf p = { 0 };
    pdfTheme(g_rpDark);
    rpBuild(&p);
    bool ok = pdfSave(&p, path, sstr("Genome report: %s", who));
    int pages = p.npages;
    pdfFree(&p);
    if (!ok) { toast("Could not save the report"); return; }
    snprintf(g_rpLast, sizeof g_rpLast, "%s", path);
    toast("Saved to Documents: %d pages", pages);
    if (getenv("STRAND_NO_OPEN")) { printf("report: %s\n", path); return; }
    char url[1400];
    snprintf(url, sizeof url, "file://%s", path);
    openUrl(url);
}

/* ---- the page ---- */
static bool rpHave(int i) { return !strcmp(g_rpSec[i].table, "meta") || dbHasTable(g_rpSec[i].table); }

static void pageReport(Rectangle r) {
    rpLoadChoices();
    float top = gxScrollBegin(&g_rpScroll, r);
    float x = r.x + 32, W = fminf(1000, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Report generator",
                         "A PDF of the results to keep, print or take to a doctor. Tick what to include; sections without the data for "
                         "them are greyed out. The file is made on this computer and saved in your Documents folder.");
    if (gxButtonEx(gxId("rp doctor"), (Rectangle){ x, y, 190, 32 }, "For a doctor", BTN_NORMAL, true,
                   "Summary, medicines, clinical findings and well-known variants"))
        for (int i = 0; i < RP__N; i++) g_rpOn[i] = g_rpSec[i].doctor;
    if (gxButtonEx(gxId("rp all"), (Rectangle){ x + 202, y, 150, 32 }, "Everything", BTN_NORMAL, true, NULL))
        for (int i = 0; i < RP__N; i++) g_rpOn[i] = true;
    if (gxButtonEx(gxId("rp none"), (Rectangle){ x + 364, y, 120, 32 }, "Clear", BTN_GHOST, true, NULL))
        for (int i = 0; i < RP__N; i++) g_rpOn[i] = false;
    y += 48;
    for (int i = 0; i < RP__N; i++) {
        bool have = rpHave(i);
        Rectangle row = { x, y, W, 36 };
        if (i & 1) drRound(row, 6, TH->zebra);
        if (have) {
            gxCheckbox(gxIdN(73000 + i), (Rectangle){ x + 12, y + 7, 280, 22 }, &g_rpOn[i], g_rpSec[i].label);
        } else {
            drRoundLine((Rectangle){ x + 12, y + 10, 16, 16 }, 4, TH->line);   /* an empty box: nothing to tick */
            txtIn(F_REG, TH->fsSmall, (Rectangle){ x + 36, y + 7, 256, 22 }, TH->faint, 0, g_rpSec[i].label);
        }
        txtIn(F_REG, TH->fsSmall, (Rectangle){ x + 300, y, W - 312, 36 }, have ? TH->dim : TH->faint, 0,
              have ? g_rpSec[i].about : "not available for this person's data");
        y += 37;
    }
    y += 12;
    gxCheckbox(gxId("rp genotypes"), (Rectangle){ x + 12, y, W - 24, 22 }, &g_rpGenotypes, "Show genotypes and variant names (for a doctor or a lab)");
    y += 32;
    gxCheckbox(gxId("rp dark"), (Rectangle){ x + 12, y, W - 24, 22 }, &g_rpDark, "Dark pages (light text on a dark background; for reading on screen rather than printing)");
    y += 44;
    int n = 0;
    for (int i = 0; i < RP__N; i++) n += g_rpOn[i] && rpHave(i);
    if (gxButtonEx(gxId("rp create"), (Rectangle){ x, y, 220, 38 }, "Create the PDF", BTN_PRIMARY, n > 0, NULL)) rpCreate();
    if (g_rpLast[0]) {
        const char *base = strrchr(g_rpLast, '/');
        txtIn(F_REG, TH->fsSmall, (Rectangle){ x + 236, y, W - 236, 38 }, TH->faint, 0, sstr("Saved in Documents: %s", base ? base + 1 : g_rpLast));
    }
    y += 70;
    gxScrollEnd(&g_rpScroll, r, y - top);
}
