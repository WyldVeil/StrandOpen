/* ========================================================================
   src/page_raw.c - the raw data and its quality (files, sequencing reads,
   alignment, calls, the chip cross-check), and the sources and notes page.
   ======================================================================== */

static Rows *g_rwFiles, *g_rwQc, *g_rwAln;
static GxScroll g_rwScroll;
static float g_rwPosQ[2][300], g_rwGc[2][101], g_rwQm[2][64];
static int g_rwPosN[2];
static bool g_rwParsed;

static int parseFloats(const char *s, float *out, int max) {
    int n = 0;
    while (s && *s && n < max) {
        while (*s && !(isdigit((unsigned char)*s) || *s == '-' || *s == '.')) s++;
        if (!*s) break;
        char *e;
        out[n++] = strtof(s, &e);
        if (e == s) break;
        s = e;
    }
    return n;
}
static void lineChart(Rectangle r, const float *v, int n, float lo, float hi, Color c, const char *yl, const char *xl) {
    drRect(r, (Color){ 22, 25, 31, 255 });
    for (int k = 0; k <= 4; k++) {
        float y = r.y + r.height * (float)k / 4;
        drHLine(r.x, y, r.width, drFade(TH->line, 0.6f));
        txt(F_REG, TH->fsTiny - 1, r.x - 26, y - 6, TH->faint, sstr("%.0f", hi - (hi - lo) * (float)k / 4));
    }
    Vector2 prev = { 0 };
    for (int i = 0; i < n; i++) {
        Vector2 p = { r.x + r.width * (float)i / (float)(n > 1 ? n - 1 : 1), r.y + r.height * (1 - clamp01((v[i] - lo) / (hi - lo))) };
        if (i) drLine(prev.x, prev.y, p.x, p.y, 2, c);
        prev = p;
    }
    if (yl) txt(F_REG, TH->fsTiny, r.x, r.y - 16, TH->faint, yl);
    if (xl) txtIn(F_REG, TH->fsTiny, (Rectangle){ r.x, r.y + r.height + 4, r.width, 14 }, TH->faint, 1, xl);
}
static void qcCard(Rectangle c, const char *label, const char *value, const char *sub, Color col) {
    gxCard(c);
    txt(F_REG, TH->fsSmall, c.x + 16, c.y + 14, TH->dim, label);
    txt(F_SEMI, TH->fsH2 + 2, c.x + 16, c.y + 36, col.a ? col : TH->text, value);
    if (sub) txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 16, c.y + 68, c.width - 30, 16 }, TH->faint, 0, sub);
}
static void pageRaw(Rectangle r) {
    if (!g_rwFiles) {
        g_rwFiles = dbHasTable("files") ? dbq("SELECT name, size, kind, about FROM files ORDER BY ord", NULL) : NULL;
        g_rwQc = dbHasTable("qc_fastq") ? dbq("SELECT mate, reads, bases, q30, mean_q, gc, n_rate, adapter, dup, pos_q, gc_hist, qmean_hist, len_mode FROM qc_fastq ORDER BY mate", NULL) : NULL;
        g_rwAln = dbHasTable("qc_align") ? dbq("SELECT key, value, note FROM qc_align ORDER BY ord", NULL) : NULL;
    }
    if (!g_rwParsed && g_rwQc) {
        g_rwParsed = true;
        for (int m = 0; m < g_rwQc->nrow && m < 2; m++) {
            g_rwPosN[m] = parseFloats(rc(g_rwQc, m, 9), g_rwPosQ[m], 300);
            parseFloats(rc(g_rwQc, m, 10), g_rwGc[m], 101);
            parseFloats(rc(g_rwQc, m, 11), g_rwQm[m], 64);
        }
    }
    float top = gxScrollBegin(&g_rwScroll, r);
    float x = r.x + 32, W = fminf(1300, r.width - 64);
    bool chip = dkChip();
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Raw data and quality",
                         sstr("What this profile was built from (%s), what the import made of it, and how good it is.", dkName()));
    float cw = (W - 48) / 4, ch = 92;
    if (chip) {
        /* a chip: what was read, what was dropped and why */
        gxSectionTitle(x, y, "THE CHIP FILE");
        y += 24;
        qcCard((Rectangle){ x, y, cw, ch }, "Markers in the file", fmtInt(atoll(dbMeta("ingest_markers_in_file"))), dkName(), BLANK);
        qcCard((Rectangle){ x + cw + 16, y, cw, ch }, "Positions read", fmtInt(atoll(dbMeta("callable_bases"))),
               sstr("%s matching the reference, %s not", fmtBig(atof(dbMeta("ingest_tested_ref"))), fmtBig(atof(dbMeta("ingest_tested_variant")))), TH->ok);
        qcCard((Rectangle){ x + 2 * (cw + 16), y, cw, ch }, "No-calls", fmtInt(atoll(dbMeta("ingest_dropped_nocall"))), "the chip could not read these", BLANK);
        qcCard((Rectangle){ x + 3 * (cw + 16), y, cw, ch }, "Insertion/deletion probes", fmtInt(atoll(dbMeta("ingest_dropped_indel"))),
               "reported as I/D without the letters: not usable", BLANK);
        y += ch + 16;
        qcCard((Rectangle){ x, y, cw, ch }, "Genome build", sstr("GRCh%s", dbMeta("source_build")),
               !strcmp(dbMeta("source_build"), "37") ? sstr("lifted to GRCh38; %s could not be", fmtInt(atoll(dbMeta("ingest_dropped_unlifted")))) : "used as it is", BLANK);
        qcCard((Rectangle){ x + cw + 16, y, cw, ch }, "Strand corrections", fmtInt(atoll(dbMeta("ingest_strand_flips"))),
               "calls reported on the other DNA strand, turned round", BLANK);
        qcCard((Rectangle){ x + 2 * (cw + 16), y, cw, ch }, "Sex chromosomes", *dbMeta("sex") ? (!strcmp(dbMeta("sex"), "male") ? "XY" : !strcmp(dbMeta("sex"), "female") ? "XX" : "?") : "-",
               sstr("Y markers read %s; X heterozygosity %s", dbMeta("ingest_y_called"), dbMeta("ingest_x_het")), BLANK);
        double het = atof(dbMeta("n_het")), hom = atof(dbMeta("n_hom"));
        qcCard((Rectangle){ x + 3 * (cw + 16), y, cw, ch }, "Het / hom ratio", sstr("%.2f", het / fmax(1, hom)), "of the positions that differ", BLANK);
        y += ch + 16;
        Rectangle c = { x, y, W, 0 };
        const char *note = "A genotyping chip reads a fixed list of positions, picked because they vary between people. That makes it good at common, "
                           "well-studied variants and poor at rare ones: in a study of nearly 50,000 people, most very rare variants reported by consumer chips "
                           "turned out to be wrong when checked by sequencing (Weedon et al. 2021, BMJ). So the rare-variant pages, read quality and "
                           "coverage need sequencing data, and for anything rare or alarming found here a clinical test may make sense before acting on it.";
        c.height = txtParaH(F_REG, TH->fsSmall, W - 40, 19, note, 8) + 36;
        gxCard(c);
        txtPara(F_REG, TH->fsSmall, c.x + 20, c.y + 18, W - 40, 19, TH->text, note, 8);
        y += c.height + 28;
    } else {
    /* calls */
    gxSectionTitle(x, y, "VARIANT CALLS");
    y += 24;
    qcCard((Rectangle){ x, y, cw, ch }, "Callable autosomes", sstr("%.1f%%", atof(dbMeta("autosome_callable_fraction")) * 100), "bases confidently called", TH->ok);
    qcCard((Rectangle){ x + cw + 16, y, cw, ch }, "Ti/Tv (passing SNVs)", dbMeta("titv"), "2.0 to 2.1 expected for a whole genome", BLANK);
    double het = atof(dbMeta("n_het")), hom = atof(dbMeta("n_hom"));
    qcCard((Rectangle){ x + 2 * (cw + 16), y, cw, ch }, "Het / hom ratio", sstr("%.2f", het / fmax(1, hom)), "1.5 to 2.0 typical for European ancestry", BLANK);
    qcCard((Rectangle){ x + 3 * (cw + 16), y, cw, ch }, "Passed filters", sstr("%.1f%%", 100.0 * atof(dbMeta("n_pass")) / fmax(1, atof(dbMeta("n_variants")))),
           sstr("%s of %s calls", fmtInt(atoll(dbMeta("n_pass"))), fmtInt(atoll(dbMeta("n_variants")))), BLANK);
    y += ch + 16;
    qcCard((Rectangle){ x, y, cw, ch }, "Sex chromosomes", *dbMeta("sex") ? (!strcmp(dbMeta("sex"), "male") ? "XY" : "XX") : "-",
           sstr("X callable %.0f%%, Y callable %.0f%%", atof(dbMeta("callable_fraction_X")) * 100, atof(dbMeta("callable_fraction_Y")) * 100), BLANK);
    qcCard((Rectangle){ x + cw + 16, y, cw, ch }, "With an rsID", sstr("%.1f%%", 100.0 * atof(dbMeta("n_rsid")) / fmax(1, atof(dbMeta("n_variants")))), "known to dbSNP", BLANK);
    qcCard((Rectangle){ x + 2 * (cw + 16), y, cw, ch }, "Not in 1000 Genomes", fmtInt(atoll(dbMeta("n_novel_kg"))), "rare, private or artefacts", BLANK);
    qcCard((Rectangle){ x + 3 * (cw + 16), y, cw, ch }, "Mean depth", atof(dbMeta("autosome_mean_depth")) > 0 ? sstr("%.1fx", atof(dbMeta("autosome_mean_depth"))) : "-",
           atof(dbMeta("autosome_mean_depth")) <= 0 ? "the file gives no depth" :
           *dbMeta("covbin_depth_source") ? "reads covering an average autosomal base" : "minimum depth per gVCF reference block", BLANK);
    y += ch + 28;
    }
    /* reads */
    if (g_rwQc && g_rwQc->nrow) {
        gxSectionTitle(x, y, "SEQUENCING READS (FASTQ)");
        y += 24;
        long long reads = 0, bases = 0;
        for (int m = 0; m < g_rwQc->nrow; m++) { reads += ri(g_rwQc, m, 1); bases += ri(g_rwQc, m, 2); }
        qcCard((Rectangle){ x, y, cw, ch }, "Read pairs", fmtBig((double)ri(g_rwQc, 0, 1)), sstr("%s reads in all", fmtInt(reads)), BLANK);
        qcCard((Rectangle){ x + cw + 16, y, cw, ch }, "Bases sequenced", sstr("%.1f Gb", (double)bases / 1e9), sstr("about %.0fx of the genome", (double)bases / 3.1e9), BLANK);
        qcCard((Rectangle){ x + 2 * (cw + 16), y, cw, ch }, "Q30 or better", sstr("%.1f%% / %.1f%%", rf(g_rwQc, 0, 3) * 100, g_rwQc->nrow > 1 ? rf(g_rwQc, 1, 3) * 100 : 0),
               "read 1 / read 2: 1 error in 1,000 or better", rf(g_rwQc, 0, 3) > 0.85 ? TH->ok : TH->warn);
        qcCard((Rectangle){ x + 3 * (cw + 16), y, cw, ch }, "Duplication (est.)", sstr("%.1f%%", rf(g_rwQc, 0, 8) * 100), sstr("adapter in %.2f%% of reads", rf(g_rwQc, 0, 7) * 100), BLANK);
        y += ch + 16;
        float hw = (W - 16) / 2;
        Rectangle q = { x, y, hw, 230 };
        gxCard(q);
        gxCardTitle(q, "Base quality along the read", "mean Phred score per position (read 1 blue, read 2 amber)");
        lineChart((Rectangle){ q.x + 44, q.y + 70, q.width - 64, 120 }, g_rwPosQ[0], g_rwPosN[0], 30, 40, TH->accent, NULL, "position in read");
        if (g_rwPosN[1]) {   /* read 2 over read 1, without a second grid */
            Rectangle a = { q.x + 44, q.y + 70, q.width - 64, 120 };
            Vector2 prev = { 0 };
            for (int i = 0; i < g_rwPosN[1]; i++) {
                Vector2 p2 = { a.x + a.width * (float)i / (float)(g_rwPosN[1] - 1), a.y + a.height * (1 - clamp01((g_rwPosQ[1][i] - 30) / 10)) };
                if (i) drLine(prev.x, prev.y, p2.x, p2.y, 2, TH->warn);
                prev = p2;
            }
        }
        Rectangle g = { x + hw + 16, y, hw, 230 };
        gxCard(g);
        gxCardTitle(g, "GC content of reads", sstr("mean %.1f%%; a single smooth peak near 40%% is what a clean human library looks like", rf(g_rwQc, 0, 5) * 100));
        float gc2[28];
        for (int k = 0; k < 28; k++) gc2[k] = g_rwGc[0][14 + 2 * k] + g_rwGc[0][15 + 2 * k];
        histogram((Rectangle){ g.x + 20, g.y + 70, g.width - 40, 140 }, gc2, 28, 14, 70, NAN, TH->accent, "GC % of each read (14 to 70, 2% bins)");
        y += 246;
    }
    /* alignment */
    if (g_rwAln && g_rwAln->nrow) {
        gxSectionTitle(x, y, "ALIGNMENT (READS ON GRCh38)");
        y += 24;
        Rectangle c = { x, y, W, 40 + 26 * (float)g_rwAln->nrow };
        gxCard(c);
        float yy = c.y + 18;
        for (int i = 0; i < g_rwAln->nrow; i++) {
            txt(F_REG, TH->fsSmall, c.x + 20, yy, TH->dim, rc(g_rwAln, i, 0));
            txt(F_SEMI, TH->fsSmall, c.x + 300, yy, TH->text, rc(g_rwAln, i, 1));
            txt(F_REG, TH->fsSmall, c.x + 470, yy, TH->faint, rc(g_rwAln, i, 2));
            yy += 26;
        }
        y += c.height + 20;
    }
    /* files */
    if (g_rwFiles) {
        gxSectionTitle(x, y, "FILES");
        y += 24;
        for (int i = 0; i < g_rwFiles->nrow; i++) {
            Rectangle row = { x, y, W, 50 };
            drRound(row, 8, (i & 1) ? TH->zebra : TH->card);
            txtIn(F_SEMI, TH->fsSmall, (Rectangle){ row.x + 16, row.y + 6, W - 260, 18 }, TH->text, 0, rc(g_rwFiles, i, 0));
            txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ row.x + 16, row.y + 26, W - 260, 18 }, TH->faint, 0, rc(g_rwFiles, i, 3));
            txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + W - 240, row.y + 6, 110, 18 }, TH->dim, 0, rc(g_rwFiles, i, 2));
            txtIn(F_SEMI, TH->fsSmall, (Rectangle){ row.x + W - 120, row.y + 6, 104, 18 }, TH->text, 2, fmtBp(atof(rc(g_rwFiles, i, 1))) + 0);
            y += 54;
        }
    }
    y += 30;
    gxScrollEnd(&g_rwScroll, r, y - top);
}

/* ---- sources and notes ---- */
static GxScroll g_abScroll;
static void pageAbout(Rectangle r) {
    float top = gxScrollBegin(&g_abScroll, r);
    float x = r.x + 32, W = fminf(1000, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Sources and notes", NULL);
    if (gxButtonEx(gxId("about export"), (Rectangle){ x, y - 8, 260, 32 }, "Summary for a doctor (HTML)", BTN_PRIMARY, true,
                   "A printable summary of medication genes, HLA type and clinical findings, written into your Documents folder"))
        exportReport();
    y += 44;
    struct { const char *name, *what; } src[] = {
        { "This person's data", sstr("%s%s, converted to GRCh38 and Ensembl chromosome names by the import.", dkName(),
                                      !strcmp(dbMeta("source_build"), "37") ? ", lifted from GRCh37 with the UCSC hg19-to-hg38 chain" : "") },
        { "Reference genome", dkReads() ? "GRCh38 primary assembly, Ensembl release 115. Reads aligned with minimap2 (short-read preset), duplicates marked with samtools."
                                        : "GRCh38 primary assembly, Ensembl release 115." },
        { "Genes", "Ensembl 115 gene models (canonical and MANE Select transcripts); consequences called with bcftools csq." },
        { "ClinVar", sstr("NCBI ClinVar (GRCh38 VCF), weekly release used at import. %s", "Classifications change: re-run the import to refresh them.") },
        { "Population frequencies", "1000 Genomes Project, 30x high-coverage panel (NYGC), 2,590 unrelated founders from 26 populations." },
        { "Missense predictions", "AlphaMissense (Cheng et al. 2023, Science), precomputed hg38 scores." },
        { "Gene function and constraint", "UniProt (reviewed human entries); gnomAD v4.1 constraint metrics; HPO gene-disease links with inheritance." },
        { "Secondary findings", "ACMG SF v3.3 gene list (Lee et al. 2025, Genetics in Medicine)." },
        { "Medications", sstr("%s, with CPIC and DPWG guidelines and FDA label annotations; PharmGKB clinical annotations.", dbMeta("pgx_tool")) },
        { "Research associations", "NHGRI-EBI GWAS Catalog (associations at p < 5e-8)." },
        { "Lineages", sstr("%s%s%s.", dbMeta("mt_tree"), *dbMeta("mt_tree") && *dbMeta("y_tree") ? "; " : "", dbMeta("y_tree")) },
        { "Traits", "A curated panel of replicated single-variant traits, each allele checked against 1000 Genomes frequencies (strand guard)." },
    };
    for (size_t i = 0; i < NELEMS(src); i++) {
        float h = txtParaH(F_REG, TH->fsSmall, W - 260, 19, src[i].what, 4);
        txt(F_SEMI, TH->fsSmall, x, y, TH->text, src[i].name);
        txtPara(F_REG, TH->fsSmall, x + 240, y, W - 260, 19, TH->dim, src[i].what, 4);
        y += fmaxf(22, h) + 12;
    }
    y += 16;
    Rectangle c = { x, y, W, 0 };
    static char note[1400];   /* longer than an sstr slot */
    snprintf(note, sizeof note, "Strand Open runs entirely on this computer. Nothing about anyone's genome leaves it: the setup downloads public reference data, "
                       "and the only outside requests happen when you click a link to a public database, carrying the public identifier you clicked.\n\n"
                       "This is a tool for curiosity and for preparing questions, not a medical test. %s Classifications in ClinVar are revised; risk "
                       "estimates from common variants are averages over populations. Anything that might matter for health deserves a clinical-grade "
                       "confirmation and a conversation with a doctor or genetic counsellor.",
                       dkChip() ? "Chip calls can be wrong, above all for rare variants, most of which turn out false when checked by sequencing."
                                : "Calls from short-read sequencing can be wrong, especially in repetitive regions and for indels.");
    c.height = txtParaH(F_REG, TH->fs, W - 40, 21, note, 30) + 40;
    gxCard(c);
    txtPara(F_REG, TH->fs, c.x + 20, c.y + 20, W - 40, 21, TH->text, note, 30);
    y += c.height + 30;
    gxSectionTitle(x, y, "KEYS");
    y += 24;
    const char *keys[][2] = { { "Ctrl+K", "search genes, rsIDs, positions, conditions" }, { "Ctrl+Q", "close (or the button at the top right)" }, { "Alt+Left / Alt+Right", "back / forward (also mouse side buttons)" },
                              { "Ctrl+plus / minus / 0", "interface size" }, { "Browser: drag, wheel", "pan, zoom around the pointer" },
                              { "Browser: double-click", "zoom in (on a gene: its page)" }, { "Browser: arrows, plus, minus", "pan and zoom" },
                              { "Browser: Shift+wheel", "scroll the reads" } };
    for (size_t i = 0; i < NELEMS(keys); i++) {
        txt(F_MONO, 12, x, y, TH->accent, keys[i][0]);
        txt(F_REG, TH->fsSmall, x + 240, y, TH->dim, keys[i][1]);
        y += 22;
    }
    y += 30;
    txt(F_REG, TH->fsTiny, x, y, TH->faint, sstr("Strand Open %s. Database built %s.", ST_VERSION, dbMeta("built")));
    y += 40;
    gxScrollEnd(&g_abScroll, r, y - top);
}
