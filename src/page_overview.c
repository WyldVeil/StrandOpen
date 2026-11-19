/* ========================================================================
   src/page_overview.c - the first page: what the genome holds at a glance,
   the karyogram, and the findings worth a look first.
   ======================================================================== */

static Rows *g_ovCsq;          /* consequence counts (one scan of the variant table, cached) */
static Rows *g_ovFind;         /* highlights (pipeline step 17) */
static bool g_ovHl;            /* g_ovFind is the highlight table, not the ClinVar-only fallback */
static bool g_ovAsked;
static GxScroll g_ovScroll;

static bool statCard(Rectangle r, const char *label, const char *value, const char *sub, Color accent) {
    bool h;
    bool clicked = gxBehave(gxId(label), r, &h, NULL);
    gxCard(r);
    if (h) { drRoundLine(r, 10, accent); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
    drRound((Rectangle){ r.x + 16, r.y + 18, 4, 18 }, 2, accent);
    txt(F_REG, TH->fsSmall, r.x + 28, r.y + 19, TH->dim, label);
    txt(F_SEMI, 28, r.x + 18, r.y + 46, TH->text, value);
    if (sub) txtIn(F_REG, TH->fsSmall, (Rectangle){ r.x + 18, r.y + 84, r.width - 30, 18 }, TH->faint, 0, sub);
    return clicked;
}

static Color hlTone(const char *t) {
    return !strcmp(t, "bad") ? TH->bad : !strcmp(t, "warn") ? TH->warn : !strcmp(t, "accent") ? TH->accent : !strcmp(t, "ok") ? TH->ok : TH->info;
}
static void ovFindings(Rectangle c) {
    gxCard(c);
    /* the card is narrow on a small window: the subtitle is cut to it, and the
       two buttons go one above the other when they do not fit side by side */
    gxCardTitle(c, "Findings to look at first", NULL);
    txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 18, c.y + 16 + (float)TH->fsH3 + 6, c.width - 36, (float)TH->fsSmall + 2 }, TH->dim, 0,
          g_ovHl ? "Health, medicines, rare variants, scores" : "ClinVar pathogenic calls, strongest first");
    const char *l1 = "All clinical findings", *l2 = "Rare and remarkable";
    float w1 = fminf(gxButtonW(l1), c.width - 36), w2 = fminf(gxButtonW(l2), c.width - 36);
    bool stack = 18 + w1 + 8 + w2 + 18 > c.width;
    float by = c.y + c.height - (stack ? 78 : 42);
    float y = c.y + 70;
    if (!g_ovFind) { txt(F_REG, TH->fsSmall, c.x + 18, y, TH->faint, "Loading..."); return; }
    /* without the highlight table (step 17) only ClinVar calls are listed: say so above the buttons */
    float rowsEnd = by - 4;
    if (!g_ovHl) {
        const char *msg = "The full list of highlights is missing. Re-import this person from the People page.";
        float nh = txtParaH(F_REG, TH->fsSmall - 1, c.width - 36, 17, msg, 4);
        rowsEnd = by - 8 - nh;
        txtPara(F_REG, TH->fsSmall - 1, c.x + 18, rowsEnd + 2, c.width - 36, 17, TH->warn, msg, 4);
        rowsEnd -= 4;
    }
    if (!g_ovFind->nrow) txt(F_REG, TH->fsSmall, c.x + 18, y, TH->dim, "Nothing stands out.");
    for (int i = 0; i < g_ovFind->nrow && y + 48 < rowsEnd; i++) {
        Rectangle row = { c.x + 10, y, c.width - 20, 46 };
        bool h;
        if (gxBehave(gxIdN(9100 + i), row, &h, NULL)) stEnter(rc(g_ovFind, i, 5));
        if (h) { drRound(row, 6, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
        Color cc = hlTone(rc(g_ovFind, i, 1));
        DrawCircleV((Vector2){ row.x + 14, row.y + 15 }, 4, cc);
        const char *badge = rc(g_ovFind, i, 4);
        float bw = gxPillW(badge, TH->fsTiny);
        txtIn(F_SEMI, TH->fs, (Rectangle){ row.x + 28, row.y + 5, row.width - 48 - bw, 20 }, TH->text, 0, rc(g_ovFind, i, 2));
        gxPill(row.x + row.width - bw - 8, row.y + 4, badge, cc, TH->fsTiny);
        txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ row.x + 28, row.y + 25, row.width - 40, 18 }, TH->faint, 0, rc(g_ovFind, i, 3));
        if (h) gxTip(gxIdN(9300 + i), row, rc(g_ovFind, i, 3));
        y += 48;
    }
    if (gxButtonEx(gxId("all findings"), (Rectangle){ c.x + 18, by, w1, 30 }, l1, BTN_GHOST, true, NULL))
        navPage(PG_HEALTH);
    if (gxButtonEx(gxId("rare findings"), stack ? (Rectangle){ c.x + 18, by + 36, w2, 30 } : (Rectangle){ c.x + 18 + w1 + 8, by, w2, 30 }, l2, BTN_GHOST, true, NULL))
        navPage(PG_RARE);
}

static void ovCsq(Rectangle c) {
    gxCard(c);
    gxCardTitle(c, "What your variants do", "Most severe effect on the gene's main transcript");
    if (!g_ovCsq) { txt(F_REG, TH->fsSmall, c.x + 18, c.y + 70, TH->faint, "Counting..."); return; }
    static const struct { const char *k, *label; int impact; } show[] = {
        { "stop_gained", "Stop gained", 3 }, { "frameshift", "Frameshift", 3 }, { "splice_donor", "Splice donor", 3 },
        { "splice_acceptor", "Splice acceptor", 3 }, { "start_lost", "Start lost", 3 }, { "stop_lost", "Stop lost", 3 },
        { "inframe_deletion", "In-frame deletion", 2 }, { "inframe_insertion", "In-frame insertion", 2 },
        { "missense", "Missense", 2 }, { "splice_region", "Splice region", 1 }, { "synonymous", "Synonymous", 1 },
        { "5_prime_utr", "5' UTR", 0 }, { "3_prime_utr", "3' UTR", 0 }, { "non_coding", "Non-coding RNA", 0 },
        { "intron", "Intron", 0 }, { "", "Intergenic", 0 },
    };
    long long mx = 1;
    for (int i = 0; i < g_ovCsq->nrow; i++) mx = ri(g_ovCsq, i, 1) > mx ? ri(g_ovCsq, i, 1) : mx;
    float y = c.y + 66;
    for (size_t k = 0; k < NELEMS(show); k++) {
        long long n = 0;
        for (int i = 0; i < g_ovCsq->nrow; i++) if (!strcmp(rc(g_ovCsq, i, 0), show[k].k)) n = ri(g_ovCsq, i, 1);
        if (!n) continue;
        float v = (float)(log10((double)n + 1) / log10((double)mx + 1));
        labelBar((Rectangle){ c.x + 18, y, c.width - 36, 22 }, show[k].label, v, 1, impactColor(show[k].impact), fmtInt(n));
        y += 24;
        if (y > c.y + c.height - 30) break;
    }
    txt(F_REG, TH->fsTiny, c.x + 18, c.y + c.height - 22, TH->faint, "Bar lengths on a log scale.");
}

static void ovMeds(Rectangle c) {
    gxCard(c);
    gxCardTitle(c, "Medicines", "Genes that change how you process some drugs");
    if (!dbHasTable("pgx_gene")) { txt(F_REG, TH->fsSmall, c.x + 18, c.y + 70, TH->faint, "Run the medications step."); return; }
    static Rows *r;
    if (!r) r = dbq("SELECT gene, diplotype, phenotype FROM pgx_gene WHERE called=1 AND phenotype!='' AND phenotype NOT LIKE 'Normal%' "
                    "AND phenotype NOT LIKE '%Uncertain%' AND phenotype NOT LIKE 'normal%' AND phenotype NOT LIKE '%negative%' ORDER BY gene", NULL);
    float y = c.y + 68;
    for (int i = 0; i < r->nrow && y < c.y + c.height - 50; i++) {
        txt(F_SEMI, TH->fs, c.x + 18, y, TH->text, rc(r, i, 0));
        txtIn(F_MONO, 12, (Rectangle){ c.x + 110, y, 150, 18 }, TH->dim, 0, rc(r, i, 1));
        txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 270, y, c.width - 290, 18 }, phenoColor(rc(r, i, 2)), 0, rc(r, i, 2));
        y += 26;
    }
    if (gxButtonEx(gxId("ov meds"), (Rectangle){ c.x + 18, c.y + c.height - 42, 140, 30 }, "All medicines", BTN_GHOST, true, NULL)) navPage(PG_PGX);
}
static void ovTraits(Rectangle c) {
    gxCard(c);
    gxCardTitle(c, "Traits at a glance", NULL);
    if (!dbHasTable("trait")) return;
    static Rows *r;
    if (!r) r = dbq("SELECT title, result FROM trait WHERE id IN ('eye','redhair','lactase','caffeine','earwax','bitter','abo','actn3','flush','secretor') ORDER BY ord", NULL);
    float y = c.y + 52, half = (c.width - 36) / 2;
    for (int i = 0; i < r->nrow; i++) {
        float x = c.x + 18 + (float)(i % 2) * half, yy = y + (float)(i / 2) * 46;
        txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ x, yy, half - 12, 16 }, TH->dim, 0, rc(r, i, 0));
        txtIn(F_SEMI, TH->fs, (Rectangle){ x, yy + 17, half - 12, 20 }, TH->text, 0, rc(r, i, 1));
    }
    if (gxButtonEx(gxId("ov traits"), (Rectangle){ c.x + 18, c.y + c.height - 42, 120, 30 }, "All traits", BTN_GHOST, true, NULL)) navPage(PG_TRAITS);
}
static void ovPrs(Rectangle c) {
    gxCard(c);
    gxCardTitle(c, "Polygenic scores furthest from average",
                sstr("Health scores; percentile among %s", *dbMeta("prs_group_label") ? dbMeta("prs_group_label") : "1000 Genomes Europeans"));
    if (!dbHasTable("prs")) { txt(F_REG, TH->fsSmall, c.x + 18, c.y + 70, TH->faint, "Appears once the polygenic step has run."); return; }
    static Rows *r;
    if (!r) r = dbq(dbColumnExists("prs", "direction") ? "SELECT trait, pct_eur, direction FROM prs WHERE category != 'Mind and behaviour' ORDER BY abs(pct_eur - 50) DESC LIMIT 8"
                    : dbColumnExists("prs", "note") ? "SELECT trait, pct_eur, 1 FROM prs WHERE category != 'Mind and behaviour' ORDER BY abs(pct_eur - 50) DESC LIMIT 8"
                                                    : "SELECT trait, pct_eur, 1 FROM prs ORDER BY abs(pct_eur - 50) DESC LIMIT 8", NULL);
    float y = c.y + 68;
    for (int i = 0; i < r->nrow && y < c.y + c.height - 50; i++) {
        double p = rf(r, i, 1);
        int dir = (int)ri(r, i, 2);
        double risk = dir < 0 ? 100 - p : p;
        Color col = dir == 0 ? (p >= 90 || p <= 10 ? TH->purple : TH->accent)
                             : risk >= 90 ? TH->bad : risk >= 75 ? TH->warn : risk <= 10 ? TH->ok : TH->accent;
        txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 18, y, c.width * 0.45f, 20 }, TH->text, 0, rc(r, i, 0));
        Rectangle b = { c.x + c.width * 0.48f, y + 6, c.width * 0.38f, 8 };
        drRound(b, 4, TH->field);
        drRect((Rectangle){ b.x + b.width / 2 - 0.5f, b.y - 3, 1, b.height + 6 }, TH->lineHi);
        DrawCircleV((Vector2){ b.x + b.width * (float)p / 100, b.y + 4 }, 6, col);
        txtIn(F_SEMI, TH->fsSmall, (Rectangle){ c.x + c.width - 92, y, 74, 20 }, col, 2, pctLabel(p));
        y += 30;
    }
    if (gxButtonEx(gxId("ov prs"), (Rectangle){ c.x + 18, c.y + c.height - 42, 150, 30 }, "All scores", BTN_GHOST, true, NULL)) navPage(PG_PRS);
}
static void ovLineage(Rectangle c) {
    gxCard(c);
    gxCardTitle(c, "Lineage and ancestry", NULL);
    float y = c.y + 52;
    const char *mt = dbMeta("mt_haplogroup"), *yh = dbMeta("y_haplogroup"), *anc = dbMeta("ancestry_summary");
    txt(F_REG, TH->fsSmall, c.x + 18, y, TH->dim, "Maternal line (mtDNA)");
    /* a line can be missing from the data: no MT calls in the file, or no Y chromosome at all */
    bool female = !strcmp(dbMeta("sex"), "female");
    if (*mt) txt(F_SEMI, TH->fsH2, c.x + 18, y + 18, TH->text, mt);
    else txt(F_REG, TH->fs, c.x + 18, y + 22, TH->faint, dbHasTable("hap_mt") ? "not in this data (no mitochondrial calls)" : "not yet computed");
    y += 56;
    txt(F_REG, TH->fsSmall, c.x + 18, y, TH->dim, "Paternal line (Y chromosome)");
    if (*yh) txt(F_SEMI, TH->fsH2, c.x + 18, y + 18, TH->text, yh);
    else txt(F_REG, TH->fs, c.x + 18, y + 22, TH->faint, female ? "none: a woman has no Y chromosome (a father's or brother's DNA shows it)"
                                                        : dbHasTable("hap_y") ? "not in this data (no Y chromosome calls)" : "not yet computed");
    y += 56;
    txt(F_REG, TH->fsSmall, c.x + 18, y, TH->dim, "Ancestry");
    txtPara(F_SEMI, TH->fs, c.x + 18, y + 18, c.width - 36, 20, TH->text, *anc ? anc : "not yet computed", 2);
}

static void pageOverview(Rectangle r) {
    if (!g_ovAsked) {
        g_ovAsked = true;
        if (dbHasTable("csq_count")) dbAsync(SLOT_PAGE, "SELECT csq, n FROM csq_count", 0, NULL);
        else dbAsync(SLOT_PAGE, "SELECT coalesce(csq,''), count(*) FROM variant GROUP BY csq", 0, NULL);
        g_ovHl = dbHasTable("highlight");
        if (g_ovHl) dbAsync(SLOT_PAGE2, "SELECT ord, tone, title, detail, badge, target FROM highlight ORDER BY ord", 0, NULL);
        else dbAsync(SLOT_PAGE2, "SELECT id, 'bad', coalesce(gene, 'intergenic') || coalesce(' ' || hgvs_p, ''), replace(clin_dn, '|', '; '), 'ClinVar', 'health:0' FROM variant "
                                 "WHERE clin IN (5,6) ORDER BY clin DESC, clin_stars DESC LIMIT 12", 0, NULL);
    }
    Rows *t;
    if ((t = dbTake(SLOT_PAGE))) { rowsFree(&g_ovCsq); g_ovCsq = t; }
    if ((t = dbTake(SLOT_PAGE2))) { rowsFree(&g_ovFind); g_ovFind = t; }

    float top = gxScrollBegin(&g_ovScroll, r);
    Rectangle rr = { r.x, top, r.width, r.height };
    char sub[600];
    bool chip = dkChip();
    if (chip)
        snprintf(sub, sizeof sub, "%s raw data, %s positions tested and converted to GRCh38; %s of them differ from the reference. "
                 "A chip reads a fixed set of positions chosen because they vary between people (about one letter in 5,000), so it sees "
                 "well-known variants but not rare ones.", dkName(), fmtInt(atoll(dbMeta("callable_bases"))), fmtInt(atoll(dbMeta("n_variants"))));
    else if (!dkGvcf())
        snprintf(sub, sizeof sub, "%s, on GRCh38. %s variant calls.%s", dkName(), fmtInt(atoll(dbMeta("n_variants"))),
                 atoi(dbMeta("ingest_ref_assumed")) ? " The file lists differences only, so positions it leaves out are taken to match the reference "
                                                      "wherever short reads can be trusted (the 1000 Genomes accessibility mask)." : "");
    else
        snprintf(sub, sizeof sub, "%s, on GRCh38. %s variant calls; %s%% of the autosomes confidently called.", dkName(),
                 fmtInt(atoll(dbMeta("n_variants"))), sstr("%.1f", atof(dbMeta("autosome_callable_fraction")) * 100));
    float y = pageHeader(rr, "Your genome", sub);
    float x0 = r.x + 32, W = r.width - 64, gap = 16;
    int ncard = W > 1200 ? 6 : 3;
    float cw = (W - gap * (float)(ncard - 1)) / (float)ncard, ch = 112;
    struct { const char *l, *v, *s; Color c; } cards[6] = {
        { "Variant calls", fmtBig(atof(dbMeta("n_variants"))), sstr("%s SNVs, %s indels", fmtBig(atof(dbMeta("n_snv"))), fmtBig(atof(dbMeta("n_indel")))), TH->accent },
        { "Zygosity", sstr("%.0f%% het", 100.0 * atof(dbMeta("n_het")) / fmax(1, atof(dbMeta("n_variants")))), sstr("%s homozygous, %s hemizygous", fmtBig(atof(dbMeta("n_hom"))), fmtBig(atof(dbMeta("n_hemi")))), TH->purple },
        { "Protein-altering", fmtInt(atoll(dbMeta("n_missense")) + atoll(dbMeta("n_high"))), sstr("%s missense, %s high impact", fmtInt(atoll(dbMeta("n_missense"))), fmtInt(atoll(dbMeta("n_high")))), TH->warn },
        { "In ClinVar", fmtInt(atoll(dbMeta("n_clinvar"))), "variants with a clinical record", TH->bad },
        { "Callable genome", sstr("%.1f%%", atof(dbMeta("autosome_callable_fraction")) * 100), "of autosomal bases called", TH->ok },
        { "Ti/Tv ratio", dbMeta("titv"), "about 2.0 to 2.1 is typical for WGS", TH->info },
    };
    if (chip) {
        cards[0].l = "Positions tested";
        cards[0].v = fmtBig(atof(dbMeta("callable_bases")));
        cards[0].s = sstr("%s carry a variant", fmtBig(atof(dbMeta("n_variants"))));
        cards[4].l = "Converted from";
        cards[4].v = sstr("GRCh%s", *dbMeta("source_build") ? dbMeta("source_build") : "38");
        cards[4].s = !strcmp(dbMeta("source_build"), "37") ? "positions lifted to GRCh38" : "already on GRCh38";
        cards[5].l = "Not read";
        cards[5].v = fmtInt(atoll(dbMeta("ingest_dropped_nocall")) + atoll(dbMeta("ingest_dropped_indel")));
        cards[5].s = "no-calls and insertion/deletion probes";
    } else if (!dkGvcf()) {
        cards[4].l = "Reference calls";
        cards[4].v = atoi(dbMeta("ingest_ref_assumed")) ? "assumed" : "not given";
        cards[4].s = atoi(dbMeta("ingest_ref_assumed")) ? "where the file lists nothing" : "unlisted positions count as untested";
    }
    for (int i = 0; i < 6; i++) {
        int col = i % ncard, row = i / ncard;
        static const StPage dest[6] = { PG_VARIANTS, PG_VARIANTS, PG_GENE, PG_HEALTH, PG_RAW, PG_RAW };
        if (statCard((Rectangle){ x0 + (cw + gap) * (float)col, y + (ch + gap) * (float)row, cw, ch }, cards[i].l, cards[i].v, cards[i].s, cards[i].c)) {
            if (i == 2) variantsShow(2, "", -1);
            else navPage(dest[i]);
        }
    }
    y += (ch + gap) * (float)((6 + ncard - 1) / ncard) + 8;

    /* karyogram + the findings column */
    float rightW = fminf(460, W * 0.36f);
    Rectangle kc = { x0, y, W - rightW - gap, 470 };
    gxCard(kc);
    gxCardTitle(kc, "Your chromosomes", dkChip() ? "Each chromosome's bands, with your variant density beside it (grey: no positions tested)"
                                                 : "Each chromosome's bands, with your variant density beside it (grey: little or no coverage)");
    int clicked = karyogram((Rectangle){ kc.x + 16, kc.y + 68, kc.width - 32, kc.height - 80 }, true);
    if (clicked) navLocus(clicked, 1, g_chrom[clicked].len);
    ovFindings((Rectangle){ x0 + kc.width + gap, y, rightW, 470 });
    y += 470 + gap;
    float half = (W - gap) / 2;
    ovCsq((Rectangle){ x0, y, half, 470 });
    ovLineage((Rectangle){ x0 + half + gap, y, half, 230 });
    ovMeds((Rectangle){ x0 + half + gap, y + 230 + gap, half, 470 - 230 - gap });
    y += 470 + gap;
    ovTraits((Rectangle){ x0, y, half, 330 });
    ovPrs((Rectangle){ x0 + half + gap, y, half, 330 });
    y += 330 + gap;
    Rectangle tips = { x0, y, W, 96 };
    gxCard(tips);
    gxCardTitle(tips, "Getting around", NULL);
    txtPara(F_REG, TH->fsSmall, tips.x + 18, tips.y + 48, tips.width - 36, 19, TH->dim,
            "Ctrl+K searches genes, rsIDs, positions, conditions, medicines and traits. Alt+Left and Alt+Right go back and forward, "
            "as do the mouse's side buttons. Ctrl+plus and Ctrl+minus change the interface size. "
            "In the genome browser, drag to pan, scroll to zoom and double-click to zoom in on a spot.", 4);
    y += 96 + 32;
    gxScrollEnd(&g_ovScroll, r, y - top);
}
