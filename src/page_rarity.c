/* ========================================================================
   src/page_rarity.c - your rarest variants, and how rare in the world.

   pipeline/18_rarity.py estimates, for each variant, the share of the
   world's people who carry it (the five 1000 Genomes continental groups
   weighted by population). Three lists: the rarest changes to proteins,
   rare two-copy genotypes, and variants rare in your own ancestry group
   but common elsewhere. Each row: what the gene does, the change, a
   rarity scale and where in the world the variant turns up.
   ======================================================================== */

static int g_ryTab;
static GxScroll g_ryScroll;
static Rows *g_ryRows[3];
static const char *const g_ryKinds[3] = { "protein", "homo", "ancestry" };

static const char *oneIn(double x) {
    if (isnan(x) || x <= 0) return "fewer than 1 in 2,590";
    if (x < 1e-5) return "fewer than 1 in 100,000";
    double n = 1.0 / x;
    if (n < 1.5) return "nearly everyone";
    if (n < 100) return sstr("about 1 in %.0f", n);
    double m = pow(10, floor(log10(n)) - 1);
    return sstr("about 1 in %s", fmtInt((long long)(round(n / m) * m)));
}
static const char *pctPeople(double x) {
    if (isnan(x) || x <= 0) return "under 0.04%";
    if (x < 0.0001) return "under 0.01%";
    if (x < 0.01) return sstr("%.2f%%", x * 100);
    return sstr("%.1f%%", x * 100);
}

/* where a variant turns up: five small bars of allele frequency */
static void ryWorld(Rectangle r, Rows *rows, int i) {
    static const char *const lab[5] = { "AFR", "EUR", "EAS", "SAS", "AMR" };
    float bw = (r.width - 4 * 6) / 5;
    for (int k = 0; k < 5; k++) {
        double f = rf(rows, i, 9 + k);
        float x = r.x + (bw + 6) * (float)k;
        Rectangle b = { x, r.y, bw, 6 };
        drRound(b, 3, TH->field);
        if (!isnan(f) && f > 0) {
            float w = fmaxf(3, bw * (float)fmin(1.0, (log10(f) + 4) / 4));
            drRound((Rectangle){ x, r.y, w, 6 }, 3, f >= 0.05 ? TH->accent : f >= 0.01 ? TH->accent2 : TH->warn);
        }
        txtIn(F_REG, TH->fsTiny - 1, (Rectangle){ x, r.y + 8, bw, 14 }, TH->faint, 1, lab[k]);
    }
}

/* a log scale from 1 in a million to everyone, with the variant's place */
static void ryScale(Rectangle r, double world, bool unseen) {
    drRound((Rectangle){ r.x, r.y + r.height / 2 - 2, r.width, 4 }, 2, TH->field);
    float t = unseen ? 0.12f : (float)clamp01((log10(fmax(world, 1e-6)) + 6) / 6);
    Color c = unseen || world < 0.001 ? TH->bad : world < 0.01 ? TH->warn : world < 0.05 ? TH->accent2 : TH->accent;
    drRound((Rectangle){ r.x, r.y + r.height / 2 - 2, r.width * t, 4 }, 2, drFade(c, 0.5f));
    DrawCircleV((Vector2){ r.x + r.width * t, r.y + r.height / 2 }, 6, c);
    DrawCircleV((Vector2){ r.x + r.width * t, r.y + r.height / 2 }, 3, TH->card);
}

static void ryRow(Rectangle row, Rows *r, int i, int kind) {
    bool hov;
    if (gxBehave(gxIdN(52000 + i), row, &hov, NULL)) lkVariantMenu(ri(r, i, 2));
    drRound(row, 8, hov ? TH->hover : TH->card);
    if (hov) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
    float x = row.x + 16, y = row.y + 12;
    txtIn(F_SEMI, TH->fsH3, (Rectangle){ x, y + 6, 48, 24 }, TH->faint, 0, sstr("%lld", ri(r, i, 1) + (kind == 2 ? 0 : 1)));
    x += 52;
    const char *gene = *rc(r, i, 3) ? rc(r, i, 3) : "Between genes";
    txt(F_SEMI, TH->fs + 1, x, y, TH->text, gene);
    float gx = x + (float)txtW(F_SEMI, TH->fs + 1, gene) + 10;
    const char *chg = *rc(r, i, 4) ? rc(r, i, 4) : csqLabel(rc(r, i, 5));
    float rightW = 330;
    txtIn(F_REG, TH->fs, (Rectangle){ gx, y + 1, row.width - (gx - row.x) - rightW - 20, 20 }, TH->dim, 0, chg);
    if (*rc(r, i, 4)) {
        const char *cs = csqLabel(rc(r, i, 5));
        float cx = gx + (float)txtW(F_REG, TH->fs, chg) + 10;
        if (cx + gxPillW(cs, TH->fsTiny) < row.x + row.width - rightW - 20) gxPill(cx, y, cs, TH->purple, TH->fsTiny);
    }
    float ty = y + 24;
    if (*rc(r, i, 6)) { txtIn(F_REG, TH->fsSmall, (Rectangle){ x, ty, row.width - (x - row.x) - rightW - 20, 18 }, TH->text, 0, rc(r, i, 6)); ty += 20; }
    txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ x, ty, row.width - (x - row.x) - rightW - 20, 18 }, TH->faint, 0, rc(r, i, 15));
    /* right: the headline number, scale and world bars */
    float rx = row.x + row.width - rightW;
    bool unseen = rnull(r, i, 14) && kind == 0;
    double world = rf(r, i, 14);
    const char *big = kind == 1 ? oneIn(world) : unseen ? "never seen" : sstr("%s of people", pctPeople(world));
    const char *sub = kind == 1 ? "expected to have two copies" : unseen ? "not in the 1000 Genomes panel" : oneIn(world);
    txtIn(F_SEMI, TH->fs + 1, (Rectangle){ rx, y - 1, rightW - 16, 20 }, unseen || world < 0.001 ? TH->bad : TH->text, 2, big);
    txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ rx, y + 20, rightW - 16, 16 }, TH->dim, 2, sub);
    ryScale((Rectangle){ rx + 20, y + 40, rightW - 36, 14 }, world, unseen);
    ryWorld((Rectangle){ rx + 20, y + 60, rightW - 36, 22 }, r, i);
}

static void ryStat(Rectangle c, const char *label, const char *value, const char *sub, Color tone) {
    gxCard(c);
    drRound((Rectangle){ c.x + 16, c.y + 18, 3, 16 }, 1.5f, tone);
    txt(F_REG, TH->fsSmall, c.x + 26, c.y + 18, TH->dim, label);
    txt(F_SEMI, 26, c.x + 18, c.y + 44, TH->text, value);
    txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 18, c.y + 82, c.width - 30, 16 }, TH->faint, 0, sub);
}

static const char *const g_ryTabs[] = { "Rarest protein changes", "Rare two-copy genotypes", "Unusual for your ancestry" };
static void pageRarity(Rectangle r) {
    /* without 1000 Genomes frequencies step 18 makes the table but leaves it empty */
    if (!dbHasTable("rarity") || !strcmp(dbMeta("rare_has_af"), "0")) {
        pageMissing(r, "Your rarest variants", dkChip()
                    ? "This page needs sequencing data. A genotyping chip reads common, well-known positions, and the rare calls it does make are "
                      "mostly wrong when checked by sequencing (Weedon et al. 2021, BMJ), so how rare your variants are cannot be judged from it."
                    : !strcmp(dbMeta("rare_has_af"), "0")
                    ? "Rarity cannot be judged without the 1000 Genomes reference. Download it from Setup on the People page, then re-import this person."
                    : "The rarity step did not run for this person. Re-import them from the People page.");
        return;
    }
    if ((unsigned)g_ryTab >= NELEMS(g_ryTabs)) g_ryTab = 0;   /* a tab index from --enter may be anything */
    for (int k = 0; k < 3; k++)
        if (!g_ryRows[k])
            g_ryRows[k] = dbq("SELECT kind, rank, vid, gene, change, csq, plain, zyg, af, af_afr, af_eur, af_eas, af_sas, af_amr, world, note "
                              "FROM rarity WHERE kind=? ORDER BY rank", g_ryKinds[k], NULL);
    float top = gxScrollBegin(&g_ryScroll, r);
    float x = r.x + 32, W = fminf(1240, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Your rarest variants",
                         "How rare your variants are around the world. Shares are estimated from the 2,590 people of 1000 Genomes, each continental group "
                         "weighted by how many people it stands for; regions the panel does not sample (the Middle East, North Africa, Central Asia, Oceania) make them rough.");
    Rows *p = g_ryRows[0];
    int unseen = 0;
    for (int i = 0; p && i < p->nrow; i++) unseen += rnull(p, i, 14);
    float cw = (W - 3 * 14) / 4;
    ryStat((Rectangle){ x, y, cw, 106 }, "Protein changes never seen", fmtInt(unseen), "not in the 1000 Genomes panel", TH->bad);
    ryStat((Rectangle){ x + (cw + 14), y, cw, 106 }, "Carried by under 1 in 1,000", fmtInt(atoll(dbMeta("rarity_under_1000"))), "of your variants that the panel has seen", TH->warn);
    ryStat((Rectangle){ x + 2 * (cw + 14), y, cw, 106 }, "Rare two-copy genotypes", fmtInt(atoll(dbMeta("rarity_homo_under_1000"))), "fewer than 1 in 1,000 expected to share them", TH->purple);
    ryStat((Rectangle){ x + 3 * (cw + 14), y, cw, 106 }, sstr("Rare in %s", *dbMeta("rarity_main_group") ? dbMeta("rarity_main_group") : "your group"),
           fmtInt(atoll(dbMeta("rarity_ancestry_n"))), "but carried by 1 in 5 or more somewhere else", TH->accent);
    y += 106 + 18;
    if (gxTabs(gxId("rarity tabs"), x - 14, y, 40, g_ryTabs, 3, &g_ryTab)) {}
    y += 40;
    drHLine(x, y, W, TH->line);
    y += 14;
    const char *const intro[3] = {
        "Changes to the protein a gene makes, rarest first. \"Never seen\" means the variant is missing from the 1000 Genomes panel, which leaves out anything "
        "found in only one of its 3,202 people, so it is probably carried by fewer than about 1 in 1,000 people. "
        "Rare is not the same as harmful; the Rare and remarkable page says how disruptive each change looks.",
        "Variants you carry on both copies of a chromosome, where few people are expected to: the chance of two copies is roughly the square of the variant's frequency. "
        "Only clean single-letter calls in or near genes are listed.",
        sstr("Variants that are rare among %s but common in another part of the world: small traces of a wider family tree, or just chance. "
             "Gene changes first, then the rest.",   /* the meta value is a place: Europe, the Americas... */
             *dbMeta("rarity_main_group") ? sstr("people with roots in %s", dbMeta("rarity_main_group")) : "people of your ancestry"),
    };
    y += txtPara(F_REG, TH->fsSmall, x, y, W, 19, TH->dim, intro[g_ryTab], 3) + 14;
    Rows *rows = g_ryRows[g_ryTab];
    for (int i = 0; rows && i < rows->nrow && i < 200; i++) {
        float h = *rc(rows, i, 6) ? 98 : 90;
        if (y + h > r.y && y < r.y + r.height) ryRow((Rectangle){ x, y, W, h }, rows, i, g_ryTab);
        y += h + 8;
    }
    y += 12;
    y += txtPara(F_REG, TH->fsTiny + 1, x, y, W, 17, TH->faint,
                 sstr("World weights: %s. Frequencies come from 1000 Genomes only; gnomAD (over 800,000 people) gives far sharper numbers and is a click away "
                      "in any variant's panel, which opens it in your browser.", dbMeta("rarity_weights")), 3) + 30;
    gxScrollEnd(&g_ryScroll, r, y - top);
}
