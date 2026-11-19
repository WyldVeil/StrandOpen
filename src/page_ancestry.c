/* ========================================================================
   src/page_ancestry.c - ancestry (composition, PCA, chromosome painting,
   closest populations) and the polygenic scores.
   ======================================================================== */

static Rows *g_anAdmix, *g_anSuper, *g_anPca[2], *g_anPaint, *g_anClose;
static int g_anPcaSet;
static GxScroll g_anScroll;
static Color superColor(const char *sp) {
    if (!strcmp(sp, "EUR")) return (Color){ 110, 178, 206, 255 };
    if (!strcmp(sp, "AFR")) return (Color){ 232, 196, 110, 255 };
    if (!strcmp(sp, "EAS")) return (Color){ 143, 192, 120, 255 };
    if (!strcmp(sp, "SAS")) return (Color){ 180, 142, 200, 255 };
    if (!strcmp(sp, "AMR")) return (Color){ 220, 104, 112, 255 };
    return TH->text;
}
static Color eurPopColor(const char *p) {
    if (!strcmp(p, "GBR")) return (Color){ 110, 178, 206, 255 };
    if (!strcmp(p, "CEU")) return (Color){ 129, 161, 230, 255 };
    if (!strcmp(p, "FIN")) return (Color){ 143, 192, 120, 255 };
    if (!strcmp(p, "IBS")) return (Color){ 232, 196, 110, 255 };
    if (!strcmp(p, "TSI")) return (Color){ 220, 104, 112, 255 };
    return TH->faint;
}
static void pcaPlot(Rectangle r, Rows *p, bool eur) {
    if (!p || !p->nrow) return;
    float mnx = 1e9f, mxx = -1e9f, mny = 1e9f, mxy = -1e9f;
    for (int i = 0; i < p->nrow; i++) {
        float x = (float)rf(p, i, 4), y = (float)rf(p, i, 5);
        mnx = fminf(mnx, x); mxx = fmaxf(mxx, x); mny = fminf(mny, y); mxy = fmaxf(mxy, y);
    }
    float px = (mxx - mnx) * 0.06f, py = (mxy - mny) * 0.06f;
    mnx -= px; mxx += px; mny -= py; mxy += py;
    drRect(r, (Color){ 22, 25, 31, 255 });
    #define PX(v) (r.x + (((v) - mnx) / (mxx - mnx)) * r.width)
    #define PY(v) (r.y + r.height - (((v) - mny) / (mxy - mny)) * r.height)
    drHLine(r.x, PY(0), r.width, drFade(TH->line, 0.6f));
    drVLine(PX(0), r.y, r.height, drFade(TH->line, 0.6f));
    int you = -1, hover = -1;
    float bd = 6;
    for (int i = 0; i < p->nrow; i++) {
        if (!strcmp(rc(p, i, 1), "YOU")) { you = i; continue; }
        float x = PX((float)rf(p, i, 4)), y = PY((float)rf(p, i, 5));
        Color c = eur ? eurPopColor(rc(p, i, 2)) : superColor(rc(p, i, 3));
        DrawCircleV((Vector2){ x, y }, 2.6f, drFade(c, 0.7f));
        float d = Vector2Distance(g_gx.mouse, (Vector2){ x, y });
        if (d < bd && CheckCollisionPointRec(g_gx.mouse, r)) { bd = d; hover = i; }
    }
    if (you >= 0) {
        float x = PX((float)rf(p, you, 4)), y = PY((float)rf(p, you, 5));
        DrawCircleV((Vector2){ x, y }, 9, drFade(WHITE, 0.25f));
        DrawCircleV((Vector2){ x, y }, 6, WHITE);
        DrawCircleV((Vector2){ x, y }, 3, TH->bad);
        txt(F_SEMI, TH->fsSmall, x + 10, y - 20, TH->text, "You");
    }
    if (hover >= 0) {
        Rows *lab = dbq("SELECT label FROM anc_pop WHERE pop=?", rc(p, hover, 2), NULL);
        gxTip(gxIdN(40000 + hover), (Rectangle){ g_gx.mouse.x - 5, g_gx.mouse.y - 5, 10, 10 }, sstr("%s (%s)\n%s", rc(lab, 0, 0), rc(p, hover, 2), rc(p, hover, 1)));
        rowsFree(&lab);
    }
    #undef PX
    #undef PY
}
/* the 22 autosomes as thin bars with marked stretches (rows: chrom, start, end); returns the height used */
static float segStrips(Rectangle c, Rows *segs, Color col, long long minLen) {
    double maxLen = 0;
    for (int k = 1; k <= 22; k++) maxLen = fmax(maxLen, (double)g_chrom[k].len);
    if (maxLen <= 0) return 0;
    float colW = (c.width - 30) / 2, rowH = 15;
    for (int k = 1; k <= 22; k++) {
        int half = (k - 1) / 11, row = (k - 1) % 11;
        float bx = c.x + half * (colW + 30) + 26, by = c.y + row * rowH, bw = colW - 26;
        txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + half * (colW + 30), by - 2, 22, 14 }, TH->faint, 2, chromName(k));
        float lw = bw * (float)((double)g_chrom[k].len / maxLen);
        drRound((Rectangle){ bx, by + 2, lw, 7 }, 3, TH->field);
        for (int i = 0; segs && i < segs->nrow; i++) {
            if (ri(segs, i, 0) != k || ri(segs, i, 2) - ri(segs, i, 1) < minLen) continue;
            float x0 = bx + bw * (float)((double)ri(segs, i, 1) / maxLen), x1 = bx + bw * (float)((double)ri(segs, i, 2) / maxLen);
            drRect((Rectangle){ x0, by + 1, fmaxf(1.5f, x1 - x0), 9 }, col);
        }
    }
    return 11 * rowH + 6;
}

static Rows *g_anRoh, *g_anArch;
static bool g_anExtraAsked;
/* parents' relatedness (runs of homozygosity) and archaic DNA, under the composition cards */
static float ancestryExtras(float x, float y, float W) {
    if (!g_anExtraAsked) {
        g_anExtraAsked = true;
        if (dbHasTable("roh")) g_anRoh = dbq("SELECT chrom, start, end FROM roh WHERE len >= 300000", NULL);
        if (dbHasTable("archaic")) g_anArch = dbq("SELECT chrom, start, end FROM archaic WHERE prob >= 0.8 AND len >= 40000", NULL);
    }
    float y0 = y;
    float cw = (W - 20) / 2;
    float h = 0;
    for (int pass = 0; pass < 2; pass++) {
        /* runs of homozygosity */
        if (*dbMeta("roh_summary")) {
            Rectangle c = { x, y, cw, h };
            float hh = 66 + txtParaH(F_REG, TH->fsSmall, cw - 40, 19, dbMeta("roh_summary"), 6) + 14 + 11 * 15 + 6 + 40;
            if (pass == 1) {
                gxCard(c);
                gxCardTitle(c, "Did your parents share recent ancestors?", "Runs where both copies of your DNA match, inherited from one ancestor on both sides");
                float ty = c.y + 66;
                ty += txtPara(F_REG, TH->fsSmall, c.x + 20, ty, cw - 40, 19, TH->text, dbMeta("roh_summary"), 6) + 14;
                segStrips((Rectangle){ c.x + 20, ty, cw - 40, 0 }, g_anRoh, TH->warn, 300000);
                txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 20, c.y + c.height - 26, cw - 40, 14 }, TH->faint, 0,
                      "Marked: runs over 0.3 Mb. Offspring of first cousins carry about 180 Mb in runs over 5 Mb.");
            }
            h = fmaxf(h, hh);
        }
        /* archaic DNA: measured, or why it cannot be */
        if (!*dbMeta("archaic_frac") && *dbMeta("archaic_unmeasurable")) {
            Rectangle c = { x + cw + 20, y, cw, h };
            float hh = 66 + txtParaH(F_REG, TH->fsSmall, cw - 40, 19, dbMeta("archaic_unmeasurable"), 8) + 24;
            if (pass == 1) {
                gxCard(c);
                gxCardTitle(c, "Neanderthal and other archaic DNA", "Stretches inherited from humans who mixed with modern humans over 40,000 years ago");
                txtPara(F_REG, TH->fsSmall, c.x + 20, c.y + 66, cw - 40, 19, TH->dim, dbMeta("archaic_unmeasurable"), 8);
            }
            h = fmaxf(h, hh);
        } else if (*dbMeta("archaic_frac")) {
            Rectangle c = { x + cw + 20, y, cw, h };
            double frac = atof(dbMeta("archaic_frac")) * 100;
            const char *t1 = sstr("About %.1f%% of your genome sits in stretches inherited from archaic humans: %s stretches, %s Mb in all, the longest %.0f kb.",
                                  frac, fmtInt(atoll(dbMeta("archaic_segments"))), sstr("%.0f", atof(dbMeta("archaic_bp")) / 1e6), atof(dbMeta("archaic_longest")) / 1e3);
            /* a fun comparison with a published average measured by a similar stretch-finding method:
               Sankararaman et al. 2014 (Nature 507:354), confidently inferred Neanderthal ancestry,
               1.15% in Europeans and 1.38% in East Asians of 1000 Genomes */
            static char cmp[400];
            static bool cmpDone;
            if (!cmpDone) {
                cmpDone = true;
                cmp[0] = 0;
                Rows *g = dbHasTable("anc_admix") ? dbq("SELECT pop, frac FROM anc_admix WHERE level='super' ORDER BY frac DESC LIMIT 1", NULL) : NULL;
                double avg = 0;
                const char *who = "";
                if (g && g->nrow && rf(g, 0, 1) >= 0.8) {
                    if (!strcmp(rc(g, 0, 0), "EUR")) { avg = 1.15; who = "Europeans"; }
                    if (!strcmp(rc(g, 0, 0), "EAS")) { avg = 1.38; who = "East Asians"; }
                }
                if (avg > 0) {
                    double rel = (frac - avg) / avg * 100;
                    if (fabs(rel) < 10) snprintf(cmp, sizeof cmp, "About average: the average for %s measured by a similar method is %.2f%% (Sankararaman et al. 2014). ", who, avg);
                    else snprintf(cmp, sizeof cmp, "That is about %.0f%% %s than the average for %s measured by a similar method, %.2f%% (Sankararaman et al. 2014). ",
                                  fabs(rel), rel > 0 ? "more" : "less", who, avg);
                }
                rowsFree(&g);
            }
            t1 = sstr("%s %s", t1, cmp);
            const char *t2 = *dbMeta("archaic_ref")
                ? dbMeta("archaic_ref")
                : "People whose ancestors lived outside Africa carry about 1 to 2 percent Neanderthal DNA. This method finds archaic stretches from patterns of variants missing in Africans, so it cannot tell Neanderthal from Denisovan stretches; in Europeans nearly all are Neanderthal.";
            float hh = 66 + txtParaH(F_REG, TH->fsSmall, cw - 40, 19, t1, 7) + 8 + txtParaH(F_REG, TH->fsSmall - 1, cw - 40, 18, t2, 6) + 14 + 11 * 15 + 6 + 40;
            if (pass == 1) {
                gxCard(c);
                gxCardTitle(c, "Neanderthal and other archaic DNA", "Stretches inherited from humans who mixed with your ancestors over 40,000 years ago");
                float ty = c.y + 66;
                ty += txtPara(F_REG, TH->fsSmall, c.x + 20, ty, cw - 40, 19, TH->text, t1, 7) + 8;
                ty += txtPara(F_REG, TH->fsSmall - 1, c.x + 20, ty, cw - 40, 18, TH->dim, t2, 6) + 14;
                segStrips((Rectangle){ c.x + 20, ty, cw - 40, 0 }, g_anArch, TH->accent, 0);
                txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 20, c.y + c.height - 26, cw - 40, 14 }, TH->faint, 0, dbMeta("archaic_method"));
            }
            h = fmaxf(h, hh);
        }
    }
    return h > 0 ? y - y0 + h + 20 : 0;
}

static void pageAncestry(Rectangle r) {
    if (!dbHasTable("anc_admix")) { pageMissing(r, "Ancestry", "Ancestry needs the 1000 Genomes reference. Download it from Setup on the People page, then re-import this person."); return; }
    if (!g_anAdmix) {
        g_anAdmix = dbq("SELECT a.pop, a.label, a.frac, a.lo, a.hi, p.super FROM anc_admix a JOIN anc_pop p ON p.pop=a.pop WHERE a.level='pop' ORDER BY a.frac DESC", NULL);
        g_anSuper = dbq("SELECT pop, label, frac, lo, hi FROM anc_admix WHERE level='super' ORDER BY frac DESC", NULL);
        g_anPca[0] = dbq("SELECT sett, sample, pop, super, pc1, pc2 FROM anc_pca WHERE sett='global'", NULL);
        g_anPca[1] = dbq("SELECT sett, sample, pop, super, pc1, pc2 FROM anc_pca WHERE sett='eur'", NULL);
        g_anPaint = dbq("SELECT chrom, start, end, a, b FROM anc_paint ORDER BY chrom, start", NULL);
        g_anClose = dbq("SELECT pop, label, score, rank FROM anc_closest ORDER BY rank LIMIT 12", NULL);
    }
    float top = gxScrollBegin(&g_anScroll, r);
    float x = r.x + 32, W = fminf(1400, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Ancestry",
                         sstr("Your genome compared with the 2,590 unrelated people of the 1000 Genomes Project, from 26 populations on five continents, "
                              "using %s independent common SNPs.", fmtInt(atoll(dbMeta("ancestry_panel_snps")))));
    /* composition */
    float leftW = fminf(520, W * 0.4f);
    Rectangle cc = { x, y, leftW, 560 };
    gxCard(cc);
    gxCardTitle(cc, "Composition", "Supervised admixture; bars show 95% bootstrap intervals");
    float v[8];
    Color cols[8];
    int n = 0;
    for (int i = 0; i < g_anSuper->nrow && n < 8; i++) { v[n] = (float)rf(g_anSuper, i, 2); cols[n] = superColor(rc(g_anSuper, i, 0)); n++; }
    Vector2 dc = { cc.x + 110, cc.y + 170 };
    donut(dc, 58, 84, v, cols, n);
    txtIn(F_SEMI, TH->fsH2, (Rectangle){ dc.x - 56, dc.y - 18, 112, 24 }, TH->text, 1, g_anSuper->nrow ? sstr("%.0f%%", rf(g_anSuper, 0, 2) * 100) : "");
    txtIn(F_REG, TH->fsTiny, (Rectangle){ dc.x - 56, dc.y + 6, 112, 16 }, TH->dim, 1, g_anSuper->nrow ? rc(g_anSuper, 0, 1) : "");
    float ly = cc.y + 92;
    for (int i = 0; i < g_anSuper->nrow; i++) {
        if (rf(g_anSuper, i, 2) < 0.001 && i > 0) continue;
        DrawCircleV((Vector2){ cc.x + 228, ly + 9 }, 5, superColor(rc(g_anSuper, i, 0)));
        txt(F_REG, TH->fsSmall, cc.x + 240, ly + 1, TH->text, rc(g_anSuper, i, 1));
        txtIn(F_SEMI, TH->fsSmall, (Rectangle){ cc.x + leftW - 130, ly, 110, 18 }, TH->text, 2,
              sstr("%.1f%%  (%.0f-%.0f)", rf(g_anSuper, i, 2) * 100, rf(g_anSuper, i, 3) * 100, rf(g_anSuper, i, 4) * 100));
        ly += 26;
    }
    float py = cc.y + 290;
    gxSectionTitle(cc.x + 18, py, "BY REFERENCE POPULATION");
    py += 22;
    for (int i = 0; i < g_anAdmix->nrow && i < 9; i++) {
        double f = rf(g_anAdmix, i, 2);
        if (f < 0.005 && i > 3) break;
        Rectangle b = { cc.x + 18, py, leftW - 36, 22 };
        txtIn(F_REG, TH->fsSmall, (Rectangle){ b.x, b.y, 220, 22 }, TH->text, 0, rc(g_anAdmix, i, 1));
        Rectangle bar = { b.x + 228, b.y + 7, b.width - 228 - 60, 8 };
        gxMeter(bar, (float)f, superColor(rc(g_anAdmix, i, 5)));
        float l0 = bar.x + bar.width * (float)rf(g_anAdmix, i, 3), l1 = bar.x + bar.width * (float)rf(g_anAdmix, i, 4);
        drRect((Rectangle){ l0, bar.y + 3, fmaxf(1, l1 - l0), 2 }, TH->text);
        txtIn(F_SEMI, TH->fsSmall, (Rectangle){ bar.x + bar.width + 4, b.y, 56, 22 }, TH->text, 2, sstr("%.0f%%", f * 100));
        py += 26;
    }
    txtPara(F_REG, TH->fsTiny + 1, cc.x + 18, cc.y + cc.height - 60, leftW - 36, 16, TH->faint,
            "1000 Genomes has five European groups (British, Utah Northwest European, Finnish, Spanish, Tuscan): splits between close neighbours "
            "such as British and Northwest European are uncertain, as the intervals show.", 4);
    /* PCA */
    Rectangle pc = { x + leftW + 20, y, W - leftW - 20, 560 };
    gxCard(pc);
    gxCardTitle(pc, g_anPcaSet ? "Within Europe" : "The world's genetic map", g_anPcaSet ? "PCA of the 525 European founders; you projected onto it"
                                                                                          : "Principal components of all 2,590 founders; you projected onto it");
    static const char *const sets[] = { "World", "Europe" };
    gxSegmented(gxId("pca set"), (Rectangle){ pc.x + pc.width - 200, pc.y + 16, 180, 30 }, sets, 2, &g_anPcaSet);
    Rectangle plot = { pc.x + 40, pc.y + 70, pc.width - 64, pc.height - 120 };
    pcaPlot(plot, g_anPca[g_anPcaSet], g_anPcaSet == 1);
    const char *ve = dbMeta(g_anPcaSet ? "pca_eur_var" : "pca_global_var");
    float v1 = 0, v2 = 0;
    sscanf(ve, "%f %f", &v1, &v2);
    txtIn(F_REG, TH->fsTiny, (Rectangle){ plot.x, plot.y + plot.height + 4, plot.width, 16 }, TH->faint, 1, sstr("PC1 (%.1f%% of variation)", v1 * 100));
    txt(F_REG, TH->fsTiny, plot.x - 32, plot.y - 16, TH->faint, sstr("PC2 (%.1f%%)", v2 * 100));
    float lx = pc.x + 40;
    if (g_anPcaSet) {
        static const char *const pops[] = { "GBR", "CEU", "FIN", "IBS", "TSI" }, *const lab[] = { "British", "NW European (Utah)", "Finnish", "Spanish", "Tuscan" };
        for (int i = 0; i < 5; i++) { DrawCircleV((Vector2){ lx, pc.y + pc.height - 22 }, 5, eurPopColor(pops[i])); txt(F_REG, TH->fsTiny + 1, lx + 9, pc.y + pc.height - 29, TH->dim, lab[i]); lx += (float)txtW(F_REG, TH->fsTiny + 1, lab[i]) + 30; }
    } else {
        static const char *const sp[] = { "EUR", "AFR", "EAS", "SAS", "AMR" }, *const lab[] = { "European", "African", "East Asian", "South Asian", "American" };
        for (int i = 0; i < 5; i++) { DrawCircleV((Vector2){ lx, pc.y + pc.height - 22 }, 5, superColor(sp[i])); txt(F_REG, TH->fsTiny + 1, lx + 9, pc.y + pc.height - 29, TH->dim, lab[i]); lx += (float)txtW(F_REG, TH->fsTiny + 1, lab[i]) + 30; }
    }
    y += 580;
    /* painting */
    Rectangle pr = { x, y, W, 120 + 22 * 24 };
    gxCard(pr);
    gxCardTitle(pr, "Your chromosomes, painted", "Each window of your genome assigned the pair of continental groups that best explains it (one per copy)");
    double maxLen = (double)g_chrom[1].len;
    float bx = pr.x + 60, bw = pr.width - 90;
    for (int c = 1; c <= 22; c++) {
        float cy = pr.y + 76 + 22 * (float)(c - 1);
        txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ pr.x + 18, cy, 34, 16 }, TH->dim, 2, chromName(c));
        float cw = (float)(bw * (double)g_chrom[c].len / maxLen);
        drRound((Rectangle){ bx, cy, cw, 16 }, 4, TH->field);
        for (int i = 0; i < g_anPaint->nrow; i++) {
            if (ri(g_anPaint, i, 0) != c) continue;
            float x0 = bx + (float)(bw * (double)ri(g_anPaint, i, 1) / maxLen), x1 = bx + (float)(bw * (double)ri(g_anPaint, i, 2) / maxLen);
            drRect((Rectangle){ x0, cy, x1 - x0, 8 }, superColor(rc(g_anPaint, i, 3)));
            drRect((Rectangle){ x0, cy + 8, x1 - x0, 8 }, superColor(rc(g_anPaint, i, 4)));
        }
    }
    y += pr.height + 20;
    /* closest */
    Rectangle cl = { x, y, W, 80 + 26 * (float)(g_anClose ? g_anClose->nrow : 0) };
    gxCard(cl);
    gxCardTitle(cl, "Closest reference populations", "How well each population's allele frequencies fit your genotypes (relative log-likelihood per SNP; 0 is the best fit)");
    float yy = cl.y + 70;
    double worst = g_anClose && g_anClose->nrow ? rf(g_anClose, g_anClose->nrow - 1, 2) : -1;
    for (int i = 0; g_anClose && i < g_anClose->nrow; i++) {
        double s = rf(g_anClose, i, 2);
        txtIn(F_REG, TH->fsSmall, (Rectangle){ cl.x + 18, yy, 260, 20 }, i < 2 ? TH->text : TH->dim, 0, sstr("%lld. %s", ri(g_anClose, i, 3), rc(g_anClose, i, 1)));
        gxMeter((Rectangle){ cl.x + 290, yy + 6, W - 400, 8 }, (float)(1 - s / (worst < 0 ? worst : -1)), i == 0 ? TH->accent : TH->accent2);
        txtIn(F_REG, TH->fsSmall, (Rectangle){ cl.x + W - 100, yy, 80, 20 }, TH->dim, 2, sstr("%.4f", s));
        yy += 26;
    }
    y += cl.height + 20;
    y += ancestryExtras(x, y, W);
    y += 10;
    gxScrollEnd(&g_anScroll, r, y - top);
}

/* ---- polygenic scores ---- */
static Rows *g_psRows, *g_psDist;
static GxScroll g_psScroll;
/* a percentile as people say it: "37th", "top 1%", "bottom 1%" */
static const char *pctLabel(double pct) {
    if (pct >= 99) return "top 1%";
    if (pct < 1) return "bottom 1%";
    int p = (int)lround(pct);
    int t = p % 100;
    const char *suf = (t >= 11 && t <= 13) ? "th" : p % 10 == 1 ? "st" : p % 10 == 2 ? "nd" : p % 10 == 3 ? "rd" : "th";
    return sstr("%d%s", p, suf);
}
static int g_psSort;
static void pagePrs(Rectangle r) {
    if (!dbHasTable("prs")) { pageMissing(r, "Polygenic scores", "Polygenic scores need the 1000 Genomes reference. Download it from Setup on the People page, then re-import this person."); return; }
    if (!g_psRows) {
        g_psRows = dbq(dbColumnExists("prs", "direction")
                       ? "SELECT id, trait, category, name, n_weights, n_used, n_called, pct_eur, z_eur, citation, method, metric, ancestry, pmid, score, note, direction FROM prs "
                         "ORDER BY category='Mind and behaviour', category, trait"
                       : dbColumnExists("prs", "note")
                       ? "SELECT id, trait, category, name, n_weights, n_used, n_called, pct_eur, z_eur, citation, method, metric, ancestry, pmid, score, note, 1 FROM prs "
                         "ORDER BY category='Mind and behaviour', category, trait"
                       : "SELECT id, trait, category, name, n_weights, n_used, n_called, pct_eur, z_eur, citation, method, metric, ancestry, pmid, score, '', 1 FROM prs "
                         "ORDER BY category, trait", NULL);
        g_psDist = dbq("SELECT id, lo, hi, bins FROM prs_dist", NULL);
    }
    float top = gxScrollBegin(&g_psScroll, r);
    float x = r.x + 32, W = fminf(1300, r.width - 64);
    static char intro[1400];   /* longer than an sstr slot */
    snprintf(intro, sizeof intro,
             "Each score adds up the small effects of up to two million common variants found in large studies. Your score is placed among %s%s. "
             "A high percentile means more of the variants that raise the trait than most people carry; for most scores here that means higher risk, "
             "but for HDL cholesterol, bone density and vitamin D a high score is the favourable end. A score shifts the odds, it does not predict. "
             "Lifestyle, family history and rare variants often matter more.%s%s", dbMeta("prs_reference"),
             !strcmp(dbMeta("prs_adjusted"), "1") ? ", adjusted for where your ancestry sits within Europe (scores drift with ancestry for mostly technical reasons)" : "",
             *dbMeta("prs_portability") ? " " : "", dbMeta("prs_portability"));
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Polygenic scores", intro);
    static const char *const sorts[] = { "By area", "Highest first" };
    gxSegmented(gxId("prs sort"), (Rectangle){ x, y - 6, 260, 30 }, sorts, 2, &g_psSort);
    y += 40;
    int order[64];
    int n = g_psRows->nrow < 64 ? g_psRows->nrow : 64;   /* the cards shown */
    for (int i = 0; i < n; i++) order[i] = i;
    if (g_psSort == 1)
        for (int i = 1; i < n; i++) { int k = order[i], j = i - 1; while (j >= 0 && rf(g_psRows, order[j], 7) < rf(g_psRows, k, 7)) { order[j + 1] = order[j]; j--; } order[j + 1] = k; }
    const char *cat = "";
    float cw = (W - 16) / 2;
    int col = 0;
    float rowH = 0;
    for (int oi = 0; oi < n; oi++) {
        int i = order[oi];
        if (g_psSort == 0 && strcmp(cat, rc(g_psRows, i, 2))) {
            if (col) { y += rowH + 16; col = 0; rowH = 0; }
            cat = rc(g_psRows, i, 2);
            gxSectionTitle(x, y + 4, cat);
            y += 28;
            if (!strcmp(cat, "Mind and behaviour")) {
                y += txtPara(F_REG, TH->fsSmall, x, y + 4, W, 19, TH->dim,
                             "Scores for mind, mood and behaviour, plus two body traits shown with them (thinness and male testosterone). "
                             "Most explain only a few percent of the differences between people, and part of that is tangled up with family and upbringing. "
                             "None says anything certain about one person, so read the note on each card. "
                             "Popular questions that genes cannot answer for one person, such as intelligence or wealth, are covered under Rare and remarkable.", 9) + 14;
            }
        }
        const char *note = rc(g_psRows, i, 15);
        float ch = 180 + (*note ? txtParaH(F_REG, TH->fsTiny + 1, cw - 36, 17, note, 6) + 8 : 0);
        if (col == 0) {   /* the row's height: the taller of its two cards */
            rowH = ch;
            if (oi + 1 < n && (g_psSort == 1 || !strcmp(rc(g_psRows, order[oi + 1], 2), rc(g_psRows, i, 2)))) {
                const char *n2 = rc(g_psRows, order[oi + 1], 15);
                rowH = fmaxf(rowH, 180 + (*n2 ? txtParaH(F_REG, TH->fsTiny + 1, cw - 36, 17, n2, 6) + 8 : 0));
            }
        }
        Rectangle c = { x + (cw + 16) * (float)col, y, cw, rowH };
        gxCard(c);
        double pct = rf(g_psRows, i, 7), z = rf(g_psRows, i, 8);
        /* disease scores: red is higher risk; the Mind and behaviour scores explain little and several have no good
           or bad end, so they all get neutral colours */
        /* colour by which end of the score is the one to watch: disease scores red at the top,
           favourable-at-the-top scores (HDL, bone density, vitamin D) the other way round,
           neutral scores (height, every behaviour score) in neutral colours */
        int dir = (int)ri(g_psRows, i, 16);
        bool mind = !strcmp(rc(g_psRows, i, 2), "Mind and behaviour");
        bool neutral = mind || dir == 0;
        double risk = dir < 0 ? 100 - pct : pct;
        Color pc = neutral ? (pct >= 90 || pct <= 10 ? TH->purple : TH->text)
                           : risk >= 90 ? TH->bad : risk >= 75 ? TH->warn : risk <= 10 ? TH->ok : TH->text;
        txtIn(F_SEMI, TH->fsH3, (Rectangle){ c.x + 18, c.y + 14, cw - 150, 22 }, TH->text, 0, rc(g_psRows, i, 1));
        txtIn(F_SEMI, 28, (Rectangle){ c.x + cw - 140, c.y + 10, 122, 34 }, pc, 2, pctLabel(pct));
        txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + cw - 160, c.y + 44, 142, 14 }, TH->faint, 2, sstr("percentile %s", *dbMeta("prs_group_short") ? dbMeta("prs_group_short") : "in Europeans"));
        /* the distribution */
        float bins[40] = { 0 };
        double lo = 0, hi = 1;
        for (int k = 0; g_psDist && k < g_psDist->nrow; k++)
            if (!strcmp(rc(g_psDist, k, 0), rc(g_psRows, i, 0))) { parseFloats(rc(g_psDist, k, 3), bins, 40); lo = rf(g_psDist, k, 1); hi = rf(g_psDist, k, 2); }
        double you = rf(g_psRows, i, 14);
        histogram((Rectangle){ c.x + 18, c.y + 62, cw - 36, 70 }, bins, 40, (float)lo, (float)hi, (float)you, neutral ? TH->accent2 : risk < 75 && risk > 10 ? TH->accent : pc, NULL);
        txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 18, c.y + 136, cw - 36, 16 }, TH->dim, 0,
              sstr("z %+.2f   %s of %s variants used, %.0f%% called   %s", z, fmtBig((double)ri(g_psRows, i, 5)), fmtBig((double)ri(g_psRows, i, 4)),
                   100.0 * (double)ri(g_psRows, i, 6) / fmax(1, (double)ri(g_psRows, i, 5)), rc(g_psRows, i, 0)));
        txtIn(F_REG, TH->fsTiny, (Rectangle){ c.x + 18, c.y + 154, cw - 36, 16 }, TH->faint, 0, *rc(g_psRows, i, 11) ? rc(g_psRows, i, 11) : rc(g_psRows, i, 10));
        if (*note) txtPara(F_REG, TH->fsTiny + 1, c.x + 18, c.y + 176, cw - 36, 17, TH->dim, note, 6);
        if (gxBoxClick(gxIdN(33500 + i), c)) { lkBegin(rc(g_psRows, i, 1)); lkPgs(rc(g_psRows, i, 0), rc(g_psRows, i, 13)); lkOpen(); }
        if (++col == 2) { col = 0; y += rowH + 16; }
    }
    if (col) y += rowH + 16;
    y += 10;
    y += txtPara(F_REG, TH->fsTiny + 1, x, y, W, 17, TH->faint, sstr("Scores from the PGS Catalog; the metric under each is how well it predicted in its evaluation. "
                 "Percentiles are within a small reference (%s), so the top and bottom few percent are approximate.", dbMeta("prs_reference")), 4) + 30;
    gxScrollEnd(&g_psScroll, r, y - top);
}
