/* ========================================================================
   src/page_diet.c - diet and supplements: what genes can and cannot say
   (pipeline step 23, table diet). Each topic shows the person's result, how
   strong the evidence is, what it means, and its sources.
   ======================================================================== */

static Rows *g_dtRows;
static GxScroll g_dtScroll;
static bool g_dtAsked;

static Color dtGrade(const char *g) {
    return !strcmp(g, "Guideline") ? TH->ok : !strcmp(g, "Strong") ? TH->accent : !strcmp(g, "Moderate") ? TH->info
         : !strcmp(g, "Weak") ? TH->dim : !strcmp(g, "Not supported") ? TH->warn : TH->faint;
}

static float dtCard(Rectangle c, Rows *r, int i, bool draw) {
    const char *topic = rc(r, i, 1), *grade = rc(r, i, 2), *res = rc(r, i, 3), *gt = rc(r, i, 4);
    const char *says = rc(r, i, 5), *means = rc(r, i, 6), *src = rc(r, i, 7), *target = rc(r, i, 8);
    bool intro = ri(r, i, 0) < 0;
    float tw = c.width - 48;
    float h = 18 + 26 + (*res ? 26 : 0) + (*gt ? 18 : 0);
    h += txtParaH(F_REG, TH->fsSmall, tw, 19, says, 8) + 8;
    h += txtParaH(F_REG, TH->fsSmall, tw - 120, 19, means, 8) + 8;
    if (*src) h += txtParaH(F_REG, TH->fsTiny + 1, tw, 16, src, 3) + 6;
    h += 12;
    if (!draw) return h;
    c.height = h;
    bool hov = false;
    if (*target && gxBehave(gxIdN(72000 + i), c, &hov, NULL)) stEnter(target);
    gxCard(c);
    if (hov) { drRound(c, 10, drFade(TH->hover, 0.6f)); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
    float x = c.x + 24, y = c.y + 18;
    txt(F_SEMI, TH->fsH3, x, y, TH->text, topic);
    if (*grade) {
        float pw = gxPillW(grade, TH->fsTiny);
        gxPill(c.x + c.width - 24 - pw, y + 2, grade, dtGrade(grade), TH->fsTiny);
    }
    y += 26;
    if (*res) { txt(F_SEMI, TH->fs, x, y, intro ? TH->text : TH->accent, sstr("You: %s", res)); y += 26; }
    if (*gt) { txtIn(F_MONO, 11, (Rectangle){ x, y, tw, 16 }, TH->faint, 0, gt); y += 18; }
    y += txtPara(F_REG, TH->fsSmall, x, y, tw, 19, TH->text, says, 8) + 8;
    txt(F_SEMI, TH->fsSmall, x, y, TH->accent2, intro ? "The evidence" : "What it means");
    y += txtPara(F_REG, TH->fsSmall, x + 120, y, tw - 120, 19, TH->dim, means, 8) + 8;
    if (*src) txtPara(F_REG, TH->fsTiny + 1, x, y, tw, 16, TH->faint, sstr("Sources: %s", src), 3);
    return h;
}

static void pageDiet(Rectangle r) {
    if (!dbHasTable("diet")) {
        pageMissing(r, "Diet and supplements", "The diet step did not run for this person. Re-import them from the People page.");
        return;
    }
    if (!g_dtAsked) {
        g_dtAsked = true;
        g_dtRows = dbq("SELECT ord, topic, grade, result, genotype, says, means, source, target FROM diet ORDER BY ord", NULL);
    }
    float top = gxScrollBegin(&g_dtScroll, r);
    float x = r.x + 32, W = fminf(1100, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Diet and supplements",
                         "What your genes can and cannot say about what you eat and take, each with how strong the evidence is: Guideline "
                         "(doctors act on it), Strong, Moderate, Weak, or Not supported (sold widely, but the evidence says otherwise). "
                         "Not medical advice: anything here that matters for you is worth a word with a doctor or dietitian.");
    for (int i = 0; g_dtRows && i < g_dtRows->nrow; i++) {
        Rectangle c = { x, y, W, 0 };
        float h = dtCard(c, g_dtRows, i, false);
        dtCard(c, g_dtRows, i, true);
        y += h + 12;
    }
    y += 30;
    gxScrollEnd(&g_dtScroll, r, y - top);
}
