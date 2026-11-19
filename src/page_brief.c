/* ========================================================================
   src/page_brief.c - your genome in brief: the findings worth knowing first,
   gathered from every other page (pipeline step 22, table brief), each in
   plain words with what it means; selecting one opens its details.
   ======================================================================== */

static Rows *g_bfRows;
static GxScroll g_bfScroll;
static bool g_bfAsked;

static Color bfTone(const char *t) {
    return !strcmp(t, "bad") ? TH->bad : !strcmp(t, "warn") ? TH->warn : !strcmp(t, "accent") ? TH->accent : !strcmp(t, "ok") ? TH->ok : TH->info;
}

/* one finding: its height when !draw, and drawn (clickable) when draw */
static float bfItem(Rectangle c, Rows *r, int i, bool draw) {
    const char *title = rc(r, i, 3), *text = rc(r, i, 4), *mean = rc(r, i, 5), *target = rc(r, i, 6);
    float tw = c.width - 56;
    float h = 18 + txtParaH(F_SEMI, TH->fs + 1, tw, 22, title, 3) + 6;
    if (*text) h += txtParaH(F_REG, TH->fsSmall, tw, 19, text, 8) + 6;
    if (*mean) h += txtParaH(F_REG, TH->fsSmall, tw - 120, 19, mean, 8) + 6;
    h += 14;
    if (!draw) return h;
    c.height = h;
    bool hov = false;
    if (*target && gxBehave(gxIdN(71000 + i), c, &hov, NULL)) stEnter(target);
    gxCard(c);
    if (hov) { drRound(c, 10, drFade(TH->hover, 0.6f)); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
    drRound((Rectangle){ c.x, c.y + 14, 4, h - 28 }, 2, bfTone(rc(r, i, 2)));
    float x = c.x + 24, y = c.y + 18;
    y += txtPara(F_SEMI, TH->fs + 1, x, y, tw, 22, TH->text, title, 3) + 6;
    if (*text) y += txtPara(F_REG, TH->fsSmall, x, y, tw, 19, TH->text, text, 8) + 6;
    if (*mean) {
        txt(F_SEMI, TH->fsSmall, x, y, TH->accent2, "What it means");
        txtPara(F_REG, TH->fsSmall, x + 120, y, tw - 120, 19, TH->dim, mean, 8);
    }
    if (*target && hov) txt(F_REG, TH->fsTiny, c.x + c.width - 90, c.y + 12, TH->faint, "details  >");
    return h;
}

static void pageBrief(Rectangle r) {
    if (!dbHasTable("brief")) {
        pageMissing(r, "In brief", "The summary step did not run for this person. Re-import them from the People page.");
        return;
    }
    if (!g_bfAsked) {
        g_bfAsked = true;
        g_bfRows = dbq("SELECT section, ord, tone, title, text, meaning, target FROM brief ORDER BY ord", NULL);
    }
    float top = gxScrollBegin(&g_bfScroll, r);
    float x = r.x + 32, W = fminf(1100, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "Your genome in brief",
                         "The findings most worth knowing, gathered from every page and put in plain words: what was found, and what it "
                         "means. Select any one for the details. Nothing here is a diagnosis; anything that matters for health deserves a "
                         "conversation with a doctor.");
    const char *sec = "";
    for (int i = 0; g_bfRows && i < g_bfRows->nrow; i++) {
        if (strcmp(sec, rc(g_bfRows, i, 0))) {
            sec = rc(g_bfRows, i, 0);
            y += i ? 14 : 0;
            char up[64];
            snprintf(up, sizeof up, "%s", sec);
            for (char *p = up; *p; p++) *p = (char)toupper((unsigned char)*p);
            gxSectionTitle(x, y, up);
            y += 26;
        }
        Rectangle c = { x, y, W, 0 };
        float h = bfItem(c, g_bfRows, i, false);
        bfItem(c, g_bfRows, i, true);
        y += h + 10;
    }
    y += 30;
    gxScrollEnd(&g_bfScroll, r, y - top);
}
