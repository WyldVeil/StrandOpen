/* ========================================================================
   src/app.c - the application shell: theme, sidebar, top bar with search
   and history, toasts, page dispatch and startup.
   ======================================================================== */

static void themeInit(void) {
    g_t = (StTheme){
        .bg = { 22, 25, 31, 255 }, .side = { 18, 20, 26, 255 }, .panel = { 30, 34, 42, 255 },
        .panel2 = { 39, 44, 54, 255 }, .card = { 28, 32, 40, 255 }, .line = { 46, 52, 64, 255 },
        .lineHi = { 66, 74, 90, 255 }, .text = { 233, 236, 242, 255 }, .dim = { 160, 168, 184, 255 },
        .faint = { 108, 117, 135, 255 }, .accent = { 110, 178, 206, 255 }, .accent2 = { 129, 161, 193, 255 },
        .accentDim = { 52, 92, 112, 255 }, .ok = { 143, 192, 120, 255 }, .warn = { 232, 196, 110, 255 },
        .bad = { 220, 104, 112, 255 }, .info = { 120, 160, 230, 255 }, .purple = { 180, 142, 200, 255 },
        .hover = { 42, 48, 60, 255 }, .press = { 52, 58, 72, 255 }, .field = { 20, 23, 29, 255 },
        .sel = { 48, 92, 120, 255 }, .selDim = { 36, 58, 74, 255 }, .shadow = { 0, 0, 0, 90 },
        .onAccent = { 14, 22, 28, 255 }, .zebra = { 26, 30, 37, 255 },
        .fs = 14, .fsSmall = 13, .fsTiny = 11, .fsH1 = 26, .fsH2 = 20, .fsH3 = 16, .rowH = 30, .pad = 16,
    };
}

/* ---- navigation history ------------------------------------------------ */
#define NAV_MAX 64
static StLoc g_navHist[NAV_MAX];
static int g_navN = 1, g_navAt = 0;
static StLoc *navCur(void) { return &g_navHist[g_navAt]; }
static void navGo(StLoc l) {
    StLoc *c = navCur();
    if (!memcmp(c, &l, sizeof l)) return;
    if (g_navAt + 1 >= NAV_MAX) { memmove(g_navHist, g_navHist + 1, sizeof(StLoc) * (NAV_MAX - 1)); g_navAt--; }
    g_navHist[++g_navAt] = l;
    g_navN = g_navAt + 1;
    pageEnter(l.page);
}
static void navPage(StPage p) {
    StLoc l = *navCur();
    l.page = p;
    navGo(l);
}
static void navGene(int geneId) {
    StLoc l = *navCur();
    l.page = PG_GENE;
    l.geneId = geneId;
    navGo(l);
}
static void navGeneSym(const char *sym) {
    Rows *r = dbq("SELECT id FROM gene WHERE symbol=? COLLATE NOCASE ORDER BY biotype!='protein_coding' LIMIT 1", sym, NULL);
    if (r->nrow) navGene((int)ri(r, 0, 0));
    else toast("No gene named %s", sym);
    rowsFree(&r);
}
static void navLocus(int chrom, long long start, long long end) {
    StLoc l = *navCur();
    l.page = PG_BROWSER;
    l.chrom = chrom;
    l.start = start;
    l.end = end;
    navGo(l);
}
static void navVariant(long long varId) {
    Rows *r = dbq("SELECT chrom, pos FROM variant WHERE id=?", sstr("%lld", varId), NULL);
    if (r->nrow) {
        StLoc l = *navCur();
        l.page = PG_BROWSER;
        l.varId = varId;
        l.chrom = (int)ri(r, 0, 0);
        long long p = ri(r, 0, 1);
        l.start = p - 60;
        l.end = p + 60;
        navGo(l);
    }
    rowsFree(&r);
}
static void navBack(void) { if (g_navAt > 0) { g_navAt--; pageEnter(navCur()->page); } }
static void navFwd(void) { if (g_navAt + 1 < g_navN) { g_navAt++; pageEnter(navCur()->page); } }

/* ---- toasts -------------------------------------------------------------- */
static char g_toast[256];
static double g_toastAt = -100;
static void toast(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(g_toast, sizeof g_toast, fmt, ap);
    va_end(ap);
    g_toastAt = (GetTime)();
}
static void toastDraw(void) {
    double age = (GetTime)() - g_toastAt;
    if (age > 3.2 || !g_toast[0]) return;
    float a = age < 0.15 ? (float)(age / 0.15) : age > 2.7 ? (float)((3.2 - age) / 0.5) : 1.0f;
    float w = (float)txtW(F_REG, TH->fs, g_toast) + 40;
    Rectangle r = { ((float)gxScreenW() - w) / 2, (float)gxScreenH() - 70, w, 40 };
    drRound(r, 10, drFade(TH->panel2, a));
    drRoundLine(r, 10, drFade(TH->lineHi, a));
    txtIn(F_REG, TH->fs, r, drFade(TH->text, a), 1, g_toast);
}
/* a web page for a database record, opened only when the person clicks */
static void openUrl(const char *url) {
    if (getenv("STRAND_NO_OPEN")) return;
    char *argv[] = { "xdg-open", (char *)url, NULL };
    Proc p = spawnReader(argv);
    if (p.pid) reap(&p, false);
    toast("Opened in your browser");
}
static void copyText(const char *s) {
    SetClipboardText(s);
    toast("Copied: %.60s", s);
}

/* ---- sidebar ------------------------------------------------------------- */
typedef struct { StPage page; const char *label; GxIconFn icon; const char *group; } NavItem;
static const NavItem g_navItems[] = {
    { PG_PEOPLE, "People", icPerson, "PEOPLE" },
    { PG_BRIEF, "In brief", icStar, "EXPLORE" },
    { PG_OVERVIEW, "Overview", icHome, NULL },
    { PG_BROWSER, "Genome browser", icHelix, NULL },
    { PG_VARIANTS, "Variants", icList, NULL },
    { PG_GENE, "Genes", icDots, NULL },
    { PG_RARE, "Rare and remarkable", icGem, NULL },
    { PG_HEALTH, "Clinical findings", icHeart, "HEALTH" },
    { PG_CARRIER, "Carrier status", icShield, NULL },
    { PG_PGX, "Medications", icPill, NULL },
    { PG_DIET, "Diet and supplements", icPill, NULL },
    { PG_PRS, "Polygenic scores", icChart, NULL },
    { PG_RARITY, "Your rarest variants", icRarity, "YOU" },
    { PG_TRAITS, "Traits", icStar, NULL },
    { PG_ANCESTRY, "Ancestry", icGlobe, NULL },
    { PG_LINEAGE, "Lineage", icTree, NULL },
    { PG_RAW, "Raw data & QC", icFile, "DATA" },
    { PG_ABOUT, "Sources & notes", icInfo, NULL },
};
static void sidebar(Rectangle r) {
    drRect(r, TH->side);
    drVLine(r.x + r.width - 1, r.y, r.height, TH->line);
    /* the mark */
    icHelix((Rectangle){ r.x + 22, r.y + 20, 18, 26 }, TH->accent);
    txt(F_SEMI, 20, r.x + 50, r.y + 21, TH->text, "Strand Open");
    txt(F_REG, TH->fsTiny, r.x + 51, r.y + 45, TH->faint, "genome explorer");
    /* every entry must stay reachable: on a short window (or a large interface
       size) the rows tighten a little, and if that is still not enough the
       list scrolls */
    int ngroup = 0;
    for (size_t i = 0; i < NELEMS(g_navItems); i++) ngroup += g_navItems[i].group != NULL;
    float need = (float)ngroup * 32 + (float)NELEMS(g_navItems) * 36, avail = r.height - 76 - 8;
    float k = need > avail ? fmaxf(0.75f, avail / need) : 1.0f;
    static GxScroll scroll;
    Rectangle list = { r.x, r.y + 76, r.width, r.height - 76 };
    float top = gxScrollBegin(&scroll, list), y = top;
    StPage cur = navCur()->page;
    for (size_t i = 0; i < NELEMS(g_navItems); i++) {
        const NavItem *it = &g_navItems[i];
        if (it->group) {
            y += 10 * k;
            gxSectionTitle(r.x + 22, y, it->group);
            y += 22 * k;
        }
        Rectangle row = { r.x + 10, y, r.width - 20, 34 * k };
        bool h;
        if (gxBehave(gxIdN(1000 + i), row, &h, NULL)) navPage(it->page);
        bool on = cur == it->page;
        if (on) drRound(row, 7, TH->selDim);
        else if (h) drRound(row, 7, TH->hover);
        if (on) drRound((Rectangle){ row.x, row.y + 8, 3, row.height - 16 }, 1.5f, TH->accent);
        if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
        it->icon((Rectangle){ row.x + 14, row.y + (row.height - 16) / 2, 16, 16 }, on ? TH->accent : TH->dim);
        txtIn(on ? F_SEMI : F_REG, TH->fs, (Rectangle){ row.x + 42, row.y, row.width - 50, row.height }, on ? TH->text : TH->dim, 0, it->label);
        y += 36 * k;
    }
    float listH = y - top + 2;
    gxScrollEnd(&scroll, list, listH);
    /* the sample, at the foot */
    const char *sample = dbMeta("sample_id");
    float fy = r.y + r.height - 58;
    if (fy > list.y + listH + 10) {
        drHLine(r.x + 16, fy - 10, r.width - 32, TH->line);
        txt(F_SEMI, TH->fsSmall, r.x + 22, fy, TH->dim, *sample ? sample : "Sample");
        txt(F_REG, TH->fsTiny, r.x + 22, fy + 20, TH->faint, sstr("%s, GRCh38", dkName()));
    }
}

/* ---- global search ---------------------------------------------------------- */
typedef struct { int kind; long long id; char title[160], sub[200], key[64]; int chrom; long long a, b; } SearchHit;
enum { SH_GENE, SH_VARIANT, SH_LOCUS, SH_CONDITION, SH_PAGE, SH_DRUG, SH_TRAIT, SH_REMARK };
static char g_search[200];
static SearchHit g_hits[24];
static int g_hitN, g_hitSel;
static bool g_searchOpen;
static char g_searchFor[200];
static void searchRun(void) {
    g_hitN = 0;
    g_hitSel = 0;
    snprintf(g_searchFor, sizeof g_searchFor, "%s", g_search);
    const char *q = g_search;
    while (*q == ' ') q++;
    if (!*q) return;
    int chrom;
    long long a, b;
    if (parseLocus(q, &chrom, &a, &b)) {
        SearchHit *h = &g_hits[g_hitN++];
        h->kind = SH_LOCUS; h->chrom = chrom; h->a = a; h->b = b;
        snprintf(h->title, sizeof h->title, "%s", a == b ? fmtPos(chrom, a) : sstr("chr%s:%s-%s", chromName(chrom), fmtInt(a), fmtInt(b)));
        snprintf(h->sub, sizeof h->sub, "Open this position in the genome browser");
    }
    if ((q[0] == 'r' || q[0] == 'R') && (q[1] == 's' || q[1] == 'S') && isdigit((unsigned char)q[2])) {
        Rows *r = dbq("SELECT id, chrom, pos, ref, alt, gene, zyg FROM variant WHERE rsid=? LIMIT 4", q + 2, NULL);
        for (int i = 0; i < r->nrow && g_hitN < 24; i++) {
            SearchHit *h = &g_hits[g_hitN++];
            h->kind = SH_VARIANT; h->id = ri(r, i, 0);
            snprintf(h->title, sizeof h->title, "rs%s  %s>%s", q + 2, rc(r, i, 3), rc(r, i, 4));
            snprintf(h->sub, sizeof h->sub, "%s %s  %s  (you: %s)", fmtPos((int)ri(r, i, 1), ri(r, i, 2)), rc(r, i, 5), "", zygName((int)ri(r, i, 6)));
        }
        if (!r->nrow) {
            SearchHit *h = &g_hits[g_hitN++];
            h->kind = SH_PAGE; h->id = -1;
            snprintf(h->title, sizeof h->title, "rs%s: not among your variants", q + 2);
            snprintf(h->sub, sizeof h->sub, "You most likely carry the reference allele here (see Traits for curated SNPs)");
        }
        rowsFree(&r);
    }
    char like[220];
    snprintf(like, sizeof like, "%s%%", q);
    Rows *r = dbq("SELECT id, symbol, name, chrom, start, nvar FROM gene WHERE symbol LIKE ? "
                  "ORDER BY length(symbol), biotype!='protein_coding' LIMIT 8", like, NULL);
    for (int i = 0; i < r->nrow && g_hitN < 24; i++) {
        SearchHit *h = &g_hits[g_hitN++];
        h->kind = SH_GENE; h->id = ri(r, i, 0);
        snprintf(h->title, sizeof h->title, "%s", rc(r, i, 1));
        snprintf(h->sub, sizeof h->sub, "%s  -  %s variants", *rc(r, i, 2) ? rc(r, i, 2) : "gene", fmtInt(ri(r, i, 5)));
    }
    rowsFree(&r);
    char mid[220];
    snprintf(mid, sizeof mid, "%%%s%%", q);
    if (strlen(q) >= 3 && g_hitN < 20) {
        char fts[240];
        snprintf(fts, sizeof fts, "\"%s\"*", q);
        for (char *c = fts + 1; *c; c++) if (*c == '"' && c[1] != '*') *c = ' ';
        r = dbq("SELECT g.id, g.symbol, g.name FROM gene_fts f JOIN gene g ON g.id=f.rowid WHERE gene_fts MATCH ? "
                "AND g.symbol NOT LIKE ? ORDER BY rank LIMIT 6", fts, like, NULL);
        for (int i = 0; i < r->nrow && g_hitN < 24; i++) {
            SearchHit *h = &g_hits[g_hitN++];
            h->kind = SH_GENE; h->id = ri(r, i, 0);
            snprintf(h->title, sizeof h->title, "%s", rc(r, i, 1));
            snprintf(h->sub, sizeof h->sub, "%s", rc(r, i, 2));
        }
        rowsFree(&r);
        if (dbHasTable("pgx_drug")) {
            r = dbq("SELECT drug, group_concat(DISTINCT gene), max(recommendation) FROM pgx_drug WHERE drug LIKE ? AND recommendation!='' GROUP BY drug LIMIT 4", like, NULL);
            for (int i = 0; i < r->nrow && g_hitN < 24; i++) {
                SearchHit *h = &g_hits[g_hitN++];
                h->kind = SH_DRUG; h->id = PG_PGX;
                snprintf(h->title, sizeof h->title, "%s", rc(r, i, 0));
                snprintf(h->sub, sizeof h->sub, "%s: %.150s", rc(r, i, 1), rc(r, i, 2));
            }
            rowsFree(&r);
        }
        if (dbHasTable("trait")) {
            r = dbq("SELECT title, result, cat FROM trait WHERE title LIKE ? OR id LIKE ? LIMIT 4", mid, like, NULL);
            for (int i = 0; i < r->nrow && g_hitN < 24; i++) {
                SearchHit *h = &g_hits[g_hitN++];
                h->kind = SH_TRAIT; h->id = !strcmp(rc(r, i, 2), "Health") ? PG_HEALTH : PG_TRAITS;
                snprintf(h->title, sizeof h->title, "%s", rc(r, i, 0));
                snprintf(h->sub, sizeof h->sub, "%s", rc(r, i, 1));
            }
            rowsFree(&r);
        }
        if (dbHasTable("remarkable")) {
            r = dbq("SELECT key, title, CASE status WHEN 2 THEN 'You carry it' WHEN 1 THEN 'You carry one copy' WHEN 0 THEN 'Not carried' ELSE 'Not called' END, rarity "
                    "FROM remarkable WHERE title LIKE ? OR gene LIKE ? OR effect LIKE ? OR cat LIKE ? ORDER BY status DESC LIMIT 4", mid, like, mid, mid, NULL);
            for (int i = 0; i < r->nrow && g_hitN < 24; i++) {
                SearchHit *h = &g_hits[g_hitN++];
                h->kind = SH_REMARK; h->id = 0;
                snprintf(h->title, sizeof h->title, "%s", rc(r, i, 1));
                snprintf(h->sub, sizeof h->sub, "%s; %s", rc(r, i, 2), rc(r, i, 3));
                snprintf(h->key, sizeof h->key, "%s", rc(r, i, 0));
            }
            rowsFree(&r);
        }
        r = dbq("SELECT id, gene, clin_dn, clin, zyg FROM variant WHERE clin>0 AND clin NOT IN (1,2) AND clin_dn LIKE ? "
                "ORDER BY clin IN (5,6) DESC, clin_stars DESC LIMIT 6", mid, NULL);
        for (int i = 0; i < r->nrow && g_hitN < 24; i++) {
            SearchHit *h = &g_hits[g_hitN++];
            h->kind = SH_CONDITION; h->id = ri(r, i, 0);
            snprintf(h->title, sizeof h->title, "%s (%s)", rc(r, i, 1), clinName((int)ri(r, i, 3)));
            snprintf(h->sub, sizeof h->sub, "%.150s", rc(r, i, 2));
        }
        rowsFree(&r);
    }
}
static void searchPick(int i) {
    if (i < 0 || i >= g_hitN) return;
    SearchHit *h = &g_hits[i];
    switch (h->kind) {
    case SH_GENE: navGene((int)h->id); break;
    case SH_VARIANT: case SH_CONDITION: navVariant(h->id); break;
    case SH_DRUG: snprintf(g_pxQ, sizeof g_pxQ, "%s", h->title); g_pxTab = 1; g_pxScroll = (GxScroll){ 0 }; g_pxSel = -1; navPage(PG_PGX); break;
    case SH_TRAIT: navPage((StPage)h->id); break;
    case SH_REMARK: rareFocus(h->key); break;
    case SH_LOCUS: {
        long long a = h->a, b = h->b;
        if (b - a < 100) { a -= 60; b += 60; }
        navLocus(h->chrom, a, b);
        break;
    }
    default: return;
    }
    g_searchOpen = false;
    g_search[0] = 0;
    g_gx.editId = 0;
}
static GxId g_searchId;
static void searchPopup(Rectangle r, void *u) {
    (void)u;
    float y = r.y + 6;
    if (!g_hitN) {
        txtIn(F_REG, TH->fsSmall, (Rectangle){ r.x + 14, y, r.width - 28, 30 }, TH->faint, 0, "No matches. Try a gene (BRCA2), an rsID (rs429358) or a position (chr19:44908684).");
        return;
    }
    static const char *const kinds[] = { "GENE", "VARIANT", "LOCUS", "CONDITION", "", "MEDICINE", "TRAIT", "REMARKABLE" };
    for (int i = 0; i < g_hitN; i++) {
        Rectangle row = { r.x + 6, y, r.width - 12, 46 };
        bool h;
        if (gxBehave(gxIdN(5000 + i), row, &h, NULL)) { searchPick(i); gxPopupClose(); return; }
        if (h) g_hitSel = i;
        if (g_hitSel == i) drRound(row, 6, TH->hover);
        Color kc = g_hits[i].kind == SH_GENE ? TH->accent : g_hits[i].kind == SH_CONDITION ? TH->warn : g_hits[i].kind == SH_DRUG ? TH->ok : g_hits[i].kind == SH_REMARK ? TH->accent2 : TH->purple;
        txt(F_SEMI, TH->fsTiny - 1, row.x + 12, row.y + 8, kc, kinds[g_hits[i].kind]);
        txtIn(F_SEMI, TH->fs, (Rectangle){ row.x + 90, row.y + 4, row.width - 100, 20 }, TH->text, 0, g_hits[i].title);
        txtIn(F_REG, TH->fsSmall - 1, (Rectangle){ row.x + 90, row.y + 24, row.width - 100, 18 }, TH->dim, 0, g_hits[i].sub);
        y += 48;
    }
}
static void topBar(Rectangle r) {
    drRect(r, TH->bg);
    drHLine(r.x, r.y + r.height - 1, r.width, TH->line);
    float x = r.x + 14;
    if (gxIconButton(gxId("nav back"), (Rectangle){ x, r.y + 12, 30, 30 }, icBack, false, "Back (Alt+Left)") ) navBack();
    x += 34;
    if (gxIconButton(gxId("nav fwd"), (Rectangle){ x, r.y + 12, 30, 30 }, icFwd, false, "Forward (Alt+Right)")) navFwd();
    x += 46;
    float sw = fminf(560, r.width - 360);
    Rectangle sr = { x, r.y + 10, sw, 34 };
    g_searchId = gxId("global search");
    int res = gxTextField(g_searchId, sr, g_search, sizeof g_search, "Search genes, rsIDs, positions, conditions   (Ctrl+K)", true);
    bool focused = g_gx.editId == g_searchId;
    if (res == 2 || (focused && strcmp(g_searchFor, g_search))) searchRun();
    if (focused && g_search[0]) {
        if (g_gx.popup != g_searchId)
            gxPopupOpen(g_searchId, sr, (Rectangle){ sr.x, sr.y + sr.height + 6, sr.width, 12 + 48 * (float)(g_hitN ? g_hitN : 1) }, searchPopup, NULL);
        g_gx.popupRect.height = 12 + 48 * (float)(g_hitN ? g_hitN : 1);
        if (IsKeyPressed(KEY_DOWN)) g_hitSel = (g_hitSel + 1) % (g_hitN ? g_hitN : 1);
        if (IsKeyPressed(KEY_UP)) g_hitSel = (g_hitSel + g_hitN - 1) % (g_hitN ? g_hitN : 1);
    } else if (g_gx.popup == g_searchId && !focused) {
        /* keep it while the mouse is on it (a click on a hit loses focus first) */
        if (!CheckCollisionPointRec(g_gx.mouse, g_gx.popupRect)) gxPopupClose();
    }
    if (res == 1 && g_hitN) { searchPick(g_hitSel); gxPopupClose(); }
    /* where we are */
    const char *crumb = "";
    for (size_t i = 0; i < NELEMS(g_navItems); i++) if (g_navItems[i].page == navCur()->page) crumb = g_navItems[i].label;
    txtIn(F_REG, TH->fsSmall, (Rectangle){ r.x + r.width - 300, r.y, 280, r.height }, TH->faint, 2, crumb);
}

/* ---- pages ----------------------------------------------------------------- */
static void pageEnter(StPage p) { (void)p; }

static void pageDraw(Rectangle r) {
    switch (navCur()->page) {
    case PG_OVERVIEW: pageOverview(r); break;
    case PG_BROWSER: pageBrowser(r); break;
    case PG_VARIANTS: pageVariants(r); break;
    case PG_GENE: pageGene(r); break;
    case PG_HEALTH: pageHealth(r); break;
    case PG_CARRIER: pageCarrier(r); break;
    case PG_PGX: pagePgx(r); break;
    case PG_TRAITS: pageTraits(r); break;
    case PG_PRS: pagePrs(r); break;
    case PG_ANCESTRY: pageAncestry(r); break;
    case PG_LINEAGE: pageLineage(r); break;
    case PG_RAW: pageRaw(r); break;
    case PG_ABOUT: pageAbout(r); break;
    case PG_RARE: pageRare(r); break;
    case PG_RARITY: pageRarity(r); break;
    case PG_PEOPLE: pagePeople(r); break;
    case PG_BRIEF: pageBrief(r); break;
    case PG_DIET: pageDiet(r); break;
    default: break;
    }
}

/* a page that has nothing yet (its pipeline step has not run) */
static void pageMissing(Rectangle r, const char *title, const char *why) {
    Rectangle c = { r.x + 32, r.y + 32, fminf(720, r.width - 64), 150 };
    gxCard(c);
    txt(F_SEMI, TH->fsH2, c.x + 24, c.y + 22, TH->text, title);
    txtPara(F_REG, TH->fs, c.x + 24, c.y + 60, c.width - 48, 21, TH->dim, why, 4);
}
/* page heading: title, subtitle; returns the y below it */
static float pageHeader(Rectangle r, const char *title, const char *sub) {
    txt(F_SEMI, TH->fsH1, r.x + 32, r.y + 26, TH->text, title);
    float y = r.y + 26 + (float)TH->fsH1 + 10;
    if (sub) { y += txtPara(F_REG, TH->fs, r.x + 32, y, fminf(980, r.width - 64), 21, TH->dim, sub, 5); }
    return y + 18;
}

/* ---- the frame --------------------------------------------------------------- */
static bool g_stReady;
static float g_stScale = 1.0f;
static void stShortcuts(void) {
    bool ctrl = IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL);
    bool alt = IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT);
    if (ctrl && IsKeyPressed(KEY_K)) { gxFocus(g_searchId, "", false); }
    if (alt && IsKeyPressed(KEY_LEFT)) navBack();
    if (alt && IsKeyPressed(KEY_RIGHT)) navFwd();
    if (IsMouseButtonPressed(3)) navBack();
    if (IsMouseButtonPressed(4)) navFwd();
    float ns = g_stScale;
    if (ctrl && (IsKeyPressed(KEY_EQUAL) || IsKeyPressed(KEY_KP_ADD))) ns = fminf(2.0f, g_stScale + 0.1f);
    if (ctrl && (IsKeyPressed(KEY_MINUS) || IsKeyPressed(KEY_KP_SUBTRACT))) ns = fmaxf(0.7f, g_stScale - 0.1f);
    if (ctrl && IsKeyPressed(KEY_ZERO)) ns = 1.0f;
    if (ns != g_stScale) {
        g_stScale = ns; g_gx.scale = ns; fontsReload(ns);
        toast("Interface size %d%%", (int)lroundf(ns * 100));
        cfgSet("so.scale", sstr("%.2f", ns));
        cfgSave();
    }
}
/* a data app left open should not draw sixty frames a second: after a few
   quiet seconds (no input, no query or read fetch in flight, no animation)
   it drops to 15 frames a second, and any input brings it straight back */
static double g_stActiveAt;
static bool g_stIdle;
static void stPace(void) {
    Vector2 md = GetMouseDelta();
    static const int keys[] = { KEY_LEFT, KEY_RIGHT, KEY_UP, KEY_DOWN, KEY_PAGE_UP, KEY_PAGE_DOWN, KEY_LEFT_CONTROL, KEY_RIGHT_CONTROL,
                                KEY_LEFT_ALT, KEY_ESCAPE, KEY_ENTER, KEY_TAB, KEY_EQUAL, KEY_MINUS };
    bool active = md.x != 0 || md.y != 0 || g_gx.wheel != 0 || g_gx.down[0] || g_gx.down[1] || g_gx.charN > 0 ||
                  g_gx.editId != 0 || readsState() == 1 || g_gx.popup != 0;
    for (size_t k = 0; k < NELEMS(keys) && !active; k++) if (IsKeyDown(keys[k])) active = true;
    for (int s = 0; s < SLOT__N && !active; s++) if (dbBusy((StSlot)s)) active = true;
    double t = (GetTime)();
    if (active) g_stActiveAt = t;
    bool idle = t - g_stActiveAt > 3.0;
    if (idle != g_stIdle) {
        g_stIdle = idle;
        (SetTargetFPS)(idle ? 15 : 0);
    }
}
static void stFrame(float dt, double now) {
    (void)dt; (void)now;
    ClearBackground(TH->bg);
    gxBeginFrame();
    if (!getenv("STRAND_NO_OPEN")) stPace();
    float W = (float)gxScreenW(), H = (float)gxScreenH();
    if (!g_stReady) {
        /* nobody imported yet (or the chosen person has no results): only the People page */
        drRect((Rectangle){ 0, 0, W, 54 }, TH->side);
        icHelix((Rectangle){ 22, 14, 18, 26 }, TH->accent);
        txt(F_SEMI, 20, 50, 16, TH->text, "Strand Open");
        Rectangle content = { 0, 54, W, H - 54 };
        gxClipPush(content);
        pagePeople(content);
        gxClipPop();
        toastDraw();
        gxEndFrame();
        return;
    }
    if (!g_gx.keysTaken) stShortcuts();
    else {
        bool ctrl = IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL);
        if (ctrl && IsKeyPressed(KEY_K)) gxFocus(g_searchId, "", false);
    }
    float sideW = 230, topH = 54;
    Rectangle content = { sideW, topH, W - sideW, H - topH };
    gxClipPush(content);
    pageDraw(content);
    gxClipPop();
    sidebar((Rectangle){ 0, 0, sideW, H });
    topBar((Rectangle){ sideW, 0, W - sideW, topH });
    toastDraw();
    gxEndFrame();
}

static void stInit(void) {
    themeInit();
    fontsInit();
    SetWindowMinSize(1100, 700);
    SetExitKey(KEY_NULL);
    /* which person: --profile, else the one open last time; nobody yet -> the People page */
    soResolveRoot();
    snprintf(g_profileSlug, sizeof g_profileSlug, "%s", *g_argProfile ? g_argProfile : cfgGet("so.profile", ""));
    snprintf(g_dataDir, sizeof g_dataDir, "%s/profiles/%s", g_root, g_profileSlug);
    char path[1200];
    snprintf(path, sizeof path, "%s/db/strand.db", g_dataDir);
    if (!*g_profileSlug || access(path, R_OK) != 0 || !dbOpen(path)) {
        g_profileSlug[0] = 0;
        g_navHist[0] = (StLoc){ .page = PG_PEOPLE };
        return;
    }
    genomeInit();
    g_stReady = true;
    g_navHist[0] = (StLoc){ .page = PG_OVERVIEW, .chrom = 7, .start = 117480000, .end = 117680000 };
    /* where the last session left off */
    float sc = (float)atof(cfgGet("so.scale", "1"));
    if (sc >= 0.7f && sc <= 2.0f && sc != 1.0f) { g_stScale = sc; g_gx.scale = sc; fontsReload(sc); }
    int pg = atoi(cfgGet("so.page", "0"));
    if (pg > 0 && pg < PG__N && pg != PG_GENE) g_navHist[0].page = (StPage)pg;
    int c; long long a, b;
    if (parseLocus(cfgGet("so.locus", ""), &c, &a, &b)) { g_navHist[0].chrom = c; g_navHist[0].start = a; g_navHist[0].end = b; }
}
static void stShutdown(void) {
    if (g_stReady) {
        StLoc *l = navCur();
        cfgSet("so.page", sstr("%d", (int)l->page));
        if (l->chrom) cfgSet("so.locus", sstr("%s:%lld-%lld", chromName(l->chrom), l->start, l->end));
        cfgSave();
    }
    dbShutdown();
}

static void stEnter(const char *name) {
    static const struct { const char *n; StPage p; } m[] = {
        { "overview", PG_OVERVIEW }, { "browser", PG_BROWSER }, { "variants", PG_VARIANTS }, { "gene", PG_GENE },
        { "health", PG_HEALTH }, { "carrier", PG_CARRIER }, { "pgx", PG_PGX }, { "traits", PG_TRAITS },
        { "prs", PG_PRS }, { "ancestry", PG_ANCESTRY }, { "lineage", PG_LINEAGE }, { "raw", PG_RAW }, { "about", PG_ABOUT }, { "rare", PG_RARE }, { "rarity", PG_RARITY },
        { "brief", PG_BRIEF }, { "people", PG_PEOPLE }, { "diet", PG_DIET },
    };
    if (!g_stReady) return;
    char page[64] = "", arg[128] = "";
    sscanf(name, "%63[^:]:%127s", page, arg);
    for (size_t i = 0; i < NELEMS(m); i++) if (!strcmp(page, m[i].n)) navPage(m[i].p);
    if (!strcmp(page, "gene") && arg[0]) navGeneSym(arg);
    if (!strcmp(page, "browser") && arg[0]) {
        int c; long long a, b;
        if (parseLocus(arg, &c, &a, &b)) navLocus(c, a, b);
        else {
            Rows *r = dbq("SELECT chrom, start, end FROM gene WHERE symbol=? LIMIT 1", arg, NULL);
            if (r->nrow) navLocus((int)ri(r, 0, 0), ri(r, 0, 1) - 2000, ri(r, 0, 2) + 2000);
            rowsFree(&r);
        }
    }
    if (!strcmp(page, "variant") && arg[0]) navVariant(atoll(arg));
    if (!strcmp(page, "export")) exportReport();
    /* a tab jump lands at the top of that tab, as a click on the tab does */
    if (!strcmp(page, "health") && arg[0]) { g_hTab = clampi(atoi(arg), 0, (int)NELEMS(g_hTabs) - 1); g_hScroll = (GxScroll){ 0 }; }
    if (!strcmp(page, "pgx") && arg[0]) { g_pxTab = clampi(atoi(arg), 0, (int)NELEMS(g_pxTabs) - 1); g_pxScroll = (GxScroll){ 0 }; g_pxSel = -1; }
    if (!strcmp(page, "rare") && arg[0]) { g_rTab = clampi(atoi(arg), 0, (int)NELEMS(g_rTabs) - 1); g_rScroll = (GxScroll){ 0 }; }
    if (!strcmp(page, "rarity") && arg[0]) { g_ryTab = clampi(atoi(arg), 0, (int)NELEMS(g_ryTabs) - 1); g_ryScroll = (GxScroll){ 0 }; }
    if (!strcmp(page, "lineage") && arg[0]) { if (!g_lnAsked) lnAsk(); g_lnLine = atoi(arg) ? 1 : 0; }
    if (!strcmp(page, "linkmenu") && arg[0]) {     /* a test hook: the link menu of a variant id or a gene */
        if (isdigit((unsigned char)arg[0])) lkVariantMenu(atoll(arg));
        else lkGeneMenu(arg);
    }
    if (!strcmp(page, "ancestry") && arg[0]) g_anPcaSet = atoi(arg);
    if (!strcmp(page, "search") && arg[0]) {
        snprintf(g_search, sizeof g_search, "%s", arg);
        for (char *c = g_search; *c; c++) if (*c == '_') *c = ' ';
        gxFocus(g_searchId ? g_searchId : (g_searchId = gxHash(0x9e3779b97f4a7c15ull, "global search", 13)), g_search, false);
        searchRun();
    }
}

static const unsigned char g_stIcon[] = {
#embed "../assets/icon/strand-256.png"
};
int main(int argc, char **argv) {
    AppDef def = {
        .title = "Strand Open",
        .iconPng = g_stIcon, .iconPngLen = sizeof g_stIcon,
        .w = 1600, .h = 1000,
        .init = stInit,
        .frame = stFrame,
        .shutdown = stShutdown,
        .enter = stEnter,
    };
    return appMain(&def, argc, argv);
}
