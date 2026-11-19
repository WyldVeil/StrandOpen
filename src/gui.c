/* ========================================================================
   src/gui.c - Strand Open's immediate-mode widget kit (gx*).

   Ids hashed from labels on an id stack, hot / active / focus, one pop-up
   at a time drawn over everything, a UI scale applied as one transform with
   the mouse mapped back through it; for a data application: a virtual
   table with sortable headers, cards, pills, meters, segmented controls and
   a scroll area with smooth wheel scrolling.
   ======================================================================== */

typedef uint64_t GxId;

typedef struct {
    Vector2 mouse;
    bool down[3], pressed[3], released[3], dbl;
    float wheel;
    double now, lastClick;
    Vector2 lastClickPos;
    GxId hot, active, focus;
    bool overUi;
    GxId idStack[32];
    int idN;
    Rectangle clip[32];
    int clipN;
    /* tooltip */
    GxId tipId;
    double tipSince;
    char tip[600];
    Vector2 tipAt;
    bool tipShown;
    /* text editing */
    GxId editId;
    char edit[512];
    int editCur, editSel, editLen;
    float editScroll;
    int chars[32], charN;
    bool keysTaken;
    /* pop-up */
    GxId popup;
    Rectangle popupAnchor, popupRect, popupRectPrev;
    bool popupOpening;
    void (*popupDraw)(Rectangle anchor, void *u);
    void *popupUser;
    int cursor;
    float scale;
    long frame;
} GxState;
static GxState g_gx = { .scale = 1.0f };

static int gxScreenW(void) { return (int)ceilf((float)(GetScreenWidth)() / g_gx.scale - 0.01f); }
static int gxScreenH(void) { return (int)ceilf((float)(GetScreenHeight)() / g_gx.scale - 0.01f); }

static GxId gxHash(GxId seed, const void *data, size_t n) {
    GxId h = seed ^ 1469598103934665603ull;
    const unsigned char *p = data;
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ull; }
    return h ? h : 1;
}
static GxId gxTop(void) { return g_gx.idN ? g_gx.idStack[g_gx.idN - 1] : 0x9e3779b97f4a7c15ull; }
static GxId gxId(const char *s) { return gxHash(gxTop(), s, strlen(s)); }
static GxId gxIdN(uint64_t v) { return gxHash(gxTop(), &v, sizeof v); }
static void gxPush(GxId id) { if (g_gx.idN < 32) g_gx.idStack[g_gx.idN++] = id; }
static void gxPushS(const char *s) { gxPush(gxId(s)); }
static void gxPop(void) { if (g_gx.idN) g_gx.idN--; }

/* ---- the frame ------------------------------------------------------ */
static void gxScaleApply(void) {
    rlDrawRenderBatchActive();
    rlMatrixMode(RL_MODELVIEW);
    rlLoadIdentity();
    if (g_gx.scale != 1.0f) rlScalef(g_gx.scale, g_gx.scale, 1.0f);
}
static void gxBeginFrame(void) {
    g_gx.frame++;
    g_gx.mouse = Vector2Scale(GetMousePosition(), 1.0f / g_gx.scale);
    for (int b = 0; b < 3; b++) {
        g_gx.down[b] = IsMouseButtonDown(b);
        g_gx.pressed[b] = IsMouseButtonPressed(b);
        g_gx.released[b] = IsMouseButtonReleased(b);
    }
    g_gx.wheel = GetMouseWheelMove();
    g_gx.now = (GetTime)();
    g_gx.dbl = false;
    if (g_gx.pressed[0]) {
        if (g_gx.now - g_gx.lastClick < 0.35 && Vector2Distance(g_gx.mouse, g_gx.lastClickPos) < 5) g_gx.dbl = true;
        g_gx.lastClick = g_gx.dbl ? 0 : g_gx.now;
        g_gx.lastClickPos = g_gx.mouse;
    }
    g_gx.hot = 0;
    g_gx.overUi = false;
    g_gx.clipN = 0;
    g_gx.idN = 0;
    g_gx.charN = 0;
    int c;
    while ((c = GetCharPressed()) > 0 && g_gx.charN < 32) g_gx.chars[g_gx.charN++] = c;
    g_gx.keysTaken = g_gx.editId != 0;
    g_gx.cursor = MOUSE_CURSOR_DEFAULT;
    g_gx.popupRectPrev = g_gx.popup ? g_gx.popupRect : (Rectangle){ 0 };
    g_gx.tipShown = false;
    gxScaleApply();
}

/* ---- clipping ------------------------------------------------------- */
static Rectangle gxRectX(Rectangle a, Rectangle b) {
    float x0 = fmaxf(a.x, b.x), y0 = fmaxf(a.y, b.y);
    float x1 = fminf(a.x + a.width, b.x + b.width), y1 = fminf(a.y + a.height, b.y + b.height);
    return (Rectangle){ x0, y0, fmaxf(0, x1 - x0), fmaxf(0, y1 - y0) };
}
static void gxClipApply(void) {
    rlDrawRenderBatchActive();
    if (!g_gx.clipN) { (EndScissorMode)(); return; }
    Rectangle r = g_gx.clip[g_gx.clipN - 1];
    float k = g_gx.scale;
    int x0 = (int)floorf(r.x * k), y0 = (int)floorf(r.y * k), x1 = (int)ceilf((r.x + r.width) * k), y1 = (int)ceilf((r.y + r.height) * k);
    (BeginScissorMode)(x0, y0, x1 - x0, y1 - y0);
}
static void gxClipPush(Rectangle r) {
    if (g_gx.clipN) r = gxRectX(r, g_gx.clip[g_gx.clipN - 1]);
    if (g_gx.clipN < 32) g_gx.clip[g_gx.clipN++] = r;
    gxClipApply();
}
static void gxClipPop(void) { if (g_gx.clipN) g_gx.clipN--; gxClipApply(); }
static bool gxInClip(Vector2 p) { return !g_gx.clipN || CheckCollisionPointRec(p, g_gx.clip[g_gx.clipN - 1]); }

/* ---- mouse ownership ------------------------------------------------ */
static bool g_gxInPopup;
static bool gxHover(Rectangle r) {
    if (!CheckCollisionPointRec(g_gx.mouse, r) || !gxInClip(g_gx.mouse)) return false;
    if (!g_gxInPopup && g_gx.popupRectPrev.width > 0 && CheckCollisionPointRec(g_gx.mouse, g_gx.popupRectPrev)) return false;
    return true;
}
/* press-and-release over the widget = a click */
static bool gxBehave(GxId id, Rectangle r, bool *hovered, bool *held) {
    bool h = gxHover(r) && (g_gx.active == 0 || g_gx.active == id);
    if (h) { g_gx.hot = id; g_gx.overUi = true; }
    bool clicked = false;
    if (h && g_gx.pressed[0]) g_gx.active = id;
    if (g_gx.active == id && g_gx.released[0]) { if (h) clicked = true; g_gx.active = 0; }
    if (hovered) *hovered = h;
    if (held) *held = g_gx.active == id && g_gx.down[0];
    return clicked;
}

/* ---- drawing primitives --------------------------------------------- */
static void drRect(Rectangle r, Color c) { if (r.width > 0 && r.height > 0) DrawRectangleRec(r, c); }
static void drRound(Rectangle r, float rad, Color c) {
    if (r.width <= 0 || r.height <= 0) return;
    float m = fminf(r.width, r.height);
    DrawRectangleRounded(r, m > 0 ? fminf(1.0f, rad * 2 / m) : 0, 8, c);
}
static void drRoundLine(Rectangle r, float rad, Color c) {
    if (r.width <= 0 || r.height <= 0) return;
    float m = fminf(r.width, r.height);
    DrawRectangleRoundedLinesEx(r, m > 0 ? fminf(1.0f, rad * 2 / m) : 0, 8, 1.0f, c);
}
static void drLine(float x0, float y0, float x1, float y1, float th, Color c) { DrawLineEx((Vector2){ x0, y0 }, (Vector2){ x1, y1 }, th, c); }
static void drHLine(float x, float y, float w, Color c) { drRect((Rectangle){ x, y, w, 1 }, c); }
static void drVLine(float x, float y, float h, Color c) { drRect((Rectangle){ x, y, 1, h }, c); }
static Color drFade(Color c, float a) { c.a = (unsigned char)(c.a * clamp01(a)); return c; }
static Color drMix(Color a, Color b, float t) {
    t = clamp01(t);
    return (Color){ (unsigned char)(a.r + (b.r - a.r) * t), (unsigned char)(a.g + (b.g - a.g) * t),
                    (unsigned char)(a.b + (b.b - a.b) * t), (unsigned char)(a.a + (b.a - a.a) * t) };
}
static Rectangle drInset(Rectangle r, float d) { return (Rectangle){ r.x + d, r.y + d, r.width - 2 * d, r.height - 2 * d }; }
static void drShadow(Rectangle r, float rad) {
    for (int i = 3; i >= 1; i--) drRound((Rectangle){ r.x - i + 1, r.y - i + 3, r.width + 2 * i - 2, r.height + 2 * i - 2 }, rad + i, drFade(TH->shadow, 0.35f));
}

/* ---- tooltips -------------------------------------------------------- */
static void gxTip(GxId id, Rectangle r, const char *text) {
    if (!text || !*text || !gxHover(r)) return;
    if (g_gx.tipId != id) { g_gx.tipId = id; g_gx.tipSince = g_gx.now; return; }
    g_gx.tipShown = true;
    if (g_gx.now - g_gx.tipSince > 0.40 && !g_gx.down[0]) {
        snprintf(g_gx.tip, sizeof g_gx.tip, "%s", text);
        g_gx.tipAt = (Vector2){ g_gx.mouse.x + 14, g_gx.mouse.y + 20 };
    }
}
static void gxTipDraw(void) {
    if (!g_gx.tipShown) { g_gx.tip[0] = 0; g_gx.tipId = 0; return; }
    if (!g_gx.tip[0]) return;
    float w = 360, lh = 18;
    float h = txtParaH(F_REG, TH->fsSmall, w, lh, g_gx.tip, 14);
    float tw = 0;
    if (h <= lh) { tw = (float)txtW(F_REG, TH->fsSmall, g_gx.tip); w = tw; }
    Rectangle r = { g_gx.tipAt.x, g_gx.tipAt.y, w + 20, h + 14 };
    if (r.x + r.width > gxScreenW() - 6) r.x = gxScreenW() - 6 - r.width;
    if (r.y + r.height > gxScreenH() - 6) r.y = g_gx.tipAt.y - r.height - 28;
    drShadow(r, 6);
    drRound(r, 6, TH->panel2);
    drRoundLine(r, 6, TH->lineHi);
    txtPara(F_REG, TH->fsSmall, r.x + 10, r.y + 7, w + 1, lh, TH->text, g_gx.tip, 14);
}

/* ---- buttons --------------------------------------------------------- */
typedef enum { BTN_NORMAL, BTN_PRIMARY, BTN_GHOST, BTN_DANGER } GxBtnKind;
static bool gxButtonEx(GxId id, Rectangle r, const char *label, GxBtnKind kind, bool enabled, const char *tip) {
    bool h, held;
    bool clicked = gxBehave(id, r, &h, &held) && enabled;
    if (h && enabled) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
    Color bg, fg = enabled ? TH->text : TH->faint;
    switch (kind) {
    case BTN_PRIMARY: bg = held ? TH->accentDim : h ? drMix(TH->accent, WHITE, 0.15f) : TH->accent; fg = TH->onAccent; break;
    case BTN_GHOST: bg = held ? TH->press : h ? TH->hover : BLANK; break;
    case BTN_DANGER: bg = held ? drMix(TH->bad, BLACK, 0.3f) : h ? drMix(TH->bad, WHITE, 0.1f) : TH->bad; fg = WHITE; break;
    default: bg = held ? TH->press : h ? TH->hover : TH->panel2; break;
    }
    if (!enabled) bg = kind == BTN_GHOST ? BLANK : TH->panel;
    drRound(r, 6, bg);
    if (kind == BTN_NORMAL) drRoundLine(r, 6, h ? TH->lineHi : TH->line);
    txtIn(F_SEMI, TH->fsSmall, r, fg, 1, label);
    if (tip) gxTip(id, r, tip);
    return clicked;
}
static bool gxButton(Rectangle r, const char *label) { return gxButtonEx(gxId(label), r, label, BTN_NORMAL, true, NULL); }
static bool gxButtonP(Rectangle r, const char *label) { return gxButtonEx(gxId(label), r, label, BTN_PRIMARY, true, NULL); }
static bool gxButtonG(Rectangle r, const char *label) { return gxButtonEx(gxId(label), r, label, BTN_GHOST, true, NULL); }
static float gxButtonW(const char *label) { return (float)txtW(F_SEMI, TH->fsSmall, label) + 26; }

/* a clickable line of text (links to genes, variants, loci) */
static bool gxLink(GxId id, float x, float y, int px, StFace f, const char *s, Color c, const char *tip) {
    float w = (float)txtW(f, px, s);
    Rectangle r = { x, y - 2, w, (float)px + 6 };
    bool h;
    bool clicked = gxBehave(id, r, &h, NULL);
    if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
    txt(f, px, x, y, h ? drMix(c, WHITE, 0.3f) : c, s);
    if (h) drHLine(x, y + (float)px + 1, w, drFade(c, 0.7f));
    if (tip) gxTip(id, r, tip);
    return clicked;
}

/* a drawn icon inside an invisible button */
typedef void (*GxIconFn)(Rectangle r, Color c);
static bool gxIconButton(GxId id, Rectangle r, GxIconFn icon, bool on, const char *tip) {
    bool h, held;
    bool clicked = gxBehave(id, r, &h, &held);
    if (on) drRound(r, 6, TH->selDim);
    else if (held) drRound(r, 6, TH->press);
    else if (h) drRound(r, 6, TH->hover);
    if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
    if (icon) icon(drInset(r, 6), on || h ? TH->text : TH->dim);
    if (tip) gxTip(id, r, tip);
    return clicked;
}

static bool gxCheckbox(GxId id, Rectangle r, bool *v, const char *label) {
    bool h;
    bool clicked = gxBehave(id, r, &h, NULL);
    Rectangle b = { r.x, r.y + (r.height - 16) / 2, 16, 16 };
    drRound(b, 4, *v ? TH->accent : (h ? TH->hover : TH->field));
    drRoundLine(b, 4, *v ? TH->accent : TH->lineHi);
    if (*v) {
        drLine(b.x + 3.5f, b.y + 8.5f, b.x + 6.5f, b.y + 11.5f, 2, TH->onAccent);
        drLine(b.x + 6.5f, b.y + 11.5f, b.x + 12.5f, b.y + 4.5f, 2, TH->onAccent);
    }
    if (label) txtIn(F_REG, TH->fsSmall, (Rectangle){ r.x + 24, r.y, r.width - 24, r.height }, h ? TH->text : TH->dim, 0, label);
    if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
    if (clicked) *v = !*v;
    return clicked;
}

/* pills and badges */
static float gxPillW(const char *s, int px) { return (float)txtW(F_SEMI, px, s) + (float)px; }
static float gxPill(float x, float y, const char *s, Color c, int px) {
    float w = gxPillW(s, px), h = (float)px + 8;
    Rectangle r = { x, y, w, h };
    drRound(r, h / 2, drFade(c, 0.18f));
    drRoundLine(r, h / 2, drFade(c, 0.45f));
    txtIn(F_SEMI, px, r, drMix(c, WHITE, 0.25f), 1, s);
    return w;
}
/* a horizontal meter 0..1 */
static void gxMeter(Rectangle r, float v, Color c) {
    drRound(r, r.height / 2, TH->field);
    if (v > 0) drRound((Rectangle){ r.x, r.y, fmaxf(r.height, r.width * clamp01(v)), r.height }, r.height / 2, c);
}

/* segmented control: returns true when the choice changed */
static bool gxSegmented(GxId id, Rectangle r, const char *const *items, int n, int *sel) {
    drRound(r, 7, TH->field);
    drRoundLine(r, 7, TH->line);
    float w = r.width / (float)n;
    bool changed = false;
    for (int i = 0; i < n; i++) {
        Rectangle c = { r.x + w * (float)i + 2, r.y + 2, w - 4, r.height - 4 };
        bool h;
        if (gxBehave(gxHash(id, &i, sizeof i), c, &h, NULL) && *sel != i) { *sel = i; changed = true; }
        if (*sel == i) drRound(c, 5, TH->panel2);
        else if (h) drRound(c, 5, TH->hover);
        if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
        txtIn(*sel == i ? F_SEMI : F_REG, TH->fsSmall, c, *sel == i ? TH->text : TH->dim, 1, items[i]);
    }
    return changed;
}

/* underline tabs: returns true when the tab changed */
static bool gxTabs(GxId id, float x, float y, float h, const char *const *items, int n, int *sel) {
    bool changed = false;
    float cx = x;
    for (int i = 0; i < n; i++) {
        float w = (float)txtW(F_SEMI, TH->fs, items[i]) + 28;
        Rectangle c = { cx, y, w, h };
        bool hv;
        if (gxBehave(gxHash(id, &i, sizeof i), c, &hv, NULL) && *sel != i) { *sel = i; changed = true; }
        if (hv) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
        txtIn(*sel == i ? F_SEMI : F_REG, TH->fs, c, *sel == i ? TH->text : hv ? TH->text : TH->dim, 1, items[i]);
        if (*sel == i) drRound((Rectangle){ cx + 10, y + h - 3, w - 20, 3 }, 1.5f, TH->accent);
        cx += w;
    }
    return changed;
}

/* ---- text field --------------------------------------------------------
   Returns 1 on Enter, 2 while the text changes (live), 0 otherwise. Esc
   clears a focused field's typing and leaves it. */
static int gxUtf8Prev(const char *s, int i) { do { i--; } while (i > 0 && ((unsigned char)s[i] & 0xC0) == 0x80); return i < 0 ? 0 : i; }
static int gxUtf8Next(const char *s, int i, int n) { do { i++; } while (i < n && ((unsigned char)s[i] & 0xC0) == 0x80); return i > n ? n : i; }
static void gxEditDelSel(void) {
    int a = g_gx.editCur < g_gx.editSel ? g_gx.editCur : g_gx.editSel, b = g_gx.editCur ^ g_gx.editSel ^ a;
    if (a == b) return;
    memmove(g_gx.edit + a, g_gx.edit + b, (size_t)(g_gx.editLen - b + 1));
    g_gx.editLen -= b - a;
    g_gx.editCur = g_gx.editSel = a;
}
static void gxEditInsert(const char *s) {
    gxEditDelSel();
    int n = (int)strlen(s);
    if (g_gx.editLen + n >= (int)sizeof g_gx.edit - 1) return;
    memmove(g_gx.edit + g_gx.editCur + n, g_gx.edit + g_gx.editCur, (size_t)(g_gx.editLen - g_gx.editCur + 1));
    memcpy(g_gx.edit + g_gx.editCur, s, (size_t)n);
    g_gx.editLen += n;
    g_gx.editCur += n;
    g_gx.editSel = g_gx.editCur;
}
static int gxCharAt(const char *s, int px, float x) {
    int n = (int)strlen(s), best = 0;
    float bestD = 1e9f;
    for (int i = 0; i <= n; i = i < n ? gxUtf8Next(s, i, n) : n + 1) {
        float d = fabsf((float)txtWn(F_REG, px, s, i) - x);
        if (d < bestD) { bestD = d; best = i; }
        if (i == n) break;
    }
    return best;
}
static void gxFocus(GxId id, const char *text, bool selectAll) {
    g_gx.editId = id;
    snprintf(g_gx.edit, sizeof g_gx.edit, "%s", text);
    g_gx.editLen = (int)strlen(g_gx.edit);
    g_gx.editCur = g_gx.editLen;
    g_gx.editSel = selectAll ? 0 : g_gx.editLen;
    g_gx.editScroll = 0;
    g_gx.keysTaken = true;
}
static void gxSearchIcon(Rectangle r, Color c);
static int gxTextField(GxId id, Rectangle r, char *buf, size_t cap, const char *hint, bool searchIcon) {
    bool h;
    gxBehave(id, r, &h, NULL);
    bool focused = g_gx.editId == id;
    if (h) g_gx.cursor = MOUSE_CURSOR_IBEAM;
    float left = searchIcon ? 30 : 10;
    if (!focused && h && g_gx.pressed[0]) {
        gxFocus(id, buf, false);
        focused = true;
        g_gx.editCur = g_gx.editSel = gxCharAt(g_gx.edit, TH->fs, g_gx.mouse.x - (r.x + left));
    }
    int result = 0;
    if (focused) {
        bool ctrl = IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL);
        bool shift = IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT);
        float x0 = r.x + left - g_gx.editScroll;
        if (g_gx.pressed[0] && h) {
            g_gx.editCur = gxCharAt(g_gx.edit, TH->fs, g_gx.mouse.x - x0);
            if (!shift) g_gx.editSel = g_gx.editCur;
            if (g_gx.dbl) { g_gx.editSel = 0; g_gx.editCur = g_gx.editLen; }
        } else if (g_gx.down[0] && g_gx.active == id) {
            g_gx.editCur = gxCharAt(g_gx.edit, TH->fs, g_gx.mouse.x - x0);
        }
        for (int i = 0; i < g_gx.charN; i++) {
            int c = g_gx.chars[i];
            if (c < 32) continue;
            int bl = 0;
            const char *u = CodepointToUTF8(c, &bl);
            char tmp[8] = { 0 };
            memcpy(tmp, u, (size_t)(bl < 7 ? bl : 7));
            gxEditInsert(tmp);
            result = 2;
        }
        if (IsKeyPressed(KEY_BACKSPACE) || IsKeyPressedRepeat(KEY_BACKSPACE)) {
            if (g_gx.editCur != g_gx.editSel) gxEditDelSel();
            else if (g_gx.editCur > 0) { g_gx.editSel = gxUtf8Prev(g_gx.edit, g_gx.editCur); gxEditDelSel(); }
            result = 2;
        }
        if (IsKeyPressed(KEY_DELETE) || IsKeyPressedRepeat(KEY_DELETE)) {
            if (g_gx.editCur != g_gx.editSel) gxEditDelSel();
            else if (g_gx.editCur < g_gx.editLen) { g_gx.editSel = gxUtf8Next(g_gx.edit, g_gx.editCur, g_gx.editLen); gxEditDelSel(); }
            result = 2;
        }
        if (IsKeyPressed(KEY_LEFT) || IsKeyPressedRepeat(KEY_LEFT)) { g_gx.editCur = gxUtf8Prev(g_gx.edit, g_gx.editCur); if (!shift) g_gx.editSel = g_gx.editCur; }
        if (IsKeyPressed(KEY_RIGHT) || IsKeyPressedRepeat(KEY_RIGHT)) { g_gx.editCur = gxUtf8Next(g_gx.edit, g_gx.editCur, g_gx.editLen); if (!shift) g_gx.editSel = g_gx.editCur; }
        if (IsKeyPressed(KEY_HOME)) { g_gx.editCur = 0; if (!shift) g_gx.editSel = 0; }
        if (IsKeyPressed(KEY_END)) { g_gx.editCur = g_gx.editLen; if (!shift) g_gx.editSel = g_gx.editCur; }
        if (ctrl && IsKeyPressed(KEY_A)) { g_gx.editSel = 0; g_gx.editCur = g_gx.editLen; }
        if (ctrl && (IsKeyPressed(KEY_C) || IsKeyPressed(KEY_X)) && g_gx.editCur != g_gx.editSel) {
            int a = g_gx.editCur < g_gx.editSel ? g_gx.editCur : g_gx.editSel, b = g_gx.editCur ^ g_gx.editSel ^ a;
            char tmp[512];
            memcpy(tmp, g_gx.edit + a, (size_t)(b - a));
            tmp[b - a] = 0;
            SetClipboardText(tmp);
            if (IsKeyPressed(KEY_X)) { gxEditDelSel(); result = 2; }
        }
        if (ctrl && IsKeyPressed(KEY_V)) {
            const char *cb = GetClipboardText();
            if (cb) {
                char tmp[512];
                snprintf(tmp, sizeof tmp, "%s", cb);
                for (char *c = tmp; *c; c++) if (*c == '\n' || *c == '\r' || *c == '\t') *c = ' ';
                gxEditInsert(tmp);
                result = 2;
            }
        }
        bool enter = IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_KP_ENTER);
        bool esc = IsKeyPressed(KEY_ESCAPE);
        bool lost = g_gx.pressed[0] && !h && !g_gxInPopup;
        if (result == 2 || enter) snprintf(buf, cap, "%s", g_gx.edit);
        if (enter) { result = 1; }
        if (esc) { g_gx.editId = 0; buf[0] = 0; result = 2; }
        else if (lost) g_gx.editId = 0;
    }
    drRound(r, 7, focused ? TH->field : (h ? drMix(TH->field, TH->hover, 0.5f) : TH->field));
    drRoundLine(r, 7, focused ? TH->accent : (h ? TH->lineHi : TH->line));
    if (searchIcon) gxSearchIcon((Rectangle){ r.x + 9, r.y + (r.height - 14) / 2, 14, 14 }, TH->faint);
    const char *shown = focused ? g_gx.edit : buf;
    Rectangle in = { r.x + left, r.y, r.width - left - 8, r.height };
    gxClipPush(in);
    float ty = r.y + (r.height - (float)TH->fs) / 2 - 1;
    if (focused) {
        float cx = (float)txtWn(F_REG, TH->fs, g_gx.edit, g_gx.editCur);
        if (cx - g_gx.editScroll > in.width - 4) g_gx.editScroll = cx - in.width + 4;
        if (cx - g_gx.editScroll < 0) g_gx.editScroll = cx;
        float x0 = in.x - g_gx.editScroll;
        if (g_gx.editSel != g_gx.editCur) {
            int a = g_gx.editCur < g_gx.editSel ? g_gx.editCur : g_gx.editSel, b = g_gx.editCur ^ g_gx.editSel ^ a;
            float xa = (float)txtWn(F_REG, TH->fs, g_gx.edit, a), xb = (float)txtWn(F_REG, TH->fs, g_gx.edit, b);
            drRect((Rectangle){ x0 + xa, r.y + 5, xb - xa, r.height - 10 }, TH->sel);
        }
        txt(F_REG, TH->fs, x0, ty, TH->text, shown);
        if (fmod(g_gx.now, 1.0) < 0.6) drRect((Rectangle){ x0 + cx, r.y + 6, 1.5f, r.height - 12 }, TH->text);
    } else if (*shown) {
        txt(F_REG, TH->fs, in.x, ty, TH->text, shown);
    } else if (hint) {
        txt(F_REG, TH->fs, in.x, ty, TH->faint, hint);
    }
    gxClipPop();
    return result;
}

/* ---- pop-ups ---------------------------------------------------------
   One pop-up at a time. The widget that owns it sets a draw callback; the
   frame draws it last, over everything, and widgets under its rectangle do
   not take the mouse. A click outside closes it. */
static void gxPopupOpen(GxId id, Rectangle anchor, Rectangle rect, void (*draw)(Rectangle, void *), void *u) {
    g_gx.popup = id;
    g_gx.popupAnchor = anchor;
    g_gx.popupRect = rect;
    g_gx.popupDraw = draw;
    g_gx.popupUser = u;
    g_gx.popupOpening = true;
}
static void gxPopupClose(void) { g_gx.popup = 0; g_gx.popupDraw = NULL; }
static void gxPopupFrame(void) {
    if (!g_gx.popup || !g_gx.popupDraw) return;
    Rectangle r = g_gx.popupRect;
    if (r.y + r.height > gxScreenH() - 8) r.y = fmaxf(8, gxScreenH() - 8 - r.height);
    if (r.x + r.width > gxScreenW() - 8) r.x = fmaxf(8, gxScreenW() - 8 - r.width);
    g_gx.popupRect = r;
    int savedClip = g_gx.clipN;
    g_gx.clipN = 0;
    gxClipApply();
    g_gxInPopup = true;
    drShadow(r, 8);
    drRound(r, 8, TH->panel2);
    drRoundLine(r, 8, TH->lineHi);
    gxClipPush(r);
    g_gx.popupDraw(r, g_gx.popupUser);
    gxClipPop();
    g_gxInPopup = false;
    g_gx.clipN = savedClip;
    if (!g_gx.popupOpening && g_gx.pressed[0] && !CheckCollisionPointRec(g_gx.mouse, r) && !CheckCollisionPointRec(g_gx.mouse, g_gx.popupAnchor))
        gxPopupClose();
    if (IsKeyPressed(KEY_ESCAPE)) gxPopupClose();
    g_gx.popupOpening = false;
}

/* combo box: a list in a pop-up */
typedef struct { const char *const *items; int n; int *sel; bool changed; GxId id; } GxComboQ;
static GxComboQ g_gxCombo;
static void gxComboDraw(Rectangle r, void *u) {
    GxComboQ *q = u;
    float rh = 28;
    for (int i = 0; i < q->n; i++) {
        Rectangle row = { r.x + 4, r.y + 4 + rh * (float)i, r.width - 8, rh };
        bool h;
        if (gxBehave(gxHash(q->id, &i, sizeof i), row, &h, NULL)) { *q->sel = i; q->changed = true; gxPopupClose(); }
        if (h || *q->sel == i) drRound(row, 5, h ? TH->hover : TH->selDim);
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 10, row.y, row.width - 20, row.height }, TH->text, 0, q->items[i]);
    }
}
static bool gxCombo(GxId id, Rectangle r, const char *const *items, int n, int *sel) {
    bool h;
    bool clicked = gxBehave(id, r, &h, NULL);
    if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
    bool open = g_gx.popup == id;
    drRound(r, 7, h || open ? TH->hover : TH->panel2);
    drRoundLine(r, 7, open ? TH->accent : TH->line);
    txtIn(F_REG, TH->fsSmall, (Rectangle){ r.x + 10, r.y, r.width - 30, r.height }, TH->text, 0, *sel >= 0 && *sel < n ? items[*sel] : "");
    float cx = r.x + r.width - 16, cy = r.y + r.height / 2;
    drLine(cx - 4, cy - 2, cx, cy + 2, 1.6f, TH->dim);
    drLine(cx, cy + 2, cx + 4, cy - 2, 1.6f, TH->dim);
    bool changed = false;
    if (g_gxCombo.id == id && g_gxCombo.changed) { changed = true; g_gxCombo.changed = false; }
    if (clicked) {
        if (open) gxPopupClose();
        else {
            g_gxCombo = (GxComboQ){ items, n, sel, false, id };
            gxPopupOpen(id, r, (Rectangle){ r.x, r.y + r.height + 4, fmaxf(r.width, 180), 8 + 28 * (float)n }, gxComboDraw, &g_gxCombo);
        }
    }
    return changed;
}

/* ---- scroll area -------------------------------------------------------
   gxScrollBegin clips to r and returns the content's top y; the caller lays
   out from there and passes the content height to gxScrollEnd, which draws
   the scroll bar and takes the wheel. */
typedef struct { float pos, target, contentH, viewH; GxId id; bool dragging; float dragOff; } GxScroll;
static float gxScrollBegin(GxScroll *s, Rectangle r) {
    s->viewH = r.height;
    float maxS = fmaxf(0, s->contentH - r.height);
    if (gxHover(r) && g_gx.wheel != 0 && !g_gx.popupRectPrev.width) s->target -= g_gx.wheel * 90;
    s->target = clampf(s->target, 0, maxS);
    s->pos += (s->target - s->pos) * 0.35f;
    if (fabsf(s->target - s->pos) < 0.5f) s->pos = s->target;
    s->pos = clampf(s->pos, 0, maxS);
    gxClipPush(r);
    return r.y - floorf(s->pos);
}
static void gxScrollEnd(GxScroll *s, Rectangle r, float contentH) {
    gxClipPop();
    s->contentH = contentH;
    if (contentH <= r.height + 1) { s->pos = s->target = 0; return; }
    float trackH = r.height - 8, barH = fmaxf(30, trackH * r.height / contentH);
    float maxS = contentH - r.height;
    float barY = r.y + 4 + (trackH - barH) * (s->pos / maxS);
    Rectangle bar = { r.x + r.width - 8, barY, 5, barH };
    Rectangle hit = { r.x + r.width - 12, r.y, 12, r.height };
    if (!s->id) s->id = gxIdN((uint64_t)(uintptr_t)s);
    bool h, held;
    gxBehave(s->id, hit, &h, &held);
    if (h && g_gx.pressed[0]) {
        s->dragging = true;
        s->dragOff = CheckCollisionPointRec(g_gx.mouse, (Rectangle){ bar.x - 4, bar.y, 13, bar.height }) ? g_gx.mouse.y - barY : barH / 2;
    }
    if (!g_gx.down[0]) s->dragging = false;
    if (s->dragging) {
        float t = (g_gx.mouse.y - s->dragOff - (r.y + 4)) / fmaxf(1, trackH - barH);
        s->target = s->pos = clampf(t, 0, 1) * maxS;
    }
    drRound(bar, 2.5f, h || s->dragging ? TH->lineHi : drFade(TH->lineHi, 0.6f));
}
static void gxScrollTo(GxScroll *s, float y) { s->target = fmaxf(0, y); }

/* ---- the table ---------------------------------------------------------
   Virtual rows: only the visible ones are drawn (the caller's cell callback
   gets the row index). Headers sort when a column allows it. */
typedef struct { const char *title; float w; int align; bool sortable; const char *tip; } GxCol;
typedef struct {
    GxScroll scroll;
    int sel;            /* selected row, -1 none */
    int sortCol;        /* -1 none */
    bool sortDesc;
    bool sortChanged;
    int hover;
} GxTable;
typedef void (*GxCellFn)(int row, int col, Rectangle cell, bool selected, void *u);
/* returns the row clicked this frame (-1 none); *dbl set on a double click */
static int gxTable(GxId id, GxTable *t, Rectangle r, const GxCol *cols, int ncol, int nrow, float rowH, GxCellFn cell, void *u, bool *dbl) {
    float hh = 32;
    float fixed = 0, flex = 0;
    for (int c = 0; c < ncol; c++) { if (cols[c].w > 0) fixed += cols[c].w; else flex += -cols[c].w; }
    float avail = r.width - 10 - fixed;
    float cw[64];
    for (int c = 0; c < ncol && c < 64; c++) cw[c] = cols[c].w > 0 ? cols[c].w : fmaxf(60, avail * (-cols[c].w) / fmaxf(flex, 1));
    /* header */
    drRect((Rectangle){ r.x, r.y, r.width, hh }, TH->panel);
    drHLine(r.x, r.y + hh - 1, r.width, TH->line);
    float x = r.x + 4;
    t->sortChanged = false;
    for (int c = 0; c < ncol; c++) {
        Rectangle hc = { x, r.y, cw[c], hh };
        bool h = false;
        if (cols[c].sortable) {
            if (gxBehave(gxHash(id, &c, sizeof c), hc, &h, NULL)) {
                if (t->sortCol == c) t->sortDesc = !t->sortDesc;
                else { t->sortCol = c; t->sortDesc = cols[c].align == 2; }
                t->sortChanged = true;
            }
            if (h) g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
        }
        Color hc_ = t->sortCol == c ? TH->text : h ? TH->text : TH->dim;
        Rectangle tr = { x + 8, r.y, cw[c] - 16 - (t->sortCol == c ? 12 : 0), hh };
        txtIn(F_SEMI, TH->fsTiny, tr, hc_, cols[c].align, cols[c].title);
        if (t->sortCol == c) {
            float ax = cols[c].align == 2 ? x + cw[c] - 14 : x + 12 + (float)txtW(F_SEMI, TH->fsTiny, cols[c].title), ay = r.y + hh / 2;
            if (cols[c].align == 2) ax = x + cw[c] - 8 - (float)txtW(F_SEMI, TH->fsTiny, cols[c].title) - 12, tr.x += 0;
            if (t->sortDesc) DrawTriangle((Vector2){ ax - 4, ay - 2 }, (Vector2){ ax, ay + 3 }, (Vector2){ ax + 4, ay - 2 }, TH->accent);
            else DrawTriangle((Vector2){ ax, ay - 3 }, (Vector2){ ax - 4, ay + 2 }, (Vector2){ ax + 4, ay + 2 }, TH->accent);
        }
        if (cols[c].tip) gxTip(gxHash(id ^ 77, &c, sizeof c), hc, cols[c].tip);
        x += cw[c];
    }
    /* body */
    Rectangle body = { r.x, r.y + hh, r.width, r.height - hh };
    t->scroll.contentH = (float)nrow * rowH;
    float top = gxScrollBegin(&t->scroll, body);
    int first = (int)floorf((body.y - top) / rowH);
    if (first < 0) first = 0;
    int last = (int)ceilf((body.y + body.height - top) / rowH);
    if (last > nrow) last = nrow;
    int clicked = -1;
    if (dbl) *dbl = false;
    t->hover = -1;
    for (int i = first; i < last; i++) {
        Rectangle row = { r.x, top + (float)i * rowH, r.width - 10, rowH };
        bool h;
        GxId rid = gxHash(id ^ 0x5151, &i, sizeof i);
        if (gxBehave(rid, row, &h, NULL)) { clicked = i; t->sel = i; }
        if (h && g_gx.dbl && dbl) { *dbl = true; clicked = i; t->sel = i; }
        if (h) t->hover = i;
        bool sel = t->sel == i;
        if (sel) drRect(row, TH->selDim);
        else if (h) drRect(row, TH->hover);
        else if (i & 1) drRect(row, TH->zebra);
        float cx = r.x + 4;
        for (int c = 0; c < ncol; c++) {
            Rectangle cr = { cx + 8, row.y, cw[c] - 16, rowH };
            cell(i, c, cr, sel, u);
            cx += cw[c];
        }
    }
    gxScrollEnd(&t->scroll, body, (float)nrow * rowH);
    return clicked;
}

/* ---- cards and headings ---------------------------------------------- */
static void gxCard(Rectangle r) {
    drRound(r, 10, TH->card);
    drRoundLine(r, 10, TH->line);
}
static void gxCardTitle(Rectangle r, const char *title, const char *sub) {
    txt(F_SEMI, TH->fsH3, r.x + 18, r.y + 16, TH->text, title);
    if (sub) txt(F_REG, TH->fsSmall, r.x + 18, r.y + 16 + (float)TH->fsH3 + 6, TH->dim, sub);
}
static void gxSectionTitle(float x, float y, const char *title) {
    txt(F_SEMI, TH->fsTiny, x, y, TH->faint, title);
}

/* the frame's end: pop-up, tooltip, cursor */
static void gxEndFrame(void) {
    g_gx.clipN = 0;
    gxClipApply();
    gxPopupFrame();
    gxTipDraw();
    if (!g_gx.down[0] && !g_gx.down[1] && !g_gx.down[2]) g_gx.active = 0;
    SetMouseCursor(g_gx.cursor);
    rlDrawRenderBatchActive();
    rlMatrixMode(RL_MODELVIEW);
    rlLoadIdentity();
}

/* ---- icons (drawn, so they scale and theme with everything else) ------- */
static void gxSearchIcon(Rectangle r, Color c) {
    float cx = r.x + r.width * 0.42f, cy = r.y + r.height * 0.42f, rad = r.width * 0.32f;
    DrawRing((Vector2){ cx, cy }, rad - 1.2f, rad + 0.4f, 0, 360, 24, c);
    drLine(cx + rad * 0.7f, cy + rad * 0.7f, r.x + r.width - 1, r.y + r.height - 1, 2, c);
}
static void icHome(Rectangle r, Color c) {
    float x = r.x, y = r.y, w = r.width, h = r.height;
    drLine(x + 1, y + h * 0.45f, x + w / 2, y + 1, 1.8f, c);
    drLine(x + w / 2, y + 1, x + w - 1, y + h * 0.45f, 1.8f, c);
    DrawRectangleLinesEx((Rectangle){ x + w * 0.18f, y + h * 0.42f, w * 0.64f, h * 0.55f }, 1.8f, c);
}
static void icHelix(Rectangle r, Color c) {
    for (int i = 0; i < 2; i++) {
        Vector2 prev = { 0 };
        for (int k = 0; k <= 16; k++) {
            float t = (float)k / 16.0f;
            float y = r.y + t * r.height;
            float x = r.x + r.width / 2 + sinf(t * 6.2832f + (float)i * 3.1416f) * r.width * 0.38f;
            Vector2 p = { x, y };
            if (k) drLine(prev.x, prev.y, p.x, p.y, 1.7f, c);
            prev = p;
        }
    }
    for (int k = 1; k < 4; k++) {
        float t = (float)k / 4.0f, y = r.y + t * r.height;
        float dx = sinf(t * 6.2832f) * r.width * 0.38f;
        drLine(r.x + r.width / 2 - dx, y, r.x + r.width / 2 + dx, y, 1.0f, drFade(c, 0.6f));
    }
}
static void icList(Rectangle r, Color c) {
    for (int i = 0; i < 3; i++) {
        float y = r.y + r.height * (0.2f + 0.3f * (float)i);
        DrawCircleV((Vector2){ r.x + 2, y }, 1.6f, c);
        drLine(r.x + 6, y, r.x + r.width, y, 1.8f, c);
    }
}
static void icHeart(Rectangle r, Color c) {
    float cx = r.x + r.width / 2, w = r.width, y = r.y;
    DrawCircleV((Vector2){ cx - w * 0.22f, y + r.height * 0.32f }, w * 0.25f, c);
    DrawCircleV((Vector2){ cx + w * 0.22f, y + r.height * 0.32f }, w * 0.25f, c);
    DrawTriangle((Vector2){ r.x + w * 0.02f, y + r.height * 0.42f }, (Vector2){ cx, r.y + r.height }, (Vector2){ r.x + w * 0.98f, y + r.height * 0.42f }, c);
}
static void icPill(Rectangle r, Color c) {
    Rectangle p = { r.x, r.y + r.height * 0.3f, r.width, r.height * 0.4f };
    DrawRectangleRoundedLinesEx(p, 1.0f, 8, 1.8f, c);
    drRect((Rectangle){ p.x + p.width / 2, p.y + 1, p.width / 2 - 2, p.height - 2 }, drFade(c, 0.6f));
}
static void icPerson(Rectangle r, Color c) {
    DrawRing((Vector2){ r.x + r.width / 2, r.y + r.height * 0.28f }, r.width * 0.16f, r.width * 0.16f + 1.8f, 0, 360, 20, c);
    DrawRing((Vector2){ r.x + r.width / 2, r.y + r.height * 1.02f }, r.width * 0.42f, r.width * 0.42f + 1.8f, 180, 360, 24, c);
}
static void icGlobe(Rectangle r, Color c) {
    Vector2 m = { r.x + r.width / 2, r.y + r.height / 2 };
    float rad = r.width / 2;
    DrawRing(m, rad - 1.6f, rad, 0, 360, 28, c);
    DrawEllipseLines((int)m.x, (int)m.y, rad * 0.45f, rad, c);
    drLine(m.x - rad, m.y, m.x + rad, m.y, 1.4f, c);
}
static void icTree(Rectangle r, Color c) {
    Vector2 top = { r.x + r.width / 2, r.y + 2 }, a = { r.x + 3, r.y + r.height - 2 }, b = { r.x + r.width - 3, r.y + r.height - 2 };
    Vector2 mid = { r.x + r.width / 2, r.y + r.height * 0.55f };
    drLine(top.x, top.y, mid.x, mid.y, 1.8f, c);
    drLine(mid.x, mid.y, a.x, a.y, 1.8f, c);
    drLine(mid.x, mid.y, b.x, b.y, 1.8f, c);
    DrawCircleV(top, 2.4f, c);
    DrawCircleV(a, 2.4f, c);
    DrawCircleV(b, 2.4f, c);
}
static void icChart(Rectangle r, Color c) {
    for (int i = 0; i < 4; i++) {
        float h = r.height * (0.3f + 0.2f * (float)((i * 7) % 4));
        drRect((Rectangle){ r.x + (float)i * r.width / 4 + 1, r.y + r.height - h, r.width / 4 - 2.5f, h }, c);
    }
}
static void icFile(Rectangle r, Color c) {
    DrawRectangleLinesEx((Rectangle){ r.x + r.width * 0.15f, r.y, r.width * 0.7f, r.height }, 1.8f, c);
    for (int i = 0; i < 3; i++) drLine(r.x + r.width * 0.3f, r.y + r.height * (0.3f + 0.2f * (float)i), r.x + r.width * 0.7f, r.y + r.height * (0.3f + 0.2f * (float)i), 1.4f, c);
}
static void icInfo(Rectangle r, Color c) {
    Vector2 m = { r.x + r.width / 2, r.y + r.height / 2 };
    DrawRing(m, r.width / 2 - 1.6f, r.width / 2, 0, 360, 28, c);
    drRect((Rectangle){ m.x - 1, m.y - 1, 2, r.height * 0.32f }, c);
    DrawCircleV((Vector2){ m.x, m.y - r.height * 0.2f }, 1.4f, c);
}
static void icShield(Rectangle r, Color c) {
    Vector2 a = { r.x + 2, r.y + 2 }, b = { r.x + r.width - 2, r.y + 2 }, d = { r.x + r.width / 2, r.y + r.height };
    drLine(a.x, a.y, b.x, b.y, 1.8f, c);
    drLine(a.x, a.y, a.x, r.y + r.height * 0.5f, 1.8f, c);
    drLine(b.x, b.y, b.x, r.y + r.height * 0.5f, 1.8f, c);
    drLine(a.x, r.y + r.height * 0.5f, d.x, d.y, 1.8f, c);
    drLine(b.x, r.y + r.height * 0.5f, d.x, d.y, 1.8f, c);
}
static void icDots(Rectangle r, Color c) {
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
        DrawCircleV((Vector2){ r.x + r.width * (0.2f + 0.3f * (float)i), r.y + r.height * (0.2f + 0.3f * (float)j) }, 1.7f, c);
}
static void icStar(Rectangle r, Color c) {
    Vector2 m = { r.x + r.width / 2, r.y + r.height / 2 + 1 };
    Vector2 p[10];
    for (int i = 0; i < 10; i++) {
        float a = -PI / 2 + (float)i * PI / 5, rr = (i & 1) ? r.width * 0.22f : r.width * 0.5f;
        p[i] = (Vector2){ m.x + cosf(a) * rr, m.y + sinf(a) * rr };
    }
    for (int i = 0; i < 10; i++) drLine(p[i].x, p[i].y, p[(i + 1) % 10].x, p[(i + 1) % 10].y, 1.6f, c);
}
static void icGem(Rectangle r, Color c) {
    float l = r.x + 1, rr = r.x + r.width - 1, t = r.y + r.height * 0.18f, m = r.y + r.height * 0.42f, b = r.y + r.height - 1;
    float a = r.x + r.width * 0.3f, d = r.x + r.width * 0.7f, cx = r.x + r.width / 2;
    drLine(a, t, d, t, 1.6f, c);
    drLine(a, t, l, m, 1.6f, c);
    drLine(d, t, rr, m, 1.6f, c);
    drLine(l, m, rr, m, 1.4f, c);
    drLine(l, m, cx, b, 1.6f, c);
    drLine(rr, m, cx, b, 1.6f, c);
    drLine(a, t, cx, m, 1.2f, c);
    drLine(d, t, cx, m, 1.2f, c);
}
static void icRarity(Rectangle r, Color c) {
    Vector2 m = { r.x + r.width / 2, r.y + r.height / 2 };
    for (int i = 0; i < 8; i++) {
        float a = (float)i * PI / 4;
        DrawCircleV((Vector2){ m.x + cosf(a) * r.width * 0.42f, m.y + sinf(a) * r.height * 0.42f }, 1.3f, drFade(c, 0.55f));
    }
    DrawCircleV(m, r.width * 0.2f, c);
}
static void icBack(Rectangle r, Color c) {
    drLine(r.x + r.width * 0.65f, r.y + 2, r.x + r.width * 0.3f, r.y + r.height / 2, 2, c);
    drLine(r.x + r.width * 0.3f, r.y + r.height / 2, r.x + r.width * 0.65f, r.y + r.height - 2, 2, c);
}
static void icFwd(Rectangle r, Color c) {
    drLine(r.x + r.width * 0.35f, r.y + 2, r.x + r.width * 0.7f, r.y + r.height / 2, 2, c);
    drLine(r.x + r.width * 0.7f, r.y + r.height / 2, r.x + r.width * 0.35f, r.y + r.height - 2, 2, c);
}
static void icClose(Rectangle r, Color c) {
    drLine(r.x + 2, r.y + 2, r.x + r.width - 2, r.y + r.height - 2, 1.8f, c);
    drLine(r.x + r.width - 2, r.y + 2, r.x + 2, r.y + r.height - 2, 1.8f, c);
}
static void icCopy(Rectangle r, Color c) {
    DrawRectangleLinesEx((Rectangle){ r.x + r.width * 0.3f, r.y, r.width * 0.7f, r.height * 0.7f }, 1.6f, c);
    DrawRectangleLinesEx((Rectangle){ r.x, r.y + r.height * 0.3f, r.width * 0.7f, r.height * 0.7f }, 1.6f, c);
}
static void icZoomIn(Rectangle r, Color c) {
    gxSearchIcon(r, c);
    float cx = r.x + r.width * 0.42f, cy = r.y + r.height * 0.42f;
    drLine(cx - 3, cy, cx + 3, cy, 1.5f, c);
    drLine(cx, cy - 3, cx, cy + 3, 1.5f, c);
}
static void icZoomOut(Rectangle r, Color c) {
    gxSearchIcon(r, c);
    float cx = r.x + r.width * 0.42f, cy = r.y + r.height * 0.42f;
    drLine(cx - 3, cy, cx + 3, cy, 1.5f, c);
}
