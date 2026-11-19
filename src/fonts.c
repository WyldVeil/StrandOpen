/* ========================================================================
   src/fonts.c - crisp text for a desktop app.

   The engine's drawTextSmart samples one 32 px atlas for every size, which
   suits a game's menus but leaves 13 px table text soft. Strand loads each
   face (Inter Regular and SemiBold, JetBrains Mono) once per pixel size it
   is asked for, rasterised at the physical size (logical size times the UI
   scale), and draws at whole pixels. The glyph set covers Latin-1, Latin
   Extended-A, Greek and the punctuation the databases use (condition names
   such as Sjogren's with its accent, alpha/beta in trait names).
   ======================================================================== */

static const unsigned char g_fontInterReg[] = {
#embed "../assets/fonts/Inter-Regular.ttf"
};
static const unsigned char g_fontInterSemi[] = {
#embed "../assets/fonts/Inter-SemiBold.ttf"
};
static const unsigned char g_fontMonoData[] = {
#embed "../assets/fonts/JetBrainsMono-Regular.ttf"
};

#define ST_FONT_SIZES 48
typedef struct { int px; Font f; bool ok; } StFontSlot;
static StFontSlot g_fonts[F__N][ST_FONT_SIZES];
static float g_fontScale = 1.0f;
static int *g_fontCps;
static int g_fontCpN;

static void fontsInit(void) {
    static const int ranges[][2] = {
        { 0x20, 0x7E }, { 0xA0, 0x17F }, { 0x391, 0x3C9 }, { 0x2010, 0x2027 }, { 0x2030, 0x203A },
        { 0x2190, 0x2193 }, { 0x2212, 0x2212 }, { 0x2264, 0x2265 }, { 0x2248, 0x2248 }, { 0x2260, 0x2260 },
        { 0x221E, 0x221E }, { 0x2713, 0x2713 }, { 0x25B2, 0x25BC }, { 0x25CF, 0x25CF }, { 0x2022, 0x2022 },
    };
    int n = 0;
    for (int i = 0; i < (int)NELEMS(ranges); i++) n += ranges[i][1] - ranges[i][0] + 1;
    g_fontCps = xalloc(sizeof(int) * (size_t)n);
    for (int i = 0; i < (int)NELEMS(ranges); i++)
        for (int c = ranges[i][0]; c <= ranges[i][1]; c++) g_fontCps[g_fontCpN++] = c;
}

static void fontsReload(float scale) {
    for (int f = 0; f < F__N; f++)
        for (int i = 0; i < ST_FONT_SIZES; i++)
            if (g_fonts[f][i].ok) { UnloadFont(g_fonts[f][i].f); g_fonts[f][i].ok = false; g_fonts[f][i].px = 0; }
    g_fontScale = scale;
}

static Font *fontAt(StFace face, int px) {
    if (px < 6) px = 6;
    if (px > 160) px = 160;
    StFontSlot *slots = g_fonts[face];
    for (int i = 0; i < ST_FONT_SIZES; i++)
        if (slots[i].ok && slots[i].px == px) return &slots[i].f;
    int free_ = -1;
    for (int i = 0; i < ST_FONT_SIZES; i++) if (!slots[i].ok) { free_ = i; break; }
    if (free_ < 0) free_ = 0, UnloadFont(slots[0].f), slots[0].ok = false;
    const unsigned char *data = face == F_REG ? g_fontInterReg : face == F_SEMI ? g_fontInterSemi : g_fontMonoData;
    int len = face == F_REG ? (int)sizeof g_fontInterReg : face == F_SEMI ? (int)sizeof g_fontInterSemi : (int)sizeof g_fontMonoData;
    int phys = (int)lroundf((float)px * g_fontScale);
    if (phys < 6) phys = 6;
    Font f = LoadFontFromMemory(".ttf", data, len, phys, g_fontCps, g_fontCpN);
    SetTextureFilter(f.texture, g_fontScale == 1.0f ? TEXTURE_FILTER_POINT : TEXTURE_FILTER_BILINEAR);
    slots[free_] = (StFontSlot){ .px = px, .f = f, .ok = true };
    return &slots[free_].f;
}

static int txtWn(StFace face, int px, const char *s, int n) {
    if (!s || n <= 0) return 0;
    Font *f = fontAt(face, px);
    float scale = (float)px / (float)f->baseSize, w = 0;
    for (int i = 0; i < n && s[i];) {
        int bytes = 0;
        int cp = GetCodepointNext(s + i, &bytes);
        if (bytes <= 0) bytes = 1;
        int gi = GetGlyphIndex(*f, cp);
        w += f->glyphs[gi].advanceX ? (float)f->glyphs[gi].advanceX * scale : f->recs[gi].width * scale;
        i += bytes;
    }
    return (int)ceilf(w);
}
static int txtW(StFace face, int px, const char *s) { return s ? txtWn(face, px, s, (int)strlen(s)) : 0; }

static void txtn(StFace face, int px, float x, float y, Color c, const char *s, int n) {
    if (!s || n <= 0 || c.a == 0) return;
    Font *f = fontAt(face, px);
    float scale = (float)px / (float)f->baseSize;
    float cx = floorf(x), cy = floorf(y);
    for (int i = 0; i < n && s[i];) {
        int bytes = 0;
        int cp = GetCodepointNext(s + i, &bytes);
        if (bytes <= 0) bytes = 1;
        int gi = GetGlyphIndex(*f, cp);
        if (cp != ' ' && cp != '\t') DrawTextCodepoint(*f, cp, (Vector2){ cx, cy }, (float)px, c);
        cx += f->glyphs[gi].advanceX ? (float)f->glyphs[gi].advanceX * scale : f->recs[gi].width * scale;
        i += bytes;
    }
}
static void txt(StFace face, int px, float x, float y, Color c, const char *s) { if (s) txtn(face, px, x, y, c, s, (int)strlen(s)); }
static void txtf(StFace face, int px, float x, float y, Color c, const char *fmt, ...) {
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    txt(face, px, x, y, c, buf);
}

/* text centred vertically in r; align 0 left, 1 centre, 2 right; an ellipsis
   when it does not fit */
static void txtIn(StFace face, int px, Rectangle r, Color c, int align, const char *s) {
    if (!s || !*s) return;
    int n = (int)strlen(s);
    int w = txtWn(face, px, s, n);
    char buf[600];
    if (w > r.width && r.width > 12) {
        int ell = txtW(face, px, "...");
        int k = n < (int)sizeof buf - 4 ? n : (int)sizeof buf - 4;
        while (k > 0) {
            k--;
            while (k > 0 && ((unsigned char)s[k] & 0xC0) == 0x80) k--;
            if (txtWn(face, px, s, k) + ell <= r.width) break;
        }
        memcpy(buf, s, (size_t)k);
        memcpy(buf + k, "...", 4);
        s = buf;
        n = k + 3;
        w = txtWn(face, px, s, n);
    }
    float x = align == 0 ? r.x : align == 1 ? r.x + (r.width - (float)w) / 2 : r.x + r.width - (float)w;
    float y = r.y + (r.height - (float)px) / 2 - (float)px * 0.08f;
    txtn(face, px, x, y, c, s, n);
}

/* word wrap: calls line(s, n) for each line no wider than w; returns lines */
static int txtWrap(StFace face, int px, float w, const char *s, int (*line)(const char *s, int n, void *u), void *u) {
    if (!s) return 0;
    int lines = 0;
    const char *p = s;
    while (*p) {
        const char *nl = strchr(p, '\n');
        int segEnd = nl ? (int)(nl - p) : (int)strlen(p);
        int start = 0;
        if (segEnd == 0) { if (line && !line(p, 0, u)) return lines + 1; lines++; }
        while (start < segEnd) {
            int best = -1, i = start;
            while (i <= segEnd) {
                if (i == segEnd || p[i] == ' ') {
                    if (txtWn(face, px, p + start, i - start) <= w) best = i;
                    else break;
                }
                i++;
            }
            if (best < 0) {   /* one word wider than the line: cut it by letters */
                int k = start + 1;
                while (k < segEnd && txtWn(face, px, p + start, k - start + 1) <= w) k++;
                while (k < segEnd && ((unsigned char)p[k] & 0xC0) == 0x80) k++;
                best = k;
            }
            lines++;
            if (line && !line(p + start, best - start, u)) return lines;
            start = best;
            while (start < segEnd && p[start] == ' ') start++;
        }
        if (!nl) break;
        p = nl + 1;
    }
    return lines;
}

typedef struct { StFace f; int px; float x, y, lh; Color c; int max, n; } StParaQ;
static int paraLine(const char *s, int n, void *u) {
    StParaQ *q = u;
    if (q->max > 0 && q->n >= q->max) return 0;
    txtn(q->f, q->px, q->x, q->y + (float)q->n * q->lh, q->c, s, n);
    q->n++;
    return 1;
}
/* a wrapped paragraph; returns its height */
static float txtPara(StFace f, int px, float x, float y, float w, float lh, Color c, const char *s, int maxLines) {
    StParaQ q = { f, px, x, y, lh, c, maxLines, 0 };
    txtWrap(f, px, w, s, paraLine, &q);
    return (float)q.n * lh;
}
static int paraCount(const char *s, int n, void *u) { (void)s; (void)n; int *k = u; (*k)++; return 1; }
static float txtParaH(StFace f, int px, float w, float lh, const char *s, int maxLines) {
    int k = 0;
    txtWrap(f, px, w, s, paraCount, &k);
    if (maxLines > 0 && k > maxLines) k = maxLines;
    return (float)k * lh;
}
