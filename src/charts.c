/* ========================================================================
   src/charts.c - the drawn figures: chromosome ideograms, the karyogram,
   histograms with a "you are here" marker, stacked bars, donuts and
   scatter plots. All take a rectangle and draw in it; none keep state
   except the karyogram's density cache.
   ======================================================================== */

static Color stainColor(const char *stain) {
    if (!strcmp(stain, "gneg")) return (Color){ 214, 218, 226, 255 };
    if (!strcmp(stain, "gpos25")) return (Color){ 172, 178, 190, 255 };
    if (!strcmp(stain, "gpos50")) return (Color){ 132, 138, 152, 255 };
    if (!strcmp(stain, "gpos75")) return (Color){ 96, 102, 116, 255 };
    if (!strcmp(stain, "gpos100")) return (Color){ 64, 68, 80, 255 };
    if (!strcmp(stain, "acen")) return (Color){ 196, 112, 112, 255 };
    if (!strcmp(stain, "gvar")) return (Color){ 120, 140, 170, 255 };
    if (!strcmp(stain, "stalk")) return (Color){ 150, 160, 180, 255 };
    return (Color){ 200, 204, 212, 255 };
}

/* a horizontal ideogram of one chromosome across r; returns the x of a position */
static float ideoX(Rectangle r, int chrom, long long pos) {
    double len = (double)g_chrom[chrom].len;
    return r.x + (float)((double)r.width * ((double)pos / (len > 0 ? len : 1)));
}
static void ideogramH(Rectangle r, int chrom) {
    double len = (double)g_chrom[chrom].len;
    if (len <= 0) return;
    for (int i = 0; i < g_bandN; i++) {
        StBand *b = &g_bands[i];
        if (b->chrom != chrom) continue;
        float x0 = r.x + (float)(r.width * (double)(b->start - 1) / len), x1 = r.x + (float)(r.width * (double)b->end / len);
        Color c = stainColor(b->stain);
        if (!strcmp(b->stain, "acen")) {
            bool pArm = b->name[0] == 'p';
            Vector2 a = { x0, r.y }, bb = { x0, r.y + r.height }, cc = { x1, r.y + r.height / 2 };
            if (!pArm) { a = (Vector2){ x1, r.y }; bb = (Vector2){ x0, r.y + r.height / 2 }; cc = (Vector2){ x1, r.y + r.height }; }
            if (pArm) DrawTriangle(a, bb, cc, c);
            else DrawTriangle(a, bb, cc, c);
            continue;
        }
        drRect((Rectangle){ x0, r.y, fmaxf(1, x1 - x0), r.height }, c);
    }
    drRoundLine(r, r.height / 2, drFade(TH->lineHi, 0.8f));
}
/* a vertical ideogram (top = pter) */
static void ideogramV(Rectangle r, int chrom, double maxLen) {
    if (maxLen <= 0) return;
    float h = (float)((double)r.height * (double)g_chrom[chrom].len / maxLen);
    Rectangle c = { r.x, r.y, r.width, h };
    for (int i = 0; i < g_bandN; i++) {
        StBand *b = &g_bands[i];
        if (b->chrom != chrom) continue;
        float y0 = c.y + (float)(h * (double)(b->start - 1) / (double)g_chrom[chrom].len);
        float y1 = c.y + (float)(h * (double)b->end / (double)g_chrom[chrom].len);
        Color col = stainColor(b->stain);
        if (!strcmp(b->stain, "acen")) {
            bool pArm = b->name[0] == 'p';
            if (pArm) DrawTriangle((Vector2){ c.x, y0 }, (Vector2){ c.x + c.width / 2, y1 }, (Vector2){ c.x + c.width, y0 }, col);
            else DrawTriangle((Vector2){ c.x + c.width / 2, y0 }, (Vector2){ c.x, y1 }, (Vector2){ c.x + c.width, y1 }, col);
            continue;
        }
        drRect((Rectangle){ c.x, y0, c.width, fmaxf(1, y1 - y0) }, col);
    }
}

/* ---- the karyogram: every chromosome, its variant density beside it ---- */
typedef struct { int n; float *dens; float *call; } KaryoChrom;
static KaryoChrom g_karyo[ST_NCHROM + 1];
static bool g_karyoLoaded;
static float g_karyoMax;
#define KARYO_MB 1000000
static void karyoLoad(void) {
    if (g_karyoLoaded) return;
    g_karyoLoaded = true;
    Rows *r = dbq("SELECT chrom, bin/100 AS mb, sum(nhet+nhom+nhemi), avg(callable) FROM covbin GROUP BY chrom, mb ORDER BY chrom, mb", NULL);
    for (int c = 1; c <= ST_NCHROM; c++) {
        int n = (int)(g_chrom[c].len / KARYO_MB) + 1;
        g_karyo[c].n = n;
        g_karyo[c].dens = xcalloc((size_t)n, sizeof(float));
        g_karyo[c].call = xcalloc((size_t)n, sizeof(float));
    }
    for (int i = 0; i < r->nrow; i++) {
        int c = (int)ri(r, i, 0), b = (int)ri(r, i, 1);
        if (c < 1 || c > ST_NCHROM || b < 0 || b >= g_karyo[c].n) continue;
        g_karyo[c].dens[b] = (float)rf(r, i, 2);
        g_karyo[c].call[b] = (float)rf(r, i, 3);
        /* a chip tests a few hundred positions per megabase: any ten count as read */
        if (dkChip()) g_karyo[c].call[b] = rf(r, i, 3) * 1e6 >= 10 ? 1.0f : 0.0f;
    }
    rowsFree(&r);
    float v[4000];
    int k = 0;
    for (int c = 1; c <= 22; c++) for (int b = 0; b < g_karyo[c].n && k < 4000; b++) if (g_karyo[c].dens[b] > 0) v[k++] = g_karyo[c].dens[b];
    /* scale by the 98th percentile so a few dense bins do not wash it out */
    for (int i = 1; i < k; i++) { float x = v[i]; int j = i - 1; while (j >= 0 && v[j] > x) { v[j + 1] = v[j]; j--; } v[j + 1] = x; }
    g_karyoMax = k ? v[(int)(k * 0.98f)] : 1;
}
static Color heatColor(float t) {
    t = clamp01(t);
    Color a = { 40, 62, 86, 255 }, b = { 110, 178, 206, 255 }, c = { 236, 214, 140, 255 };
    return t < 0.6f ? drMix(a, b, t / 0.6f) : drMix(b, c, (t - 0.6f) / 0.4f);
}
/* returns the chromosome clicked (0 none) */
static int karyogram(Rectangle r, bool showY) {
    karyoLoad();
    int n = showY ? 24 : 23;
    double maxLen = (double)g_chrom[1].len;
    float colW = r.width / (float)n;
    float chrW = fminf(14, colW * 0.34f);
    float top = r.y + 4, h = r.height - 30;
    int clicked = 0;
    for (int i = 0; i < n; i++) {
        int c = i + 1;
        float cx = r.x + colW * (float)i + colW / 2;
        Rectangle chr = { cx - chrW - 2, top, chrW, h };
        float ch = (float)(h * (double)g_chrom[c].len / maxLen);
        /* density strip */
        KaryoChrom *k = &g_karyo[c];
        float sx = cx + 2;
        for (int b = 0; b < k->n; b++) {
            float y0 = top + (float)(h * (double)b * KARYO_MB / maxLen);
            float y1 = top + (float)(h * (double)(b + 1) * KARYO_MB / maxLen);
            if (y1 > top + ch) y1 = top + ch;
            Color col = k->call[b] < 0.3f ? (Color){ 52, 56, 66, 255 } : heatColor(k->dens[b] / fmaxf(g_karyoMax, 1));
            drRect((Rectangle){ sx, y0, chrW * 0.8f, fmaxf(1, y1 - y0 + 0.5f) }, col);
        }
        ideogramV(chr, c, maxLen);
        Rectangle hit = { r.x + colW * (float)i, top - 2, colW, ch + 28 };
        bool hv;
        if (gxBehave(gxIdN(7000 + c), hit, &hv, NULL)) clicked = c;
        if (hv) {
            g_gx.cursor = MOUSE_CURSOR_POINTING_HAND;
            drRoundLine((Rectangle){ hit.x + 2, hit.y, hit.width - 4, hit.height }, 6, TH->accent);
            float my = clampf(g_gx.mouse.y, top, top + ch);
            long long pos = (long long)((my - top) / h * maxLen);
            const char *band = "";
            for (int j = 0; j < g_bandN; j++) if (g_bands[j].chrom == c && pos >= g_bands[j].start && pos <= g_bands[j].end) band = g_bands[j].name;
            gxTip(gxIdN(7100 + c), hit, sstr("Chromosome %s, band %s%s (%s)\nClick to open it in the genome browser.", chromName(c), chromName(c), band, fmtBp((double)pos)));
        }
        txtIn(F_SEMI, TH->fsSmall, (Rectangle){ hit.x, top + ch + 6, hit.width, 18 }, hv ? TH->text : TH->dim, 1, chromName(c));
    }
    return clicked;
}

/* ---- a distribution with the sample marked ---------------------------------
   bins: counts; lo/hi: the x range; you: the sample's value (NAN = none) */
static void histogram(Rectangle r, const float *bins, int n, float lo, float hi, float you, Color col, const char *xlabel) {
    float mx = 0;
    for (int i = 0; i < n; i++) mx = fmaxf(mx, bins[i]);
    if (mx <= 0) mx = 1;
    Rectangle a = { r.x, r.y, r.width, r.height - 22 };
    float bw = a.width / (float)n;
    for (int i = 0; i < n; i++) {
        float h = a.height * bins[i] / mx;
        float bx = a.x + bw * (float)i;
        bool mine = !isnan(you) && you >= lo + (hi - lo) * (float)i / (float)n && you < lo + (hi - lo) * (float)(i + 1) / (float)n;
        drRect((Rectangle){ bx + 0.5f, a.y + a.height - h, bw - 1, h }, mine ? col : drFade(col, 0.35f));
    }
    drHLine(a.x, a.y + a.height, a.width, TH->line);
    if (!isnan(you)) {
        float x = a.x + a.width * clamp01((you - lo) / (hi - lo));
        drRect((Rectangle){ x - 1, a.y - 4, 2, a.height + 6 }, TH->text);
        DrawTriangle((Vector2){ x - 5, a.y - 10 }, (Vector2){ x, a.y - 3 }, (Vector2){ x + 5, a.y - 10 }, TH->text);
    }
    if (xlabel) txtIn(F_REG, TH->fsTiny, (Rectangle){ a.x, a.y + a.height + 4, a.width, 16 }, TH->faint, 1, xlabel);
}

/* a horizontal stacked bar of fractions (sum <= 1) */
static void stackBar(Rectangle r, const float *v, const Color *cols, int n) {
    drRound(r, r.height / 2, TH->field);
    float x = r.x;
    gxClipPush(r);
    for (int i = 0; i < n; i++) {
        float w = r.width * v[i];
        if (w <= 0) continue;
        drRect((Rectangle){ x, r.y, w + 0.5f, r.height }, cols[i]);
        x += w;
    }
    gxClipPop();
}

/* a donut; values need not sum to one */
static void donut(Vector2 c, float r0, float r1, const float *v, const Color *cols, int n) {
    float tot = 0;
    for (int i = 0; i < n; i++) tot += v[i];
    if (tot <= 0) { DrawRing(c, r0, r1, 0, 360, 64, TH->field); return; }
    float a = -90;
    for (int i = 0; i < n; i++) {
        float s = 360 * v[i] / tot;
        if (s > 0.2f) DrawRing(c, r0, r1, a, a + s - (n > 1 ? 0.8f : 0), 64, cols[i]);
        a += s;
    }
}

/* a horizontal value bar with a label and number on one row */
static void labelBar(Rectangle r, const char *label, float v, float vmax, Color c, const char *num) {
    float lw = fminf(220, r.width * 0.38f);
    txtIn(F_REG, TH->fsSmall, (Rectangle){ r.x, r.y, lw - 8, r.height }, TH->text, 0, label);
    Rectangle b = { r.x + lw, r.y + r.height / 2 - 5, r.width - lw - 64, 10 };
    gxMeter(b, vmax > 0 ? v / vmax : 0, c);
    txtIn(F_SEMI, TH->fsSmall, (Rectangle){ b.x + b.width + 8, r.y, 56, r.height }, TH->text, 2, num);
}
