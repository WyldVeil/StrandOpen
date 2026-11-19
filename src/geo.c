/* ========================================================================
   src/geo.c - maps.

   Land outlines from Natural Earth 1:50m (public domain; assets/geo/
   land50.bin: int32 rings, int32 points, then per ring int32 n and n
   float lon/lat pairs). A view is an equirectangular box fitted to a
   rectangle with the right shape at its middle latitude. The land is
   filled by scanlines into a texture once per view and size, so drawing a
   map each frame costs one textured quad plus the points.
   ======================================================================== */

static const unsigned char g_landBin[] = {
#embed "../assets/geo/land50.bin"
};

typedef struct { float lon0, lon1, lat0, lat1; Rectangle r; } GeoView;

static struct {
    int n;
    const float *pts[2048];
    int len[2048];
    float bx0[2048], bx1[2048], by0[2048], by1[2048];
} g_land;

static int32_t geoI32(size_t off) {
    int32_t v;
    memcpy(&v, g_landBin + off, 4);
    return v;
}
static void geoLoad(void) {
    if (g_land.n) return;
    int rings = geoI32(0);
    size_t off = 8;
    for (int i = 0; i < rings && g_land.n < 2048 && off + 4 <= sizeof g_landBin; i++) {
        int n = geoI32(off);
        off += 4;
        if (n <= 0 || off + (size_t)n * 8 > sizeof g_landBin) break;
        const float *p = (const float *)(g_landBin + off);
        float x0 = 1e9f, x1 = -1e9f, y0 = 1e9f, y1 = -1e9f;
        for (int k = 0; k < n; k++) {
            x0 = fminf(x0, p[2 * k]); x1 = fmaxf(x1, p[2 * k]);
            y0 = fminf(y0, p[2 * k + 1]); y1 = fmaxf(y1, p[2 * k + 1]);
        }
        int j = g_land.n++;
        g_land.pts[j] = p; g_land.len[j] = n;
        g_land.bx0[j] = x0; g_land.bx1[j] = x1; g_land.by0[j] = y0; g_land.by1[j] = y1;
        off += (size_t)n * 8;
    }
}

/* a view of the box lon0..lon1 x lat0..lat1 (plus pad degrees) widened to the
   rectangle's shape */
static GeoView geoFit(Rectangle r, float lon0, float lon1, float lat0, float lat1, float pad) {
    lon0 -= pad; lon1 += pad; lat0 -= pad; lat1 += pad;
    lat0 = fmaxf(lat0, -80); lat1 = fminf(lat1, 84);
    float k = cosf((lat0 + lat1) * 0.5f * DEG2RAD);
    if (k < 0.2f) k = 0.2f;
    float w = (lon1 - lon0) * k, h = lat1 - lat0, ar = r.width / fmaxf(1, r.height);
    if (w / fmaxf(h, 1e-3f) < ar) { float nw = h * ar / k, c = (lon0 + lon1) * 0.5f; lon0 = c - nw / 2; lon1 = c + nw / 2; }
    else { float nh = w / ar, c = (lat0 + lat1) * 0.5f; lat0 = c - nh / 2; lat1 = c + nh / 2; }
    return (GeoView){ lon0, lon1, lat0, lat1, r };
}
static Vector2 geoXY(const GeoView *v, float lon, float lat) {
    return (Vector2){ v->r.x + (lon - v->lon0) / (v->lon1 - v->lon0) * v->r.width,
                      v->r.y + (v->lat1 - lat) / (v->lat1 - v->lat0) * v->r.height };
}

static struct { Texture2D tex; float lon0, lon1, lat0, lat1; int w, h; Color col; bool ok; } g_geoTex;
static int geoCmpF(const void *a, const void *b) { float x = *(const float *)a, y = *(const float *)b; return (x > y) - (x < y); }

/* the land, filled, for this view */
static void geoLand(const GeoView *v, Color land) {
    geoLoad();
    int w = (int)(v->r.width * g_gx.scale), h = (int)(v->r.height * g_gx.scale);
    if (w < 8 || h < 8) return;
    if (!g_geoTex.ok || g_geoTex.w != w || g_geoTex.h != h || g_geoTex.lon0 != v->lon0 || g_geoTex.lon1 != v->lon1 ||
        g_geoTex.lat0 != v->lat0 || g_geoTex.lat1 != v->lat1 || memcmp(&g_geoTex.col, &land, sizeof land)) {
        if (g_geoTex.ok) UnloadTexture(g_geoTex.tex);
        Image img = GenImageColor(w, h, BLANK);
        Color *px = img.data;
        float xs[512];
        for (int y = 0; y < h; y++) {
            float lat = v->lat1 - ((float)y + 0.5f) / (float)h * (v->lat1 - v->lat0);
            for (int i = 0; i < g_land.n; i++) {
                if (lat < g_land.by0[i] || lat > g_land.by1[i] || g_land.bx1[i] < v->lon0 || g_land.bx0[i] > v->lon1) continue;
                const float *p = g_land.pts[i];
                int n = g_land.len[i], m = 0;
                for (int k = 0; k < n && m < 512; k++) {
                    const float *a = p + 2 * k, *b = p + 2 * ((k + 1) % n);
                    if ((a[1] > lat) != (b[1] > lat)) xs[m++] = a[0] + (lat - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
                }
                if (m < 2) continue;
                qsort(xs, (size_t)m, sizeof xs[0], geoCmpF);
                for (int k = 0; k + 1 < m; k += 2) {
                    int x0 = (int)ceilf((xs[k] - v->lon0) / (v->lon1 - v->lon0) * (float)w - 0.5f);
                    int x1 = (int)floorf((xs[k + 1] - v->lon0) / (v->lon1 - v->lon0) * (float)w - 0.5f);
                    if (x0 < 0) x0 = 0;
                    if (x1 >= w) x1 = w - 1;
                    for (int x = x0; x <= x1; x++) px[(size_t)y * (size_t)w + (size_t)x] = land;
                }
            }
        }
        g_geoTex.tex = LoadTextureFromImage(img);
        SetTextureFilter(g_geoTex.tex, TEXTURE_FILTER_BILINEAR);
        UnloadImage(img);
        g_geoTex.ok = true; g_geoTex.w = w; g_geoTex.h = h; g_geoTex.col = land;
        g_geoTex.lon0 = v->lon0; g_geoTex.lon1 = v->lon1; g_geoTex.lat0 = v->lat0; g_geoTex.lat1 = v->lat1;
    }
    DrawTexturePro(g_geoTex.tex, (Rectangle){ 0, 0, (float)w, (float)h }, v->r, (Vector2){ 0, 0 }, 0, WHITE);
}
