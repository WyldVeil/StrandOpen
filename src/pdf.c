/* ========================================================================
   src/pdf.c - a small PDF writer for the report generator.

   A4 pages, the standard Helvetica, Helvetica-Bold and Courier fonts (every
   PDF reader has them, so nothing is embedded), headings, wrapped
   paragraphs, simple tables with wrapping cells, a footer with page
   numbers. Text is WinAnsi: UTF-8 input is mapped to Latin-1 and anything
   beyond it becomes '?'. Pages are kept in memory and the file is written
   in one go at the end, with its cross-reference table.
   ======================================================================== */

#define PDF_W 595.0f
#define PDF_H 842.0f
#define PDF_M 48.0f
#define PDF_BOTTOM 70.0f      /* above the footer */

typedef struct { char *s; size_t n, cap; } PdfBuf;
typedef struct {
    PdfBuf *pages;
    int npages, cap;
    float y;                  /* the next line's top, in PDF units from the bottom */
    char footer[240];
    const char *const *thead; /* the open table's header, repeated when a row starts a new page */
    const float *tw;
    int tn;
} Pdf;

/* Helvetica and Helvetica-Bold advance widths for ASCII 32-126, in 1/1000 em (Adobe AFM) */
static const unsigned short g_pdfHelv[95] = {
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556, 556, 556, 556, 556, 556, 556, 556,
    278, 278, 584, 584, 584, 556, 1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667, 778, 722, 667,
    611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556, 333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833,
    556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584 };
static const unsigned short g_pdfHelvB[95] = {
    278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556, 556, 556, 556, 556, 556, 556, 556,
    333, 333, 584, 584, 584, 611, 975, 722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611, 833, 722, 778, 667, 778, 722, 667,
    611, 722, 667, 944, 667, 667, 611, 333, 278, 333, 584, 556, 333, 556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889,
    611, 611, 611, 611, 389, 556, 333, 611, 556, 778, 556, 556, 500, 389, 280, 389, 584 };

enum { PDF_REG = 1, PDF_BOLD = 2, PDF_MONO = 3 };

static void pbPut(PdfBuf *b, const char *s, size_t n) {
    if (b->n + n + 1 > b->cap) {
        b->cap = (b->n + n + 1) * 2 + 4096;
        b->s = xrealloc(b->s, b->cap);
    }
    memcpy(b->s + b->n, s, n);
    b->n += n;
    b->s[b->n] = 0;
}
static void pbf(PdfBuf *b, const char *fmt, ...) __attribute__((format(printf, 2, 3)));
static void pbf(PdfBuf *b, const char *fmt, ...) {
    char tmp[4096];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(tmp, sizeof tmp, fmt, ap);
    va_end(ap);
    if (n > 0) pbPut(b, tmp, (size_t)n < sizeof tmp ? (size_t)n : sizeof tmp - 1);
}

/* UTF-8 to Latin-1, one byte per character ('?' beyond it) */
static size_t pdfLatin1(const char *s, unsigned char *out, size_t cap) {
    size_t n = 0;
    const unsigned char *p = (const unsigned char *)s;
    while (*p && n + 1 < cap) {
        unsigned cp = *p;
        if (cp < 0x80) p++;
        else if ((cp & 0xE0) == 0xC0 && p[1]) { cp = ((cp & 0x1F) << 6) | (p[1] & 0x3F); p += 2; }
        else if ((cp & 0xF0) == 0xE0 && p[1] && p[2]) { cp = ((cp & 0x0F) << 12) | ((p[1] & 0x3F) << 6) | (p[2] & 0x3F); p += 3; }
        else if ((cp & 0xF8) == 0xF0 && p[1] && p[2] && p[3]) { cp = 0x10000; p += 4; }
        else { p++; cp = '?'; }
        if (cp == 0x2013 || cp == 0x2014) cp = '-';
        else if (cp == 0x2018 || cp == 0x2019) cp = '\'';
        else if (cp == 0x201C || cp == 0x201D) cp = '"';
        else if (cp == 0x2026) cp = '.';
        out[n++] = (unsigned char)(cp < 256 ? cp : '?');
    }
    out[n] = 0;
    return n;
}

static float pdfCharW(int font, unsigned char c) {
    if (font == PDF_MONO) return 600;
    if (c >= 32 && c <= 126) return font == PDF_BOLD ? g_pdfHelvB[c - 32] : g_pdfHelv[c - 32];
    return 556;
}
static float pdfTextW(int font, float size, const unsigned char *s, size_t n) {
    float w = 0;
    for (size_t i = 0; i < n; i++) w += pdfCharW(font, s[i]);
    return w * size / 1000.0f;
}

/* the colours: light pages by default, dark pages (light text on near-black) on request */
static Color PDF_INK = { 30, 34, 42, 255 }, PDF_DIM = { 96, 104, 118, 255 }, PDF_HEAD = { 22, 88, 120, 255 };
static bool g_pdfDark;
static const char *g_pdfRule = "0.80 0.84 0.88", *g_pdfBand = "0.92 0.94 0.96", *g_pdfLine = "0.88 0.90 0.93";
static void pdfTheme(bool dark) {
    g_pdfDark = dark;
    if (dark) {
        PDF_INK = (Color){ 226, 230, 236, 255 }; PDF_DIM = (Color){ 150, 158, 172, 255 }; PDF_HEAD = (Color){ 120, 190, 226, 255 };
        g_pdfRule = "0.30 0.34 0.40"; g_pdfBand = "0.16 0.18 0.22"; g_pdfLine = "0.22 0.25 0.30";
    } else {
        PDF_INK = (Color){ 30, 34, 42, 255 }; PDF_DIM = (Color){ 96, 104, 118, 255 }; PDF_HEAD = (Color){ 22, 88, 120, 255 };
        g_pdfRule = "0.80 0.84 0.88"; g_pdfBand = "0.92 0.94 0.96"; g_pdfLine = "0.88 0.90 0.93";
    }
}

static void pdfNewPage(Pdf *p) {
    if (p->npages == p->cap) {
        p->cap = p->cap ? p->cap * 2 : 8;
        p->pages = xrealloc(p->pages, sizeof *p->pages * (size_t)p->cap);
    }
    p->pages[p->npages++] = (PdfBuf){ 0 };
    p->y = PDF_H - PDF_M;
    if (g_pdfDark) pbf(&p->pages[p->npages - 1], "0.086 0.094 0.114 rg 0 0 %.0f %.0f re f\n", PDF_W, PDF_H);
}
static PdfBuf *pdfCur(Pdf *p) { if (!p->npages) pdfNewPage(p); return &p->pages[p->npages - 1]; }
static void pdfEnsure(Pdf *p, float h) { if (!p->npages || p->y - h < PDF_BOTTOM) pdfNewPage(p); }

/* one run of text at a baseline */
static void pdfRun(Pdf *p, int font, float size, float x, float y, Color c, const unsigned char *s, size_t n) {
    PdfBuf *b = pdfCur(p);
    pbf(b, "%.3f %.3f %.3f rg BT /F%d %.1f Tf %.2f %.2f Td (", c.r / 255.0f, c.g / 255.0f, c.b / 255.0f, font, size, x, y);
    for (size_t i = 0; i < n; i++) {
        unsigned char ch = s[i];
        if (ch == '(' || ch == ')' || ch == '\\') { pbPut(b, "\\", 1); pbPut(b, (const char *)&ch, 1); }
        else if (ch < 32) pbPut(b, " ", 1);
        else if (ch > 126) pbf(b, "\\%03o", ch);
        else pbPut(b, (const char *)&ch, 1);
    }
    pbPut(b, ") Tj ET\n", 8);
}

/* wrap text to a width: the start and length of each line (at most max lines) */
static int pdfWrap(int font, float size, const unsigned char *s, size_t n, float width, size_t *st, size_t *ln, int max) {
    int k = 0;
    size_t i = 0;
    while (i < n && k < max) {
        while (i < n && s[i] == ' ') i++;
        if (i >= n) break;
        size_t lineStart = i, lastBreak = 0;
        float w = 0;
        size_t j = i;
        for (; j < n; j++) {
            if (s[j] == '\n') break;
            float cw = pdfCharW(font, s[j]) * size / 1000.0f;
            if (w + cw > width && j > lineStart) break;
            w += cw;
            if (s[j] == ' ') lastBreak = j;
        }
        size_t end = j;
        if (j < n && s[j] != '\n' && lastBreak > lineStart) end = lastBreak;   /* break at the last space */
        st[k] = lineStart;
        ln[k] = end - lineStart;
        k++;
        i = (j < n && s[j] == '\n' && end == j) ? j + 1 : end;
    }
    return k;
}

/* a wrapped paragraph at x with width; returns its height (drawn when draw) */
static float pdfParaAt(Pdf *p, int font, float size, float x, float y, float width, Color c, const char *text, bool draw) {
    static unsigned char buf[16384];
    size_t n = pdfLatin1(text ? text : "", buf, sizeof buf);
    size_t st[200], ln[200];
    int k = pdfWrap(font, size, buf, n, width, st, ln, 200);
    float lh = size * 1.32f;
    if (draw) for (int i = 0; i < k; i++) pdfRun(p, font, size, x, y - size - lh * (float)i, c, buf + st[i], ln[i]);
    return lh * (float)(k ? k : 1);
}


static void pdfPara(Pdf *p, const char *text, float size, Color c, int font) {
    float w = PDF_W - 2 * PDF_M;
    float h = pdfParaAt(p, font, size, PDF_M, 0, w, c, text, false);
    pdfEnsure(p, h);
    pdfParaAt(p, font, size, PDF_M, p->y, w, c, text, true);
    p->y -= h + size * 0.45f;
}
static void pdfHeading(Pdf *p, const char *text, int level) {
    float size = level == 1 ? 15 : 11.5f;
    pdfEnsure(p, size * 4);             /* never a heading alone at the foot of a page */
    p->y -= level == 1 ? 10 : 6;
    float h = pdfParaAt(p, PDF_BOLD, size, PDF_M, p->y, PDF_W - 2 * PDF_M, level == 1 ? PDF_HEAD : PDF_INK, text, true);
    p->y -= h + 4;
    if (level == 1) {
        pbf(pdfCur(p), "%s RG 0.6 w %.2f %.2f m %.2f %.2f l S\n", g_pdfRule, PDF_M, p->y + 2, PDF_W - PDF_M, p->y + 2);
        p->y -= 6;
    }
}

/* tables: cells wrapped to their column widths, a shaded header band that stays with the first row and is
   repeated at the top of each page the table runs onto */
static float pdfRowH(Pdf *p, const char *const *cells, const float *widths, int ncol, bool header, int monoCol) {
    float size = header ? 8.5f : 8.8f, pad = 4, h = 0;
    for (int i = 0; i < ncol; i++) {
        int font = header ? PDF_BOLD : (i == monoCol ? PDF_MONO : PDF_REG);
        float hh = pdfParaAt(p, font, size, 0, 0, widths[i] - 2 * pad, PDF_INK, cells[i] ? cells[i] : "", false);
        if (hh > h) h = hh;
    }
    return h + 2 * pad;
}
static void pdfRowDraw(Pdf *p, const char *const *cells, const float *widths, int ncol, bool header, int monoCol, float h) {
    float size = header ? 8.5f : 8.8f, pad = 4;
    PdfBuf *b = pdfCur(p);
    if (header) pbf(b, "%s rg %.2f %.2f %.2f %.2f re f\n", g_pdfBand, PDF_M, p->y - h, PDF_W - 2 * PDF_M, h);
    else pbf(b, "%s RG 0.4 w %.2f %.2f m %.2f %.2f l S\n", g_pdfLine, PDF_M, p->y - h, PDF_W - PDF_M, p->y - h);
    float x = PDF_M;
    for (int i = 0; i < ncol; i++) {
        int font = header ? PDF_BOLD : (i == monoCol ? PDF_MONO : PDF_REG);
        pdfParaAt(p, font, size, x + pad, p->y - pad, widths[i] - 2 * pad, PDF_INK, cells[i] ? cells[i] : "", true);
        x += widths[i];
    }
    p->y -= h;
}
/* the header, kept on the same page as the first row (first may be NULL) */
static void pdfTableBegin(Pdf *p, const char *const *head, const float *widths, int ncol, const char *const *first, int monoCol) {
    float hh = pdfRowH(p, head, widths, ncol, true, -1), fh = first ? pdfRowH(p, first, widths, ncol, false, monoCol) : 0;
    pdfEnsure(p, hh + fh);
    pdfRowDraw(p, head, widths, ncol, true, -1, hh);
    p->thead = head;
    p->tw = widths;
    p->tn = ncol;
}
static void pdfRow(Pdf *p, const char *const *cells, const float *widths, int ncol, int monoCol) {
    float h = pdfRowH(p, cells, widths, ncol, false, monoCol);
    if (p->npages && p->y - h >= PDF_BOTTOM) { pdfRowDraw(p, cells, widths, ncol, false, monoCol, h); return; }
    pdfNewPage(p);
    if (p->thead) pdfRowDraw(p, p->thead, p->tw, p->tn, true, -1, pdfRowH(p, p->thead, p->tw, p->tn, true, -1));
    pdfRowDraw(p, cells, widths, ncol, false, monoCol, h);
}
static void pdfTableEnd(Pdf *p) { p->thead = NULL; p->tw = NULL; p->tn = 0; }

/* the whole file: catalog, pages, three standard fonts, info, xref */
static bool pdfSave(Pdf *p, const char *path, const char *title) {
    FILE *f = fopen(path, "wb");
    if (!f) return false;
    int n = p->npages;
    long *off = xcalloc((size_t)(6 + 2 * n + 1), sizeof(long));
    fputs("%PDF-1.4\n%\xe2\xe3\xcf\xd3\n", f);
    int obj = 0;
#define PDF_OBJ() (off[++obj] = ftell(f), fprintf(f, "%d 0 obj\n", obj))
    PDF_OBJ(); fputs("<< /Type /Catalog /Pages 2 0 R >>\nendobj\n", f);
    PDF_OBJ(); fputs("<< /Type /Pages /Kids [", f);
    for (int i = 0; i < n; i++) fprintf(f, " %d 0 R", 7 + 2 * i);
    fprintf(f, " ] /Count %d >>\nendobj\n", n);
    PDF_OBJ(); fputs("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>\nendobj\n", f);
    PDF_OBJ(); fputs("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>\nendobj\n", f);
    PDF_OBJ(); fputs("<< /Type /Font /Subtype /Type1 /BaseFont /Courier /Encoding /WinAnsiEncoding >>\nendobj\n", f);
    PDF_OBJ();
    {
        unsigned char t[400];
        pdfLatin1(title, t, sizeof t);
        fputs("<< /Producer (Strand Open) /Title (", f);
        for (unsigned char *c = t; *c; c++) { if (*c == '(' || *c == ')' || *c == '\\') fputc('\\', f); fputc(*c, f); }
        fputs(") >>\nendobj\n", f);
    }
    for (int i = 0; i < n; i++) {
        /* the footer, now that the page count is known */
        PdfBuf *b = &p->pages[i];
        unsigned char ft[300];
        size_t fl = pdfLatin1(sstr("%s   page %d of %d", p->footer, i + 1, n), ft, sizeof ft);
        pbf(b, "%s RG 0.5 w %.2f %.2f m %.2f %.2f l S\n", g_pdfRule, PDF_M, 46.0f, PDF_W - PDF_M, 46.0f);
        {
            int save = p->npages;
            p->npages = i + 1;
            pdfRun(p, PDF_REG, 7.5f, PDF_M, 34, PDF_DIM, ft, fl);
            p->npages = save;
        }
        PDF_OBJ();
        fprintf(f, "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.0f %.0f] /Resources << /Font << /F1 3 0 R /F2 4 0 R /F3 5 0 R >> >> "
                   "/Contents %d 0 R >>\nendobj\n", PDF_W, PDF_H, obj + 1);
        PDF_OBJ();
        fprintf(f, "<< /Length %zu >>\nstream\n", b->n);
        if (b->n) fwrite(b->s, 1, b->n, f);
        fputs("\nendstream\nendobj\n", f);
    }
#undef PDF_OBJ
    long xref = ftell(f);
    fprintf(f, "xref\n0 %d\n0000000000 65535 f \n", obj + 1);
    for (int i = 1; i <= obj; i++) fprintf(f, "%010ld 00000 n \n", off[i]);
    fprintf(f, "trailer\n<< /Size %d /Root 1 0 R /Info 6 0 R >>\nstartxref\n%ld\n%%%%EOF\n", obj + 1, xref);
    free(off);
    bool ok = !ferror(f);
    if (fclose(f) != 0) ok = false;
    return ok;
}

static void pdfFree(Pdf *p) {
    for (int i = 0; i < p->npages; i++) free(p->pages[i].s);
    free(p->pages);
    *p = (Pdf){ 0 };
}
