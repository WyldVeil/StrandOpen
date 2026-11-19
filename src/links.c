/* ========================================================================
   src/links.c - the link menu.

   Every box that stands for a variant, gene, score, medicine gene or
   haplogroup opens the same small menu when clicked: where to see it in
   Strand first, then the public sites that know about it (SNPedia, dbSNP,
   ClinVar, gnomAD, GeneCards, OMIM...), opened in the default browser.

       lkBegin("FLG p.Ser761fs");
       lkVariant(variantId);        // Strand views + variant sites
       lkGene("FLG");               // gene sites
       lkOpen();                    // the menu, at the mouse

   gxBoxClick() makes a whole box clickable; call it after the box's own
   buttons and links so they keep their clicks.
   ======================================================================== */

typedef enum { LK_HEAD, LK_URL, LK_DETAILS, LK_READS, LK_GENE, LK_ENTER, LK_LOCUS, LK_FN } LkKind;
typedef struct { LkKind kind; char label[64]; char sub[48]; char arg[400]; long long id; void (*fn)(long long); } LkItem;
/* fullH: the menu's height with every row shown; scroll: how far a menu
   taller than the window has been scrolled */
static struct { LkItem it[40]; int n; char title[120]; bool hasWeb; float fullH, scroll; } g_lk;
static const GxId LK_POPUP = 0x51a7d00dull;
static void lkOpen(void);

static void lkBegin(const char *title) {
    g_lk.n = 0;
    g_lk.hasWeb = false;
    snprintf(g_lk.title, sizeof g_lk.title, "%s", title ? title : "");
}
static void lkAdd(LkKind k, const char *label, const char *sub, const char *arg, long long id) {
    if (g_lk.n >= (int)NELEMS(g_lk.it)) return;
    LkItem *t = &g_lk.it[g_lk.n++];
    t->kind = k;
    t->id = id;
    snprintf(t->label, sizeof t->label, "%s", label);
    snprintf(t->sub, sizeof t->sub, "%s", sub ? sub : "");
    snprintf(t->arg, sizeof t->arg, "%s", arg ? arg : "");
    t->fn = NULL;
}
static void lkHead(const char *s) { lkAdd(LK_HEAD, s, NULL, NULL, 0); }
static void lkUrl(const char *label, const char *site, const char *url) {
    if (!g_lk.hasWeb) { lkHead("On the web"); g_lk.hasWeb = true; }
    lkAdd(LK_URL, label, site, url, 0);
}
static void lkEnter(const char *label, const char *target) { lkAdd(LK_ENTER, label, "Strand Open", target, 0); }
/* an action inside Strand */
static void lkFn(const char *label, void (*fn)(long long), long long arg) {
    lkAdd(LK_FN, label, "Strand Open", NULL, arg);
    g_lk.it[g_lk.n - 1].fn = fn;
}

/* the sites that know a dbSNP id */
static void lkRsid(long long rs) {
    if (rs <= 0) return;
    lkUrl(sstr("rs%lld on SNPedia", rs), "snpedia.com", sstr("https://www.snpedia.com/index.php/Rs%lld", rs));
    lkUrl(sstr("rs%lld in dbSNP", rs), "ncbi.nlm.nih.gov", sstr("https://www.ncbi.nlm.nih.gov/snp/rs%lld", rs));
    lkUrl("Studies mentioning it", "pubmed.ncbi.nlm.nih.gov", sstr("https://pubmed.ncbi.nlm.nih.gov/?term=rs%lld", rs));
    lkUrl("GWAS Catalog associations", "ebi.ac.uk/gwas", sstr("https://www.ebi.ac.uk/gwas/variants/rs%lld", rs));
    lkUrl("Ensembl variant page", "ensembl.org", sstr("https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=rs%lld", rs));
}
static void lkLocus(int chrom, long long a, long long b) {
    if (chrom <= 0) return;
    lkUrl("UCSC Genome Browser", "genome.ucsc.edu", sstr("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr%s:%lld-%lld", chromName(chrom), a, b));
}

/* a variant of the person: Strand views, then its public records */
static void lkVariant(long long vid) {
    Rows *v = dbq("SELECT chrom, pos, ref, alt, rsid, clin_id, gene, hgvs_p FROM variant WHERE id=?", sstr("%lld", vid), NULL);
    if (!v->nrow) { rowsFree(&v); return; }
    int chrom = (int)ri(v, 0, 0);
    long long pos = ri(v, 0, 1), rs = ri(v, 0, 4), clin = ri(v, 0, 5);
    lkHead("In Strand Open");
    lkAdd(LK_DETAILS, "Variant details", "Strand Open", rc(v, 0, 6), vid);
    lkAdd(LK_READS, dkReads() ? "See the reads" : "Show in the genome browser", "genome browser", NULL, vid);
    if (*rc(v, 0, 6)) lkAdd(LK_GENE, sstr("%s gene page", rc(v, 0, 6)), "Strand Open", rc(v, 0, 6), 0);
    lkRsid(rs);
    if (clin > 0) lkUrl("ClinVar record", "ncbi.nlm.nih.gov/clinvar", sstr("https://www.ncbi.nlm.nih.gov/clinvar/variation/%lld/", clin));
    lkUrl("gnomAD frequencies (800,000 people)", "gnomad.broadinstitute.org",
          sstr("https://gnomad.broadinstitute.org/variant/%s-%lld-%s-%s?dataset=gnomad_r4", chromName(chrom), pos, rc(v, 0, 2), rc(v, 0, 3)));
    lkLocus(chrom, pos - 50, pos + 50);
    rowsFree(&v);
}

/* a gene: its Strand page and the reference sites */
static void lkGene(const char *sym) {
    if (!sym || !*sym) return;
    Rows *g = dbq("SELECT id, ensg, hgnc, uniprot, (SELECT count(*) FROM gene_disease d WHERE d.gene=gene.id) FROM gene "
                  "WHERE symbol=? ORDER BY biotype='protein_coding' DESC LIMIT 1", sym, NULL);
    bool strandPage = false;
    for (int i = 0; i < g_lk.n; i++) if (g_lk.it[i].kind == LK_GENE) strandPage = true;
    if (!strandPage && g->nrow) {
        /* the gene page belongs with the Strand views, so when web links are
           already in the menu it goes in just before their heading */
        bool head = false;
        int web = g_lk.n;
        for (int i = 0; i < g_lk.n; i++) {
            if (g_lk.it[i].kind != LK_HEAD) continue;
            if (!strcmp(g_lk.it[i].label, "In Strand")) head = true;
            if (!strcmp(g_lk.it[i].label, "On the web") && web == g_lk.n) web = i;
        }
        int n0 = g_lk.n;
        if (n0 + (head ? 1 : 2) <= (int)NELEMS(g_lk.it)) {
            if (!head) lkHead("In Strand Open");
            lkAdd(LK_GENE, sstr("%s gene page", sym), "Strand Open", sym, 0);
        }
        int k = g_lk.n - n0;
        if (k > 0 && web < n0) {
            LkItem add[2];
            memcpy(add, &g_lk.it[n0], (size_t)k * sizeof add[0]);
            memmove(&g_lk.it[web + k], &g_lk.it[web], (size_t)(n0 - web) * sizeof add[0]);
            memcpy(&g_lk.it[web], add, (size_t)k * sizeof add[0]);
        }
    }
    lkUrl(sstr("%s on GeneCards", sym), "genecards.org", sstr("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", sym));
    lkUrl(sstr("%s in NCBI Gene", sym), "ncbi.nlm.nih.gov/gene", sstr("https://www.ncbi.nlm.nih.gov/gene/?term=%s%%5Bsym%%5D+AND+human%%5Borgn%%5D", sym));
    lkUrl("Conditions in OMIM", "omim.org", sstr("https://www.omim.org/search?search=%s", sym));
    if (g->nrow && *rc(g, 0, 3)) lkUrl("Protein in UniProt", "uniprot.org", sstr("https://www.uniprot.org/uniprotkb/%s/entry", rc(g, 0, 3)));
    if (g->nrow && ri(g, 0, 4) > 0) {
        char low[40];
        snprintf(low, sizeof low, "%s", sym);
        for (char *p = low; *p; p++) *p = (char)tolower((unsigned char)*p);
        lkUrl("Plain-language guide (MedlinePlus)", "medlineplus.gov", sstr("https://medlineplus.gov/genetics/gene/%s/", low));
    }
    lkUrl("Where it is active (GTEx)", "gtexportal.org", sstr("https://gtexportal.org/home/gene/%s", sym));
    if (g->nrow && *rc(g, 0, 1)) lkUrl("Variation in gnomAD", "gnomad.broadinstitute.org", sstr("https://gnomad.broadinstitute.org/gene/%s?dataset=gnomad_r4", rc(g, 0, 1)));
    lkUrl("Wikipedia", "wikipedia.org", sstr("https://en.wikipedia.org/w/index.php?search=%s+gene", sym));
    rowsFree(&g);
}
static void lkPubmed(const char *pmid, const char *label) {
    if (!pmid || !*pmid) return;
    char first[24] = "";
    sscanf(pmid, "%23[0-9]", first);
    if (first[0]) lkUrl(label ? label : "The study on PubMed", "pubmed.ncbi.nlm.nih.gov", sstr("https://pubmed.ncbi.nlm.nih.gov/%s/", first));
}
static void lkPgs(const char *pgsId, const char *pmid) {
    lkUrl(sstr("%s in the PGS Catalog", pgsId), "pgscatalog.org", sstr("https://www.pgscatalog.org/score/%s/", pgsId));
    lkPubmed(pmid, "The study behind the score");
}
/* a medicine gene: the guideline sites */
static void lkPgxGene(const char *gene) {
    lkUrl(sstr("%s guidelines (CPIC)", gene), "cpicpgx.org", sstr("https://cpicpgx.org/genes-drugs/?gene=%s", gene));
    lkUrl(sstr("%s on ClinPGx (PharmGKB)", gene), "clinpgx.org", sstr("https://www.clinpgx.org/search?query=%s", gene));
    if (!strncmp(gene, "CYP", 3) || !strcmp(gene, "NUDT15") || !strcmp(gene, "SLCO1B1") || !strcmp(gene, "DPYD") || !strcmp(gene, "TPMT") || !strcmp(gene, "UGT1A1"))
        lkUrl(sstr("%s star alleles (PharmVar)", gene), "pharmvar.org", sstr("https://www.pharmvar.org/gene/%s", gene));
}
static void lkDrug(const char *drug) {
    lkUrl(sstr("%s on ClinPGx (PharmGKB)", drug), "clinpgx.org", sstr("https://www.clinpgx.org/search?query=%s", drug));
    lkUrl(sstr("%s on MedlinePlus", drug), "medlineplus.gov", sstr("https://vsearch.nlm.nih.gov/vivisimo/cgi-bin/query-meta?v%%3Aproject=medlineplus&query=%s", drug));
}

/* a stretch of a chromosome (structural variants, deletions) */
static void lkRegion(int chrom, long long a, long long b) {
    lkHead("In Strand Open");
    lkAdd(LK_LOCUS, "Show in the genome browser", "Strand Open", sstr("%lld %lld", a, b), chrom);
    lkLocus(chrom, a, b);
    lkUrl("Structural variants in gnomAD", "gnomad.broadinstitute.org",
          sstr("https://gnomad.broadinstitute.org/region/%s-%lld-%lld?dataset=gnomad_sv_r4", chromName(chrom), a, b));
    lkUrl("DECIPHER (patients with changes here)", "deciphergenomics.org",
          sstr("https://www.deciphergenomics.org/browser#q/grch38:%s:%lld-%lld", chromName(chrom), a, b));
}
/* a trait or well-known variant card: its rsIDs (from the genotype text) and gene (in the title's brackets) */
static void lkTrait(const char *title, const char *genotypes) {
    char gene[40] = "";
    const char *lp = strrchr(title, '(');
    if (lp) sscanf(lp + 1, "%39[A-Za-z0-9-]", gene);
    if (gene[0]) {
        Rows *g = dbq("SELECT 1 FROM gene WHERE symbol=?", gene, NULL);
        if (!g->nrow) gene[0] = 0;
        rowsFree(&g);
    }
    if (gene[0]) { lkHead("In Strand Open"); lkAdd(LK_GENE, sstr("%s gene page", gene), "Strand Open", gene, 0); }
    int n = 0;
    for (const char *p = genotypes; p && (p = strstr(p, "rs")) && n < 4; p += 2) {
        long long rs = atoll(p + 2);
        if (rs <= 0) continue;
        lkUrl(sstr("rs%lld on SNPedia", rs), "snpedia.com", sstr("https://www.snpedia.com/index.php/Rs%lld", rs));
        lkUrl(sstr("rs%lld in dbSNP", rs), "ncbi.nlm.nih.gov", sstr("https://www.ncbi.nlm.nih.gov/snp/rs%lld", rs));
        n++;
    }
    if (gene[0]) {
        lkUrl(sstr("%s on GeneCards", gene), "genecards.org", sstr("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", gene));
        lkUrl(sstr("%s in NCBI Gene", gene), "ncbi.nlm.nih.gov/gene", sstr("https://www.ncbi.nlm.nih.gov/gene/?term=%s%%5Bsym%%5D+AND+human%%5Borgn%%5D", gene));
    }
}
/* the usual menus, built and opened in one call */
static void lkVariantMenu(long long vid) {
    Rows *v = dbq("SELECT gene, hgvs_p, rsid, chrom, pos FROM variant WHERE id=?", sstr("%lld", vid), NULL);
    if (!v->nrow) { rowsFree(&v); return; }
    const char *gene = rc(v, 0, 0);
    lkBegin(*rc(v, 0, 1) ? sstr("%s %s", gene, rc(v, 0, 1)) : ri(v, 0, 2) > 0 ? sstr("%s rs%lld", gene, ri(v, 0, 2))
                                                                  : sstr("%s %s", gene, fmtPos((int)ri(v, 0, 3), ri(v, 0, 4))));
    char g[40];
    snprintf(g, sizeof g, "%s", gene);
    rowsFree(&v);
    lkVariant(vid);
    lkGene(g);
    lkOpen();
}
static void lkGeneMenu(const char *sym) {
    lkBegin(sym);
    lkGene(sym);
    lkOpen();
}

/* ---- the menu ---- */
static void lkDraw(Rectangle r, void *u) {
    (void)u;
    /* a menu taller than the window is capped to it and its rows scroll under
       the title; checked every frame, as the window or the interface size can
       change while the menu is open */
    float cap = fminf(g_lk.fullH, (float)gxScreenH() - 16);
    if (g_gx.popupRect.height != cap) g_gx.popupRect.height = cap;
    float maxS = fmaxf(0, g_lk.fullH - r.height);
    if (g_gx.wheel != 0 && CheckCollisionPointRec(g_gx.mouse, r)) g_lk.scroll -= g_gx.wheel * 52;
    g_lk.scroll = clampf(g_lk.scroll, 0, maxS);
    txtIn(F_SEMI, TH->fsSmall, (Rectangle){ r.x + 14, r.y + 10, r.width - 28, 18 }, TH->text, 0, g_lk.title);
    Rectangle body = { r.x, r.y + 36, r.width, r.height - 36 };
    gxClipPush(body);
    /* the second press of the double-click that opened the menu must not pick
       the entry that has just appeared under the mouse */
    bool dbl = g_gx.dbl;
    float y = body.y - g_lk.scroll;
    for (int i = 0; i < g_lk.n; i++) {
        LkItem *t = &g_lk.it[i];
        if (t->kind == LK_HEAD) {
            y += 4;
            txt(F_SEMI, TH->fsTiny, r.x + 14, y, TH->faint, t->label);
            y += 18;
            continue;
        }
        Rectangle row = { r.x + 6, y, r.width - 12, 26 };
        if (row.y + row.height <= body.y || row.y >= body.y + body.height) { y += 26; continue; }   /* scrolled out of view */
        bool h;
        if (dbl) h = gxHover(row);
        else if (gxBehave(gxHash(LK_POPUP, &i, sizeof i), row, &h, NULL)) {
            LkItem it = *t;
            gxClipPop();
            gxPopupClose();
            switch (it.kind) {
            case LK_URL: openUrl(it.arg); break;
            case LK_DETAILS: variantsShow(0, it.arg, it.id); break;
            case LK_READS: navVariant(it.id); break;
            case LK_GENE: navGeneSym(it.arg); break;
            case LK_ENTER: stEnter(it.arg); break;
            case LK_LOCUS: { long long a = 0, b = 0; sscanf(it.arg, "%lld %lld", &a, &b); long long pad = (b - a) / 5 + 200; navLocus((int)it.id, a - pad, b + pad); break; }
            case LK_FN: if (it.fn) it.fn(it.id); break;
            default: break;
            }
            return;
        }
        if (h) { drRound(row, 5, TH->hover); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
        txtIn(F_REG, TH->fsSmall, (Rectangle){ row.x + 8, row.y, row.width - 180, row.height }, TH->text, 0, t->label);
        txtIn(F_REG, TH->fsTiny, (Rectangle){ row.x + row.width - 172, row.y, 146, row.height }, TH->faint, 2, t->sub);
        if (t->kind == LK_URL) {   /* an arrow out of a box: leaves the app */
            float ax = row.x + row.width - 12, ay = row.y + 9;
            drLine(ax - 5, ay + 7, ax + 2, ay, 1.3f, TH->dim);
            drLine(ax - 2, ay, ax + 2, ay, 1.3f, TH->dim);
            drLine(ax + 2, ay, ax + 2, ay + 4, 1.3f, TH->dim);
        }
        y += 26;
    }
    gxClipPop();
    if (maxS > 0) {   /* how far down the list the view is */
        float barH = fmaxf(24, body.height * body.height / (g_lk.fullH - 36));
        drRound((Rectangle){ r.x + r.width - 6, body.y + (body.height - barH) * (g_lk.scroll / maxS), 3, barH }, 1.5f, TH->lineHi);
    }
}
static void lkOpen(void) {
    if (!g_lk.n) return;
    float h = 46;
    for (int i = 0; i < g_lk.n; i++) h += g_lk.it[i].kind == LK_HEAD ? 22 : 26;
    g_lk.fullH = h;
    g_lk.scroll = 0;
    float W = (float)gxScreenW(), H = (float)gxScreenH(), w = 400, hh = fminf(h, H - 16);
    Vector2 m = g_gx.mouse;
    /* beside the click rather than over it: to the left of the mouse when it
       does not fit to the right, and above when it fits there but not below */
    float x = m.x + 4 + w > W - 8 ? fmaxf(8, m.x - 4 - w) : m.x + 4;
    float y = m.y + 4 + hh > H - 8 && m.y - 4 - hh >= 8 ? m.y - 4 - hh : m.y + 4;
    gxPopupOpen(LK_POPUP, (Rectangle){ m.x - 2, m.y - 2, 4, 4 }, (Rectangle){ x, y, w, hh }, lkDraw, NULL);
}

/* a whole box as one button (after its own widgets); a faint outline on hover */
static bool gxBoxClick(GxId id, Rectangle r) {
    bool h;
    bool c = gxBehave(id, r, &h, NULL);
    if (h) { drRoundLine(r, 10, TH->lineHi); g_gx.cursor = MOUSE_CURSOR_POINTING_HAND; }
    return c;
}
