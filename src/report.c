/* ========================================================================
   src/report.c - a printable summary for a doctor or pharmacist.

   One HTML file (inline styles, prints on A4) with what a clinician would
   want first: the medication genes and the guidance that differs from
   standard dosing, the HLA alleles that matter for drug reactions,
   ClinVar pathogenic findings read against inheritance, carrier status,
   the well-known single variants, and the lineages. Written only when the
   person asks, into their Documents folder (XDG_DOCUMENTS_DIR when set),
   through the engine's atomic commit.
   ======================================================================== */

static void htmlEsc(FILE *f, const char *s) {
    for (; s && *s; s++) {
        switch (*s) {
        case '<': fputs("&lt;", f); break;
        case '>': fputs("&gt;", f); break;
        case '&': fputs("&amp;", f); break;
        case '"': fputs("&quot;", f); break;
        case '|': fputs("; ", f); break;
        case '\n': fputs("<br>", f); break;
        default: fputc(*s, f);
        }
    }
}
static void htmlRows(FILE *f, const char *sql, const char *const *heads, int ncol) {
    Rows *r = dbq(sql, NULL);
    fputs("<table><tr>", f);
    for (int c = 0; c < ncol; c++) { fputs("<th>", f); htmlEsc(f, heads[c]); fputs("</th>", f); }
    fputs("</tr>", f);
    for (int i = 0; i < r->nrow; i++) {
        fputs("<tr>", f);
        for (int c = 0; c < ncol && c < r->ncol; c++) { fputs("<td>", f); htmlEsc(f, rc(r, i, c)); fputs("</td>", f); }
        fputs("</tr>", f);
    }
    if (!r->nrow) fprintf(f, "<tr><td colspan=\"%d\" class=\"none\">None found.</td></tr>", ncol);
    fputs("</table>", f);
    rowsFree(&r);
}
static bool reportDocsDir(char *out, size_t n) {
    const char *x = getenv("XDG_DOCUMENTS_DIR"), *home = getenv("HOME");
    if (x && *x) snprintf(out, n, "%s", x);
    else if (home && *home) snprintf(out, n, "%s/Documents", home);
    else return false;
    mkdir(out, 0755);
    struct stat st;
    return stat(out, &st) == 0 && S_ISDIR(st.st_mode);
}
static void exportReport(void) {
    char dir[1024], path[1200], tmp[1300];
    if (!reportDocsDir(dir, sizeof dir)) { toast("No Documents folder to write into"); return; }
    time_t now = time(NULL);
    struct tm tmv;
    localtime_r(&now, &tmv);
    char day[32];
    strftime(day, sizeof day, "%Y-%m-%d", &tmv);
    snprintf(path, sizeof path, "%s/Strand Open summary %s.html", dir, day);
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { toast("Could not write the summary"); return; }
    fputs("<!doctype html><html><head><meta charset=\"utf-8\"><title>Genome summary</title><style>"
          "body{font-family:Inter,Helvetica,Arial,sans-serif;max-width:900px;margin:32px auto;color:#1c2330;font-size:13px;line-height:1.45}"
          "h1{font-size:24px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 8px;border-bottom:2px solid #6eb2ce;padding-bottom:4px}"
          "table{border-collapse:collapse;width:100%;margin:6px 0}th,td{text-align:left;padding:5px 8px;border-bottom:1px solid #dde2ea;vertical-align:top}"
          "th{background:#eef2f6;font-size:11px;text-transform:uppercase;letter-spacing:.04em}.none{color:#888}.note{color:#555;font-size:12px}"
          ".meta{color:#555}@media print{body{margin:12mm}h2{break-after:avoid}}</style></head><body>", f);
    fputs("<h1>Genome summary</h1><div class=\"meta\">", f);
    fputs("Sample ", f); htmlEsc(f, dbMeta("sample_id"));
    fputs(" &middot; ", f); htmlEsc(f, dkName());
    fprintf(f, ", GRCh38 &middot; prepared %s with Strand Open", day);
    fprintf(f, "</div><p class=\"note\">For discussion with a doctor, pharmacist or genetic counsellor. Research-grade data from %s: "
            "confirm anything that would change care with a clinical test.</p>", dkChip() ? "a consumer genotyping chip" : "consumer sequencing");

    fputs("<h2>Medicines: pharmacogenes</h2>", f);
    { static const char *const h[] = { "Gene", "Result", "Phenotype" };
      htmlRows(f, "SELECT gene, diplotype, phenotype FROM pgx_gene WHERE called=1 ORDER BY gene", h, 3); }
    fprintf(f, "<p class=\"note\">Star alleles by %s; CYP2D6 by Cyrius from the aligned reads (%s); HLA by T1K.</p>", dbMeta("pgx_tool"), dbMeta("cyp2d6_cyrius"));
    fputs("<h2>Medicines: guidance that differs from standard use</h2>", f);
    { static const char *const h[] = { "Medicine", "Source", "Your genotype", "Guidance" };
      htmlRows(f, "SELECT drug, CASE WHEN source LIKE 'CPIC%' THEN 'CPIC' WHEN source LIKE 'DPWG%' THEN 'DPWG' ELSE source END, phenotype, recommendation "
                  "FROM pgx_drug WHERE recommendation!='' AND source NOT LIKE 'FDA%' "
                  "AND (phenotype LIKE '%Poor%' OR phenotype LIKE '%Intermediate%' OR phenotype LIKE '%Rapid%' OR phenotype LIKE '%Decreased%' "
                  "OR phenotype LIKE '%Increased%' OR phenotype LIKE '%positive%' OR phenotype LIKE '%Deficien%') "
                  "AND recommendation NOT LIKE '%does not provide a recommendation%' AND recommendation NOT LIKE 'Initiate therapy with recommended starting dose%' "
                  "AND recommendation NOT LIKE 'Prescribe desired starting dose%' AND recommendation NOT LIKE '%per standard dosing%' "
                  "AND recommendation NOT LIKE 'No recommendation%' AND recommendation NOT LIKE '%no reason to avoid%' "
                  "AND recommendation NOT LIKE 'Initiate standard%' AND recommendation NOT LIKE '%no action is needed%' AND recommendation NOT LIKE '%no action is required%' "
                  "AND recommendation NOT LIKE '%use at standard dose%' AND recommendation NOT LIKE 'No action recommended%' AND recommendation NOT LIKE 'Initiate therapy with standard%' "
                  "GROUP BY drug, source ORDER BY drug", h, 4); }
    fputs("<h2>HLA type</h2>", f);
    { static const char *const h[] = { "Gene", "Allele 1", "Allele 2" };
      htmlRows(f, "SELECT gene, a1, a2 FROM hla WHERE gene IN ('HLA-A','HLA-B','HLA-C','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DPB1')", h, 3); }

    fputs("<h2>ClinVar pathogenic and likely pathogenic variants</h2>", f);
    { static const char *const h[] = { "Gene", "Variant", "Copies", "ClinVar", "Conditions", "Inheritance" };
      htmlRows(f, "SELECT v.gene, coalesce(v.hgvs_p, 'chr'||(CASE v.chrom WHEN 23 THEN 'X' WHEN 24 THEN 'Y' WHEN 25 THEN 'MT' ELSE v.chrom END)||':'||v.pos||' '||v.ref||'>'||v.alt), "
                  "CASE v.zyg WHEN 1 THEN 'heterozygous' WHEN 2 THEN 'homozygous' WHEN 3 THEN 'hemizygous' ELSE 'het (2 alts)' END, v.clin_sig, v.clin_dn, "
                  "(SELECT group_concat(DISTINCT d.inheritance) FROM gene_disease d WHERE d.gene=v.gene_id) FROM variant v WHERE v.clin IN (5,6) ORDER BY v.gene", h, 6); }
    fputs("<p class=\"note\">Screened: every ClinVar-classified variant in the genome. A single copy in a recessive (AR) gene means carrier status.</p>", f);
    fputs("<h2>Well-known variants</h2>", f);
    { static const char *const h[] = { "Variant", "Result", "Genotype" };
      htmlRows(f, "SELECT title, result, genotype FROM trait WHERE cat='Health' ORDER BY ord", h, 3); }
    fputs("<h2>Lineage and ancestry</h2><table>", f);
    fputs("<tr><th>Maternal line (mtDNA)</th><td>", f); htmlEsc(f, dbMeta("mt_haplogroup")); fputs("</td></tr>", f);
    fputs("<tr><th>Paternal line (Y)</th><td>", f); htmlEsc(f, dbMeta("y_haplogroup")); fputs("</td></tr>", f);
    fputs("<tr><th>Ancestry</th><td>", f); htmlEsc(f, *dbMeta("ancestry_summary") ? dbMeta("ancestry_summary") : "-"); fputs("</td></tr></table>", f);
    fputs("<p class=\"note\">Sources: ClinVar (NCBI), PharmCAT with CPIC and DPWG guidelines, PharmGKB, ACMG SF v3.3, GWAS Catalog, gnomAD, 1000 Genomes, PhyloTree, YFull.</p>", f);
    fputs("</body></html>\n", f);
    if (!commitFile(f, tmp, path)) { toast("Could not save the summary"); return; }
    toast("Saved to Documents: Strand summary %s.html", day);
    if (getenv("STRAND_NO_OPEN")) { printf("report: %s\n", path); return; }
    char url[1300];
    snprintf(url, sizeof url, "file://%s", path);
    char *argv[] = { "xdg-open", url, NULL };
    Proc p = spawnReader(argv);
    if (p.pid) reap(&p, false);
}
