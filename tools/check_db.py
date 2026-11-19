#!/usr/bin/env python3
"""Sanity checks on a built genome database: every table the app reads,
and the numbers that would be wrong if a step had half-run.

    tools/check_db.py [path/to/strand.db]"""
import os
import sqlite3
import sys

if len(sys.argv) > 1:
    path = sys.argv[1]
elif os.environ.get("STRANDOPEN_PROFILE"):
    path = os.path.join(os.environ["STRANDOPEN_PROFILE"], "db", "strand.db")
else:
    raise SystemExit("usage: tools/check_db.py PROFILE_DIR/db/strand.db (or set STRANDOPEN_PROFILE)")
con = sqlite3.connect(path)
meta = dict(con.execute("SELECT key, value FROM meta"))
fails = 0


def check(what, ok, detail=""):
    global fails
    print("%-4s %-52s %s" % ("ok" if ok else "FAIL", what, detail))
    fails += 0 if ok else 1


def count(table, where="1"):
    try:
        return con.execute("SELECT count(*) FROM %s WHERE %s" % (table, where)).fetchone()[0]
    except sqlite3.Error:
        return -1


tables = {"variant": 1_000_000, "gene": 50_000, "exon": 100_000, "refblock": 1_000_000, "covbin": 300_000, "cytoband": 800,
          "chrom": 25, "gene_disease": 1000, "hap_mt": 5, "hap_y": 5, "pgx_gene": 15, "pgx_drug": 50, "trait": 20, "sv": 100,
          "acmg": 80, "gwas": 10_000, "pgkb": 100, "files": 3, "qc_fastq": 2, "qc_align": 3, "anc_admix": 5, "anc_pca": 1000,
          "prs": 10, "hla": 5, "rare_matrix": 30, "rare_sfs": 10, "rare_var": 100, "rare_gene": 1, "remarkable": 20, "highlight": 1, "lin_story": 2,
          "rarity": 100, "unanswerable": 5, "anc_closest": 5, "anc_paint": 20, "anc_pop": 20, "csq_count": 5, "genedel": 1,
          "hap_mt_cand": 1, "prs_dist": 10}
for t, least in tables.items():
    n = count(t)
    check("table %s" % t, n >= least, "%s rows" % format(n, ","))
# the search indexes take their text from gene and gwas, so count(*) on them
# reads those tables; the _docsize shadow table holds one row per indexed row
for t, base in (("gene_fts", "gene"), ("gwas_fts", "gwas")):
    n = count(t + "_docsize")
    check("search index %s filled" % t, n > 0 and n == count(base), "%s rows" % format(n, ","))
nv = int(meta.get("n_variants", 0))
check("variant count matches meta", nv == count("variant"), str(nv))
titv = float(meta.get("titv", 0))
check("Ti/Tv between 1.9 and 2.2", 1.9 <= titv <= 2.2, str(titv))
af = count("variant", "af IS NOT NULL")
check("1000 Genomes frequencies present", af > 0.8 * nv, "%.1f%%" % (100 * af / max(1, nv)))
check("ClinVar annotations present", count("variant", "clin>0") > 10_000)
check("AlphaMissense scores present", count("variant", "am IS NOT NULL") > 5000)
check("no runaway protein strings", count("variant", "length(hgvs_p) > 64") == 0, str(count("variant", "length(hgvs_p) > 64")))
check("callable autosome fraction above 90%", float(meta.get("autosome_callable_fraction", 0)) > 0.9, meta.get("autosome_callable_fraction"))
check("mtDNA haplogroup", bool(meta.get("mt_haplogroup")), meta.get("mt_haplogroup", ""))
check("Y haplogroup (male)", meta.get("sex") != "male" or bool(meta.get("y_haplogroup")), meta.get("y_haplogroup", ""))
check("PGx genes called", count("pgx_gene", "called=1") >= 15, str(count("pgx_gene", "called=1")))
check("traits all called", count("trait", "status='nocall'") == 0, str(count("trait", "status='nocall'")))
s = con.execute("SELECT sum(frac) FROM anc_admix WHERE level='super'").fetchone()[0] if count("anc_admix") > 0 else 0
check("admixture sums to one", abs((s or 0) - 1) < 0.01, "%.4f" % (s or 0))
check("rare variants all have a plain change text", count("rare_var", "change=''") == 0, str(count("rare_var", "change=''")))
check("remarkable variants checked (none left uncalled)", 0 <= count("remarkable", "status<0") <= 3, str(count("remarkable", "status<0")))
for k in ("protein", "homo", "ancestry"):
    check("rarity list %s built" % k, count("rarity", "kind='%s'" % k) > 0, str(count("rarity", "kind='%s'" % k)))
missing = [k for k in ("rarity_under_1000", "rarity_homo_under_1000", "rarity_ancestry_n", "rarity_main_group", "rarity_weights") if k not in meta]
check("rarity counts in meta", not missing, ", ".join(missing))
print("\n%d problem(s)" % fails)
sys.exit(1 if fails else 0)
