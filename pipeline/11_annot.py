#!/usr/bin/env python3
"""Step 11: research and clinical annotations.

- acmg: the ACMG SF v3.3 genes (data/acmg_sf_v3.3.tsv, from the policy
  statement's Table 1), the genes a clinical genome reports incidental
  findings in.
- gwas: GWAS Catalog associations at genome-wide significance (p < 5e-8)
  where the sample's genotype is known, with how many copies of the reported
  risk allele the sample carries. A risk allele that is neither the reference
  base nor an allele seen at the site (sample or 1000 Genomes) is left out:
  its strand cannot be trusted.
- pgkb: PharmGKB clinical annotations (levels 1A to 3) for rsID variants,
  with the text for the sample's own genotype.

Tables: acmg, gwas, pgkb."""
import csv
import re
import io
import subprocess
import zipfile
from strandlib import *
from genotype import Genotyper

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
con = db()
cur = con.cursor()
G = Genotyper(con)
csv.field_size_limit(10 ** 9)

_fai = {}
with open(FA + ".fai") as f:
    for line in f:
        p = line.split("\t")
        _fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
_fa = open(FA, "rb")


def refbase(chrom, p):
    ln, off, lb, lw = _fai[chrom]
    z = p - 1
    _fa.seek(off + z // lb * lw + z % lb)
    return _fa.read(1).decode().upper()


# ---- ACMG SF
cur.executescript("DROP TABLE IF EXISTS acmg; CREATE TABLE acmg(gene TEXT, category TEXT, disease TEXT, inheritance TEXT, report TEXT);")
with open(os.path.join(ROOT, "data", "acmg_sf_v3.3.tsv")) as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) == 5:
            cur.execute("INSERT INTO acmg VALUES(?,?,?,?,?)", p)
log("ACMG SF genes", cur.execute("SELECT count(DISTINCT gene) FROM acmg").fetchone()[0])

# 1000 Genomes alleles, when step 4 has run, to vouch for risk alleles
kgalt = {}

# ---- GWAS Catalog
cur.executescript("""DROP TABLE IF EXISTS gwas;
CREATE TABLE gwas(rsid TEXT, chrom INTEGER, pos INTEGER, risk TEXT, copies INTEGER, ploidy INTEGER, trait TEXT, mapped TEXT,
  p REAL, mlog REAL, effect TEXT, ci TEXT, raf TEXT, pmid TEXT, author TEXT, year TEXT, gene TEXT, study TEXT, n TEXT);""")
zf = zipfile.ZipFile(os.path.join(RAW, "gwas-catalog-associations_ontology-annotated-full.zip"))
name = [n for n in zf.namelist() if n.endswith(".tsv")][0]
cand = []      # (record fields..., chrom, pos, risk, ref, sample alleles)
with zf.open(name) as fh:
    rd = csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8", errors="replace"), delimiter="\t")
    for r in rd:
        try:
            mlog = float(r.get("PVALUE_MLOG") or 0)
        except ValueError:
            continue
        if mlog < 7.30103:
            continue
        sra = r.get("STRONGEST SNP-RISK ALLELE", "")
        if ";" in sra or " x " in sra or "-" not in sra:
            continue
        rs, risk = sra.rsplit("-", 1)
        risk = risk.strip()
        if risk not in ("A", "C", "G", "T"):
            continue
        c = (r.get("CHR_ID") or "").strip()
        if c not in CHROM_IDX or ";" in (r.get("CHR_POS") or ""):
            continue
        try:
            pos = int(r["CHR_POS"])
        except (ValueError, KeyError):
            continue
        cand.append((r, rs, risk, c, pos, mlog))
log("GWAS candidate associations", len(cand))
# the 1000 Genomes alleles at those positions vouch for risk alleles that are
# neither the reference nor carried (the sample is then 0 copies)
kgalt = {}
if kg_site_files():
    for line in kg_query([(x[3], x[4]) for x in cand], "%CHROM\t%POS\t%REF\t%ALT\n"):
        q = line.split("\t")
        if len(q) == 4:
            kgalt.setdefault((q[0], int(q[1])), set()).update([q[2], q[3]])
    log("1000 Genomes alleles found at", len(kgalt), "GWAS positions")
rows = []
seen = set()
kept = skipped = 0
for r, rs, risk, c, pos, mlog in cand:
    ci = CHROM_IDX[c]
    ref = refbase(c, pos)
    al = G.alleles(ci, pos, ref)
    if not al:
        skipped += 1
        continue
    if risk != ref and risk not in al and risk not in kgalt.get((c, pos), ()):
        skipped += 1      # strand or allele cannot be confirmed
        continue
    copies = sum(1 for a in al if a == risk)
    key = (rs, r.get("MAPPED_TRAIT") or r.get("DISEASE/TRAIT"), r.get("PUBMEDID"))
    if key in seen:
        continue
    seen.add(key)
    try:
        pv = float(r.get("P-VALUE") or "nan")
    except ValueError:
        pv = None
    rows.append((rs, ci, pos, risk, copies, len(al), r.get("DISEASE/TRAIT", ""), r.get("MAPPED_TRAIT", ""), pv, mlog,
                 r.get("OR or BETA", ""), r.get("95% CI (TEXT)", ""), r.get("RISK ALLELE FREQUENCY", ""), r.get("PUBMEDID", ""),
                 r.get("FIRST AUTHOR", ""), (r.get("DATE") or "")[:4], r.get("MAPPED_GENE", ""), r.get("STUDY", "")[:300],
                 r.get("INITIAL SAMPLE SIZE", "")[:200]))
    kept += 1
    if len(rows) >= 50000:
        cur.executemany("INSERT INTO gwas VALUES(" + ",".join("?" * 19) + ")", rows)
        rows.clear()
if rows:
    cur.executemany("INSERT INTO gwas VALUES(" + ",".join("?" * 19) + ")", rows)
cur.execute("CREATE INDEX gwas_rs ON gwas(rsid)")
cur.execute("CREATE INDEX gwas_pos ON gwas(chrom, pos)")
cur.execute("DROP TABLE IF EXISTS gwas_fts")
cur.execute("CREATE VIRTUAL TABLE gwas_fts USING fts5(trait, mapped, gene, content='gwas', content_rowid='rowid')")
cur.execute("INSERT INTO gwas_fts(rowid, trait, mapped, gene) SELECT rowid, trait, mapped, gene FROM gwas")
log("GWAS associations kept", kept, "skipped (no call or unverifiable strand)", skipped)

# ---- PharmGKB clinical annotations
locs = {}
with zipfile.ZipFile(os.path.join(RAW, "variants.zip")) as z:
    with z.open("variants.tsv") as fh:
        for r in csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8", errors="replace"), delimiter="\t"):
            loc = r.get("Location") or ""
            m = re.match(r"NC_0000(\d\d)\.\d+:(\d+)", loc)
            if m and r["Variant Name"].startswith("rs"):
                n = int(m.group(1))
                c = str(n) if n <= 22 else "X" if n == 23 else "Y"
                locs[r["Variant Name"]] = (c, int(m.group(2)))
import re as _re
cur.executescript("""DROP TABLE IF EXISTS pgkb;
CREATE TABLE pgkb(id TEXT, rsid TEXT, gene TEXT, level TEXT, category TEXT, drugs TEXT, phenotypes TEXT,
  genotype TEXT, text TEXT, function TEXT, url TEXT);""")
with zipfile.ZipFile(os.path.join(RAW, "clinicalAnnotations.zip")) as z:
    anns = {}
    with z.open("clinical_annotations.tsv") as fh:
        for r in csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8", errors="replace"), delimiter="\t"):
            anns[r["Clinical Annotation ID"]] = r
    alle = {}
    with z.open("clinical_ann_alleles.tsv") as fh:
        for r in csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8", errors="replace"), delimiter="\t"):
            alle.setdefault(r["Clinical Annotation ID"], []).append(r)
n = 0
for aid, a in anns.items():
    var = a.get("Variant/Haplotypes", "")
    lvl = a.get("Level of Evidence", "")
    if not var.startswith("rs") or lvl not in ("1A", "1B", "2A", "2B", "3"):
        continue
    if var not in locs:
        continue
    c, pos = locs[var]
    ref = refbase(c, pos)
    al = G.alleles(CHROM_IDX[c], pos, ref)
    if not al:
        continue
    mine = "".join(sorted(al))
    text = func = gt = ""
    for x in alle.get(aid, []):
        g = x.get("Genotype/Allele", "").strip()
        if "".join(sorted(g)) == mine or (len(al) == 1 and g == al[0]):
            text, func, gt = x.get("Annotation Text", ""), x.get("Allele Function", ""), g
            break
    if not text:
        continue
    cur.execute("INSERT INTO pgkb VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                (aid, var, a.get("Gene", ""), lvl, a.get("Phenotype Category", ""), a.get("Drug(s)", ""), a.get("Phenotype(s)", ""),
                 gt, text, func, a.get("URL", "")))
    n += 1
log("PharmGKB clinical annotations matched to your genotype", n)
con.commit()
