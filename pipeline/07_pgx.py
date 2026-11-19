#!/usr/bin/env python3
"""Step 7: pharmacogenomics with PharmCAT (CPIC and DPWG guidance).

PharmCAT matches star alleles from a VCF holding every position its allele
definitions use. That VCF is written here from the sample's own calls: the
variant table where a variant sits, 0/0 where a gVCF reference block covers
the base, ./. where nothing was called (PharmCAT then says what it could
not decide). CYP2D6 copy number and the HLA genes cannot be called from
these SNV and indel calls; PharmCAT reports them as not called.

Tables: pgx_gene (one row per gene: diplotype, phenotype, activity score,
calls), pgx_drug (each drug's guidance from CPIC and DPWG)."""
import glob
import json
import subprocess
import re
import shutil
from strandlib import *
from genotype import Genotyper

PC = os.path.join(RAW, "pharmcat")
JAR = glob.glob(os.path.join(PC, "pharmcat-*-all.jar"))[0]
POSV = glob.glob(os.path.join(PC, "pharmcat_positions_*.vcf"))[0]
OUT = os.path.join(WORK, "pgx")
shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(OUT)

con = db()
cur = con.cursor()
G = Genotyper(con)
male = G.male

hdr, body = [], []
nref = nvar = nmiss = 0
with open(POSV) as f:
    for line in f:
        if line.startswith("##"):
            hdr.append(line.rstrip("\n"))
            continue
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        c = p[0][3:] if p[0].startswith("chr") else p[0]
        if c == "M":
            c = "MT"
        if c not in CHROM_IDX:
            continue
        ci = CHROM_IDX[c]
        pos = int(p[1])
        ref = p[3]
        alts = [a for a in p[4].split(",") if a and a != "."]
        hap = male and c in ("X", "Y")
        r = G.at(ci, pos)
        gt = "./." if not hap else "."
        if r[0] == "ref":
            gt = "0" if hap else "0/0"
            nref += 1
        elif r[0] == "var":
            idx = []
            zyg = None
            for vref, valt, vz, ps, dp in r[1]:
                # the same event spelled with more reference context (PharmCAT
                # writes the UGT1A1 TA repeat as CAT>CATAT; the calls have C>CAT)
                if len(ref) > len(vref) and ref.startswith(vref):
                    sfx = ref[len(vref):]
                    vref, valt = vref + sfx, valt + sfx
                if vref == ref and valt in alts:
                    idx.append(alts.index(valt) + 1)
                    zyg = vz
                elif vref == ref:
                    alts.append(valt)
                    idx.append(len(alts))
                    zyg = vz
            if idx:
                nvar += 1
                if hap or zyg == 3:
                    gt = str(idx[0])
                elif zyg == 2:
                    gt = "%d/%d" % (idx[0], idx[0])
                elif len(idx) >= 2:
                    gt = "%d/%d" % (idx[0], idx[1])
                else:
                    gt = "0/%d" % idx[0]
            else:
                # a variant here with another representation (an indel overlapping
                # a defined SNV): leave it uncalled rather than guess
                gt = "./." if not hap else "."
                nmiss += 1
        else:
            nmiss += 1
        body.append("\t".join([p[0], p[1], p[2], ref, ",".join(alts) if alts else ".", ".", "PASS", p[7], "GT", gt]))
vcf = os.path.join(OUT, "strand.pgx.vcf")
with open(vcf, "w") as f:
    for h in hdr:
        f.write(h + "\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    for b in body:
        f.write(b + "\n")
log("PharmCAT input: %d reference, %d variant, %d uncalled positions" % (nref, nvar, nmiss))

# Outside calls: what the VCF cannot carry.
#  CYP2D6 from Cyrius (copy number and hybrids read from the CRAM, step 13a)
#  MT-RNR1 from the mitochondrial calls, against PharmCAT's own allele list
outside = []
cy = os.path.join(WORK, "cyp2d6", "strand.tsv")
if os.path.exists(cy):
    for line in open(cy):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3 and p[0] != "Sample" and p[1] not in ("None", "") and p[2] == "PASS":
            outside.append("CYP2D6\t%s" % p[1])
            set_meta(con, "cyp2d6_cyrius", "%s (%s)" % (p[1], p[2]))
# HLA from T1K (reads of the MHC region typed against IPD-IMGT/HLA)
hla_rows = []
ht = os.path.join(WORK, "hla", "strand_genotype.tsv")
if os.path.exists(ht):
    for line in open(ht):
        p = line.rstrip("\n").split("\t")
        if len(p) < 8 or not p[0].startswith("HLA-"):
            continue
        def two_field(a):
            a = a.split(",")[0].strip()
            if "*" not in a:
                return None
            g, rest = a.split("*", 1)
            return "*" + ":".join(rest.split(":")[:2])
        a1, a2 = two_field(p[2]), two_field(p[5]) if p[5] not in (".", "") else None
        q1 = float(p[4]) if p[4] not in (".", "") else 0
        q2 = float(p[7]) if len(p) > 7 and p[7] not in (".", "") else 0
        if a1 and not a2:
            a2 = a1   # T1K reports one allele when both copies are the same
        hla_rows.append((p[0], a1 or "", a2 or "", q1, q2))
        if p[0] in ("HLA-A", "HLA-B") and a1 and q1 >= 10 and (q2 >= 10 or a2 == a1):
            outside.append("%s\t%s/%s" % (p[0], a1, a2))
cur.executescript("DROP TABLE IF EXISTS hla; CREATE TABLE hla(gene TEXT, a1 TEXT, a2 TEXT, q1 REAL, q2 REAL);")
cur.executemany("INSERT INTO hla VALUES(?,?,?,?,?)", hla_rows)
import zipfile
with zipfile.ZipFile(JAR) as z:
    mt = json.loads(z.read("org/pharmgkb/pharmcat/phenotype/MT_RNR1.json").decode("utf-8", "replace"))
known = mt["haplotypes"]
mine = []
for pos, ref, alt in cur.execute("SELECT pos, ref, alt FROM variant WHERE chrom=25 AND pos BETWEEN 648 AND 1601 AND length(ref)=1 AND length(alt)=1"):
    name = "m.%d%s>%s" % (pos, ref, alt)
    if name in known:
        mine.append(name)
rank = lambda n: 2 if "increased" in known[n].lower() else 1 if "uncertain" in known[n].lower() else 0
pick = max(mine, key=rank) if mine and max(rank(n) for n in mine) > 0 else (mine[0] if mine else "Reference")
outside.append("MT-RNR1\t%s" % pick)
ofile = os.path.join(OUT, "outside.tsv")
open(ofile, "w").write("\n".join(outside) + "\n")
log("outside calls:", "; ".join(outside))
with open(os.path.join(OUT, "run.log"), "w") as runlog:     # an argument list: paths with spaces are safe
    rc = subprocess.run(["java", "-jar", JAR, "-vcf", vcf, "-o", OUT, "-bf", "strand", "-po", ofile,
                         "-reporterJson", "-matcher", "-phenotyper", "-reporter"], stdout=runlog, stderr=subprocess.STDOUT).returncode
if rc != 0:
    print(open(os.path.join(OUT, "run.log")).read()[-3000:])
    raise SystemExit("PharmCAT failed")
rep = json.load(open(os.path.join(OUT, "strand.report.json")))
log("report keys", list(rep.keys()))

cur.executescript("""
DROP TABLE IF EXISTS pgx_gene; DROP TABLE IF EXISTS pgx_drug;
CREATE TABLE pgx_gene(gene TEXT, source TEXT, diplotype TEXT, phenotype TEXT, activity TEXT, called INTEGER,
  alleles TEXT, functions TEXT, notes TEXT, missing INTEGER, uncallable TEXT);
CREATE TABLE pgx_drug(drug TEXT, gene TEXT, source TEXT, phenotype TEXT, recommendation TEXT, implication TEXT,
  strength TEXT, guideline TEXT, url TEXT, matched INTEGER, population TEXT);
""")
# genes
genes = rep.get("genes", {})
for gsym, g in genes.items():
    if not isinstance(g, dict):
        continue
    src = g.get("alleleDefinitionSource") or ""
    if True:
        dips = g.get("sourceDiplotypes") or g.get("recommendationDiplotypes") or []
        names, phen, act, fn = [], [], [], []
        for d in dips:
            a1 = (d.get("allele1") or {})
            a2 = (d.get("allele2") or {})
            label = d.get("label") or "/".join(x for x in [a1.get("name"), a2.get("name")] if x)
            names.append(label)
            phen += d.get("phenotypes") or []
            if d.get("activityScore") is not None:
                act.append(str(d.get("activityScore")))
            fn.append(" / ".join(x for x in [a1.get("function"), a2.get("function")] if x))
        called = 1 if g.get("callSource") not in (None, "NONE") and names and not any("Unknown" in n or "No Result" in n for n in names) else 0
        missing = len(g.get("uncalledHaplotypes") or [])
        notes = "; ".join(m.get("message", "") if isinstance(m, dict) else str(m) for m in (g.get("messages") or []))[:1500]
        cur.execute("INSERT INTO pgx_gene VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                    (gsym, src, " or ".join(dict.fromkeys(names)), "; ".join(dict.fromkeys(phen)), ", ".join(dict.fromkeys(act)),
                     called, json.dumps(dips)[:20000], " | ".join(fn), notes, missing,
                     ", ".join(g.get("uncalledHaplotypes") or [])[:2000]))
def clean(t):
    """PharmCAT's DPWG texts carry HTML lists: plain text with dashes"""
    t = re.sub(r"<li>\s*", "\n- ", t or "")
    t = re.sub(r"<[^>]+>", "", t)
    t = t.replace("&nbsp;", " ").replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
    return re.sub(r"\n\s*\n+", "\n", t).strip()


# drugs
drugs = rep.get("drugs", {})
nd = 0
for src, dmap in drugs.items():
    for dname, d in dmap.items():
        for gl in d.get("guidelines", []):
            gurl = (gl.get("url") or "")
            gname = gl.get("name") or ""
            for ann in gl.get("annotations", []):
                genes_in = []
                phenos = []
                for gp in ann.get("genotypes", []) or []:
                    for dip in gp.get("diplotypes", []) or []:
                        gsym = dip.get("gene") or ""
                        genes_in.append(gsym)
                        phenos += [("%s: %s" % (gsym, p)) for p in (dip.get("phenotypes") or [])]
                impl = clean("; ".join(ann.get("implications") or []))
                cur.execute("INSERT INTO pgx_drug VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                            (dname, ", ".join(dict.fromkeys(genes_in)), src, "; ".join(dict.fromkeys(phenos)),
                             clean(ann.get("drugRecommendation") or ""), impl, ann.get("classification") or "",
                             gname, gurl, 1 if ann.get("drugRecommendation") else 0, ann.get("population") or ""))
                nd += 1
# HLA-A*31:01 has a tag SNP in Europeans (rs1061235, chr6:29945521 A>T):
# the read-based HLA call is worth checking against it
tag = cur.execute("SELECT zyg, ad_ref, ad_alt FROM variant WHERE chrom=6 AND pos=29945521 AND ref='A' AND alt='T'").fetchone()
row = cur.execute("SELECT diplotype FROM pgx_gene WHERE gene='HLA-A'").fetchone()
if row and "*31:01" in (row[0] or ""):
    if tag:
        note = ("HLA-A*31:01 typed from the reads (T1K) and supported by its European tag SNP rs1061235, which you carry "
                "(%s, %d of %d reads)." % ("one copy" if tag[0] in (1, 4) else "two copies", tag[2] or 0, (tag[1] or 0) + (tag[2] or 0)))
    else:
        note = "HLA-A*31:01 typed from the reads (T1K); its European tag SNP rs1061235 was not seen, so treat this call with caution."
    cur.execute("UPDATE pgx_gene SET notes = ? || CASE WHEN notes!='' THEN ' ' || notes ELSE '' END WHERE gene='HLA-A'", (note,))
log("pgx genes", cur.execute("SELECT count(*) FROM pgx_gene").fetchone()[0], "drug annotations", nd)
set_meta(con, "pgx_tool", "PharmCAT " + os.path.basename(JAR).split("-")[1] + " (CPIC and DPWG)")
con.commit()
