#!/usr/bin/env python3
"""Step 20: runs of homozygosity, or how closely related your parents were.

Where both copies of a stretch of a chromosome came from one shared
ancestor, every position in it is homozygous. Short runs (under a couple of
megabases) are in everyone: they come from ancestors shared hundreds of
generations back. Long runs (5 to 10 Mb and more) come from ancestors
shared within the last few generations, so their total length says how
closely related the parents were (offspring of first cousins carry about a
sixteenth of the genome, roughly 180 Mb, in such runs).

bcftools roh (an HMM over the called variants with the allele frequencies
of your own 1000 Genomes continental group) finds the runs.

Table roh(chrom, start, end, len, nmark, qual); meta roh_* and roh_summary.
"""
import gzip
import subprocess
from strandlib import *

con = db()
cur = con.cursor()
src = os.path.join(WORK, "variants.af.vcf.gz")
if not os.path.exists(src):
    log("no frequency-annotated variants (step 5): skipped")
    raise SystemExit(0)
main = "EUR"
try:
    r = cur.execute("SELECT pop FROM anc_admix WHERE level='super' ORDER BY frac DESC LIMIT 1").fetchone()
    if r and r[0] in ("EUR", "AFR", "EAS", "SAS", "AMR"):
        main = r[0]
except Exception:
    pass
tag = "KG" + main
out = os.path.join(WORK, "roh")
os.makedirs(out, exist_ok=True)
res = os.path.join(out, "roh.txt.gz")
subprocess.run(["bcftools", "roh", "--AF-tag", tag, "-G30", "--threads", "8", "-Orz", "-o", res,
                "-r", ",".join(str(i) for i in range(1, 23)), src], check=True, capture_output=True)

cur.executescript("DROP TABLE IF EXISTS roh; CREATE TABLE roh(chrom INTEGER, start INTEGER, end INTEGER, len INTEGER, nmark INTEGER, qual REAL);")
runs = []
with gzip.open(res, "rt") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] != "RG" or p[2] not in CHROM_IDX:
            continue
        c, a, b, n, q = CHROM_IDX[p[2]], int(p[3]), int(p[4]), int(p[6]), float(p[7])
        runs.append((c, a, b, b - a + 1, n, q))
cur.executemany("INSERT INTO roh VALUES(?,?,?,?,?,?)", runs)


def total(minlen):
    return sum(r[3] for r in runs if r[3] >= minlen)


n1, t1 = sum(1 for r in runs if r[3] >= 1_000_000), total(1_000_000)
n5, t5 = sum(1 for r in runs if r[3] >= 5_000_000), total(5_000_000)
longest = max((r[3] for r in runs), default=0)
mb = lambda x: x / 1e6
if t5 >= 150e6:
    text = ("Long runs of matching DNA add up to about %.0f Mb, the pattern expected when parents are related about as closely as first cousins, "
            "or come from a small community where everyone is distantly related." % mb(t5))
elif t5 >= 40e6:
    text = ("Long runs add up to about %.0f Mb: consistent with parents who are distant relatives (around second cousins), "
            "or with ancestry from a small or isolated community." % mb(t5))
elif n5:
    text = ("A few long runs (%d over 5 Mb, %.0f Mb in all) point to an ancestor shared by both parents several generations back, "
            "which is common and usually unknown to families." % (n5, mb(t5)))
else:
    text = ("No long runs of matching DNA: no sign that your parents were related. Your longest run is %.1f Mb; short runs like these "
            "are in everyone and come from ancestors shared hundreds of generations ago." % mb(longest))
for k, v in (("roh_n_1mb", n1), ("roh_total_1mb", t1), ("roh_n_5mb", n5), ("roh_total_5mb", t5), ("roh_longest", longest)):
    set_meta(con, k, str(v))
set_meta(con, "roh_summary", text)
set_meta(con, "roh_af_group", main)
con.commit()
log(text)
