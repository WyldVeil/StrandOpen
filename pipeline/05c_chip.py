#!/usr/bin/env python3
"""Step 5c (chip data only): rare chip calls are set aside as unconfirmed.

Genotyping chips are accurate for common variants and unreliable for rare
ones: in UK Biobank (nearly 50,000 people with both chip and exome data)
most variants rarer than 1 in 100,000 that the chips reported were not
there when sequenced (Weedon et al. 2021, BMJ 372:n214). So a chip call of a
variant rarer than 1 in 1,000 (or missing from the reference panel) is
marked FILTER=ChipRare and pass=0, and every step after this one reads its
position as untested. A ClinVar class it carries moves to 100 + class
("unconfirmed chip call"), so it stays visible, labelled, but out of every
list of findings.

Frequencies: 1000 Genomes (step 5) when the population reference is
installed; otherwise the ExAC, ESP and 1000 Genomes frequencies ClinVar
publishes with each variant, for ClinVar-listed variants only (others are
left alone: without frequencies nothing says they are rare).

Meta: chip_rare (calls set aside), chip_rare_clinvar (of them, ClinVar-listed
as pathogenic or likely pathogenic)."""
import subprocess
from strandlib import *

con = db()
cur = con.cursor()
cm = cur.execute("SELECT value FROM meta WHERE key='call_model'").fetchone()
if not cm or cm[0] != "chip":
    log("not chip data: nothing to do")
    raise SystemExit(0)
RARE = 0.001
has_kg = cur.execute("SELECT 1 FROM variant WHERE af IS NOT NULL LIMIT 1").fetchone() is not None
# chromosomes 1-22 and X only: the 1000 Genomes panel has no Y or mitochondrial
# frequencies, and the lineage callers read many markers at once, so a few
# wrong Y or MT calls cannot move them
if has_kg:
    rare = [r[0] for r in cur.execute("SELECT id FROM variant WHERE zyg>0 AND chrom<=23 AND (af IS NULL OR af < ?)", (RARE,))]
else:
    # ClinVar's own frequency fields, for the ClinVar-listed calls
    rows = cur.execute("SELECT id, chrom, pos, ref, alt FROM variant WHERE zyg>0 AND chrom<=23 AND clin_id IS NOT NULL").fetchall()
    names = {i + 1: c for i, c in enumerate(CHROMS)}
    reg = os.path.join(WORK, "chip_clinvar_regions.tsv")
    with open(reg, "w") as f:
        for vid, c, p, r, a in sorted(rows, key=lambda x: (x[1], x[2])):
            f.write("%s\t%d\n" % (names[c], p))
    freq = {}
    q = subprocess.run(["bcftools", "query", "-R", reg, "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_EXAC\t%INFO/AF_TGP\t%INFO/AF_ESP\n",
                        os.path.join(RAW, "clinvar.vcf.gz")], capture_output=True, text=True)
    for line in q.stdout.splitlines():
        c, p, r, a, *fs = line.split("\t")
        vals = [float(x) for x in fs if x not in (".", "")]
        freq[(c, int(p), r, a)] = max(vals) if vals else None
    rare = []
    for vid, c, p, r, a in rows:
        f = freq.get((names[c], p, r, a))
        if f is None or f < RARE:
            rare.append(vid)
log("chip calls set aside as rare (unconfirmed): %d (frequencies from %s)" % (len(rare), "1000 Genomes" if has_kg else "ClinVar"))
for i in range(0, len(rare), 50000):
    chunk = rare[i:i + 50000]
    q = ",".join("?" * len(chunk))
    cur.execute("UPDATE variant SET filter='ChipRare', pass=0 WHERE id IN (%s)" % q, chunk)
    cur.execute("UPDATE variant SET clin=100+clin WHERE clin BETWEEN 1 AND 11 AND id IN (%s)" % q, chunk)
n_plp = cur.execute("SELECT count(*) FROM variant WHERE clin IN (105, 106)").fetchone()[0]
set_meta(con, "chip_rare", len(rare))
set_meta(con, "chip_rare_clinvar", n_plp)
con.commit()
log("of them ClinVar pathogenic or likely pathogenic:", n_plp)
