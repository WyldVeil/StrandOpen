#!/usr/bin/env python3
"""Step 5b: mapping quality of every call (from the normalised variant VCF
of step 1): variant.mq, the average mapping quality of the reads, and
variant.mq0f, the fraction of reads that could have come from elsewhere in
the genome (mapping quality 0). Low values mark duplicated DNA where a
call is often a mis-mapping artefact; the rare-variant steps use them."""
from strandlib import *

con = db()
cur = con.cursor()
src = os.path.join(WORK, "variants.annot.vcf.gz")
cols = [r[1] for r in cur.execute("PRAGMA table_info(variant)")]
for c, t in (("mq", "REAL"), ("mq0f", "REAL")):
    if c not in cols:
        cur.execute("ALTER TABLE variant ADD COLUMN %s %s" % (c, t))
cur.executescript("DROP TABLE IF EXISTS temp.mq; CREATE TEMP TABLE mq(chrom INTEGER, pos INTEGER, ref TEXT, alt TEXT, mq REAL, mq0f REAL);")
batch = []
n = 0
for line in pipe_lines("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/MQ\\t%%INFO/MQ0F\\n' '%s'" % src):
    c, pos, ref, alt, mq, mq0f = line.rstrip("\n").split("\t")
    ci = CHROM_IDX.get(c)
    if not ci:
        continue
    batch.append((ci, int(pos), ref, alt, None if mq == "." else float(mq), None if mq0f == "." else float(mq0f)))
    if len(batch) >= 200000:
        cur.executemany("INSERT INTO mq VALUES(?,?,?,?,?,?)", batch)
        n += len(batch)
        batch = []
cur.executemany("INSERT INTO mq VALUES(?,?,?,?,?,?)", batch)
n += len(batch)
log("records", n)
cur.execute("CREATE INDEX temp.mq_k ON mq(chrom, pos, ref, alt)")
cur.execute("""UPDATE variant SET mq = m.mq, mq0f = m.mq0f FROM mq m
               WHERE m.chrom = variant.chrom AND m.pos = variant.pos AND m.ref = variant.ref AND m.alt = variant.alt""")
log("matched", cur.execute("SELECT count(*) FROM variant WHERE mq IS NOT NULL").fetchone()[0], "of", cur.execute("SELECT count(*) FROM variant").fetchone()[0])

# rsIDs for clinically known calls that came without one (most indels): ClinVar
# records the dbSNP id of every variation it holds
cv = os.path.join(RAW, "clinvar.vcf.gz")
if os.path.exists(cv):
    rs = {}
    for line in pipe_lines("bcftools query -f '%%ID\\t%%INFO/RS\\n' '%s'" % cv):
        a, b = line.rstrip("\n").split("\t")
        if b not in (".", "") and a.isdigit():
            rs[int(a)] = int(b.split(",")[0].split("|")[0])
    upd = [(rs[c], vid) for vid, c in cur.execute("SELECT id, clin_id FROM variant WHERE (rsid IS NULL OR rsid=0) AND clin_id > 0").fetchall() if c in rs]
    cur.executemany("UPDATE variant SET rsid=? WHERE id=?", upd)
    log("rsIDs filled from ClinVar", len(upd))
con.commit()
