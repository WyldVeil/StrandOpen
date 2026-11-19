#!/usr/bin/env python3
"""Step 12: the sample, its files and the read-level quality summary.

- sex from the data: share of heterozygous calls on chromosome X outside
  the pseudoautosomal regions, and how much of Y could be called
- files: every source file with its size and what it holds
- qc_fastq: tools/fqstat's summary of each FASTQ (work/qc/r1.json, r2.json)

Tables: files, qc_fastq; meta keys sample_id, seq_kind, sex, built."""
import json
import time
from strandlib import *

con = db()
cur = con.cursor()

def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path, encoding="utf-8", errors="replace"):
            if "=" in line:
                k, v = line.rstrip("\n").split("=", 1)
                out[k] = v
    return out


prof = kv(os.path.join(DATA, "profile.kv"))
ing = kv(os.path.join(CANON, "ingest.kv"))
meta = dict(cur.execute("SELECT key, value FROM meta").fetchall())

# sex: what the person said, else what ingest inferred from the data
sex = prof.get("sex", "") if prof.get("sex", "") in ("male", "female") else ing.get("sex", meta.get("sex", "unknown"))
set_meta(con, "sex", sex)
log("sex", sex)
set_meta(con, "sample_id", prof.get("name", "This person"))
fmt = ing.get("format", "")
if fmt == "chip":
    kind = "%s genotyping chip" % (ing.get("vendor") or "Consumer")
    model = "chip"
elif fmt == "gvcf":
    kind = "Whole genome (gVCF)"
    model = "gvcf"
else:
    kind = "Variant calls (VCF)" if ing.get("coverage") != "genome" else "Whole genome (VCF)"
    model = "vcf"
set_meta(con, "seq_kind", kind)
set_meta(con, "call_model", model)
set_meta(con, "source_build", ing.get("build", ""))
for k, v in ing.items():                      # what the import read and dropped, for the raw-data page
    set_meta(con, "ingest_" + k, v)
set_meta(con, "built", time.strftime("%Y-%m-%d %H:%M"))

# files: what the person gave, and what the import made
cur.executescript("DROP TABLE IF EXISTS files; CREATE TABLE files(name TEXT, size INTEGER, kind TEXT, about TEXT, ord INTEGER);")
k = 0
for key in sorted(x for x in prof if x.startswith("input.")):
    p = prof[key]
    if os.path.exists(p):
        cur.execute("INSERT INTO files VALUES(?,?,?,?,?)", (os.path.basename(p), os.path.getsize(p), "Input",
                                                          "The file this profile was made from (%s)" % kind, k))
        k += 1
for p, ftype, about in [(GVCF, "gVCF", "Your calls converted to GRCh38 (made by the import)"),
                        (DBPATH, "SQLite", "This person's genome database")]:
    if os.path.exists(p):
        # the database is built as next.db and renamed when the import succeeds
        cur.execute("INSERT INTO files VALUES(?,?,?,?,?)", ("strand.db" if p == DBPATH else os.path.basename(p), os.path.getsize(p), ftype, about, k))
        k += 1

# FASTQ quality
cur.executescript("""DROP TABLE IF EXISTS qc_fastq;
CREATE TABLE qc_fastq(mate INTEGER, reads INTEGER, bases INTEGER, q30 REAL, mean_q REAL, gc REAL, n_rate REAL, adapter REAL, dup REAL,
  pos_q TEXT, gc_hist TEXT, qmean_hist TEXT, len_mode INTEGER);""")
for m in (1, 2):
    p = os.path.join(WORK, "qc", "r%d.json" % m)
    if not os.path.exists(p):
        continue
    j = json.load(open(p))
    qh = j["q_hist"]
    tot = sum(qh)
    q30 = sum(qh[30:]) / max(1, tot)
    meanq = sum(i * v for i, v in enumerate(qh)) / max(1, tot)
    gh = j["gc_hist"]
    gc = sum(i * v for i, v in enumerate(gh)) / max(1, sum(gh)) / 100
    lh = j["len_hist"]
    cur.execute("INSERT INTO qc_fastq VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (m, j["reads"], j["bases"], q30, meanq, gc, j["n_bases"] / max(1, j["bases"]), j["adapter_reads"] / max(1, j["adapter_checked"]),
                 j["dup_hits"] / max(1, j["dup_sampled"]), json.dumps(j["pos_mean_q"]), json.dumps(gh), json.dumps(j["qmean_hist"]), lh.index(max(lh))))
    log("read %d: %d reads, %.2f Gb, Q30 %.1f%%, mean Q %.1f, GC %.1f%%" % (m, j["reads"], j["bases"] / 1e9, q30 * 100, meanq, gc * 100))
con.commit()
