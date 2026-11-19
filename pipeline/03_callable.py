#!/usr/bin/env python3
"""Step 3: what the sequencing could see.

The gVCF holds reference blocks (stretches called as matching the reference,
with their minimum depth) between the variant records. Those blocks are
what turns "no variant here" into "reference genotype here"; a position in
no block and with no variant was not called.

Tables: refblock (callable reference stretches), covbin (10 kb bins: mean
depth, callable share, variant counts), meta keys for genome-wide totals."""
import numpy as np
from strandlib import *

BIN = 10000
con = db()
cur = con.cursor()
# what kind of data this is decides how later steps read it (chip: no depth)
_ing = {}
_p = os.path.join(CANON, "ingest.kv")
if os.path.exists(_p):
    for _l in open(_p):
        if "=" in _l:
            _k, _v = _l.rstrip("\n").split("=", 1)
            _ing[_k] = _v
set_meta(con, "call_model", {"chip": "chip", "gvcf": "gvcf"}.get(_ing.get("format", ""), "vcf"))
set_meta(con, "source_format", _ing.get("format", ""))
set_meta(con, "source_vendor", _ing.get("vendor", ""))
set_meta(con, "source_build", _ing.get("build", ""))
cur.executescript("""
DROP TABLE IF EXISTS refblock; DROP TABLE IF EXISTS covbin;
CREATE TABLE refblock(chrom INTEGER, start INTEGER, end INTEGER, dp INTEGER, pass INTEGER);
CREATE TABLE covbin(chrom INTEGER, bin INTEGER, depth REAL, callable REAL, nhet INTEGER, nhom INTEGER,
  nhemi INTEGER, nlow INTEGER, PRIMARY KEY(chrom, bin)) WITHOUT ROWID;
""")
lens = chrom_lengths()
nb = {c: lens[c] // BIN + 1 for c in CHROMS}
dpsum = {c: np.zeros(nb[c]) for c in CHROMS}
callsum = {c: np.zeros(nb[c]) for c in CHROMS}
vc = {c: np.zeros((nb[c], 4), dtype=np.int64) for c in CHROMS}


def add_span(c, s, e, dp):
    """bases s..e (1-based inclusive) at depth dp into the bins"""
    a, b = s - 1, e
    while a < b:
        bi = a // BIN
        hi = min(b, (bi + 1) * BIN)
        n = hi - a
        dpsum[c][bi] += n * dp
        callsum[c][bi] += n
        a = hi


rows = []
nblock = 0
called = 0
called_pass = 0
for line in opentext(GVCF):
    if line[0] == "#":
        continue
    p = line.split("\t", 10)
    c = p[0]
    if c not in CHROM_IDX:
        continue
    pos = int(p[1])
    vals = p[9].rstrip("\n").split(":")
    gt = vals[0]
    fmt = p[8].split(":")
    if p[4] == ".":
        if gt not in ("0/0", "0", "0|0"):
            continue
        end = pos
        info = p[7]
        i = info.find("END=")
        if i >= 0:
            j = info.find(";", i)
            end = int(info[i + 4:j if j >= 0 else None])
        dp = 0
        i = info.find("MinDP=")
        if i >= 0:
            j = info.find(";", i)
            dp = int(info[i + 6:j if j >= 0 else None])
        elif "DP" in fmt:
            try:
                dp = int(vals[fmt.index("DP")])
            except ValueError:
                dp = 0
        ok = 1 if p[6] == "PASS" else 0
        rows.append((CHROM_IDX[c], pos, end, dp, ok))
        add_span(c, pos, end, dp)
        called += end - pos + 1
        if ok:
            called_pass += end - pos + 1
        nblock += 1
        if len(rows) >= 200000:
            cur.executemany("INSERT INTO refblock VALUES(?,?,?,?,?)", rows)
            rows.clear()
    else:
        if gt in ("./.", ".", "0/0", "0"):
            continue
        dp = 0
        if "DP" in fmt:
            try:
                dp = int(vals[fmt.index("DP")])
            except ValueError:
                dp = 0
        rl = len(p[3])
        add_span(c, pos, pos + rl - 1, dp)
        called += rl
        bi = (pos - 1) // BIN
        if p[6] != "PASS":
            vc[c][bi, 3] += 1
        if "/" in gt or "|" in gt:
            a = gt.replace("|", "/").split("/")
            vc[c][bi, 0 if a[0] != a[1] else 1] += 1
        else:
            vc[c][bi, 2] += 1
if rows:
    cur.executemany("INSERT INTO refblock VALUES(?,?,?,?,?)", rows)
log("reference blocks", nblock, "called bases", called)

brow = []
for c in CHROMS:
    cs = callsum[c]
    with np.errstate(invalid="ignore", divide="ignore"):
        dep = np.where(cs > 0, dpsum[c] / np.maximum(cs, 1), 0.0)
    for bi in range(nb[c]):
        span = min(BIN, lens[c] - bi * BIN)
        brow.append((CHROM_IDX[c], bi, float(dep[bi]), float(cs[bi] / span) if span > 0 else 0.0,
                     int(vc[c][bi, 0]), int(vc[c][bi, 1]), int(vc[c][bi, 2]), int(vc[c][bi, 3])))
cur.executemany("INSERT INTO covbin VALUES(?,?,?,?,?,?,?,?)", brow)
cur.execute("CREATE INDEX refblock_pos ON refblock(chrom, start)")

auto = [c for c in CHROMS if c not in ("X", "Y", "MT")]
tot = sum(lens[c] for c in auto)
cal = sum(callsum[c].sum() for c in auto)
mean_dp = sum(dpsum[c].sum() for c in auto) / max(cal, 1)
set_meta(con, "callable_bases", int(called))
set_meta(con, "callable_pass_bases", int(called_pass))
set_meta(con, "autosome_callable_fraction", "%.4f" % (cal / tot))
set_meta(con, "autosome_mean_depth", "%.2f" % mean_dp)
set_meta(con, "covbin_size", BIN)
for c in ("X", "Y"):
    set_meta(con, "callable_fraction_" + c, "%.4f" % (callsum[c].sum() / lens[c]))
    set_meta(con, "mean_depth_" + c, "%.2f" % (dpsum[c].sum() / max(callsum[c].sum(), 1)))
con.commit()
log("autosomes callable %.3f, mean depth %.1f" % (cal / tot, mean_dp))
