"""Shared helpers for the Strand Open import steps (standard library + numpy)."""
import gzip
import os
import sqlite3
import subprocess
import sys
import time

# One person's data lives in a profile folder; reference data is shared by
# everyone and only read here.
#   STRANDOPEN_PROFILE  the person's folder: canon/ (their calls, from
#                       ingest.py), work/ (step outputs), db/strand.db
#   STRANDOPEN_REF      reference downloads (GRCh38.fa, raw/ annotations)
#   STRANDOPEN_KG       1000 Genomes products (sites/, pgen/, founders lists)
PROFILE = os.environ.get("STRANDOPEN_PROFILE", "")
if not PROFILE:
    raise SystemExit("STRANDOPEN_PROFILE is not set: run the import through pipeline/run.py")
DATA = PROFILE
REF = os.environ.get("STRANDOPEN_REF", os.path.join(os.path.dirname(os.path.dirname(PROFILE)), "shared", "ref"))
KG = os.environ.get("STRANDOPEN_KG", os.path.join(os.path.dirname(os.path.dirname(PROFILE)), "shared", "1kg"))
RAW = os.path.join(REF, "raw")
WORK = os.path.join(DATA, "work")
DBDIR = os.path.join(DATA, "db")
DBPATH = os.environ.get("STRANDOPEN_DBPATH", os.path.join(DBDIR, "strand.db"))
FA = os.path.join(REF, "GRCh38.fa")
CANON = os.path.join(DATA, "canon")
GVCF = os.path.join(CANON, "calls.g.vcf.gz")
CRAM = os.path.join(DATA, "align", "reads.cram")
os.makedirs(WORK, exist_ok=True)

NPROC = os.cpu_count() or 4
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
CHROM_IDX = {c: i + 1 for i, c in enumerate(CHROMS)}
CHROM_LEN = {}


def chrom_lengths():
    if not CHROM_LEN:
        with open(FA + ".fai") as f:
            for line in f:
                p = line.split("\t")
                if p[0] in CHROM_IDX:
                    CHROM_LEN[p[0]] = int(p[1])
    return CHROM_LEN


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, file=sys.stderr, flush=True)


def opentext(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def db():
    os.makedirs(DBDIR, exist_ok=True)
    con = sqlite3.connect(DBPATH)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA temp_store=MEMORY")
    con.execute("PRAGMA cache_size=-1000000")
    con.execute("CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT)")
    return con


def set_meta(con, key, value):
    con.execute("INSERT OR REPLACE INTO meta VALUES(?,?)", (key, str(value)))


def pipe_lines(cmd):
    """Run a shell pipeline and yield its stdout lines (text)."""
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
                         stderr=subprocess.DEVNULL)
    assert p.stdout is not None
    for line in p.stdout:
        yield line
    p.wait()
    if p.returncode not in (0, None):
        raise RuntimeError("command failed (%d): %s" % (p.returncode, cmd))


def unquote_gff(s):
    return (s.replace("%3B", ";").replace("%3D", "=").replace("%2C", ",")
             .replace("%26", "&").replace("%25", "%"))


# ---- the 1000 Genomes site files (step 4: one per chromosome) ----
KG_SITES = os.path.join(KG, "sites")


def kg_site_files():
    """{chrom: path} for the per-chromosome site files that exist"""
    out = {}
    for c in CHROMS:
        p = os.path.join(KG_SITES, "chr%s.bcf" % c)
        if os.path.exists(p) and os.path.exists(p + ".csi"):
            out[c] = p
    return out


def kg_query(positions, fmt):
    """bcftools query on the site files at the given (chrom, pos) list;
    yields output lines. One query per chromosome, run in parallel."""
    import concurrent.futures
    files = kg_site_files()
    by = {}
    for c, p in positions:
        if c in files:
            by.setdefault(c, set()).add(int(p))

    def one(c):
        rf = os.path.join(WORK, "kgq_%s_%d.tsv" % (c, os.getpid()))
        with open(rf, "w") as f:
            for p in sorted(by[c]):
                f.write("%s\t%d\n" % (c, p))
        r = subprocess.run(["bcftools", "query", "-R", rf, "-f", fmt, files[c]], capture_output=True, text=True)
        os.remove(rf)
        return r.stdout.splitlines()
    with concurrent.futures.ThreadPoolExecutor(8) as ex:
        for lines in ex.map(one, sorted(by)):
            yield from lines
