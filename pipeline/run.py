#!/usr/bin/env python3
"""Strand Open's import: from one person's raw file to their genome database.

    run.py --profile PROFILE_DIR [--steps id,id,...]

The profile folder holds profile.kv (name, sex, input.1=PATH, ...) and gets
canon/, work/, db/ and run/. Reference data is shared between people and
only read: <root>/shared.kv names where it lives (ref=, kg=), so one machine
can reuse folders that already exist.

Progress for the app: run/progress.kv (state, step, title, done, of,
skipped) is replaced after every step; run/log.txt keeps every step's
output. The database is built as db/next.db and moved over db/strand.db
only when every essential step succeeded, so a failed import never spoils
a working profile.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path, encoding="utf-8", errors="replace"):
            if "=" in line and not line.startswith("#"):
                k, v = line.rstrip("\n").split("=", 1)
                out[k.strip()] = v.strip()
    return out


# id, script, title, needs, essential. needs: "*" any data; "seq" sequencing
# (a VCF or gVCF, not a chip); "genome" the whole genome called (a gVCF, or a
# whole-genome VCF whose unlisted positions are taken as reference); "kg" the
# 1000 Genomes reference; "sv" a structural-variant file
STEPS = [
    ("ingest", None, "Reading your file", "*", True),
    ("01", "01_variants.sh", "Normalising variants, consequences, ClinVar", "*", True),
    ("02", "02_genes.py", "Genes, transcripts and conditions", "*", True),
    ("03", "03_callable.py", "Which positions were tested", "*", True),
    ("05", "05_variants.py", "The variant table and frequencies", "*", True),
    ("05b", "05b_mapq.py", "Call quality and rsIDs", "seq", False),
    ("05c", "05c_chip.py", "Setting aside rare chip calls", "chip", True),
    ("sex", None, "Chromosomal sex", "*", True),
    ("06", "06_haplogroups.py", "Maternal and paternal lines", "*", False),
    ("07", "07_pgx.py", "Medicines (PharmCAT)", "*", False),
    ("08", "08_traits.py", "Traits and well-known variants", "*", False),
    ("11", "11_annot.py", "ACMG genes, GWAS Catalog, PharmGKB", "*", False),
    ("12", "12_meta.py", "Sample summary", "*", True),
    ("13", "13_ancestry.py", "Ancestry", "kg", False),
    ("15", "15_prs.py", "Polygenic scores", "kg", False),
    ("16", "16_rare.py", "Rare and remarkable variants", "*", False),
    ("18", "18_rarity.py", "How rare your variants are", "seq", False),
    ("19", "19_lineage.py", "Where your lineages come from", "*", False),
    ("20", "20_roh.py", "Parents' relatedness", "*", False),
    ("21", "21_archaic.py", "Neanderthal DNA", "genome", False),
    ("17", "17_highlights.py", "Findings to look at first", "*", False),
    ("23", "23_diet.py", "Diet and supplements", "*", False),
    ("22", "22_brief.py", "Your genome in brief", "*", False),
]
WHY = {
    "seq": "needs sequencing data: rare chip calls are too often wrong to report",
    "genome": "needs the whole genome called (sequencing, not a chip or an exome)",
    "sv": "needs a structural-variant file from whole-genome sequencing",
    "kg": "needs the 1000 Genomes reference (Setup)",
    "reads": "needs sequencing reads",
}


class Progress:
    def __init__(self, pdir, of):
        self.dir = os.path.join(pdir, "run")
        os.makedirs(self.dir, exist_ok=True)
        self.st = {"state": "running", "step": "", "title": "", "done": "0", "of": str(of), "skipped": "",
                   "started": time.strftime("%Y-%m-%d %H:%M:%S"), "pid": str(os.getpid())}
        self.write()

    def write(self):
        tmp = os.path.join(self.dir, "progress.kv.tmp")
        with open(tmp, "w") as f:
            for k, v in self.st.items():
                f.write("%s=%s\n" % (k, str(v).replace("\n", " ")))
        os.replace(tmp, os.path.join(self.dir, "progress.kv"))

    def set(self, **kw):
        self.st.update({k: str(v) for k, v in kw.items()})
        self.write()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--steps", default="")
    a = ap.parse_args()
    pdir = os.path.abspath(a.profile)
    root = os.path.dirname(os.path.dirname(pdir))
    prof = kv(os.path.join(pdir, "profile.kv"))
    shared = kv(os.path.join(root, "shared.kv"))
    inputs = [prof[k] for k in sorted(prof) if k.startswith("input.")]
    if not inputs:
        raise SystemExit("the profile names no input file")
    env = dict(os.environ)
    env["STRANDOPEN_PROFILE"] = pdir
    env["STRANDOPEN_REF"] = shared.get("ref", os.path.join(root, "shared", "ref"))
    env["STRANDOPEN_KG"] = shared.get("kg", os.path.join(root, "shared", "1kg"))
    env["STRANDOPEN_SHARED"] = os.path.join(root, "shared")
    db_next = os.path.join(pdir, "db", "next.db")
    env["STRANDOPEN_DBPATH"] = db_next
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if shared.get("archaic"):
        env["STRANDOPEN_ARCHAIC"] = shared["archaic"]
    # the bundled tools first, then any extra tool folders shared.kv names (path=a:b)
    extra = [d for d in shared.get("path", "").split(os.pathsep) if d]
    for d in reversed([os.path.join(APP, "runtime", "bin")] + extra):
        if os.path.isdir(d):
            env["PATH"] = d + os.pathsep + env.get("PATH", "")
    os.makedirs(os.path.join(pdir, "db"), exist_ok=True)
    want = set(x for x in a.steps.split(",") if x)
    final = os.path.join(pdir, "db", "strand.db")
    if not want:                                      # a full run starts clean
        if os.path.exists(db_next):
            os.remove(db_next)
    elif not os.path.exists(db_next) and os.path.exists(final):
        shutil.copyfile(final, db_next)               # chosen steps rework a copy of the finished database
    steps = [s for s in STEPS if not want or s[0] in want]
    prog = Progress(pdir, len(steps))
    log = open(os.path.join(pdir, "run", "log.txt"), "a")
    log.write("\n==== import started %s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
    kind = None
    ing = {}
    skipped, failed = [], []
    has_kg = os.path.isdir(os.path.join(env["STRANDOPEN_KG"], "pgen"))
    for i, (sid, script, title, needs, essential) in enumerate(steps):
        prog.set(step=sid, title=title, done=i)
        if kind is None and sid != "ingest":
            ing = kv(os.path.join(pdir, "canon", "ingest.kv"))
            kind = ing.get("format", "")
        # does this kind of data support the step?
        if needs != "*" and sid != "ingest":
            ok = {"kg": has_kg, "seq": kind in ("vcf", "gvcf"), "chip": kind == "chip",
                  "genome": kind == "gvcf" or (kind == "vcf" and ing.get("coverage") == "genome"),
                  "sv": os.path.exists(os.path.join(pdir, "canon", "sv.vcf.gz"))}.get(needs, False)
            if not ok:
                if needs != "chip":            # a chip-only step is not something other data lacks
                    skipped.append("%s: %s" % (title, WHY.get(needs, "not available for this data")))
                log.write("-- skip %s (%s)\n" % (sid, WHY.get(needs, "")))
                continue
        t0 = time.time()
        if sid == "ingest":
            cmd = [sys.executable, "-B", os.path.join(HERE, "ingest.py"), "--profile", pdir, inputs[0]]
            if prof.get("sample"):
                cmd[4:4] = ["--sample", prof["sample"]]
        elif sid == "sex":
            # the chromosomal sex every later step depends on, before any of them reads chromosome X
            ing = kv(os.path.join(pdir, "canon", "ingest.kv"))
            sex = prof.get("sex") if prof.get("sex") in ("male", "female") else ing.get("sex", "")
            if not sex or sex == "unknown":
                r = subprocess.run("bcftools view -H -r X:2781480-155701382 -i 'GT=\"het\"' '%s' | wc -l; bcftools view -H -r X:2781480-155701382 -i 'GT=\"alt\"' '%s' | wc -l"
                                   % ((os.path.join(pdir, "canon", "calls.g.vcf.gz"),) * 2), shell=True, capture_output=True, text=True)
                try:
                    het, alt = (int(x) for x in r.stdout.split())
                    sex = "male" if het < 0.05 * max(1, alt) else "female"
                except ValueError:
                    sex = "unknown"
            cmd = [sys.executable, "-B", "-c", "import sys; sys.path.insert(0, %r); from strandlib import *; c=db(); set_meta(c, 'sex', %r); c.commit()" % (HERE, sex)]
        else:
            cmd = ["bash", os.path.join(HERE, script)] if script.endswith(".sh") else [sys.executable, "-B", os.path.join(HERE, script)]
        log.write("-- %s %s\n" % (sid, title))
        log.flush()
        # each step runs in the person's work folder: the program folder may be read-only, and
        # tools such as PharmCAT write logs where they start
        work = os.path.join(pdir, "work")
        os.makedirs(work, exist_ok=True)
        r = subprocess.run(cmd, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
        log.write("-- %s finished in %.0f s with code %d\n" % (sid, time.time() - t0, r.returncode))
        log.flush()
        if r.returncode != 0:
            if essential:
                prog.set(state="failed", step=sid, title=title, error=last_error(pdir))
                raise SystemExit(1)
            failed.append(title)
    # the new database replaces the old one only now
    if os.path.exists(db_next):
        subprocess.run([sys.executable, "-c", "import sqlite3; c=sqlite3.connect(%r); c.execute('ANALYZE'); c.commit()" % db_next])
        os.replace(db_next, final)
    prog.set(state="done_with_problems" if failed else "done", done=len(steps), step="", title="",
             skipped=" | ".join(skipped), failed=" | ".join(failed), finished=time.strftime("%Y-%m-%d %H:%M:%S"))


def last_error(pdir):
    """the last meaningful line of the log, for the app to show"""
    try:
        lines = open(os.path.join(pdir, "run", "log.txt"), encoding="utf-8", errors="replace").read().splitlines()
    except OSError:
        return ""
    for line in reversed(lines):
        if line.strip() and not line.startswith("--"):
            return line.strip()[-300:]
    return ""


if __name__ == "__main__":
    main()
