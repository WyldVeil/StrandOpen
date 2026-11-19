#!/usr/bin/env python3
"""Setup, once for everyone: the 1000 Genomes reference (NYGC 30x, GRCh38,
3,202 people) turned into what the import steps read. The Windows edition
runs this; on Linux pipeline/04_1kg.sh does the same work in the shell, and
the two produce the same files. Founders only (no children of trios). Per
chromosome:
  sites/chrN.bcf  every biallelic variant with its frequency overall and in
                  each of the five continental groups (rarity, strand checks)
  pgen/chrN.*     genotypes of variants with a minor allele frequency of at
                  least 0.5% (PCA, admixture, polygenic score distributions)
Chromosomes are cut into 10 Mb pieces processed in parallel, then joined.
  kg_build.py            every chromosome
  kg_build.py 21 22      just these
Paths come from the environment pipeline/setup.py sets (as for env.sh)."""
import concurrent.futures
import glob
import os
import subprocess
import sys
import threading
import time

PROFILE = os.environ["STRANDOPEN_PROFILE"]
ROOT = os.path.dirname(os.path.dirname(PROFILE))
REF = os.environ.get("STRANDOPEN_REF", os.path.join(ROOT, "shared", "ref"))
P = os.environ.get("STRANDOPEN_KG", os.path.join(ROOT, "shared", "1kg"))
RAW, KG, FA = os.path.join(REF, "raw"), os.path.join(REF, "1kg"), os.path.join(REF, "GRCh38.fa")
THREADS = int(os.environ.get("THREADS") or os.cpu_count() or 4)
KEEP = "^INFO/AF,INFO/AC,INFO/AN,INFO/AF_EUR,INFO/AF_AFR,INFO/AF_EAS,INFO/AF_SAS,INFO/AF_AMR"
NULL = subprocess.DEVNULL

for d in ("sites", "pgen", "chunks"):
    os.makedirs(os.path.join(P, d), exist_ok=True)
os.chdir(P)

# the sample lists: founders, their continental group, population and sex
founders = []
for n, line in enumerate(open(os.path.join(RAW, "20130606_g1k_3202_samples_ped_population.txt"))):
    f = line.split()
    if n > 0 and len(f) >= 7 and f[2] == "0" and f[3] == "0":      # awk 'NR>1 && $3=="0" && $4=="0"'
        founders.append(f)
with open("founders.txt", "w", newline="\n") as o:
    o.writelines("%s\n" % f[1] for f in founders)
with open("groups.txt", "w", newline="\n") as o:
    o.writelines("%s\t%s\n" % (f[1], f[6]) for f in founders)
with open("founders.pop.tsv", "w", newline="\n") as o:
    o.writelines("%s\t%s\t%s\t%s\t%s\n" % (f[1], f[1], f[5], f[6], f[4]) for f in founders)
with open("founders.sex", "w", newline="\n") as o:
    o.write("#IID\tSEX\n")
    o.writelines("%s\t%s\n" % (f[1], f[4]) for f in founders)
with open("rename.txt", "w", newline="\n") as o:
    o.writelines("chr%s %s\n" % (c, c) for c in [str(i) for i in range(1, 23)] + ["X"])


def panel(c):
    v = os.path.join(KG, "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" % c)
    if os.path.exists(v) and os.path.getsize(v) > 0:
        return v
    return os.path.join(KG, "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz" % c)


def chain(cmds, stdin=None, stdout=None, quiet=()):
    """start a pipeline of argument lists; quiet: indexes whose stderr is discarded"""
    procs = []
    for i, cmd in enumerate(cmds):
        last = i == len(cmds) - 1
        p = subprocess.Popen(cmd, stdin=procs[-1].stdout if procs else stdin,
                             stdout=stdout if last else subprocess.PIPE, stderr=NULL if i in quiet else None)
        if procs:
            procs[-1].stdout.close()      # the next process owns the read end now
        procs.append(p)
    return procs


def tee(src, dests):
    """copy one stream into several (bash's tee >(...))"""
    try:
        while True:
            b = src.read(1 << 20)
            if not b:
                break
            for d in dests:
                d.write(b)
    finally:
        for d in dests:
            try:
                d.close()
            except OSError:
                pass


def piece(c, a, b):
    o = os.path.join("chunks", "chr%s.%d" % (c, a))
    if os.path.exists(o + ".done"):
        return
    head = chain([["bcftools", "view", "-r", "chr%s:%d-%d" % (c, a, b), "-S", "founders.txt", "--force-samples",
                   "-m2", "-M2", "-v", "snps,indels", "-Ou", panel(c)],
                  ["bcftools", "view", "-i", "POS>=%d" % a, "-Ou"],
                  ["bcftools", "annotate", "--rename-chrs", "rename.txt", "-x", "INFO", "-Ou"],
                  ["bcftools", "+fill-tags", "-Ou", "--", "-S", "groups.txt", "-t", "AF,AC,AN"],
                  ["bcftools", "annotate", "-x", KEEP, "-Ou"]], stdout=subprocess.PIPE, quiet=(0, 3))
    sites = subprocess.Popen(["bcftools", "view", "-G", "-Ob", "-o", o + ".sites.bcf"], stdin=subprocess.PIPE)
    common = chain([["bcftools", "view", "-q", "0.005:minor", "-Ou"],
                    ["bcftools", "annotate", "-x", "INFO", "-Ob", "-o", o + ".common.bcf"]], stdin=subprocess.PIPE)
    t = threading.Thread(target=tee, args=(head[-1].stdout, [sites.stdin, common[0].stdin]))
    t.start()
    t.join()
    head[-1].stdout.close()
    rcs = [p.wait() for p in head + [sites] + common]
    if any(rcs):
        raise SystemExit("chr%s:%d failed (%s)" % (c, a, rcs))
    time.sleep(1)
    open(o + ".done", "w").close()


lens = {}
for line in open(FA + ".fai"):
    f = line.split("\t")
    lens[f[0]] = int(f[1])
CS = sys.argv[1:] or [str(i) for i in range(1, 23)] + ["X"]


def finished(c):
    pg, ix = os.path.join("pgen", "chr%s.pgen" % c), os.path.join("sites", "chr%s.bcf.csi" % c)
    return all(os.path.exists(x) and os.path.getsize(x) > 0 for x in (pg, ix))


jobs = [(c, a, a + 9999999) for c in CS if not finished(c) for a in range(1, lens[c] + 1, 10000000)]
with concurrent.futures.ThreadPoolExecutor(int(os.environ.get("KGJOBS") or THREADS)) as ex:
    for f in [ex.submit(piece, *j) for j in jobs]:
        f.result()


def start_of(path):   # chunks/chr22.10000001.sites.bcf -> 10000001 (sort -t. -k2 -n)
    return int(os.path.basename(path).split(".")[1])


for c in CS:
    if finished(c):
        continue
    for kind in ("sites", "common"):
        parts = sorted(glob.glob(os.path.join("chunks", "chr%s.*.%s.bcf" % (c, kind))), key=start_of)
        with open(os.path.join("chunks", "chr%s.%s.list" % (c, kind)), "w", newline="\n") as o:
            o.writelines(p.replace(os.sep, "/") + "\n" for p in parts)
    # --naive: the pieces share one header, so their blocks are joined as they are
    sbcf = "sites/chr%s.bcf" % c
    subprocess.run(["bcftools", "concat", "--naive", "-f", "chunks/chr%s.sites.list" % c, "-Ob", "-o", sbcf], check=True)
    subprocess.run(["bcftools", "index", "-f", sbcf], check=True)
    cbcf = "pgen/chr%s.common.bcf" % c
    subprocess.run(["bcftools", "concat", "--naive", "-f", "chunks/chr%s.common.list" % c, "-Ob", "-o", cbcf], check=True)
    sex = ["--update-sex", "founders.sex", "--split-par", "hg38"] if c == "X" else []
    with open("pgen/chr%s.log" % c, "w") as log:
        subprocess.run(["plink2", "--bcf", cbcf, "--keep", "founders.txt", "--max-alleles", "2",
                        "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "60", "missing"] + sex +
                       ["--threads", "4", "--memory", "4000", "--make-pgen", "--out", "pgen/chr%s" % c],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    os.remove(cbcf)
    for x in glob.glob(os.path.join("chunks", "chr%s.*" % c)):
        os.remove(x)
    print("chr%s done" % c, flush=True)

with open(os.path.join("pgen", "merge.txt"), "w", newline="\n") as o:
    o.writelines("pgen/%s\n" % os.path.basename(x)[:-5] for x in sorted(glob.glob(os.path.join("pgen", "*.pgen"))))
with open("sites.list", "w", newline="\n") as o:
    o.writelines("sites/%s\n" % os.path.basename(x) for x in sorted(glob.glob(os.path.join("sites", "*.bcf"))))
print("1000 Genomes reference ready")
