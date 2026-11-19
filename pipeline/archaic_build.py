#!/usr/bin/env python3
"""Setup, once for everyone: the African outgroup for the archaic-DNA step.
The Windows edition runs this; on Linux pipeline/setup_archaic.sh does the
same work in the shell, and the two produce the same files. Variants
polymorphic among the 1000 Genomes founders from five African populations
(YRI, ESN, MSL, LWK, GWD), in hmmix's outgroup format, one file per
chromosome in <1kg>/archaic. Needs the phased panel (<ref>/1kg) and its
sample sheet. Paths come from the environment pipeline/setup.py sets.
  archaic_build.py            chromosomes 1 to 22
  archaic_build.py 21 22      just these"""
import concurrent.futures
import os
import subprocess
import sys

PROFILE = os.environ["STRANDOPEN_PROFILE"]
ROOT = os.path.dirname(os.path.dirname(PROFILE))
REF = os.environ.get("STRANDOPEN_REF", os.path.join(ROOT, "shared", "ref"))
KGOUT = os.environ.get("STRANDOPEN_KG", os.path.join(ROOT, "shared", "1kg"))
RAW, KG = os.path.join(REF, "raw"), os.path.join(REF, "1kg")
OUT = os.path.join(KGOUT, "archaic")
os.makedirs(OUT, exist_ok=True)
os.chdir(OUT)

AFRICAN = {"YRI", "ESN", "MSL", "LWK", "GWD"}
with open("outgroup_samples.txt", "w", newline="\n") as o:
    for n, line in enumerate(open(os.path.join(RAW, "20130606_g1k_3202_samples_ped_population.txt"))):
        f = line.split()
        if n > 0 and len(f) >= 7 and f[2] == "0" and f[3] == "0" and f[5] in AFRICAN:
            o.write(f[1] + "\n")
print("outgroup founders: %d" % sum(1 for _ in open("outgroup_samples.txt")), flush=True)


def panel(c):
    v = os.path.join(KG, "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" % c)
    if os.path.exists(v) and os.path.getsize(v) > 0:
        return v
    return os.path.join(KG, "1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz" % c)


def one(c):
    if os.path.exists("og_chr%s.tsv" % c) and os.path.getsize("og_chr%s.tsv" % c) > 0:
        return
    cmds = [["bcftools", "view", "-S", "outgroup_samples.txt", "--force-samples", "-v", "snps", "-Ou", panel(c)],
            ["bcftools", "norm", "-m", "-any", "-Ou"],
            ["bcftools", "view", "-v", "snps", "-Ou"],
            ["bcftools", "+fill-tags", "-Ou", "--", "-t", "AC,AN"],
            ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n"]]
    procs = []
    for i, cmd in enumerate(cmds):
        p = subprocess.Popen(cmd, stdin=procs[-1].stdout if procs else None, stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL if i < 2 else None)
        if procs:
            procs[-1].stdout.close()
        procs.append(p)
    tmp = "og_chr%s.tmp" % c
    with open(tmp, "w", newline="\n") as o:
        for line in procs[-1].stdout:
            f = line.decode().rstrip("\n").split("\t")
            try:
                ac, an = int(f[4]), int(f[5])
            except (ValueError, IndexError):
                continue
            if 0 < ac < an:          # polymorphic in the outgroup
                anc = f[2] if an - ac >= ac else f[3]
                o.write("%s\t%s\t%s:%d\t%s:%d\t%s\n" % (f[0], f[1], f[2], an - ac, f[3], ac, anc))
    procs[-1].stdout.close()
    rcs = [p.wait() for p in procs]
    if any(rcs):
        raise SystemExit("chr%s failed (%s)" % (c, rcs))
    os.replace(tmp, "og_chr%s.tsv" % c)
    print("chr%s done" % c, flush=True)


with concurrent.futures.ThreadPoolExecutor(int(os.environ.get("JOBS") or 8)) as ex:
    for f in [ex.submit(one, c) for c in (sys.argv[1:] or [str(i) for i in range(1, 23)])]:
        f.result()
print("archaic outgroup ready")
