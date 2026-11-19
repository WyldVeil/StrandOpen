#!/usr/bin/env python3
"""Step 21: DNA inherited from archaic humans (Neanderthals, Denisovans).

People outside Africa carry stretches of DNA from Neanderthals (about 1 to 2
percent of the genome) and, in some populations, Denisovans. Such stretches
are recognisable without an archaic genome: they hold clusters of variants
that are absent from sub-Saharan Africans, whose ancestors did not mix with
Neanderthals in the same way. hmmix (Skov et al. 2018, PLoS Genetics) fits a
two-state hidden Markov model to the density of those variants in 1 kb
windows and decodes each stretch as human or archaic.

Inputs: the outgroup (variants polymorphic among 529 African 1000 Genomes
founders: YRI, ESN, MSL, LWK, GWD; built once for everyone by
pipeline/setup_archaic.sh into the shared 1000 Genomes folder), the
person's called variants, and where their genome was called as weights.

Table archaic(chrom, start, end, len, prob, nvar); meta archaic_*.
"""
import glob
import json
import shutil
import subprocess
import sys
from strandlib import *

con = db()
cur = con.cursor()
A = os.path.join(WORK, "archaic")
os.makedirs(A, exist_ok=True)
OG = os.environ.get("STRANDOPEN_ARCHAIC", os.path.join(KG, "archaic"))
# hmmix: its command on PATH (the Linux bundle's launcher), or else its module in this
# Python (the Windows bundle installs it there, with no command of its own)
HM = [shutil.which("hmmix")] if shutil.which("hmmix") else []
if not HM:
    import importlib.util
    if importlib.util.find_spec("hmm_functions") and importlib.util.find_spec("main"):
        HM = [sys.executable, "-c", "import sys; from main import main; sys.argv[0] = 'hmmix'; sys.exit(main())"]
og_parts = sorted(glob.glob(os.path.join(OG, "og_chr*.tsv")))
if not HM or len(og_parts) < 22:
    log("hmmix or the African outgroup is missing (Setup builds it): skipped")
    raise SystemExit(0)
sample = "SAMPLE"


def run(args, **kw):
    log("hmmix", args[0])
    r = subprocess.run(HM + args, cwd=A, capture_output=True, text=True, **kw)
    if r.returncode:
        raise SystemExit("hmmix %s failed: %s" % (args[0], r.stderr[-600:]))


# ---- outgroup file in hmmix's format (sorted by chromosome then position)
og = os.path.join(A, "outgroup.txt")      # each person's own copy: the shared folder may be read-only
if not os.path.exists(og) or os.path.getmtime(og) < max(os.path.getmtime(p) for p in og_parts):
    with open(og, "w") as out:
        out.write("chrom\tpos\tref_allele_info\talt_allele_info\tancestral_base\n")
        for c in range(1, 23):
            with open(os.path.join(OG, "og_chr%d.tsv" % c)) as f:
                rows = sorted((int(l.split("\t", 2)[1]), l) for l in f)
            out.writelines(l for _, l in rows)

# ---- weights: where your genome was called (gVCF reference blocks that passed, and variant sites)
wb = os.path.join(A, "weights.bed")
iv = {}
for c, a, b in cur.execute("SELECT chrom, start, end FROM refblock WHERE pass=1 AND chrom<=22"):
    iv.setdefault(c, []).append((a - 1, b))
for c, a, b in cur.execute("SELECT chrom, pos, end FROM variant WHERE pass=1 AND chrom<=22"):
    iv.setdefault(c, []).append((a - 1, b))
# ... inside the 1000 Genomes strict accessibility mask (where short reads map
# reliably): elsewhere mis-mapped reads pile up variants the African panel
# never recorded, and the model mistakes those clusters for archaic DNA
mask = {}
mf = os.path.join(RAW, "1kg_mask", "strict.bed")
if os.path.exists(mf):
    for line in open(mf):
        p = line.split("\t")
        if p[0].startswith("chr") and p[0][3:].isdigit():
            mask.setdefault(int(p[0][3:]), []).append((int(p[1]), int(p[2])))
else:
    log("no strict accessibility mask in ref/raw/1kg_mask: results will be noisier")


def merge(ivs):
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def intersect(x, y):
    i = j = 0
    out = []
    while i < len(x) and j < len(y):
        a, b = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
        if a < b:
            out.append((a, b))
        if x[i][1] < y[j][1]:
            i += 1
        else:
            j += 1
    return out


callable_bp = 0
with open(wb, "w") as f:
    for c in sorted(iv):
        merged = merge(iv[c])
        if mask:
            merged = intersect(merged, merge(mask.get(c, [])))
        for a, b in merged:
            f.write("chr%d\t%d\t%d\n" % (c, a, b))
            callable_bp += b - a
log("callable autosomal bases", callable_bp)

# ---- your variants, autosomes, chr names (as in the outgroup)
vcf = os.path.join(A, "ingroup.vcf.gz")
ren = os.path.join(A, "chr.txt")
with open(ren, "w") as f:
    f.write("".join("%d chr%d\n" % (c, c) for c in range(1, 23)))
src = os.path.join(WORK, "variants.annot.vcf.gz")
name = os.path.join(A, "ingroup_name.txt")    # the sample's name in the outgroup's terms (a file, not <(...): no bash needed)
with open(name, "w", newline="\n") as f:
    f.write(sample + "\n")
subprocess.run("bcftools view -f PASS -v snps -r %s '%s' -Ou | bcftools annotate --rename-chrs '%s' -Ou | bcftools reheader -s '%s' "
               "| bcftools view -Oz -o '%s' && bcftools index -f '%s'" % (",".join(str(c) for c in range(1, 23)), src, ren, name, vcf, vcf),
               shell=True, check=True, executable="/bin/bash", capture_output=True)
json.dump({"ingroup": [sample], "outgroup": []}, open(os.path.join(A, "ingroup.json"), "w"))

run(["mutation_rate", "-outgroup=outgroup.txt", "-weights=weights.bed", "-window_size=1000000", "-out=mutationrate.bed"])
# your observations, as hmmix's create_ingroup writes them without an ancestral
# genome: heterozygous SNVs inside the callable regions and absent from the
# outgroup (written here with bcftools; hmmix's own route needs vcftools and
# a list of several VCF files)
pos = os.path.join(A, "og_pos.tsv")
if not os.path.exists(pos) or os.path.getmtime(pos) < os.path.getmtime(og):
    subprocess.run("tail -n +2 '%s' | cut -f1,2 > '%s'" % (og, pos), shell=True, check=True)
obs = "obs.%s.txt" % sample
# quality screens, each only where the file carries the field: a provider's
# VCF may have no QUAL, mapping quality or depth (1000 Genomes calls have none)


def carried(expr):
    out = subprocess.run("bcftools query -f '%s\\n' '%s' 2>/dev/null | head -2000" % (expr, vcf), shell=True,
                         capture_output=True, text=True, executable="/bin/bash").stdout.split()
    return sum(1 for x in out if x not in (".", "")) > len(out) / 2


screens = [(e, f) for e, f in (("%QUAL", "QUAL>=30"), ("%INFO/MQ", "INFO/MQ>=50"), ("%INFO/MQ0F", "INFO/MQ0F<0.05"),
                                ("[%GQ]", "FMT/GQ>=30"), ("[%DP]", "FMT/DP>=10 && FMT/DP<=80")) if carried(e)]
cond = " && ".join(f for _, f in screens)
log("observation screens:", cond or "none (the file has no quality fields)")
set_meta(con, "archaic_screens", cond or "none")
incl = ("-i '" + cond + "' ") if cond else ""
subprocess.run("(printf 'chrom\\tpos\\tancestral_base\\tgenotype\\n'; "
               "bcftools view -v snps -g het " + incl + "-T '%s' '%s' -Ou "
               "| bcftools view -T '^%s' -Ou "
               "| bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t[%%TGT]\\n' | sed 's#[/|]##') > '%s'"
               % (wb, vcf, pos, os.path.join(A, obs)), shell=True, check=True, executable="/bin/bash")
nobs = sum(1 for _ in open(os.path.join(A, obs))) - 1
log("observations", nobs)
# The outgroup is sub-Saharan African, so for someone of mainly sub-Saharan
# African ancestry the model compares them with their own populations: their
# variants are mostly in the outgroup and the result means nothing. Then, or
# with too few variants to fit the model, the step says so instead.
main_grp = ""
try:
    r = cur.execute("SELECT pop FROM anc_admix WHERE level='super' ORDER BY frac DESC LIMIT 1").fetchone()
    main_grp = r[0] if r else ""
except Exception:
    pass
why = ""
if main_grp == "AFR":
    why = ("This cannot be measured for you: the method finds archaic stretches by comparing a genome with sub-Saharan African ones, and your "
           "ancestry is mainly sub-Saharan African, so it would be comparing you with your own populations. Africans carry little archaic DNA "
           "from Neanderthals (some came back with migrations into Africa, and other archaic groups contributed too), but measuring it needs "
           "a different reference.")
elif nobs < 5000:
    why = "Too few usable variants in this file to fit the model (%d); a whole-genome file with quality fields gives enough." % nobs
if why:
    cur.executescript("DROP TABLE IF EXISTS archaic; CREATE TABLE archaic(chrom INTEGER, start INTEGER, end INTEGER, len INTEGER, prob REAL, nvar INTEGER);")
    for k in ("archaic_segments", "archaic_bp", "archaic_frac", "archaic_longest"):
        set_meta(con, k, "")
    set_meta(con, "archaic_unmeasurable", why)
    con.commit()
    log("archaic DNA not measured:", why[:80])
    raise SystemExit(0)
set_meta(con, "archaic_unmeasurable", "")
run(["train", "-obs=%s" % obs, "-weights=weights.bed", "-mutrates=mutationrate.bed", "-out=trained.json"])
run(["decode", "-obs=%s" % obs, "-weights=weights.bed", "-mutrates=mutationrate.bed", "-param=trained.json", "-out=decoded"])

# ---- the archaic stretches
dec = sorted(glob.glob(os.path.join(A, "decoded*.txt")), key=os.path.getmtime)[-1]
cur.executescript("DROP TABLE IF EXISTS archaic; CREATE TABLE archaic(chrom INTEGER, start INTEGER, end INTEGER, len INTEGER, prob REAL, nvar INTEGER);")
rows = []
with open(dec) as f:
    hdr = f.readline().split()
    ix = {h: i for i, h in enumerate(hdr)}
    for line in f:
        p = line.split()
        if p[ix["state"]] != "Archaic":
            continue
        c = int(p[ix["chrom"]].replace("chr", ""))
        a, b = int(p[ix["start"]]), int(p[ix["end"]])
        prob = float(p[ix["mean_prob"]]) if "mean_prob" in ix else 1.0
        nv = int(p[ix["snps"]]) if "snps" in ix else 0
        rows.append((c, a, b, b - a, prob, nv))
cur.executemany("INSERT INTO archaic VALUES(?,?,?,?,?,?)", rows)
confident = [r for r in rows if r[4] >= 0.8 and r[3] >= 40000]
total = sum(r[3] for r in confident)
frac = total / max(1, callable_bp)
for k, v in (("archaic_segments", len(confident)), ("archaic_bp", total), ("archaic_callable_bp", callable_bp), ("archaic_frac", "%.5f" % frac),
             ("archaic_longest", max((r[3] for r in confident), default=0)),
             ("archaic_method", "hmmix 0.8 (Skov 2018) against 529 African 1000 Genomes founders")):
    set_meta(con, k, str(v))
set_meta(con, "archaic_ref", "Genome-wide studies put Neanderthal ancestry at about 2 percent in people of European and Asian ancestry, a little more in East "
         "Asians, and close to none in Africans. Stretch-finding methods like this one report a little less, because short or very old stretches are "
         "hard to detect. It cannot tell Neanderthal from Denisovan stretches; in Europeans nearly all are Neanderthal.")
con.commit()
log("archaic stretches %d covering %.1f Mb = %.2f%% of the called genome" % (len(confident), total / 1e6, 100 * frac))
