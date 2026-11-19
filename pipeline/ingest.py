#!/usr/bin/env python3
"""Ingest: turn a person's raw DNA file into the canonical gVCF every
later step reads (GRCh38, chromosomes named 1..22, X, Y, MT).

Accepted inputs
  chip raw data   23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, LivingDNA,
                  tellmeGen and similar exports (tab or comma separated, plain
                  or zipped), GRCh37 or GRCh38 positions
  VCF / gVCF      from any sequencing provider, GRCh37 or GRCh38, "chr1" or
                  "1" names; one sample (the first, or --sample NAME)

What comes out (in the profile's canon/ folder)
  calls.g.vcf.gz(.tbi)  variants as normal records; for a chip, every tested
                  position that matches the reference becomes a one-base
                  reference block (ALT ".", GT 0/0, END=POS), so "tested and
                  normal" stays distinct from "never tested"
  ingest.kv       what was read: format, vendor, build and how it was decided,
                  markers, no-calls, lifted, dropped and why, strand flips,
                  inferred sex and the numbers behind it

usage: ingest.py --profile DIR [--sample NAME] FILE
"""
import argparse
import bisect
import gzip
import io
import os
import re
import shutil
import subprocess
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from strandlib import FA, RAW, CHROMS, log  # noqa: E402

COMP = str.maketrans("ACGT", "TGCA")
PRIMARY = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
# X pseudoautosomal regions on GRCh38 (diploid in men)
PAR38 = [(10001, 2781479), (155701383, 156030895)]


# ---------------------------------------------------------------- reference
class Fasta:
    """single bases from the indexed GRCh38 FASTA (Ensembl names)"""
    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        for line in open(path + ".fai"):
            p = line.split("\t")
            self.idx[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))

    def base(self, c, pos):
        if c not in self.idx:
            return "N"
        ln, off, lb, lby = self.idx[c]
        if pos < 1 or pos > ln:
            return "N"
        i = pos - 1
        self.f.seek(off + (i // lb) * lby + i % lb)
        return self.f.read(1).decode().upper()

    def seq(self, c, pos, n):
        """n bases from pos (1-based), upper case"""
        return "".join(self.base(c, pos + i) for i in range(n))

    def chrom(self, c):
        """a whole chromosome as one upper-case string (index 0 = position 1)"""
        ln, off, lb, lby = self.idx[c]
        self.f.seek(off)
        nlines = (ln + lb - 1) // lb
        raw = self.f.read(nlines * lby)
        return raw.replace(b"\n", b"").replace(b"\r", b"")[:ln].decode().upper()


# ---------------------------------------------------------------- liftover
class Lifter:
    """UCSC chain file (GRCh37 -> GRCh38), primary chromosomes only.
    lift(chrom, pos) -> (chrom, pos, strand) or None; 1-based positions."""
    def __init__(self, path):
        blocks = {}
        with gzip.open(path, "rt") as f:
            keep, t, q, tpos, qpos, qstrand, qsize = False, "", "", 0, 0, "+", 0
            for line in f:
                p = line.split()
                if not p:
                    continue
                if p[0] == "chain":
                    tname, tpos = p[2], int(p[5])
                    qname, qsize, qstrand, qpos = p[7], int(p[8]), p[9], int(p[10])
                    t = tname[3:] if tname.startswith("chr") else tname
                    q = qname[3:] if qname.startswith("chr") else qname
                    keep = t in PRIMARY and q in PRIMARY
                    continue
                size = int(p[0])
                if keep:
                    blocks.setdefault(t, []).append((tpos, tpos + size, q, qpos, qstrand, qsize))
                tpos += size
                qpos += size
                if len(p) == 3:
                    tpos += int(p[1])
                    qpos += int(p[2])
        self.blocks = {}
        for t, bl in blocks.items():
            bl.sort()
            self.blocks[t] = (bl, [b[0] for b in bl])

    def lift(self, c, pos):
        if c not in self.blocks:
            return None
        bl, starts = self.blocks[c]
        i = bisect.bisect_right(starts, pos - 1) - 1
        # chains rarely overlap; look back a few blocks in case they do
        for j in range(i, max(-1, i - 8), -1):
            s, e, q, qs, strand, qsize = bl[j]
            if s <= pos - 1 < e:
                off = pos - 1 - s
                if strand == "+":
                    return (q, qs + off + 1, "+")
                return (q, qsize - (qs + off), "-")
        return None


# ---------------------------------------------------------------- chip files
def open_text(path):
    """a text stream for a plain, gzipped or zipped file (the biggest text member of a zip)"""
    if zipfile.is_zipfile(path):
        z = zipfile.ZipFile(path)
        members = [i for i in z.infolist() if not i.is_dir() and i.file_size > 0]
        if not members:
            raise SystemExit("ZIP_EMPTY: the zip file holds no data file")
        m = max(members, key=lambda i: i.file_size)
        if m.file_size > 2_000_000_000 or (m.compress_size and m.file_size / m.compress_size > 200):
            raise SystemExit("ZIP_SUSPICIOUS: the zip file's contents are implausibly large")
        if m.flag_bits & 0x1:
            raise SystemExit("ZIP_ENCRYPTED: the zip file is password protected; extract it first")
        return io.TextIOWrapper(z.open(m), encoding="utf-8", errors="replace"), m.filename
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace"), os.path.basename(path)
    return open(path, encoding="utf-8", errors="replace"), os.path.basename(path)


VENDORS = [("23andme", "23andMe"), ("ancestrydna", "AncestryDNA"), ("ancestry.com", "AncestryDNA"), ("myheritage", "MyHeritage"),
           ("family tree dna", "FamilyTreeDNA"), ("familytreedna", "FamilyTreeDNA"), ("ftdna", "FamilyTreeDNA"), ("living dna", "LivingDNA"),
           ("livingdna", "LivingDNA"), ("tellmegen", "tellmeGen"), ("nebula", "Nebula"), ("dante", "Dante Labs")]


def detect(path):
    """('vcf'|'chip'|'unknown', vendor, header comment text)"""
    f, name = open_text(path)
    head = []
    for i, line in enumerate(f):
        head.append(line.rstrip("\n"))
        if i >= 60:
            break
    f.close()
    text = "\n".join(head).lower()
    if head and head[0].startswith("##fileformat=VCF"):
        return "vcf", "", text
    vendor = next((v for k, v in VENDORS if k in text or k in name.lower()), "")
    first = next((l for l in head if l.strip() and not l.startswith("#")), "")
    low = first.lower().replace('"', "")
    if re.search(r"\b(rsid|snp|name)\b", low) and re.search(r"chrom", low):
        return "chip", vendor or "chip", text
    if re.match(r"^(rs|i)\d+[\t,]", low):
        return "chip", vendor or "chip", text
    return "unknown", vendor, text


def read_chip(path):
    """yield (rsid, chrom, pos, alleles) with chrom in 1..22/X/Y/MT, alleles a
    string like 'AG', 'A' (haploid) or '' for a no-call; 'I'/'D' kept as is"""
    f, _ = open_text(path)
    cols = None
    delim = "\t"
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        if cols is None:
            delim = "\t" if "\t" in line else ","
            p = [x.strip().strip('"').lower() for x in line.rstrip("\n").split(delim)]
            if re.match(r"^(rs|i)\d+$", p[0]):
                cols = {"rsid": 0, "chrom": 1, "pos": 2}
                if len(p) >= 5:
                    cols.update(a1=3, a2=4)
                else:
                    cols.update(gt=3)
            else:
                def find(*names):
                    for n in names:
                        for i, h in enumerate(p):
                            if h == n or h.startswith(n):
                                return i
                    return None
                cols = {"rsid": find("rsid", "snp", "name", "marker"), "chrom": find("chromosome", "chrom", "chr"),
                        "pos": find("position", "pos"), "gt": find("genotype", "result"), "a1": find("allele1"), "a2": find("allele2")}
                continue
        p = [x.strip().strip('"') for x in line.rstrip("\n").split(delim)]
        try:
            c = p[cols["chrom"]].upper().replace("CHR", "")
            pos = int(p[cols["pos"]])
        except (ValueError, IndexError, TypeError):
            continue
        c = {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT", "XY": "X", "0": ""}.get(c, c)
        if c not in PRIMARY:
            continue
        if cols.get("gt") is not None:
            g = p[cols["gt"]].upper() if cols["gt"] < len(p) else ""
        else:
            g = (p[cols["a1"]] + p[cols["a2"]]).upper() if cols["a2"] is not None and cols["a2"] < len(p) else ""
        g = g.replace("0", "").replace("-", "").replace(" ", "")
        yield (p[cols["rsid"]] if cols["rsid"] is not None else "."), c, pos, g
    f.close()


def in_par(c, pos):
    return c == "X" and any(a <= pos <= b for a, b in PAR38)


def kg_alts(markers):
    """{(chrom, pos): set of ALT alleles} from the 1000 Genomes site files, for strand checks"""
    from strandlib import kg_query
    out = {}
    for line in kg_query([(c, p) for c, p in markers], "%CHROM\t%POS\t%REF\t%ALT\n"):
        c, p, r, a = line.split("\t")
        out.setdefault((c, int(p)), set()).update(a.split(","))
    return out


def chip_to_gvcf(path, outdir, vendor, header):
    fa = Fasta(FA)
    rows = list(read_chip(path))
    if len(rows) < 1000:
        raise SystemExit("TOO_FEW_MARKERS: only %d genotyped markers were found in the file" % len(rows))
    st = {"format": "chip", "vendor": vendor, "markers_in_file": len(rows)}
    m = re.search(r"build\s*(\d+)|grch(\d+)|hg(\d+)", header)
    claimed = next((g for g in (m.groups() if m else ()) if g), "")
    # ---- which build: does the reference base match homozygous calls?
    lifter = Lifter(os.path.join(RAW, "chains", "hg19ToHg38.over.chain.gz"))
    # heterozygous calls: in the right build the reference base is one of the
    # two alleles (on either strand) almost always; in a wrong build about half the time
    hets = [r for r in rows if len(r[3]) == 2 and r[3][0] != r[3][1] and all(a in "ACGT" for a in r[3]) and r[1] in PRIMARY[:22]]
    sample = hets[::max(1, len(hets) // 4000)][:4000]
    def match(build):
        ok = n = 0
        for _, c, pos, g in sample:
            if build == "38":
                c2, p2, strand = c, pos, "+"
            else:
                L = lifter.lift(c, pos)
                if not L:
                    continue
                c2, p2, strand = L
            b = fa.base(c2, p2)
            n += 1
            # chips report the forward strand of their build; a reversed chain block flips it
            ok += b in (g if strand == "+" else g.translate(COMP))
        return ok / max(1, n)
    m38, m37 = match("38"), match("37")
    build = "37" if m37 >= m38 else "38"
    st.update(build=build, build_claimed=claimed, build_match_37="%.3f" % m37, build_match_38="%.3f" % m38)
    if max(m37, m38) < 0.90:
        raise SystemExit("UNKNOWN_BUILD: the positions match neither GRCh37 nor GRCh38 (%.0f%% / %.0f%%); an NCBI36 file?" % (100 * m37, 100 * m38))
    log("chip:", vendor, "markers", len(rows), "build", build, "(GRCh37 match %.3f, GRCh38 match %.3f)" % (m37, m38))
    # ---- sex from the data: Y calls and X heterozygosity
    ycalled = sum(1 for r in rows if r[1] == "Y" and r[3])
    ytotal = sum(1 for r in rows if r[1] == "Y")
    xhet = sum(1 for r in rows if r[1] == "X" and len(r[3]) == 2 and r[3][0] != r[3][1])
    xcalled = sum(1 for r in rows if r[1] == "X" and r[3])
    yfrac = ycalled / max(1, ytotal)
    xh = xhet / max(1, xcalled)
    sex = "male" if (ytotal == 0 or yfrac > 0.5) and xh < 0.03 else "female" if xh > 0.05 or (ytotal and yfrac < 0.2) else "unknown"
    if ytotal == 0 and xh < 0.03:
        sex = "male"
    st.update(sex=sex, y_called="%d/%d" % (ycalled, ytotal), x_het="%.4f" % xh)
    # ---- lift, place, orient
    recs = {}
    drop = {"nocall": 0, "indel": 0, "unlifted": 0, "ref_n": 0, "conflict": 0, "odd_alleles": 0}
    for rsid, c, pos, g in rows:
        if not g:
            drop["nocall"] += 1
            continue
        if any(ch not in "ACGT" for ch in g):
            drop["indel"] += 1
            continue
        if build == "37" and c != "MT":            # chip MT positions are already on rCRS, as GRCh38's MT
            L = lifter.lift(c, pos)
            if not L or L[0] != c:
                drop["unlifted"] += 1
                continue
            c2, p2, s = L
            if s == "-":
                g = g.translate(COMP)
        else:
            c2, p2 = c, pos
        ref = fa.base(c2, p2)
        if ref not in "ACGT":
            drop["ref_n"] += 1
            continue
        key = (c2, p2)
        if key in recs and recs[key][1] != g:
            recs[key] = None
            drop["conflict"] += 1
            continue
        if key not in recs:
            recs[key] = (rsid, g, ref)
    recs = {k: v for k, v in recs.items() if v}
    # strand: a call whose alleles are neither the reference nor a known
    # alternative allele, but whose complement is, was reported on the other strand
    odd = [k for k, (rs, g, ref) in recs.items() if any(a != ref for a in g)]
    alts = kg_alts(odd) if odd else {}
    flips = 0
    for k in odd:
        rs, g, ref = recs[k]
        known = {ref} | alts.get(k, set())
        if all(a in known for a in g):
            continue
        cg = g.translate(COMP)
        if alts.get(k) and all(a in known for a in cg):
            recs[k] = (rs, cg, ref)
            flips += 1
    st["strand_flips"] = flips
    # ---- write the canonical gVCF
    male = sex == "male"
    order = {c: i for i, c in enumerate(PRIMARY)}
    out = os.path.join(outdir, "calls.g.vcf")
    nref = nvar = 0
    with open(out, "w") as f:
        f.write("##fileformat=VCFv4.2\n##source=Strand Open ingest (chip)\n##reference=GRCh38\n")
        for c in PRIMARY:
            if c in fa.idx:
                f.write("##contig=<ID=%s,length=%d>\n" % (c, fa.idx[c][0]))
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End of a reference block">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FILTER=<ID=PASS,Description="Called by the chip">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for (c, p), (rs, g, ref) in sorted(recs.items(), key=lambda t: (order[t[0][0]], t[0][1])):
            haploid = c in ("Y", "MT") or (c == "X" and male and not in_par(c, p))
            if haploid and len(g) == 2:
                if g[0] != g[1]:
                    drop["odd_alleles"] += 1           # a heterozygous call where a man has one copy
                    continue
                g = g[0]
            if all(a == ref for a in g):
                f.write("%s\t%d\t%s\t%s\t.\t.\tPASS\tEND=%d\tGT\t%s\n" % (c, p, rs, ref, p, "0" if haploid else "0/0"))
                nref += 1
                continue
            alts_here = []
            for a in g:
                if a != ref and a not in alts_here:
                    alts_here.append(a)
            codes = ["0" if a == ref else str(alts_here.index(a) + 1) for a in g]
            gt = codes[0] if haploid else "/".join(sorted(codes, key=int))
            f.write("%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\tGT\t%s\n" % (c, p, rs, ref, ",".join(alts_here), gt))
            nvar += 1
    subprocess.run(["bgzip", "-f", out], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", out + ".gz"], check=True)
    st.update(tested_ref=nref, tested_variant=nvar, **{"dropped_" + k: v for k, v in drop.items()})
    log("chip: %d reference sites, %d variant sites; dropped %s; %d strand flips; sex %s (Y %s, X het %s)" % (
        nref, nvar, drop, flips, sex, st["y_called"], st["x_het"]))
    return st


# ---------------------------------------------------------------- VCF / gVCF
def symbolic(a):
    """an allele that is not a DNA sequence: <NON_REF>, <*>, <DEL>, breakends, '*' and '.'"""
    return a in (".", "*") or a.startswith("<") or "[" in a or "]" in a
END_RE = re.compile(r"(?:^|;)END=(\d+)")
MINDP_RE = re.compile(r"(?:^|;)MinDP=(\d+)")
KEEP_INFO = ("DP", "MQ", "MQ0F", "MQ0")
RENAME = {}
for _c in PRIMARY:
    RENAME[_c] = RENAME["chr" + _c] = _c
RENAME.update({"chrM": "MT", "M": "MT", "chrMT": "MT"})


def unzip_vcf(path, outdir):
    """a VCF inside a zip is extracted next to the profile's other files"""
    if not zipfile.is_zipfile(path):
        return path
    z = zipfile.ZipFile(path)
    members = [i for i in z.infolist() if not i.is_dir() and re.search(r"\.(vcf|vcf\.gz|bcf)$", i.filename, re.I)]
    if not members:
        raise SystemExit("ZIP_NO_VCF: the zip file holds no VCF")
    m = max(members, key=lambda i: i.file_size)
    if m.flag_bits & 0x1:
        raise SystemExit("ZIP_ENCRYPTED: the zip file is password protected; extract it first")
    out = os.path.join(outdir, "input" + (".vcf.gz" if m.filename.lower().endswith(".gz") else ".bcf" if m.filename.lower().endswith(".bcf") else ".vcf"))
    with z.open(m) as src, open(out, "wb") as dst:
        while True:
            buf = src.read(1 << 20)
            if not buf:
                break
            dst.write(buf)
    return out


def vcf_build(path, sample, hdr, fa, lifter):
    """'37' or '38' from REF alleles of the first SNVs (and the header's chr1 length)"""
    q = subprocess.run("bcftools view --no-version -s '%s' -v snps -H '%s' 2>/dev/null | head -5000" % (sample, path),
                       shell=True, capture_output=True, text=True).stdout.splitlines()
    m38 = m37 = n = 0
    for line in q:
        p = line.split("\t")
        c = RENAME.get(p[0])
        if c not in PRIMARY[:22] or len(p[3]) != 1:
            continue
        n += 1
        m38 += fa.base(c, int(p[1])) == p[3].upper()
        L = lifter.lift(c, int(p[1]))
        if L:
            m37 += fa.base(L[0], L[1]) == (p[3].upper() if L[2] == "+" else p[3].upper().translate(COMP))
    m38, m37 = m38 / max(1, n), m37 / max(1, n)
    m = re.search(r"##contig=<ID=(?:chr)?1,length=(\d+)", hdr)
    said = {"249250621": "37", "248956422": "38"}.get(m.group(1), "") if m else ""
    if n < 50:
        build = said or "38"
    else:
        build = "38" if m38 >= m37 else "37"
        if max(m37, m38) < 0.9:
            raise SystemExit("UNKNOWN_BUILD: the REF alleles match neither GRCh37 nor GRCh38 (%.0f%% / %.0f%%)" % (100 * m37, 100 * m38))
    return build, m37, m38, said


def vcf_to_gvcf(path, outdir, sample):
    """one sample of a VCF or gVCF as the canonical gVCF.

    gVCF reference blocks come in several dialects (Illumina and bcftools:
    ALT '.', GATK and DRAGEN: '<NON_REF>', DeepVariant: '<*>'; depth in INFO
    MinDP or FORMAT MIN_DP or DP); all become ALT '.', END and MinDP.
    Symbolic alleles leave variant records (with their AD entries), and only
    a few INFO fields stay, so no per-allele field can disagree with ALT.

    A VCF without blocks lists differences only. When it covers a whole
    genome (a million variants or more), the positions it does not list are
    taken as matching the reference inside the 1000 Genomes strict
    accessibility mask, on each chromosome it has calls for; an exome or
    panel VCF gets no such assumption. GRCh37 files are lifted to GRCh38;
    their blocks cannot be lifted reliably, so they are read as a VCF.
    The sample column is renamed SAMPLE: the file's own name may be a
    person's name or kit number."""
    fa = Fasta(FA)
    lifter = Lifter(os.path.join(RAW, "chains", "hg19ToHg38.over.chain.gz"))
    orig = path
    path = unzip_vcf(path, outdir)
    extracted = path != orig
    hdr = subprocess.run(["bcftools", "view", "-h", path], capture_output=True, text=True)
    if hdr.returncode != 0:
        raise SystemExit("BAD_VCF: bcftools cannot read this file: " + hdr.stderr.strip()[-200:])
    hdr = hdr.stdout
    samples = hdr.strip().split("\n")[-1].split("\t")[9:]
    if not samples:
        raise SystemExit("NO_SAMPLE: the VCF has no genotype column")
    if sample and sample not in samples:
        raise SystemExit("SAMPLE_NOT_FOUND: that sample is not in the VCF (it has %d samples)" % len(samples))
    sample = sample or samples[0]
    build, m37, m38, said = vcf_build(path, sample, hdr, fa, lifter)
    # gVCF dialects announce themselves in the header
    gvcf_hdr = any(k in hdr for k in ("ID=NON_REF", "<NON_REF>", "GVCFBlock", "ID=MinDP", "BLOCKAVG", "##ALT=<ID=*,"))
    lift = build == "37"
    st = {"vendor": "", "samples_in_file": len(samples), "build": build, "build_header": said,
          "build_match_37": "%.3f" % m37, "build_match_38": "%.3f" % m38}
    log("vcf: build", build, "(GRCh37 match %.3f, GRCh38 match %.3f, header says %s)" % (m37, m38, said or "nothing"))
    tmp = os.path.join(outdir, "tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    outs, last, unsorted = {}, {}, set()
    cnt = {"blocks": 0, "variants": 0, "nocall": 0, "other_contig": 0, "unlifted": 0, "ref_mismatch": 0, "blocks_not_lifted": 0}
    xhet = xhom = yalt = 0
    nvar_chrom = {}

    def put(c, pos, text):
        f = outs.get(c)
        if f is None:
            f = outs[c] = open(os.path.join(tmp, c), "w")
        if pos < last.get(c, 0):
            unsorted.add(c)
        last[c] = pos
        f.write(text)

    # outside a gVCF only the sample's own variant sites matter (-c1 after
    # subsetting: a many-sample file lists everyone's)
    cmd = ["bcftools", "view", "--no-version", "-s", sample] + ([] if gvcf_hdr else ["-c1"]) + ["-Ov", path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    seqc = [None, ""]

    def refseq(c, pos, n):
        if seqc[0] != c:
            seqc[0], seqc[1] = c, fa.chrom(c)
        return seqc[1][pos - 1:pos - 1 + n]
    for line in proc.stdout:
        if line[0] == "#":
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 10:
            continue
        c = RENAME.get(p[0])
        if c is None:
            cnt["other_contig"] += 1
            continue
        pos = int(p[1])
        fmt = p[8].split(":")
        vals = p[9].split(":")
        fv = dict(zip(fmt, vals))
        gt = fv.get("GT", "")
        alleles = re.split(r"[/|]", gt) if gt else []
        alts = p[4].split(",")
        flt = "PASS" if p[6] in (".", "PASS") else p[6]
        idx, k = {0: 0}, 0
        for i, a in enumerate(alts, 1):
            if not symbolic(a):
                k += 1
                idx[i] = k
        if k == 0:
            # a reference block, or one reference site
            if not alleles or any(a != "0" for a in alleles):
                cnt["nocall"] += 1
                continue
            if lift:
                cnt["blocks_not_lifted"] += 1
                continue
            m = END_RE.search(p[7])
            end = int(m.group(1)) if m else pos + len(p[3]) - 1
            m = MINDP_RE.search(p[7])
            dp = m.group(1) if m else fv.get("MIN_DP", fv.get("DP", "0"))
            dp = dp if dp.isdigit() else "0"
            put(c, pos, "%d\t.\t%s\t.\t.\t%s\tEND=%d;MinDP=%s\tGT:DP\t%s:%s\n" % (pos, p[3][0].upper(), flt, end, dp, gt, dp))
            cnt["blocks"] += 1
            continue
        new = []
        for a in alleles:
            j = idx.get(int(a)) if a.isdigit() else None
            new.append("." if j is None else str(j))
        if not new or all(a == "." for a in new):
            cnt["nocall"] += 1
            continue
        sep = "|" if "|" in gt else "/"
        if all(a in ("0", ".") for a in new):
            if all(a == "0" for a in new) and not lift and gvcf_hdr:
                dp = fv.get("DP", "0")
                dp = dp if dp.isdigit() else "0"
                put(c, pos, "%d\t.\t%s\t.\t.\t%s\tEND=%d;MinDP=%s\tGT:DP\t%s:%s\n" % (pos, p[3][0].upper(), flt, pos, dp, sep.join(new), dp))
                cnt["blocks"] += 1
            continue
        ref = p[3].upper()
        real = [a.upper() for a in alts if not symbolic(a)]
        if lift:
            L = lifter.lift(c, pos) if c != "MT" else ("MT", pos, "+")
            if not L or L[0] != c:
                cnt["unlifted"] += 1
                continue
            if L[2] == "-":
                if len(ref) != 1 or any(len(a) != 1 for a in real):
                    cnt["unlifted"] += 1          # indels on reversed blocks would need re-anchoring
                    continue
                ref, real = ref.translate(COMP), [a.translate(COMP) for a in real]
            pos = L[1]
        if refseq(c, pos, len(ref)) != ref:
            cnt["ref_mismatch"] += 1
            continue
        out_fmt, out_val = ["GT"], [sep.join(new)]
        if fv.get("DP", ".").isdigit():
            out_fmt.append("DP")
            out_val.append(fv["DP"])
        ad = fv.get("AD", "")
        if ad and ad != ".":
            parts = ad.split(",")
            if len(parts) == len(alts) + 1:
                keep = [parts[0]] + [parts[i] for i in range(1, len(alts) + 1) if i in idx]
                out_fmt.append("AD")
                out_val.append(",".join(keep))
        if fv.get("GQ", ".").replace(".", "", 1).isdigit():
            out_fmt.append("GQ")
            out_val.append(str(int(float(fv["GQ"]))))
        info = ";".join(kv for kv in p[7].split(";") if kv.split("=", 1)[0] in KEEP_INFO) or "."
        put(c, pos, "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pos, p[2], ref, ",".join(real), p[5], flt, info, ":".join(out_fmt), ":".join(out_val)))
        cnt["variants"] += 1
        nvar_chrom[c] = nvar_chrom.get(c, 0) + 1
        if c == "X" and 2781480 <= pos <= 155701382 and len(new) == 2 and "." not in new:
            if new[0] != new[1]:
                xhet += 1
            else:
                xhom += 1
        elif c == "Y":
            yalt += 1
    proc.wait()
    seqc[1] = ""
    for f in outs.values():
        f.close()
    if proc.returncode != 0:
        raise SystemExit("BAD_VCF: bcftools stopped reading the file (code %d)" % proc.returncode)
    if cnt["variants"] < 100:
        raise SystemExit("TOO_FEW_VARIANTS: only %d usable variant calls were found" % cnt["variants"])
    xh = xhet / max(1, xhet + xhom)
    if xhet + xhom < 50:
        sex = "unknown"
    elif xh < 0.15 and yalt >= 20:
        sex = "male"
    elif xh > 0.25 and yalt < 20:
        sex = "female"
    else:
        sex = "male" if xh < 0.15 else "female" if xh > 0.25 else "unknown"
    gvcf = gvcf_hdr and cnt["blocks"] > 1000
    wgs = cnt["variants"] >= 1_000_000
    assume = (not gvcf) and wgs
    st.update(format="gvcf" if gvcf else "vcf", sex=sex, x_het="%.4f" % xh, y_variants=yalt, ref_assumed=int(assume),
              coverage="genome" if (gvcf or wgs) else "partial", **{k: v for k, v in cnt.items()})
    log("vcf: %d variants, %d reference blocks; sex %s (X het %.3f, Y variants %d); %s" % (
        cnt["variants"], cnt["blocks"], sex, xh, yalt, "reference assumed inside the accessibility mask" if assume else "no reference assumed"))
    # ---- pass 2: chromosomes in order, sorted, with assumed reference where it applies
    mask = {}
    if assume:
        mp = os.path.join(RAW, "1kg_mask", "strict.bed")
        if not os.path.exists(mp):
            raise SystemExit("NO_MASK: the 1000 Genomes accessibility mask is missing (Setup downloads it)")
        with open(mp) as f:
            for line in f:
                q = line.split("\t", 3)
                c = RENAME.get(q[0])
                if c:
                    mask.setdefault(c, []).append((int(q[1]) + 1, int(q[2])))
    male = sex == "male"
    body = os.path.join(outdir, "calls.g.vcf.gz")
    bg = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=open(body, "wb"), text=True)
    w = bg.stdin
    w.write("##fileformat=VCFv4.2\n##source=Strand Open ingest (%s, GRCh%s%s)\n##reference=GRCh38\n" % (
        st["format"], build, ", lifted" if lift else ""))
    for c in PRIMARY:
        if c in fa.idx:
            w.write("##contig=<ID=%s,length=%d>\n" % (c, fa.idx[c][0]))
    w.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End of a reference block">\n'
            '##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Lowest depth in the block (0: not measured)">\n'
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
            '##INFO=<ID=MQ,Number=1,Type=Float,Description="Root-mean-square mapping quality">\n'
            '##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Share of reads with mapping quality 0">\n'
            '##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Reads with mapping quality 0">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Reads per allele">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
    filters = set()
    for c in outs:
        with open(os.path.join(tmp, c)) as f:
            for line in f:
                ft = line.split("\t", 6)[5]
                if ft != "PASS":
                    filters.update(ft.split(";"))
    for ft in sorted(filters):
        w.write('##FILTER=<ID=%s,Description="From the original file">\n' % ft)
    w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    nassumed = 0
    for c in PRIMARY:
        if c not in outs:
            continue
        fp = os.path.join(tmp, c)
        if c in unsorted:
            lines = open(fp).readlines()
            lines.sort(key=lambda l: int(l.split("\t", 1)[0]))
            src = iter(lines)
        else:
            src = open(fp)
        ivs = mask.get(c, [])
        if c == "MT" and assume:
            ivs = [(1, fa.idx["MT"][0])]
        if c == "Y" and (not male or nvar_chrom.get("Y", 0) < 20):
            ivs = []
        if not ivs:
            for line in src:
                w.write(c + "\t" + line)
            continue
        seq = fa.chrom(c)
        hap = male and c in ("X", "Y", "MT")
        if c == "MT":
            hap = True

        def block(a, b):
            g = "0" if hap and not (c == "X" and (a <= 2781479 or b >= 155701383)) else "0/0"
            w.write("%s\t%d\t.\t%s\t.\t.\tPASS\tEND=%d;MinDP=0\tGT:DP\t%s:0\n" % (c, a, seq[a - 1], b, g))

        nxt = None

        def advance():
            line = next(src, None)
            if line is None:
                return None
            q = line.split("\t", 4)
            return (int(q[0]), len(q[2]), line)
        nxt = advance()
        # done: the last position any written record covers. A stretch is
        # assumed to match the reference only beyond it, so nothing is
        # assumed inside a deletion, and records stay in position order.
        done = 0
        for s0, e in ivs:
            a = max(s0, done + 1)
            while a <= e:
                while nxt and nxt[0] < a:
                    w.write(c + "\t" + nxt[2])
                    done = max(done, nxt[0] + nxt[1] - 1)
                    a = max(a, done + 1)
                    nxt = advance()
                if a > e:
                    break
                if nxt and nxt[0] <= e:
                    if nxt[0] > a:
                        block(a, nxt[0] - 1)
                        nassumed += 1
                    w.write(c + "\t" + nxt[2])
                    done = max(done, nxt[0] - 1, nxt[0] + nxt[1] - 1)
                    a = max(a, done + 1)
                    nxt = advance()
                else:
                    block(a, e)
                    nassumed += 1
                    done = e
                    a = e + 1
        while nxt:
            w.write(c + "\t" + nxt[2])
            nxt = advance()
        del seq
    w.close()
    bg.wait()
    shutil.rmtree(tmp, ignore_errors=True)
    if extracted:
        os.remove(path)              # our copy extracted from the zip, never the person's own file
    st["assumed_blocks"] = nassumed
    subprocess.run(["tabix", "-f", "-p", "vcf", body], check=True)
    return st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--sample", default="")
    ap.add_argument("file")
    a = ap.parse_args()
    outdir = os.path.join(a.profile, "canon")
    os.makedirs(outdir, exist_ok=True)
    kind, vendor, header = detect(a.file)
    if kind == "chip":
        st = chip_to_gvcf(a.file, outdir, vendor, header)
    elif kind == "vcf":
        st = vcf_to_gvcf(a.file, outdir, a.sample)
    else:
        raise SystemExit("UNKNOWN_FORMAT: this file is not a DNA raw data export or VCF that Strand Open can read")
    with open(os.path.join(outdir, "ingest.kv"), "w") as f:
        for k, v in st.items():
            f.write("%s=%s\n" % (k, v))
    log("ingest done:", st)


if __name__ == "__main__":
    main()
