#!/usr/bin/env python3
"""Step 16: rare variants, in general and the remarkable ones.

Tables:
  rare_matrix  pass-filter variants by effect group and frequency class
  rare_sfs     histogram of the 1000 Genomes frequency of your variants
  rare_var     rare protein-altering variants with a plain description of
               the change, what the gene does and how much a broken copy
               tends to matter
  rare_gene    genes with two rare protein-altering hits (or one in both
               copies, or a damaging one on the only copy of a man's X),
               with read-backed phase where the reads reach both
  remarkable   the catalog of remarkable variants (remarkable_a/b.py)
               checked against the genome

"Not in 1000 Genomes" only counts confident heterozygous calls: a
homozygous call missing from the panel is nearly always a common variant
at a site the panel filtered out, not something new.

Needs variant.mq/mq0f from step 5b (run first here if missing). Without
1000 Genomes frequencies in the variant table (step 5 run before step 4
finished), meta rare_has_af is '0' and the frequency tables stay empty."""
import re
import subprocess
import sys
from strandlib import *
from genotype import Genotyper

HERE = os.path.dirname(os.path.abspath(__file__))
con = db()
cur = con.cursor()
# chip data: rare calls from genotyping arrays are mostly wrong (Weedon et al.
# 2021, BMJ 372:n214: under 20% of very rare chip calls confirmed by
# sequencing), so only the remarkable-variant catalogue runs, on the chip's
# own tested positions
_cm = cur.execute("SELECT value FROM meta WHERE key='call_model'").fetchone()
CHIP = bool(_cm) and _cm[0] == "chip"
# a lone re-run of step 5 rebuilds the variant table without step 5b's columns
if not CHIP and ("mq" not in [r[1] for r in cur.execute("PRAGMA table_info(variant)")] or
                 not cur.execute("SELECT 1 FROM variant WHERE mq IS NOT NULL LIMIT 1").fetchone()):
    log("no mapping quality in the variant table: running step 5b first")
    con.commit()
    subprocess.run([sys.executable, os.path.join(HERE, "05b_mapq.py")], check=True)
# every frequency-based result below would call common variants unseen
has_af = cur.execute("SELECT 1 FROM variant WHERE af IS NOT NULL LIMIT 1").fetchone() is not None
rare_ok = has_af and not CHIP
set_meta(con, "rare_has_af", "1" if rare_ok else "0")
if not has_af:
    log("no 1000 Genomes frequencies in the variant table (run steps 4 and 5): rare-variant tables left empty")

LOF = ("stop_gained", "frameshift", "splice_donor", "splice_acceptor", "start_lost")
GROUPS = [("lof", "Protein-truncating", LOF + ("stop_lost",)), ("missense", "Missense", ("missense",)),
          ("inframe", "In-frame insertion or deletion", ("inframe_deletion", "inframe_insertion", "inframe_altering")),
          ("splice", "Near a splice site", ("splice_region",)), ("syn", "Synonymous", ("synonymous", "stop_retained", "start_retained", "coding_sequence")),
          ("utr", "UTR (untranslated ends)", ("5_prime_utr", "3_prime_utr")), ("nc", "Non-coding gene", ("non_coding",)),
          ("other", "Intron or between genes", ("intron", ""))]
# confident heterozygote; a file without depth or read counts (some providers' VCFs) skips those screens
CONF_HET = ("(zyg=1 AND (dp IS NULL OR dp BETWEEN 15 AND 70) AND (ad_alt IS NULL OR (ad_alt*4 >= (ad_ref+ad_alt) "
            "AND ad_alt*4 <= 3*(ad_ref+ad_alt))))")
FCLASS = ["af IS NULL AND " + CONF_HET, "af < 0.001", "af >= 0.001 AND af < 0.01", "af >= 0.01 AND af < 0.05", "af >= 0.05"]

# ---- matrix
cur.executescript("DROP TABLE IF EXISTS rare_matrix; CREATE TABLE rare_matrix(grp TEXT, label TEXT, ord INTEGER, fclass INTEGER, n INTEGER, nhom INTEGER);")
csq_grp = {}
for gi, (g, label, cs) in enumerate(GROUPS):
    for c in cs:
        csq_grp[c] = gi
counts = {}
for csq, fc, n, nh in cur.execute(
        "SELECT coalesce(csq,''), CASE WHEN " + " WHEN ".join("%s THEN %d" % (f, i) for i, f in enumerate(FCLASS)) +
        " ELSE -1 END AS fc, count(*), sum(zyg=2) FROM variant WHERE pass=1 GROUP BY 1, 2").fetchall() if rare_ok else []:
    if fc < 0:
        continue
    gi = csq_grp.get(csq, len(GROUPS) - 1)
    k = (gi, fc)
    a = counts.get(k, (0, 0))
    counts[k] = (a[0] + n, a[1] + (nh or 0))
for (gi, fc), (n, nh) in counts.items():
    cur.execute("INSERT INTO rare_matrix VALUES(?,?,?,?,?,?)", (GROUPS[gi][0], GROUPS[gi][1], gi, fc, n, nh))
# homozygous only: X and Y in a man are single calls (zyg 3), and the panel
# has no Y or mitochondrial sites at all
unmatched_hom = cur.execute("SELECT count(*) FROM variant WHERE pass=1 AND af IS NULL AND zyg=2 AND chrom<=23").fetchone()[0] if rare_ok else 0
unmatched_low = cur.execute("SELECT count(*) FROM variant WHERE pass=1 AND af IS NULL AND zyg=1 AND NOT " + CONF_HET).fetchone()[0] if rare_ok else 0
set_meta(con, "rare_unmatched_hom", str(unmatched_hom))
set_meta(con, "rare_unmatched_lowconf", str(unmatched_low))
log("matrix done")

# ---- frequency spectrum: log10 bins from 1/5180 to 1
cur.executescript("DROP TABLE IF EXISTS rare_sfs; CREATE TABLE rare_sfs(lo REAL, hi REAL, n INTEGER, nhet INTEGER);")
# the panel leaves out singletons, so nothing sits below 2 in 5,180
edges = [10 ** (-3.5 + 0.25 * i) for i in range(15)]
edges[0] = 0.0
edges[-1] = 1.0001
for i in range(14 if rare_ok else 0):
    n, nh = cur.execute("SELECT count(*), sum(zyg=1) FROM variant WHERE pass=1 AND af >= ? AND af < ?", (edges[i], edges[i + 1])).fetchone()
    cur.execute("INSERT INTO rare_sfs VALUES(?,?,?,?)", (edges[i], min(edges[i + 1], 1.0), n, nh or 0))
log("spectrum done")

# ---- plain gene descriptions (data/gene_plain*.tsv, written by hand from UniProt)
plain = {}
for fn in sorted(os.listdir(os.path.join(HERE, "..", "data"))):
    if fn.startswith("gene_plain") and fn.endswith(".tsv"):
        for line in open(os.path.join(HERE, "..", "data", fn), encoding="ascii"):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1]:
                plain[p[0]] = p[1]
log("plain descriptions", len(plain))

AA = {"Ala": "alanine", "Arg": "arginine", "Asn": "asparagine", "Asp": "aspartic acid", "Cys": "cysteine", "Gln": "glutamine",
      "Glu": "glutamic acid", "Gly": "glycine", "His": "histidine", "Ile": "isoleucine", "Leu": "leucine", "Lys": "lysine",
      "Met": "methionine", "Phe": "phenylalanine", "Pro": "proline", "Ser": "serine", "Thr": "threonine", "Trp": "tryptophan",
      "Tyr": "tyrosine", "Val": "valine", "Ter": "a stop", "Sec": "selenocysteine"}
HGVS_POS = re.compile(r"p\.[A-Z][a-z]{2}(\d+)")
MIS = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def change_text(csq, hgvs, am, ref, alt):
    m = HGVS_POS.match(hgvs or "")
    at = " at position %s of the protein" % m.group(1) if m else ""
    if csq == "stop_gained":
        return "Puts a stop signal in the gene%s, so this copy's protein is cut short; usually the cell destroys such messages, leaving this copy silent." % at
    if csq == "frameshift":
        return "Adds or removes letters so the gene is read out of step from here%s; the rest of the protein is garbled and this copy usually makes none." % at
    if csq in ("splice_donor", "splice_acceptor"):
        return "Hits one of the two letters that mark where a piece of the gene is cut out of the message; the message is likely spliced wrongly and this copy often stops working."
    if csq == "start_lost":
        return "Removes the start signal, so this copy may make little or no protein (unless another start nearby takes over)."
    if csq == "stop_lost":
        return "Removes the stop signal, so the protein runs on past its normal end."
    if csq == "missense":
        mm = MIS.match(hgvs or "")
        s = ("Swaps one building block: %s becomes %s at position %s." % (AA.get(mm.group(1), mm.group(1)), AA.get(mm.group(3), mm.group(3)), mm.group(2))
             if mm else "Changes one building block of the protein.")
        if am is not None:
            s += (" AlphaMissense predicts the change likely disrupts the protein." if am > 0.564 else
                  " AlphaMissense cannot tell whether this harms the protein." if am >= 0.34 else
                  " AlphaMissense predicts the protein copes with it.")
        return s
    if csq.startswith("inframe"):
        d = abs(len(alt) - len(ref)) // 3
        return "%s %d building block%s without shifting the rest; the effect depends on how important that stretch is." % (
            "Removes" if len(ref) > len(alt) else "Adds", max(d, 1), "" if d == 1 else "s")
    return "Changes the protein."


def gene_means(loeuf, inh, lof):
    out = []
    if loeuf is not None and lof:
        if loeuf < 0.35:
            out.append("Healthy people very rarely carry a broken copy of this gene, so losing one copy is likely to matter somehow.")
        elif loeuf < 0.6:
            out.append("Broken copies of this gene are uncommon in healthy people.")
        elif loeuf > 1.0:
            out.append("Plenty of healthy people carry a broken copy of this gene, so one copy is usually harmless.")
    if inh:
        ad, ar = "AD" in inh, "AR" in inh
        if ar and not ad:
            out.append("Its known conditions are recessive: they need both copies changed, so one copy would make you a carrier.")
        elif ad and not ar:
            out.append("Some of its known conditions need only one changed copy.")
        elif ad and ar:
            out.append("It has both one-copy and two-copy conditions.")
        elif "XL" in inh:
            out.append("Its known conditions are X-linked.")
    return " ".join(out)


def tandem_repeat(seq, centers, maxk=70):
    """(period, span) of the longest tandem repeat that covers any of
    seq[centers], allowing some divergence between copies (a variant is often
    exactly where copies differ); (0, 0) if none. A call inside one is
    unreliable: reads from a repeat with a different number of copies
    misalign. Read in units of a multiple of its period, a repeat tolerates
    more mismatches and runs a few letters further, so the shortest period
    that covers nearly as much is the one given: (GCA)n has units of 3."""
    found = []
    L = len(seq)
    for center in centers:
        for k in range(1, maxk + 1):
            ident = 0.9 if k < 10 else 0.8
            start, end = max(0, center - 2 * maxk), min(L - k, center + 2 * maxk)
            m = [seq[i] == seq[i + k] for i in range(start, end)]
            c = center - start
            if c < 0 or c >= len(m):
                continue
            lo, hi = max(0, c - 10), min(len(m) - 1, c + 10)
            good = sum(m[lo:hi + 1])
            if good < ident * (hi - lo + 1):
                continue
            while True:
                grew = False
                if lo > 0 and good + m[lo - 1] >= ident * (hi - lo + 2):
                    lo -= 1
                    good += m[lo]
                    grew = True
                if hi < len(m) - 1 and good + m[hi + 1] >= ident * (hi - lo + 2):
                    hi += 1
                    good += m[hi]
                    grew = True
                if not grew:
                    break
            span = hi - lo + 1 + k
            if span >= max(2.5 * k, 24):
                found.append((k, span))
    if not found:
        return (0, 0)
    k, span = max(found, key=lambda t: t[1])
    return min(((p, s) for p, s in found if k % p == 0 and s + k >= span), key=lambda t: (t[0], -t[1]))


def exact_run(seq, i, k):
    """length of the exact run of a k-letter unit that covers seq[i]"""
    lo = hi = i
    while lo > 0 and seq[lo - 1] == seq[lo - 1 + k]:
        lo -= 1
    while hi + k < len(seq) and seq[hi + k] == seq[hi]:
        hi += 1
    return hi - lo + k


def repeats_at(sites):
    """sites: list of (chrom index, pos, ref, alt) -> list of (period, span)"""
    F = 140
    out = []
    for i in range(0, len(sites), 500):
        chunk = sites[i:i + 500]
        regs = ["%s:%d-%d" % (CHROMS[c - 1], max(1, p - F), p + F) for c, p, ref, alt in chunk]
        txt = subprocess.run(["samtools", "faidx", FA] + regs, capture_output=True, text=True).stdout
        seqs, cur_s = [], None
        for line in txt.splitlines():
            if line.startswith(">"):
                if cur_s is not None:
                    seqs.append(cur_s)
                cur_s = ""
            else:
                cur_s += line.strip().upper()
        seqs.append(cur_s or "")
        for (c, p, ref, alt), sq in zip(chunk, seqs):
            at = min(F, p - 1)
            if len(ref) == len(alt):
                out.append(tandem_repeat(sq, [at]))
                continue
            # a left-aligned insertion or deletion is written at the letter
            # before the stretch it changes, so look from there to a window
            # wholly past it; a short run of one or two letters repeated (A)n
            # or (CA)n is too short for the test above, but a common place
            # for reads to slip
            rp = tandem_repeat(sq, range(at, at + min(len(ref), len(alt)) + 11))
            if not rp[0]:
                for k in (1, 2):
                    n = max((exact_run(sq, j, k) for j in range(at, min(at + 3, len(sq) - k))), default=0)
                    if n >= 8:
                        rp = (k, n)
                        break
            out.append(rp)
    return out


# ---- rare protein-altering variants
cur.executescript("""DROP TABLE IF EXISTS rare_var;
CREATE TABLE rare_var(vid INTEGER, gene_id INTEGER, symbol TEXT, csq TEXT, hgvs TEXT, zyg INTEGER, af REAL, am REAL, loeuf REAL,
  clin INTEGER, verdict INTEGER, plain TEXT, change TEXT, means TEXT, disease TEXT, inheritance TEXT, chrom INTEGER, pos INTEGER, ad_ref INTEGER, ad_alt INTEGER, suspect INTEGER);""")
gdis = {}
for gid, dis, inh in cur.execute("SELECT gene, disease, inheritance FROM gene_disease").fetchall():
    d = gdis.setdefault(gid, [[], set()])
    if dis and dis not in d[0] and len(d[0]) < 3:
        d[0].append(dis)
    for t in (inh or "").split(","):
        if t.strip():
            d[1].add(t.strip())
rows = cur.execute("""SELECT v.id, v.gene_id, v.gene, v.csq, v.hgvs_p, v.zyg, v.af, v.am, g.loeuf, v.clin, v.ref, v.alt, v.chrom, v.pos, v.ad_ref, v.ad_alt,
    coalesce(v.mq, 60), coalesce(v.mq0f, 0)
  FROM variant v LEFT JOIN gene g ON g.id=v.gene_id
  WHERE v.pass=1 AND v.impact>=2 AND (v.ad_alt IS NULL OR v.ad_alt*5 >= (v.ad_ref+v.ad_alt)) AND v.clin NOT IN (1,2)
    AND ((v.af IS NOT NULL AND v.af < 0.01) OR v.af IS NULL)""").fetchall() if rare_ok else []
# a gene where much of what was called is missing from the panel usually
# sits in duplicated DNA where reads land in the wrong copy (PABPC1, NBPF,
# PRAMEF): its rare calls are mostly mapping artefacts
hard = {}
for gid, n, u in cur.execute("SELECT gene_id, count(*), avg(af IS NULL) FROM variant WHERE gene_id IS NOT NULL GROUP BY gene_id"):
    if n >= 3 and u >= 0.4:
        hard[gid] = u
# gene families with close copies elsewhere in the genome (segmental
# duplications), or DNA that varies unusually between people (HLA, KIR,
# LILR), where short reads can land in the wrong place whatever the panel
# shows. Not the zinc-finger (ZNF) genes: most map uniquely, and the other
# screens here catch those that do not
SEGDUP = re.compile(r"^(NBPF|PRAMEF|GOLGA[68]|POTE[A-Z]|ANKRD(20|30|36)|NPIP[AB]|TBC1D3|SPATA31|FAM90A|USP17L|RGPD|CT45|GAGE|XAGE|SPANX|"
                    r"TP53TG3|FRG2|DUX4|CDY|TSPY|SMN[12]|NAIP|SERF1|GTF2IRD2|NCF1|CCDC144|LRRC37A|ARHGAP11|CHRFAM7A|HERC2P|NOTCH2NL|SRGAP2|"
                    r"FAM72|NOMO[123]|PDPR|OR[0-9]+[A-Z]+[0-9]+|KRTAP|LILR|KIR[23]D|HLA-|MUC[0-9]|CROCC|PSG[0-9]|SIGLEC1[0-9]|"
                    r"FCGR[23]|AMY[12]|CYP2D[67]|CYP4F|UGT2B|PMS2P|SULT1A[34]|TCEB3C|CTAGE|MAGEA|SSX[0-9])")
family = {gid for gid, sym in cur.execute("SELECT id, symbol FROM gene").fetchall() if sym and SEGDUP.match(sym) and gid not in hard}
HARD_NOTE = " Much of what was called in this gene is missing from the reference panel, a sign of duplicated DNA where reads land in the wrong copy: treat this call with suspicion."
FAMILY_NOTE = (" This gene belongs to a family of closely similar genes, or sits in DNA that varies unusually between people, "
               "where short reads can land in the wrong place: treat this call with suspicion.")
REP = dict(zip([r[0] for r in rows], repeats_at([(r[12], r[13], r[10], r[11]) for r in rows])))
# a cluster of calls missing from the panel within a few dozen letters is the
# signature of reads misaligned as a block
import bisect
unseen_pos = {}
if rare_ok:
    for c, p in cur.execute("SELECT chrom, pos FROM variant WHERE af IS NULL ORDER BY chrom, pos"):
        unseen_pos.setdefault(c, []).append(p)


def unseen_near(c, p, w=60):
    a = unseen_pos.get(c, [])
    return bisect.bisect_right(a, p + w) - bisect.bisect_left(a, p - w) - 1
log("rare calls inside tandem repeats", sum(1 for r in rows if REP[r[0]][0] and r[6] is None))
nv = 0
for vid, gid, sym, csq, hgvs, zyg, af, am, loeuf, clin, ref, alt, chrom, pos, adr, ada, mq, mq0f in rows:
    rep = af is None and (REP[vid][0] > 0 or unseen_near(chrom, pos) >= 2)
    tot = (adr or 0) + (ada or 0)
    mismap = mq < 40 or mq0f > 0.1
    dp_ok = 15 <= tot <= 70
    lowconf = af is None and (zyg != 1 or not dp_ok or ada * 4 < tot or ada * 4 > 3 * tot)
    if af is None and zyg >= 2:
        # a homozygous call missing from the panel: a filtered common site. So
        # are most single-copy calls missing from it (X in a man: they go
        # missing about as often as homozygous ones), and the panel has no Y
        # or mitochondrial sites
        continue
    suspect = gid in hard or gid in family or lowconf or mismap or rep
    note = ""
    if lowconf:
        if ada * 4 < tot:
            note = " Only %d of %d reads show it" % (ada, tot)
        elif ada * 4 > 3 * tot:
            note = " %d of %d reads show it: too many for one copy, too few for two" % (ada, tot)
        else:
            note = " The read depth here is unusual (%d reads)" % tot
        note += ", and it is not in the reference panel: it may be a sequencing artefact."
    if mismap:
        note += " The reads here map poorly (average mapping quality %d of 60): they may belong to a similar stretch of DNA elsewhere, so the call is doubtful." % mq
    elif rep:
        note += (" It sits in a stretch of DNA repeated several times in a row (units of %d letter%s), where reads misalign easily; "
                 "a change found only here and in no one else is often an artefact." % (REP[vid][0], "" if REP[vid][0] == 1 else "s") if REP[vid][0] else
                 " Several other changes nobody else carries sit within a few dozen letters of it, the usual sign of reads misaligned as a block: "
                 "it is probably an artefact.")
    elif gid in hard:
        note += HARD_NOTE
    elif gid in family:
        note += FAMILY_NOTE
    lof = csq in LOF
    if lof:
        verdict = 3 if (loeuf is not None and loeuf < 0.6) else 2
    elif csq == "missense":
        verdict = 3 if (am or 0) > 0.564 else 2 if (am or 0) >= 0.34 or am is None else 1
    elif csq.startswith("inframe") or csq == "stop_lost":
        verdict = 2
    else:
        verdict = 1
    if suspect:
        verdict = min(verdict, 1)
    d = gdis.get(gid, [[], set()])
    inh = ",".join(sorted(d[1]))
    cur.execute("INSERT INTO rare_var VALUES(" + ",".join("?" * 21) + ")",
                (vid, gid, sym, csq, hgvs, zyg, af, am, loeuf, clin, verdict, plain.get(sym, ""), change_text(csq, hgvs, am, ref, alt),
                 (gene_means(loeuf, inh, lof) + note).strip(),
                 "; ".join(d[0]), inh, chrom, pos, adr, ada, int(suspect)))
    nv += 1
cur.execute("CREATE INDEX rare_var_v ON rare_var(verdict, af)")
log("rare protein-altering variants", nv)


# ---- genes with two hits, phase from the reads where possible
def read_base(seq, cigar, start, pos):
    """base of a read at a reference position, None if not covered"""
    rp, qp = start, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            if rp <= pos < rp + n:
                return seq[qp + pos - rp]
            rp += n
            qp += n
        elif op in "IS":
            qp += n
        elif op in "DN":
            if rp <= pos < rp + n:
                return "-"
            rp += n
    return None


def phase(chrom, a, b):
    """a, b: (pos, ref, alt) SNVs; returns (phase, fragments spanning both):
    'cis' or 'trans' when 4 or more span both and 80% of them agree, 'mixed'
    when they do not agree, 'few' under 4, 'far' more than 700 letters apart
    (no read spans both), 'noreads' without a readable alignment file"""
    if abs(a[0] - b[0]) > 700:
        return "far", 0
    cram = os.path.join(DATA, "align", "sample.cram")
    if not os.path.exists(cram):
        return "noreads", 0
    reg = "%s:%d-%d" % (CHROMS[chrom - 1], min(a[0], b[0]), max(a[0], b[0]))
    frag = {}
    try:
        out = subprocess.run(["samtools", "view", "-T", FA, "-F", "0xF04", "-q", "20", cram, reg], capture_output=True, text=True, timeout=60, check=True).stdout
    except Exception:
        return "noreads", 0
    for line in out.splitlines():
        q = line.split("\t")
        f = frag.setdefault(q[0], {})
        for k, (p, rf, al) in enumerate((a, b)):
            base = read_base(q[9], q[5], int(q[3]), p)
            if base == al:
                f[k] = 1
            elif base == rf:
                f[k] = 0
    cis = trans = 0
    for f in frag.values():
        if 0 in f and 1 in f:
            if f[0] == 1 and f[1] == 1:
                cis += 1
            elif f[0] != f[1]:
                trans += 1
    n = cis + trans
    if n < 4:
        return "few", n
    if cis >= 0.8 * n:
        return "cis", n
    if trans >= 0.8 * n:
        return "trans", n
    return "mixed", n


# phase: 'hom' a change on both copies, 'hemi' on the only copy (X in a man),
# 'multi' more than two changes, 'notsnv' one is not a single-letter change,
# else from the reads (see phase()); not attempted for the first four
cur.executescript("""DROP TABLE IF EXISTS rare_gene;
CREATE TABLE rare_gene(gene_id INTEGER, symbol TEXT, n INTEGER, nhom INTEGER, vids TEXT, changes TEXT, phase TEXT, phase_reads INTEGER,
  plain TEXT, disease TEXT, inheritance TEXT, loeuf REAL, worst INTEGER, suspect INTEGER, ndmg INTEGER, nhemi INTEGER);""")
bygene = {}
for r in cur.execute("SELECT vid, gene_id, symbol, csq, hgvs, zyg, verdict, chrom, pos, suspect FROM rare_var WHERE NOT (suspect AND af IS NULL AND (ad_alt*4 < ad_ref+ad_alt OR ad_ref+ad_alt < 15 OR ad_ref+ad_alt > 70)) ORDER BY pos").fetchall():
    bygene.setdefault(r[1], []).append(r)


def leaves_none(v):
    """one change that leaves no untouched copy: on both copies, or a
    damaging one on the only copy"""
    return v[5] in (2, 4) or (v[5] == 3 and v[6] >= 2)


def n_events(vs):
    """changes a few letters apart are usually one complex change written as several"""
    return sum(1 for i, v in enumerate(vs) if i == 0 or v[8] - vs[i - 1][8] > 15)


for gid, vs in bygene.items():
    nhom = sum(1 for v in vs if v[5] in (2, 4))
    nhemi = sum(1 for v in vs if v[5] == 3)
    events = n_events(vs)
    if events < 2 and not any(leaves_none(v) for v in vs):
        continue
    # the gene is doubtful when its calls that passed every screen no longer make a pair
    clean = [v for v in vs if not v[9]]
    suspect = gid in hard or events >= 5 or (n_events(clean) < 2 and not any(leaves_none(v) for v in clean))
    ph, nreads = "", 0
    if nhom:
        ph = "hom"
    elif nhemi:
        ph = "hemi"
    elif len(vs) > 2:
        ph = "multi"
    else:
        refalt = {r[0]: (r[1], r[2]) for r in cur.execute("SELECT id, ref, alt FROM variant WHERE id IN (?,?)", (vs[0][0], vs[1][0]))}
        a, b = refalt[vs[0][0]], refalt[vs[1][0]]
        if len(a[0]) == len(a[1]) == 1 and len(b[0]) == len(b[1]) == 1:
            ph, nreads = phase(vs[0][7], (vs[0][8], a[0], a[1]), (vs[1][8], b[0], b[1]))
        else:
            ph = "far" if vs[1][8] - vs[0][8] > 700 else "notsnv"
    d = gdis.get(gid, [[], set()])
    lo = cur.execute("SELECT loeuf FROM gene WHERE id=?", (gid,)).fetchone()
    word = {3: "likely disruptive", 2: "possibly disruptive", 1: "probably harmless"}
    chg = ", ".join("%s (%s%s)" % (v[4] or v[3], "doubtful call" if v[9] else word[v[6]],
                                   ", both copies" if v[5] in (2, 4) else ", your only copy" if v[5] == 3 else "")
                    for v in vs[:6]) + (", ..." if len(vs) > 6 else "")
    cur.execute("INSERT INTO rare_gene VALUES(" + ",".join("?" * 16) + ")",
                (gid, vs[0][2], len(vs), nhom, ",".join(str(v[0]) for v in vs), chg, ph, nreads,
                 plain.get(vs[0][2], ""), "; ".join(d[0]), ",".join(sorted(d[1])), lo[0] if lo else None, max(v[6] for v in vs), int(suspect),
                 sum(1 for v in vs if v[6] >= 2), nhemi))
log("genes with two or more hits", cur.execute("SELECT count(*) FROM rare_gene").fetchone()[0])

# ---- remarkable variants
ENTRIES = []
for mod in sorted(f[:-3] for f in os.listdir(HERE) if f.startswith("remarkable_") and f.endswith(".py")):
    try:
        ENTRIES += __import__(mod).ENTRIES
    except ImportError:
        log("no", mod)
G = Genotyper(con)
cur.executescript("""DROP TABLE IF EXISTS remarkable;
CREATE TABLE remarkable(key TEXT, cat TEXT, gene TEXT, title TEXT, status INTEGER, copies INTEGER, found TEXT, vid INTEGER,
  effect TEXT, absent TEXT, rarity TEXT, source TEXT, pmid TEXT, conf TEXT, need TEXT, af REAL, kind TEXT, rsids TEXT);""")
# status: 2 carries (enough copies for the effect), 1 carries one copy of a two-copy effect, 0 does not carry, -1 not called
# entries whose listed alleles are never found on the same copy, so one of
# each is two copies (APOL1 G1 and G2: Genovese 2010, and none of the 1000
# Genomes founders carrying both has them together)
IN_TRANS = {"b_apol1_g1g2"}
for e in ENTRIES:
    # copies: the most copies seen of any one listed allele. Two different
    # alleles, one copy each, may sit on the same copy (a haplotype such as
    # PER3 P415A/H417R) or on different ones: phase unknown, so they only
    # count as one copy (and say so) unless an allele is itself homozygous.
    # A call on a man's only X counts as two: there is no other copy.
    status, copies, found, vid, af, nhet = 0, 0, [], None, None, 0
    need2 = e.get("copies") == "two"
    chip_rare_hit = False
    if e["kind"] == "alleles":
        anycall = False
        for chrom, pos, ref, alt, label in e["alleles"]:
            ci = CHROM_IDX[str(chrom)]
            r = cur.execute("SELECT id, zyg, af, coalesce(filter,'') FROM variant WHERE chrom=? AND pos=? AND ref=? AND alt=?", (ci, pos, ref, alt)).fetchone()
            if r and r[3] == "ChipRare":
                chip_rare_hit = True
                r = None
            if r:
                anycall = True
                copies = max(copies, 1 if r[1] == 1 else 2)
                nhet += r[1] == 1
                found.append("%s (%s)" % (label, "two copies" if r[1] == 2 else "one copy" if r[1] == 1 else "single copy"))
                if vid is None or (af is not None and r[2] is not None and r[2] < af):
                    vid, af = r[0], r[2]
            else:
                g = G.at(ci, pos)
                if g[0] in ("ref", "var"):
                    anycall = True
        if not anycall:
            status = -1
        if chip_rare_hit and not copies:
            status = -1
            found.append("the chip reported it, but chip calls this rare are usually wrong, so it counts as untested")
    else:
        gid = cur.execute("SELECT id FROM gene WHERE symbol=? AND biotype='protein_coding'", (e["gene"],)).fetchone()
        if not rare_ok:
            status = -1   # rare_var is empty without frequencies or on a chip: not checked
        elif gid:
            # the quality-checked rare calls only (mapping, repeat and cluster screens of
            # rare_var), and rarer than 1 in 1,000: the gene studies counted ultra-rare
            # variants, and sequencing-artefact hotspots in genes such as SETD1A sit at
            # a few per thousand in 1000 Genomes
            for r in cur.execute("SELECT vid, hgvs, csq, zyg, af FROM rare_var WHERE gene_id=? AND suspect=0 AND csq IN (%s) "
                                 "AND (af IS NULL OR af < 0.001)" % ",".join("'%s'" % c for c in LOF), (gid[0],)).fetchall():
                copies = max(copies, 1 if r[3] == 1 else 2)
                nhet += r[3] == 1
                found.append("%s, %s (%s)" % (r[1] or "", r[2].replace("_", " "), "two copies" if r[3] in (2, 4) else "single copy" if r[3] == 3 else "one copy"))
                vid, af = r[0], r[4]
    if need2 and copies == 1 and nhet >= 2:
        if e.get("trans") or e["key"] in IN_TRANS:
            copies = 2
            found.append("one of each: they are never found on the same copy, so this counts as two copies")
        else:
            found.append("two different variants: they would need to be on different copies, which short reads cannot tell here")
    if copies:
        status = 2 if (copies >= 2 or not need2) else 1
        if CHIP and (af is None or af < 0.001):
            found.append("a chip call of a very rare variant: these can be wrong, so a clinical test may make sense before acting on it")
    rsids = ",".join(dict.fromkeys(re.findall(r"rs\d+", " ".join(a[4] for a in e["alleles"]))))
    cur.execute("INSERT INTO remarkable VALUES(" + ",".join("?" * 18) + ")",
                (e["key"], e["cat"], e["gene"], e["title"], status, copies, "; ".join(found), vid, e["effect"], e["absent"], e["rarity"],
                 e["source"], e.get("pmid", ""), e["conf"], e.get("copies", "one"), af, e["kind"], rsids))
log("remarkable entries", len(ENTRIES), "carried", cur.execute("SELECT count(*) FROM remarkable WHERE status>0").fetchone()[0])
# ---- questions genes cannot answer (yet): data/unanswerable.tsv
cur.executescript("DROP TABLE IF EXISTS unanswerable; CREATE TABLE unanswerable(ord INTEGER, title TEXT, answer TEXT, explanation TEXT, sources TEXT, pmids TEXT);")
ua = os.path.join(HERE, "..", "data", "unanswerable.tsv")
if os.path.exists(ua):
    for i, line in enumerate(open(ua, encoding="ascii").read().splitlines()[1:]):
        f = line.split("\t")
        if len(f) >= 5:
            cur.execute("INSERT INTO unanswerable VALUES(?,?,?,?,?,?)", (i, f[0], f[1], f[2], f[3], f[4]))
con.commit()
log("done")
