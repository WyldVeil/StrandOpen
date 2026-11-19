#!/usr/bin/env python3
"""Step 5: the variant table.

Reads work/variants.annot.vcf.gz (step 1), adds 1000 Genomes frequencies
(the per-chromosome site files of step 4, when present) and AlphaMissense scores
for missense changes, and writes one row per variant allele.

Tables: variant (+ indexes), gene.nvar / gene.ncoding filled in."""
import re
import subprocess
from strandlib import *

ANNOT = os.path.join(WORK, "variants.annot.vcf.gz")
WITHAF = os.path.join(WORK, "variants.af.vcf.gz")

CSQ_RANK = {  # higher = more severe
    "stop_gained": 30, "frameshift": 29, "splice_acceptor": 28, "splice_donor": 28, "start_lost": 27,
    "stop_lost": 26, "inframe_deletion": 22, "inframe_insertion": 22, "inframe_altering": 21,
    "missense": 20, "splice_region": 14, "coding_sequence": 12, "feature_elongation": 12,
    "start_retained": 11, "stop_retained": 11, "synonymous": 10, "5_prime_utr": 7, "3_prime_utr": 6,
    "non_coding": 4, "intron": 2,
}
IMPACT = {r: (3 if r >= 26 else 2 if r >= 20 else 1 if r >= 10 else 0) for r in range(0, 40)}
AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
       "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
       "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa"}

CLIN = [  # (pattern in the first CLNSIG term, code, label)
    ("pathogenic/likely_pathogenic", 6), ("pathogenic", 6), ("likely_pathogenic", 5),
    ("conflicting", 4), ("uncertain", 3), ("likely_benign", 2), ("benign/likely_benign", 1), ("benign", 1),
    ("established_risk_allele", 7), ("likely_risk_allele", 7), ("risk_factor", 7),
    ("drug_response", 8), ("protective", 9), ("association", 10), ("affects", 10),
]
# flags for secondary terms (bit values)
CLINFLAG = {"risk_factor": 1, "drug_response": 2, "protective": 4, "association": 8, "affects": 16,
            "established_risk_allele": 1, "likely_risk_allele": 1}
STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_conflicting_classifications": 1, "criteria_provided,_conflicting_interpretations": 1,
         "criteria_provided,_single_submitter": 1}


def clin_code(sig):
    if not sig or sig == ".":
        return 0, 0
    s = sig.lower()
    first = re.split(r"[|,]", s)[0]
    code = 11
    for pat, c in CLIN:
        if first.startswith(pat):
            code = c
            break
    if first in ("not_provided", "no_classification_for_the_single_variant", "other"):
        code = 11
    flags = 0
    for term in re.split(r"[|,/]", s)[1:]:
        flags |= CLINFLAG.get(term.strip("_"), 0)
    return code, flags


PLACEHOLDER_CONDITIONS = {"not provided", "not specified", "see cases"}


def cond_names(dn):
    """ClinVar's condition names, without its placeholders ('not provided', 'see cases')"""
    if not dn or dn == ".":
        return None
    names = [x for x in dn.replace("_", " ").split("|") if x.strip() and x.strip().lower() not in PLACEHOLDER_CONDITIONS]
    return "|".join(names) or None


def stars(rev):
    if not rev or rev == ".":
        return 0
    return STARS.get(rev.lower().replace("|", ","), 0)


def a3(s):
    return "".join(AA3.get(x, x) for x in s)


def hgvs_p(aa, csq):
    """bcftools '117R>117H' -> 'p.Arg117His'. Frameshifts and indels come
    from bcftools as whole altered peptides; they are cut to the usual short
    HGVS forms (p.Ser761fs, p.Leu12_Ala14del, p.Gly5_Ser6insAla)."""
    m = re.match(r"^(\d+)([A-Z*]*)(?:>(\d+)([A-Z*]*))?$", aa)
    if not m:
        return aa[:60]
    pos, a, _pos2, b = m.groups()
    pos = int(pos)
    if not a:
        return None
    if b is None or a == b:
        return "p.%s%d=" % (a3(a[0]), pos)
    if csq == "frameshift":
        return "p.%s%dfs" % (a3(a[0]), pos)
    if len(a) == 1 and len(b) == 1:
        return "p.%s%d%s" % (a3(a), pos, a3(b))
    # trim residues shared at both ends, keeping the flanks for insertions
    k = 0
    while k < len(a) and k < len(b) and a[k] == b[k]:
        k += 1
    j = 0
    while j < len(a) - k and j < len(b) - k and a[-1 - j] == b[-1 - j]:
        j += 1
    da, db = a[k:len(a) - j], b[k:len(b) - j]
    if not da and db and 0 < k < len(a):
        ins = a3(db) if len(db) <= 6 else "%d" % len(db)
        return "p.%s%d_%s%dins%s" % (a3(a[k - 1]), pos + k - 1, a3(a[k]), pos + k, ins)
    if len(b) < len(a) and a.startswith(b):
        d = a[len(b):]
        s = pos + len(b)
        return ("p.%s%ddel" % (a3(d), s)) if len(d) == 1 else ("p.%s%d_%s%ddel" % (a3(d[0]), s, a3(d[-1]), s + len(d) - 1))
    if len(b) > len(a) and b.startswith(a):
        ins = b[len(a):]
        s = pos + len(a) - 1
        txt = a3(ins) if len(ins) <= 6 else "%d" % len(ins)
        return "p.%s%d_%dins%s" % (a3(a[-1]), s, s + 1, txt)
    if b.startswith("*") or (b and b[0] == "*"):
        return "p.%s%d*" % (a3(a[0]), pos)
    return ("p.%s%d_%ddelins%s" % (a3(a[0]), pos, pos + len(a) - 1, a3(b) if len(b) <= 6 else "%daa" % len(b)))[:60]


con = db()
cur = con.cursor()
canon = set(r[0] for r in cur.execute("SELECT tx FROM gene WHERE tx!=''"))
gid = {}
for i, sym, chrom in cur.execute("SELECT id, symbol, chrom FROM gene"):
    gid.setdefault((sym, chrom), i)
log("canonical transcripts", len(canon))

src = ANNOT
SITEF = kg_site_files()
have_af = len(SITEF) >= 22
if have_af:
    if not os.path.exists(WITHAF) or os.path.getmtime(WITHAF) < max(os.path.getmtime(f) for f in SITEF.values()):
        log("adding 1000 Genomes frequencies, one chromosome at a time")
        import concurrent.futures
        tmpd = os.path.join(WORK, "af_tmp")
        os.makedirs(tmpd, exist_ok=True)
        cols = "INFO/KGAF:=INFO/AF,INFO/KGEUR:=INFO/AF_EUR,INFO/KGAFR:=INFO/AF_AFR,INFO/KGEAS:=INFO/AF_EAS," \
               "INFO/KGSAS:=INFO/AF_SAS,INFO/KGAMR:=INFO/AF_AMR,INFO/KGAC:=INFO/AC"

        def ann(c):
            out = os.path.join(tmpd, "%s.vcf.gz" % c)
            if c in SITEF:
                cmd = ["bcftools", "annotate", "-r", c, "-a", SITEF[c], "-c", cols, "-Oz", "-o", out, ANNOT]
            else:
                cmd = ["bcftools", "view", "-r", c, "-Oz", "-o", out, ANNOT]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError("annotate %s: %s" % (c, r.stderr[-400:]))
            return out
        with concurrent.futures.ThreadPoolExecutor(8) as ex:
            parts = list(ex.map(ann, CHROMS))
        r = subprocess.run(["bcftools", "concat", "-Oz", "-o", WITHAF] + parts, capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit("concat failed: " + r.stderr[-400:])
        subprocess.run(["bcftools", "index", "-f", "-t", WITHAF], check=True)
        for f in parts:
            os.remove(f)
    src = WITHAF
afcols = "\t%INFO/KGAF\t%INFO/KGEUR\t%INFO/KGAFR\t%INFO/KGEAS\t%INFO/KGSAS\t%INFO/KGAMR" if have_af else ""
# per-sample fields that a chip or a minimal VCF does not carry are written as '.'
_hdr = subprocess.run(["bcftools", "view", "-h", src], capture_output=True, text=True).stdout
_fmt = lambda tag: "[%%%s]" % tag if ("##FORMAT=<ID=%s," % tag) in _hdr else "."
fmt = ("%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/BCSQ\t%INFO/CLNVID\t%INFO/CLNSIG\t%INFO/CLNDN\t"
       "%INFO/CLNREVSTAT\t[%GT]\t" + _fmt("DP") + "\t" + _fmt("AD") + "\t" + _fmt("GQ") + afcols + "\n")

cur.executescript("""
DROP TABLE IF EXISTS variant;
CREATE TABLE variant(id INTEGER PRIMARY KEY, chrom INTEGER, pos INTEGER, end INTEGER, ref TEXT, alt TEXT,
  rsid INTEGER, kind INTEGER, zyg INTEGER, filter TEXT, pass INTEGER, qual REAL, dp INTEGER, ad_ref INTEGER,
  ad_alt INTEGER, gq INTEGER, gene TEXT, gene_id INTEGER, csq TEXT, csq_rank INTEGER, impact INTEGER,
  tx TEXT, aa TEXT, hgvs_p TEXT, dna TEXT, clin INTEGER, clin_sig TEXT, clin_dn TEXT, clin_stars INTEGER,
  clin_id INTEGER, clin_flags INTEGER, af REAL, af_eur REAL, af_afr REAL, af_eas REAL, af_sas REAL,
  af_amr REAL, am REAL, am_class TEXT);
""")


def num(s, f=float):
    if s in (".", "", None):
        return None
    try:
        return f(s)
    except ValueError:
        return None


rows = []
nrow = 0
AMKEYS = {}   # (chr, pos, ref, alt) of missense SNVs -> variant id
INS = "INSERT INTO variant VALUES(" + ",".join("?" * 39) + ")"
for line in pipe_lines("bcftools query -f '%s' %s" % (fmt.replace("'", "'\\''"), src)):
    p = line.rstrip("\n").split("\t")
    chrom = p[0]
    if chrom not in CHROM_IDX:
        continue
    pos = int(p[1])
    ref, alt = p[3], p[4]
    if alt in ("*", "."):
        continue
    rsid = int(p[2][2:]) if p[2].startswith("rs") and p[2][2:].split(";")[0].isdigit() else 0
    if len(ref) == 1 and len(alt) == 1:
        kind = 0
    elif len(ref) < len(alt) and alt.startswith(ref):
        kind = 1
    elif len(ref) > len(alt) and ref.startswith(alt):
        kind = 2
    else:
        kind = 3
    gt = p[12].replace("|", "/")
    if "/" in gt:
        a = gt.split("/")
        if a[0] == a[1]:
            zyg = 2
        elif a[0] not in ("0", ".") and a[1] not in ("0", "."):
            zyg = 4
        else:
            zyg = 1
    else:
        zyg = 3
    ad = p[14].split(",")
    adr = num(ad[0], int) if len(ad) > 0 else None
    ada = num(ad[1], int) if len(ad) > 1 else None
    # consequence: the most severe on a canonical transcript, else the most severe anywhere
    best = None
    bestc = None
    if p[7] not in (".", ""):
        for e in p[7].split(","):
            if e.startswith("@"):
                continue
            f = e.split("|")
            terms = f[0].lstrip("*").split("&")
            r = max(CSQ_RANK.get(t, 1) for t in terms)
            t = max(terms, key=lambda t: CSQ_RANK.get(t, 1))
            cand = (r, t, f)
            if best is None or r > best[0]:
                best = cand
            if len(f) > 2 and f[2] in canon and (bestc is None or r > bestc[0]):
                bestc = cand
    ch = bestc or best
    gene = csq = tx = aa = hp = dna = None
    rank = 0
    if ch:
        rank, csq, f = ch
        gene = f[1] if len(f) > 1 else None
        tx = f[2] if len(f) > 2 else None
        aa = f[5] if len(f) > 5 and f[5] else None
        dna = f[6] if len(f) > 6 and f[6] else None
        hp = hgvs_p(aa, csq) if aa else None
        if aa and len(aa) > 60:
            aa = aa[:57] + "..."
        if dna and len(dna) > 60:
            dna = dna[:57] + "..." 
    code, flags = clin_code(p[9])
    clin_id = num(p[8], int)
    nrow += 1
    row = [nrow, CHROM_IDX[chrom], pos, pos + len(ref) - 1, ref, alt, rsid, kind, zyg, p[6], 1 if p[6] == "PASS" else 0,
           num(p[5]), num(p[13], int), adr, ada, num(p[15], int), gene, gid.get((gene, CHROM_IDX[chrom])) if gene else None,
           csq, rank, IMPACT[rank], tx, aa, hp, dna, code, None if p[9] == "." else p[9],
           cond_names(p[10]), stars(p[11]) if code else None, clin_id, flags]
    if have_af:
        row += [num(x) for x in p[16:22]]
    else:
        row += [None] * 6
    row += [None, None]
    if csq == "missense" and kind == 0:
        AMKEYS[("chr" + chrom, pos, ref, alt)] = nrow
    rows.append(row)
    if len(rows) >= 100000:
        cur.executemany(INS, rows)
        rows.clear()
if rows:
    cur.executemany(INS, rows)
    rows.clear()
log("variants inserted", nrow, "missense SNVs", len(AMKEYS))

# AlphaMissense: one pass over the 71M precomputed scores, keeping the sample's
am = os.path.join(RAW, "AlphaMissense_hg38.tsv.gz")
if os.path.exists(am):
    for line in pipe_lines("zcat " + am):
        if line[0] == "#":
            continue
        p = line.split("\t", 10)
        k = (p[0], int(p[1]), p[2], p[3])
        vid = AMKEYS.get(k)
        if vid:
            rows.append((float(p[8]), p[9].strip(), vid))
    cur.executemany("UPDATE variant SET am=?, am_class=? WHERE id=?", rows)
    log("AlphaMissense scores", len(rows))
    rows.clear()

log("indexing")
cur.executescript("""
CREATE INDEX var_pos ON variant(chrom, pos);
CREATE INDEX var_rsid ON variant(rsid) WHERE rsid > 0;
CREATE INDEX var_gene ON variant(gene_id);
CREATE INDEX var_clin ON variant(clin) WHERE clin > 0;
CREATE INDEX var_impact ON variant(impact) WHERE impact >= 2;
ALTER TABLE gene ADD COLUMN nvar INTEGER;
ALTER TABLE gene ADD COLUMN ncoding INTEGER;
""") if "nvar" not in [r[1] for r in cur.execute("PRAGMA table_info(gene)")] else cur.executescript("""
CREATE INDEX var_pos ON variant(chrom, pos);
CREATE INDEX var_rsid ON variant(rsid) WHERE rsid > 0;
CREATE INDEX var_gene ON variant(gene_id);
CREATE INDEX var_clin ON variant(clin) WHERE clin > 0;
CREATE INDEX var_impact ON variant(impact) WHERE impact >= 2;
""")
cur.execute("UPDATE gene SET nvar = 0, ncoding = 0")
cur.execute("""UPDATE gene SET nvar = c.n, ncoding = c.k FROM
  (SELECT gene_id, count(*) n, sum(impact >= 1 AND csq != 'splice_region') k FROM variant WHERE gene_id IS NOT NULL GROUP BY gene_id) c
  WHERE gene.id = c.gene_id""")
for k, q in [("n_variants", "SELECT count(*) FROM variant"),
             ("n_snv", "SELECT count(*) FROM variant WHERE kind=0"),
             ("n_indel", "SELECT count(*) FROM variant WHERE kind IN (1,2)"),
             ("n_het", "SELECT count(*) FROM variant WHERE zyg IN (1,4)"),
             ("n_hom", "SELECT count(*) FROM variant WHERE zyg=2"),
             ("n_hemi", "SELECT count(*) FROM variant WHERE zyg=3"),
             ("n_pass", "SELECT count(*) FROM variant WHERE pass=1"),
             ("n_rsid", "SELECT count(*) FROM variant WHERE rsid>0"),
             ("n_coding", "SELECT count(*) FROM variant WHERE impact>=1"),
             ("n_high", "SELECT count(*) FROM variant WHERE impact=3"),
             ("n_missense", "SELECT count(*) FROM variant WHERE csq='missense'"),
             ("n_clinvar", "SELECT count(*) FROM variant WHERE clin>0"),
             ("n_novel_kg", "SELECT count(*) FROM variant WHERE af IS NULL")]:
    set_meta(con, k, cur.execute(q).fetchone()[0])
# transitions / transversions over PASS SNVs
ti = cur.execute("""SELECT count(*) FROM variant WHERE kind=0 AND pass=1 AND
  ((ref='A' AND alt='G') OR (ref='G' AND alt='A') OR (ref='C' AND alt='T') OR (ref='T' AND alt='C'))""").fetchone()[0]
tv = cur.execute("SELECT count(*) FROM variant WHERE kind=0 AND pass=1").fetchone()[0] - ti
set_meta(con, "titv", "%.3f" % (ti / max(tv, 1)))
cur.execute("DROP TABLE IF EXISTS csq_count")
cur.execute("CREATE TABLE csq_count AS SELECT coalesce(csq, '') AS csq, count(*) AS n FROM variant GROUP BY csq")
set_meta(con, "frequencies_source", "1000 Genomes 30x (NYGC, 2590 founders)" if have_af else "")
con.commit()
log("done")
