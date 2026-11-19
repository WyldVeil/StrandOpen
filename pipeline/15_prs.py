#!/usr/bin/env python3
"""Step 15: polygenic scores (PGS Catalog, harmonised to GRCh38).

Each score's weights are matched to the 1000 Genomes variants of step 4 and
turned into ALT-allele weights (an effect allele that is the reference
allele flips the sign; a constant shift changes nobody's rank). plink2
scores the 2,590 founders; the sample is scored on exactly the same
variants with its own calls, a missing call taking the founders' mean
dosage, as plink2 does. The sample's place is its percentile among the
European founders (525 people) and among everyone.

Tables: prs (one row per score), prs_dist (the European distribution as a
histogram)."""
import csv
import gzip
import re
import json
import subprocess
import numpy as np
from strandlib import *
from genotype import Genotyper

PG = os.path.join(RAW, "pgs")
P = KG
OUT = os.path.join(WORK, "prs")
os.makedirs(OUT, exist_ok=True)
# the chosen scores live with the program (data/), their weight files in the shared reference
DATAD = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
ids = [l.split("\t")[0].strip() for l in open(os.path.join(DATAD, "pgs_selected.tsv")) if l.strip()]
# the mind and behaviour scores: ID, short label, a plain note on what the score can and cannot say
MIND = {}
mf = os.path.join(DATAD, "pgs_selected_mind.tsv")
if os.path.exists(mf):
    for l in open(mf, encoding="ascii"):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 3 and f[0] not in ids:
            MIND[f[0]] = (f[1], f[2])
            ids.append(f[0])
meta = {r["Polygenic Score (PGS) ID"]: r for r in csv.DictReader(open(os.path.join(RAW, "pgs_all_metadata_scores.csv"), encoding="utf-8"))}
perf = {}
for r in csv.DictReader(open(os.path.join(RAW, "pgs_all_metadata_performance_metrics.csv"), encoding="utf-8")):
    perf.setdefault(r.get("Evaluated Score"), []).append(r)
pss_anc = {}
pssf = os.path.join(RAW, "pgs_all_metadata_evaluation_sample_sets.csv")
if os.path.exists(pssf):
    for r in csv.DictReader(open(pssf, encoding="utf-8")):
        pss_anc.setdefault(r.get("PGS Sample Set (PSS)"), set()).add(r.get("Broad Ancestry Category", ""))
# Scores whose published weights run the wrong way round: both GSCAN 2022
# scores (Saunders 2022) as harmonised in the PGS Catalog. Evidence, from the
# 1000 Genomes Europeans after removing ancestry: smoking initiation
# correlated -0.47 with the risk-taking score, +0.21 with education and +0.40
# with age at first sex (all three known genetic correlations have the
# opposite sign); drinks per week correlated -0.56 with the risk-taking score
# (which is built partly from drinking), and its KLB rs11940694 A allele
# raises the score although A lowers drinking (Schumann 2016, PNAS). After the
# reversal, an independent drinks-per-week score (PGS000203, GSCAN 2019)
# correlates +0.86 with PGS003376.
FLIP = {"PGS003360", "PGS003376"}
# which end of a score is the one to watch: +1 higher means more risk (most
# disease scores), -1 higher is favourable, 0 neither (body measures, and
# every mind and behaviour score)
PROTECTIVE = {"PGS003768", "PGS002315", "PGS000882"}   # HDL, bone mineral density, vitamin D
NEUTRAL = {"PGS002977"}                                 # height
CAT = {"Coronary": "Heart and circulation", "Atrial": "Heart and circulation", "Systolic": "Heart and circulation", "LDL": "Heart and circulation",
       "HDL": "Heart and circulation", "triglycerides": "Heart and circulation", "Venous": "Heart and circulation", "stroke": "Heart and circulation",
       "Type 2": "Metabolism", "Body mass": "Metabolism", "Gout": "Metabolism", "vitamin": "Metabolism", "Hypothyroid": "Metabolism",
       "kidney": "Metabolism", "Height": "Body", "Bone": "Body", "Prostate": "Cancer", "Colorectal": "Cancer", "Lung": "Cancer",
       "Melanoma": "Cancer", "Basal": "Cancer", "Alzheimer": "Brain", "Multiple sclerosis": "Immune", "bowel": "Immune", "Asthma": "Immune",
       "Psoriasis": "Immune", "Celiac": "Immune", "Type 1": "Immune", "rheumatoid": "Immune", "Glaucoma": "Eyes", "macular": "Eyes"}


NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE]-?\d+)?")


def tidy_metric(label, value):
    """'Other Metric(s)' 'R\u00b2 = 0.1316966' -> 'R-squared 0.132'; numbers to 3 significant digits"""
    v = value.replace("\u00b2", "2").replace("R^2", "R2").replace("R2 =", "R-squared").replace("R2=", "R-squared ").replace("R2", "R-squared")
    v = v.replace("\u2206", "change in ").replace("  ", " ")
    if label.startswith("Other"):
        label = ""
    v = NUM.sub(lambda m: "%.3g" % float(m.group()), v)
    v = "".join(ch for ch in v if ord(ch) < 128)
    return ("%s %s" % (label, v)).strip()


def same_trait(a, b):
    a, b = (a or "").lower().strip(), (b or "").lower().strip()
    return bool(a and b) and (a == b or a in b or b in a)


def category(trait):
    for k, v in CAT.items():
        if k.lower() in trait.lower():
            return v
    return "Other"


# ---- the weights of every score, split by chromosome into work files; each
# chromosome is then read alone, so memory stays near a gigabyte (all of them
# at once took over 12 GB)
K = len(ids)
nvar = [0] * K
SPLIT = os.path.join(OUT, "weights_by_chrom")
os.makedirs(SPLIT, exist_ok=True)
wout = {c: open(os.path.join(SPLIT, "chr%s.tsv" % c), "w") for c in [str(i) for i in range(1, 23)]}
for k, sid in enumerate(ids):
    path = os.path.join(PG, "%s_hmPOS_GRCh38.txt.gz" % sid)
    if not os.path.exists(path):
        log("missing", sid)
        continue
    with gzip.open(path, "rt") as f:
        hdr = None
        for line in f:
            if line.startswith("#"):
                continue
            if hdr is None:
                hdr = {h: i for i, h in enumerate(line.rstrip("\n").split("\t"))}
                continue
            p = line.rstrip("\n").split("\t")
            try:
                c = p[hdr["hm_chr"]]
                pos = int(p[hdr["hm_pos"]])
                w = float(p[hdr["effect_weight"]]) * (-1.0 if sid in FLIP else 1.0)
            except (ValueError, KeyError, IndexError):
                continue
            ea = p[hdr["effect_allele"]].upper()
            oa = ""
            if "other_allele" in hdr and p[hdr["other_allele"]]:
                oa = p[hdr["other_allele"]].upper()
            elif "hm_inferOtherAllele" in hdr:
                oa = p[hdr["hm_inferOtherAllele"]].upper()
            if c not in CHROM_IDX or c in ("X", "Y", "MT"):
                continue
            wout[c].write("%d\t%d\t%s\t%s\t%r\n" % (pos, k, ea, oa, w))
            nvar[k] += 1
    log(sid, "weights", nvar[k])
for f in wout.values():
    f.close()


def weights_of(c):
    """{position: [(score index, effect allele, other allele, weight)]} for one chromosome"""
    out = {}
    with open(os.path.join(SPLIT, "chr%s.tsv" % c)) as f:
        for line in f:
            pos, k, ea, oa, w = line.rstrip("\n").split("\t")
            out.setdefault(int(pos), []).append((int(k), ea, oa, float(w)))
    return out

G = Genotyper(con := db())
cur = con.cursor()
fam = [l.split() for l in open(os.path.join(P, "pgen", "chr1.psam")) if not l.startswith("#")]
samples = [x[0] if len(x) == 1 else x[1] if fam and len(fam[0]) > 2 else x[0] for x in fam]
psam_hdr = open(os.path.join(P, "pgen", "chr1.psam")).readline()
popmap = {}
for line in open(os.path.join(P, "founders.pop.tsv")):
    q = line.split()
    popmap[q[0]] = (q[2], q[3])

total_ref = None
user = np.zeros(K)
matched = np.zeros(K, dtype=int)
called = np.zeros(K, dtype=int)
for c in [str(i) for i in range(1, 23)]:
    need_c = weights_of(c)
    if not need_c:
        continue
    pv = os.path.join(P, "pgen", "chr%s.pvar" % c)
    rows = []          # (id, alt, weights[K])
    ci = CHROM_IDX[c]
    with open(pv) as f:
        for line in f:
            if line.startswith("#"):
                continue
            q = line.split("\t", 5)
            pos = int(q[1])
            ws = need_c.get(pos)
            if not ws:
                continue
            ref, alt = q[3], q[4].strip()
            wv = np.zeros(K)
            hit = False
            for k, ea, oa, w in ws:
                if ea == alt and (not oa or oa == ref):
                    wv[k] += w
                    hit = True
                elif ea == ref and (not oa or oa == alt):
                    wv[k] -= w
                    hit = True
            if not hit:
                continue
            # only variants this person was actually called at take part, on
            # both sides: a chip tests a fraction of a score's variants, and
            # filling the rest with population averages would pull every
            # chip result towards the middle
            r = G.at(ci, pos)
            d = None
            if r[0] == "ref":
                d = 0.0
            elif r[0] == "var":
                for vref, valt, zyg, ps, dp in r[1]:
                    if vref == ref and valt == alt:
                        d = 2.0 if zyg == 2 else 1.0
                        break
                else:
                    # called here, carrying a different alternative allele
                    if all(vref == ref for vref, *_ in r[1]):
                        d = 0.0
            nz = wv != 0
            matched[nz] += 1
            if d is None:
                continue
            called[nz] += 1
            user += wv * d
            rows.append((q[2], ref, alt, pos, wv))
    if not rows:
        continue
    wf = os.path.join(OUT, "w%s.tsv" % c)
    with open(wf, "w") as f:
        f.write("ID\tA1\t" + "\t".join("S%d" % k for k in range(K)) + "\n")
        for vid, ref, alt, pos, wv in rows:
            f.write("%s\t%s\t%s\n" % (vid, alt, "\t".join("%.6g" % x for x in wv)))
    pf = os.path.join(P, "pgen", "chr%s" % c)
    subprocess.run(f"plink2 --pfile {pf} --score {wf} 1 2 header-read cols=scoresums,denom no-mean-imputation --score-col-nums 3-{K + 2} "
                   f"--threads {NPROC} --out {OUT}/s{c} > {OUT}/s{c}.log 2>&1", shell=True, check=True)
    # reference scores
    with open(f"{OUT}/s{c}.sscore") as f:
        h = f.readline().lstrip("#").split()
        cols = [h.index("S%d_SUM" % k) for k in range(K)]
        idcol = h.index("IID")
        sc = []
        names = []
        for line in f:
            q = line.split()
            names.append(q[idcol])
            sc.append([float(q[j]) for j in cols])
    sc = np.array(sc)
    total_ref = sc if total_ref is None else total_ref + sc
    log("chr", c, "variants used", len(rows))
ref_names = names
# Ancestry adjustment. Scores drift with ancestry even inside Europe (a
# British genome scores systematically differently from a Finnish or
# Spanish one for reasons that are mostly technical), so each score is
# regressed on the European principal components (step 13: the founders and
# you projected onto the same axes) and you are placed by the residual.
pcs = {}
try:
    # references by their out-of-sample (leave-one-out) coordinates, you by
    # your projection: both projected onto a PCA that did not include them
    for sett, smp, a1, a2, a3, a4 in cur.execute("SELECT sett, sample, pc1, pc2, pc3, pc4 FROM anc_pca "
                                                 "WHERE (sett='eur_oos') OR (sett='eur' AND sample='YOU')"):
        pcs[smp] = (a1, a2, a3, a4)
except Exception:
    pcs = {}
np.save(os.path.join(OUT, "ref_scores.npy"), total_ref)
with open(os.path.join(OUT, "ref_scores.ids"), "w") as f:
    f.write("\n".join(ids) + "\n" + "#samples\n" + "\n".join(names) + "\n")
eur = np.array([popmap.get(s, ("", ""))[1] == "EUR" for s in ref_names])
# The reference group follows the person's ancestry: a score's distribution
# differs between continents (allele frequencies and linkage differ, and most
# scores come from European studies), so a Yoruba genome placed among
# Europeans reads as an extreme for no reason. One continental group making
# up 80 percent or more of the ancestry gives that group's founders;
# otherwise everyone, and the results are marked as rough.
GROUP_LABEL = {"EUR": ("Europeans", "in Europeans"), "AFR": ("people of African ancestry", "in African ancestry"),
               "EAS": ("people of East Asian ancestry", "in East Asian ancestry"), "SAS": ("people of South Asian ancestry", "in South Asian ancestry"),
               "AMR": ("people of the 1000 Genomes American populations", "in American populations"), "ALL": ("everyone", "in everyone")}
group = "EUR"
try:
    r = cur.execute("SELECT pop, frac FROM anc_admix WHERE level='super' ORDER BY frac DESC LIMIT 1").fetchone()
    if r:
        group = r[0] if r[1] >= 0.8 and r[0] in GROUP_LABEL else "ALL"
except Exception:
    pass
sel = eur if group == "EUR" else np.ones(len(ref_names), bool) if group == "ALL" else \
    np.array([popmap.get(s, ("", ""))[1] == group for s in ref_names])
log("reference group", group, int(sel.sum()), "founders")
adjust = group == "EUR" and "YOU" in pcs and all(n in pcs for n, e in zip(ref_names, eur) if e)
if adjust:
    Xe = np.array([[1.0] + list(pcs[n]) for n, e in zip(ref_names, eur) if e])
    xu = np.array([1.0] + list(pcs["YOU"]))
log("ancestry adjustment", "on" if adjust else "off (no European PCA from step 13)")
cur.executescript("""DROP TABLE IF EXISTS prs; DROP TABLE IF EXISTS prs_dist;
CREATE TABLE prs(id TEXT, trait TEXT, category TEXT, name TEXT, n_weights INTEGER, n_used INTEGER, n_called INTEGER, score REAL,
  pct_eur REAL, pct_all REAL, z_eur REAL, mean_eur REAL, sd_eur REAL, pmid TEXT, citation TEXT, method TEXT, metric TEXT, ancestry TEXT, note TEXT,
  pct_eur_raw REAL, flipped INTEGER, direction INTEGER);
CREATE TABLE prs_dist(id TEXT, lo REAL, hi REAL, bins TEXT);""")
for k, sid in enumerate(ids):
    if called[k] == 0:
        continue
    m = meta.get(sid, {})
    ref_e = total_ref[sel, k]
    you = user[k]
    raw_pe = float((ref_e < you).mean() * 100 + (ref_e == you).mean() * 50)
    if adjust:
        beta = np.linalg.lstsq(Xe, ref_e, rcond=None)[0]
        base = float(ref_e.mean())
        ref_e = ref_e - Xe @ beta + base          # what each founder would score at the average European position
        you = float(you - xu @ beta + base)
    mu, sd = float(ref_e.mean()), float(ref_e.std())
    pe = float((ref_e < you).mean() * 100 + (ref_e == you).mean() * 50)
    pa = float((total_ref[:, k] < user[k]).mean() * 100)
    z = (you - mu) / sd if sd > 0 else 0.0
    lo, hi = float(ref_e.min()), float(ref_e.max())
    span = hi - lo
    lo2, hi2 = min(lo, you) - span * 0.05, max(hi, you) + span * 0.05
    hist, _ = np.histogram(ref_e, bins=40, range=(lo2, hi2))
    # the evaluation to quote: of the score's own trait only (the catalog also
    # lists tests of a score against other traits), one in Europeans if there
    # is one, an odds or hazard ratio or AUROC before anything else; for a
    # score whose weights were reversed here, only a metric that carries no
    # sign (the published sign refers to the uncorrected file)
    metric = ""
    best = None
    own = m.get("Reported Trait", "")
    cols = ("Odds Ratio (OR)", "Hazard Ratio (HR)", "Area Under the Receiver-Operating Characteristic Curve (AUROC)",
            "Concordance Statistic (C-index)", "Beta", "Other Metric(s)")
    if sid in FLIP:
        cols = ("Area Under the Receiver-Operating Characteristic Curve (AUROC)", "Concordance Statistic (C-index)", "Other Metric(s)")
    for pr in perf.get(sid, []):
        if not same_trait(pr.get("Reported Trait", ""), own):
            continue
        anc = pss_anc.get(pr.get("PGS Sample Set (PSS)"), set())
        is_eur = any("European" in a for a in anc) and len(anc) == 1
        for rank, col in enumerate(cols):
            v = pr.get(col, "")
            if v:
                key = (0 if is_eur else 1, rank)
                if best is None or key < best[0]:
                    label = {"Area Under the Receiver-Operating Characteristic Curve (AUROC)": "AUROC",
                             "Concordance Statistic (C-index)": "C-index"}.get(col, col.split(" (")[0])
                    best = (key, "%s%s" % (tidy_metric(label, v[:90]), " in Europeans" if is_eur else ""))
                break
    if best:
        metric = best[1]
    cite = "%s" % (m.get("Publication (doi)", ""))
    trait = MIND[sid][0] if sid in MIND else m.get("Reported Trait", sid)
    cat = "Mind and behaviour" if sid in MIND else category(trait)
    note = MIND[sid][1] if sid in MIND else ""
    if sid in FLIP:
        note = "Weights reversed here: the published file runs the wrong way round (checked against independent scores). " + note
    direction = 0 if (sid in MIND or sid in NEUTRAL) else -1 if sid in PROTECTIVE else 1
    cur.execute("INSERT INTO prs VALUES(" + ",".join("?" * 22) + ")",
                (sid, trait, cat, m.get("PGS Name", ""), nvar[k], int(matched[k]), int(called[k]), float(you), pe, pa, float(z),
                 mu, sd, m.get("Publication (PMID)", ""), cite, m.get("PGS Development Method", ""), metric,
                 m.get("Ancestry Distribution (%) - Source of Variant Associations (GWAS)", ""), note, raw_pe, int(sid in FLIP), direction))
    cur.execute("INSERT INTO prs_dist VALUES(?,?,?,?)", (sid, lo2, hi2, json.dumps([int(x) for x in hist])))
    log("%s %-35s used %d/%d  called %.1f%%  percentile %.0f (z %+.2f; before ancestry adjustment %.0f)" % (
        sid, trait[:35], matched[k], nvar[k], 100 * called[k] / max(1, matched[k]), pe, z, raw_pe))
set_meta(con, "prs_reference", "the %d %s of 1000 Genomes" % (int(sel.sum()), GROUP_LABEL[group][0]) if group != "ALL"
         else "all %d people of 1000 Genomes (your ancestry is mixed, so treat these as rough)" % len(ref_names))
set_meta(con, "prs_group", group)
set_meta(con, "prs_group_label", GROUP_LABEL[group][0])
set_meta(con, "prs_group_short", GROUP_LABEL[group][1])
set_meta(con, "prs_portability", "" if group == "EUR" else
         "Most of these scores were built from studies of people of European ancestry, and they predict noticeably less well for "
         "other ancestries (Martin et al. 2019, Nature Genetics), so read them as rougher than they look.")
set_meta(con, "prs_adjusted", "1" if adjust else "0")
con.commit()
