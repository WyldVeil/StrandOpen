#!/usr/bin/env python3
"""Step 17: the overview's highlights, the few results most worth a look.

Drawn from every other step, strongest first: medicine safety (HLA risk
alleles), ClinVar pathogenic variants, drug-metabolism genes that change a
dose, well-known health variants, remarkable variants carried, outlying
polygenic scores and the count of rare disruptive changes. Each row names the app page it opens
(`target`, the same names as --enter)."""
import re
from strandlib import *

con = db()
cur = con.cursor()


def has(t):
    return cur.execute("SELECT 1 FROM sqlite_master WHERE name=?", (t,)).fetchone() is not None


rows = []   # (tone, title, detail, badge, target)
# guideline advice that keeps the usual dose, or gives none for this type
SAME = re.compile(r"(initiate (therapy )?(with )?(the )?(recommended|standard|normal)|initiate standard|initiate with a dose of|"
                  r"use (\S+ )?(label recommended|per standard|at standard)|if considering \S+,? use at standard|no action|"
                  r"no recommendation|no dose adjustment|no reason to|there is no need|the guideline does not provide|"
                  r"may (alter|affect) systemic|in the event of side effects)")
# what to do once the usual dose fails or causes side effects is not a change
FALLBACK = re.compile(r"\bif |\bin the case of|\bin the event of")
CHANGE = re.compile("|".join((
    r"\bavoid\b(?! (co-medication|concomitant|moderate|strong|potent)\b)", r"contra-?indicated", r"not recommended",
    r"do not (use|select|prescribe)", r"\b(choose|select|selecting|prescribe|consider) (an |a |another )?(alternat|other|clinically appropriate)",
    r"alternative (agent|treatment|drug|medic|therap|antidepressant|statin)", r"an alternative (can|should|may) be",
    r"not (predominantly )?metaboli[sz]ed by", r"aromatase inhibitor",
    r"\d+ ?% (of (the )?(standard|normal|usual|recommended|starting|initial)|less than|lower than)",
    r"\d+(\.\d+)? (to \d+(\.\d+)? )?times (the )?(standard|normal|higher|recommended)", r"\d(\.\d+)?x higher",
    r"(half|double) (of )?the (standard|normal|usual)",
    r"reduc(e|ed|tion of|tion in|ing) (the )?(recommended |normal |standard )?(starting |initial |maintenance |daily )?dos",
    r"(lower|lowest|higher|reduced|decreased) (recommended )?(starting |initial |maintenance |daily )?dose(?! compared to poor)",
    r"increase (the )?(starting |initial )?(daily )?dose", r"smaller steps", r"do not exceed", r"doses exceeding", r"no more than",
    r"prescribe (\u2264|<) ?\d", r"limit (the )?dose", r"use the maximum dose")))


def changes(rec):
    """True when the advice for this type is a different dose or another
    medicine from the start: not the usual dose, no advice, or a change
    only once the usual dose fails"""
    t = (rec or "").replace("&quot;", "\"").lower().lstrip("\"'.-: \n")
    if SAME.match(t):
        return False
    m = FALLBACK.search(t, 1)
    return CHANGE.search(t[:m.start()] if m else t) is not None


def actionable_drugs(gene):
    out = []
    for drug, rec in cur.execute("SELECT drug, recommendation FROM pgx_drug WHERE gene=? AND matched=1", (gene,)):
        if changes(rec) and drug not in out:
            out.append(drug)
    return out


if has("pgx_gene"):
    HLA = {"*31:01": ("carbamazepine", "a much higher chance of serious skin reactions to carbamazepine; guidelines advise another medicine where possible"),
           "*15:02": ("carbamazepine and phenytoin", "a much higher chance of serious skin reactions to carbamazepine and phenytoin"),
           "*57:01": ("abacavir", "a high chance (about half of carriers) of a serious hypersensitivity reaction to abacavir; guidelines say do not use it"),
           "*58:01": ("allopurinol", "a much higher chance of severe skin reactions to allopurinol")}
    for gene, dip, ph in cur.execute("SELECT gene, diplotype, phenotype FROM pgx_gene WHERE gene IN ('HLA-A','HLA-B') AND called=1").fetchall():
        for al, (drug, what) in HLA.items():
            if (al + " positive") in (ph or ""):
                rows.append(("bad", "%s%s: %s warning" % (gene, al, drug), "You carry %s%s, which brings %s." % (gene, al, what), "Medicine", "pgx:0"))
    KEY = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "DPYD", "TPMT", "NUDT15", "UGT1A1", "VKORC1", "G6PD", "CYP2B6", "NAT2"]
    for gene in KEY:
        r = cur.execute("SELECT diplotype, phenotype FROM pgx_gene WHERE gene=? AND called=1", (gene,)).fetchone()
        if not r or not r[1]:
            continue
        ph = r[1]
        if ph.lower().startswith(("normal", "indeterminate", "n/a")) or "uncertain" in ph.lower():
            continue
        drugs = actionable_drugs(gene)
        if not drugs:
            continue
        rows.append(("warn", "%s: %s" % (gene, ph.lower()), "%s. Guidelines suggest a different dose or choice for %d medicine%s: %s%s." % (
            r[0], len(drugs), "" if len(drugs) == 1 else "s", ", ".join(drugs[:4]),
            " and more" if len(drugs) > 4 else ""), "Medicine", "pgx:0"))

_cm = cur.execute("SELECT value FROM meta WHERE key='call_model'").fetchone()
CHIP = bool(_cm) and _cm[0] == "chip"
for vid, gene, hgvs, zyg, dn in cur.execute("SELECT id, gene, hgvs_p, zyg, clin_dn FROM variant WHERE clin IN (5,6) ORDER BY clin DESC, clin_stars DESC").fetchall():
    cond = (dn or "").replace("|", ";").split(";")
    cond = sorted({c.strip() for c in cond if c.strip() and c.strip().lower() not in ("not provided", "not specified", "see cases")
                   and "related disorder" not in c.lower() and "inborn genetic" not in c.lower()}, key=len)
    row = ("bad", "%s %s" % (gene or "intergenic", hgvs or ""), "%s of a variant ClinVar calls pathogenic%s." % (
        "One copy" if zyg == 1 else "Two copies" if zyg == 2 else "A single copy", ", linked to %s" % cond[0].strip().lower() if cond else "")
        + (" A chip call: rare ones can be wrong, so a clinical test may make sense before acting on it." if CHIP else ""),
        "Health", "health:0")
    rows.insert(sum(1 for r in rows if r[0] == "bad"), row)   # after the HLA rows and earlier ClinVar rows, before any warning

if CHIP:
    n = cur.execute("SELECT count(*) FROM variant WHERE clin IN (105,106)").fetchone()[0]
    if n:
        rows.append(("info", "%d rare chip call%s set aside" % (n, "" if n == 1 else "s"),
                     "The chip reported %s ClinVar lists as disease-causing, but chip calls this rare are usually wrong, so %s set aside as "
                     "unconfirmed. A clinical test may make sense before reading anything into %s." % (
                         "a rare variant" if n == 1 else "%d rare variants" % n, "it is" if n == 1 else "they are", "it" if n == 1 else "them"),
                     "Chip", "health:0"))

if has("trait"):
    # a carrier of recessive HFE C282Y, PiMS, PiSS and an uncalled result have no real effect on risk
    for tid, title, result, detail in cur.execute("SELECT id, title, result, detail FROM trait WHERE cat='Health' "
                                                  "AND id IN ('apoe','fvl','pt','hfe282','a1at','lpa')").fetchall():
        result = result or ""
        bad = (tid == "apoe" and "e4" in result) or (tid in ("fvl", "pt", "lpa") and not result.lower().startswith(("not", "no "))) or \
              (tid == "hfe282" and result.startswith("Homozygous")) or (tid == "a1at" and result.startswith("Pi") and "Z" in result)
        if bad:
            rows.append(("warn", "%s: %s" % (title, result), "%s The clinical page explains it." % detail if detail else
                         "A well-known variant with a real effect on health risk; the clinical page explains it.", "Health", "health:3"))

if has("remarkable"):
    # strongest evidence first, then the rarest; disputed claims are left to the Remarkable page
    for title, gene, st, af in cur.execute("SELECT title, gene, status, af FROM remarkable WHERE status=2 AND (af IS NULL OR af < 0.05) "
                                           "AND coalesce(conf, '') != 'disputed' ORDER BY CASE conf WHEN 'strong' THEN 0 "
                                           "WHEN 'moderate' THEN 1 ELSE 2 END, coalesce(af, 0)").fetchall()[:3]:
        how = "rare %s variant" % gene if af is None or af < 0.01 else "%s variant (about %d in 100 people carry it)" % (gene, max(1, round(200 * af * (1 - af) + 100 * af * af)))
        rows.append(("accent", title, "You carry this %s%s." % (how, "" if st == 2 else "; one copy, the full effect needs two"), "Remarkable", "rare:0"))

if has("prs"):
    cols = [r[1] for r in cur.execute("PRAGMA table_info(prs)")]
    health = " AND category != 'Mind and behaviour'" if "note" in cols else ""
    who = (cur.execute("SELECT value FROM meta WHERE key='prs_group_label'").fetchone() or ["Europeans"])[0]
    for trait, pct in cur.execute("SELECT trait, pct_eur FROM prs WHERE (pct_eur >= 95 OR pct_eur <= 5)" + health +
                                  " ORDER BY abs(pct_eur - 50) DESC").fetchall()[:3]:
        hi = pct >= 50
        rows.append(("info", "%s: %s polygenic score" % (trait, "high" if hi else "low"),
                     ("Higher than %d%% of %s in the reference. A score shifts the odds; it does not decide them." % (min(99, round(pct)), who)) if hi else
                     ("Lower than %d%% of %s in the reference." % (min(99, round(100 - pct)), who)), "Score", "prs"))

if has("rare_var"):
    n = cur.execute("SELECT count(*) FROM rare_var WHERE verdict=3 AND suspect=0").fetchone()[0]
    if n:
        rows.append(("info", "%d rare changes likely disrupt a protein" % n, "Each is explained in plain words. Everyone carries some; nearly all are harmless.", "Rare", "rare:1"))

cur.executescript("DROP TABLE IF EXISTS highlight; CREATE TABLE highlight(ord INTEGER, tone TEXT, title TEXT, detail TEXT, badge TEXT, target TEXT);")
for i, r in enumerate(rows):
    cur.execute("INSERT INTO highlight VALUES(?,?,?,?,?,?)", (i,) + r)
    log(r[0], r[1], "|", r[2])
con.commit()
