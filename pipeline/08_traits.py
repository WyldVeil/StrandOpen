#!/usr/bin/env python3
"""Step 8: the curated traits and notable variants (traits_def.py).

rsIDs are placed on GRCh38 from data/trait_positions.json (resolved once
from Ensembl's REST service and shipped, so nothing is looked up online).
Each effect allele is checked against Ensembl's alleles for the rsID and,
when step 4 has run, against its 1000 Genomes European frequency: a
mismatch stops the step instead of reporting a flipped strand.

Table: trait."""
import json
import subprocess
from strandlib import *
from genotype import Genotyper
from traits_def import TRAITS

con = db()
cur = con.cursor()
G = Genotyper(con)
# the panel's rsIDs on GRCh38, resolved once from Ensembl (and NCBI dbSNP) and
# shipped with the program, so an import sends nothing anywhere
pos = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "trait_positions.json")))
missing = [r for t in TRAITS for r in t["rs"] if r not in pos]
if missing:
    log("not placed:", missing)

_fai = {}
with open(FA + ".fai") as f:
    for line in f:
        p = line.split("\t")
        _fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
_fa = open(FA, "rb")


def refbase(chrom, p):
    ln, off, lb, lw = _fai[chrom]
    z = p - 1
    _fa.seek(off + z // lb * lw + z % lb)
    return _fa.read(1).decode().upper()


# European frequencies from 1000 Genomes (step 4), for the strand check
kg = {}
if kg_site_files():
    pts = [(v["chrom"], v["start"]) for v in pos.values() if v["start"] == v["end"]]
    for line in kg_query(pts, "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF_EUR\n"):
        c, p, ref, alt, af = line.split("\t")
        try:
            kg[(c, int(p), alt)] = float(af)
            kg[(c, int(p), ref)] = 1 - float(af)
        except ValueError:
            pass

problems = []


def copies(rs, eff):
    """(effect allele copies, allele text, ploidy) or None when not called"""
    v = pos.get(rs)
    if not v:
        return None
    c = v["chrom"]
    ci = CHROM_IDX[c]
    if eff == "del":   # a deletion (CCR5 delta-32): look for it
        L = v["end"] - v["start"] + 1
        rows = cur.execute("SELECT zyg, ref, alt FROM variant WHERE chrom=? AND pos BETWEEN ? AND ? AND length(ref)-length(alt) BETWEEN ? AND ?",
                           (ci, v["start"] - 40, v["end"] + 5, L - 2, L + 2)).fetchall()
        if rows:
            z = rows[0][0]
            return (2 if z == 2 else 1), ("del/del" if z == 2 else "ins/del"), 2
        r = G.at(ci, v["start"])
        return (0, "ins/ins", 2) if r[0] == "ref" else None
    ref = refbase(c, v["start"])
    al = G.alleles(ci, v["start"], ref)
    if not al:
        return None
    if eff not in v["alleles"].split("/"):
        problems.append("%s: effect allele %s is not one of Ensembl's %s" % (rs, eff, v["alleles"]))
    return sum(1 for a in al if a == eff), "".join(sorted(al)), len(al)


cur.executescript("""
DROP TABLE IF EXISTS trait;
CREATE TABLE trait(id TEXT, cat TEXT, title TEXT, result TEXT, detail TEXT, conf TEXT, genotype TEXT,
  snps TEXT, about TEXT, source TEXT, status TEXT, copies INTEGER, ord INTEGER);
""")
for k, t in enumerate(TRAITS):
    calls = []
    for rs, eff, eur in zip(t["rs"], t["eff"], t["eur"]):
        c = copies(rs, eff)
        v = pos.get(rs, {})
        af = kg.get((v.get("chrom"), v.get("start"), eff))
        if af is not None and eur is not None and abs(af - eur) > 0.2 and abs(af - (1 - eur)) < abs(af - eur):
            problems.append("%s: %s is at %.2f in 1000G Europeans, expected about %.2f (strand?)" % (rs, eff, af, eur))
        calls.append(dict(rs=rs, eff=eff, call=c, chrom=v.get("chrom"), pos=v.get("start"), af_eur=af))
    gtxt = ", ".join("%s %s" % (x["rs"], x["call"][1] if x["call"] else "no call") for x in calls)
    status = "ok" if all(x["call"] for x in calls) else "partial" if any(x["call"] for x in calls) else "nocall"
    result = detail = ""
    n = None
    m = t["model"]
    if status == "nocall":
        result, detail = "Not called", "None of the positions this needs were confidently called."
    elif m in ("add", "dom", "rec"):
        n, txt, ploidy = calls[0]["call"]
        if ploidy == 1:
            idx = 2 if n == 1 else 0
        else:
            idx = n
        detail = t["texts"][idx]
        result = ["No copies", "One copy", "Two copies"][idx] + " of the %s allele" % t["eff"][0] if t["eff"][0] != "del" else \
            ["No delta-32", "One delta-32 copy", "Two delta-32 copies"][idx]
        if t.get("res"):
            result = t["res"][idx]
        if t["id"] == "lactase":
            result = "Likely lactose intolerant" if n == 0 else "Likely lactose tolerant"
        elif t["id"] == "eye":
            result = ["Brown eyes likely", "Brown, hazel or green likely", "Blue eyes likely"][idx]
        elif t["id"] == "earwax":
            result = "Dry earwax" if n == 2 else "Wet earwax"
        elif t["id"] == "actn3":
            result = ["RR (working)", "RX", "XX (no alpha-actinin-3)"][idx]
        elif t["id"] == "secretor":
            result = "Non-secretor" if n == 2 else "Secretor"
        elif t["id"] == "flush":
            result = "No flush variant" if n == 0 else "Alcohol flush"
        elif t["id"] == "bitter":
            result = ["Non-taster", "Taster", "Strong taster"][idx]
        elif t["id"] == "caffeine":
            result = ["Slow metaboliser", "Intermediate", "Fast metaboliser"][idx]
        elif t["id"] == "comt":
            result = ["Val/Val", "Val/Met", "Met/Met"][idx]
        elif t["id"] == "mthfr":
            result = ["677 CC", "677 CT", "677 TT"][idx]
    elif m == "mc1r":
        n = sum(x["call"][0] for x in calls if x["call"])
        result = ["No strong red-hair alleles", "One red-hair allele (carrier)", "Two red-hair alleles"][min(n, 2)]
        detail = ["Red hair is unlikely from MC1R.",
                  "Red hair is unlikely, but freckling and sunburn are more likely, and you could pass the allele on: two carriers can have a red-haired child.",
                  "Red or auburn hair (or a red beard) is likely, with fair skin that burns easily."][min(n, 2)]
    elif m == "apoe":
        c4 = calls[0]["call"][0] if calls[0]["call"] else None
        c2 = calls[1]["call"][0] if calls[1]["call"] else None
        if c4 is None or c2 is None:
            result, detail = "Not called", "Both APOE positions are needed."
        else:
            e3 = max(0, 2 - c4 - c2)
            al = ["e2"] * c2 + ["e3"] * e3 + ["e4"] * c4
            result = "/".join(al[:2])
            n = c4
            detail = {0: "No e4: average or lower Alzheimer's risk from APOE.",
                      1: "One e4: about two to three times the average risk of late-onset Alzheimer's disease; most people with one e4 never develop it.",
                      2: "Two e4: substantially raised risk of late-onset Alzheimer's disease (roughly eight to twelve times); worth discussing with a doctor or genetic counsellor."}[c4]
            if c2 == 2:
                detail += " Two e2: a risk of type III hyperlipoproteinaemia (it needs other triggers)."
    elif m == "abo":
        o = calls[0]["call"][0] if calls[0]["call"] else None
        b = calls[1]["call"][0] if calls[1]["call"] else 0
        if o is None:
            result, detail = "Not called", ""
        else:
            a = max(0, 2 - o - b)
            result = "AB" if a and b else "A" if a else "B" if b else "O"
            al = ["A"] * a + ["B"] * b + ["O"] * o
            detail = "Estimated ABO alleles: %s. Tag-SNP estimates are usually right but not always (rare alleles, A2): a blood test is definitive." % "/".join(al[:2])
    elif m == "a1at":
        z = calls[0]["call"][0] if calls[0]["call"] else 0
        s = calls[1]["call"][0] if calls[1]["call"] else 0
        g = "M" * max(0, 2 - z - s) + "S" * s + "Z" * z
        result = "Pi" + g
        detail = {"MM": "Normal alpha-1 antitrypsin.", "MS": "Carrier of S: levels slightly lower, no known risk.",
                  "MZ": "Carrier of Z: levels lower; avoid smoking, which adds lung risk; a mildly raised risk of liver disease.",
                  "SS": "Mildly lower levels, usually no disease.", "SZ": "Lower levels: a raised emphysema risk with smoking.",
                  "ZZ": "Alpha-1 antitrypsin deficiency: a clear risk of emphysema and liver disease; see a doctor."}.get(g, "")
    elif m == "lpa":
        a = calls[0]["call"][0] if calls[0]["call"] else 0
        b = calls[1]["call"][0] if calls[1]["call"] else 0
        n = a + b
        result = "Raised Lp(a) likely" if n else "No high-Lp(a) markers"
        detail = ("You carry %d high-Lp(a) marker allele(s). People with these usually have Lp(a) in the top fifth; a blood test confirms it."
                  % n) if n else "Neither common high-Lp(a) marker. Lp(a) can still be high for other genetic reasons (kringle IV repeats)."
    elif m == "celiac":
        nc = [x["call"][0] if x["call"] else None for x in calls]
        dq25, dq8 = (nc[0] or 0) > 0, (nc[1] or 0) > 0
        # DQ2.2 needs rs2395182 T and rs7775228 C; a G at rs4713586 means that haplotype may be DQ4 instead (phase unknown)
        dq22_maybe = (nc[2] or 0) > 0 and (nc[3] or 0) > 0
        dq22 = dq22_maybe and nc[4] == 0
        dq7 = (nc[5] or 0) > 0
        missing = any(x is None for x in nc[:5])
        if dq25 or dq8:
            result = " and ".join(x for x, y in (("DQ2.5", dq25), ("DQ8", dq8)) if y) + " present"
            detail = ("You carry an HLA type most people with coeliac disease have (as do about a third of all Europeans); "
                      "it is necessary but far from sufficient: most carriers never develop it.")
        elif dq22:
            result = "DQ2.2 present"
            detail = ("No DQ2.5 or DQ8, but DQ2.2, which carries a lower risk and is found in a small share of people with coeliac "
                      "disease. Coeliac disease is unlikely but not ruled out.")
        elif dq22_maybe:
            result = "DQ2.2 possible"
            detail = ("No DQ2.5 or DQ8. The tags cannot tell DQ2.2 from DQ4 here, so coeliac disease is unlikely but not ruled out.")
        elif missing:
            result = "Not fully tested"
            detail = "Some of the tag positions were not read, so coeliac disease cannot be ruled out from these genes."
        else:
            result = "None of the coeliac HLA types"
            seqd = (cur.execute("SELECT value FROM meta WHERE key='call_model'").fetchone() or ["gvcf"])[0] != "chip"
            detail = ("Coeliac disease is very unlikely without DQ2.5, DQ8 or DQ2.2: doctors use a negative HLA test to rule it out."
                      + (" (You carry DQ7, which on its own adds little risk.)" if dq7 else "")
                      + (" Sequencing pipelines differ in how they read the HLA genes, so from sequencing this is a little less certain "
                         "than from a dedicated HLA test." if seqd else ""))
    cur.execute("INSERT INTO trait VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (t["id"], t["cat"], t["title"], result, detail, t["conf"], gtxt, json.dumps(calls), t["about"], t["source"], status, n, k))
con.commit()
for p in problems:
    log("CHECK:", p)
if any("not one of" in p for p in problems):
    raise SystemExit("allele check failed; fix traits_def.py")
if problems and any("strand" in p for p in problems):
    raise SystemExit("strand check failed; fix traits_def.py")
log("traits", len(TRAITS), "1000G checks" if kg else "(no 1000G frequencies yet: strand check skipped)")
for r in cur.execute("SELECT title, result, genotype FROM trait ORDER BY ord"):
    print("  %-45s %-32s %s" % r)
