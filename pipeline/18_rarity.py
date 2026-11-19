#!/usr/bin/env python3
"""Step 18: how rare your variants are in the world.

The only population data on this machine is 1000 Genomes: 2,590 unrelated
people in five continental groups. A world figure is estimated by weighting
each group by the share of the world's people it stands for (rough, and
the Middle East, North Africa and Central Asia are split between the
nearest groups). For each variant:
  carriers   share of people with at least one copy, 1 - (1 - p)^2 per group
  genotype   share with your genotype (two copies: p^2)
A variant missing from the panel was carried by at most one of the 2,590
(the panel leaves out singletons), so it is "fewer than 1 in 2,590".

Table rarity(kind, ...):
  protein   your rarest changes to proteins (quality-checked)
  homo      rare two-copy genotypes: variants in or near genes you carry
            twice that fewer than 1 in 1,000 people are expected to
  ancestry  common somewhere in the world, rare in your own ancestry group

Needs step 16's rare_var and step 5b's mapping quality (variant.mq), and
runs those steps first when they are missing. Without 1000 Genomes
frequencies in the variant table it records meta rare_has_af = '0' and
leaves the table empty.
"""
import math
import subprocess
import sys
from strandlib import *

HERE = os.path.dirname(os.path.abspath(__file__))


def run_step(script):
    subprocess.run([sys.executable, os.path.join(HERE, script)], check=True)


# a lone re-run of step 5 drops variant.mq (added by 5b) and renumbers the
# variants, so rare_var is rebuilt then too
con = db()
has_mq = "mq" in [r[1] for r in con.execute("PRAGMA table_info(variant)")] and \
    con.execute("SELECT 1 FROM variant WHERE mq IS NOT NULL LIMIT 1").fetchone() is not None
has_rare = con.execute("SELECT 1 FROM sqlite_master WHERE name='rare_var'").fetchone() is not None
con.close()
if not has_mq:
    log("variant.mq is missing (step 5 ran without step 5b): running steps 5b and 16 first")
    run_step("05b_mapq.py")
    run_step("16_rare.py")
elif not has_rare:
    log("rare_var is missing: running step 16 first")
    run_step("16_rare.py")

con = db()
cur = con.cursor()
POPS = ("afr", "eur", "eas", "sas", "amr")
# billions: Africa + diaspora; Europe + European-descended + half of MENA;
# East and Southeast Asia; South and Central Asia + half of MENA; admixed Americas
PEOPLE = {"afr": 1.35, "eur": 1.55, "eas": 2.30, "sas": 2.35, "amr": 0.40}
TOT = sum(PEOPLE.values())
W = {k: v / TOT for k, v in PEOPLE.items()}
NAMES = {"afr": "Africa", "eur": "Europe", "eas": "East Asia", "sas": "South Asia", "amr": "the Americas"}

# the person's main ancestry group
anc = cur.execute("SELECT pop, frac FROM anc_admix WHERE level='super' ORDER BY frac DESC LIMIT 1").fetchone() if \
    cur.execute("SELECT 1 FROM sqlite_master WHERE name='anc_admix'").fetchone() else None
MAIN = {"EUR": "eur", "AFR": "afr", "EAS": "eas", "SAS": "sas", "AMR": "amr"}.get((anc[0] if anc else "EUR").upper(), "eur")
log("main ancestry group", MAIN)

plain = {}
for fn in sorted(os.listdir(os.path.join(HERE, "..", "data"))):
    if fn.startswith("gene_plain") and fn.endswith(".tsv"):
        for line in open(os.path.join(HERE, "..", "data", fn), encoding="ascii"):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1]:
                plain[p[0]] = p[1]
gname = dict(cur.execute("SELECT symbol, name FROM gene WHERE biotype='protein_coding'").fetchall())


def carriers(fr):
    return sum(W[k] * (1 - (1 - (fr[k] or 0)) ** 2) for k in POPS)


def genotype2(fr):
    return sum(W[k] * (fr[k] or 0) ** 2 for k in POPS)


def where_common(fr):
    k = max(POPS, key=lambda q: fr[q] or 0)
    return k


cur.executescript("""DROP TABLE IF EXISTS rarity;
CREATE TABLE rarity(kind TEXT, rank INTEGER, vid INTEGER, gene TEXT, change TEXT, csq TEXT, plain TEXT, zyg INTEGER, af REAL,
  af_afr REAL, af_eur REAL, af_eas REAL, af_sas REAL, af_amr REAL, world REAL, note TEXT);""")
set_meta(con, "rarity_main_group", NAMES[MAIN])
set_meta(con, "rarity_weights", ", ".join("%s %.0f%%" % (NAMES[k], 100 * W[k]) for k in POPS))

# step 5 builds the variant table without frequencies when the 1000 Genomes
# site files (step 4) are missing: every list here would then call common
# variants "never seen", so they stay empty
has_af = cur.execute("SELECT 1 FROM variant WHERE af IS NOT NULL LIMIT 1").fetchone() is not None
set_meta(con, "rare_has_af", "1" if has_af else "0")
if not has_af:
    cur.execute("DELETE FROM meta WHERE key IN ('rarity_under_1000', 'rarity_homo_under_1000', 'rarity_ancestry_n')")
    con.commit()
    log("no 1000 Genomes frequencies in the variant table (run steps 4 and 5): lists skipped")
    sys.exit(0)


def put(kind, rank, r, world, note):
    vid, gene, hgvs, csq, zyg, af, fr = r
    cur.execute("INSERT INTO rarity VALUES(" + ",".join("?" * 16) + ")",
                (kind, rank, vid, gene or "", hgvs or "", csq or "", plain.get(gene or "", "") or gname.get(gene or "", ""), zyg, af,
                 fr["afr"], fr["eur"], fr["eas"], fr["sas"], fr["amr"], world, note))


def rows(sql, args=()):
    for vid, gene, hgvs, csq, zyg, af, a1, a2, a3, a4, a5 in cur.execute(sql, args).fetchall():
        yield (vid, gene, hgvs, csq, zyg, af, {"afr": a1, "eur": a2, "eas": a3, "sas": a4, "amr": a5})


COLS = "v.id, v.gene, v.hgvs_p, v.csq, v.zyg, v.af, v.af_afr, v.af_eur, v.af_eas, v.af_sas, v.af_amr"

# ---- rarest protein changes: the quality-checked rare protein changes of step 16
prot = []
change = {}
for r in rows("SELECT " + COLS + " FROM rare_var rv JOIN variant v ON v.id=rv.vid WHERE rv.suspect=0"):
    world = carriers(r[6]) if r[5] is not None else None
    prot.append((world if world is not None else -1.0, r))
for vid, chg in cur.execute("SELECT vid, change FROM rare_var WHERE suspect=0"):
    change[vid] = chg
prot.sort(key=lambda t: t[0])
for i, (world, r) in enumerate(prot):
    note = change.get(r[0], "")
    if world >= 0:
        k = where_common(r[6])
        note += " Most common in %s (%s of people there carry it)." % (NAMES[k], "%.2g%%" % (100 * (1 - (1 - r[6][k]) ** 2)))
    put("protein", i, r, None if world < 0 else world, note)
log("protein", len(prot))

# ---- rare two-copy genotypes (variants called twice that under 1 in 1,000 people have twice)
homo = []
# clean single-letter calls only: indel lengths in repeats are where callers stumble
for r in rows("SELECT " + COLS + " FROM variant v WHERE v.zyg=2 AND v.kind=0 AND v.pass=1 AND v.af >= 0.001 AND v.af < 0.15 "
              "AND coalesce(v.dp, 30) BETWEEN 10 AND 80 AND coalesce(v.gq, 99) >= 40 AND coalesce(v.ad_ref, 0) * 10 <= coalesce(v.dp, 30) "
              "AND coalesce(v.mq, 60) >= 50 AND coalesce(v.mq0f, 0) < 0.05"):
    homo.append((genotype2(r[6]), r))
CODING = ("missense", "stop_gained", "frameshift", "splice_donor", "splice_acceptor", "start_lost", "stop_lost", "inframe_deletion",
          "inframe_insertion", "synonymous", "splice_region")


def tier(r):
    csq = r[3] or ""
    return 0 if csq in CODING else 1 if csq in ("5_prime_utr", "3_prime_utr") else 2 if r[1] else 3


# the card counts the same set the list is drawn from: rare, and in or near a gene
homo = [t for t in homo if t[0] < 0.001 and tier(t[1]) < 3]
homo.sort(key=lambda t: (tier(t[1]), t[0]))
nhomo_rare = len(homo)


def one_in(x):
    if x < 1e-5:
        return "fewer than 1 in 100,000"
    n = 1 / x
    n = round(n, -int(math.floor(math.log10(n))) + 1) if n >= 100 else round(n)
    return "about 1 in %s" % format(int(n), ",")


seen = {}
i = 0
for g, r in homo:
    key = r[1] or ("%d" % r[0])
    if key in seen:
        seen[key] += 1
        continue
    seen[key] = 0
    put("homo", i, r, g, "%s people are expected to have two copies" % one_in(g).capitalize())
    i += 1
    if i >= 200:
        break
set_meta(con, "rarity_homo_under_1000", str(nhomo_rare))
log("two-copy genotypes under 1 in 1,000, in or near genes:", nhomo_rare)

# ---- unusual for your ancestry: rare in your main group, common elsewhere
anc_rows = []
for r in rows("SELECT " + COLS + " FROM variant v WHERE v.pass=1 AND v.af_%s < 0.01 AND max(v.af_afr, v.af_eas, v.af_sas, v.af_amr, v.af_eur) >= 0.2 "
              "AND (v.ad_alt IS NULL OR v.ad_alt*5 >= (v.ad_ref+v.ad_alt)) AND coalesce(v.mq, 60) >= 50 AND coalesce(v.mq0f, 0) < 0.05 "
              "AND coalesce(v.gq, 99) >= 30" % MAIN):
    fr = r[6]
    k = max((q for q in POPS if q != MAIN), key=lambda q: fr[q] or 0)
    coding = r[3] in ("missense", "stop_gained", "frameshift", "splice_donor", "splice_acceptor", "synonymous", "inframe_deletion", "inframe_insertion")
    anc_rows.append(((0 if coding else 1), -(fr[k] - fr[MAIN]), k, r))
anc_rows.sort(key=lambda t: (t[0], t[1]))
seen = set()
i = 0
for c, d, k, r in anc_rows:
    if r[1] and r[1] in seen:
        continue
    seen.add(r[1])
    i += 1
    if i > 200:
        break
    put("ancestry", i, r, carriers(r[6]),
        "%s of people in %s carry it, %s in %s" % ("%.0f%%" % (100 * (1 - (1 - r[6][k]) ** 2)), NAMES[k],
                                                   "%.2g%%" % (100 * (1 - (1 - r[6][MAIN]) ** 2)) if r[6][MAIN] else "none seen", NAMES[MAIN]))
set_meta(con, "rarity_ancestry_n", str(len(anc_rows)))
log("unusual for ancestry", len(anc_rows))

# ---- headline counts
n_world_1000 = 0
for r in rows("SELECT " + COLS + " FROM variant v WHERE v.pass=1 AND v.af IS NOT NULL AND v.af < 0.002"):
    if carriers(r[6]) < 0.001:
        n_world_1000 += 1
set_meta(con, "rarity_under_1000", str(n_world_1000))
log("variants carried by under 1 in 1,000 people:", n_world_1000)
con.commit()
log("done")
