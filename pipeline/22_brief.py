#!/usr/bin/env python3
"""Step 22: your genome in brief.

The findings worth knowing first, gathered from every other step and written
in plain words, each with what it means and where the details are. It reads
only what the earlier steps stored (run it last); a section whose step did
not run is left out. Wording rules: measured, never alarmist, never advice
to act alone; a score or a rare change "can" mean something, "may" matter.

Table brief(section, ord, tone, title, text, meaning, target)."""
import re
from strandlib import *

con = db()
cur = con.cursor()
meta = dict(cur.execute("SELECT key, value FROM meta").fetchall())


def has(t):
    return cur.execute("SELECT 1 FROM sqlite_master WHERE name=?", (t,)).fetchone() is not None


def first_sentences(text, limit=260):
    """the opening of a longer text, cut at a sentence end"""
    text = (text or "").strip()
    out = ""
    for s in re.split(r"(?<=[.!?])\s+", text):
        if out and len(out) + len(s) > limit:
            break
        out = (out + " " + s).strip()
    return out or text[:limit]


def num(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def pct(x):
    return "%.0f%%" % (100 * x) if x >= 0.01 else "under 1%"


items = []


def add(section, tone, title, text, meaning, target):
    items.append((section, len(items), tone, title, text, meaning, target))


chip = meta.get("call_model") == "chip"
sex = meta.get("sex", "")

# ---------------------------------------------------------------- remarkable and rare
if has("remarkable"):
    rows = cur.execute("SELECT title, found, effect, rarity, af, cat FROM remarkable WHERE status=2 ORDER BY coalesce(af, 1)").fetchall()
    shown = 0
    for title, found, effect, rarity, af, cat in rows:
        if shown >= 3 or (af is not None and af > 0.05 and shown):
            break
        add("Rare and remarkable", "accent", title,
            "You carry it: %s. How common: %s." % ((found or "").split(";")[0].strip(), (rarity or "").rstrip(".")),
            first_sentences(effect), "rare:0")
        shown += 1
if has("rarity"):
    n1000 = int(num(meta.get("rarity_under_1000")))
    if n1000:
        add("Rare and remarkable", "info", "%s of your variants are carried by fewer than 1 in 1,000 people" % format(n1000, ","),
            "Estimated from the 1000 Genomes frequencies on every continent, weighted by how many people live there.",
            "Everyone carries thousands of rare variants, and nearly all are harmless: they are what makes a genome one person's own. "
            "The rarest of yours are listed with a plain description of each gene.", "rarity:0")
    r = cur.execute("SELECT gene, change, plain, note FROM rarity WHERE kind='homo' ORDER BY rank LIMIT 1").fetchone()
    if r:
        add("Rare and remarkable", "info", "Your rarest two-copy genotype: %s %s" % (r[0], r[1]),
            "%s. %s" % ((r[3] or "").rstrip("."), ("The gene: " + r[2].rstrip(".") + ".") if r[2] else ""),
            "Having the same rare letter on both copies happens when each parent passed on the same rare variant. Most such changes, "
            "like this one, have no known effect.", "rarity:1")

# ---------------------------------------------------------------- medicines
if has("highlight"):
    for tone, title, detail in cur.execute("SELECT tone, title, detail FROM highlight WHERE target LIKE 'pgx%' ORDER BY ord").fetchall():
        add("Medicines", tone, title, detail,
            "Worth telling a doctor or pharmacist before any of these medicines is prescribed. The Medications page lists each medicine "
            "and what the guidelines (CPIC, DPWG) advise.", "pgx:0")

# ---------------------------------------------------------------- health
# inheritance of the condition itself where it can be matched (HBB has a rare
# dominant condition but sickle cell disease is recessive), else of the gene;
# the wording follows the Clinical findings page
PLACEHOLDER = {"not provided", "not specified", "see cases"}
STOP = {"disease", "disorder", "syndrome", "type", "deficiency", "susceptibility", "to", "of", "the", "and", "with", "due", "hereditary", "familial"}
gdis = {}
if has("gene_disease"):
    for gid, dis, inh in cur.execute("SELECT gene, disease, inheritance FROM gene_disease"):
        gdis.setdefault(gid, []).append(((dis or "").lower(), {x.strip().upper() for x in (inh or "").split(",") if x.strip()}))


def words(t):
    return set(re.findall(r"[a-z0-9]+", t.lower())) - STOP


def inheritance(gid, conds):
    rows = gdis.get(gid, [])
    got = set()
    for c in conds:
        cw = words(c)
        for dis, inh in rows:
            dw = words(dis)
            if cw and dw and len(cw & dw) / min(len(cw), len(dw)) >= 0.6:
                got |= inh
    if not got:
        for dis, inh in rows:
            got |= inh
    return got


def meaning_for(inh, zyg):
    ad, ar = "AD" in inh, "AR" in inh
    xl = any(x.startswith("XL") for x in inh)
    if zyg in (2, 3):
        if zyg == 3 and xl:
            return "One copy on your single X: for an X-linked condition this can be enough to cause it. Worth discussing with a doctor."
        if ar:
            return "Two copies: for a recessive condition this can be enough to cause it. Worth discussing with a doctor or genetic counsellor."
        return "Two copies of a variant linked to a condition. Worth discussing with a doctor or genetic counsellor."
    if ar and not ad:
        return ("You are a carrier: with one copy of a recessive condition you are usually unaffected. It can matter for children only "
                "if a partner carries a change in the same gene, which a carrier test can check.")
    if ad and not ar:
        return ("A dominant condition, where one copy can be enough, though many people with one copy have mild signs or none. Worth "
                "discussing with a doctor or genetic counsellor; a clinical test may make sense to confirm it.")
    if ad and ar:
        return "The gene is linked to both dominant and recessive conditions, so one copy may matter depending on which. A genetic counsellor can say."
    return "What it means depends on the condition and on family history; a genetic counsellor can help read it."


if has("variant"):
    for gene, gid, hgvs, zyg, dn in cur.execute("SELECT gene, gene_id, hgvs_p, zyg, clin_dn FROM variant WHERE clin IN (5,6) "
                                                 "ORDER BY clin DESC, clin_stars DESC").fetchall():
        conds = [c.strip().rstrip(".") for c in (dn or "").replace("|", ";").split(";") if c.strip() and c.strip().lower() not in PLACEHOLDER]
        # name the condition by the shortest name of at least two words (single words are often protein names)
        named = [c for c in conds if len(c.split()) >= 2] or conds
        cond = min(named, key=len) if named else ""
        if cond.isupper():      # ClinVar sometimes shouts: HEMOGLOBIN S -> Hemoglobin S
            w = [x if len(x) <= 2 else x.lower() for x in cond.split()]
            cond = " ".join([w[0][:1].upper() + w[0][1:]] + w[1:])
        inh = inheritance(gid, conds)
        carrier = zyg == 1 and "AR" in inh and "AD" not in inh
        add("Health", "warn" if carrier else "bad", ("%s %s" % (gene or "", hgvs or "")).strip(),
            "%s of a variant ClinVar lists as disease-causing%s." % ("One copy" if zyg == 1 else "Two copies" if zyg == 2 else "A single copy",
                                                                 ", linked to " + cond if cond else ""),
            meaning_for(inh, zyg), "health:0")
    n_set = int(num(meta.get("chip_rare_clinvar")))
    if chip and n_set:
        add("Health", "info", "%d rare chip call%s set aside" % (n_set, "" if n_set == 1 else "s"),
            "The chip reported rare variants ClinVar lists as disease-causing.",
            "Chip calls this rare are usually wrong (most were not there when the same people were sequenced, Weedon et al. 2021), so "
            "they are kept apart as unconfirmed. A clinical test may make sense before reading anything into them.", "health:0")
if has("prs"):
    ext = cur.execute("SELECT trait, category, pct_eur FROM prs WHERE category != 'Mind and behaviour' AND (pct_eur >= 95 OR pct_eur <= 5) "
                      "ORDER BY abs(pct_eur - 50) DESC").fetchall()
    immune = [t for t, c, p in ext if c == "Immune"]
    for trait, cat, p in ext[:5]:   # the most extreme first
        hi = p >= 50
        text = "Your polygenic score is %s than %d%% of %s in the reference." % (
            "higher" if hi else "lower", min(99, round(p if hi else 100 - p)), meta.get("prs_group_label", "Europeans"))
        if meta.get("prs_portability"):
            text += " " + meta["prs_portability"]
        meaning = ("A polygenic score adds up thousands of small effects. It shifts the odds a little, and most people with a high "
                   "score never develop the condition (and many with low scores do).")
        if cat == "Immune" and len(immune) >= 2:
            meaning += " Immune scores lean on the same HLA region, so they tend to move together."
        add("Health", "info", "%s: %s score" % (trait, "high" if hi else "low"), text, meaning, "prs")

# ---------------------------------------------------------------- diet: only what a guideline or strong evidence acts on
if has("diet"):
    for topic, grade, result, says, means in cur.execute("SELECT topic, grade, result, says, means FROM diet WHERE ord >= 0 ORDER BY ord").fetchall():
        act = ((topic.startswith("Iron") and ("Two copies" in result or "C282Y and H63D" in result)) or topic.startswith("Fava")
               or (topic.startswith("Alcohol") and result != "No flush variant"))
        if act:
            add("Diet", "warn", "%s: %s" % (topic, result), says, means, "diet")

# ---------------------------------------------------------------- ancestry
if has("anc_admix"):
    top = cur.execute("SELECT label, frac FROM anc_admix WHERE level='pop' AND frac >= 0.03 ORDER BY frac DESC LIMIT 5").fetchall()
    sup = cur.execute("SELECT label, frac FROM anc_admix WHERE level='super' AND frac >= 0.01 ORDER BY frac DESC").fetchall()
    if top:
        add("Ancestry", "accent", ", ".join("%s %s" % (pct(f), l) for l, f in sup) if sup else "Your ancestry",
            "Closest to: " + "; ".join("%s %s" % (l, pct(f)) for l, f in top) + ".",
            "Compared with 26 populations of the 1000 Genomes Project. Neighbouring groups (%s, say) share much of their history, so "
            "the split among them is uncertain; the continental picture is solid." % {
                "European": "British and Northwest European", "African": "Yoruba and Esan", "East Asian": "Han Chinese from the north and the south",
                "South Asian": "Gujarati and Telugu"}.get(sup[0][0] if sup else "", "neighbouring populations"), "ancestry:0")
if meta.get("roh_summary"):
    long5 = num(meta.get("roh_total_5mb"))
    add("Ancestry", "info", "Did your parents share recent ancestors? %s" % ("Possibly" if long5 >= 40e6 else "No sign of it"),
        meta["roh_summary"],
        "When both parents descend from the same recent ancestor, a child inherits long stretches where both copies match. Short "
        "matching stretches are in everyone; only long ones point to related parents.", "ancestry:0")
if meta.get("archaic_frac"):
    f = num(meta["archaic_frac"]) * 100
    # a fun comparison with the published average from a similar stretch-finding method
    # (Sankararaman et al. 2014, Nature 507:354: 1.15% in Europeans, 1.38% in East Asians)
    title, cmp = "About %.1f%% Neanderthal" % f, ""
    if has("anc_admix"):
        r = cur.execute("SELECT pop, frac FROM anc_admix WHERE level='super' ORDER BY frac DESC LIMIT 1").fetchone()
        avg, who = {"EUR": (1.15, "Europeans"), "EAS": (1.38, "East Asians")}.get(r[0] if r and r[1] >= 0.8 else "", (0, ""))
        if avg:
            rel = (f - avg) / avg * 100
            if abs(rel) < 10:
                title += ": about average"
                cmp = " The average for %s measured by a similar method is %.2f%% (Sankararaman et al. 2014)." % (who, avg)
            else:
                title += ": %.0f%% %s than the average %s" % (abs(rel), "more" if rel > 0 else "less", who[:-1])
                cmp = " The average for %s measured by a similar method is %.2f%% (Sankararaman et al. 2014)." % (who, avg)
    add("Ancestry", "accent", title,
        "%s stretches of your genome came from archaic humans who mixed with your ancestors over 40,000 years ago.%s" %
        (format(int(num(meta.get("archaic_segments"))), ","), cmp),
        "Everyone whose ancestors lived outside Africa carries some; genome-wide studies put it at about 2 percent, and stretch-finding "
        "methods like this one report a little less. Some of these stretches carry variants that still affect immunity, skin and "
        "hair today.", "ancestry:0")
elif meta.get("archaic_unmeasurable"):
    add("Ancestry", "info", "Neanderthal DNA: not measured", meta["archaic_unmeasurable"], "", "ancestry:0")

# ---------------------------------------------------------------- the two lines
if has("lin_story"):
    for line, hg_key, label, tgt in (("Y", "y_haplogroup", "Paternal line", "lineage:0"), ("MT", "mt_haplogroup", "Maternal line", "lineage:1")):
        hg = meta.get(hg_key, "")
        if not hg:
            continue
        rows = [t for (t,) in cur.execute("SELECT text FROM lin_story WHERE line=? ORDER BY ord", (line,))]
        if not rows:
            continue
        anc = [t for t in rows if t.startswith("The ancient people closest")]
        clade = [t for t in rows if not t.startswith("The ancient people closest") and not t.startswith("Your own branch")]
        add("Your two lines", "accent", "%s: %s" % (label, hg), first_sentences(clade[-1], 300) if clade else "",
            first_sentences(anc[0], 320) if anc else "Ancient genomes on your branches are shown on a map on the Lineage page.", tgt)
if sex in ("male", "female"):
    if sex == "male":
        add("Your two lines", "info", "50.8% from your mother, 49.2% from your father",
            "Your mother gave you her X chromosome (156 million letters), your father his much smaller Y (57 million); every other "
            "chromosome is one copy from each. All your mitochondria came from your mother.",
            "So a man's DNA is a little more his mother's. Your X carries DNA from both your mother's parents and none from your father's "
            "side.", "lineage:0")
    else:
        add("Your two lines", "info", "Half from each parent",
            "Each parent gave you one copy of every chromosome, an X each. All your mitochondria came from your mother.",
            "The X from your father is the one he got from his mother, so it comes entirely from your paternal grandmother.", "lineage:0")

# ---------------------------------------------------------------- the data itself
kind = meta.get("seq_kind", "")
if chip:
    add("About your data", "info", kind or "Genotyping chip",
        "A chip reads a fixed set of about %s positions, chosen because they vary between people." % format(int(num(meta.get("callable_bases"))), ","),
        "Good for well-known, common variants; rare variants need sequencing to find reliably, so those pages say so instead of "
        "guessing. %s" % ("%s rare chip calls were set aside as unconfirmed." % format(int(num(meta.get("chip_rare"))), ",")
                          if num(meta.get("chip_rare")) else ""), "raw")
elif kind:
    note = ("Positions your file does not list are taken to match the reference where short reads can be trusted."
            if meta.get("ingest_ref_assumed") == "1" else "")
    add("About your data", "info", kind, ("%s variant calls on GRCh38. %s" % (format(int(num(meta.get("n_variants"))), ","), note)).strip(),
        "Sequencing reads nearly every letter, so rare variants can be found, though single calls can still be wrong, especially in "
        "repetitive DNA.", "raw")

cur.executescript("DROP TABLE IF EXISTS brief; CREATE TABLE brief(section TEXT, ord INTEGER, tone TEXT, title TEXT, text TEXT, meaning TEXT, target TEXT);")
cur.executemany("INSERT INTO brief VALUES(?,?,?,?,?,?,?)", items)
con.commit()
log("brief: %d items in %d sections" % (len(items), len({i[0] for i in items})))
