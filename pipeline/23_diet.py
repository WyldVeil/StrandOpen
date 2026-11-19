#!/usr/bin/env python3
"""Step 23: diet and supplements, what genes can and cannot say.

Built from the trait results (step 8) and the medicines step (G6PD, step 7):
for each gene and nutrient pair a genome can speak to, the person's result,
how strong the evidence is, and what it means, including the popular claims
the evidence does not support. Genes rarely change what anyone should eat or
take; a blood test of the nutrient itself usually tells more.

Grades: Guideline (a clinical guideline acts on it), Strong, Moderate, Weak,
Not supported (marketed, but the evidence says otherwise).

Table diet(ord, topic, grade, result, genotype, says, means, source, target).
"""
from strandlib import *

con = db()
cur = con.cursor()


def has(t):
    return cur.execute("SELECT 1 FROM sqlite_master WHERE name=?", (t,)).fetchone() is not None


meta = dict(cur.execute("SELECT key, value FROM meta").fetchall())
tr = {}
if has("trait"):
    for tid, result, genotype, copies, status in cur.execute("SELECT id, result, genotype, copies, status FROM trait"):
        tr[tid] = (result or "", genotype or "", copies, status or "")
main = ""
if has("anc_admix"):
    r = cur.execute("SELECT pop, frac FROM anc_admix WHERE level='super' ORDER BY frac DESC LIMIT 1").fetchone()
    main = r[0] if r and r[1] >= 0.8 else ""
rows = []


def add(topic, grade, result, genotype, says, means, source, target="traits"):
    rows.append((len(rows), topic, grade, result, genotype, says, means, source, target))


def called(tid):
    return tid in tr and tr[tid][3] != "nocall"


# ---- iron
if called("hfe282"):
    c282 = tr["hfe282"][2] or 0
    h63 = (tr["hfe63"][2] or 0) if called("hfe63") else 0
    gt = "; ".join(x for x in (tr["hfe282"][1], tr["hfe63"][1] if "hfe63" in tr else "") if x)
    if c282 >= 2:
        add("Iron (HFE haemochromatosis)", "Guideline", "Two copies of C282Y", gt,
            "Two copies of the main haemochromatosis variant. Many people with two copies never become ill, but some slowly build up "
            "iron over the years, men far more often than women: in one large study about 28 in 100 men and 1 in 100 women with two "
            "copies developed iron-related disease.",
            "A blood test (ferritin and transferrin saturation) shows whether iron is building up, and regular blood removal treats it if "
            "it is. European liver guidelines advise people with haemochromatosis against iron and vitamin C supplements, and say a "
            "special diet is not needed. Worth raising with a doctor.",
            "Allen et al. 2008, New England Journal of Medicine; EASL Clinical Practice Guidelines on haemochromatosis, 2022")
    elif c282 == 1 and h63 >= 1:
        add("Iron (HFE haemochromatosis)", "Guideline", "C282Y and H63D, one copy each", gt,
            "One copy each of C282Y and H63D. This mild combination only rarely leads to iron overload, usually when something else, such "
            "as heavy drinking or liver disease, adds to it.",
            "Nothing to change on genes alone. If a blood test ever shows high iron, this combination is part of the explanation.",
            "EASL Clinical Practice Guidelines on haemochromatosis, 2022")
    elif c282 == 1 or h63 >= 1:
        add("Iron (HFE haemochromatosis)", "Strong", "Carrier", gt,
            "One copy of an HFE variant (%s). Carriers handle iron normally; it matters only for children, if a partner also carries one." %
            ("C282Y" if c282 else "H63D" if h63 == 1 else "H63D, two copies"),
            "No reason to avoid iron or change diet." if h63 < 2 else
            "Two copies of the mild H63D variant rarely cause iron overload; no reason to avoid iron or change diet.",
            "EASL Clinical Practice Guidelines on haemochromatosis, 2022")
    else:
        add("Iron (HFE haemochromatosis)", "Strong", "No haemochromatosis variants", gt,
            "Neither of the common HFE variants that cause iron overload.",
            "Nothing in these genes affects how much iron you should eat or take; iron needs depend on age, sex, diet and blood loss.",
            "EASL Clinical Practice Guidelines on haemochromatosis, 2022")

# ---- G6PD (from the medicines step)
if has("pgx_gene"):
    r = cur.execute("SELECT phenotype, diplotype FROM pgx_gene WHERE gene='G6PD' AND called=1").fetchone()
    if r and r[0] and "normal" not in r[0].lower() and "not" not in r[0].lower():
        add("Fava beans and some remedies (G6PD)", "Guideline", r[0], r[1] or "",
            "A G6PD variant lowers an enzyme that protects red blood cells, so some foods, medicines and chemicals can make them break "
            "down.",
            "Fava (broad) beans are the classic trigger; so are some medicines (the Medications page lists them, with CPIC's advice), and "
            "henna and naphthalene mothballs. Worth telling a doctor or pharmacist.",
            "CPIC guideline for G6PD", "pgx:0")

# ---- lactose
if called("lactase"):
    res, gt = tr["lactase"][0], tr["lactase"][1]
    persistent = (tr["lactase"][2] or 0) > 0
    non_eur = main not in ("EUR", "")
    add("Milk (lactase)", "Strong" if not non_eur else "Moderate", res, gt,
        ("You carry the European lactase persistence variant, which keeps the enzyme that digests milk sugar working into adulthood."
         if persistent else
         "You do not carry the European lactase persistence variant, so your lactase activity probably fell after childhood, as it "
         "does in most of the world's adults.")
        + (" Other persistence variants, common in parts of Africa and the Middle East, are not checked here, so this is less certain "
           "for you." if non_eur else ""),
        ("Nothing to change." if persistent else
         "That does not mean you must avoid milk: most adults with low lactase can handle about 12 g of lactose at a time (a glass of "
         "milk), more with food, and hard cheese and yoghurt contain little. Symptoms, not genes, are the guide."),
        "Enattah et al. 2002, Nature Genetics; EFSA Panel on Dietetic Products 2010 (lactose thresholds)")

# ---- gluten
if called("celiac"):
    res, gt = tr["celiac"][0], tr["celiac"][1]
    neg = res.startswith("None")
    add("Gluten (coeliac HLA types)", "Guideline", res, gt,
        ("You carry none of the HLA types that coeliac disease almost always needs, which makes it very unlikely." if neg else
         "You carry an HLA type found in nearly everyone with coeliac disease, and in about a third of people of European ancestry, "
         "most of whom never develop it."),
        ("Doctors use a negative HLA result like this to rule coeliac disease out." if neg else
         "No reason to avoid gluten because of this. If you have symptoms, the blood tests for coeliac disease only work while you "
         "are still eating gluten, so get tested before cutting it out."),
        "Monsuur et al. 2008, PLoS One; NICE guideline NG20 (coeliac disease)")

# ---- alcohol
if called("flush"):
    carrier = (tr["flush"][2] or 0) > 0
    if carrier:
        add("Alcohol (ALDH2)", "Strong", tr["flush"][0], tr["flush"][1],
            "You carry ALDH2*2, which slows the clearing of acetaldehyde, the toxic product of alcohol: it causes the face to flush.",
            "For people with this variant, regular drinking raises the risk of cancer of the oesophagus several times over, so it is a "
            "good reason to drink less.", "Brooks et al. 2009, PLoS Medicine")
    else:
        add("Alcohol (ALDH2)", "Strong", "No flush variant", tr["flush"][1],
            "You do not carry ALDH2*2, the alcohol flush variant common in East Asia.",
            "Alcohol's usual risks apply; nothing in this gene adds to them.", "Brooks et al. 2009, PLoS Medicine")

# ---- folate
if called("mthfr"):
    tt = "TT" in tr["mthfr"][0] or (tr["mthfr"][2] or 0) >= 2
    add("Folate (MTHFR)", "Not supported", tr["mthfr"][0], tr["mthfr"][1],
        ("Two copies of MTHFR 677T, carried by about 1 in 10 Europeans, slightly raise homocysteine when folate intake is low."
         if tt else "MTHFR C677T is one of the most marketed gene tests, but it says little about health."),
        "Ordinary folic acid works for everyone, whatever their MTHFR type, and the American College of Medical Genetics advises "
        "against MTHFR testing. There is no good evidence that 'methylated' folate supplements or special diets help people with these "
        "variants.", "ACMG practice guideline 2013, reaffirmed 2020; US CDC on MTHFR and folic acid")

# ---- caffeine
if called("caffeine"):
    add("Coffee (CYP1A2)", "Not supported", tr["caffeine"][0], tr["caffeine"][1],
        "This variant changes how quickly caffeine is cleared, which can affect how long its effect lasts.",
        "A 2006 study linked slow clearance plus heavy coffee drinking to heart attacks, but much larger later studies found coffee "
        "drinkers' health did not differ by these genes. Let how coffee makes you feel, and sleep, be the guide.",
        "Cornelis et al. 2006, JAMA; Loftfield et al. 2018, JAMA Internal Medicine; Zhou and Hypponen 2019, American Journal of Clinical Nutrition")

# ---- vitamin D
if called("vitd"):
    add("Vitamin D (GC)", "Weak", tr["vitd"][0], tr["vitd"][1],
        "Common GC variants shift blood vitamin D a little on average, but they cannot tell what your own level is.",
        "The advice is the same whatever the genes: in the UK, everyone is advised to consider 10 micrograms a day in autumn and "
        "winter. A blood test shows your actual level.", "Wang et al. 2010, Lancet (GC variants); NHS advice on vitamin D")

# ---- omega-3
if called("omega3"):
    add("Omega-3 (FADS1)", "Weak", tr["omega3"][0], tr["omega3"][1],
        "FADS1 variants change how readily the body makes the long omega-3 and omega-6 fats from plant sources; they have a clear "
        "effect on blood fat levels.",
        "Whether that changes health is unclear. Oily fish supplies the long omega-3 fats (EPA and DHA) directly, whatever your genes.",
        "Lemaitre et al. 2011, PLoS Genetics (CHARGE consortium)")

# ---- the general picture, first on the page
rows.insert(0, (-1, "What genes can and cannot say about diet", "", "", "",
                "Only a few genes change what anyone should eat or take, and most are listed below with how strong the evidence is. "
                "Blood tests of the nutrient itself (iron, vitamin D, B12, folate) tell far more than any of these variants.",
                "In the largest trial of gene-based diet advice (Food4Me, 1,600 people in seven countries), advice that used genes did no "
                "better than personalised advice without them.", "Celis-Morales et al. 2017, International Journal of Epidemiology", ""))
cur.executescript("DROP TABLE IF EXISTS diet; CREATE TABLE diet(ord INTEGER, topic TEXT, grade TEXT, result TEXT, genotype TEXT, says TEXT, "
                  "means TEXT, source TEXT, target TEXT);")
cur.executemany("INSERT INTO diet VALUES(?,?,?,?,?,?,?,?,?)", rows)
con.commit()
log("diet: %d topics" % (len(rows) - 1))
