"""Remarkable variants, part B: protective and "natural experiment" variants
in heart, metabolism, liver, immunity, infection and kidney.

Every allele is GRCh38, 1-based VCF style, left-normalised against the
reference with bcftools norm, and its REF checked against the FASTA. Common
alleles were cross-checked against the 1000 Genomes sites files.
"""

ENTRIES = [
    dict(
        key="b_pcsk9_lof",
        cat="Heart and lipids",
        gene="PCSK9",
        title="Naturally low LDL cholesterol (PCSK9 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="PCSK9 normally clears LDL receptors from the liver, so a broken copy leaves more "
               "receptors to pull LDL cholesterol out of the blood. In a 15-year US study of Black "
               "Americans, carriers of p.Tyr142Ter or p.Cys679Ter had about 28 percent lower LDL and "
               "about 88 percent lower risk of coronary heart disease, though that estimate rested on "
               "few cases. The finding helped lead to the PCSK9-blocking cholesterol drugs.",
        absent="No protein-truncating change was found in PCSK9.",
        rarity="rare: about 1 in 40 people of African ancestry, much rarer in others",
        source="Cohen 2005, Nat Genet; Cohen 2006, NEJM",
        pmid="16554528",
        conf="strong",
    ),
    dict(
        key="b_pcsk9_r46l",
        cat="Heart and lipids",
        gene="PCSK9",
        title="Slightly lower LDL cholesterol (PCSK9 p.Arg46Leu)",
        kind="alleles",
        alleles=[
            ("1", 55039974, "G", "T", "p.Arg46Leu (rs11591147)"),
        ],
        copies="one",
        effect="This change partly weakens PCSK9. Carriers have on average about 10 to 15 percent lower "
               "LDL cholesterol and a lower risk of coronary heart disease; the first study found about "
               "half the risk, and larger later studies found a smaller but real reduction.",
        absent="You do not carry the PCSK9 p.Arg46Leu variant.",
        rarity="uncommon: about 1 in 25 Europeans carry one copy",
        source="Cohen 2006, NEJM; Benn 2010, J Am Coll Cardiol",
        pmid="16554528",
        conf="strong",
    ),
    dict(
        key="b_angptl3_lof",
        cat="Heart and lipids",
        gene="ANGPTL3",
        title="Low blood fats (ANGPTL3 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="ANGPTL3 slows the enzymes that clear fat from the blood. People with one broken copy "
               "have on average about 17 to 27 percent lower triglycerides, about 10 percent lower LDL "
               "cholesterol and, in two large studies, about 35 to 40 percent lower odds of coronary "
               "artery disease. People with two broken copies have very low levels of all the main blood"
               " fats and appear healthy.",
        absent="No protein-truncating change was found in ANGPTL3.",
        rarity="rare: about 1 in 300 people",
        source="Dewey 2017, NEJM; Stitziel 2017, J Am Coll Cardiol",
        pmid="28538136",
        conf="strong",
    ),
    dict(
        key="b_apoc3_lof",
        cat="Heart and lipids",
        gene="APOC3",
        title="Low triglycerides (APOC3 loss of function)",
        kind="alleles",
        alleles=[
            ("11", 116830637, "C", "T", "p.Arg19Ter (rs76353203)"),
            ("11", 116830638, "G", "A", "IVS2+1G>A splice (rs138326449)"),
            ("11", 116830897, "G", "T", "IVS3+1G>T splice (rs140621530)"),
            ("11", 116830844, "G", "A", "p.Ala43Thr (rs147210663)"),
        ],
        copies="one",
        effect="APOC3 slows the breakdown of fat-carrying particles in the blood. Carriers of one of "
               "these four variants have about 40 percent lower triglycerides and about 40 percent lower"
               " risk of coronary heart disease. Drugs that block APOC3 are now in use for very high "
               "triglycerides.",
        absent="You do not carry any of the four protective APOC3 variants checked here.",
        rarity="rare: about 1 in 150 people carry one of them",
        source="TG and HDL Working Group (Crosby) 2014, NEJM; Jorgensen 2014, NEJM",
        pmid="24941081",
        conf="strong",
    ),
    dict(
        key="b_angptl4_e40k",
        cat="Heart and lipids",
        gene="ANGPTL4",
        title="Lower triglycerides (ANGPTL4 p.Glu40Lys)",
        kind="alleles",
        alleles=[
            ("19", 8364439, "G", "A", "p.Glu40Lys (rs116843064)"),
        ],
        copies="one",
        effect="This change disables ANGPTL4, a brake on the enzyme that clears fat from the blood. "
               "Carriers have lower triglycerides, higher HDL cholesterol and about 20 percent lower "
               "odds of coronary artery disease.",
        absent="You do not carry the ANGPTL4 p.Glu40Lys variant.",
        rarity="uncommon: about 1 in 25 Europeans carry one copy, rarer elsewhere and essentially absent"
               " in East Asians",
        source="Dewey 2016, NEJM; Stitziel 2016, NEJM",
        pmid="26933753",
        conf="strong",
    ),
    dict(
        key="b_asgr1_lof",
        cat="Heart and lipids",
        gene="ASGR1",
        title="Lower heart disease risk (ASGR1, Icelandic variants)",
        kind="alleles",
        alleles=[
            ("17", 7176936, "CCCCCAGCCCCAG", "C", "del12 in intron 4, causes frameshift"),
            ("17", 7174259, "C", "T", "p.Trp158Ter (c.473G>A)"),
            ("17", 7174258, "C", "T", "p.Trp158Ter (c.474G>A)"),
        ],
        copies="one",
        effect="These variants knock out one copy of ASGR1, a liver receptor. In studies of Icelanders "
               "and other Europeans, carriers had lower non-HDL cholesterol and about 34 percent lower "
               "risk of coronary artery disease, a larger drop than the cholesterol change alone would "
               "explain.",
        absent="You do not carry the known ASGR1 loss-of-function variants.",
        rarity="rare: about 1 in 120 Icelanders and roughly 1 in 80 other Europeans carry the main "
               "variant, rarer elsewhere",
        source="Nioi 2016, NEJM",
        pmid="27192541",
        conf="strong",
    ),
    dict(
        key="b_cetp_lof",
        cat="Heart and lipids",
        gene="CETP",
        title="High HDL cholesterol (CETP loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="CETP moves cholesterol from HDL to LDL particles. Carriers of a broken copy have much "
               "higher HDL, somewhat lower LDL and triglycerides, and about 30 percent lower odds of "
               "coronary heart disease. One broken copy is far more common in parts of Japan than "
               "elsewhere; two broken copies, rare even there, give very high HDL.",
        absent="No protein-truncating change was found in CETP.",
        rarity="rare: well under 1 in 100 people outside Japan",
        source="Nomura 2017, Circ Res",
        pmid="28506971",
        conf="moderate",
    ),
    dict(
        key="b_npc1l1_lof",
        cat="Heart and lipids",
        gene="NPC1L1",
        title="Less cholesterol absorbed from food (NPC1L1 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="NPC1L1 absorbs cholesterol in the gut and is the target of the drug ezetimibe. Carriers "
               "of a broken copy had LDL about 12 mg/dL (0.3 mmol/L) lower and about half the risk of "
               "coronary heart disease in the study that found them.",
        absent="No protein-truncating change was found in NPC1L1.",
        rarity="rare: about 1 in 650 people",
        source="Stitziel 2014, NEJM",
        pmid="25390462",
        conf="moderate",
    ),
    dict(
        key="b_lpl_s447x",
        cat="Heart and lipids",
        gene="LPL",
        title="Better fat clearance (LPL S447X)",
        kind="alleles",
        alleles=[
            ("8", 19962213, "C", "G", "p.Ser474Ter, S447X in mature LPL (rs328)"),
        ],
        copies="one",
        effect="Unusually, this stop change makes the fat-clearing enzyme LPL work slightly better. "
               "Carriers have somewhat lower triglycerides, slightly higher HDL and a modestly lower "
               "risk of coronary heart disease, roughly 10 to 20 percent.",
        absent="You do not carry the LPL S447X variant (p.Ser474Ter).",
        rarity="common: about 1 in 4 Europeans carry one copy",
        source="Rip 2006, Arterioscler Thromb Vasc Biol",
        pmid="16574898",
        conf="moderate",
    ),
    dict(
        key="b_b4galt1_n352s",
        cat="Heart and lipids",
        gene="B4GALT1",
        title="Lower LDL and fibrinogen (B4GALT1, Amish variant)",
        kind="alleles",
        alleles=[
            ("9", 33113783, "T", "C", "p.Asn352Ser (rs551564683)"),
        ],
        copies="one",
        effect="This change alters how the body adds sugar chains to proteins. Each copy lowers LDL "
               "cholesterol by about 14 mg/dL (0.36 mmol/L) and also lowers fibrinogen, a clotting "
               "protein. Whether it lowers heart disease is not yet proven.",
        absent="You do not carry the B4GALT1 p.Asn352Ser variant.",
        rarity="about 1 in 8 Old Order Amish carry it; extremely rare elsewhere",
        source="Montasser 2021, Science",
        pmid="34855475",
        conf="moderate",
    ),
    dict(
        key="b_slc30a8_lof",
        cat="Metabolism and weight",
        gene="SLC30A8",
        title="Protection from type 2 diabetes (SLC30A8 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="SLC30A8 loads zinc into insulin-making cells. People with one broken copy (for example "
               "p.Arg138Ter or p.Lys34SerfsTer50) had about 65 percent lower risk of type 2 diabetes, "
               "even though the gene was expected to do the opposite.",
        absent="No protein-truncating change was found in SLC30A8.",
        rarity="rare: under 1 in 500 in most groups, more common in western Finland and Iceland",
        source="Flannick 2014, Nat Genet",
        pmid="24584071",
        conf="strong",
    ),
    dict(
        key="b_gpr75_lof",
        cat="Metabolism and weight",
        gene="GPR75",
        title="Lower body weight (GPR75 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="Carriers of a broken GPR75 copy weighed on average about 5 kg less (BMI about 1.8 lower)"
               " and had about half the odds of obesity in a study of 640,000 people. Mice without the "
               "gene resist weight gain on a rich diet.",
        absent="No protein-truncating change was found in GPR75.",
        rarity="rare: about 1 in 3,000 people",
        source="Akbari 2021, Science",
        pmid="34210852",
        conf="moderate",
    ),
    dict(
        key="b_inhbe_lof",
        cat="Metabolism and weight",
        gene="INHBE",
        title="Healthier fat distribution (INHBE loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="INHBE makes a liver hormone that tells fat tissue to store fat. Carriers of a broken "
               "copy have a smaller waist-to-hip ratio for their weight, less belly fat, and in some "
               "analyses a lower risk of type 2 diabetes.",
        absent="No protein-truncating change was found in INHBE.",
        rarity="rare: well under 1 in 100 people",
        source="Akbari 2022, Nat Commun; Deaton 2022, Nat Commun",
        pmid="35999217",
        conf="moderate",
    ),
    dict(
        key="b_mc4r_gof",
        cat="Metabolism and weight",
        gene="MC4R",
        title="Slightly lower weight (MC4R gain of function)",
        kind="alleles",
        alleles=[
            ("18", 60372043, "C", "T", "p.Val103Ile (rs2229616)"),
            ("18", 60371599, "T", "G", "p.Ile251Leu (rs52820871)"),
        ],
        copies="one",
        effect="MC4R is the brain's main fullness switch. These two changes make it signal slightly more"
               " strongly. Carriers are on average a little lighter (BMI about 0.3 to 0.5 lower) with "
               "lower odds of obesity, and in one large study also lower risk of type 2 diabetes and "
               "heart disease.",
        absent="You do not carry either MC4R gain-of-function variant.",
        rarity="uncommon: about 1 in 30 people carry p.Val103Ile, about 1 in 60 Europeans p.Ile251Leu",
        source="Lotta 2019, Cell",
        pmid="31002796",
        conf="moderate",
    ),
    dict(
        key="b_mc4r_lof",
        cat="Metabolism and weight",
        gene="MC4R",
        title="Stronger appetite from early childhood (MC4R loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="A broken copy of MC4R weakens the brain's fullness signal. Carriers tend to feel "
               "hungrier, gain weight from childhood and are often several BMI units heavier, though "
               "lifestyle still matters. It is the most common single-gene cause of obesity, and a "
               "specific drug (setmelanotide) acts on this pathway.",
        absent="No protein-truncating change was found in MC4R.",
        rarity="rare: truncating MC4R changes are well under 1 in 1,000 people",
        source="Farooqi 2003, NEJM",
        pmid="12646665",
        conf="strong",
    ),
    dict(
        key="b_crebrf_r457q",
        cat="Metabolism and weight",
        gene="CREBRF",
        title="Heavier but less diabetes (CREBRF, Samoan variant)",
        kind="alleles",
        alleles=[
            ("5", 173108771, "G", "A", "p.Arg457Gln (rs373863828)"),
        ],
        copies="one",
        effect="In Samoans, each copy raises BMI by about 1.4 units, yet carriers have about 40 percent "
               "lower odds of type 2 diabetes. It is thought to help cells store energy, a possible "
               "advantage during past food shortages.",
        absent="You do not carry the CREBRF p.Arg457Gln variant.",
        rarity="common in Polynesians (about 1 in 4 gene copies in Samoans, 1 in 10 to 1 in 4 in other "
               "groups); a few percent or less in Micronesians and Melanesians; almost absent elsewhere",
        source="Minster 2016, Nat Genet",
        pmid="27455349",
        conf="strong",
    ),
    dict(
        key="b_pparg_p12a",
        cat="Metabolism and weight",
        gene="PPARG",
        title="Slightly lower diabetes risk (PPARG p.Pro12Ala)",
        kind="alleles",
        alleles=[
            ("3", 12351626, "C", "G", "p.Pro12Ala (rs1801282)"),
        ],
        copies="one",
        effect="PPARG controls fat cell growth and is the target of the glitazone diabetes drugs. The "
               "alanine form is linked to better insulin sensitivity and roughly 15 to 20 percent lower "
               "risk of type 2 diabetes.",
        absent="You do not carry the PPARG p.Pro12Ala variant.",
        rarity="common: about 1 in 4 Europeans carry one copy",
        source="Altshuler 2000, Nat Genet",
        pmid="10973253",
        conf="strong",
    ),
    dict(
        key="b_hsd17b13_ta",
        cat="Liver",
        gene="HSD17B13",
        title="Protection from liver disease (HSD17B13 splice variant)",
        kind="alleles",
        alleles=[
            ("4", 87310240, "T", "TA", "rs72613567 TA insertion"),
        ],
        copies="one",
        effect="This insertion disrupts splicing and shuts off HSD17B13, an enzyme in liver fat "
               "droplets. One copy lowered the risk of alcohol-related liver disease by about 40 percent"
               " and fatty liver disease by about 17 percent; two copies gave larger reductions. Several"
               " drugs now try to copy this effect.",
        absent="You do not carry the protective HSD17B13 variant.",
        rarity="common: about 4 in 10 Europeans carry at least one copy",
        source="Abul-Husn 2018, NEJM",
        pmid="29562163",
        conf="strong",
    ),
    dict(
        key="b_cideb_lof",
        cat="Liver",
        gene="CIDEB",
        title="Protection from liver disease (CIDEB loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="CIDEB helps fat droplets in the liver grow. In a large study, people with a rare "
               "truncating or missense change in one copy had about 33 percent lower odds of liver "
               "disease of any cause and about half the odds of cirrhosis. Truncating changes alone "
               "showed a similar fall in liver enzymes, and the protective link was strongest in people "
               "with obesity or type 2 diabetes.",
        absent="No protein-truncating change was found in CIDEB.",
        rarity="rare: about 1 in 700 people",
        source="Verweij 2022, NEJM",
        pmid="35939579",
        conf="moderate",
    ),
    dict(
        key="b_ccr5_delta32",
        cat="Immunity and infection",
        gene="CCR5",
        title="HIV resistance (CCR5 delta32)",
        kind="alleles",
        alleles=[
            ("3", 46373452, "TACAGTCAGTATCAATTCTGGAAGAATTTCCAG", "T", "delta32 (rs333)"),
        ],
        copies="two",
        effect="Most HIV strains enter white blood cells through the CCR5 receptor. With two copies of "
               "this 32-letter deletion there is no CCR5 on the cell surface and strong resistance to "
               "the common form of HIV; one copy slows the disease somewhat. Two copies may raise the "
               "risk from West Nile virus.",
        absent="You do not carry the CCR5 delta32 deletion.",
        rarity="about 1 in 5 Europeans carry one copy and 1 in 100 have two; rare outside Europe",
        source="Samson 1996, Nature; Liu 1996, Cell",
        pmid="8751444",
        conf="strong",
    ),
    dict(
        key="b_fut2_nonsecretor",
        cat="Immunity and infection",
        gene="FUT2",
        title="Norovirus resistance (FUT2 non-secretor)",
        kind="alleles",
        alleles=[
            ("19", 48703417, "G", "A", "p.Trp154Ter (rs601338)"),
        ],
        copies="two",
        effect="With two copies you are a non-secretor: your blood group sugars are not released into "
               "saliva and the gut lining. Non-secretors are largely resistant to the most common "
               "norovirus strains and some rotavirus strains, but have a slightly higher risk of a few "
               "other conditions such as Crohn's disease.",
        absent="You do not have two copies of FUT2 p.Trp154Ter, the usual non-secretor variant in "
               "Europeans and Africans, so you are probably a secretor. Other FUT2 variants, the main "
               "ones in East Asia, are not checked here.",
        rarity="common: about 1 in 5 Europeans and Africans have two copies",
        source="Lindesmith 2003, Nat Med",
        pmid="12692541",
        conf="strong",
    ),
    dict(
        key="b_ackr1_duffy_null",
        cat="Immunity and infection",
        gene="ACKR1",
        title="Vivax malaria resistance (Duffy-null)",
        kind="alleles",
        alleles=[
            ("1", 159204893, "T", "C", "c.-67T>C GATA box (rs2814778)"),
        ],
        copies="two",
        effect="With two copies, red blood cells lack the Duffy protein that the malaria parasite "
               "Plasmodium vivax uses to get in, giving strong resistance to vivax malaria. It also "
               "gives a naturally lower neutrophil count on blood tests, which is harmless but can be "
               "mistaken for a problem.",
        absent="You do not have the Duffy-null genotype.",
        rarity="nearly universal in West and Central Africa; rare in Europeans and East Asians",
        source="Tournamille 1995, Nat Genet; Miller 1976, NEJM",
        pmid="7663520",
        conf="strong",
    ),
    dict(
        key="b_tyk2_p1104a",
        cat="Immunity and infection",
        gene="TYK2",
        title="Fewer autoimmune diseases, more tuberculosis risk (TYK2 p.Pro1104Ala)",
        kind="alleles",
        alleles=[
            ("19", 10352442, "G", "C", "p.Pro1104Ala (rs34536443)"),
        ],
        copies="one",
        effect="This change partly dampens TYK2, an immune signalling enzyme. One copy lowers the risk "
               "of several autoimmune diseases such as rheumatoid arthritis, lupus, psoriasis and "
               "multiple sclerosis. Two copies markedly raise the risk of tuberculosis if exposed. "
               "Findings like this helped mark TYK2 as a drug target, and a TYK2-blocking medicine is "
               "now used to treat psoriasis.",
        absent="You do not carry the TYK2 p.Pro1104Ala variant.",
        rarity="uncommon: about 1 in 12 Europeans carry one copy, about 1 in 600 have two",
        source="Dendrou 2016, Sci Transl Med; Boisson-Dupuis 2018, Sci Immunol",
        pmid="27807284",
        conf="strong",
    ),
    dict(
        key="b_ifih1_t1d",
        cat="Immunity and infection",
        gene="IFIH1",
        title="Lower type 1 diabetes risk (IFIH1 rare variants)",
        kind="alleles",
        alleles=[
            ("2", 162268127, "T", "C", "p.Ile923Val (rs35667974)"),
            ("2", 162277580, "C", "A", "p.Glu627Ter (rs35744605)"),
            ("2", 162279995, "C", "G", "IVS8+1G>C splice (rs35337543)"),
            ("2", 162268086, "C", "T", "IVS14+1G>A splice (rs35732034)"),
        ],
        copies="one",
        effect="IFIH1 senses viral RNA and switches on antiviral defences. These four variants weaken "
               "it, and each lowers the risk of type 1 diabetes by roughly a quarter to a half, probably"
               " by reducing the immune overreaction that can follow viral infections.",
        absent="You do not carry any of the four protective IFIH1 variants.",
        rarity="uncommon: each is carried by about 1 in 25 to 1 in 125 Europeans",
        source="Nejentsev 2009, Science",
        pmid="19264985",
        conf="strong",
    ),
    dict(
        key="b_hbb_sickle_trait",
        cat="Immunity and infection",
        gene="HBB",
        title="Malaria protection (sickle cell trait, HBB)",
        kind="alleles",
        alleles=[
            ("11", 5227002, "T", "A", "p.Glu7Val, HbS (rs334)"),
        ],
        copies="one",
        effect="One copy (sickle cell trait) protects strongly against severe malaria, by about 90 "
               "percent. Trait carriers are usually healthy, with small extra risks under extreme "
               "exertion or dehydration. Two copies cause sickle cell disease, so a child of two "
               "carriers has a 1 in 4 chance of it.",
        absent="You do not carry the sickle cell (HbS) variant.",
        rarity="about 1 in 13 Black Americans carry one copy; up to 1 in 4 people in parts of Africa",
        source="Allison 1954, BMJ; Taylor 2012, Lancet Infect Dis",
        pmid="22445352",
        conf="strong",
    ),
    dict(
        key="b_il23r_r381q",
        cat="Immunity and infection",
        gene="IL23R",
        title="Protection from Crohn's disease and psoriasis (IL23R p.Arg381Gln)",
        kind="alleles",
        alleles=[
            ("1", 67240275, "G", "A", "p.Arg381Gln (rs11209026)"),
        ],
        copies="one",
        effect="This change weakens the receptor for IL-23, a driver of gut and skin inflammation. In "
               "people of European ancestry, carriers have less than half the usual odds of Crohn's "
               "disease and somewhat lower odds of ulcerative colitis, psoriasis and ankylosing "
               "spondylitis. Drugs that block IL-23 are now used to treat Crohn's disease, ulcerative "
               "colitis and psoriasis.",
        absent="You do not carry the IL23R p.Arg381Gln variant.",
        rarity="uncommon: about 1 in 8 Europeans carry one copy",
        source="Duerr 2006, Science",
        pmid="17068223",
        conf="strong",
    ),
    dict(
        key="b_card9_splice",
        cat="Immunity and infection",
        gene="CARD9",
        title="Protection from inflammatory bowel disease (CARD9 splice variant)",
        kind="alleles",
        alleles=[
            ("9", 136365140, "C", "G", "c.IVS11+1G>C splice (rs141992399)"),
        ],
        copies="one",
        effect="This splice change makes a shortened CARD9 protein that blunts one inflammation signal. "
               "Carriers had roughly one third the risk of Crohn's disease and ulcerative colitis.",
        absent="You do not carry the protective CARD9 splice variant.",
        rarity="rare: about 1 in 125 Europeans",
        source="Rivas 2011, Nat Genet",
        pmid="21983784",
        conf="strong",
    ),
    dict(
        key="b_prnp_g127v",
        cat="Immunity and infection",
        gene="PRNP",
        title="Prion disease resistance (PRNP p.Gly127Val, kuru)",
        kind="alleles",
        alleles=[
            ("20", 4699600, "G", "T", "p.Gly127Val"),
        ],
        copies="one",
        effect="Found only in the kuru region of Papua New Guinea, this change became common among the "
               "Fore during the kuru epidemic, a prion disease spread by funeral practices. No kuru "
               "patient tested carried it, and families with it had less kuru. Mice making both variant "
               "and normal protein, like a carrier, resisted kuru and classical CJD, though variant CJD "
               "infected some; mice with only the variant resisted every strain tested.",
        absent="You do not carry the PRNP p.Gly127Val variant.",
        rarity="found only in the kuru-affected region of Papua New Guinea",
        source="Mead 2009, NEJM; Asante 2015, Nature",
        pmid="19923577",
        conf="moderate",
    ),
    dict(
        key="b_apol1_g1g2",
        cat="Kidney",
        gene="APOL1",
        title="Sleeping sickness protection, kidney risk (APOL1 G1/G2)",
        kind="alleles",
        alleles=[
            ("22", 36265860, "A", "G", "G1: p.Ser342Gly (rs73885319)"),
            ("22", 36265995, "AATAATT", "A", "G2: p.Asn388_Tyr389del (rs71785313)"),
        ],
        copies="two",
        effect="One copy of G1 or G2 helps the blood kill the parasite behind the East African form of "
               "sleeping sickness in laboratory tests, and G2 was linked to protection in a Ugandan "
               "study; the West African form, which causes most cases, largely resists it. Two copies "
               "(any mix of G1 and G2) raise the risk of some kidney diseases several-fold, though most "
               "people with two copies never develop kidney disease.",
        absent="You do not carry the APOL1 G1 or G2 variants.",
        rarity="about 1 in 8 African Americans and about 1 in 4 people in Ghana and Nigeria have two "
               "copies; very rare without recent African ancestry",
        source="Genovese 2010, Science",
        pmid="20647424",
        conf="strong",
    ),
]
