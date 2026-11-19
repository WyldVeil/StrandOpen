"""Remarkable variants, part D: famous cases found missing by a completeness
review of parts A to C (ancient-human DNA, natural experiments, landmark
risk and protective variants). Each entry was written and then checked
again independently: REF against GRCh38.fa, position and alleles against
dbSNP, frequencies against the local 1000 Genomes sites, and every claim
and PMID against the papers.

Where the reference genome itself carries the allele with the effect
(OAS1, FOXO3, IFNL4), the entry describes the other allele, so the app's
"carries the ALT allele" test still reads the right way.

kind 'alleles'   the person carries any of the listed alleles
kind 'gene_lof'  any rare protein-truncating variant in the gene
copies 'one'     the effect shows with one copy
       'two'     the effect needs two copies; one copy means a carrier
"""

ENTRIES = [
    dict(
        key="d_epas1_tibetan_denisovan",
        cat="Ancient ancestry",
        gene="EPAS1",
        title="Tibetan high-altitude marker on Denisovan-derived DNA (EPAS1)",
        kind="alleles",
        alleles=[
            ("2", 46360880, "C", "G", "rs150877473 derived G (Tibetan EPAS1 tag)"),
        ],
        copies="one",
        effect="This is the EPAS1 change most enriched in Tibetans (about 87 percent of their copies, "
               "against about 1 percent in other East Asians) and lies on DNA inherited from Denisovans "
               "or a related archaic group. Above 4,000 m, Tibetan carriers have lower red cell counts "
               "and haemoglobin, which may spare them the overly thick blood of chronic mountain "
               "sickness. No effect at sea level is known.",
        absent="You do not carry the G version of rs150877473, the main marker of the Tibetan high-"
               "altitude form of EPAS1; this single marker cannot rule out other pieces of the "
               "Denisovan-derived stretch.",
        rarity="rare outside the Himalayas: about 1 in 50 East Asians carry one copy, essentially absent"
               " in Europeans and Africans",
        source="Huerta-Sanchez 2014, Nature; Yi 2010, Science; Beall 2010, PNAS; Hackinger 2016, Hum Genet",
        pmid="25043035",
        conf="strong",
    ),
    dict(
        key="d_tbc1d4_arg684ter",
        cat="Metabolism and weight",
        gene="TBC1D4",
        title="Greenlandic muscle insulin resistance (TBC1D4 p.Arg684Ter)",
        kind="alleles",
        alleles=[
            ("13", 75324385, "G", "A", "p.Arg684Ter (rs61736969)"),
        ],
        copies="two",
        effect="This stop change removes the muscle form of TBC1D4, a switch that helps insulin move "
               "glucose into muscle. In Greenlanders, people with two copies had blood sugar about 3.8 "
               "mmol/L higher two hours after a sugary drink and about 10 times the odds of type 2 "
               "diabetes, though fasting sugar was slightly lower. One copy gave only a small rise.",
        absent="You do not carry two copies of the TBC1D4 p.Arg684Ter variant.",
        rarity="common in Greenlanders (about 17 percent of gene copies); very rare elsewhere",
        source="Moltke 2014, Nature",
        pmid="25043022",
        conf="strong",
    ),
    dict(
        key="d_cpt1a_arctic_p479l",
        cat="Metabolism and weight",
        gene="CPT1A",
        title="Arctic fat-burning variant (CPT1A p.Pro479Leu)",
        kind="alleles",
        alleles=[
            ("11", 68780662, "G", "A", "p.Pro479Leu, c.1436C>T (rs80356779)"),
        ],
        copies="one",
        effect="This change weakens the liver enzyme that lets fat be burned for fuel when fasting. It "
               "rose to about 70 to 80 percent of gene copies in Inuit and Northeast Siberians in one of"
               " the strongest selective sweeps known in humans. In Greenlanders each copy went with "
               "about 2 cm less height; in Alaska Native infants, two copies went with about 2.5 times "
               "the odds of death, perhaps partly via infections.",
        absent="You do not carry the CPT1A Arctic variant p.Pro479Leu; you have the ancestral version "
               "found in almost everyone outside the Arctic.",
        rarity="very rare outside the Arctic (absent from 1000 Genomes); about 70 to 80 percent of gene "
               "copies in Inuit and Northeast Siberians",
        source="Clemente 2014, Am J Hum Genet; Skotte 2017, Circ Cardiovasc Genet; Gessner 2016, Genet Med",
        pmid="25449608",
        conf="moderate",
    ),
    dict(
        key="d_lztfl1_covid_neanderthal",
        cat="Immunity and infection",
        gene="LZTFL1",
        title="Higher risk of severe COVID-19 (Neanderthal chromosome 3 segment)",
        kind="alleles",
        alleles=[
            ("3", 45867532, "A", "G", "Neanderthal core haplotype tag (rs35044562)"),
        ],
        copies="one",
        effect="This stretch of chromosome 3, inherited from Neanderthals, is the strongest common "
               "genetic risk factor found for severe COVID-19. In early-pandemic studies each copy "
               "raised the odds of hospital admission or respiratory failure by about 1.6 to 1.8 times, "
               "more so in people aged 60 or under. The risk version seems to boost LZTFL1, a gene that "
               "shapes how cells lining the lungs respond to the virus.",
        absent="You do not carry the Neanderthal-derived chromosome 3 segment tagged here, though age, "
               "other health conditions and vaccination still affect COVID-19 risk.",
        rarity="common in some groups: about 1 in 6 Europeans and half of South Asians carry a copy, "
               "rare in East Asians and Africans",
        source="Zeberg 2020, Nature; Ellinghaus 2020, NEJM; Nakanishi 2021, J Clin Invest; Downes 2021, "
               "Nat Genet",
        pmid="32998156",
        conf="strong",
    ),
    dict(
        key="d_oas1_splice_aa",
        cat="Immunity and infection",
        gene="OAS1",
        title="No Neanderthal-like OAS1 protection against severe COVID-19",
        kind="alleles",
        alleles=[
            ("12", 112919388, "G", "A", "OAS1 splice acceptor, common A version (rs10774671)"),
        ],
        copies="two",
        effect="OAS1 helps cells destroy viral RNA. Its G version at this splice site makes a more "
               "active form of the enzyme; outside Africa G mostly came from Neanderthals. Two copies of"
               " the common A version mean no G. Each copy of the Neanderthal stretch carrying G was "
               "linked to about 22 percent lower risk of needing intensive care for COVID-19, and each G"
               " copy to 6 to 8 percent lower odds of hospital admission.",
        absent="You do not carry the common A version, so both your copies are the G version, linked to "
               "a modestly lower chance of severe COVID-19.",
        rarity="common: about 4 in 10 Europeans, over half of East Asians and about 1 in 8 people of "
               "African ancestry have two A copies",
        source="Zeberg 2021, PNAS; Zhou 2021, Nat Med; Huffman 2022, Nat Genet",
        pmid="33593941",
        conf="strong",
    ),
    dict(
        key="d_ifnl4_hcv_clearance",
        cat="Immunity and infection",
        gene="IFNL4",
        title="Less likely to clear hepatitis C unaided (IFNL4 rs12979860 T)",
        kind="alleles",
        alleles=[
            ("19", 39248147, "C", "T", "rs12979860 T (IFNL4 intron 1; the classic IL28B marker)"),
        ],
        copies="one",
        effect="This marker travels with a nearby IFNL4 change that lets liver cells make an extra "
               "interferon, which oddly goes with poorer control of hepatitis C. If infected, people "
               "with a T copy are less likely to clear the virus unaided: CC people had about 2.6 to 3.1"
               " times the odds of clearing it. Older interferon treatment cured CC people at least "
               "twice as often as TT; modern tablets cure most people regardless.",
        absent="You carry two copies of the C version, linked to a better chance of clearing hepatitis C"
               " unaided if ever infected.",
        rarity="common: about half of Europeans carry a T, most people of African ancestry, about 1 in 7"
               " East Asians",
        source="Thomas 2009, Nature; Ge 2009, Nature; Prokunina-Olsson 2013, Nat Genet",
        pmid="19759533",
        conf="strong",
    ),
    dict(
        key="d_hbg2_xmni",
        cat="Blood and oxygen",
        gene="HBG2",
        title="More fetal haemoglobin (HBG2 XmnI site)",
        kind="alleles",
        alleles=[
            ("11", 5254939, "G", "A", "HBG2 -158 C>T, XmnI site (rs7482144)"),
        ],
        copies="one",
        effect="The XmnI change just upstream of a fetal haemoglobin gene goes with that gene staying a "
               "little more active in adults, a small, harmless rise in healthy people. In sickle cell "
               "disease and beta thalassaemia, where fetal haemoglobin eases the illness, carriers tend "
               "to make more and have milder disease. Along with variants in BCL11A and near MYB, it "
               "helps explain 20 to 50 percent of the variation in levels.",
        absent="You do not carry the XmnI variant; fetal haemoglobin levels also depend on other genes "
               "such as BCL11A.",
        rarity="common: about 1 in 3 people carry at least one copy, about half of Europeans",
        source="Lettre 2008, PNAS; Thein 2009, Hum Mol Genet",
        pmid="18667698",
        conf="strong",
    ),
    dict(
        key="d_trem2_r47h",
        cat="Brain and memory",
        gene="TREM2",
        title="Raised Alzheimer's risk (TREM2 p.Arg47His)",
        kind="alleles",
        alleles=[
            ("6", 41161514, "C", "T", "p.Arg47His (rs75932628)"),
        ],
        copies="one",
        effect="TREM2 helps the brain's immune cells respond to and clear damage, and this change "
               "weakens it. In Iceland carriers had about 3 times the odds of late-onset Alzheimer's "
               "disease, similar to one copy of APOE e4, and a second study found about 5 times the "
               "odds. Many carriers never develop the disease, and absolute risk depends heavily on age,"
               " APOE type and other factors.",
        absent="You do not carry TREM2 p.Arg47His; your Alzheimer's risk is shaped by other genes, age "
               "and lifestyle.",
        rarity="rare: about 1 in 100 Europeans carry one copy; very rare in East Asians and Africans",
        source="Jonsson 2013, NEJM; Guerreiro 2013, NEJM",
        pmid="23150908",
        conf="strong",
    ),
    dict(
        key="d_apoa1_milano",
        cat="Heart and lipids",
        gene="APOA1",
        title="Very low HDL without the expected artery damage (apoA-I Milano)",
        kind="alleles",
        alleles=[
            ("11", 116836023, "G", "A", "p.Arg197Cys, mature Arg173Cys, apoA-I Milano (rs28931573)"),
        ],
        copies="one",
        effect="This change lets apoA-I, the main protein of HDL, form unusual paired molecules. "
               "Carriers have very low HDL cholesterol and raised triglycerides, yet in a 2001 study of "
               "21 carriers their neck artery walls were no thicker than relatives' (0.63 vs 0.64 mm), "
               "unlike two groups with low HDL from other causes (about 0.87 mm). Lab-made infusions "
               "looked promising in 2003 but did not help in a 2018 trial.",
        absent="You do not carry the apoA-I Milano variant, which is expected: it is known almost only "
               "in one Italian family.",
        rarity="extremely rare: known almost only in descendants of one family from Limone sul Garda, "
               "Italy",
        source="Franceschini 1980, J Clin Invest; Sirtori 2001, Circulation; Nissen 2003, JAMA; Nicholls"
               " 2018, JAMA Cardiol",
        pmid="7430351",
        conf="moderate",
    ),
    dict(
        key="d_foxo3_rs2802292",
        cat="Longevity",
        gene="FOXO3",
        title="Without the longevity-linked FOXO3 version (rs2802292 T/T)",
        kind="alleles",
        alleles=[
            ("6", 108587315, "G", "T", "rs2802292 T (non-longevity allele)"),
        ],
        copies="two",
        effect="People with two copies of the T version at this FOXO3 marker lack the G allele linked to"
               " reaching very old age. In Japanese-American men in Hawaii, those with one G had about "
               "1.9 times and those with two G about 2.75 times the odds of reaching 95 rather than "
               "dying before 81, compared with T/T men. Later studies found smaller effects, clearest in"
               " men; many T/T people live long lives.",
        absent="You carry at least one copy of the G allele, the FOXO3 version linked to longevity in "
               "several studies, though its effect on any one person's lifespan is small.",
        rarity="common: about 1 in 3 Europeans, half of East Asians and 1 in 30 people of African "
               "ancestry are T/T",
        source="Willcox 2008, PNAS; Flachsbart 2009, PNAS; Bao 2014, Asian J Androl",
        pmid="18765803",
        conf="moderate",
    ),
    dict(
        key="d_tcf7l2_rs7903146",
        cat="Metabolism and weight",
        gene="TCF7L2",
        title="Higher risk of type 2 diabetes (TCF7L2 rs7903146 T)",
        kind="alleles",
        alleles=[
            ("10", 112998590, "C", "T", "intronic T allele (rs7903146)"),
        ],
        copies="one",
        effect="Each copy of the T allele raises the odds of type 2 diabetes by about 45 percent; people"
               " with two copies have about two to two and a half times the odds of those with none. It "
               "seems to act mainly by cutting insulin release from the pancreas, not by causing insulin"
               " resistance. It is the strongest common genetic risk factor for the disease in "
               "Europeans, but weight, activity and age matter more.",
        absent="You do not carry the higher-risk T allele at rs7903146, which removes one common genetic"
               " push towards type 2 diabetes but does not by itself make you low risk.",
        rarity="very common: about half of Europeans and South Asians and about 45 percent of Africans "
               "carry at least one copy; about 1 in 20 East Asians",
        source="Grant 2006, Nat Genet; Cauchi 2007, J Mol Med; Florez 2006, NEJM",
        pmid="16415884",
        conf="strong",
    ),
    dict(
        key="d_erap2_black_death",
        cat="Famous claims",
        gene="ERAP2",
        title="The ERAP2 version thought to have fared worse in the Black Death (rs2549794 T)",
        kind="alleles",
        alleles=[
            ("5", 96908845, "C", "T", "rs2549794 T (lower ERAP2 expression)"),
        ],
        copies="one",
        effect="This T allele lowers ERAP2 levels. A 2022 ancient DNA study found it fell after the 1348"
               " Black Death in London and Denmark, estimating that people with two C copies were about "
               "40 percent likelier to survive than those with two T copies; a 2025 reanalysis found no "
               "clear change. Each T copy today brings slightly higher odds of severe pneumonia and "
               "about 14 percent lower odds of Crohn's disease.",
        absent="You do not carry the T allele here, so both your copies are probably the C version, "
               "which the 2022 study proposed was favoured during the Black Death.",
        rarity="very common: about 2 in 3 copies worldwide are T (about 6 in 10 in Europeans), so most "
               "people carry at least one",
        source="Klunk 2022, Nature; Barton 2025, Nature; Hamilton 2023, Am J Hum Genet",
        pmid="36261521",
        conf="disputed",
    ),
    dict(
        key="d_or7d4_wm",
        cat="Pain and senses",
        gene="OR7D4",
        title="A musky steroid smells weaker and milder (OR7D4 WM form)",
        kind="alleles",
        alleles=[
            ("19", 9214576, "G", "A", "p.Arg88Trp (rs61729907)"),
            ("19", 9214440, "G", "A", "p.Thr133Met (rs5020278)"),
        ],
        copies="one",
        effect="These two changes, almost always inherited together, largely disable OR7D4, a smell "
               "receptor tuned to androstenone, a musky steroid found in sweat and in meat from male "
               "pigs. In several studies, carriers were on average less sensitive to androstenone and "
               "found it less unpleasant than people with two working copies. It is only one influence: "
               "androstenone perception is about 30 percent heritable.",
        absent="You do not carry the common weakened (WM) form of OR7D4; people without it tend on "
               "average to find androstenone stronger and more unpleasant, though perception still "
               "varies a lot.",
        rarity="common: about 1 in 3 people of European or Asian ancestry carry at least one copy, about"
               " 1 in 9 of African ancestry",
        source="Keller 2007, Nature; Knaapila 2012, Chem Senses; Lunde 2012, PLoS One",
        pmid="17873857",
        conf="strong",
    ),
    dict(
        key="d_abcg2_q141k",
        cat="Metabolism and weight",
        gene="ABCG2",
        title="Higher uric acid and gout risk (ABCG2 p.Gln141Lys)",
        kind="alleles",
        alleles=[
            ("4", 88131171, "G", "T", "p.Gln141Lys (rs2231142)"),
        ],
        copies="one",
        effect="This change weakens ABCG2, a pump that moves uric acid out of the body through the "
               "kidneys and gut; in lab tests the altered pump moved about half as much. Carriers have "
               "modestly higher blood uric acid, and each copy raises the odds of gout by about 1.7 "
               "times. At least 1 in 10 gout cases in white Americans has been put down to this variant.",
        absent="You do not carry the ABCG2 p.Gln141Lys variant, though other genes, diet and kidney "
               "function also shape uric acid levels and gout risk.",
        rarity="common: about 1 in 5 Europeans and about half of East Asians carry at least one copy, "
               "rare in Africans",
        source="Woodward 2009, PNAS; Ichida 2012, Nat Commun",
        pmid="19506252",
        conf="strong",
    ),
    dict(
        key="d_9p21_cad_rs1333049",
        cat="Heart and lipids",
        gene="CDKN2B-AS1 (9p21)",
        title="Modestly higher coronary heart disease risk (9p21 locus, rs1333049)",
        kind="alleles",
        alleles=[
            ("9", 22125504, "G", "C", "9p21 risk allele C (rs1333049)"),
        ],
        copies="one",
        effect="Reported by several groups in 2007, this was the first coronary heart disease region "
               "found by genome-wide scans and remains among the strongest common signals. In one of "
               "those studies each copy of C raised the odds by about 36 percent; larger later studies "
               "put it at 20 to 30 percent. It acts independently of cholesterol and blood pressure, "
               "which stay the main risk factors you can change.",
        absent="You do not carry the 9p21 risk allele C at rs1333049, so you lack this one common, "
               "modest risk factor; it says nothing about the many others.",
        rarity="very common: about 3 in 4 Europeans and Asians carry at least one copy, about 4 in 10 "
               "Africans",
        source="Samani 2007, NEJM; Helgadottir 2007, Science; McPherson 2007, Science; Schunkert 2008, "
               "Circulation",
        pmid="17634449",
        conf="strong",
    ),
]
