"""Remarkable variants, part A: rare or unusual alleles with striking effects
outside heart and metabolism (sleep, pain, bone, blood, brain, appearance).

Every allele is (chrom, pos, ref, alt, label) on GRCh38, VCF style (1-based),
plus strand, left-normalised against the reference. Each REF was checked
against GRCh38.fa with samtools faidx, each position against dbSNP (NCBI
SPDI) and each protein change against the Ensembl 115 MANE Select CDS.

kind 'alleles'   the person carries any of the listed alleles
kind 'gene_lof'  any rare protein-truncating variant in the gene
copies 'one'     the effect shows with one copy
       'two'     the effect needs two copies; one copy means a carrier
"""

ENTRIES = [
    dict(
        key="a_dec2_p384r",
        cat="Sleep",
        gene="BHLHE41",
        title="Natural short sleeper (DEC2 P384R)",
        kind="alleles",
        alleles=[
            ("12", 26122364, "G", "C", "p.Pro384Arg (rs121912617)"),
        ],
        copies="one",
        effect="The first gene found for natural short sleep. A mother and daughter carrying it slept "
               "about 6 hours a night, against about 8 for relatives without it, and felt fully rested. "
               "Mice given the same change sleep less, and after being kept awake they make up less of "
               "the lost sleep than normal mice.",
        absent="You do not carry this variant, so your sleep need is not shortened by it.",
        rarity="found in a single family; essentially absent from population databases",
        source="He 2009, Science",
        pmid="19679812",
        conf="moderate",
    ),
    dict(
        key="a_dec2_y362h",
        cat="Sleep",
        gene="BHLHE41",
        title="Short sleep and resilience to sleep loss (DEC2 Y362H)",
        kind="alleles",
        alleles=[
            ("12", 26122431, "A", "G", "p.Tyr362His (rs1591838266)"),
        ],
        copies="one",
        effect="Found in one of a pair of non-identical twins. Over a week at home the carrier twin "
               "slept about an hour less a night on average than his twin, though in the sleep lab the "
               "gap was smaller and varied by night. During 38 hours without sleep he had fewer "
               "attention lapses, and afterwards slept about 90 minutes less. In cell tests the change "
               "weakens DEC2's grip on its target genes, like the P384R variant.",
        absent="You do not carry this variant.",
        rarity="found in one person in a twin study",
        source="Pellegrino 2014, Sleep",
        pmid="25083013",
        conf="limited",
    ),
    dict(
        key="a_adrb1_a187v",
        cat="Sleep",
        gene="ADRB1",
        title="Natural short sleeper (ADRB1 A187V)",
        kind="alleles",
        alleles=[
            ("10", 114044692, "C", "T", "p.Ala187Val (rs776439595)"),
        ],
        copies="one",
        effect="Carriers in one family slept about 6 hours a night or less with no daytime tiredness. "
               "Mice with the same change slept about an hour less a day; the variant makes the beta-1 "
               "adrenergic receptor less stable, which keeps wake-promoting brain cells more active.",
        absent="You do not carry this variant.",
        rarity="very rare in population databases",
        source="Shi 2019, Neuron",
        pmid="31473062",
        conf="moderate",
    ),
    dict(
        key="a_npsr1_y206h",
        cat="Sleep",
        gene="NPSR1",
        title="Natural short sleeper (NPSR1 Y206H)",
        kind="alleles",
        alleles=[
            ("7", 34827538, "T", "C", "p.Tyr206His (rs1406844918)"),
        ],
        copies="one",
        effect="Found in a family of natural short sleepers who need far less sleep than average. The "
               "change makes the neuropeptide S receptor more responsive; mice carrying it slept less "
               "and, unusually, kept their memory performance after sleep deprivation.",
        absent="You do not carry this variant.",
        rarity="found in a single family",
        source="Xing 2019, Science Translational Medicine",
        pmid="31619542",
        conf="limited",
    ),
    dict(
        key="a_grm1_nss",
        cat="Sleep",
        gene="GRM1",
        title="Natural short sleeper (GRM1)",
        kind="alleles",
        alleles=[
            ("6", 146352435, "T", "G", "p.Ser458Ala (rs151255685)"),
            ("6", 146426563, "A", "T", "p.Arg889Trp, mGluR1b isoform (rs768023437)"),
        ],
        copies="one",
        effect="Two different changes in the glutamate receptor mGluR1 were found in two unrelated "
               "families of people who naturally sleep 6 hours or less. Both reduce receptor signalling,"
               " and mice carrying either one slept about 25 minutes less a day.",
        absent="You do not carry either of these GRM1 variants.",
        rarity="rare: roughly 1 in 2,000 Europeans carry p.Ser458Ala and 1 in 30,000 carry p.Arg889Trp; "
               "rarer elsewhere",
        source="Shi 2021, Current Biology",
        pmid="33065013",
        conf="limited",
    ),
    dict(
        key="a_sik3_n783y",
        cat="Sleep",
        gene="SIK3",
        title="Natural short sleeper (SIK3 N783Y)",
        kind="alleles",
        alleles=[
            ("11", 116861352, "T", "A", "p.Asn783Tyr (rs767407831)"),
        ],
        copies="one",
        effect="Found in a person who naturally sleeps far less than average. The change lowers the "
               "activity of SIK3, a kinase that tracks sleep pressure; mice with two copies slept about "
               "30 minutes less a day and slept more deeply.",
        absent="You do not carry this variant.",
        rarity="very rare; about 1 in 8,000 people",
        source="Chen 2025, PNAS",
        pmid="40324078",
        conf="limited",
    ),
    dict(
        key="a_cry1_dsps",
        cat="Sleep",
        gene="CRY1",
        title="Night owl body clock (CRY1 delayed sleep phase)",
        kind="alleles",
        alleles=[
            ("12", 106992962, "T", "G", "c.1657+3A>C, exon 11 skipping (rs184039278)"),
        ],
        copies="one",
        effect="This splice change makes a CRY1 protein missing a small piece that holds the clock's "
               "brake on for longer, lengthening the body's daily rhythm. In the families studied, "
               "carriers slept about 40 minutes later on average, some with broken sleep. In a large UK "
               "population study the shift was only about 5 minutes, with a slightly higher chance of "
               "being an evening person. One copy is enough.",
        absent="You do not carry this variant, so it is not pushing your body clock later.",
        rarity="about 1 in 115 people of European ancestry; about 1 in 15 people of Ashkenazi Jewish "
               "ancestry",
        source="Patke 2017, Cell",
        pmid="28388406",
        conf="moderate",
    ),
    dict(
        key="a_per2_s662g",
        cat="Sleep",
        gene="PER2",
        title="Extreme early bird (PER2 S662G)",
        kind="alleles",
        alleles=[
            ("2", 238257003, "T", "C", "p.Ser662Gly (rs121908635)"),
        ],
        copies="one",
        effect="The first human body-clock mutation found. Carriers have familial advanced sleep phase: "
               "they get sleepy around 7 to 8 pm and wake around 3 to 4 am, about 4 hours earlier than "
               "usual. The change removes a phosphorylation site and shortens the clock's cycle.",
        absent="You do not carry this variant.",
        rarity="found in a single family",
        source="Toh 2001, Science",
        pmid="11232563",
        conf="moderate",
    ),
    dict(
        key="a_csnk1d_fasp",
        cat="Sleep",
        gene="CSNK1D",
        title="Extreme early bird (CSNK1D)",
        kind="alleles",
        alleles=[
            ("17", 82265743, "T", "C", "p.Thr44Ala (rs104894561)"),
            ("17", 82265736, "T", "C", "p.His46Arg (rs397514693)"),
        ],
        copies="one",
        effect="Changes in casein kinase 1 delta, the enzyme that tags PER2, cause familial advanced "
               "sleep phase: a strong drive to sleep and wake hours earlier than most people. Migraine "
               "was also common in the carrier families, and mice with Thr44Ala are more sensitive to "
               "migraine triggers.",
        absent="You do not carry these CSNK1D variants.",
        rarity="each found in a single family",
        source="Xu 2005, Nature; Brennan 2013, Science Translational Medicine",
        pmid="15800623",
        conf="moderate",
    ),
    dict(
        key="a_cry2_a260t",
        cat="Sleep",
        gene="CRY2",
        title="Extreme early bird (CRY2 A260T)",
        kind="alleles",
        alleles=[
            ("11", 45867648, "G", "A", "p.Ala260Thr (rs201220841)"),
        ],
        copies="one",
        effect="Found in a family with familial advanced sleep phase, sleeping and waking several hours "
               "earlier than usual. The change makes CRY2 easier to break down, which speeds up the "
               "clock; mice carrying it showed the same early shift.",
        absent="You do not carry this variant.",
        rarity="found in a single family; very rare in databases",
        source="Hirano 2016, eLife",
        pmid="27529127",
        conf="limited",
    ),
    dict(
        key="a_per3_fasp",
        cat="Sleep",
        gene="PER3",
        title="Early bird and seasonal mood (PER3 P415A/H417R)",
        kind="alleles",
        alleles=[
            ("1", 7809893, "C", "G", "p.Pro415Ala (rs150812083)"),
            ("1", 7809900, "A", "G", "p.His417Arg (rs139315125)"),
        ],
        copies="one",
        effect="Two changes that travel together on the same copy of PER3. In the one family where they "
               "were found, all three carriers had early sleep timing and higher scores for low mood and"
               " seasonal mood swings, mostly worse in winter. Mice given the pair showed shifted clock "
               "timing and some low-mood-like behaviour, mainly under short days. The pair is fairly "
               "common, and its effect outside that family is untested.",
        absent="You do not carry these PER3 variants.",
        rarity="about 1 in 50 people of European ancestry; rarer elsewhere",
        source="Zhang 2016, PNAS",
        pmid="26903630",
        conf="limited",
    ),
    dict(
        key="a_timeless_r1081x",
        cat="Sleep",
        gene="TIMELESS",
        title="Extreme early bird (TIMELESS R1081X)",
        kind="alleles",
        alleles=[
            ("12", 56418347, "G", "A", "p.Arg1081Ter (rs1465092391)"),
        ],
        copies="one",
        effect="A truncating change found in a small family with advanced sleep phase: the two carriers "
               "were strong morning types who woke about 2 hours earlier than a relative without it, "
               "though they fell asleep at about the same time. Mice with one copy of the matching "
               "change shifted sleep and waking about 30 minutes earlier and responded differently to "
               "light at night, but their body clock ran at a normal speed.",
        absent="You do not carry this variant.",
        rarity="described in a single family; very rare in public databases",
        source="Kurien 2019, PNAS",
        pmid="31138685",
        conf="limited",
    ),
    dict(
        key="a_scn9a_lof",
        cat="Pain and senses",
        gene="SCN9A",
        title="Congenital insensitivity to pain (SCN9A)",
        kind="gene_lof",
        alleles=[],
        copies="two",
        effect="SCN9A makes Nav1.7, the sodium channel that pain nerves need to fire. People with two "
               "broken copies have never felt physical pain and usually cannot smell. One broken copy, "
               "as in a carrier, has no known effect on pain.",
        absent="No rare truncating variant in SCN9A was found.",
        rarity="carriers are uncommon; complete pain insensitivity is about 1 in a million",
        source="Cox 2006, Nature",
        pmid="17167479",
        conf="strong",
    ),
    dict(
        key="a_ngf_r221w",
        cat="Pain and senses",
        gene="NGF",
        title="Loss of deep pain (NGF, Norrbotten)",
        kind="alleles",
        alleles=[
            ("1", 115286135, "G", "A", "p.Arg221Trp, R100W in mature NGF (rs11466112)"),
        ],
        copies="two",
        effect="With two copies, people lose deep pain from bones and joints and often have painless "
               "fractures and joint damage, yet keep normal intelligence and sweating. Carriers of one "
               "copy are mostly unaffected, though some develop mild nerve changes later in life.",
        absent="You do not carry this variant.",
        rarity="very rare; affected families are from Norrbotten in northern Sweden, and most carriers "
               "seen elsewhere are Finnish",
        source="Einarsdottir 2004, Human Molecular Genetics",
        pmid="14976160",
        conf="moderate",
    ),
    dict(
        key="a_zfhx2_r1913k",
        cat="Pain and senses",
        gene="ZFHX2",
        title="Reduced pain sensitivity (ZFHX2, Marsili)",
        kind="alleles",
        alleles=[
            ("14", 23524204, "C", "T", "p.Arg1913Lys (rs1555344723)"),
        ],
        copies="one",
        effect="An Italian family with this change feels heat and injury far less than most people; "
               "several members broke bones or burned themselves without much pain. Mice with the same "
               "change also have altered pain responses.",
        absent="You do not carry this variant.",
        rarity="very rare; described in one Italian family, and carried by only 3 of about 690,000 "
               "people in a large population database",
        source="Habib 2018, Brain",
        pmid="29253101",
        conf="limited",
    ),
    dict(
        key="a_lrp5_hbm",
        cat="Muscle and bone",
        gene="LRP5",
        title="Very dense bones (LRP5 high bone mass)",
        kind="alleles",
        alleles=[
            ("11", 68357673, "G", "T", "p.Gly171Val (rs121908668)"),
            ("11", 68357672, "G", "C", "p.Gly171Arg (rs121908669)"),
            ("11", 68357801, "G", "A", "p.Ala214Thr (rs121908671)"),
            ("11", 68357802, "C", "T", "p.Ala214Val (rs121908672)"),
            ("11", 68363784, "G", "A", "p.Ala242Thr (rs121908670)"),
            ("11", 68363818, "C", "T", "p.Thr253Ile (rs121908673)"),
        ],
        copies="one",
        effect="These changes stop natural brakes such as sclerostin and DKK1 from damping LRP5 "
               "signalling, so bone keeps being built. Carriers have very dense, strong bones and rarely"
               " break them. In two of the families, X-rays after one member's car crash showed no "
               "fractures but unusually dense bone. Some have a wide jaw or a bony lump on the palate.",
        absent="You do not carry any of these LRP5 high-bone-mass variants.",
        rarity="found in a handful of families worldwide",
        source="Boyden 2002, NEJM; Van Wesenbeeck 2003, Am J Hum Genet",
        pmid="12015390",
        conf="strong",
    ),
    dict(
        key="a_sost_lof",
        cat="Muscle and bone",
        gene="SOST",
        title="Denser bones (SOST, sclerostin)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="Sclerostin is the body's brake on bone building. Carriers of one broken copy tend to "
               "have higher bone density than others of their age, from normal to clearly raised, with "
               "no known harm. Two broken copies cause sclerosteosis, an overgrowth of bone. The "
               "osteoporosis drug romosozumab copies this effect by blocking sclerostin.",
        absent="No rare truncating variant in SOST was found.",
        rarity="very rare; most known carriers are in Afrikaner families in South Africa",
        source="Gardner 2005, J Clin Endocrinol Metab; Brunkow 2001, Am J Hum Genet",
        pmid="16189254",
        conf="moderate",
    ),
    dict(
        key="a_mstn_lof",
        cat="Muscle and bone",
        gene="MSTN",
        title="Extra muscle (myostatin, MSTN)",
        kind="gene_lof",
        alleles=[],
        copies="two",
        effect="Myostatin limits muscle growth. A boy with two broken copies was born with unusually "
               "large muscles and very little body fat, like double-muscled cattle and whippets. His "
               "mother, with one copy, was a strong athlete, but the effect of a single copy has not "
               "been measured in people.",
        absent="No rare truncating variant in MSTN was found.",
        rarity="extremely rare; two broken copies are known in very few people",
        source="Schuelke 2004, NEJM",
        pmid="15215484",
        conf="limited",
    ),
    dict(
        key="a_stc2_r44l",
        cat="Muscle and bone",
        gene="STC2",
        title="Extra height (STC2 R44L)",
        kind="alleles",
        alleles=[
            ("5", 173328063, "C", "A", "p.Arg44Leu (rs148833559)"),
        ],
        copies="one",
        effect="One of the largest single-variant effects on adult height in the general population. "
               "Carriers are on average about 2 cm taller; the change weakens STC2, which normally holds"
               " back the growth-promoting enzyme PAPP-A.",
        absent="You do not carry this variant.",
        rarity="about 1 in 500 people of European ancestry",
        source="Marouli 2017, Nature",
        pmid="28146470",
        conf="strong",
    ),
    dict(
        key="a_epor_mantyranta",
        cat="Blood and oxygen",
        gene="EPOR",
        title="Extra red blood cells (EPOR, Mantyranta)",
        kind="alleles",
        alleles=[
            ("19", 11378195, "C", "T", "p.Trp439Ter (rs121917830)"),
        ],
        copies="one",
        effect="This change cuts off the tail of the erythropoietin receptor that normally switches it "
               "off, so carriers make more red blood cells: haemoglobin is often 20 to 50 percent above "
               "normal, with no harm reported. It runs in the family of Eero Mantyranta, the Finnish "
               "skier who won three Olympic gold medals.",
        absent="You do not carry this variant.",
        rarity="found in one Finnish family",
        source="de la Chapelle 1993, PNAS",
        pmid="8506290",
        conf="strong",
    ),
    dict(
        key="a_piezo1_e756del",
        cat="Blood and oxygen",
        gene="PIEZO1",
        title="Slower-closing red cell channel (PIEZO1 E756del)",
        kind="alleles",
        alleles=[
            ("16", 88733964, "GTCC", "G", "p.Glu756del (rs59446030)"),
        ],
        copies="one",
        effect="This deletion makes the PIEZO1 channel in red cells slower to close when the cell is "
               "squeezed. One lab study found carriers' red cells were slightly dehydrated and less "
               "easily infected by the malaria parasite, which may explain why the variant is common in "
               "Africa. A later study found normal water content and normal parasite growth in carriers'"
               " cells, and studies of severe malaria in patients disagree.",
        absent="You do not carry this variant.",
        rarity="common in people of African ancestry (about 1 in 3 carry it); rare elsewhere",
        source="Ma 2018, Cell; Nguetse 2020, PNAS",
        pmid="29576450",
        conf="moderate",
    ),
    dict(
        key="a_egln1_tibetan",
        cat="Blood and oxygen",
        gene="EGLN1",
        title="High-altitude adaptation (EGLN1, Tibetan)",
        kind="alleles",
        alleles=[
            ("1", 231421877, "G", "C", "p.Asp4Glu (rs186996510)"),
        ],
        copies="one",
        effect="The Tibetan form of the oxygen sensor PHD2 carries this change and a second one (C127S, "
               "common worldwide on its own). In cell studies the pair dampens the usual surge in red "
               "blood cell production in thin air, and Tibetans, most of whom carry it, tend to have "
               "lower haemoglobin at altitude than Han Chinese at the same height. Its effect within "
               "Tibetans is modest, and how it alters PHD2 is still debated.",
        absent="You do not carry the Tibetan EGLN1 variant.",
        rarity="carried by most Tibetans; about 1 in 20 East Asians; rare elsewhere",
        source="Lorenzo 2014, Nature Genetics",
        pmid="25129147",
        conf="moderate",
    ),
    dict(
        key="a_app_a673t",
        cat="Brain and memory",
        gene="APP",
        title="Protection from Alzheimer's (APP Icelandic)",
        kind="alleles",
        alleles=[
            ("21", 25897620, "C", "T", "p.Ala673Thr (rs63750847)"),
        ],
        copies="one",
        effect="This change sits next to the cut site that releases amyloid beta and cuts its production"
               " by about 40 percent. Carriers are about 5 times less likely to develop Alzheimer's "
               "disease and show less memory decline in old age.",
        absent="You do not carry this variant; your Alzheimer's risk is shaped by other genes and life "
               "factors.",
        rarity="about 1 in 120 Icelanders, Swedes and Finns and 1 in 250 Norwegians; very rare "
               "elsewhere, Denmark included",
        source="Jonsson 2012, Nature",
        pmid="22801501",
        conf="strong",
    ),
    dict(
        key="a_apoe_christchurch",
        cat="Brain and memory",
        gene="APOE",
        title="Alzheimer's resistance (APOE Christchurch)",
        kind="alleles",
        alleles=[
            ("19", 44908756, "C", "A", "p.Arg154Ser, R136S in mature APOE (rs121918393)"),
        ],
        copies="one",
        effect="A Colombian woman with two copies stayed free of memory problems until her seventies, "
               "about 30 years later than expected, despite carrying a gene change that causes "
               "Alzheimer's in middle age. In her family, relatives with that gene change and one copy "
               "of this variant developed memory problems about 5 years later (median age 52, against 47"
               " without it). Its effect in people without that gene change is unknown.",
        absent="You do not carry this variant.",
        rarity="rare; well under 1 percent of people",
        source="Arboleda-Velasquez 2019, Nature Medicine; Quiroz 2024, NEJM",
        pmid="31686034",
        conf="moderate",
    ),
    dict(
        key="a_apoe_jacksonville",
        cat="Brain and memory",
        gene="APOE",
        title="Protection from Alzheimer's (APOE Jacksonville V236E)",
        kind="alleles",
        alleles=[
            ("19", 44909057, "T", "A", "p.Val254Glu, V236E in mature APOE (rs199768005)"),
        ],
        copies="one",
        effect="Carriers have roughly a third of the usual risk of Alzheimer's disease. The change makes"
               " APOE less prone to clumping and less likely to help amyloid build up in the brain.",
        absent="You do not carry this variant.",
        rarity="rare; roughly 1 in 600 people of European ancestry, rarer elsewhere",
        source="Medway 2014, Mol Neurodegener; Le Guen 2022, JAMA Neurology",
        pmid="24607147",
        conf="moderate",
    ),
    dict(
        key="a_apoe_r251g",
        cat="Brain and memory",
        gene="APOE",
        title="Protection from Alzheimer's (APOE R251G)",
        kind="alleles",
        alleles=[
            ("19", 44909101, "C", "G", "p.Arg269Gly, R251G in mature APOE (rs267606661)"),
        ],
        copies="one",
        effect="Found on the APOE e4 background. People who carry e4 together with R251G have roughly "
               "half the Alzheimer's risk of other e4 carriers, cancelling much of the usual e4 effect.",
        absent="You do not carry this variant.",
        rarity="rare; roughly 1 in 1,000 people",
        source="Le Guen 2022, JAMA Neurology",
        pmid="35639372",
        conf="moderate",
    ),
    dict(
        key="a_plcg2_p522r",
        cat="Brain and memory",
        gene="PLCG2",
        title="Protection from dementia (PLCG2 P522R)",
        kind="alleles",
        alleles=[
            ("16", 81908423, "C", "G", "p.Pro522Arg (rs72824905)"),
        ],
        copies="one",
        effect="Slightly boosts an enzyme in the brain's immune cells. Carriers have about a third lower"
               " risk of Alzheimer's disease, lower risk of two other dementias, and were more common "
               "among people who reached 100 with good memory.",
        absent="You do not carry this variant.",
        rarity="about 1 in 50 people of European ancestry; rarer elsewhere",
        source="Sims 2017, Nature Genetics; van der Lee 2019, Acta Neuropathologica",
        pmid="28714976",
        conf="strong",
    ),
    dict(
        key="a_reln_colbos",
        cat="Brain and memory",
        gene="RELN",
        title="Alzheimer's resistance (Reelin COLBOS)",
        kind="alleles",
        alleles=[
            ("7", 103472855, "T", "C", "p.His3447Arg (rs201731543)"),
        ],
        copies="one",
        effect="A Colombian man with this change, who also carried his family's early-onset Alzheimer's "
               "mutation, stayed cognitively intact until 67 and developed dementia at 72, over 20 years"
               " later than usual for that mutation. His sister, who had both, was delayed less, by 12 "
               "to 14 years. The change makes Reelin signal more strongly, which may have helped protect"
               " his brain despite heavy amyloid build-up.",
        absent="You do not carry this variant.",
        rarity="very rare",
        source="Lopera 2023, Nature Medicine",
        pmid="37188781",
        conf="limited",
    ),
    dict(
        key="a_kl_vs",
        cat="Longevity",
        gene="KL",
        title="Klotho KL-VS (common)",
        kind="alleles",
        alleles=[
            ("13", 33054001, "T", "G", "p.Phe352Val (rs9536314)"),
        ],
        copies="one",
        effect="A common klotho variant. People with one copy have slightly higher klotho levels, have "
               "scored a little better on thinking tests and, among APOE e4 carriers aged 60 to 80, show"
               " less amyloid and lower Alzheimer's risk. Two copies do not seem to share this benefit.",
        absent="You do not carry the KL-VS variant; most people do not.",
        rarity="common; about 1 in 4 to 1 in 3 people of European or African ancestry carry one copy; "
               "rare in East Asians",
        source="Arking 2002, PNAS; Dubal 2014, Cell Reports; Belloy 2020, JAMA Neurology",
        pmid="24813892",
        conf="moderate",
    ),
    dict(
        key="a_tyrp1_r93c",
        cat="Appearance",
        gene="TYRP1",
        title="Melanesian blond hair (TYRP1 R93C)",
        kind="alleles",
        alleles=[
            ("9", 12694273, "C", "T", "p.Arg93Cys (rs387907171)"),
        ],
        copies="two",
        effect="The main cause of the natural blond hair seen in about 5 to 10 percent of Solomon "
               "Islanders, who have dark skin. Two copies usually give blond hair; one copy lightens "
               "hair only slightly on average. It arose independently of the genes behind European blond"
               " hair.",
        absent="You do not carry this variant; it is almost unknown outside Oceania.",
        rarity="about 1 in 4 copies of the gene in the Solomon Islands; almost absent elsewhere",
        source="Kenny 2012, Science",
        pmid="22556244",
        conf="strong",
    ),
    dict(
        key="a_smarcad1_noprints",
        cat="Appearance",
        gene="SMARCAD1",
        title="No fingerprints (adermatoglyphia)",
        kind="alleles",
        alleles=[
            ("4", 94253672, "G", "T", "c.1281+665G>T, skin-specific isoform splice site (rs1057519613)"),
            ("4", 94253672, "G", "A", "c.1281+665G>A (rs1057519613)"),
            ("4", 94253673, "T", "C", "c.1281+666T>C (rs1114167276)"),
            ("4", 94253676, "G", "C", "c.1281+669G>C (rs1114167277)"),
        ],
        copies="one",
        effect="These changes break a short, skin-only form of SMARCAD1. Carriers are born without "
               "fingerprints and often sweat less from their hands. In some families the same changes "
               "also bring blisters and tiny cysts at birth (Basan syndrome), or hardened, thinned skin "
               "on the hands and feet and a higher skin cancer risk (Huriez syndrome). It is nicknamed "
               "'immigration delay disease' because fingerprint scanners fail.",
        absent="You do not carry these variants.",
        rarity="found in a few families",
        source="Nousbeck 2011, Am J Hum Genet",
        pmid="21820097",
        conf="moderate",
    ),
    dict(
        key="a_fgf5_lof",
        cat="Appearance",
        gene="FGF5",
        title="Very long eyelashes (FGF5)",
        kind="gene_lof",
        alleles=[],
        copies="two",
        effect="FGF5 tells hair follicles when to stop growing. Two broken copies cause familial "
               "trichomegaly: eyelashes and body hair that keep growing far longer than normal, as in "
               "long-haired cats and angora mice. One copy has no known visible effect.",
        absent="No rare truncating variant in FGF5 was found.",
        rarity="two broken copies are known in a few families",
        source="Higgins 2014, PNAS",
        pmid="24989505",
        conf="limited",
    ),
]
