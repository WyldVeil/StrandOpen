"""Remarkable variants, part C: mind, behaviour, habits and body size,
including some famous single-gene claims that did not hold up.

Every allele is GRCh38, 1-based VCF style, and its REF checked against the
FASTA; positions come from dbSNP and protein changes were checked against
the Ensembl 115 MANE Select CDS. Common alleles were cross-checked against
the 1000 Genomes sites files.

conf 'disputed' marks a claim that was widely reported but that large or
pooled studies did not confirm; the effect text says both. Repeat
polymorphisms (DRD4 7R, 5-HTTLPR, the MAOA promoter repeat) cannot be read
from these variant calls and are left out.
"""

ENTRIES = [
    dict(
        key="c_setd1a_lof",
        cat="Mind and mood",
        gene="SETD1A",
        title="Raised schizophrenia risk (SETD1A loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="SETD1A marks the packaging around DNA to switch genes on. Broken copies are among the "
               "strongest single genetic risk factors for schizophrenia known, raising the odds roughly "
               "20-fold in the largest study, and many carriers also have learning difficulties. "
               "Schizophrenia affects under 1 in 100 people, so many carriers never develop it.",
        absent="No protein-truncating change was found in SETD1A.",
        rarity="extremely rare: about 1 in 30,000 people",
        source="Singh 2016, Nat Neurosci; Singh 2022 (SCHEMA), Nature",
        pmid="26974950",
        conf="strong",
    ),
    dict(
        key="c_cul1_lof",
        cat="Mind and mood",
        gene="CUL1",
        title="Raised schizophrenia risk (CUL1 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="CUL1 is the backbone of a protein complex that tags other proteins for disposal. In an "
               "exome study of 24,000 people with schizophrenia and 97,000 without, broken copies turned"
               " up in 8 patients and 1 control: an estimated 36-fold rise in odds, but with very wide "
               "uncertainty. A carrier may well never develop the condition.",
        absent="No protein-truncating change was found in CUL1.",
        rarity="extremely rare: roughly 1 in 10,000 to 100,000 people",
        source="Singh 2022 (SCHEMA), Nature",
        pmid="35396579",
        conf="moderate",
    ),
    dict(
        key="c_xpo7_lof",
        cat="Mind and mood",
        gene="XPO7",
        title="Raised schizophrenia risk (XPO7 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="XPO7 ferries proteins out of the cell nucleus. In an exome study of 24,000 people with "
               "schizophrenia and 97,000 without, broken copies turned up in 12 patients and 1 control: "
               "an estimated 50-fold rise in odds, but with very wide uncertainty. A carrier may well "
               "never develop the condition.",
        absent="No protein-truncating change was found in XPO7.",
        rarity="extremely rare: roughly 1 in 10,000 to 100,000 people",
        source="Singh 2022 (SCHEMA), Nature",
        pmid="35396579",
        conf="moderate",
    ),
    dict(
        key="c_grin2a_lof",
        cat="Mind and mood",
        gene="GRIN2A",
        title="Raised schizophrenia risk (GRIN2A loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="GRIN2A makes part of the NMDA receptor, a key learning switch between nerve cells. "
               "Broken copies turned up in 9 people with schizophrenia and 2 controls in the largest "
               "exome study (about 18 times the odds), and common variants near the gene also affect "
               "risk. Broken copies can also cause childhood epilepsy with speech problems, though many "
               "carriers are unaffected.",
        absent="No protein-truncating change was found in GRIN2A.",
        rarity="very rare: roughly 1 in 5,000 to 50,000 people",
        source="Singh 2022 (SCHEMA), Nature",
        pmid="35396579",
        conf="moderate",
    ),
    dict(
        key="c_akap11_lof",
        cat="Mind and mood",
        gene="AKAP11",
        title="Raised bipolar disorder risk (AKAP11 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="AKAP11 anchors a signalling enzyme inside nerve cells and binds GSK3B, the suspected "
               "target of the mood stabiliser lithium. Broken copies were found in 33 of 38,000 people "
               "with bipolar disorder or schizophrenia but only 13 of 106,000 controls, about 7 times "
               "the odds. Most carriers would still not be expected to develop either condition.",
        absent="No protein-truncating change was found in AKAP11.",
        rarity="very rare: about 1 in 5,000 to 8,000 people",
        source="Palmer 2022, Nat Genet",
        pmid="35410376",
        conf="moderate",
    ),
    dict(
        key="c_cacna1c_rs1006737",
        cat="Mind and mood",
        gene="CACNA1C",
        title="Slightly higher bipolar and schizophrenia risk (CACNA1C)",
        kind="alleles",
        alleles=[
            ("12", 2236129, "G", "A", "rs1006737, intron"),
        ],
        copies="one",
        effect="CACNA1C makes a calcium channel found in heart and brain cells, the target of some blood"
               " pressure drugs. This was one of the first genome-wide findings for bipolar disorder and"
               " is also linked to schizophrenia, but each copy raises the odds by only about 10 to 18 "
               "percent. Smaller studies linking it to differences in brain activity have given mixed "
               "results.",
        absent="You do not carry the CACNA1C rs1006737 risk version.",
        rarity="common: about half of Europeans carry at least one copy",
        source="Ferreira 2008, Nat Genet; Ripke 2013, Nat Genet",
        pmid="18711365",
        conf="moderate",
    ),
    dict(
        key="c_maoa_lof",
        cat="Mind and mood",
        gene="MAOA",
        title="Brunner syndrome (MAOA loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="MAOA breaks down serotonin, noradrenaline and dopamine. In a Dutch family studied in "
               "1993, men with no working MAOA had mild learning disability and bouts of impulsive "
               "aggression. The gene is on the X chromosome, so one copy matters in men; women with one "
               "changed copy have usually been unaffected. The 'warrior gene' is a different, common "
               "repeat this test cannot read; it predicts little for any one person.",
        absent="No protein-truncating change was found in MAOA.",
        rarity="very rare: only a few affected families are known",
        source="Brunner 1993, Science; Piton 2014, Eur J Hum Genet",
        pmid="8211186",
        conf="moderate",
    ),
    dict(
        key="c_ecm1_fear",
        cat="Mind and mood",
        gene="ECM1",
        title="Loss of fear in Urbach-Wiethe disease (ECM1)",
        kind="gene_lof",
        alleles=[],
        copies="two",
        effect="Two broken copies cause lipoid proteinosis (Urbach-Wiethe disease), with a hoarse voice "
               "from infancy and thickened skin. In some patients it slowly turns both amygdalae, the "
               "brain's fear centres, to calcium. The best-studied of them, known as S.M., felt almost "
               "no fear of snakes, spiders or a haunted house. One broken copy has no known effect.",
        absent="No protein-truncating change was found in ECM1.",
        rarity="two copies are very rare; in Namaqualand, South Africa, about 1 in 9 people carry one",
        source="Hamada 2002, Hum Mol Genet; Feinstein 2011, Curr Biol",
        pmid="21167712",
        conf="limited",
    ),
    dict(
        key="c_foxp2_lof",
        cat="Mind and mood",
        gene="FOXP2",
        title="Developmental speech disorder (FOXP2 loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="One broken copy of FOXP2 causes childhood apraxia of speech: trouble planning the mouth "
               "movements of speech, often with wider language difficulties. Found in the KE family in "
               "2001, it was hailed as 'the language gene', but it is one of many genes involved, and a "
               "claimed recent evolutionary sweep in humans was not confirmed in 2018.",
        absent="No protein-truncating change was found in FOXP2.",
        rarity="very rare: a few dozen affected families are known",
        source="Lai 2001, Nature; Atkinson 2018, Cell",
        pmid="11586359",
        conf="strong",
    ),
    dict(
        key="c_gnptab_e1200k",
        cat="Mind and mood",
        gene="GNPTAB",
        title="Persistent stuttering (GNPTAB p.Glu1200Lys)",
        kind="alleles",
        alleles=[
            ("12", 101753376, "C", "T", "p.Glu1200Lys (rs137853825)"),
        ],
        copies="one",
        effect="GNPTAB helps tag enzymes for delivery to the lysosome, the cell's recycling unit. This "
               "change ran in about 1 in 10 Pakistani families with inherited stuttering, and mice given"
               " it make altered pup calls. Many carriers never stutter. Fully broken copies of the gene"
               " cause a severe childhood disease, but this mild change does not.",
        absent="You do not carry the GNPTAB p.Glu1200Lys variant.",
        rarity="found mainly in South Asia, where roughly 1 to 5 in 100 people carry it; very rare "
               "elsewhere",
        source="Kang 2010, NEJM; Barnes 2016, Curr Biol",
        pmid="20147709",
        conf="limited",
    ),
    dict(
        key="c_aldh2_rs671",
        cat="Habits and addiction",
        gene="ALDH2",
        title="Alcohol flush (ALDH2*2, p.Glu504Lys)",
        kind="alleles",
        alleles=[
            ("12", 111803962, "G", "A", "p.Glu504Lys, ALDH2*2 (rs671)"),
        ],
        copies="one",
        effect="This change all but disables the enzyme that clears acetaldehyde, a toxic product of "
               "alcohol. One copy causes flushing, a racing heart and nausea after a drink, and carriers"
               " tend to drink less, with a much lower risk of alcohol dependence. Carriers who drink "
               "regularly anyway have a much higher risk of oesophageal cancer.",
        absent="You do not carry the ALDH2*2 flush variant.",
        rarity="common in East Asia (about 1 in 3 people carry it); very rare elsewhere",
        source="Luczak 2006, Psychol Bull; Brooks 2009, PLoS Med",
        pmid="19320537",
        conf="strong",
    ),
    dict(
        key="c_chrna5_d398n",
        cat="Habits and addiction",
        gene="CHRNA5",
        title="Heavier smoking (CHRNA5 p.Asp398Asn)",
        kind="alleles",
        alleles=[
            ("15", 78590583, "G", "A", "p.Asp398Asn (rs16969968)"),
        ],
        copies="one",
        effect="This change makes a nicotine receptor in the brain less responsive. In people who smoke,"
               " each copy adds about one cigarette a day on average and raises the risk of nicotine "
               "dependence; each copy is also linked to about 30 percent higher odds of lung cancer. It "
               "does not affect whether people start smoking.",
        absent="You do not carry the CHRNA5 p.Asp398Asn variant.",
        rarity="common: about 6 in 10 Europeans carry at least one copy; rarer in Africa and East Asia",
        source="Thorgeirsson 2008, Nature; TAG Consortium 2010, Nat Genet",
        pmid="18385739",
        conf="strong",
    ),
    dict(
        key="c_oprm1_n40d",
        cat="Habits and addiction",
        gene="OPRM1",
        title="Slightly lower opioid addiction risk (OPRM1 p.Asn40Asp)",
        kind="alleles",
        alleles=[
            ("6", 154039662, "A", "G", "p.Asn40Asp, A118G (rs1799971)"),
        ],
        copies="one",
        effect="OPRM1 makes the mu-opioid receptor, the main target of morphine and other opioid "
               "painkillers. After decades of mixed results, very large studies of opioid-exposed people"
               " of European ancestry confirmed that each copy lowers the odds of opioid use disorder by"
               " about 6 to 11 percent, a small shift next to other factors. No link was seen in people "
               "of African ancestry, where it is rare, and other groups are less studied.",
        absent="You do not carry the OPRM1 p.Asn40Asp variant.",
        rarity="common: about 3 in 10 Europeans and 6 in 10 East Asians carry at least one copy",
        source="Zhou 2020, JAMA Psychiatry; Kember 2022, Nat Neurosci",
        pmid="32492095",
        conf="moderate",
    ),
    dict(
        key="c_cadm2_cautious",
        cat="Habits and addiction",
        gene="CADM2",
        title="Slightly more cautious (CADM2 risk-taking locus)",
        kind="alleles",
        alleles=[
            ("3", 85550770, "C", "G", "rs6783138 G, intron"),
        ],
        copies="one",
        effect="CADM2 helps nerve cells connect, and it is the strongest single genetic signal for self-"
               "reported risk taking in studies of over a million people; the region is also linked to "
               "speeding, drinking, cannabis use and body weight. This version, the more common one in "
               "Europeans, goes with very slightly lower risk tolerance, explaining well under 0.1 "
               "percent of differences between people.",
        absent="You carry two copies of the other, slightly more risk-tolerant CADM2 version.",
        rarity="very common: about 85 percent of Europeans carry at least one copy",
        source="Karlsson Linner 2019, Nat Genet",
        pmid="30643258",
        conf="moderate",
    ),
    dict(
        key="c_fto_rs9939609",
        cat="Body size",
        gene="FTO",
        title="Slightly higher weight (FTO rs9939609)",
        kind="alleles",
        alleles=[
            ("16", 53786615, "T", "A", "rs9939609, intron"),
        ],
        copies="one",
        effect="The first common variant found for body weight. Each copy adds on average about 1 to 1.5"
               " kg and raises the odds of obesity by about 20 to 30 percent, and carriers tend to eat "
               "slightly more. In physically active adults the effect is about a quarter smaller. It "
               "seems to act through nearby genes (IRX3 and IRX5) rather than FTO itself.",
        absent="You do not carry the FTO rs9939609 weight-raising version.",
        rarity="common: about 2 in 3 Europeans carry at least one copy",
        source="Frayling 2007, Science; Kilpelainen 2011, PLoS Med",
        pmid="17434869",
        conf="strong",
    ),
    dict(
        key="c_fbn1_e1297g",
        cat="Body size",
        gene="FBN1",
        title="Shorter stature (FBN1 p.Glu1297Gly, Peruvian variant)",
        kind="alleles",
        alleles=[
            ("15", 48481729, "T", "C", "p.Glu1297Gly (rs200342067)"),
        ],
        copies="one",
        effect="Found almost only in people with Native American ancestry, most often in Peru. Each copy"
               " reduces adult height by about 2.2 cm (4.4 cm with two), the largest effect known for a "
               "common height variant. FBN1 builds the elastic fibres of connective tissue; other "
               "changes in it cause Marfan syndrome, which often makes people tall. The variant shows "
               "signs of natural selection.",
        absent="You do not carry the FBN1 p.Glu1297Gly variant.",
        rarity="about 1 in 10 Peruvians carry one copy; almost absent elsewhere",
        source="Asgari 2020, Nature",
        pmid="32499652",
        conf="moderate",
    ),
    dict(
        key="c_ghr_laron_e180",
        cat="Body size",
        gene="GHR",
        title="Laron syndrome: very short, little cancer or diabetes (GHR)",
        kind="alleles",
        alleles=[
            ("5", 42699978, "A", "G", "c.594A>G, E180 splice (rs121909360)"),
        ],
        copies="two",
        effect="Two copies disable the growth hormone receptor, causing Laron syndrome: very short "
               "stature (often under about 1.3 m as adults) with normal intelligence. In a 22-year study"
               " in Ecuador, people with two copies had one non-fatal cancer and no diabetes, compared "
               "with 17 and 5 percent of controls. One copy has little or no effect on height.",
        absent="You do not carry the Ecuadorian GHR E180 splice variant.",
        rarity="found mainly in southern Ecuador, where about 100 people with two copies are known; very"
               " rare elsewhere",
        source="Guevara-Aguirre 2011, Sci Transl Med; Rosenbloom 1990, NEJM",
        pmid="21325617",
        conf="strong",
    ),
    dict(
        key="c_acan_lof",
        cat="Body size",
        gene="ACAN",
        title="Short stature and early arthritis (ACAN loss of function)",
        kind="gene_lof",
        alleles=[],
        copies="one",
        effect="Aggrecan is the main cushioning protein of cartilage, including the growth plates where "
               "bones lengthen. In families with one broken copy, adults were a median of about 18 to 20"
               " cm shorter than average (from near normal to very short), as growth stopped early, and "
               "early-onset osteoarthritis and spinal disc problems were common.",
        absent="No protein-truncating change was found in ACAN.",
        rarity="rare: roughly 1 in 3,000 people or fewer",
        source="Nilsson 2014; Gkourogianni 2017, J Clin Endocrinol Metab",
        pmid="27870580",
        conf="strong",
    ),
    dict(
        key="c_comt_v158m",
        cat="Famous claims",
        gene="COMT",
        title="The 'worrier' COMT version (Val158Met, two Met copies)",
        kind="alleles",
        alleles=[
            ("22", 19963748, "G", "A", "p.Val158Met (rs4680)"),
        ],
        copies="two",
        effect="The Met version makes a less stable COMT enzyme, so dopamine lingers longer in the front"
               " of the brain; that part is solid. It was then claimed that Met carriers are 'worriers' "
               "with sharper working memory but more anxiety and pain, and Val carriers 'warriors' who "
               "cope better under stress. Large pooled studies found little or no effect on thinking or "
               "depression.",
        absent="You have two copies of the Val version, the so-called 'warrior' type.",
        rarity="very common: about 1 in 4 Europeans have two Met copies",
        source="Egan 2001, PNAS; Barnett 2008, Biol Psychiatry; Border 2019",
        pmid="18339359",
        conf="disputed",
    ),
    dict(
        key="c_bdnf_v66m",
        cat="Famous claims",
        gene="BDNF",
        title="BDNF Val66Met and memory",
        kind="alleles",
        alleles=[
            ("11", 27658369, "C", "T", "p.Val66Met (rs6265)"),
        ],
        copies="one",
        effect="The Met version reduces how much BDNF, a growth factor for nerve cells, is released when"
               " they fire. A 2003 study linked it to poorer memory and many later ones to mood, but a "
               "pooled analysis of 7,000 people found no effect on memory or thinking, and depression "
               "studies of over 60,000 people found none. It does slightly lower body weight and the "
               "chance of taking up smoking.",
        absent="You do not carry the BDNF Met version.",
        rarity="common: about 1 in 3 Europeans and 3 in 4 East Asians carry at least one copy",
        source="Egan 2003, Cell; Mandelman 2012, Genes Brain Behav; Border 2019",
        pmid="21980924",
        conf="disputed",
    ),
    dict(
        key="c_oxtr_rs53576",
        cat="Famous claims",
        gene="OXTR",
        title="The 'empathy' version of the oxytocin receptor (OXTR rs53576)",
        kind="alleles",
        alleles=[
            ("3", 8762685, "A", "G", "rs53576 G, intron"),
        ],
        copies="two",
        effect="A small 2009 study found that people with two G copies scored higher on empathy tests "
               "and showed a smaller stress response than people with an A, and the variant became known"
               " as an empathy or 'sociability gene'. A pooled analysis of over 17,000 people found no "
               "effect on personality, social behaviour, mental health or autism.",
        absent="You have no G copy at OXTR rs53576.",
        rarity="very common: about 4 in 10 Europeans have two G copies",
        source="Rodrigues 2009, PNAS; Bakermans-Kranenburg 2014, Psychiatr Genet",
        pmid="23921259",
        conf="disputed",
    ),
    dict(
        key="c_kibra_rs17070145",
        cat="Famous claims",
        gene="WWC1",
        title="The KIBRA 'memory gene' (WWC1 rs17070145)",
        kind="alleles",
        alleles=[
            ("5", 168418786, "C", "T", "rs17070145 T, intron"),
        ],
        copies="one",
        effect="A 2006 genome scan reported that people carrying the T version remembered word lists "
               "better, and KIBRA was hailed as a memory gene. Later results were mixed: a pooled "
               "analysis of about 9,000 people found an effect so small it explains about 0.5 percent of"
               " differences in memory, and large genome-wide studies of memory and thinking have not "
               "singled it out.",
        absent="You do not carry the KIBRA T version.",
        rarity="common: about 6 in 10 Europeans carry at least one T copy",
        source="Papassotiropoulos 2006, Science; Milnik 2012, Am J Med Genet B",
        pmid="23065961",
        conf="disputed",
    ),
    dict(
        key="c_drd2_taq1a",
        cat="Famous claims",
        gene="ANKK1",
        title="The dopamine 'reward gene' (DRD2 Taq1A A1)",
        kind="alleles",
        alleles=[
            ("11", 113400106, "G", "A", "Taq1A A1, ANKK1 p.Glu713Lys (rs1800497)"),
        ],
        copies="one",
        effect="In 1990 the A1 version was found in 77 percent of brain samples from people with alcohol"
               " dependence against 28 percent of others, launching 'reward deficiency' theories. It "
               "lies beside the dopamine D2 receptor gene, and carriers may have slightly fewer D2 "
               "receptors. Pooled studies found at most a small effect with signs of publication bias, "
               "and genome-wide studies of problem drinking flag other variants nearby.",
        absent="You do not carry the Taq1A A1 version.",
        rarity="common: about 1 in 3 Europeans carry at least one A1 copy",
        source="Munafo 2007, Mol Psychiatry; Blum 1990, JAMA; Zhou 2020, Nat Neurosci",
        pmid="17453061",
        conf="disputed",
    ),
]
