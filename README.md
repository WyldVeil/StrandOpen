# Strand Open

A genome explorer for Linux that works from the raw DNA file you already
have. Give it a consumer chip export (23andMe, AncestryDNA, MyHeritage,
FamilyTreeDNA, LivingDNA and similar) or a VCF/gVCF from whole-genome
sequencing, and it turns it into something you can browse: every variant,
every gene, what ClinVar and the other big databases say, how your genes
affect medicines and traits, ancestry and both family lines, the ancient
people who shared them, and polygenic scores.

Everything runs on your computer. The one-time setup downloads public
reference data; nothing about anyone's genome is sent anywhere.

## Screenshots

The person in these screenshots is a fictional example.

**People**: add anyone's raw DNA file; imports run in the background

![People](docs/screenshots/people.png)

**In brief**: the findings worth knowing first, each with what it means

![In brief](docs/screenshots/brief.png)

**Overview**: counts, chromosomes, findings, lines and scores

![Overview](docs/screenshots/overview.png)

**Genome browser**: genes, variants and coverage at any zoom

![Genome browser](docs/screenshots/browser.png)

**Variants**: every call, filtered by effect, ClinVar class, zygosity or frequency

![Variants](docs/screenshots/variants.png)

**Genes**: the genes worth a look, the ACMG list and the pharmacogenes

![Genes](docs/screenshots/genes.png)

**Rare and remarkable**: notable variants in plain words

![Rare and remarkable](docs/screenshots/remarkable.png)

**Clinical findings**: ClinVar, risk factors, ACMG genes and well-known variants

![Clinical findings](docs/screenshots/clinical.png)

**Carrier status**: recessive conditions carried in one copy

![Carrier status](docs/screenshots/carrier.png)

**Medications**: PharmCAT results with CPIC and DPWG guidance

![Medications](docs/screenshots/medications.png)

**Diet and supplements**: what genes can and cannot say, graded by evidence

![Diet and supplements](docs/screenshots/diet.png)

**Polygenic scores**: placed among people of the same ancestry

![Polygenic scores](docs/screenshots/scores.png)

**Your rarest variants**: how rare each variant is around the world

![Your rarest variants](docs/screenshots/rarest.png)

**Traits**: what well-studied variants say about looks, taste and body

![Traits](docs/screenshots/traits.png)

**Ancestry**: composition, world and European maps, painted chromosomes

![Ancestry](docs/screenshots/ancestry.png)

**Lineage**: where the maternal line comes from, and ancient people on it

![Lineage](docs/screenshots/lineage.png)

**Raw data and quality**: what the file contained and how good it is

![Raw data and quality](docs/screenshots/raw.png)

**Sources and notes**: every data source, and a summary for a doctor

![Sources and notes](docs/screenshots/sources.png)

## Getting started

Download `StrandOpen-<version>-linux-x86_64.tar.gz` from the Releases page,
unpack it anywhere and run `./strandopen` inside (`./add-to-menu.sh` adds it
to your application menu). Everything it needs is in the folder: nothing to
install.

Needs: 64-bit Linux with a desktop (X11, or Wayland with XWayland), 8 GB of
memory, and disk space for the reference data (about 7 GB for the
essentials, 55 GB free while the population reference is built, 19 GB
afterwards).

1. Start Strand Open. The People page opens.
2. Download the reference data (one click). The essentials are about 3 GB
   and cover health, traits, medicines and lineages. The population
   reference (about 33 GB more to download, 12 GB kept) adds ancestry,
   polygenic scores, rarity and Neanderthal DNA. Downloads resume if they
   are interrupted, and carry on if you close the window.
3. Add a person: a name, and their raw data file (drop it on the window or
   choose it). The import runs in the background, usually 10 to 30 minutes.
4. Open them. Each person's results are kept apart; add as many as you like.

Your data folder defaults to `~/.local/share/StrandOpen`; "Change folder" on
the People page starts using another one, for example on a bigger disk (what
is already in the old folder stays there).

## What it reads

| file | notes |
|---|---|
| 23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, LivingDNA, tellmeGen | plain, gzipped or zipped; GRCh37 or GRCh38 positions, detected from the data |
| VCF or gVCF | any provider; GATK, DRAGEN, DeepVariant and bcftools gVCF styles; `chr1` or `1`; the first sample in the file |

GRCh37 files are lifted to GRCh38 with the UCSC chain, strand flips are
checked against 1000 Genomes, and chromosomal sex is read from the data
(you can also state it). A VCF that lists only differences is taken as
matching the reference elsewhere when it covers a whole genome, inside the
1000 Genomes accessibility mask; an exome or panel VCF gets no such
assumption.

## What you get, by kind of data

| | chip | sequencing (VCF/gVCF) |
|---|---|---|
| in brief: the headline findings, each with what it means | yes | yes |
| ClinVar findings, ACMG genes, well-known variants, traits | yes | yes |
| medicines (PharmCAT, CPIC, DPWG, FDA labels, PharmGKB) | yes, for the positions tested | yes |
| diet and supplements: what genes can and cannot say, graded by evidence | yes | yes |
| what each parent gave you (inheritance by sex, X and Y, grandparents) | yes | yes |
| maternal and paternal lines, ancient people on your branches | yes | yes, if the file has MT and Y |
| ancestry, chromosome painting, parents' relatedness | yes | yes |
| polygenic scores | yes, over the variants the chip read | yes |
| remarkable variants (short sleep, low LDL, infection resistance...) | the ones on the chip | yes |
| rare variants, your rarest variants | no: rare chip calls are mostly wrong (Weedon et al. 2021, BMJ) | yes |
| Neanderthal DNA | no | whole genomes |

Chip calls and single sequencing calls can be wrong. This is a tool for
curiosity and for preparing questions, not a medical test: anything that
might matter for health deserves clinical confirmation and a conversation
with a doctor or genetic counsellor.

## Building from source

    ./build.sh                  # build/strandopen (gcc 15+ for #embed)

The program is C (raylib 6.2 and SQLite, both vendored); the import is
Python 3 with numpy, driving bcftools, samtools, plink2, Java (PharmCAT) and
hmmix. Release bundles ship all of these in `runtime/`, so nothing needs
installing; when running from source they are found on your PATH.

    pipeline/setup.py --root DIR [--population]   # the reference data
    pipeline/run.py --profile DIR/profiles/NAME   # one person's import

Headless screenshots for development (never on the real desktop):

    STRANDOPEN_ROOT=DIR tools/shot.sh out.png people 40

## Data sources

GRCh38 and gene models (Ensembl 115), ClinVar (NCBI), AlphaMissense (Cheng
et al. 2023), gnomAD v4.1 constraint, UniProt, the Human Phenotype Ontology,
PharmCAT and PharmGKB/ClinPGx, CPIC and DPWG guidelines, the NHGRI-EBI GWAS
Catalog, the PGS Catalog, the 1000 Genomes 30x panel (NYGC), PhyloTree 17
(rCRS forensic update), the YFull tree with YBrowse positions, the Allen
Ancient DNA Resource, UCSC chain files and bands, and Natural Earth. Each is
downloaded from its publisher at setup and stays under its own licence.

## Licence

Strand Open's own code is MIT licensed (LICENSE). Bundled tools keep their
own licences, listed in `runtime/LICENSES`.
