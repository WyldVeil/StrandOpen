# Paths shared by every shell step, set by pipeline/run.py:
#   STRANDOPEN_PROFILE  the person's folder (canon/, work/, db/)
#   STRANDOPEN_REF      shared reference downloads; STRANDOPEN_KG 1000 Genomes products
: "${STRANDOPEN_PROFILE:?run the import through pipeline/run.py}"
PROFILE="$STRANDOPEN_PROFILE"
REF="${STRANDOPEN_REF:-$(dirname "$(dirname "$PROFILE")")/shared/ref}"
KGOUT="${STRANDOPEN_KG:-$(dirname "$(dirname "$PROFILE")")/shared/1kg}"
RAW="$REF/raw"
FA="$REF/GRCh38.fa"
KG="$REF/1kg"
WORK="$PROFILE/work"
DB="$PROFILE/db"
GVCF="$PROFILE/canon/calls.g.vcf.gz"
CRAM="$PROFILE/align/reads.cram"
SVVCF="$PROFILE/canon/sv.vcf.gz"
CNVVCF="$PROFILE/canon/cnv.vcf.gz"
CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
THREADS="${THREADS:-$(nproc)}"
mkdir -p "$WORK" "$DB"
