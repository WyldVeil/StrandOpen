#!/bin/bash
# Setup, once for everyone: the African outgroup for the archaic-DNA step.
# Variants polymorphic among the 1000 Genomes founders from five African
# populations (YRI, ESN, MSL, LWK, GWD), in hmmix's outgroup format, one file
# per chromosome in $KGOUT/archaic. Needs the phased panel ($REF/1kg) and its
# sample sheet; takes about half an hour on 8 cores.
set -euo pipefail
. "$(dirname "$0")/env.sh"
OUT="$KGOUT/archaic"
mkdir -p "$OUT"
cd "$OUT"
PED="$RAW/20130606_g1k_3202_samples_ped_population.txt"
awk 'NR>1 && $3=="0" && $4=="0" && ($6=="YRI"||$6=="ESN"||$6=="MSL"||$6=="LWK"||$6=="GWD") {print $2}' "$PED" > outgroup_samples.txt
echo "outgroup founders: $(wc -l < outgroup_samples.txt)"
one() {
    local c=$1 in="$KG/1kGP_high_coverage_Illumina.chr$1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    [ -s "$in" ] || in="$KG/1kGP_high_coverage_Illumina.chr$1.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    [ -s "og_chr$c.tsv" ] && return 0
    bcftools view -S outgroup_samples.txt --force-samples -v snps -Ou "$in" 2>/dev/null \
      | bcftools norm -m -any -Ou 2>/dev/null \
      | bcftools view -v snps -Ou \
      | bcftools +fill-tags -Ou -- -t AC,AN \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' \
      | awk -F'\t' '$5>0 && $5<$6 {anc = ($6-$5 >= $5) ? $3 : $4; print $1"\t"$2"\t"$3":"($6-$5)"\t"$4":"$5"\t"anc}' > "og_chr$c.tmp"
    mv "og_chr$c.tmp" "og_chr$c.tsv"
    echo "chr$c done"
}
export -f one
export KG
seq 1 22 | xargs -P "${JOBS:-8}" -I{} bash -c 'one {}'
echo "archaic outgroup ready"
