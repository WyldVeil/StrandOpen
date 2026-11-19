#!/bin/bash
# Step 1: the sample's variant calls, cleaned and annotated.
#   - keep the primary chromosomes, drop reference-only records
#   - split multi-allelic records and left-align against GRCh38
#   - consequence calls from Ensembl gene models (bcftools csq)
#   - ClinVar significance, conditions, review status and variation id
# Output: work/variants.annot.vcf.gz (+ index)
set -euo pipefail
. "$(dirname "$0")/env.sh"
cd "$WORK"

# Ensembl GFF3 for csq: primary chromosomes only
if [ ! -s genes.gff3.gz ]; then
    zcat "$RAW/Homo_sapiens.GRCh38.115.gff3.gz" \
      | awk -F'\t' 'BEGIN{ok["1"];for(i=2;i<=22;i++)ok[i""];ok["X"];ok["Y"];ok["MT"]} /^#/ || ($1 in ok)' \
      | bgzip -@4 > genes.gff3.gz
fi

# ClinVar as a tab file, so its variation id lands in INFO/CLNVID and the
# sample's own rsID column is left alone
if [ ! -s clinvar.tab.gz ]; then
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\t%INFO/CLNSIGCONF\n' \
        "$RAW/clinvar.vcf.gz" | awk -F'\t' '$4!="."' | bgzip -@4 > clinvar.tab.gz
    tabix -f -s1 -b2 -e2 clinvar.tab.gz
fi
cat > clinvar.hdr <<'EOF'
##INFO=<ID=CLNVID,Number=1,Type=String,Description="ClinVar variation id">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar clinical significance">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar preferred condition names">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status">
##INFO=<ID=CLNGENE,Number=.,Type=String,Description="ClinVar gene info">
##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="ClinVar conflicting interpretations">
EOF

REGIONS=$(echo $CHROMS | tr ' ' ',')
bcftools view -t "$REGIONS" -i 'GT="alt"' -Ou "$GVCF" \
  | bcftools norm -m-any -f "$FA" -c w -Ou 2> norm.log \
  | bcftools view -i 'GT="alt"' -Ou \
  | bcftools csq -f "$FA" -g genes.gff3.gz -p a --ncsq 64 -Ou 2> csq.log \
  | bcftools annotate -a clinvar.tab.gz -h clinvar.hdr \
        -c CHROM,POS,REF,ALT,INFO/CLNVID,INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNGENE,INFO/CLNSIGCONF \
        -Oz -o variants.annot.vcf.gz
bcftools index -f -t variants.annot.vcf.gz
echo "variants: $(bcftools index -n variants.annot.vcf.gz)"
