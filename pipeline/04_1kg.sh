#!/bin/bash
# Setup, once for everyone: the 1000 Genomes reference (NYGC 30x, GRCh38,
# 3,202 people) turned into what the import steps read. Founders only (no
# children of trios). Per chromosome:
#   sites/chrN.bcf  every biallelic variant with its frequency overall and in
#                   each of the five continental groups (rarity, strand checks)
#   pgen/chrN.*     genotypes of variants with a minor allele frequency of at
#                   least 0.5% (PCA, admixture, polygenic score distributions)
# Only the frequencies the steps read are kept, and the genotype files carry
# no annotation column: about 10 GB in all instead of nearly 60.
# Chromosomes are cut into 10 Mb pieces processed in parallel (one stream per
# chromosome is single-threaded), then joined.
#   pipeline/04_1kg.sh            every chromosome
#   pipeline/04_1kg.sh 21 22      just these
set -euo pipefail
. "$(dirname "$0")/env.sh"
P="$KGOUT"
mkdir -p "$P/sites" "$P/pgen" "$P/chunks"
cd "$P"
PED="$RAW/20130606_g1k_3202_samples_ped_population.txt"
awk 'NR>1 && $3=="0" && $4=="0" {print $2}' "$PED" > founders.txt
awk 'NR>1 && $3=="0" && $4=="0" {print $2"\t"$7}' "$PED" > groups.txt
awk 'NR>1 && $3=="0" && $4=="0" {print $2"\t"$2"\t"$6"\t"$7"\t"$5}' "$PED" > founders.pop.tsv
awk 'NR>1 && $3=="0" && $4=="0" {print $2"\t"$5}' "$PED" | sed '1i #IID\tSEX' > founders.sex
for c in $(seq 1 22) X; do echo "chr$c $c"; done > rename.txt
KEEP='^INFO/AF,INFO/AC,INFO/AN,INFO/AF_EUR,INFO/AF_AFR,INFO/AF_EAS,INFO/AF_SAS,INFO/AF_AMR'

piece() {   # piece <chrom> <start> <end>
    local c=$1 a=$2 b=$3 in="$KG/1kGP_high_coverage_Illumina.chr$1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    [ -s "$in" ] || in="$KG/1kGP_high_coverage_Illumina.chr$1.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    local o="chunks/chr$c.$a"
    [ -s "$o.done" ] && return 0
    bcftools view -r "chr$c:$a-$b" -S founders.txt --force-samples -m2 -M2 -v snps,indels -Ou "$in" 2>/dev/null \
      | bcftools view -i "POS>=$a" -Ou \
      | bcftools annotate --rename-chrs rename.txt -x INFO -Ou \
      | bcftools +fill-tags -Ou -- -S groups.txt -t AF,AC,AN 2>/dev/null \
      | bcftools annotate -x "$KEEP" -Ou \
      | tee >(bcftools view -G -Ob -o "$o.sites.bcf") \
      | bcftools view -q 0.005:minor -Ou | bcftools annotate -x INFO -Ob -o "$o.common.bcf"
    sleep 1
    touch "$o.done"
}
export -f piece
export KG KEEP
CS=$(if [ $# -gt 0 ]; then echo "$@"; else echo $(seq 1 22) X; fi)
for c in $CS; do
    [ -s "pgen/chr$c.pgen" ] && [ -s "sites/chr$c.bcf.csi" ] && continue
    len=$(awk -v c="$c" '$1==c{print $2}' "$FA.fai")
    for ((a = 1; a <= len; a += 10000000)); do echo "$c $a $((a + 9999999))"; done
done | xargs -P "${KGJOBS:-$THREADS}" -L 1 bash -c 'piece $0 $1 $2'
for c in $CS; do
    [ -s "pgen/chr$c.pgen" ] && [ -s "sites/chr$c.bcf.csi" ] && continue
    ls chunks/chr$c.*.sites.bcf | sort -t. -k2 -n > chunks/chr$c.sites.list
    ls chunks/chr$c.*.common.bcf | sort -t. -k2 -n > chunks/chr$c.common.list
    # --naive: the pieces share one header, so their blocks are joined as
    # they are (bcftools 1.25's re-encoding concat corrupted these records)
    bcftools concat --naive -f chunks/chr$c.sites.list -Ob -o "sites/chr$c.bcf" && bcftools index -f "sites/chr$c.bcf"
    bcftools concat --naive -f chunks/chr$c.common.list -Ob -o "pgen/chr$c.common.bcf"
    plink2 --bcf "pgen/chr$c.common.bcf" --keep founders.txt --max-alleles 2 \
        --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 60 missing $( [ "$c" = X ] && echo "--update-sex founders.sex --split-par hg38" ) \
        --threads 4 --memory 4000 --make-pgen --out "pgen/chr$c" > "pgen/chr$c.log" 2>&1
    rm -f "pgen/chr$c.common.bcf" chunks/chr$c.*
    echo "chr$c done"
done
ls pgen/*.pgen | sed 's/\.pgen$//' > pgen/merge.txt
ls sites/*.bcf > sites.list
echo "1000 Genomes reference ready"
