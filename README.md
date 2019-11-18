# KRG-vcf
korean reference genome variant table to vcf file

# Overview
The KRGDB (http://coda.nih.go.kr/coda/KRGDB/index.jsp) provides some important statistics of genetic profile from 1722 Korean individuals.
but they provides variant informations via only table format and allele frequency of many indels need to correction.
there are many differences between MAF from "Merged All" and sum of "1st phase and "2nd phase"

# Features
- This repository provides vcf format file of KRGDB statistics table
- calculate MAF using sum of allele counts from 1st phase and 2nd phase (using home-made python script, KRG_phase_merge.py)
- realign reference and position at non-vcf indel allele (ex. AGACG/AGA -> ACG/A, +2 position)
- check reference and filter out errors using bcftools --check-ref with human_g1k_v37.fasta
- add rsid using SnpSift and dbsnp151

# Limits
- I don't validate all potential errors at realign alleles
- I uploaded only chr22.vcf.gz for sample because of size limitation. you can make vcf file using my python scripts.

# add rsid to merged vcf

for chrNo in {1..22}
do
  bgzip -c KRG_merge_song_chr${chrNo}.vcf > KRG_merge_song_chr${chrNo}.vcf.gz
  tabix -p vcf KRG_merge_song_chr${chrNo}.vcf.gz
  bcftools norm --check-ref wx -f human_g1k_v37.fasta -Oz -o KRG_merged_chr${chrNo}.checkref.vcf.gz KRG_merge_song_chr${chrNo}.vcf.gz
  bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' KRG_merge_song_chr${chrNo}.checkref.vcf.gz | bcftools norm --rm-dup none | bcftools annotate -x ID -Ov | SnpSift -Xmx32g annotate -id -noInfo -a -v dbsnp_151.hg19.vcf.gz - | bcftools annotate --set-id +'%CHROM:%POS:%REF:%ALT' -Oz > KRG_merged_chr${chrNo}.checkref.SnpSift.vcf.gz
done

for chrNo in X Y
do
  bgzip -c KRG_merge_song_chr${chrNo}.vcf > KRG_merge_song_chr${chrNo}.vcf.gz
  tabix -p vcf KRG_merge_song_chr${chrNo}.vcf.gz
  bcftools norm --check-ref wx -f human_g1k_v37.fasta -Oz -o KRG_merged_chr${chrNo}.checkref.vcf.gz KRG_merge_song_chr${chrNo}.vcf.gz
  bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' KRG_merge_song_chr${chrNo}.checkref.vcf.gz | bcftools norm --rm-dup none | bcftools annotate -x ID -Ov | SnpSift -Xmx32g annotate -id -noInfo -a -v dbsnp_151.hg19.vcf.gz - | bcftools annotate --set-id +'%CHROM:%POS:%REF:%ALT' -Oz > KRG_merged_chr${chrNo}.checkref.SnpSift.vcf.gz
done
