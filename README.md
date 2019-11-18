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

# Limits
- I don't validate all potential errors at realign alleles
