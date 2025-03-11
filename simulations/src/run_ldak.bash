#!/usr/bin/env bash

mkdir -p ldak; cd ldak

## generate the .fam file from the pedigree output
## .fam: one row per ind; 6 columns (id, fam id, mom id, dad id, sex, phenotype)
awk 'BEGIN {OFS="\t"} NR>1 {print $1, "1", $3, $2, $4, $7}' ../pedigree.txt > simu.fam

## generate the .bim file from the SNP_INFO output
## .bim: one row per snp; 6 columns (chr, name, genetic distance, physical distance, A1 allele, A2 allele)
LC_NUMERIC=C awk -F'\t' 'NR==2 {divisor = $3} NR > 1 {print $2, $1, $3, $3/divisor, $4, $5}' OFS='\t' ../SNP_INFO.txt > simu.bim

## generate the .bed file from the genotypes
## first need to generate a file with the following information: SNP genotypes (0, 1, 2, NA) for the A1 allele (row = snp, column = ind)
## and then convert it to .deb using ldak
awk 'NR > 1 { $1=""; print substr($0,2) }' ../genotype.txt | tr -s ' ' '\t' | \
awk '
BEGIN { 
    FS = OFS = "\t" 
}
{
    for (i = 1; i <= NF; i++) {
        a[i] = a[i] $i "\t"
    }
}
END {
    for (i = 1; i in a; i++) {
        print substr(a[i], 1, length(a[i]) - 1)
    }
}
' > simu.txt

ldak6 --make-bed simu2 --gen simu.txt --bim simu.bim --fam simu.fam --gen-skip 0 --gen-headers 0 --gen-probs 1 --threshold 1

## calculate the weights
ldak6 --cut-weights simu_freq --bfile simu2
ldak6 --calc-weights-all simu_freq --bfile simu2 

## calculate the GRM
ldak6 --calc-kins-direct LDAKgmat_freq --bfile simu2 --weights simu_freq/weights.short --power -.25 --kinship-raw YES
