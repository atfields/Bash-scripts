#!/bin/bash

### Create quality graphics for a given VCF file ###

if [ ! -s popmap ]; then
echo "No popmap file present"
exit 1
fi

var1=$(echo $1 | sed 's/.vcf//g'| sed 's/.recode//g')
date=`date +%Y-%m-%d`

## Usage charts.sh VCF_file ##

# Making data files #

mkdir $var1.$date
cp $1 $var1.$date
cd $var1.$date
vcftools --vcf ../$1 --out out --missing-indv
vcftools --vcf ../$1 --out out --missing-site
vcftools --vcf ../$1 --out out --het
vcftools --vcf ../$1 --out out --depth
vcftools --vcf ../$1 --out out --site-mean-depth
vcftools --vcf ../$1 --out out --freq
vcftools --vcf ../$1 --out out --site-quality
vcftools --vcf ../$1 --out out --relatedness
vcftools --vcf ../$1 --out out --site-depth
vcftools --vcf ../$1 --out out --singletons
vcftools --vcf ../$1 --out out --geno-depth
vcftools --vcf ../$1 --out opt --012
echo "CHROM POS N_ALLELES N_CHR FREQ1 FREQ2 FREQ3 FREQ4 FREQ5" > out.fin.frq & echo "$(tail -n +2 out.frq)" >> out.fin.frq
cp ../popmap .
cp ~/bin/plotting.r .
cp ~/bin/opt.script .
cp ~/bin/ddRAD_compare.R .

# Running Rscripts
Rscript plotting.r 2>> error.Rout &
Rscript opt.script 2>> error.Rout &
Rscript ddRAD_compare.R 2>> error.Rout
wait

# Clean-up
rm plotting.r
rm Rplots.pdf
rm opt.script
rm ddRAD_compare.R
