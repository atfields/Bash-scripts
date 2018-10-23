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
cut -f 1 out.imiss > tmp
tail -n +2 tmp | sed -r 's/([0-9]+)_/\1\t/g' | sed -r 's/([0-9]+)./\1\t/g' > tmp2
sed -i '1s/^/Sample\tLib\tIndex\tRow\n/' tmp2
paste tmp tmp2 > out.samples
sed -i '1s/^/CHROM POS N_ALLELES N_CHR FREQ1 FREQ2 FREQ3 FREQ4 FREQ5\n/' out.frq
#echo "CHROM POS N_ALLELES N_CHR FREQ1 FREQ2 FREQ3 FREQ4 FREQ5" > out.fin.frq & echo "$(tail -n +2 out.frq)" >> out.fin.frq
cp ../popmap .
cp ~/bin/plotting.r .
cp ~/bin/opt.script .

# Running Rscripts
Rscript plotting.r 2>> error.Rout &
Rscript opt.script 2>> error.Rout &
wait

# Clean-up
rm tmp tmp2
rm plotting.r
rm Rplots.pdf
rm opt.script
