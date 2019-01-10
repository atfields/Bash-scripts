#!/bin/bash
#Script to optimize the mapping of a sample to a given reference through mapping it multiple times with different parameters

#Usage "Opti_mapping.sh <Fwd_fastq.file.gz> <ref.list> <max_A> <max_B> <max_O> <#_processors>"
#The reference list is a list of the names of the references to map

if [ "$#" -lt 2 ]; then
echo "USAGE: Opti_mapping.sh <Fwd_fastq.file.gz> <ref.list> <max_A> <max_B> <max_O> <#_processors>"
exit 1
fi

#Preparing folders
if [ ! -d "bamfiles" ]
then
    mkdir bamfiles
fi

if [ ! -d "logfiles" ]
then
    mkdir logfiles
fi

if [ ! -d "bedfiles" ]
then
    mkdir bedfiles
fi

if [ ! -d "bedcov_files" ]
then
    mkdir bedcov_files
fi

if [ ! -d "summaries/csv" ]; then
    mkdir -p summaries/csv
fi

if [ ! -d "samfiles" ]; then
    mkdir samfiles
fi

#Variables
var=$(echo $1 | sed -e 's/.R1.fq.gz//g')
FASTQ=$1
REF=$2
ALPHA=${3:-10}
BETA=${4:-10}
GAMMA=${5:-10}
DELTA=${6:-10}
LOG=logfiles/
BAM=bamfiles/
BED=bedfiles/
COV=bedcov_files/
SAM=samfiles/
CSV=summaries/csv/

echo "A = "$ALPHA
echo "B = "$BETA
echo "O = "$GAMMA
echo "# of Processors = "$DELTA

#Summary setup
echo "A B       O       Reference      Total   QC      Mapped  Paired #Q60" > $var.bwa_mapping_summary.txt

#Preparing references
echo "preparing references"
for ref in $(cat $REF)
do
samtools faidx $ref
bwa index $ref &>> index.log
done

#The Loop
for i in $(seq 1 $ALPHA)
do
for j in $(seq 1 $BETA)
do
for k in $(seq 1 $GAMMA)
do
for ref in $(cat $REF)
do
echo "mapping $var $i$j$k $ref"

#Mapping
bwa mem $ref $var.R1.fq.gz $var.R2.fq.gz -L 20,5 -t $DELTA -a -M -T 10 -A $i -B $j -O $k -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> $LOG/bwa.$var.$ref.$i$j$k.log | mawk '$6 !~/[2-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -@$DELTA -q 1 -SbT $ref - > $BAM/$var.$ref.$i$j$k.bam 2>$LOG/$var.$ref.$i$j$k.bam.log

#Making Samfile
samtools view $BAM/$var.$ref.$i.$j.$k.bam > $SAM/$var.$ref.$i.$j.$k.sam

#Mapping Mapping Stats
cut  --output-delimiter="," -f 3,4,5 $SAM/$var.$ref.$i.$j.$k.sam > part.1
cut -f 6,6 $SAM/$var.$ref.$i.$j.$k.sam | awk -F"[A-Z]" '{print $1}'> part.2
for a in A C G T; do cut -f 10,10 $SAM/$var.$ref.$i.$j.$k.sam | awk -F"$a" '{print NF-1}' > part.3$a; done
for b in A C G T; do awk -v val='MD' -F" " '{for (c=1; c<=NF; c++) if ($c ~ val) {print $c} }' $SAM/$var.$ref.$i.$j.$k.sam | awk -F"$b" '{print NF-1}'> part.4$b; done
awk -v val='NM:i' -F" " '{for (c=1; c<=NF; c++) if ($c ~ val) {print $c} }' $SAM/$var.$ref.$i.$j.$k.sam | sed 's/NM:i://' > part.5
awk -v val='AS:i' -F" " '{for (c=1; c<=NF; c++) if ($c ~ val) {print $c} }' $SAM/$var.$ref.$i.$j.$k.sam | sed 's/AS:i://' > part.6
paste -d "," part.1 part.2 part.3A part.3C part.3G part.3T part.4A part.4C part.4G part.4T part.5 part.6 > data.csv
echo -e "Ref_contig,start_bp,Map_Q,1st_CIGAR_#,Seq_A,Seq_C,Seq_G,Seq_T,Ref_A,Ref_C,Ref_G,Ref_T,Mismatch,Align_score" | cat - data.csv > $CSV/$var.$ref.$i.$j.$k.summary.csv

#Making Bedfile
samtools sort -@$DELTA $BAM/$var.$ref.$i$j$k.bam -o $BAM/$var.$ref.$i$j$k.bam
samtools index $BAM/$var.$ref.$i$j$k.bam
wait
bamToBed -i $BAM/$var.$ref.$i$j$k.bam | bedtools merge -i - > $BED/$var.$ref.$i$j$k.bed

echo "evaluating $var $i$j$k $ref"
#Summarizing Reads mapped
TOT=$(zgrep @ $var.R1.fq.gz| wc -l)
QC=$(samtools flagstat $BAM/$var.$ref.$i$j$k.bam | grep total  | grep -E -o '[0-9]+' | head -1)
MAP=$(samtools flagstat $BAM/$var.$ref.$i$j$k.bam | grep mapped | grep -E -o '[0-9]+' | head -1)
PAIR=$(samtools flagstat $BAM/$var.$ref.$i$j$k.bam | grep properly | grep -E -o '[0-9]+' | head -1)
Q=$(samtools view $BAM/$var.$ref.$i$j$k.bam | cut -f 5,5 | grep 60 - | wc -l)
echo "$i        $j      $k      $ref      $TOT    $QC     $MAP    $PAIR    $Q" >> $var.bwa_mapping_summary.txt

#Coverage per contig
samtools bedcov  $BED/$var.$ref.$i.$j.$k.bed  $BAM/$var.$ref.$i.$j.$k.bam > $COV/$var.$ref.$i.$j.$k.bedcov
awk '{$5 = ($3 -$2)}1' $COV/$var.$ref.$i.$j.$k.bedcov | awk '{$6 = $4 / $5}1' | awk '!($4="")'> temp
echo $'Contig Start_bp End_bp Length Avg_coverage' | cat - temp > $COV/$var.$ref.$i.$j.$k.bedcov

#Coverage summaries
tail -n +2 bedcov_files/$var.$ref.$i.$j.$k.bedcov | awk '{print$4}' | sort -n | awk -v i=$i -v j=$j -v k=$k -v var=$var 'BEGIN {c=0; sum=0;} $1 ~ /^[0-9]*(\.[0-9]*)?$/ {a[c++] = $1; sum += $1;} END { ave = sum / c; if( (c % 2) == 1 ) {median = a[ int(c/2) ];} else {median = ( a[c/2] + a[c/2-1] ) / 2;} OFS="\t"; print var, i, j, k, ave, median, a[0], a[c-1];}'  >>  $ref.bedcov_len_summary.txt
tail -n +2 bedcov_files/$var.$ref.$i.$j.$k.bedcov | awk '{print$5}' | sort -n | awk -v i=$i -v j=$j -v k=$k -v var=$var 'BEGIN {c=0; sum=0;} $1 ~ /^[0-9]*(\.[0-9]*)?$/ {a[c++] = $1; sum += $1;} END { ave = sum / c; if( (c % 2) == 1 ) {median = a[ int(c/2) ];} else {median = ( a[c/2] + a[c/2-1] ) / 2;} OFS="\t"; print var, i, j, k, ave, median, a[0], a[c-1];}'  >>  $ref.bedcov_cov_summary.txt

done
done
done
done
echo "It is finished"
