fastqc <file.fastq>
unzip <file_fastqc.zip>
awk '/Overrepresented sequences/,/END/' <file_fastqc>/fastqc_data.txt | tail -n +3 | head -n -1 | cut -f 1 > tmp
awk ' BEGIN{print ">Overrepresented sequence"}{print;} NR % 1 == 0 { print ">Overrepresented sequence"; }' tmp | head -n -1 > tmp2
cat /usr/local/bin/TruSeq2-PE.fa tmp2 > adator_trim.fa
rm tmp tmp2
java -jar /usr/local/bin/trimmomatic-0.33.jar PE -threads <Num_Threads> -phred33 <sample>.F.fq.gz <sample>.R.fq.gz <sample>.R1.fq.gz <sample>.unpairedF.fq.gz <sample>.R2.fq.gz <sample>.unpairedR.fq.gz ILLUMINACLIP:./adator_trim.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:10 2>&1 | tee -a <sample>.trim.log
