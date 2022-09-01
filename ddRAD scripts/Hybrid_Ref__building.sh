### Making a script for mapping to the reference ###
#Getting the directory setup
cd /home/afields/Workspace/King_Mack/kdye/Ref_70
mkdir -p data/logfiles reference
cd data
cat /home/shared/King_Mack/sample_groups/King_Mack_70.txt | while read i; do cp -s /home/shared/King_Mack/hiseq/$i.F.fq.gz /home/shared/King_Mack/hiseq/$i.R.fq.gz .; done

# Setting Variables
K1_MAX=2
K2_MAX=2
DATA_DIR="/home/afields/Workspace/King_Mack/kdye/Ref_70/data"
NUMProc=10

# Parameters and functions from dDocent
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
AWK4='{for(i=0;i<$1;i++)print}'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
SED1='s/^[ \t]*//'
SED2='s/\s/\t/g'
FRL=$(gunzip -c ${NAMES[0]}.F.fq.gz | mawk '{ print length() | "sort -rn" }' | head -1)

special_uniq(){
    mawk -v x=$1 '$1 >= x' $2  |cut -f2 | sed -e 's/NNNNNNNNNN/     /g' | cut -f1 | uniq
}
export -f special_uniq

pmerge(){
    num=$( echo $1 | sed -e 's/^0*//g')
    if [ "$num" -le 100 ]; then
        j=$num
        k=$(($num -1))
    else
        num=$(($num - 99))
        j=$(python -c "print ("$num" * 100)")
        k=$(python -c "print ("$j" - 100)")
    fi

    mawk -v x="$j" -v y="$k" '$5 <= x && $5 > y'  rbdiv.out | cut -f-5 > rbdiv.out.$1
	 

    if [ -s "rbdiv.out.$1" ]; then
        rainbow merge -o rbasm.out.$1 -a -i rbdiv.out.$1 -r 2 -N10000 -R10000 -l 20 -f 0.75
    fi
}
export -f pmerge

#Finding overlapping reads
cd $DATA_DIR
ls *.F.fq.gz | sed 's/.F.fq.gz//g' > namelist
cat namelist | while read i; do
echo "Processing $i"
pearRM -f $i.F.fq.gz -r $i.R.fq.gz -m 290 -n 97 -o $i -j 20 1> logfiles/${i}_pear.log
done

#Gathering unique sequences
cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.assembled.fastq | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs"
cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.unassembled.forward.fastq | mawk '$AWK2' > {}.forward"
cat namelist | parallel --no-notice -j $NUMProc "mawk '$AWK1' {}.unassembled.reverse.fastq | mawk '$AWK2' > {}.reverse"
cat namelist | parallel --no-notice -j $NUMProc "paste -d '-' {}.forward {}.reverse | mawk '$AWK3'| sed -e 's/-/NNNNNNNNNN/' | perl -e '$PERLT' >> {}.uniq.seqs"
rm *.forward
rm *.reverse

cd ..
#if [ ! -f "uniq.seqs" ]; then cat *.uniq.seqs > uniq.seqs

#Preparing the directories
for i in $(seq 1 $K1_MAX); do
for j in $(seq 1 $K1_MAX); do
if [ $i -eq 1 ] && [ $j -eq 1 ]; then continue; fi
for k in $(seq 0.8 0.02 0.98); do
#Making a directory for each combination
mkdir ref_${i}_${j}_${k};
#Moving into the directory
cd ref_${i}_${j}_${k};
#Soft copying the data into the directory
cp -s $DATA_DIR/*.uniq.seqs .;
cd ..;
done
done
done

#Preparing the references
for i in $(seq 1 $K1_MAX); do
    for j in $(seq 1 $K1_MAX); do
        if [ $i -eq 1 && $j -eq 1 ]; then continue; fi
	    for k in $(seq 0.8 0.02 0.98); do
    		cd ref_${i}_${j}_${k}
          #Create a data file with the number of unique sequences and the number of occurrences(straight from ReferenceOpt.sh or dDocent)
            parallel --no-notice mawk -v x=$i \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' | mawk -v x=$j '$1 >= x' > uniq.k.$i.c.$j.seqs
			cut -f2 uniq.k.$i.c.$j.seqs > totaluniqseq
			mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' totaluniqseq > uniq.full.fasta
			LENGTH=$(mawk '!/>/' uniq.full.fasta  | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
			LENGTH=$(($LENGTH * 3 / 4))
			awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+""\n"gensub(/./, "I", "g", $2)}' uniq.full.fasta > uniq.fq
			fastp -i uniq.fq -o uniq.fq1 -w $NUMProc -Q &> assemble.trim.log
			mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' uniq.fq1 | paste - - | sort -k1,1 -V | tr "\t" "\n" > uniq.fasta
            mawk '!/>/' uniq.fasta > totaluniqseq
            rm uniq.fq*
          #CDHit used to cluster reads
            sed -e 's/NNNNNNNNNN/   /g' totaluniqseq | cut -f1 | sort --parallel=$NUMProc -S 2G | uniq | mawk '{c= c + 1; print ">dDocent_Contig_" c "\n" $1}' > uniq.F.fasta
            CDHIT=$(python -c "print (max("$k" - 0.1,0.8))")
            cd-hit-est -i uniq.F.fasta -o xxx -c $CDHIT -T $NUMProc -M 0 -g 1 -d 100 &>cdhit.log
            mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed -e 's/[>dDocent_Contig_,...]//g' | sort -g -k1 -S 2G --parallel=$NUMProc > sort.contig.cluster.ids
            paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
#            paste sort.contig.cluster.ids <(mawk '!/>/' uniq.F.fasta) > contig.cluster.totaluniqseq
#            sed -e 's/NNNNNNNNNN/   /g' totaluniqseq | sort --parallel=$NUMProc -k1 -S 20G | mawk '{print $0 "\t" NR}'  > totaluniqseq.CN
#            join -t $'\t' -1 3 -2 1 contig.cluster.Funiq totaluniqseq.CN -o 2.3,1.2,2.1,2.2 > contig.cluster.totaluniqseq
          #CD-hit output is converted to rainbow format
            sort -k2,2 -g contig.cluster.totaluniqseq -S 2G --parallel=$NUMProc | sed -e 's/NNNNNNNNNN/    /g' > rcluster
		  #Running rainbow for assembly
            rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10
            CLUST=(`tail -1 rbdiv.out | cut -f5`)
            CLUST1=$(( $CLUST / 100 + 1))
            CLUST2=$(( $CLUST1 + 100 ))

#Input File Format: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>[\t<pre_cluster_id:int>]
#            seq -w 1 $CLUST2 | parallel --no-notice -j $NUMProc --env pmerge pmerge {}
            rainbow merge -o rbasm.out -a -i rbdiv.out -r 2 -N10000 -R10000 -l 20 -f 0.75
#            cat rbasm.out.[0-9]* > rbasm.out
#            rm rbasm.out.[0-9]* rbdiv.out.[0-9]*
			
          #This AWK code replaces rainbow's contig selection perl script
            LENGTH=$(cut -f3 rbdiv.out |mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')

            LENGTH=$(( $LENGTH * 11 / 10 ))

            cat rbasm.out <(echo "E") |sed -e 's/[0-9]*:[0-9]*://g' | mawk -v mlen=$LENGTH  '{
                if (NR == 1) e=$2;
                else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_A_Contig_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
                else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
                else if ($1 ~/C/) clus=$2;
                else if ($1 ~/L/) len=$2;
                else if ($1 ~/S/) seq=$2;
                else if ($1 ~/N/) freq=$2;
                else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
                else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/ && $0 ~/^R 0/ && len <= mlen) {seq1 = seq; fclus=clus;lenf=len}
                else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/ && $0 ~!/^R 0/ && len > mlen) {seq1 = seq; fclus=clus; len1=len}
                else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/ && $0 ~!/^R 0/ && len <= mlen) {seq1 = seq; fclus=clus; lenf=len}
                else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
                }' > rainbow.fasta

            seqtk seq -r rainbow.fasta > rainbow.RC.fasta
            mv rainbow.RC.fasta rainbow.fasta

          #The rainbow assembly is checked for overlap between newly assembled Forward and Reverse reads using the software PEAR
            grep -A1 "dDocent_A_Contig_" rainbow.fasta | mawk '!/^--/' | sed -e 's/dDocent_A_Contig_/dDocent_Contig_/g' > rainbow.asm.fasta
            grep -A1 "dDocent_Contig_" rainbow.fasta | mawk '!/^--/' > rainbow.n.fasta

            sed -e 's/NNNNNNNNNN/   /g' rainbow.asm.fasta | cut -f1 | seqtk seq -F I - > ref.F.fq
            sed -e 's/NNNNNNNNNN/   /g' rainbow.asm.fasta | cut -f2 | seqtk seq -F I - > ref.R.fq

            seqtk seq -r ref.R.fq > ref.RC.fq
            mv ref.RC.fq ref.R.fq
            LENGTH=$(mawk '!/>/' rainbow.fasta | mawk '(NR==1||length<shortest){shortest=length} END {print shortest}')
            LENGTH=$(( $LENGTH * 5 / 4))

            pearRM -f ref.F.fq -r ref.R.fq -o overlap -p 0.001 -j $NUMProc -n $LENGTH &>kopt.log
            rm ref.F.fq ref.R.fq

            if[ -s overlap ]; then
            mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.assembled.fastq > overlap.fasta
            mawk '/>/' overlap.fasta > overlap.loci.names
            mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.forward.fastq > other.F
            mawk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' overlap.unassembled.reverse.fastq > other.R
            paste other.F other.R | mawk '{if ($1 ~ />/) print $1; else print $0}' | sed -e 's/     /NNNNNNNNNN/g' > other.FR

            cat other.FR overlap.fasta rainbow.n.fasta > totalover.fasta
			else 
            cat other.FR overlap.fasta rainbow.n.fasta > totalover.fasta
			fi

            paste <(mawk '{if (NR % 2) print $0}' totalover.fasta) <(mawk '{if (NR % 2 == 0) print $0}' totalover.fasta) | sort -V | sed -e 's/     /\'$'\n/g' > totalover.s.fasta
            cat totalover.s.fasta | tr "\t" "\n" > totalover.fasta
            rm *.F *.R
			cd-hit-est -i totalover.fasta -o reference.fasta.original -M 0 -T 0 -c $k &>cdhit2.log
            sed -e 's/^C/NC/g' -e 's/^A/NA/g' -e 's/^G/NG/g' -e 's/^T/NT/g' -e 's/T$/TN/g' -e 's/A$/AN/g' -e 's/C$/CN/g' -e 's/G$/GN/g' reference.fasta.original > reference.fasta
			samtools faidx reference.fasta &> index.log
            bwa index reference.fasta >> index.log 2>&1
	    	cd ..
		done
    done
done

#Mapping reads
for i in $(seq 1 $K1_MAX); do
    for j in $(seq 1 $K1_MAX); do
        if [ $i -eq 1 && $j -eq 1 ]; then continue; fi
	    for k in $(seq 0.8 0.02 0.98); do
    		cd ref_${i}_${j}_${k}
		  # Adding fastq files
		    cp -s $DATA_DIR/*.R1.fq.gz $DATA_DIR/*.R2.fq.gz .
		  # Adding the dDocent configuration files
		    cp ~/bin/dDocent.config .
		  # Altering the dDocent configuration files
		    awk -v VAR=$i '$0 ~ "K1" {print $0; getline; $0=VAR; print $0; next}{print $0}' dDocent.config
		    awk -v VAR=$j '$0 ~ "K2" {print $0; getline; $0=VAR; print $0; next}{print $0}' dDocent.config
		    awk -v VAR=$k '$0 ~ "Clustering_Similarity" {print $0; getline; $0=VAR; print $0; next}{print $0}' dDocent.config
		    awk -v VAR=$NUMProc '$0 ~ "Processors" {print $0; getline; $0=VAR; print $0; next}{print $0}' dDocent.config
		  # Running dDocent Mapping
		    dDocent dDocent.config
	    	cd ..
		done
    done
done

#Getting statistics
echo -e "File\tK1\tK2\tc-value\tTotal\tQC\tMapped\tPaired" > bwa_mapping_summary.txt

for i in $(seq 1 $K1_MAX); do
    for j in $(seq 1 $K1_MAX); do
        if [ $i -eq 1 && $j -eq 1 ]; then continue; fi
	    for k in $(seq 0.8 0.02 0.98); do
    		cd ref_${i}_${j}_${k}
			    cat namelist | while read l; do
                  #Summarizing Reads mapped
                    TOT=$(zgrep $l.R1.fq.gz| wc -l)
                    QC=$(samtools flagstat -@NUMProc $l.bam | grep total  | grep -E -o '[0-9]+' | head -1)
                    MAP=$(samtools flagstat -@NUMProc $l.bam | grep mapped | grep -E -o '[0-9]+' | head -1)
                    PAIR=$(samtools flagstat -@NUMProc $l.bam | grep properly | grep -E -o '[0-9]+' | head -1)
                    echo -e "$l\t$i\t$j\t$k\t$TOT\t$QC\t$MAP\t$PAIR" >> ../bwa_mapping_summary.txt
			    done
	    	cd ..
		done
    done
done
