#!/bin/bash

#This script is to provide a quick way of changing the population structure in a genepop file using a popmap file

if [[ $# -ne 3 ]]; then
echo "Usage: genepop_newpop.sh file.gen popmap output.gen"
exit 1
fi

INPUT=$1
POP=$2
OUTPUT=$3

awk '$0 ~ "POP"{exit 1}{print $0}' $INPUT > $OUTPUT
awk '{FS=" |\t"} {print $2}' $POP | sort | uniq > tmp
sed -i '/^\s*$/d' tmp
cat tmp | while read i; do
echo "Processing $i";
grep -w $i $POP | awk '{FS=" |\t"} {print $1}' > $i.names;
grep -wf $i.names $INPUT > $i.gen;
awk 'BEGIN {print "POP"}{print $0}' $i.gen > tmp2;
mv tmp2 $i.gen;
cat $OUTPUT $i.gen > tmp2;
mv tmp2 $OUTPUT;
done

