#!/bin/bash

# type in path to text file with list of samples
SAMPLE_LIST=$1

# type in output file name
OUT=$2


# number of threads used in samtools
NUM_THREADS=$3

T0=$(date +%T)
echo "Start:"
echo $T0

touch $OUT

# write sequencing depth at every (?) position to text file
samtools depth -H -@ $NUM_THREADS -f $SAMPLE_LIST > $OUT


# reformat the header of the file
HEADER_START="CHROM\tPOS\t"
HEAD=$(head -1 $OUT)

new_str=""
length=${#HEAD}
for ((i = 0; i < length; i++)); do
	pattern="${HEAD:i:7}"
	if [ "$pattern" == "_reads/" ]
	then
    		sample="${HEAD:i+7:3}"
		if [ ${sample:2:1} == "_" ]
		then
			sample="${sample:0:2}"
		fi
		new_string+=$sample
		new_string+='\t'
	fi
done

NEW_HEADER_TAIL=$new_string
NEW_HEADER="$HEADER_START$NEW_HEADER_TAIL"

# replace old header with new one
sed -i "1s/.*/$NEW_HEADER/" $OUT

T1=$(date +%T)
echo "Finished:"
echo $T1


