#!/bin/bash

# type in path to dir above folders with bam files
path=$1

#conda activate samtools

cd $path

for d in */
do
	if [ "$d" != "subsets/" ]
	then
		for file in $d*dedup.bam
		do
			FILENAME="$file"
			FILENAME=${FILENAME%.bam*}
			echo $FILENAME
			# create random subset with size of 0.0001 of original
			samtools view -bo "$FILENAME"_subset.bam -@8 -s 123.0001 "$file"
		done
	fi

done
