#!/bin/bash

# type in path to dir above folders with bam files
path=$1

#conda activate samtools

cd $path

# loops through directories
for d in */
do
	if [ "$d" != "subsets/" ]
	then
		# finds subset file in directory
		for file in $d*_subset.bam
		do
			echo $file
			# create random subset with size of 0.0001 of original
			samtools index -@ 8 $file
			#samtools view -bo "$FILENAME"_subset.bam -@8 -s 123.0001 "$file"
		done
	fi

done


