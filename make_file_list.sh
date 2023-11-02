#!/bin/bash

# type in path to dir above folders with bam files
path=$1

# output file
out=$2

touch $out

cd $path

for d in */
do
	if [ "$d" != "subsets/" ]
	then
		for file in $d*dedup.bam
		do
			FILENAME="$file"
			echo $path$FILENAME >> $out
		done
	fi

done
