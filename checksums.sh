#!/bin/bash


# type in path to dir above folders with bam files
path=$1

cd $path

for d in */
do
	if [ "$d" != "subsets/" ]
	then
		for file in $d*dedup.bam
		do
			md5sum $file
		done
	fi

done
