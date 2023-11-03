#!/bin/bash

# type in path to text file with list of samples
SAMPLE_LIST=$1

# type in output file name
OUT=$2

# number of samples
NUM_SAMPLES=$3
((MAX_COL=$NUM_SAMPLES+2))


# number of threads used in samtools
NUM_THREADS=$4

T0=$(date +%T)
echo $T0

touch $OUT

# write sequencing depth at every (?) position to text file
samtools depth -H -@ $NUM_THREADS -f $SAMPLE_LIST > $OUT

T1=$(date +%T)
echo $T1

echo "done"

