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

touch $OUT

# get average coverage for "typical" chromosome 3 = NC_031770.1
samtools depth -r NC_031770.1 -@ $NUM_THREADS -f $SAMPLE_LIST | awk '{n++; for(i=3;i<='"$MAX_COL"';i++) sum[i]+=$i;}END{for(i=3;i<='"$MAX_COL"';i++) print "NC_031770.1","\011",sum[i]/n}' > $OUT

# get average coverage for Z chromosome (ZZ => male, WZ => female)
samtools depth -r NC_031799.1 -@ $NUM_THREADS -f $SAMPLE_LIST | awk '{n++; for(i=3;i<='"$MAX_COL"';i++) sum[i]+=$i;}END{for(i=3;i<='"$MAX_COL"';i++) print "NC_031799.1","\011",sum[i]/n}' >> $OUT

echo "done"

