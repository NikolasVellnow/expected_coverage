#!/bin/bash

# type in path to dir above folders with bam files
path=$1

touch test.txt

# get average coverage for "typical" chromosome 3 = NC_031770.1
samtools depth -r NC_031770.1 -f $1 | awk '{n++; for(i=3;i<=5;i++) sum[i]+=$i;}END{for(i=3;i<=5;i++) print "NC_031770.1","\011",sum[i]/n}' >> test.txt

# get average coverage for Z chromosome (ZZ => male, WZ => female)
samtools depth -r NC_031799.1 -f $1 | awk '{n++; for(i=3;i<=5;i++) sum[i]+=$i;}END{for(i=3;i<=5;i++) print "NC_031799.1","\011",sum[i]/n}' >> test.txt

echo "done"

# samtools depth -r NC_031799.1 -f $1 | head -100