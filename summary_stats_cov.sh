#!/bin/bash

INPUT_FILE=$1

# number of samples
NUM_SAMPLES=$2
((MAX_COL=$NUM_SAMPLES+2))

# calculating mean, standard deviation, and coefficient of variation for each sample (row=sample)
zcat $INPUT_FILE | head -10 | awk '{n++; for(i=3;i<='"$MAX_COL"';i++) sum[i]+=$i;}END{for(i=3;i<='"$MAX_COL"';i++) print sum[i]/n,"\011",sum[i]/n}'

#zcat $INPUT_FILE | head -10 | awk '{n++; for(i=3;i<='"$MAX_COL"';i++) sum[i]+=$i;}END{for(i=3;i<='"$MAX_COL"';i++) print sum[i]/n,"\011",sum[i]/n}'

#zcat $INPUT_FILE | head -10 | awk '{n++; for(i=3;i<='"$MAX_COL"';i++) sum[i]+=$i;}END{for(i=3;i<='"$MAX_COL"';i++) print sum[i]/n,"\011",sum[i]/n}'
