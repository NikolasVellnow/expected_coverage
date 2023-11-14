#!/bin/bash -l
#SBATCH --partition=med
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=07:58:00 
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=samtools_depth_job
#SBATCH --mail-user=nikolas.vellnow@tu-dortmund.de
#SBATCH --mail-type=All



# type in path to text file with list of samples
SAMPLE_LIST=$1

# type in output file name
OUT=$2


# number of threads used in samtools
NUM_THREADS=32


OUT_PATH=/scratch/mnikvell/samtools_depth_job${SLURM_JOBID}/

# create directories in scratch-dir
rm -rf /scratch/mnikvell/samtools_depth_job${SLURM_JOBID}/
mkdir -p /scratch/mnikvell/samtools_depth_job${SLURM_JOBID}/


T0=$(date +%T)
echo "Start:"
echo $T0

conda activate samtools

touch $OUT_PATH$OUT

# write sequencing depth at every (?) position to text file
samtools depth -H -@ $NUM_THREADS -f $SAMPLE_LIST > $OUT_PATH$OUT


# reformat the header of the file
HEADER_START="CHROM\tPOS\t"
HEAD=$(head -1 $OUT_PATH$OUT)

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
sed -i "1s/.*/$NEW_HEADER/" $OUT_PATH$OUT

gzip $OUT_PATH$OUT

# copy outputs back to
cp -a "${OUT_PATH}." /"work/mnikvell/data/mapped_reads/"
rm -rf /scratch/mnikvell/samtools_depth_job${SLURM_JOBID}/


T1=$(date +%T)
echo "Finished:"
echo $T1


conda deactivate

