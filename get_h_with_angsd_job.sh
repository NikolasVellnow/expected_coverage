#!/bin/bash -l
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:20:00 
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=gen_get_h_with_angsd_job
#SBATCH --mail-user=nikolas.vellnow@tu-dortmund.de
#SBATCH --mail-type=All




NUM_THREADS=20

# type in path to bam file
DATA_FILE_PATH=$1

# type in output file name
OUT=$2

FILE_PATH="$(dirname -- $DATA_FILE_PATH)"
DATA_FILE="$(basename -- $DATA_FILE_PATH)"

GENOME_PATH=/home/mnikvell/Desktop/work/data/genomes/refseq/vertebrate_other/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna

OUT_PATH=/scratch/mnikvell/get_h_with_angsd_job${SLURM_JOBID}/

# create directories in scratch-dir
rm -rf /scratch/mnikvell/get_h_with_angsd_job${SLURM_JOBID}/
mkdir -p /scratch/mnikvell/get_h_with_angsd_job${SLURM_JOBID}/


T0=$(date +%T)
echo "Start copying data file:"
echo $T0

# copy data file
cp -a "${DATA_FILE_PATH}" "${OUT_PATH}"

T1=$(date +%T)
echo "Finished copying data file:"
echo $T1

# check content of scratch dir
ls "${OUT_PATH}"


conda activate angsd

cd "${OUT_PATH}"

angsd -i "${DATA_FILE}" -doSaf 1 -anc "${GENOME_PATH}" -r NC_031779.1:5000-7000 -P "${NUM_THREADS}"

realSFS angsdput.saf.idx >est.ml

rm "${DATA_FILE}"

# copy outputs back to
cp -a "${OUT_PATH}." "${FILE_PATH}"
rm -rf "${OUT_PATH}"


conda deactivate

