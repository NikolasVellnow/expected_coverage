#!/bin/bash -l
#SBATCH --partition=med
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=07:59:00 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=gen_cov_hist_job
#SBATCH --mail-user=nikolas.vellnow@tu-dortmund.de
#SBATCH --mail-type=All



# type in path to python script
PY_SCRIPT=$1

# type in path to data file with coverage values for all positions
DATA_FILE=$2

# type in output file name
OUT=$3

OUT_PATH=/scratch/mnikvell/gen_cov_hist_job${SLURM_JOBID}/

# create directories in scratch-dir
rm -rf /scratch/mnikvell/gen_cov_hist_job${SLURM_JOBID}/
mkdir -p /scratch/mnikvell/gen_cov_hist_job${SLURM_JOBID}/


conda activate samtools

python3 "${PY_SCRIPT}" -i "${DATA_FILE}" -o "${OUT}"


# copy outputs back to
cp -a "${OUT_PATH}." /"work/mnikvell/data/mapped_reads/"
rm -rf /scratch/mnikvell/gen_cov_hist_job${SLURM_JOBID}/


conda deactivate

