#!/bin/bash -l
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=gen_cov_hist_job
#SBATCH --mail-user=nikolas.vellnow@tu-dortmund.de
#SBATCH --mail-type=All



# type in path to python script
PY_SCRIPT=$1

# type in path to data file with coverage values for all positions
DATA_FILE_PATH=$2

# type in output file name
OUT=$3

DATA_FILE="$(basename -- $DATA_FILE_PATH)"

OUT_PATH=/scratch/mnikvell/gen_cov_hist_job${SLURM_JOBID}/

# create directories in scratch-dir
rm -rf /scratch/mnikvell/gen_cov_hist_job${SLURM_JOBID}/
mkdir -p /scratch/mnikvell/gen_cov_hist_job${SLURM_JOBID}/

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

conda activate samtools

python3 "${PY_SCRIPT}" -i "${OUT_PATH}${DATA_FILE}" -o "${OUT}"


# copy outputs back to
cp -a "${OUT_PATH}." /"work/mnikvell/data/mapped_reads/"
rm -rf "${OUT_PATH}"


conda deactivate

