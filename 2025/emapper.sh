#!/bin/bash
#
#SBATCH --job-name=emapper
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=02:00:00

# Exit the slurm script if a command fails
set -e

# Activate the eggNOG mapper environment (replace version according to your needs)
module load conda
conda activate eggnog-mapper-2.1.12

QUERY=$1

# Cache the eggNOG database and store the path of the cache in the environment variable EGGNOG_DATA_DIR
export EGGNOG_DATA_DIR=$(lisc_localcache /lisc/scratch/mirror/eggnog-mapper/2.1.12)

# Run eggNOG mapper and always specify a local temporary folder for the tmp_dir parameter.
emapper.py -i "$QUERY" -o query.emapper.txt --cpu "$SLURM_CPUS_PER_TASK" --temp_dir="$TMPDIR"

# If we reached this point, the database search succeeded. We clean up resources.
rm -rf "$TMPDIR"