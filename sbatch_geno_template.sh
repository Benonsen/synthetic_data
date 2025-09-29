#!/bin/bash

#SBATCH --array=1-22
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=8
#SBATCH --time 24:00:00
#SBATCH --partition="batch-mc2"
#SBATCH --job-name=amr_1000_geno
#SBATCH --output=R-%x.%j.out

n=$SLURM_ARRAY_TASK_ID

CONFIG=data/amr/config_1000 # prefix for config file

# generate a config for each chromosome
cp ${CONFIG}.yaml ${CONFIG}$n.yaml
sed -i 's/${chr}'"/$n/g" ${CONFIG}$n.yaml

# generate data for each chromosome
singularity exec --bind /nvme01/compgenomics/bsteinegger/syntetic_data/data:/data --env JULIA_DEPOT_PATH=/data/.julia  synthetic-genetic-data_latest.sif generate_geno 8 ${CONFIG}$n.yaml