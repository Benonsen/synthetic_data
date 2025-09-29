#!/bin/bash
#SBATCH --job-name=amr_1000_pheno
#SBATCH --output=R-%x.%j.out
#SBATCH --time=2-24:00:00
#SBATCH --mem=32G
#SBATCH --partition="batch-mc2"

CONFIG=data/amr/config_1000 # prefix for config file


# replace the chromosome wildcard with all
cp ${CONFIG}.yaml ${CONFIG}_pheno$n.yaml
sed -i 's/${chr}'"/all/g" ${CONFIG}_pheno$n.yaml

# generate phenotype data
singularity exec --bind /nvme01/compgenomics/bsteinegger/syntetic_data/data:/data --env JULIA_DEPOT_PATH=/data/.julia synthetic-genetic-data_latest.sif generate_pheno ${CONFIG}_pheno.yaml
singularity exec --bind /nvme01/compgenomics/bsteinegger/syntetic_data/data:/data --env JULIA_DEPOT_PATH=/data/.julia synthetic-genetic-data_latest.sif validate  ${CONFIG}_pheno.yaml
singularity exec --bind /nvme01/compgenomics/bsteinegger/syntetic_data/data:/data --env JULIA_DEPOT_PATH=/data/.julia synthetic-genetic-data_latest.sif convert  ${CONFIG}_pheno.yaml