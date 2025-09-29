#!/bin/bash

# Populations
POPS=("amr" "afr" "eur" "csa" "eas")

# Sample sizes
SIZES=(1000 10000 100000 1000000 10000000)

# Template files (your originals with amr+1000 placeholders)
PHENO_TEMPLATE="sbatch_pheno_template.sh"
GENO_TEMPLATE="sbatch_geno_template.sh"

for pop in "${POPS[@]}"; do
  for size in "${SIZES[@]}"; do
    pheno_out="${pop}_${size}_pheno.sbatch"
    geno_out="${pop}_${size}_geno.sbatch"

    # Generate sbatch scripts
    sed -e "s/amr/${pop}/g" \
        -e "s/1000/${size}/g" \
        "$PHENO_TEMPLATE" > "$pheno_out"

    sed -e "s/amr/${pop}/g" \
        -e "s/1000/${size}/g" \
        "$GENO_TEMPLATE" > "$geno_out"

    # Submit geno script first (array job)
    geno_jobid=$(sbatch "$geno_out" | awk '{print $4}')
    echo "Submitted GENO job $geno_jobid ($geno_out)"

    # Submit pheno script, waiting until the geno job finishes
    pheno_jobid=$(sbatch --dependency=afterok:${geno_jobid} "$pheno_out" | awk '{print $4}')
    echo "Submitted PHENO job $pheno_jobid ($pheno_out), depends on $geno_jobid"
  done
done


