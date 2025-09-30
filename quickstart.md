# How to use this tool?

For a more comprehensive guide, see the **README** file.

---

## Step 1: Pull the container

Pull the container (Docker or Singularity) using this command:

```bash
singularity pull docker://benonsen/synthetic-genetic-data
```

If some singularity out of memory issues appear, set the cache and tmp dir before pulling the container.
```bash
export SINGULARITY_TMPDIR=<path-to-tmp>/tmp
export SINGULARITY_CACHEDIR=<path-to-tmp>/tmp/cache
```

---

## Step 2: Initialization

Run `init` and install all the needed packages using this command. Also, bind the data folder to make it persistent; move the `config.yaml` file into `data/config.yaml` beforehand.

If you encounter any Julia-related errors, add this after the bind parameter to the command:

```
--env JULIA_DEPOT_PATH=/data/.julia
```

Command:

```bash
singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif init
```

---

## Step 3: Fetch reference data

Fetch all the reference data. This will download the reference dataset:

```bash
singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif fetch
```

---

## Step 4: Generate data

To generate data, run:

```bash
singularity exec --bind data/:/data/ containers/synthetic-genetic-data_latest.sif generate_geno <number-of-threads> data/config.yaml
```

In practice, we do not use this single command; instead, we use a script that automatically submits all the jobs with the correct dependencies. I will move the files to the correct directory.

---

## File structure

**What files do we have?**

* `/data/inputs/processed/1KG+HGDP/`:
  Reference dataset from 1KG. This path contains 5 folders: `EUR`, `CSA`, `EAS`, `AMR`, `AFR`.
  Each folder contains the `CausalList` file for a trait.

  * Trait 1 corresponds to `CausalList1`
  * Trait 2 corresponds to `CausalList2`
  * etc.
    First, specify the SNP and then the effect size.
    I provided a script (see create_causalList.py) to convert GWAS Catalog associations into a `CausalList` file automatically. As input provide the downloaded file from GWASCatalog and the output file.

* `/data/`:
  Also contains 5 folders (`eur`, `csa`, `eas`, `amr`, `afr`) with the corresponding `config.yaml` files.

  * For each sample size, a new config file was created.
  * Make sure the path for the `CausalList` points to the correct folder and the correct population is specified. (I already checked it; for this case, it is fine.)
  * The `chromosome` parameter can be `"all"` or an integer between `1` and `22`. Since we are operating in a cluster, we set it to `${chr}` and replace it with a number between `1` and `22` to generate in parallel.

After generating the genotype, we generate the phenotype, validate, and convert (all sequentially). The script automatically replaces the `${chr}` placeholder with `"all"`.

---

## Step 5: Submit jobs

After running the `init` and `fetch` commands, just run the `submit_jobs.sh` script to submit all the jobs.

**Important:**
This script uses two underlying template scripts and replaces the sample size and population identifier. Therefore, please check the paths in those scripts before running `submit_jobs.sh`.

The template scripts are:

* `sbatch_geno_template.sh`
* `sbatch_pheno_template.sh` (this one generates phenotypes, validates, and converts).

---

## Output

* Results can be found here:

  ```
  /data/outputs/test/{population}/{sample size}
  ```

* VCF files are in the same path but in a dedicated `VCF` folder.

* Evaluation results can be found here:

  ```
  /data/outputs/test/{population}/{sample size}/evaluation
  ```

This folder contains multiple `sumstats` files: one for each trait and each chromosome.

* For summary statistics of a trait across all chromosomes, see:

  ```
  test_chr-1.traitX.all.sumstat
  ```

* This path also includes all **Manhattan plots**, **QQ plots**, and intermediate results (`pca`, `glm.linear` files, etc.).

