"""Executes the pipeline for evaluating synthetic data quality
"""

include("../utils/reference_data.jl")
include("metrics/eval_aats.jl")
include("metrics/eval_kinship_detail.jl")
include("metrics/eval_ld_corr.jl")
include("metrics/eval_ld_decay.jl")
include("metrics/eval_maf.jl")
include("metrics/eval_pca.jl")
include("metrics/eval_gwas.jl")
 

function run_kinship_evaluation(ibsfile_real, ibsfile_synt, ibsfile_cross)
    run_kinship(ibsfile_real, ibsfile_synt, ibsfile_cross)
end


function run_aats_evaluation(ibsfile_cross)
    run_aats(ibsfile_cross)
end


function run_ld_corr_evaluation(real_data_prefix, synt_data_prefix, eval_dir, plink_path)
    run_ld_corr(real_data_prefix, synt_data_prefix, eval_dir, plink_path)
end


function run_ld_decay_evaluation(synt_data_prefix, real_data_prefix, plink_path, mapthin_path, eval_dir, bp_to_cm_map, chromosome)
    run_ld_decay(synt_data_prefix, real_data_prefix, plink_path, mapthin_path, eval_dir, bp_to_cm_map, chromosome)
end


function run_maf_evaluation(real_maf_file, synt_maf_file)
    run_maf(real_maf_file, synt_maf_file)
end


function run_pca_evaluation(real_data_pca_prefix, synt_data_pca_prefix, pcaproj_file, real_data_prefix, synt_data_prefix, eval_dir)
    real_data_pop = string(real_data_prefix, ".sample")
    synt_data_pop = string(synt_data_prefix, ".sample")
    run_pca(real_data_pca_prefix, synt_data_pca_prefix, pcaproj_file, real_data_pop, synt_data_pop, eval_dir)
end


function run_gwas_evaluation(ntraits, plink2, covar, synt_data_prefix, eval_dir, chromosome)
    run_gwas(ntraits, plink2, covar, synt_data_prefix, eval_dir, chromosome)
end


function run_mia_evaluation(orig_prefix, gen_prefix, out_dir, balancing, pca_components, train_fraction, shuffle, random_state)
    @info "Running MIA (Python)"
    py_script = joinpath(@__DIR__, "metrics", "eval_mia.py")
    println(`python3 $py_script --orig_prefix $orig_prefix --gen_prefix $gen_prefix --out_dir $out_dir --balancing $balancing --pca_components $pca_components --train_fraction $train_fraction --shuffle $shuffle --random_state $random_state`)
    
    run(`python3 $py_script --orig_prefix $orig_prefix --gen_prefix $gen_prefix --out_dir $out_dir --balancing $balancing --pca_components $pca_components --train_fraction $train_fraction --shuffle $shuffle --random_state $random_state`)
end


"""Computations for MAF using PLINK
"""
function run_maf_tools(plink, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for MAF"
    reffile_out = @sprintf("%s.ref.maf", outdir)
    synfile_out = @sprintf("%s.syn.maf", outdir)
    run(`$plink --bfile $reffile_prefix --freq --out $reffile_out`)
    run(`$plink --bfile $synfile_prefix --freq --out $synfile_out`)
    real_maffile = @sprintf("%s.frq", reffile_out)
    syn_maffile = @sprintf("%s.frq", synfile_out)
    return real_maffile, syn_maffile
end


"""Computations for PCA using PLINK
"""
function run_pca_tools(plink2, king, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for PCA"
    reffile_out = @sprintf("%s.ref.pca", outdir)
    synfile_out = @sprintf("%s.syn.pca", outdir)
    run(`$plink2 --bfile $reffile_prefix --freq counts --pca approx allele-wts --out $reffile_out`)
    run(`$plink2 --bfile $synfile_prefix --freq counts --pca approx allele-wts --out $synfile_out`)

    # For PCA projection
    projfile_prefix = @sprintf("%s_proj", outdir)
    projfile_out = @sprintf("%spc.txt", projfile_prefix)
    run(`$king -b $reffile_prefix.bed,$synfile_prefix.bed --pca --projection --prefix $projfile_prefix`)

    return reffile_out, synfile_out, projfile_out
end


"""Computations for relatedness using KING
"""
function run_relatedness_tools(king, reffile_prefix, synfile_prefix, outdir)
    @info "Running external tools for relatedness"
    reffile_bed = @sprintf("%s.bed", reffile_prefix)
    synfile_bed = @sprintf("%s.bed", synfile_prefix)
    reffile_out = @sprintf("%s.ref.king", outdir)
    synfile_out = @sprintf("%s.syn.king", outdir)
    crossfile_out = @sprintf("%s.cross.king", outdir)
    run(`$king -b $reffile_bed --ibs --prefix $reffile_out`)
    run(`$king -b $synfile_bed --ibs --prefix $synfile_out`)
    run(`$king -b $reffile_bed,$synfile_bed --ibs --prefix $crossfile_out`)
    real_ibsfile = @sprintf("%s.ibs0", reffile_out)
    syn_ibsfile = @sprintf("%s.ibs0", synfile_out)
    cross_ibsfile = @sprintf("%s.ibs0", crossfile_out)
    return real_ibsfile, syn_ibsfile, cross_ibsfile
end


"""Runs computations using external tools such as PLINK and KING, based on the selection of evaluation metrics
"""
function run_external_tools(metrics, reffile_prefix, synfile_prefix, filepaths)
    external_files = Dict()
    if metrics["maf"]
        external_files["real_maffile"], external_files["syn_maffile"] = run_maf_tools(filepaths.plink, reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["pca"] || metrics["gwas"]
        external_files["real_pcafile"], external_files["syn_pcafile"], external_files["pcaproj_file"] = run_pca_tools(
            filepaths.plink2, 
            filepaths.king,
            reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["kinship"] || metrics["aats"]
        external_files["real_ibsfile"], external_files["syn_ibsfile"], external_files["cross_ibsfile"] = run_relatedness_tools(filepaths.king, reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    return external_files
end


"""Executes evaluation for the metrics specified in the configuration file
"""
function run_pipeline(options, chromosome, superpopulation, metrics, filepaths, genomic_metadata, reffile_prefix, nsamples_ref, synfile_prefix, external_files)
    if metrics["aats"]
        run_aats_evaluation(external_files["cross_ibsfile"])
    end
    if metrics["kinship"]
        run_kinship_evaluation(external_files["real_ibsfile"], external_files["syn_ibsfile"], external_files["cross_ibsfile"])
    end
    if metrics["ld_corr"]
        run_ld_corr_evaluation(reffile_prefix, synfile_prefix, filepaths.evaluation_output, filepaths.plink)
    end
    if metrics["ld_decay"]
        bp_to_cm_map = create_bp_cm_ref(filepaths.genetic_distfile)
        run_ld_decay_evaluation(synfile_prefix, reffile_prefix, filepaths.plink, filepaths.mapthin, filepaths.evaluation_output, bp_to_cm_map, chromosome)
    end
    if metrics["maf"]
        run_maf_evaluation(external_files["real_maffile"],  external_files["syn_maffile"])
    end
    if metrics["pca"]
        run_pca_evaluation(
            external_files["real_pcafile"], external_files["syn_pcafile"], external_files["pcaproj_file"],
            reffile_prefix, synfile_prefix, filepaths.evaluation_output)
    end
    if metrics["mia"]
        mia_cfg = haskey(options["evaluation"], "mia") ? options["evaluation"]["mia"] : Dict{String,Any}()
        balancing = haskey(mia_cfg, "balancing") ? mia_cfg["balancing"] : "undersample"
        pca_components = haskey(mia_cfg, "pca_components") ? mia_cfg["pca_components"] : 20
        train_fraction = haskey(mia_cfg, "train_fraction") ? mia_cfg["train_fraction"] : 0.7
        shuffle = haskey(mia_cfg, "shuffle") ? mia_cfg["shuffle"] : true
        out_dir = filepaths.evaluation_output
        random_state = options["global_parameters"]["random_seed"]
        run_mia_evaluation(reffile_prefix, synfile_prefix, out_dir, balancing, pca_components, train_fraction, shuffle, random_state)
    end
    if metrics["gwas"]
        run_gwas_evaluation(options["phenotype_data"]["nTrait"], filepaths.plink2, @sprintf("%s.eigenvec", external_files["syn_pcafile"]), synfile_prefix, filepaths.evaluation_output, chromosome)
    end
    

end


"""Entry point to running the evaluation pipeline for genotype data

Note that the evaluation pipeline assumes that the synthetic data you want 
to evaluate has already been geneerated, using the setup specified in the
configuration file. It is therefore recommended to run the pipeline with the 
--genotype and --evaluation flags together, so that the program generates 
the data and then immediately evaluates it using the correct settings.
"""
function run_evaluation(options)
    chromosome = parse_chromosome(options)
    superpopulation = parse_superpopulation(options)

    metrics = options["evaluation"]["metrics"]
    gwas = metrics["gwas"] 
    if chromosome == "all"      
        covar_paths = String[]
        for chromosome_i in 1:22
            # Prepare per-chromosome filepaths and metadata
            filepaths = parse_filepaths(options, chromosome_i, superpopulation)
            genomic_metadata = parse_genomic_metadata(options, superpopulation, filepaths)

            reffile_prefix, nsamples_ref = create_reference_dataset(filepaths.vcf_input_processed, filepaths.popfile_processed, genomic_metadata.population_weights, filepaths.plink, filepaths.reference_dir, chromosome_i)
            synfile_prefix =  filepaths.synthetic_data_prefix

            # Ensure PCA is computed per chromosome if GWAS is requested
            local_metrics = copy(metrics)
            local_metrics["gwas"] = false
            if gwas
                local_metrics["pca"] = true
            end

            external_files = run_external_tools(local_metrics, reffile_prefix, synfile_prefix, filepaths)
            run_pipeline(options, chromosome_i, superpopulation, local_metrics, filepaths, genomic_metadata, reffile_prefix, nsamples_ref, synfile_prefix, external_files)

            # Collect per-chromosome covariate path (PCA eigenvectors)
            if gwas
                push!(covar_paths, @sprintf("%s.eigenvec", external_files["syn_pcafile"]))
            end
        end
        # If GWAS was originally requested, run it across all chromosomes now
        if gwas
            filepaths_all = parse_filepaths(options, 1, superpopulation)
            synfile_prefix_all = filepaths_all.synthetic_data_prefix
            run_gwas_evaluation(options["phenotype_data"]["nTrait"], filepaths_all.plink2, covar_paths, synfile_prefix_all, filepaths_all.evaluation_output, "all")
        end
    else
        filepaths = parse_filepaths(options, chromosome, superpopulation)
        genomic_metadata = parse_genomic_metadata(options, superpopulation, filepaths)

        reffile_prefix, nsamples_ref = create_reference_dataset(filepaths.vcf_input_processed, filepaths.popfile_processed, genomic_metadata.population_weights, filepaths.plink, filepaths.reference_dir, chromosome)
        synfile_prefix =  filepaths.synthetic_data_prefix

        external_files = run_external_tools(metrics, reffile_prefix, synfile_prefix, filepaths)
        run_pipeline(options, chromosome, superpopulation, metrics, filepaths, genomic_metadata, reffile_prefix, nsamples_ref, synfile_prefix, external_files)
    end
end
