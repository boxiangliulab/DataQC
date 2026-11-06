def helpMessage() {
    log.info"""
    =======================================================
     DataQC v${workflow.manifest.version}
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run DataQC.nf \
        --bfile EstBB_HT12v3\
        --expfile EstBB_HT12v3_exp.txt\
        --gte gte_EstBB_HT12v3.txt\
        --exp_platform HT12v3\
        --cohort_name EstBB_HT12v3\
        --genome_build GRCh37
        --outdir EstBB_HT12v3_PreImputationQCd\
        -profile slurm\
        -resume

    Mandatory arguments:
      --cohort_name                 Name of the cohort.
      --genome_build                Genome build of the cohort. Either hg18, GRCh36, hg19, GRCh37, hg38 or GRCh38.
      --bfile                       Path to unimputed genotype files in plink bed/bim/fam format (without extensions bed/bim/fam).
      --vcf                         Path to a vcf file.
      --fam                         Path to a plink fam file. This is especially helpful for sex annotation of samples in VCF files.
      --expfile                     Path to the un-preprocessed gene expression matrix (genes/probes in the rows, samples in the columns). Can be from RNA-seq experiment or from array. NB! For Affymetrix arrays (AffyU219, AffyExon) we assume that standard preprocessing and normalisation is already done.
      --gte                         Genotype-to-expression linking file. Tab-delimited, no header. First column: sample ID for genotype data. Second column: corresponding sample ID for gene expression data. Can be used to filter samples from the analysis.
      --exp_platform                Indicator indicating the gene expression platform. HT12v3, HT12v4, HuRef8, RNAseq, AffyU219, AffyHumanExon, RNAseq_HGNC.
      --outdir                      Path to the output directory.
      --GenOutThresh                "Outlierness" score threshold for excluding ethnic outliers. Defaults to 0.4 but it should be adjusted according to visual inspection.
      --GenSdThresh                 Threshold for declaring samples outliers based on genetic PC1 and PC2. Defaults to 3 SD from the mean of PC1 and PC2 but should be adjusted according to visual inspection.
      --ExpSdThresh                 Standard deviation threshold for excluding gene expression outliers. By default, samples away by 3 SDs from the mean of PC1 are removed.
      --ContaminationSlope          Slope (in degrees) that is used to draw a line separating females and males based on sex-specific genes. Must be between 0 (a completely horizontal line), and 90 (a completely vertical line). 45 by default
      --ContaminationArea           Area that marks likely contaminated samples based on sex chromosome gene expression. Must be an angle between 0 and 90. The angle represents the total area around the y = x function (45 degrees, or any other slope depending on the --ContaminationSlope setting)
      --gen_qc_steps                Either generic, array-based, QC or including also WGS specific QC (only valid with VCF datasets). 'Array' (default) or 'WGS' (Generic + WGS qc).

    Optional arguments
      --InclusionList               File with sample IDs to restrict to the analysis. Useful for keeping in the inclusion list of the samples. By default, all samples are kept.
      --ExclusionList               File with sample IDs to remove from the analysis. Useful for removing the ancestry outliers or restricting the genotype data to one superpopulation. Samples are also removed from the inclusion list. By default, all samples are kept.
      --AdditionalCovariates        File with additional cohort-specific covariates. First column name SampleID is the sample ID. Following columns are named by covariates.  Categorical covariates need to be text-based (e.g. batch1, batch2, etc). 
      --preselected_sex_check_vars  Path to a plink ranges file that defines which variants to use for the check-sex command. Use this when the automatic selection does not yield satisfactory results.
      --plink_executable            Path to plink executable. By default this is automatically downloaded from internet. Use this setting when you have to work offline.
      --plink2_executable           Path to plink2 executable. By default this is automatically downloaded from internet. Use this setting when you have to work offline.
      --reference_1000g_folder      Path to 1000g reference folder. By default this is automatically downloaded from internet. Use this setting when you have to work offline.
      --chain_path                  Path to folder containing hg19ToHg38 and hg38ToHg19 chain files. By default these are automatically downloaded from internet. Use this setting when you have to work offline and your build is hg38.

    """.stripIndent()
}

// Define location of Report_template.Rmd
params.report_template = "$baseDir/bin/Report_template.Rmd"

// Define set of accepted genome builds:
def genome_builds_accepted = ['hg18', 'GRCh36', 'hg19', 'GRCh37', 'hg38', 'GRCh38']
def genotyping_platforms_accepted = ['Array', 'WGS']

// Define input channels

// Take care of the following options:
// --bfile <single_plink_dataset>
// --vcf <single_vcf_file>
// --vcf <vcf_files_through_globbing>

params.vcf = ''
params.bfile = ''
params.fam = ''

params.plink_executable = ''
params.plink2_executable = ''
params.reference_1000g_folder = ''
params.chain_path = ''

if (params.vcf != '') {

  Channel
      .fromPath(params.vcf, checkIfExists: true).view()
      .ifEmpty { exit 1, "Input vcf files not found!" }
      .into { vcffile_ch; vcffile_ch2; vcf_contig_counter }

  Channel.empty()
    .set { bfile_ch }

} else {

  Channel.empty()
    .into { vcffile_ch; vcffile_ch2; vcf_contig_counter }

  Channel
    .from(params.bfile)
    .ifEmpty { exit 1, "Input plink prefix not found!" }
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { bfile_ch }

}

vcf_contig_counter.count().view().set { vcf_contig_count }

if (params.fam != '') {

  Channel
    .fromPath(params.fam, checkIfExists: true)
    .set { fam_annot_ch }

} else {

  Channel.empty()
    .set { fam_annot_ch }

}

// Splicing/expression input channels (combined_matrix_ch, expression_matrix_ch, junction_site_ch,
// ir_site_ch) are created further below, after the params block, from a single --splice_pheno_dir.

Channel
    .from(params.gte)
    .map { study -> [file("${study}")]}
    .ifEmpty { exit 1, "Input GTE file not found!" }
    .into { gte_ch_gen; gte_ch_exp; gte_ch_exp2 }

Channel
    .fromPath(params.report_template)
    .ifEmpty { exit 1, "Input report not found!" }
    .set { report_ch }

if (params.plink_executable != '') {
  Channel
    .fromPath(params.plink_executable)
    .ifEmpty('EMPTY')
    .set { plink_executable_ch }
} else {
  Channel.empty().set {plink_executable_ch}
}

if (params.plink2_executable != '') {
  Channel
    .fromPath(params.plink2_executable)
    .ifEmpty('EMPTY')
    .set { plink2_executable_ch }
} else {
  Channel.empty().set {plink2_executable_ch}
}
if (params.reference_1000g_folder != '') {
  Channel
    .fromPath(params.reference_1000g_folder)
    .ifEmpty('EMPTY')
    .set { reference_1000g_ch }
} else {
  Channel.empty().set {reference_1000g_ch}
}
if (params.chain_path != '') {
  Channel
    .fromPath(params.chain_path)
    .ifEmpty('EMPTY')
    .set { chain_path_ch }
} else {
  Channel.empty().set {chain_path_ch}
}


params.GenOutThresh = 0.4
params.GenSdThresh = 3
params.ExpSdThresh = 4
params.ContaminationSlope = 45
params.ContaminationArea = 30
params.exp_platform = ''
params.cohort_name = ''

// Splicing QC + covariate params (the splicing analog of the expression-QC inputs)
params.combined_matrix = ''
params.sd_threshold = 4
params.int_transform = true
params.num_gen_pcs = 10
params.num_splicing_pcs = 50
params.n_top_features = 20000      // TOP-FEATURES variant: splicing PCA/MDS use the N most-variable junctions
params.max_pca_missingness = 0.10  // PRE-FILTER: splicing-PCA basis keeps features with per-site missingness < this
params.min_pca_meancount   = 5     // PRE-FILTER: splicing-PCA basis keeps features with mean raw read depth/donor >= this
params.num_expression_pcs = 100   // eQTLGen uses up to 100 expression PCs in the covariate file
// Optional cohort-specific covariates folded into CovariatePCs.Splicing/Expression.txt (MakeCovariateFile.R behaviour).
// Sentinel 1000G_pops.txt = none supplied (matches the reference default).
params.AdditionalCovariates = "$baseDir/data/1000G_pops.txt"
params.outdir = ''
params.genome_build = 'hg19'

// By default define random non-colliding file names in data folder. If default, these are ignored by corresponding script.
params.InclusionList = "$baseDir/data/EmpiricalProbeMatching_AffyHumanExon.txt"
params.ExclusionList = "$baseDir/data/EmpiricalProbeMatching_AffyU219.txt"

// ---- Stage-0 splicing + expression inputs ----
// A single --splice_pheno_dir (the 0_SplicingPhenotype/output folder) derives everything by cohort_name:
//   combined_matrix/<cohort>_combined_splice_site.txt.gz   (QC'd splice input)
//   expression_matrix/<cohort>_gene_counts.txt.gz          (raw gene counts, for expression QC + 4_ permuted/LD path)
//   03_phenotype/<cohort>.site  +  04_IR/<cohort>_IR.site  (raw counts, for the splice read-count/missingness summary)
// Any of these can be overridden explicitly via --combined_matrix / --expression_matrix / --junction_site / --ir_site.
params.splice_pheno_dir    = ''
params.expression_matrix   = ''
params.junction_site       = ''
params.ir_site             = ''
params.expr_remove_outliers = true    // ACTIVE expression QC: remove expr PCA/MDS/sex/contamination outliers, then
                                       // re-align the splicing set to the intersection (geno n splice n expr) via FINALIZE_SPLICING

def _spd = params.splice_pheno_dir
def combined_matrix_path   = (params.combined_matrix   && params.combined_matrix   != '') ? params.combined_matrix   : "${_spd}/combined_matrix/${params.cohort_name}_combined_splice_site.txt.gz"
def expression_matrix_path = (params.expression_matrix && params.expression_matrix != '') ? params.expression_matrix : "${_spd}/expression_matrix/${params.cohort_name}_gene_counts.txt.gz"
def junction_site_path     = (params.junction_site     && params.junction_site     != '') ? params.junction_site     : "${_spd}/03_phenotype/${params.cohort_name}.site"
def ir_site_path           = (params.ir_site           && params.ir_site           != '') ? params.ir_site           : "${_spd}/04_IR/${params.cohort_name}_IR.site"

Channel.fromPath(combined_matrix_path)
    .ifEmpty { exit 1, "Combined splice matrix not found: ${combined_matrix_path}" }
    .set { combined_matrix_ch }
Channel.fromPath(expression_matrix_path)
    .ifEmpty { exit 1, "Expression matrix not found: ${expression_matrix_path}" }
    .set { expression_matrix_ch }
// .site files are large and optional: if absent, fall back to a sentinel so the splice QC still runs;
// compute_site_readcount_summary.py simply skips a file that has no include:total fields.
Channel.fromPath(junction_site_path).ifEmpty { file("${baseDir}/data/1000G_pops.txt") }.set { junction_site_ch }
Channel.fromPath(ir_site_path).ifEmpty { file("${baseDir}/data/1000G_pops.txt") }.set { ir_site_ch }
params.AdditionalCovariates = "$baseDir/data/1000G_pops.txt"

params.preselected_sex_check_vars = ''
params.gen_qc_steps = 'Array'

if ((params.gen_qc_steps in genotyping_platforms_accepted) == false) {
  exit 1, "[Pipeline error] Genotype QC steps $params.gen_qc_steps not one of: $genotyping_platforms_accepted \n"
}

if ((params.genome_build in genome_builds_accepted) == false) {
  exit 1, "[Pipeline error] Genome build $params.genome_build not in accepted genome builds: $genome_builds_accepted \n"
}

// Header log info
log.info """=======================================================
DataQC v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']            = 'DataQC'
summary['Pipeline Version']         = workflow.manifest.version
summary['PLINK bfile']              = params.bfile
summary['VCF dataset']              = params.vcf
summary['FAM file']                 = params.fam
summary['Gen QC steps']             = params.gen_qc_steps
summary['Genome Build']             = params.genome_build
summary['S threshold']              = params.GenOutThresh
summary['Gen SD threshold']         = params.GenSdThresh
summary['Exp SD threshold']         = params.ExpSdThresh
summary['Contamination slope']       = params.ContaminationSlope
summary['Contamination area']       = params.ContaminationArea
summary['Expression matrix']        = params.expfile
summary['GTE file']                 = params.gte
summary['Max Memory']               = params.max_memory
summary['Max CPUs']                 = params.max_cpus
summary['Max Time']                 = params.max_time
summary['Cohort name']              = params.cohort_name
if(params.AdditionalCovariates!="$baseDir/data/1000G_pops.txt") summary['Additional covariates'] = params.AdditionalCovariates
if(params.InclusionList!="$baseDir/data/EmpiricalProbeMatching_AffyHumanExon.txt") summary['Inclusion list'] = params.InclusionList
if(params.ExclusionList!="$baseDir/data/EmpiricalProbeMatching_AffyHumanExon.txt") summary['Exclusion list'] = params.ExclusionList
if(params.preselected_sex_check_vars) summary['Pruned variants for sex check'] = params.preselected_sex_check_vars
summary['Expression platform']      = params.exp_platform
summary['Plink executable']         = params.plink_executable
summary['Plink 2 executable']       = params.plink2_executable
summary['Reference 1000G folder']   = params.reference_1000g_folder
summary['Chain folder']             = params.chain_path
summary['Output dir']               = params.outdir
summary['Working dir']              = workflow.workDir
summary['Container Engine']         = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']             = "$HOME"
summary['Current user']             = "$USER"
summary['Current path']             = "$PWD"
summary['Working dir']              = workflow.workDir
summary['Script dir']               = workflow.projectDir
summary['Config Profile']           = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="

process ListChromosomes {
    tag {ListChromosomes}

    input:
      file(input_vcf) from ( vcf_contig_count.value > 1 ? Channel.empty() : vcffile_ch)

    output:
      file("chromosomes.txt") into vcf_chromosome_list_file
      file("${input_vcf}.csi") into input_vcf_all_index

    when:
      vcf_contig_count.value == 1

    script:
      """
      bcftools index ${input_vcf}
      tabix -l ${input_vcf} > "chromosomes.txt"
      """
}

vcf_chromosome_list_file.splitText().map{ line -> line.trim() }.view().set { chromosome_channel }

process SplitVcf {

    tag {SplitVcf}

    input:
      file(input_vcf) from vcffile_ch2.collect()
      file(index_file) from input_vcf_all_index.collect()
      val(chr) from chromosome_channel

    output:
      tuple val(chr), file("split_${chr}.vcf.gz") into split_vcf
      tuple val(chr), file("split_${chr}.vcf.gz") into split_vcf_to_clean

    script:
      """
      bcftools view ${input_vcf} -r ${chr} -Oz -o "split_${chr}.vcf.gz"
      """
}

process ExpandVcfChannel {
    tag {ExpandVcfChannel}

    input:
      file(input_vcf) from ( vcf_contig_count.value > 1 ? vcffile_ch : Channel.empty())

    output:
      tuple env(chr), file(input_vcf) into ext_vcf_ch

    when:
      vcf_contig_count.value > 1

    script:
      """
      bcftools index ${input_vcf}
      tabix -l ${input_vcf} > chr.txt
      chr=`cat chr.txt | tr -d '\n'`
      """      
}

process WgsNorm {

    tag {WgsNorm}

    input:
      tuple val(chr), file(input_vcf) from ( vcf_contig_count.value == 1 ? split_vcf : ext_vcf_ch )

    output:
      tuple val(chr), file("norm.vcf.gz") into vcf_normalised
      tuple val(chr), file("norm.vcf.gz") into vcf_normalised_to_clean
      set val(chr), val(1) into clean_split_vcf_signal

    script:
      """
      echo "chromosome ${chr}"
      bcftools norm -m -any ${input_vcf} -Oz -o "norm.vcf.gz"
      """
}

process WgsQC {

    tag {WgsQC}

    input:
      tuple val(chr), file(input_vcf) from ( params.gen_qc_steps == 'WGS' ? vcf_normalised : Channel.empty() )

    output:
      tuple val(chr), file("norm-filtered.vcf.gz") into vcf_wgs_qced
      tuple val(chr), file("norm-filtered.vcf.gz") into vcf_wgs_qced_to_clean

      file("VCFFilterSummaryStats.txt") into wgs_qc_stats
      file("VCFFilterSettings.txt") into wgs_qc_settings

    when:
      params.gen_qc_steps == 'WGS'

    script:
      if (chr in ["X", "Y", "chrX", "chrY"]) 
      """
      python3 $baseDir/bin/custom_vcf_filter.py --input ${input_vcf} --hardy_weinberg_equilibrium 0 --output norm \
      | tee custom_vcf_filter.log
      
      python3 $baseDir/bin/print_WGS_VCF_filter_overview.py \
        --workdir . --chr ${chr} \
        --vcf_file_format "norm.vcf.gz"
      
      """
      else
      """
      python3 $baseDir/bin/custom_vcf_filter.py --input ${input_vcf} --output norm \
      | tee custom_vcf_filter.log
      
      python3 $baseDir/bin/print_WGS_VCF_filter_overview.py \
        --workdir .  --chr ${chr} \
        --vcf_file_format "norm.vcf.gz"

      """
}

process VcfToPlink {

    tag {VcfToPlink}

    input:
      tuple val(chr), file(vcf) from ( params.gen_qc_steps == 'WGS' ? vcf_wgs_qced : vcf_normalised )

    output:
      file("${chr}_converted_vcf.bed") into vcf_to_plink_bed_ch
      file("${chr}_converted_vcf.bim") into vcf_to_plink_bim_ch
      file("${chr}_converted_vcf.fam") into vcf_to_plink_fam_ch
      val("${chr}_converted_vcf") into vcf_to_plink_prefix_ch
      set val(chr), val(1) into clean_wgs_norm_signal
      set val(chr), val(1) into clean_wgs_qc_signal

    script:
      """
      # Make plink file
      plink2 --vcf ${vcf} --const-fid --split-x 'hg38' --make-bed --out "${chr}_converted_vcf"
      """
}

process MergePlink {

    tag {MergePlink}

    input:
      file(bed) from vcf_to_plink_bed_ch.collect()
      file(bim) from vcf_to_plink_bim_ch.collect()
      file(fam) from vcf_to_plink_fam_ch.collect()
      file(mergelist) from vcf_to_plink_prefix_ch.collectFile(name: 'mergelist.txt', newLine: true)
      
    output:
      tuple file("chrAll.bed"), file("chrAll.bim"), file("chrAll.fam") into plink_merged

    script:
      """
      plink2 --merge-list ${mergelist} --make-bed --out "chrAll"
      """
}

split_vcf_to_clean.mix(clean_split_vcf_signal).groupTuple(size: 2).view().set {clean_split_vcf_ready}

process CleanSplitVcf {
    tag {CleanSplitVcf}

    input:
        set val(chr), val(files_list) from clean_split_vcf_ready

    script:
    """
    # clean_work_files.sh "${files_list[0]}"
    """
}

vcf_normalised_to_clean.mix(clean_wgs_norm_signal).groupTuple(size: 2).view().set { clean_wgs_norm_ready }

process CleanWgsNorm {
    tag {CleanWgsNorm}

    input:
        set val(chr), val(files_list) from clean_wgs_norm_ready

    script:
    """
    # clean_work_files.sh "${files_list[0]}"
    """
}

vcf_wgs_qced_to_clean.mix(clean_wgs_qc_signal).groupTuple(size: 2).view().set { clean_wgs_qc_ready }

process CleanWgsQc {
    tag {CleanWgsQc}

    input:
        set val(chr), val(files_list) from clean_wgs_qc_ready

    script:
    """
    # clean_work_files.sh "${files_list[0]}"
    """
}

process GenotypeQC {

    tag {GenotypeQC}

    input:
      set file(bfile), file(bim), file(fam) from ( vcf_contig_count.value > 0 ? plink_merged : bfile_ch )
      file(fam_annot) from fam_annot_ch.ifEmpty { 'EMPTY' }
      file gte from gte_ch_gen
      val s_stat from params.GenOutThresh
      val sd_thresh from params.GenSdThresh
      val optional_pruned_variants_sex_check from params.preselected_sex_check_vars
      path ExclusionList from params.ExclusionList
      path InclusionList from params.InclusionList
      val genome_build from params.genome_build
      file(plink_executable) from plink_executable_ch.ifEmpty { 'EMPTY' }
      file(plink2_executable) from plink2_executable_ch.ifEmpty { 'EMPTY' }
      file(reference_1000g_folder) from reference_1000g_ch.ifEmpty { 'EMPTY' }
      file(chain_path) from chain_path_ch.ifEmpty { 'EMPTY' }

    output:
      path ('outputfolder_gen') into output_ch_genotypes
      file 'outputfolder_gen/gen_data_QCd/SexCheck.txt' into (sexcheck, sexcheck_expr)
      file 'outputfolder_gen/gen_data_QCd/*fam' into (sample_qc, sample_qc_expr)
      file '1000Gref.afreq.gz' into ref_allele_frequencies
      file 'target.afreq.gz' into target_allele_frequencies

    script:
    if (params.reference_1000g_folder == '')
      reference_1000g_prefix_arg = "--ref_1000g data/1000G_phase3_common_norel"
    else
      reference_1000g_prefix_arg = "--ref_1000g $reference_1000g_folder/1000G_phase3_common_norel"

    fam_arg = (params.fam != '') ? "--fam $fam_annot" : ""
    plink_arg = (params.plink_executable != '') ? "--plink_executable $plink_executable" : ""
    plink2_arg = (params.plink2_executable != '') ? "--plink2_executable $plink2_executable" : ""
    chain_path_arg = (params.chain_path != '') ? "--chain_path $chain_path" : ""

    """
    Rscript --vanilla $baseDir/bin/GenQcAndPosAssign.R  \
    --target_bed ${bfile} \
    $fam_arg \
    --genome_build ${genome_build} \
    --gen_exp ${gte} \
    --sample_list $baseDir/data/unrelated_reference_samples_ids.txt \
    --pops $baseDir/data/1000G_pops.txt \
    --S_threshold ${s_stat} \
    --SD_threshold ${sd_thresh} \
    --inclusion_list "${InclusionList}" \
    --exclusion_list "${ExclusionList}" \
    --output outputfolder_gen \
    --pruned_variants_sex_check "${optional_pruned_variants_sex_check}" \
    --liftover_path $baseDir/bin/liftOver \
    $plink_arg \
    $plink2_arg \
    $reference_1000g_prefix_arg \
    $chain_path_arg
    """
    
}

process SPLICING_SAMPLE_QC {

    tag {SPLICING_SAMPLE_QC}

    input:
      file combined_matrix from combined_matrix_ch
      file gte from gte_ch_exp
      file geno_filter from sample_qc
      file junction_site from junction_site_ch
      file ir_site from ir_site_ch

    output:
      path ('outputfolder_splice') into (output_ch_splicing_for_finalize, output_ch_splicing_for_expr)

    script:
    """
    mkdir -p outputfolder_splice outputfolder_splice/splice_data_summary

    # Analysed feature ids = the combined-matrix row ids (junction_<id> / IR_<id>). Used to restrict
    # the per-site read-count summary to the analysed set AND to feed the splicing-PCA pre-filter.
    zcat -f ${combined_matrix} | cut -f1 | tail -n +2 > feature_ids.txt

    # Splicing-specific summary FIRST: per-site RAW read-count depth + missingness (junction & IR), from
    # the Stage-0 .site files, restricted to the analysed feature set. Runs BEFORE ProcessSplicing so the
    # splicing-PCA pre-filter can select a well-measured basis. (Skips gracefully if a .site is absent.)
    python3 $baseDir/bin/compute_site_readcount_summary.py \\
      --junction_site ${junction_site} \\
      --ir_site ${ir_site} \\
      --keep_sites feature_ids.txt \\
      --output_dir outputfolder_splice/splice_data_summary

    # Splicing QC (ProcessSplicing.R = ProcessExpression.R analog): PCA/MDS outlier QC + per-feature
    # rank-INT. geno_filter (gen_data_QCd/*_ToImputation.fam) restricts QC to the genotyped samples,
    # exactly like ProcessExpression.R does for expression. This QC defines the analysed sample set.
    # --site_summary_dir enables the PRE-FILTER: the top-features PCA basis is restricted to features with
    # missingness < max_pca_missingness AND mean_count >= min_pca_meancount (full matrix still written out).
    Rscript --vanilla $baseDir/bin/ProcessSplicing.R \\
      --psi_matrix ${combined_matrix} \\
      --gte ${gte} \\
      --geno_filter ${geno_filter} \\
      --sd ${params.sd_threshold} \\
      --int_transform ${params.int_transform} \\
      --n_top_features ${params.n_top_features} \\
      --site_summary_dir outputfolder_splice/splice_data_summary \\
      --max_pca_missingness ${params.max_pca_missingness} \\
      --min_pca_meancount ${params.min_pca_meancount} \\
      --output_dir outputfolder_splice

    # Transpose features x samples -> samples x features (eQTLGen/HASE format)
    python3 $baseDir/bin/transpose_matrix.py \\
      --input outputfolder_splice/splice_data_QCd/splice_data_preprocessed.txt \\
      --output outputfolder_splice/splice_data_QCd/splice_data_preprocessed_T.txt \\
      --junction_list outputfolder_splice/splice_data_QCd/junction_list.txt
    mv outputfolder_splice/splice_data_QCd/splice_data_preprocessed_T.txt outputfolder_splice/splice_data_QCd/splice_data_preprocessed.txt

    md5sum outputfolder_splice/splice_data_QCd/* > outputfolder_splice/splice_data_QCd/md5sums.txt
    """
}

process EXPRESSION_SAMPLE_QC {

    tag {EXPRESSION_SAMPLE_QC}

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      file expression_matrix from expression_matrix_ch
      file gte from gte_ch_exp2
      path output_splice from output_ch_splicing_for_expr
      file sexcheck from sexcheck_expr
      file geno_filter from sample_qc_expr

    output:
      path ('outputfolder_exp') into (output_ch_expression, output_ch_expression_for_finalize)

    script:
    """
    mkdir -p outputfolder_exp

    # Expression QC: CPM-filter + TMM + per-gene INT (needed for the 4_ permuted/LD path), aligned to the
    # SPLICING-QC'd sample set; iterative PCA outlier detection is DIAGNOSTIC (reported, not removed,
    # unless --expr_remove_outliers true). Reads the QC'd splice sample list via cut -f1 (never loads the
    # wide splice matrix). Normalization + outlier logic mirror eQTLGen ProcessExpression.R.
    Rscript --vanilla $baseDir/bin/ProcessExpression.R \\
      --expression_matrix ${expression_matrix} \\
      --gte ${gte} \\
      --sample_list ${output_splice}/splice_data_QCd/splice_data_preprocessed.txt \\
      --sd ${params.sd_threshold} \\
      --remove_outliers ${params.expr_remove_outliers} \\
      --sex_info ${sexcheck} \\
      --emp_probe_mapping $baseDir/data/EmpiricalProbeMatching_RNAseq.txt \\
      --contamination_area ${params.ContaminationArea} \\
      --contamination_slope ${params.ContaminationSlope} \\
      --geno_filter ${geno_filter} \\
      --output outputfolder_exp

    md5sum outputfolder_exp/exp_data_QCd/* > outputfolder_exp/exp_data_QCd/md5sums.txt
    """
}

process FINALIZE_SPLICING {

    tag {FINALIZE_SPLICING}

    input:
      path output_splice from output_ch_splicing_for_finalize
      path output_exp    from output_ch_expression_for_finalize

    output:
      path ('outputfolder_splice') into output_ch_splicing

    script:
    """
    # ACTIVE expression QC (mod 1): re-align the splicing phenotype + covariate PCs to the final analysed
    # set E = genotype n splice n expression-QC (exp_data_summary/final_samples.txt). Rank-INT is monotone,
    # so re-ranking the already-INT'd phenotype restricted to E is BIT-EXACT INT-on-E (no raw-matrix
    # re-read). If expression QC removed nothing, this is an idempotent no-op.
    # Replace the staged input symlink with a writable copy so edits never touch upstream work files.
    cp -rL outputfolder_splice _splice_final
    rm outputfolder_splice
    mv _splice_final outputfolder_splice

    cd outputfolder_splice/splice_data_QCd
    # The on-disk phenotype is WIDE (samples x up-to-~670k features); data.table::fread segfaults on that
    # many columns and pandas.read_csv OOMs (per-column overhead). transpose_matrix_np.py reads the WIDE
    # with numpy (compact) and writes the TALL, so R can fread the tall (few columns) and re-INT per feature;
    # the TALL->WIDE step below uses pandas (the output is tall-read, few columns, safe). junction_list regen'd.
    python3 $baseDir/bin/transpose_matrix_np.py \\
      --input splice_data_preprocessed.txt --output _tall.txt --junction_list _samples_tmp.txt
    cd ../..

    Rscript --vanilla $baseDir/bin/finalize_splice_intersection.R \\
      --tall outputfolder_splice/splice_data_QCd/_tall.txt \\
      --splice_dir outputfolder_splice \\
      --final_samples ${output_exp}/exp_data_summary/final_samples.txt \\
      --num_splicing_pcs ${params.num_splicing_pcs}

    cd outputfolder_splice/splice_data_QCd
    python3 $baseDir/bin/transpose_matrix.py \\
      --input _tall.txt --output splice_data_preprocessed.txt --junction_list junction_list.txt
    rm -f _tall.txt _samples_tmp.txt
    md5sum splice_data_preprocessed.txt > md5sums.txt
    cd ../..
    """
}

wgs_qc_stats.ifEmpty("EMPTY").collectFile(name: 'wgs_qc_table_combined.txt', skip: 1, newLine: false, keepHeader: true).set { wgs_qc_stats_file_ch }
wgs_qc_settings.ifEmpty("EMPTY").collectFile(name: 'wgs_qc_settings_combined.txt', skip: 1, newLine: false, keepHeader: true).set { wgs_qc_settings_file_ch }

process RenderReport {

    tag {RenderReport}

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      path output_gen from output_ch_genotypes
      path output_splice from output_ch_splicing
      path output_expression from output_ch_expression
      path report from report_ch
      path ref_af from ref_allele_frequencies
      path target_af from target_allele_frequencies
      path sexcheck from sexcheck
      val num_gen_pcs from params.num_gen_pcs
      val num_splicing_pcs from params.num_splicing_pcs
      path additional_covariates from params.AdditionalCovariates
      file wgs_qc_stats_collected from wgs_qc_stats_file_ch
      file wgs_qc_settings_collected from wgs_qc_settings_file_ch

    output:
      path ('outputfolder_gen/*') into output_ch2
      path ('outputfolder_splice/*') into output_ch3
      path ('Report_DataQc*') into report_ch2
      path ('CovariatePCs.Splicing.txt') into splicing_covariates
      path ('CovariatePCs.Expression.txt') into expression_covariates
      path ('wgs_qc_table_combined.txt') into wgs_qc_stats_output_ch
      path ('wgs_qc_settings_combined.txt') into wgs_qc_settings_output_ch

    script:
    """
    # Two covariate files, one per Stage-4 path (MakeCovariateFile.R analog):
    #   CovariatePCs.Splicing.txt   = genotype PCs + SPLICING PCs + GenSex (+ additional)  -> empirical/meta-analysis (splicing)
    #   CovariatePCs.Expression.txt = genotype PCs + EXPRESSION PCs + GenSex (+ additional) -> permuted/LD (expression)
    python3 $baseDir/bin/generate_covariate_pcs.py \\
      --genotype_pcs outputfolder_gen/gen_PCs/GenotypePCs.txt \\
      --splicing_pcs outputfolder_splice/splice_PCs/splicing_PCs.txt \\
      --pc_prefix SplicingPC \\
      --gen_sex ${sexcheck} \\
      --num_gen_pcs ${num_gen_pcs} \\
      --num_splicing_pcs ${num_splicing_pcs} \\
      --additional_covariates ${additional_covariates} \\
      --output CovariatePCs.Splicing.txt

    python3 $baseDir/bin/generate_covariate_pcs.py \\
      --genotype_pcs outputfolder_gen/gen_PCs/GenotypePCs.txt \\
      --splicing_pcs outputfolder_exp/exp_PCs/exp_PCs.txt \\
      --pc_prefix ExpPC \\
      --gen_sex ${sexcheck} \\
      --num_gen_pcs ${num_gen_pcs} \\
      --num_splicing_pcs ${params.num_expression_pcs} \\
      --additional_covariates ${additional_covariates} \\
      --output CovariatePCs.Expression.txt

    # Single DataQC report: genotype QC + splicing QC + expression QC (merged; expression PCA/MDS/scree
    # are now a section of this report rather than a separate Report_ExpressionDiagnostic).
    cp -L ${report} notebook.Rmd
    R -e 'library(rmarkdown);rmarkdown::render("notebook.Rmd", "html_document",
    output_file = "Report_DataQc_${params.cohort_name}.html",
    params = list(
    dataqc_version = "${workflow.manifest.version}",
    dataset_name = "${params.cohort_name}",
    N = "outputfolder_splice/splice_data_QCd/splice_data_preprocessed.txt",
    N_expr = "outputfolder_exp/exp_data_summary/summary_table.txt",
    S = ${params.GenOutThresh},
    SD = ${params.GenSdThresh},
    SD_splice = ${params.sd_threshold},
    SD_exp = ${params.sd_threshold}))'
    """
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
