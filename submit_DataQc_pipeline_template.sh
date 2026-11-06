#!/bin/bash
# =============================================================================
# 1_DataQC — genotype QC + splicing QC + expression QC + covariates + report (one pipeline).
# Input: the combined splice matrix from 0_SplicingPhenotype. Output (-> ../output/):
#   CovariatePCs.Splicing.txt, CovariatePCs.Expression.txt, Report_DataQc_<cohort>.html
#   (share the report with central analysts and await review BEFORE 2_Imputation).
# Run from inside 1_DataQC/DataQC/ :   qsub submit_DataQc_pipeline_template.sh
# =============================================================================
#PBS -P 11003054
#PBS -q normal
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=48:00:00
#PBS -m bea
#PBS -N 1_dataqc

#### Load an environment with Nextflow, Singularity and Java (>=11, <21). ####

## Singapore cluster:
#source ~/.bashrc
#conda activate nf
#cd "$PBS_O_WORKDIR"   # go to the submission dir if your scheduler doesn't (e.g. PBS); not needed for SLURM

## UT HPC:
#module load java-1.8.0_40
#module load singularity/3.5.3
#module load squashfs/4.4

export SINGULARITY_CACHEDIR=../../singularitycache
export NXF_HOME=../../nextflowcache

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f

# Define paths
nextflow_path=../../tools

# =========================================================
# USER CONFIGURATION REQUIRED
# =========================================================

# Genotype Data
#input: fill bfile_path (PLINK, default) OR uncomment input_vcf (VCF/WGS)
bfile_path=[full path to your input genotype files without .bed/.bim/.fam extension]
#input_vcf=[path to input genotype VCF, bgzipped .vcf.gz with unique variant IDs]
gen_qc_steps=[WGS for whole-genome-sequencing genotypes | Array for SNP-array genotypes]
genome_build=hg38

# Other data
cohort_name=[name of the cohort]
gte_path=[full path to your genotype-to-expression file]
# (pre-filled defaults; change only if needed). splice_pheno_dir = the Stage-0 output folder; DataQC derives
# the combined splice matrix, gene-count matrix and .site files from it by cohort_name (so define cohort_name above).
splice_pheno_dir=../../0_SplicingPhenotype/output
expression_matrix=../../0_SplicingPhenotype/output/expression_matrix/${cohort_name}_gene_counts.txt.gz
output_path=../output

# Profiles: keep 'singularity'; change 'pbs' to your scheduler (slurm, sge, local, ...) or site profile.
NEXTFLOW_PROFILES="pbs,singularity"
# =========================================================


# Additional settings and optional arguments for the command

# --GenOutThresh [numeric threshold]
# --GenSdThresh [numeric threshold]
# --ExpSdThresh [numeric threshold]
# --ContaminationArea [number between 0 and 90, default 30]
# --ContaminationSlope [number between 0 and 90, default 45]
# --InclusionList [file with the list of samples to restrict the analysis]
# --ExclusionList [file with the list of samples to remove from the analysis]
# --preselected_sex_check_vars "data/Affy6_pruned_chrX_variant_positions.txt"
# --AdditionalCovariates [file with additional covariates. First column should be `SampleID`]
# --fam [PLINK .fam file. Takes precedence over .fam file supplied with `--bfile`]
# --plink_executable [path to plink executable (PLINK v1.90b6.26 64-bit)]
# --plink2_executable [path to plink2 executable (PLINK v2.00a3.7LM 64-bit Intel)]
# --reference_1000g_folder [path to folder with 1000G reference data]
# --chain_path [folder with hg19->hg38 and hg38->hg19 chain files]
# --combined_matrix / --expression_matrix / --junction_site / --ir_site [override individual Stage-0 inputs instead of --splice_pheno_dir]
# --expr_remove_outliers [default true = active expression QC (remove outliers, intersect sets); false = diagnostic only]

# Select the genotype input flag from whichever variable you filled above
if   [ -n "$input_vcf" ];  then geno_arg="--vcf $input_vcf"
elif [ -n "$bfile_path" ]; then geno_arg="--bfile $bfile_path"
else echo "ERROR: set either input_vcf (VCF) or bfile_path (PLINK) above"; exit 1; fi

NXF_VER=21.10.6 ${nextflow_path}/nextflow run DataQC.nf \
  $geno_arg --gen_qc_steps "$gen_qc_steps" --genome_build "$genome_build" \
  --splice_pheno_dir "$splice_pheno_dir" --expression_matrix "$expression_matrix" --gte ${gte_path} --cohort_name "$cohort_name" \
  --chain_path deps/chain \
  --outdir ${output_path} \
  -profile ${NEXTFLOW_PROFILES} -resume \
  -with-report logs/dataqc_report.html -with-trace logs/dataqc_trace.txt

echo ""
echo "==================================================================="
echo " DataQC complete. SHARE: ${output_path}/Report_DataQc_${cohort_name}.html"
echo " Await central review, THEN run 2_Imputation."
echo "==================================================================="


