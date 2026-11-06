# eQTLGen sQTL (splicing) data QC pipeline

Automatic data quality check and processing for unimputed genotype data, **splicing phenotypes** and gene
expression. This is a splicing (sQTL) adaptation of the eQTLGen [DataQC](https://github.com/eQTLGen/DataQC)
pipeline: the genotype QC is unchanged, a **splicing QC** step replaces the expression QC as the primary
phenotype, and expression QC is retained (now active) — it detects and removes sample outliers (expression
PCA/MDS outliers and sex/contamination mismatches) and, together with the genotype and splicing QC, defines
the final shared analysed sample set.

Performs the following main steps:

- Genotypes
  - WGS-specific variant QC filtering (if input is WGS).
    - WGS-specific variant-level QC (VQSR thresholds, PASS labels, inbreeding coefficients, etc.).
  - Standard variant QC filtering (call-rate>0.95, Hardy-Weinberg P>1e-6, MAF>0.01).
  - Individual-level missingness filter <0.05.
  - Comparison of reported and genetic check, removal of mismatched samples.
  - Removal of samples with unclear genetic sex.
  - Removal of samples with excess heterozygosity (+/-3SD from the mean).
  - Removal of related samples (3d degree relatives). From each pair of related samples, one is kept in the data.
  - Visualisation of samples in the genetic reference space, instructions how to remove or split the data in case of multi-ancestry samples.
  - Removal of in-sample genetic outliers.
  - Calculates 10 first genetic principal components (PCs), used in analyses as covariates to correct for population stratification.
- Splicing (primary sQTL phenotype):
  - Aligns sample IDs between genotype and splicing samples (via the gte file), restricting to genotype-QC'd samples.
  - QCs the combined splice-site (PSI) matrix assembled in Stage 0 (`0_SplicingPhenotype`).
  - Iterative PCA / MDS outlier detection to remove sample outliers (run on the most-variable junctions, so standard PCA fits in modest memory — no Gram/dual PCA needed).
  - **Pre-filters the PCA basis** to well-measured features (per-site missingness < 0.10 **and** mean read depth per donor >= 5) before selecting the top (default 20,000) most-variable junctions, so the covariate PCs are driven by robust phenotypes rather than low-depth / high-missingness artifacts.
  - Per-feature rank inverse-normal transform (INT) of the full junction matrix (the phenotype).
  - Calculates 50 splicing PCs, used in analyses as covariates to correct for unknown variation.
  - Computes per-site read-depth & missingness summaries (junction + intron-retention), used for *post-hoc* QC in the meta-analysis.
- Gene expression (sample-level QC):
  - Aligns to the splicing-QC'd samples; CPM-filters, TMM-normalises and applies per-gene inverse normal transformation.
  - **Active QC**: iteratively removes expression PCA / MDS outlier samples.
  - The final analysed sample set is the intersection **genotype &cap; splicing &cap; expression**; the splicing phenotype and covariate PCs are re-aligned to this set.
  - Calculates up to 100 expression PCs and writes out the processed expression matrix.
- Additional steps:
  - Removes samples whose genetic sex does not match with the expression of sex-specific genes.
  - Reorders the genotype samples into random order.
  - Organises all the QCd data into the standard folder format.
  - Provides a single commented `html` QC report (genotype + splicing + expression) which should be used to get an overview of the quality of the data.

## Usage information

### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java >=11 installed.
- HPC has Singularity installed and running.

### Setup of the pipeline

You can either clone it by using git (if available in HPC):

`git clone https://github.com/boxiangliulab/DataQC.git`

Or just download this from the gitlab/github download link and unzip.

### Input files

- Unimputed genotype file in plink `.bed` format (https://www.cog-genomics.org/plink/1.9/input#bed). Genome build has to be in **hg19/GRCh37 (default)** or **hg38/GRCh38**. It is advisable that .fam file also includes observed sex for all samples (format: males=1, females=2), so that pipeline does extra check on that. However, if this information is not available for all samples, pipeline just skips this check. Input path has to be without `.bed/.bim/.fam` extension.
- We also allow unimputed `.vcf` datasets as input, specifically for WGS datasets for which we can do WGS-specific variant level QC. If your dataset concerns a VCF dataset, replace the `--bfile` argument with the `--vcf` argument. The `--vcf` argument expects a full path to a `.vcf` dataset. Pathname extension using globbing is allowed (using `*` or `?`), but the path should be provided without pathway extension.
- **Splicing phenotypes** produced by Stage 0 (`0_SplicingPhenotype`): the combined splice-site (PSI) matrix (`combined_matrix/<cohort>_combined_splice_site.txt.gz`) and the raw `.site` read-count files (`03_phenotype/<cohort>.site`, `04_IR/<cohort>_IR.site`). These are located automatically from `--splice_pheno_dir` by cohort name (or override individually with `--combined_matrix` / `--junction_site` / `--ir_site`).
- **Gene expression matrix**: raw, unprocessed gene counts (RNA-seq; genes in rows with stable ENSEMBL gene IDs, samples in columns), used for the expression-based sample QC. Provided via `--expression_matrix` (defaults to the Stage-0 location by cohort name).
- Genotype-to-expression linking file (gte). Tab-delimited, no header, 2 columns: sample ID in genotype data, corresponding sample ID in the RNA (splicing/expression) data.

### Required inputs

`--cohort_name`                 Name of the cohort.

`--bfile`                       Path to the unimputed genotype files in plink bed/bim/fam format (without extensions). Use `--vcf` instead for a VCF/WGS dataset.

`--genome_build`                Genome build of the cohort. Either hg19, GRCh37, hg38 or GRCh38. Defaults to hg19.

`--gen_qc_steps`                `Array` (default; generic SNP-array QC) or `WGS` (generic + WGS-specific QC; only valid with VCF datasets).

`--splice_pheno_dir`            Path to the Stage-0 (`0_SplicingPhenotype`) output folder. The combined splice matrix, gene-count matrix and raw `.site` files are derived from it by cohort name.

`--expression_matrix`           Path to the raw gene-count matrix (RNA-seq) used for the expression-based sample QC.

`--gte`                         Genotype-to-expression linking file. Tab-delimited, no header. First column: sample ID for genotype data. Second column: corresponding sample ID for gene expression data. Can be used to filter samples from the analysis.

`--outdir`                      Path to the output directory.

### Additional settings

There are six arguments which can be used to adjust certain outlier detection thresholds. These should be adjusted after initial run with the default settings and after investigating the diagnostic plots in the `Report_DataQc_[cohort name].html`. Then the pipeline should be re-run with adjusted settings.

`--GenOutThresh`        Threshold for a genotype sample genetic outlier, based on the LOF "outlierness" metric. Default 0.4.

`--GenSdThresh`         Threshold for a genotype sample genetic outlier, based on deviation from the mean of the first two genetic PCs. Default 3 SD.

`--sd_threshold`        SD threshold for splicing / expression PCA/MDS sample-outlier detection. Default 4.

`--ExpSdThresh`         SD threshold for expression sample outliers. Default 4.

`--ContaminationArea`   Threshold (degrees) for declaring samples contaminated from the X-vs-Y chromosome expression plot. Default 30.

`--ContaminationSlope`  Slope (degrees, 0-90) of the line separating males/females in that plot. Default 45.

Optional arguments:

`--expr_remove_outliers` `true` (default) = **active** expression QC: expression PCA/MDS and sex/contamination outliers are removed and the splicing set is re-aligned to the genotype &cap; splicing &cap; expression intersection. `false` = diagnostic only (flagged, not removed; splicing QC alone defines the set).

`--max_pca_missingness` / `--min_pca_meancount`  Splicing-PCA pre-filter cutoffs (defaults 0.10 and 5): a feature enters the top-features PCA basis only if its per-site missingness is below `max_pca_missingness` and its mean read depth per donor is at least `min_pca_meancount`.

`--num_gen_pcs` / `--num_splicing_pcs` / `--num_expression_pcs`  Number of PCs written to the covariate files. Defaults 10 / 50 / 100.

`--int_transform`       Apply the per-feature rank-INT to the splicing phenotype. Default true.

`--AdditionalCovariates` Tab-separated file with additional covariates (first column header `SampleID`, following columns named covariates; categorical values as text, one-hot encoded for you). Must cover every QC-passing sample, no NAs.

`--InclusionList` / `--ExclusionList`  Files with genotype IDs to keep / remove (one per row, no header). On overlap, the intersection is kept.

`--preselected_sex_check_vars`  plink ranges file defining which variants to use for the sex check.

`--fam`                 PLINK `.fam` file; takes precedence over the `.fam` supplied with `--bfile`.

`--plink_executable` / `--plink2_executable` / `--reference_1000g_folder`  Local paths to plink / plink2 / the 1000G reference; by default these are downloaded automatically. Use when working offline.

`--chain_path`          Folder with the `hg19ToHg38`/`hg38ToHg19` liftOver chain files. Only needed for **hg38/GRCh38** cohorts (to align the hg19 1000G reference); hg19 cohorts never use it. By default the chain is fetched automatically, but a local copy (e.g. `deps/chain`) is more robust where the download is unreachable.

### Running the data QC command

Modify the submit template `submit_DataQc_pipeline_template.sh` with your input paths. The template defaults
to a PBS scheduler; commented examples for loading the environment (conda / modules) and for other schedulers
are included. Only a few variables need filling (cohort name, genotype input, gte, build, Stage-0 folder);
the PC counts and thresholds come from the pipeline defaults above and can be overridden by adding the flag.

The core command is:

```bash
NXF_VER=21.10.6 ${nextflow_path}/nextflow run DataQC.nf \
  $geno_arg --cohort_name "$cohort_name" --genome_build "$genome_build" --gen_qc_steps "$gen_qc_steps" \
  --splice_pheno_dir "$splice_pheno_dir" --expression_matrix "$expression_matrix" --gte ${gte_path} \
  --outdir ${output_path} \
  -profile ${NEXTFLOW_PROFILES} -resume \
  -with-report logs/dataqc_report.html -with-trace logs/dataqc_trace.txt
```

where `$geno_arg` is `--bfile <path>` (default) or `--vcf <path>`. Select the HPC scheduler via the profile:

- Slurm: `-profile slurm,singularity`
- PBS/TORQUE: `-profile pbs,singularity`
- SGE: `-profile sge,singularity`

Submit the script (e.g. `qsub submit_DataQc_pipeline_template.sh`). This initiates the pipeline, builds the
analysis environment (Singularity) and submits the steps in the correct order; a separate `work` directory
holds all interim files. On a crash (e.g. walltime) just resubmit — Nextflow resumes from incomplete steps.

### Monitoring and debugging

- Monitoring:
  - Monitor the scheduler log file (e.g. `1_dataqc.o[JobID]` on PBS) and check if all the steps finish without error. Trick: command `watch tail -n 20 1_dataqc.o[JobID]` helps you to interactively monitor the status of the jobs.
  - Use `qstat -u [YourUserName]` to see if individual tasks are in the queue.
- If the pipeline crashes (e.g. due to walltime), you can just resubmit the same script after the fixes. Nextflow does not rerun completed steps and continues only from the steps which had not completed.
- When the work has finished, download and check the job report. This file is automatically written when the pipeline is run with `-with-report`, for potential errors or warnings. E.g. `logs/dataqc_report.html`.
- When you need to do some debugging, then you can use the last section of aforementioned report to figure out in which subfolder from `work` folder the actual step was run. You can then navigate to this folder and investigate the following hidden files:
  - `.command.sh`: script which was submitted.
  - `.command.log`: log file for seeing the analysis outputs/errors.
  - `.command.err`: file which lists the errors, if any.

### Output

Pipeline makes the following output (most relevant files outlined):

```
|--output
    |--outputfolder_gen
        |--gen_data_QCd/*_ToImputation.{bed,bim,fam,log}
        |--gen_PCs/GenotypePCs.txt
        |--gen_data_summary/ ...
        |--gen_plots/ ...
    |--outputfolder_splice
        |--splice_data_QCd/splice_data_preprocessed.txt   (samples x junctions, rank-INT phenotype)
        |--splice_PCs/splicing_PCs.txt
        |--splice_data_summary/ (summary_table.txt, site_readcount_summary_*.tsv, ...)
        |--splice_plots/ ...
    |--outputfolder_exp
        |--exp_data_QCd/exp_data_preprocessed.txt         (processed expression matrix)
        |--exp_PCs/exp_PCs.txt
        |--exp_data_summary/ (raw_gene_summary.txt, processed_gene_summary.txt, ...)
        |--exp_plots/ ...
    |--CovariatePCs.Splicing.txt      (genotype PCs + splicing PCs + GenSex)
    |--CovariatePCs.Expression.txt    (genotype PCs + expression PCs + GenSex)
    |--Report_DataQc_[cohort name].html
    |--pipeline_info/ ...
```

#### Steps to take

1. Investigate `Report_DataQc_[cohort name].html`, fix any genotype / splicing / expression issues according to the plots and instructions (adjust the threshold arguments above and re-run if needed).

2. If data had quality issues, re-run once resolved and re-check the report.

When all issues are solved:

3. `output/outputfolder_gen/gen_data_QCd/*_ToImputation.{bed,bim,fam}` are the QCd genotype files for the imputation pipeline (`https://gitlab.com/eqtlgen-group/eqtlgen-imputation-pipeline`).

4. The whole `output` folder is the input for the per-cohort preparations / encoding pipeline, which uses the QCd splicing phenotype + covariate files to run the sQTL encoding and partial-derivative calculation, and extracts QC files for sharing with the central site — most importantly `Report_DataQc_[cohort name].html` and the per-site/per-gene summary statistics.

## Genotype QC only

To only run genotype QC without the splicing / expression QC, use `DataQCgenoOnly.nf` (inherited from upstream; not adapted/validated for the splicing report).

## Acknowledgements

Genotype QC and covariate preparations make extensive use of the [bigsnpr package](https://privefl.github.io/bigsnpr/) and [plink 2](https://www.cog-genomics.org/plink/2.0/).

Gene expression processing makes use of [preprocesscore](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) R packages.

### Citation

[Privé, F., Luu, K., Blum, M. G. B., McGrath, J. J., &#38; Vilhjálmsson, B. J. (2020). Efficient toolkit implementing best practices for principal component analysis of population genetic data. <i>Bioinformatics</i>, <i>36</i>(16), 4449–4457. https://doi.org/10.1093/BIOINFORMATICS/BTAA520](https://academic.oup.com/bioinformatics/article/36/16/4449/5838185)

[Chang, C. C., Chow, C. C., Tellier, L. C. A. M., Vattikuti, S., Purcell, S. M., &#38; Lee, J. J. (2015). Second-generation PLINK: Rising to the challenge of larger and richer datasets. GigaScience, 4(1). https://doi.org/10.1186/s13742-015-0047-8](https://academic.oup.com/gigascience/article/4/1/s13742-015-0047-8/2707533)

[Bolstad B (2021). preprocessCore: A collection of pre-processing functions. R package version 1.56.0, https://github.com/bmbolstad/preprocessCore.](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html)

[Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140. doi: 10.1093/bioinformatics/btp616.](10.1093/bioinformatics/btp616)

### Contacts

For this splicing (sQTL) adaptation: Tian Chi, e1101920@u.nus.edu.

For the original eQTLGen DataQC pipeline: urmo.vosa at gmail.com.
