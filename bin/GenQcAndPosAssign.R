#!/usr/bin/env Rscript

library(bigreadr)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(data.table)
library(optparse)
library(patchwork)
library(stringr)
library(rmarkdown)
library(Cairo)
library(igraph)

system_verbose <- function(..., verbose) {
  system(..., ignore.stdout = !verbose, ignore.stderr = !verbose)
}

#' Relationship-based pruning
#'
#' Quality Control based on KING-robust kinship estimator. More information can
#' be found at \url{https://www.cog-genomics.org/plink/2.0/distance#king_cutoff}.
#'
#' @param plink2.path Path to the executable of PLINK 2.
#' @inheritParams snp_plinkIBDQC
#' @param thr.king  Note that KING kinship coefficients are scaled such that
#'   duplicate samples have kinship 0.5, not 1. First-degree relations
#'   (parent-child, full siblings) correspond to ~0.25, second-degree relations
#'   correspond to ~0.125, etc. It is conventional to use a cutoff of ~0.354
#'   (2^-1.5, the geometric mean of 0.5 and 0.25) to screen for monozygotic
#'   twins and duplicate samples, ~0.177 (2^-2.5) to remove first-degree
#'   relations as well, and ~0.0884 (2^-3.5, **default**) to remove
#'   second-degree relations as well, etc.
#' @param extra.options Other options to be passed to PLINK2 as a string.
#' @param make.bed Whether to create new bed/bim/fam files (default).
#'   Otherwise, returns a table with coefficients of related pairs.
#'
#' @return See parameter `make-bed`.
#' @export
#'
#' @inherit snp_plinkQC references
#' @references
#' Manichaikul, Ani, Josyf C. Mychaleckyj, Stephen S. Rich, Kathy Daly,
#' Michele Sale, and Wei-Min Chen. "Robust relationship inference in genome-wide
#' association studies." Bioinformatics 26, no. 22 (2010): 2867-2873.
#'
#' @seealso [download_plink2] [snp_plinkQC]
#'
#' @examples
#' \dontrun{
#'
#' bedfile <- system.file("extdata", "example.bed", package = "bigsnpr")
#' plink2 <- download_plink2(AVX2 = FALSE)
#'
#' bedfile2 <- snp_plinkKINGQC(plink2, bedfile,
#'                             bedfile.out = tempfile(fileext = ".bed"),
#'                             ncores = 2)
#'
#' df_rel <- snp_plinkKINGQC(plink2, bedfile, make.bed = FALSE, ncores = 2)
#' str(df_rel)
#' }
#'
snp_plinkKINGQC <- function(plink2.path,
                            bedfile.in,
                            bedfile.out = NULL,
                            thr.king = 2^-3.5,
                            make.bed = TRUE,
                            ncores = 1,
                            extra.options = "",
                            verbose = TRUE) {

  # check PLINK version
  v <- system(paste(plink2.path, "--version"), intern = TRUE)
  if (substr(v, 1, 8) != "PLINK v2")
    stop2("This requires PLINK v2; got '%s' instead.", v)

  # get file without extension
  prefix.in <- sub_bed(bedfile.in)

  # get possibly new file
  if (make.bed) {

    if (is.null(bedfile.out)) bedfile.out <- paste0(prefix.in, "_norel.bed")
    assert_noexist(bedfile.out)

    # compute KING-robust kinship coefficients and filter
    system_verbose(
      paste(
        plink2.path,
        "--bfile", prefix.in,
        "--make-bed --king-cutoff", thr.king,
        "--out", sub_bed(bedfile.out),
        "--threads", ncores,
        extra.options
      ),
      verbose = verbose
    )

    bedfile.out

  } else {

    prefix.out <- tempfile()

    # compute table of KING-robust kinship coefficients
    system_verbose(
      paste(
        plink2.path,
        "--bfile", prefix.in,
        "--make-king-table --king-table-filter", thr.king,
        "--out", prefix.out,
        "--threads", ncores,
        extra.options
      ),
      verbose = verbose
    )

    rel_df <- bigreadr::fread2(paste0(prefix.out, ".kin0"), header = TRUE, keepLeadingZeros = TRUE,
                               colClasses = list(character = c(1,2,3,4)),
                               nThread = ncores)
    names(rel_df) <- sub("^#(.*)$", "\\1", names(rel_df))
    rel_df

  }
}

print(system_verbose)

# Function
read_fam <- function(path) {
  NAMES.FAM <- c("family.ID", "sample.ID", "paternal.ID",
                 "maternal.ID", "sex", "affection")

  pattern <- "\\.bed$"
  
  if (!grepl(pattern, path)) {
    famfile <- paste0(path, ".fam")
  } else {
    famfile <- sub(pattern, ".fam", path)
  }

  fam <- bigreadr::fread2(famfile, col.names = NAMES.FAM, keepLeadingZeros = TRUE,
                          colClasses = list(character = c(1,2)), nThread = 1)

  return(fam)
}

# Function modified from bigsnpr to work with offline chain files
snp_modifyBuild2 <- function(info_snp,
                             liftOver,
                             from = "hg18",
                             to = "hg19",
                             check_reverse = TRUE,
                             chain_path = NULL) {

  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop2("Expecting variables 'chr' and 'pos' in input 'info_snp'.")

  # Make sure liftOver is executable
  # comment out, not working in HPC
  #liftOver <- make_executable(normalizePath(liftOver))

  # Need BED UCSC file for liftOver
  info_BED <- with(info_snp, data.frame(
    # sub("^0", "", c("01", 1, 22, "X")) -> "1"  "1"  "22" "X"
    chrom = paste0("chr", sub("^0", "", chr)),
    start = pos - 1L, end = pos,
    id = seq_along(pos)))

  BED <- tempfile(fileext = ".BED")
  bigreadr::fwrite2(stats::na.omit(info_BED),
                    BED, col.names = FALSE, sep = " ", scipen = 50)

  # Need chain file
  # url <- paste0("ftp://hgdownload.cse.ucsc.edu/goldenPath/", from, "/liftOver/",
  #                from, "To", tools::toTitleCase(to), ".over.chain.gz")
  chain <- tempfile(fileext = ".over.chain.gz")
  chain_file <- paste0(chain_path, "/", from, "To", tools::toTitleCase(to), ".over.chain.gz")
  message(chain_file)
  #utils::download.file(url, destfile = chain, quiet = TRUE)
  file.copy(chain_file, chain)

  # Run liftOver (usage: liftOver oldFile map.chain newFile unMapped)
  lifted <- tempfile(fileext = ".BED")

  system2(liftOver, c(BED, chain, lifted, tempfile(fileext = ".txt")))
  message("Liftover done")
  # Read the ones lifter + some QC
  new_pos <- bigreadr::fread2(lifted, nThread = 1)
  is_bad <- vctrs::vec_duplicate_detect(new_pos$V4) |
    (new_pos$V1 != info_BED$chrom[new_pos$V4])
  new_pos <- new_pos[which(!is_bad), ]

  pos0 <- info_snp$pos
  info_snp$pos <- NA_integer_
  info_snp$pos[new_pos$V4] <- new_pos$V3

  if (check_reverse) {
    pos2 <- Recall(info_snp, liftOver, from = to, to = from, check_reverse = FALSE, chain_path = chain_path)$pos
    info_snp$pos[pos2 != pos0] <- NA_integer_
  }

  bigassertr::message2("%d variants have not been mapped.", sum(is.na(info_snp$pos)))

  info_snp
}

# Argument parser
option_list <- list(
    make_option(c("-t", "--target_bed"), type = "character",
    help = "Name of the target genotype file (bed/bim/fam format). Required file extension: .bed."),
    make_option(c("-f", "--fam"), type = "character", default = NULL,
    help = "Path to a separate fam file. Has priority over fam associated with --target_bed"),
    make_option(c("-g", "--gen_exp"), type = "character",
    help = "Tab-delimited genotype-to-expression sample ID linking file."),
    make_option(c("-s", "--sample_list"), type = "character",
    help = "Path to the file listing unrelated samples for reference data (tab-delimited .txt)."),
    make_option(c("-p", "--pops"), type = "character",
    help = "Path to the file indicating the population for each sample in reference data."),
    make_option(c("-a", "--pruned_variants_sex_check"), type = "character",
    help = "Path to a file with pruned X-chromosome variants to use in the sex check"),
    make_option(c("-o", "--output"), type = "character", help = "Folder with all the output files."),
    make_option(c("-S", "--S_threshold"), default = 0.4,
    help = paste0("Numeric threshold to declare samples outliers, based on the genotype PCs. ", 
                  "Defaults to 0.4 but should always be visually checked and changed, if needed.")),
    make_option(c("-d", "--SD_threshold"), default = 0.4,
    help = paste0("Numeric threshold to declare samples outliers, based on the genotype PCs. ", 
                  "Defaults to 0.4 but should always be visually checked and changed, if needed.")),
    make_option(c("-i", "--inclusion_list"), type = "character",
    help = "Path to the file with sample IDs to include."),
    make_option(c("-e", "--exclusion_list"), type = "character",
    help = "Path to the file with sample IDs to exclude. This also removes samples from inclusion list."),
    make_option(c("-b", "--genome_build"), type = "character",
    help = "Genome build of the target genotype file."),
    make_option(c("--liftover_path"), type = "character",
    help = "Liftover executable."),
    make_option(c("--plink_executable"), type = "character", default = NULL,
                help = "Plink executable."),
    make_option(c("--plink2_executable"), type = "character", default = NULL,
                help = "Plink2 executable."),
    make_option(c("--ref_1000g"), type = "character", default = NULL,
                help = "reference 1000g prefix."),
    make_option(c("--chain_path"), type = "character", default = NULL,
                help = "Folder with liftOver chain files.")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Remove the check of parallel blas
options(bigstatsr.check.parallel.blas = FALSE)

# Report settings
print(args$target_bed)
print(args$fam)
print(args$genome_build)
print(args$gen_exp)
print(args$sample_list)
print(args$pops)
print(args$pruned_variants_sex_check)
print(args$output)
print(args$S_threshold)
print(args$SD_threshold)
print(args$exclusion_list)
print(args$liftover_path)
print(args$plink_executable)
print(args$plink2_executable)
print(args$chain_path)

if (!is.numeric(args$S_threshold) | !is.numeric(args$SD_threshold) | !is.numeric(args$SD_threshold)){
  message("Some of the QC thresholds is not numeric!")
  stop()
}

# Map genome builds to build codes.
build_code <- "b37"
ucsc_code <- "hg19"
variant_format <- r"(@:#[b37]\$r,\$a)"

if (args$genome_build %in% c("hg19", "GRCh37")) {
  message("Using genome build hg19/GRCh37")
} else if (args$genome_build %in% c("hg38", "GRCh38")) {
  message("Using genome build hg38/GRCh38")
  variant_format <- r"(@:#[b38]\$r,\$a)"
  build_code <- "b38"
  ucsc_code <- "hg38"
} else if (args$genome_build %in% c("hg18", "GRCh36")) {
  message("Using genome build hg18/GRCh36")
  variant_format <- r"(@:#[b36]\$r,\$a)"
  build_code <- "b36"
  ucsc_code <- "hg18"
} else {
  stop(sprintf("Genome build %s is not recognized as an available genome build!", args$genome_build))
}

bed_simplepath <- stringr::str_replace(args$target_bed, ".bed", "")

# Make output folder structure
dir.create(args$output)
dir.create(paste0(args$output, "/gen_plots"))
dir.create(paste0(args$output, "/gen_data_QCd"))
dir.create(paste0(args$output, "/gen_PCs"))
dir.create(paste0(args$output, "/gen_data_summary"))

# Download plink executables
make_executable <- function(exe) {
  Sys.chmod(exe, mode = (file.info(exe)$mode | "111"))
}

PLINK <- args$plink_executable
PLINK2 <- args$plink2_executable

if (is.null(PLINK2) || PLINK2 == "" || !file.exists(PLINK2)) {
  message(sprintf("PLINK 2 executable empty, or not found at %s.", PLINK2))
  message("Attempting to download PLINK 2 executable")

  dir.create("plink")

  # Download plink 2 executable
  utils::download.file("https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250915.zip",
                       destfile = "plink/plink2.zip", verbose = TRUE)
  PLINK2 <- utils::unzip("plink/plink2.zip",
                        files = "plink2",
                        exdir = "plink")
} else {
  PLINK2 <- normalizePath(PLINK2) 
  message(sprintf("PLINK 2 executable found at %s.", PLINK2))
}

make_executable(PLINK2)

if (is.null(PLINK) || PLINK == "" || !file.exists(PLINK)) {
  message(sprintf("PLINK 1.9 executable empty, or not found at %s.", PLINK))
  message("Attempting to download PLINK 1.9 executable")

  dir.create("plink")
  
  # Download plink 1.9 executable
  utils::download.file("https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip",
  destfile = "plink/plink.zip", verbose = TRUE)
  PLINK <- utils::unzip("plink/plink.zip",
                          files = "plink",
                          exdir = "plink")
} else {
  PLINK <- normalizePath(PLINK)
  message(sprintf("PLINK 1.9 executable found at %s.", PLINK))
}

make_executable(PLINK)

ref_1000g_prefix <- "data"
if (!is.null(args$ref_1000g) && args$ref_1000g != "") {
  if (endsWith(args$ref_1000g, "1000G_phase3_common_norel"))
  ref_1000g_prefix <- args$ref_1000g
}

if (file.exists(paste0(ref_1000g_prefix, ".bed"))
  & file.exists(paste0(ref_1000g_prefix, ".bim"))
  & file.exists(paste0(ref_1000g_prefix, ".fam"))) {
  message(paste0("found 1000G reference at ", ref_1000g_prefix, "'.<bim/bed/fam>'."))
} else {
  # Download subsetted 1000G reference
  message(paste0("1000G reference does not exist at ", ref_1000g_prefix, "'.<bim/bed/fam>'."))
  message("Attempting to download the 1000G reference data")
  bedfile <- download_1000G(dirname(ref_1000g_prefix))
}

## Chain files for LiftOver
chain_path <- args$chain_path
if (file.exists(paste0(chain_path, "/hg19ToHg38.over.chain.gz")) & 
file.exists(paste0(chain_path, "/hg38ToHg19.over.chain.gz"))){message(paste0("Found liftOver chain files at ", chain_path))}

## Calculate AFs for reference data
system(paste0(PLINK2, " --bfile ", ref_1000g_prefix, " --threads 4 --freq 'cols=+pos' --out 1000Gref"))

if ("hg19" != ucsc_code) {
  target_frequencies <- fread("1000Gref.afreq", sep="\t", data.table=F, header=T,
                 col.names=c("chr", "pos", "ID", "REF", "ALT", "PROVISIONAL_REF", "ALT_FREQS", "OBS_CT"))

if (!is.null(args$chain_path) && args$chain_path != "") {
  target_frequencies_mapped <- snp_modifyBuild2(
    target_frequencies, file.path(".", R.utils::getRelativePath(args$liftover_path)),
    from = "hg19", to = ucsc_code, chain_path = chain_path)
}else{
  target_frequencies_mapped <- snp_modifyBuild(
    target_frequencies, file.path(".", R.utils::getRelativePath(args$liftover_path)),
    from = "hg19", to = ucsc_code)
}
  colnames(target_frequencies_mapped)[1:2] <- c("#CHROM", "POS")

  fwrite(target_frequencies_mapped[!is.na(target_frequencies_mapped$POS),], "1000Gref.afreq.gz", col.names=T, row.names=F, quote=F, sep="\t")

  system("rm 1000Gref.afreq")
} else {
  system("gzip 1000Gref.afreq --force")
}

# Target data
## Original file
message("Read in target data.")
target_bed <- bed(args$target_bed)
target_bed$.fam <- read_fam(args$target_bed)

## Calculate AFs for target data
system(paste0(PLINK2, " --bfile ", bed_simplepath, " --threads 4 --freq 'cols=+pos' --out target"))
system("gzip target.afreq --force")

# eQTL samples

if(args$gen_exp == "NA"){
	gte <- as.data.table(data.frame(V1 = target_bed$.fam$`sample.ID`))
} else {
	gte <- fread(args$gen_exp, sep = "\t", header = FALSE,
             keepLeadingZeros = TRUE,
             colClasses = "character")
}

print(target_bed$ncol)
print(target_bed$nrow)
print(nrow(gte[gte$V1 %in% target_bed$.fam$`sample.ID`, ]))


print(str(gte))

summary_table <- data.frame(stage = "Raw file", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = target_bed$nrow,
Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% target_bed$.fam$`sample.ID`, ]))

if(nrow(gte[gte$V1 %in% target_bed$.fam$`sample.ID`, ]) < 100){stop("Less than 100 samples are in genotype-to-expression file!")}

# Prepare and normalise fam file
#
fam <- target_bed$fam
fam$`family.ID` <- '0'

if (any(duplicated(fam$`sample.ID`))) {
  stop(sprintf("Error! samples in PLINK fam file are not unique. Exiting"))
}

if (!is.null(args$fam) && args$fam != "") {
  new_fam <- fread(args$fam, data.table = FALSE, header = FALSE, col.names = colnames(fam),
                   keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
  new_fam$`family.ID` <- '0'

  # Check if all sample ids in new fam are unique
  if (any(duplicated(new_fam$IID))) {
    stop(sprintf("Error! samples in fam file '%s' are not unique. Exiting", args$fam))
  }
  # Check if all sample ids from plink fam are in new fam
  if (!all(fam$`sample.ID` %in% new_fam$`sample.ID`)) {
    stop(sprintf("Error! samples in PLINK fam file are not all in '%s'. Exiting", args$fam))
  }
  # Check if all sample ids from new fam are in plink fam
  if (!all(new_fam$`sample.ID` %in% fam$`sample.ID`)) {
    stop(sprintf("Error! samples in fam file '%s' are not all in PLINK fam file. Exiting", args$fam))
  }

  fam <- new_fam[order(match(new_fam$`sample.ID`, fam$`sample.ID`)),]
}

# Write normalized fam
fwrite(fam, "fam_normalized.fam", col.names=F, row.names=F, quote=F, sep="\t")

## If specified, keep in only samples which are in the sample whitelist
if (args$inclusion_list != "" & args$inclusion_list != "EmpiricalProbeMatching_AffyHumanExon.txt"){
  inc_list <- fread(args$inclusion_list, header = FALSE,
                    keepLeadingZeros = TRUE, colClasses = "character")
  samples_to_include <- fam[fam$`sample.ID` %in% inc_list$V1, ]
  message("Sample inclusion filter active!")

  temp_QC <- data.frame(stage = "Samples in inclusion list",
                        Nr_of_SNPs = target_bed$ncol,
                        Nr_of_samples = nrow(samples_to_include),
                        Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% samples_to_include$`sample.ID`, ]))
  summary_table <- rbind(summary_table, temp_QC)
}

## Keep in only samples which are present in genotype-to-expression file AND additional up to 5000 samples (better phasing)
samples_to_include_gte <- fam[fam$`sample.ID` %in% gte$V1, ]

if (exists("samples_to_include")){
  print(table(samples_to_include_gte$`sample.ID` %in% samples_to_include$`sample.ID`))
  samples_to_include_gte <- samples_to_include_gte[samples_to_include_gte$`sample.ID` %in% samples_to_include$`sample.ID`, ]
  fam <- fam[fam$`sample.ID` %in% samples_to_include$`sample.ID`, ]
}

# Here add up to 5000 samples which are not already included
if (nrow(fam[!fam$`sample.ID` %in% samples_to_include_gte$`sample.ID`, ]) > 0){
  set.seed(123)
  add_samples <- sample(fam[!fam$`sample.ID` %in% samples_to_include_gte$`sample.ID`, ]$`sample.ID`,
                        min(5000, nrow(fam[!fam$`sample.ID` %in% samples_to_include_gte$`sample.ID`, ])))
  set.seed(NULL)
  fam2 <- fam[fam$`sample.ID` %in% add_samples, ]
  samples_to_include_temp <- rbind(samples_to_include_gte, fam2)
} else {samples_to_include_temp <- samples_to_include_gte}

if (exists("samples_to_include") && nrow(samples_to_include) > 0){
  samples_to_include <- samples_to_include[samples_to_include$`sample.ID` %in% samples_to_include_temp$`sample.ID`, ]
  print(nrow(samples_to_include))
} else {samples_to_include <- samples_to_include_temp}

temp_QC <- data.frame(stage = "Samples in genotype-to-expression file + 5000", Nr_of_SNPs = target_bed$ncol,
Nr_of_samples = nrow(samples_to_include),
Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% samples_to_include$`sample.ID`, ]))
summary_table <- rbind(summary_table, temp_QC)

# Remove samples which are in the exclusion list
if (args$exclusion_list != "" & args$exclusion_list != "EmpiricalProbeMatching_AffyU219.txt"){
exc_list <- fread(args$exclusion_list, header = FALSE,
                  keepLeadingZeros = TRUE, colClasses = "character")
samples_to_include <- samples_to_include[!samples_to_include$`sample.ID` %in% exc_list$V1, ]
message("Sample exclusion filter active!")
}

fwrite(data.table(`#FID` = '0', `IID` = samples_to_include$`sample.ID`), "SamplesToInclude.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

temp_QC <- data.frame(stage = "Samples after removing exclusion list", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = nrow(samples_to_include),
Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% samples_to_include$`sample.ID`, ]))
summary_table <- rbind(summary_table, temp_QC)

# Remove samples not in GTE + 5k samples
system(paste0(PLINK2, " --bfile ", bed_simplepath, " --fam fam_normalized.fam",
" --output-chr 26 --keep SamplesToInclude.txt --geno 0.05 --make-bed --threads 4 --out ", bed_simplepath, "_filtered"))

# Do a first pass over variants to remove the bulk of highly missed variants
# List the missingness per variant
#system(paste0(PLINK, " --bfile ",
#              paste0(bed_simplepath, "_filtered"),
#              " --threads 4 --missing --out initial_pass_missingness"))

#initial_pass_missing_variants <- fread("initial_pass_missingness.lmiss") %>%
#  filter(FMISS < 0.05) %>%
#  pull(SNP)

#fwrite(initial_pass_missing_variants, "variants_callrate_95.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Do SNP and sample missingness QC on raw genotype bed
message("Do SNP and genotype QC.")
snp_plinkQC(
  plink.path = PLINK2,
  prefix.in = paste0(bed_simplepath, "_filtered"),
  prefix.out = paste0(bed_simplepath, "_QC"),
  file.type = "--bfile",
  maf = 0.01,
  geno = 0.05,
  mind = 0.05,
  hwe = 1e-6,
  autosome.only = FALSE,
  extra.options = paste0("--output-chr 26 --not-chr 0 25-26 --set-all-var-ids ", variant_format, " --new-id-max-allele-len 20 truncate --threads 4"),
  verbose = TRUE
)

qc_bim <- fread(paste0(bed_simplepath, "_QC.bim"),
                data.table = FALSE, keepLeadingZeros = TRUE)
consecutive_runs <- unlist(lapply(rle(qc_bim[,2])$lengths, seq_len))
consequtive_runs_values <- qc_bim[consecutive_runs != 1, 2]
qc_bim[consecutive_runs != 1, 2] <- paste(qc_bim[consecutive_runs != 1, 2], consecutive_runs[consecutive_runs != 1], sep = "_")

fwrite(qc_bim, paste0(bed_simplepath, "_QC.bim"), sep="\t", row.names=F, col.names=F)

# Read in reference and target genotype data
ref_bed <- bed(paste0(ref_1000g_prefix, ".bed"))
# Read in QCd target genotype data
target_bed <- bed(paste0(bed_simplepath, "_QC.bed"))
target_bed$.fam <- read_fam(paste0(bed_simplepath, "_QC"))

temp_QC <- data.frame(stage = "SNP QC: HWE P>1e-6; MAF>0.01; GENO<0.05; MIND<0.05", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = target_bed$nrow,
Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% target_bed$.fam$`sample.ID`, ]))

summary_table <- rbind(summary_table, temp_QC)

## Assert that all IIDs are unique
if (any(duplicated(target_bed$fam$`sample.ID`))) {
  stop("Individual sample IDs should be unique. Exiting...")
}

sex_check_data_set_chromosomes <- unique(target_bed$map$chromosome)

sex_check_out_path <- paste0(args$output, "/gen_data_QCd/SexCheck.txt")
sex_check_removed_out_path <- paste0(args$output, "/gen_data_QCd/SexCheckFailed.txt")
sex_check_samples <- target_bed$fam

if (23 %in% sex_check_data_set_chromosomes) {

  # Do sex check
  message("Do sex check.")

  # Split x if needed

  pruned_variants_sex_check <- args$pruned_variants_sex_check

  if (!is.null(pruned_variants_sex_check)
    && pruned_variants_sex_check != "") {
    
    if (!file.exists(pruned_variants_sex_check)) {
      stop(sprintf("file '%s' does not exist"))
    }

    message("Using predefined pruned variants for sex-check:")
    message(pruned_variants_sex_check)

    if (ucsc_code != "hg19") {
      variants_sex_check <- fread(
        pruned_variants_sex_check, sep = " ", data.table = FALSE, header = FALSE,
        col.names = c("chr", "pos", "pos.end", "id"))

      variants_sex_check$chr <- "X"

    if (!is.null(args$chain_path) && args$chain_path != "") {
        variants_sex_check_new <- snp_modifyBuild2(
        variants_sex_check, file.path(".", R.utils::getRelativePath(args$liftover_path)),
        from = "hg19", to = ucsc_code, chain_path = chain_path)
    }else{
        variants_sex_check_new <- snp_modifyBuild(
        variants_sex_check, file.path(".", R.utils::getRelativePath(args$liftover_path)),
        from = "hg19", to = ucsc_code)
    }

      variants_sex_check_new$chr <- "23"
      variants_sex_check_new$pos.end <- variants_sex_check_new$pos

      fwrite(variants_sex_check_new[!is.na(variants_sex_check_new$pos),], "mapped_sex_check_variants.txt", col.names=F, row.names=F, quote=F, sep=" ")

      system(paste0(
        PLINK, " --bfile ", bed_simplepath, "_QC", " --extract range mapped_sex_check_variants.txt",
        " --maf 0.05 --make-bed --out ", bed_simplepath, "_split"))

    } else {
 
      system(paste0(
        PLINK, " --bfile ", bed_simplepath, "_QC", " --extract range ", pruned_variants_sex_check,
        " --maf 0.05 --make-bed --out ", bed_simplepath, "_split"))
 
   }
  } else {

    message("Not using predefined pruned variants for sex-check")

    system(paste0(
      PLINK, " --bfile ", bed_simplepath, "_QC",
      " --chr X --maf 0.05 --split-x ", build_code, " no-fail --make-bed --out ", bed_simplepath, "_split"))

  }

  ## Pruning
  system(paste0(PLINK2, " --bfile ", bed_simplepath, "_split",
                " --rm-dup 'exclude-mismatch' --indep-pairwise 20000 200 0.2 --out check_sex_x --threads 4"))

  ## Sex check
  system(paste0(PLINK, " --bfile ", bed_simplepath, "_split --extract check_sex_x.prune.in --check-sex --threads 4"))

  ## If there is sex info in the fam file for all samples then remove samples which fail the sex check or genotype-based F is >0.2 & < 0.8
  sexcheck <- fread("plink.sexcheck", keepLeadingZeros = TRUE,
                    colClasses = list(character = c(1,2)))
  ## Annotate samples who have clear sex

  sexcheck$F_PASS <- !(sexcheck$F > 0.2 & sexcheck$F < 0.8)
  temp_QC <- data.frame(stage = "Sex check (0.2<F<0.8)",
                        Nr_of_SNPs = target_bed$ncol,
                        Nr_of_samples = sum(sexcheck$F_PASS),
                        Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% sexcheck[sexcheck$F_PASS == TRUE, ]$IID, ]))
  summary_table <- rbind(summary_table, temp_QC)

  sexcheck$MATCH_PASS <- case_when(sexcheck$PEDSEX == 0 ~ T,
                                   sexcheck$STATUS == "PROBLEM" ~ F,
                                   TRUE ~ T)

  sexcheck$PASS <- sexcheck$MATCH_PASS & sexcheck$F_PASS

  if (any(sexcheck$PEDSEX %in% c(1, 2))) {

    temp_QC <- data.frame(stage = "Sex check (reported and genetic sex mismatch)",
                          Nr_of_SNPs = target_bed$ncol,
                          Nr_of_samples = sum(sexcheck$PASS),
                          Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% sexcheck[sexcheck$PASS == TRUE, ]$IID, ]))
    summary_table <- rbind(summary_table, temp_QC)

  } else {
    message("No sex info in the .fam file.")
  }

  sex_cols <- c("0" = "black", "1" = "orange", "2" = "blue")

  p <- ggplot(sexcheck, aes(x = F, fill = factor(PEDSEX))) +
    geom_histogram(position="stack", color = "black", alpha = 0.5) +
    scale_fill_manual(values = sex_cols, breaks = c("0", "1", "2"), labels = c("Unknown", "Male", "Female"), name = "Reported sex") +
    geom_vline(xintercept = c(0.2, 0.8), colour = "red", linetype = 2) + theme_bw()

  ggsave(paste0(args$output, "/gen_plots/SexCheck.png"), p, type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
  ggsave(paste0(args$output, "/gen_plots/SexCheck.pdf"), p, height = 7 / 2, width = 9, units = "in", dpi = 300)

  fwrite(sexcheck, sex_check_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  fwrite(sexcheck[!sexcheck$PASS,], sex_check_removed_out_path, sep = "\t", quote = FALSE, row.names = FALSE)

} else {
  warning("No X chromosome present. Skipping sex-check...")

  sexcheck <- sex_check_samples[,c(1,2,5)]
  colnames(sexcheck) <- c("FID", "IID", "PEDSEX")
  sexcheck$PEDSEX_COPY <- sexcheck$PEDSEX
  sexcheck$STATUS <- NA_character_
  sexcheck$F <- NA_real_

  fwrite(sexcheck, sex_check_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  file.create(sex_check_removed_out_path)
}

# Remove sex chromosomes
snp_plinkQC(
  plink.path = PLINK2,
  prefix.in = paste0(bed_simplepath, "_QC"),
  prefix.out = paste0(bed_simplepath, "_QC", "_QC"),
  file.type = "--bfile",
  maf = 0.01,
  geno = 0.05,
  mind = 0.05,
  hwe = 1e-6,
  autosome.only = TRUE,
  extra.options = paste0("--output-chr 26 --remove ", sex_check_removed_out_path, " --threads 4"),
  verbose = TRUE
)

# replace the previous QC version
system(paste0("rm ", bed_simplepath, "_QC.*"))
system(paste0("mv ", bed_simplepath, "_QC_QC.bed ", bed_simplepath, "_QC.bed"))
system(paste0("mv ", bed_simplepath, "_QC_QC.bim ", bed_simplepath, "_QC.bim"))
system(paste0("mv ", bed_simplepath, "_QC_QC.fam ", bed_simplepath, "_QC.fam"))

# Read in again QCd target genotype data
target_bed <- bed(paste0(bed_simplepath, "_QC.bed"))
target_bed$.fam <- read_fam(paste0(bed_simplepath, "_QC"))

temp_QC <- data.frame(stage = "Removed X/Y", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = nrow(target_bed$fam),
Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% target_bed$.fam$`sample.ID`, ]))
summary_table <- rbind(summary_table, temp_QC)

# Do heterozygosity check
message("Do heterozygosity check.")

# Get path where to write heterozygosity failed samples to
het_failed_samples_out_path <- paste0(args$output, "/gen_data_QCd/HeterozygosityFailed.txt")

# Prune variants
system(paste0(PLINK2, " --bfile ", bed_simplepath, "_QC --rm-dup 'exclude-mismatch' --indep-pairwise 50 1 0.2 --threads 4"))

system(paste0(PLINK2, " --bfile ", bed_simplepath, "_QC --extract plink2.prune.in --het --threads 4"))
het <- fread("plink2.het", header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
het$het_rate <- (het$OBS_CT - het$`O(HOM)`) / het$OBS_CT

het_fail_samples <- het[het$het_rate < mean(het$het_rate) - 3 * sd(het$het_rate) | het$het_rate > mean(het$het_rate) + 3 * sd(het$het_rate), ]

print(str(het_fail_samples))

# Get the indices of those samples that passed heterozygozity check
indices_of_het_failed_samples <- match(het_fail_samples$IID, target_bed$fam$`sample.ID`)
indices_of_het_passed_samples <- rows_along(target_bed)
if (length(indices_of_het_failed_samples) > 0) {
  indices_of_het_passed_samples <- rows_along(target_bed)[-indices_of_het_failed_samples]
}

print("het_failed_samples:")
print(indices_of_het_failed_samples)
print("het_passed_samples:")
print(indices_of_het_passed_samples)

fwrite(het_fail_samples, het_failed_samples_out_path, sep = "\t", quote = FALSE, row.names = FALSE)

het_s <- data.frame(ID = target_bed$.fam$`sample.ID`, FAMID = target_bed$.fam$`family.ID`)
het_s <- het_s[!het_s$ID %in% het_fail_samples$IID, ]

temp_QC <- data.frame(stage = "Excess heterozygosity (mean+/-3SD)", Nr_of_SNPs = target_bed$ncol,
Nr_of_samples = length(indices_of_het_passed_samples),
Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% het_s$ID, ]))

summary_table <- rbind(summary_table, temp_QC)

p <- ggplot(het, aes(x = het_rate)) + geom_histogram(color = "#000000", fill = "#000000", alpha = 0.5) +
xlab("Heterozygosity rate") +
geom_vline(xintercept = c(mean(het$het_rate), mean(het$het_rate) + 3 * sd(het$het_rate), mean(het$het_rate) - 3 * sd(het$het_rate)), linetype = 2, colour = "red") +
theme_bw()

ggsave(paste0(args$output, "/gen_plots/HetCheck.png"), type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/HetCheck.pdf"), height = 7 / 2, width = 9, units = "in", dpi = 300)

# Project the data on QCd 1000G reference
message("Projecting samples to 1000G reference.")
unrelated_ref_samples <- fread(args$sample_list, keepLeadingZeros = TRUE, colClasses = 'character')
unrelated_ref_samples <- as.numeric(unrelated_ref_samples$ind.row)

if (ucsc_code != "hg19" && !is.null(args$chain_path) && args$chain_path != "") {

message("Using offline version of PCA sample projection function.")

map_new <- setNames(target_bed$map[-3], c("chr", "rsid", "pos", "a1", "a0"))

map_new_lifted <- snp_modifyBuild2(map_new, 
liftOver = R.utils::getRelativePath(args$liftover_path), 
from = ucsc_code,
to = "hg19", 
chain_path = chain_path)

lifted_bim <- data.table(chr = map_new_lifted$chr,
rsid = map_new_lifted$rsid,
seq = 0,
pos = map_new_lifted$pos,
a1 = map_new_lifted$a1,
a0 = map_new_lifted$a0)

# Drop variants whose ID collides (from --set-all-var-ids truncating long indel alleles to the same
# chr:pos[build]ref,alt string). plink --update-map errors on duplicate IDs in the update file, and
# bed_projectPCA below joins by position anyway, so removing the few colliding variants is safe.
dup_ids <- unique(lifted_bim$rsid[duplicated(lifted_bim$rsid) | duplicated(lifted_bim$rsid, fromLast = TRUE)])
writeLines(dup_ids, "pca_dup_ids.txt")
fwrite(lifted_bim[!is.na(lifted_bim$pos) & !(lifted_bim$rsid %in% dup_ids), ], "lifted_map.bim", sep = "\t", col.names = FALSE, row.names = FALSE)

# Remove the ID-colliding variants in a SEPARATE pass first: plink1.9 validates ID-uniqueness for
# --update-map independently of --exclude, so exclude+update in one command still errors "Duplicate ID".
system(paste0(PLINK,  " --bfile ",  bed_simplepath, "_QC --exclude pca_dup_ids.txt --make-bed --out temp_for_PCA_nodup"))
system(paste0(PLINK,  " --bfile temp_for_PCA_nodup --update-chr lifted_map.bim 1 2 --update-map lifted_map.bim 4 2 --make-bed --out temp_for_PCA"))

proj_PCA <- bed_projectPCA(
  obj.bed.ref = ref_bed,
  ind.row.ref = unrelated_ref_samples,
  obj.bed.new = bed("temp_for_PCA.bed"),
  ind.row.new = indices_of_het_passed_samples,
  k = 10,
  strand_flip = TRUE,
  join_by_pos = TRUE,
  match.min.prop = 0.01,
  build.new = "hg19",
  build.ref = "hg19",
  liftOver = R.utils::getRelativePath(args$liftover_path),
  verbose = TRUE,
  ncores = 4
)

system("rm temp_for_PCA*")

} else {
proj_PCA <- bed_projectPCA(
  obj.bed.ref = ref_bed,
  ind.row.ref = unrelated_ref_samples,
  obj.bed.new = target_bed,
  ind.row.new = indices_of_het_passed_samples,
  k = 10,
  strand_flip = TRUE,
  join_by_pos = TRUE,
  match.min.prop = 0.01,
  build.new = ucsc_code,
  build.ref = "hg19",
  liftOver = R.utils::getRelativePath(args$liftover_path),
  verbose = TRUE,
  ncores = 4
)
}

## Visualise PCs
abi <- as.data.frame(proj_PCA$OADP_proj)
colnames(abi) <- paste0("PC", 1:10)

PCs_ref <- predict(proj_PCA$obj.svd.ref)
abi2 <- as.data.frame(PCs_ref)
colnames(abi2) <- paste0("PC", 1:10)

abi2$sample <- ref_bed$fam$`sample.ID`[unrelated_ref_samples]
abi2 <- abi2[, c(11, 1:10)]

pops <- fread(args$pops, keepLeadingZeros = TRUE, colClasses = list(character = c(2, 6, 7)))
pops <- pops[, c(2, 6, 7)]
abi2 <- merge(abi2, pops, by.x = "sample", by.y = "SampleID")
abi2 <- abi2[, c(1, 12, 13, 2:11)]

abi <- data.frame(sample = target_bed$fam$`sample.ID`[indices_of_het_passed_samples],
                  Population = "Target", Superpopulation = "Target", abi)

abi$type <- "Target"
abi2$type <- "1000G"

combined <- rbind(abi, abi2)

combined$Superpopulation <- factor(combined$Superpopulation, levels = c("Target", "EUR", "EAS", "AMR", "SAS", "AFR"))

p00 <- ggplot(combined, aes(x = PC1, y = PC2, alpha = type)) +
geom_point() +
theme_bw() +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0)) +
ggtitle("Target sample projections\nin 1000G PC space")

combined_h <- combined[combined$Superpopulation == "Target", ]

p0 <- ggplot(combined_h, aes(x = PC1, y = PC2)) +
geom_point() +
theme_bw() +
ggtitle("Target sample projections\nzoomed in")

p1 <- ggplot(combined, aes(x = PC1, y = PC2, colour = Superpopulation, alpha = type)) +
geom_point() +
theme_bw() +
scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p2 <- ggplot(combined, aes(x = PC3, y = PC4, colour = Superpopulation, alpha = type)) +
geom_point() + theme_bw() +
scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p3 <- ggplot(combined, aes(x = PC5, y = PC6, colour = Superpopulation, alpha = type)) +
geom_point() + theme_bw() +
scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p4 <- ggplot(combined, aes(x = PC7, y = PC8, colour = Superpopulation, alpha = type)) +
geom_point() + theme_bw() +
scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p5 <- ggplot(combined, aes(x = PC9, y = PC10, colour = Superpopulation, alpha = type)) +
geom_point() + theme_bw() +
scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p <- p00 + p0 + p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 4)

ggsave(paste0(args$output, "/gen_plots/SamplesPCsProjectedTo1000G.png"), type = "cairo", height = 20, width = 9.5 * 1.6, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/SamplesPCsProjectedTo1000G.pdf"), height = 20, width = 9.5 * 1.6, units = "in", dpi = 300)
fwrite(abi[, -c(2, 3, ncol(abi))], paste0(args$output, "/gen_data_summary/1000G_PC_projections.txt"), sep = "\t", quote = FALSE )

## Assign each sample to the superpopulation
message("Assign each sample to 1000G superpopulation.")
### Calculate distance of each sample to all samples per each population
target_samples <- abi[, -c(2, 3, ncol(abi))]

#### Use 3 PCs
target_samples <- target_samples[, c(1:4)]
rownames(target_samples) <- target_samples$sample
target_samples <- target_samples[, -1]

population_assign_res <- data.frame(sample = rownames(target_samples), abi = rownames(target_samples))

#### EUR
for(population in c("EUR", "EAS", "AMR", "SAS", "AFR")){
abi_e <- abi2[abi2$Superpopulation == population, ]
head(abi_e)

sup_pop_samples <- abi_e[, -c(2, 3, ncol(abi_e))]

sup_pop_samples <- sup_pop_samples[, c(1:4)]
rownames(sup_pop_samples) <- sup_pop_samples$sample
sup_pop_samples <- sup_pop_samples[, -1]

head(sup_pop_samples)

comb <- rbind(target_samples, sup_pop_samples)
head(comb)
distance <- as.matrix(dist(comb, method = "euclidean"))
head(distance)
distance <- distance[c(1:nrow(target_samples)), -c(1:nrow(target_samples))]
head(distance)

head(rowMeans(distance))

distance <- data.frame(sample = rownames(target_samples), MeanDistance = rowMeans(distance))
colnames(distance)[2] <- population

population_assign_res <- cbind(population_assign_res, distance[, -1])

print(paste("distance:", population))

}

colnames(population_assign_res)[3:ncol(population_assign_res)] <- c("EUR", "EAS", "AMR", "SAS", "AFR")
fwrite(population_assign_res[, -1], paste0(args$output, "/gen_data_summary/PopAssignResults.txt"), sep = "\t", quote = FALSE )

# Find related samples
message("Find related samples.")
related <- snp_plinkKINGQC(
  plink2.path = PLINK2,
  bedfile.in = paste0(bed_simplepath, "_QC.bed"),
  thr.king = 2^-4.5,
  make.bed = FALSE,
  ncores = 4,
  extra.options = paste0("--remove ", het_failed_samples_out_path)
)

# Filter in only related individuals from genotype-to-expression file

related <- related[related$IID1 %in% gte$V1 & related$IID2 %in% gte$V1, ]

fwrite(related, "related.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Remove samples that are related to each other
related$IID1 <- as.character(related$IID1)
related$IID2 <- as.character(related$IID2)
# First, get the total list of all samples with some relatedness above a predefined threshold (see above in plink call)
related_individuals <- unique(c(related$IID1, related$IID2))

# If there are related samples, find the samples that should be removed so that the maximum set of samples remains,
# but that also guarantees that no relatedness remains.

# First make an empty vector with samples to remove.
samples_to_remove_due_to_relatedness <- c()

# If there are related individuals, remove these in the following step.
if (length(related_individuals) > 0) {

  # Define a graph wherein each relation depicts an edge between vertices (samples)
  relatedness_graph <- graph_from_edgelist(
    as.matrix(related[,c("IID1", "IID2")]),
    directed = F)

  relatedness_graph <- simplify(
    relatedness_graph,
    remove.multiple = TRUE,
    remove.loops = FALSE,
    edge.attr.comb = igraph_opt("edge.attr.comb")
  )

  # For final version: do not write out, here are original sample IDs
  # pdf(paste0(args$output, "/gen_plots/relatedness.pdf"))
  # plot(relatedness_graph)
  # dev.off()

  # Now, get a list of samples that should be removed due to relatedness
  # We get this through a greedy algorithm trying to find a large possible set of unrelated samples.
  # This is a heuristic solution since the problem is really hard.
  while (length(V(relatedness_graph)) > 1) {

    # Get the degrees (how many edges does each vertex have)
    degrees_named <- degree(relatedness_graph)

    # Get the vertex with the least amount of degrees (edges)
    least_vertex_samples <- names(degrees_named)[min(degrees_named) == degrees_named]

    # Prioritize vertices which are in genotype-to-expression file
    if (length(least_vertex_samples[least_vertex_samples %in% gte$V1]) > 0){
      # if there are multiple related sample IDs from GTE, then take just first
      curr_vertex <- least_vertex_samples[least_vertex_samples %in% gte$V1][1]
    } else {
      # if there are multiple related sample IDs (not in GTEs), then take just first
      curr_vertex <- least_vertex_samples[1]
    }

    # Get all vertices that have an edge with curr_vertex
    related_vertices <- names(relatedness_graph[curr_vertex][relatedness_graph[curr_vertex] > 0])

    # Add these vertexes to the list of vertices to remove
    samples_to_remove_due_to_relatedness <- c(samples_to_remove_due_to_relatedness, related_vertices)

    # Remove the vertices to remove
    relatedness_graph <- delete_vertices(relatedness_graph, c(curr_vertex, related_vertices))
  }

  # Get the indices of those samples that should be removed.
  indices_of_relatedness_failed <- match(
    samples_to_remove_due_to_relatedness,
    target_bed$fam$`sample.ID`)

  # Remove these indices from the indices that remained after the previous check.
  indices_of_passed_samples <- indices_of_het_passed_samples[
    (!indices_of_het_passed_samples %in% indices_of_relatedness_failed)]

  print("relatedness_failed_samples:")
  print(indices_of_relatedness_failed)
  print("relatedness_passed_samples:")
  print(indices_of_passed_samples)

} else {
  # No relatedness observed, proceeding with all samples that passed the previous check.
  indices_of_passed_samples <- indices_of_het_passed_samples
}

temp_QC <- data.frame(stage = "Relatedness for sQTL samples: thr. KING>2^-4.5", Nr_of_SNPs = target_bed$ncol,
Nr_of_samples = length(indices_of_passed_samples),
Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% het_s[!het_s$ID %in% samples_to_remove_due_to_relatedness, ]$ID, ])
)
summary_table <- rbind(summary_table, temp_QC)

### Do PCA on target data
message("Find genetic outliers.")
message("Find genetic outliers: do PCA on QCd target data.")
### PCA
target_pca <- bed_autoSVD(target_bed, ind.row = indices_of_passed_samples, k = 10, ncores = 4)

### Find outlier samples
prob <- bigutilsr::prob_dist(target_pca$u, ncores = 4)
S <- prob$dist.self / sqrt(prob$dist.nn)

# Put threshold for outlier samples, this is by default 0.4!
Sthresh <- args$S_threshold

p <- ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = 0:5 / 5, limits = c(0, NA)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)") +
  geom_vline(aes(xintercept = Sthresh), colour = "red", linetype = 2)

ggsave(paste0(args$output, "/gen_plots/PC_dist_outliers_S.png"), type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/PC_dist_outliers_S.pdf"), height = 7 / 2, width = 9, units = "in", dpi = 300)

# Visualise PCs, outline individual outlier samples
PCs <- predict(target_pca)

PCs <- as.data.frame(PCs)
colnames(PCs) <- paste0("PC", 1:10)
PCs$S <- S

PCs$outlier_ind <- "no"

if (any(PCs$S > Sthresh)) {
  PCs[PCs$S > Sthresh, ]$outlier_ind <- "yes"
}
PCs$sd_outlier <- "no"
sd_outlier_selection <- ((PCs$PC1 > mean(PCs$PC1) + args$SD_threshold * sd(PCs$PC1)
  | PCs$PC1 < mean(PCs$PC1) - args$SD_threshold * sd(PCs$PC1))
  | (PCs$PC2 > mean(PCs$PC2) + args$SD_threshold * sd(PCs$PC2)
  | PCs$PC2 < mean(PCs$PC2) - args$SD_threshold * sd(PCs$PC2)))
if (any(sd_outlier_selection)) {
  PCs[sd_outlier_selection, ]$sd_outlier <- "yes"
}

PCs$outlier <- "no"
if (nrow(PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "no", ]) > 0){
PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "no", ]$outlier <- "S outlier"}
if(nrow(PCs[PCs$outlier_ind == "no" & PCs$sd_outlier == "yes", ]) > 0){
PCs[PCs$outlier_ind == "no" & PCs$sd_outlier == "yes", ]$outlier <- "SD outlier"}
if(nrow(PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "yes", ]) > 0){
PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "yes", ]$outlier <- "S and SD outlier"
}
# For first 2 PCs also remove samples which deviate from the mean

p1 <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick")) +
geom_vline(xintercept = c(mean(PCs$PC1) + 3 * sd(PCs$PC1), mean(PCs$PC1) - 3 * sd(PCs$PC1)), colour = "firebrick", linetype = 2) +
geom_hline(yintercept = c(mean(PCs$PC2) + 3 * sd(PCs$PC2), mean(PCs$PC2) - 3 * sd(PCs$PC2)), colour = "firebrick", linetype = 2)
p2 <- ggplot(PCs, aes(x = PC3, y = PC4, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p3 <- ggplot(PCs, aes(x = PC5, y = PC6, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p4 <- ggplot(PCs, aes(x = PC7, y = PC8, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p5 <- ggplot(PCs, aes(x = PC9, y = PC10, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))

p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

ggsave(paste0(args$output, "/gen_plots/PCA_outliers.png"), type = "cairo", height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/PCA_outliers.pdf"), height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)

# Filter out related samples and outlier samples, write out QCd data
message("Filter out related samples and outlier samples, write out QCd data.")
indices_of_passed_samples <- indices_of_passed_samples[PCs$outlier == "no"]
samples_to_include <- data.frame(family.ID = target_bed$.fam$`family.ID`[indices_of_passed_samples], sample.IDD2 = target_bed$.fam$sample.ID[indices_of_passed_samples])

temp_QC <- data.frame(stage = paste0("Outlier samples: thr. S>", Sthresh, " PC1/PC2 SD deviation thresh ", args$SD_threshold), Nr_of_SNPs = target_bed$ncol,
Nr_of_samples = nrow(samples_to_include),
Nr_of_eQTL_samples = nrow(gte[gte$V1 %in% samples_to_include$`sample.IDD2`, ]))
summary_table <- rbind(summary_table, temp_QC)

fwrite(data.table::data.table(samples_to_include), "SamplesToInclude.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# Remove samples
system(paste0(PLINK2, " --bfile ", bed_simplepath, "_QC --output-chr 26 --keep SamplesToInclude.txt --make-bed --threads 4 --out ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation"))

# Reorder the samples and write out the sample file
message("Shuffle sample order.")
rows <- sample(nrow(samples_to_include))
samples_to_include2 <- samples_to_include[rows, ]
fwrite(samples_to_include2, "ShuffledSampleOrder.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

system(paste0(PLINK2, " -bfile ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation ",
"--indiv-sort f ShuffledSampleOrder.txt ",
"--make-bed ",
"--out ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation_temp --threads 4"))

# Do final SNP QC (for MAF, etc filters on filtered SNPs)
# Remove unfiltered samples
system(paste0("rm ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.*"))

message("Final SNP QC.")

snp_plinkQC(
  plink.path = PLINK2,
  prefix.in = paste0(args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation_temp"),
  prefix.out = paste0(args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation"),
  file.type = "--bfile",
  maf = 0.01,
  geno = 0.05,
  mind = 0.05,
  hwe = 1e-6,
  autosome.only = TRUE,
  extra.options = "--output-chr 26 --threads 4",
  verbose = TRUE
)

system(paste0("rm ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation_temp*"))

# Final rerun PCA on QCd data
message("Final PCA on QCd data.")
bed_qc <- bed(paste0(args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.bed"))
bed_qc$.fam <- read_fam(paste0(args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation"))

target_pca_qcd <- bed_autoSVD(bed_qc, k = 10, ncores = 4)

# Visualise loadings
message("Plot PC post QC loadings.")
plot(target_pca_qcd, type = "loadings", loadings = 1:10, coeff = 0.6)
ggsave(paste0(args$output, "/gen_plots/Target_PCs_postQC_Loadings.png"), type = "cairo", height = (5 * 7) * 0.7, width = (5 * 7) * 0.7, units = "in", dpi = 300)

PCsQ <- predict(target_pca_qcd)
PCsQ <- as.data.frame(PCsQ)

colnames(PCsQ) <- paste0("PC", 1:10)
rownames(PCsQ) <- bed_qc$fam$sample.ID

# Visualise
p1 <- ggplot(PCsQ, aes(x = PC1, y = PC2)) + theme_bw() + geom_point(alpha = 0.5)
p2 <- ggplot(PCsQ, aes(x = PC3, y = PC4)) + theme_bw() + geom_point(alpha = 0.5)
p3 <- ggplot(PCsQ, aes(x = PC5, y = PC6)) + theme_bw() + geom_point(alpha = 0.5)
p4 <- ggplot(PCsQ, aes(x = PC7, y = PC8)) + theme_bw() + geom_point(alpha = 0.5)
p5 <- ggplot(PCsQ, aes(x = PC9, y = PC10)) + theme_bw() + geom_point(alpha = 0.5)
p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

ggsave(paste0(args$output, "/gen_plots/Target_PCs_postQC.png"), type = "cairo", height = 10 * 1.5, width = 9 * 1.3, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/Target_PCs_postQC.pdf"), height = 10 * 1.5, width = 9 * 1.3, units = "in", dpi = 300)

# Write out
fwrite(PCsQ, paste0(args$output, "/gen_PCs/GenotypePCs.txt"), row.names = TRUE, sep = "\t", quote = FALSE)

# Write out scree plots

#p <- plot(target_pca_qcd)

singlar_value <- data.frame(
  PC = paste0("PC", 1:10), 
  sv = target_pca_qcd$d
  )
singlar_value$PC <- factor(singlar_value$PC, levels = as.character(singlar_value$PC))
message("Plot scree plot.")

p <- ggplot(singlar_value, aes(x = PC, y = sv)) + 
geom_bar(stat = "identity") + 
theme_bw() + 
ylab("Singular value")

ggsave(paste0(args$output, "/gen_plots/Target_PCs_scree_postQC.png"), type = "cairo", height = 5, width = 9, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/Target_PCs_scree_postQC.pdf"), height = 5, width = 9, units = "in", dpi = 300)


# Count samples in overlapping with GTE
final_samples <- fread(paste0(args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.fam"), header = FALSE,
                       keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))

temp_QC <- data.frame(stage = "QCd samples overlapping with genotype-to-expression file and SNP QC filters on full dataset",
Nr_of_SNPs = bed_qc$ncol,
Nr_of_samples = nrow(final_samples),
Nr_of_eQTL_samples = nrow(final_samples[final_samples$V2 %in% gte$V1, ]))
summary_table <- rbind(summary_table, temp_QC)

# Write out final summary
message("Write out final sample summary table.")
colnames(summary_table) <- c("Stage", "Nr. of SNPs", "Nr. of genotype samples", "Nr. of sQTL samples")
fwrite(summary_table, paste0(args$output, "/gen_data_summary/summary_table.txt"), sep = "\t", quote = FALSE)

system("rm *.bed", wait = TRUE, intern = FALSE)
system("rm *.bim", wait = TRUE, intern = FALSE)
system("rm *.fam", wait = TRUE, intern = FALSE)
system("rm *.id", wait = TRUE, intern = FALSE)
system("rm *.log", wait = TRUE, intern = FALSE)
system("rm *.hh", wait = TRUE, intern = FALSE)
message("Temporary genotype files cleaned up.")
