#!/usr/bin/env Rscript

# ProcessSplicing.R
# Splicing QC with iterative PCA/MDS outlier detection
# Similar to eQTLGen ProcessExpression.R but adapted for splicing PSI data

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)  # For combining plots
  library(MASS)  # For isoMDS
  library(edgeR) # For data handling
})

# Command line arguments
option_list <- list(
  make_option(c("--psi_matrix"), type = "character", default = NULL,
    help = "Path to PSI matrix (junctions in rows, samples in columns)"),
  make_option(c("--gte"), type = "character", default = NULL,
    help = "Genotype-to-expression linking file (GenotypeID, RnaID)"),
  make_option(c("--geno_filter"), type = "character", default = NULL,
    help = "File with genotype-filtered sample IDs"),
  make_option(c("--sd"), type = "double", default = 4,
    help = "Standard deviation threshold for removing splicing outliers. Default: 4 SDs from mean of PC1/PC2"),
  make_option(c("--output_dir"), type = "character", default = ".",
    help = "Output directory base (default: current directory)"),
  make_option(c("--int_transform"), type = "logical", default = TRUE,
    help = "Apply per-phenotype (per-junction) rank inverse-normal transform on the QC'd matrix, matching eQTLGen ProcessExpression.R (forces each junction to N(0,1) across samples). Default: TRUE. Pass --int_transform FALSE to keep the upstream per-sample qqnorm only."),
  make_option(c("--n_top_features"), type = "integer", default = 20000,
    help = "TOP-FEATURES variant: compute the splicing PCA/MDS (sample-outlier detection + covariate PCs) on only the N most-variable junctions (by raw per-junction variance), so standard prcomp fits in modest memory and NO Gram/dual PCA is needed. The FULL junction matrix is still written out as the phenotype. Default: 20000. Set <=0 to use all junctions (reverts to full prcomp)."),
  make_option(c("--site_summary_dir"), type = "character", default = NULL,
    help = "PRE-FILTER (optional): directory holding site_readcount_summary_{junction,IR}.tsv (from compute_site_readcount_summary.py). When supplied, the TOP-FEATURES PCA basis is restricted to features that are well-measured across donors (see --max_pca_missingness / --min_pca_meancount) BEFORE ranking by variance. Only the PCA/MDS basis is affected; the FULL matrix is still written out. If omitted/missing, the pre-filter is skipped (rank by variance over all QC'd features)."),
  make_option(c("--max_pca_missingness"), type = "double", default = 0.10,
    help = "PRE-FILTER: keep a feature as a PCA-basis candidate only if its per-site missingness (fraction of donors with zero raw reads, from --site_summary_dir) is < this. Default: 0.10."),
  make_option(c("--min_pca_meancount"), type = "double", default = 5,
    help = "PRE-FILTER: keep a feature as a PCA-basis candidate only if its mean raw read depth per donor (mean_count, from --site_summary_dir) is >= this. Default: 5.")
)

args <- parse_args(OptionParser(option_list = option_list))

# Per-phenotype rank inverse-normal transform (identical to eQTLGen DataQC ProcessExpression.R).
# Applied per ROW (junction) across samples; NA-safe.
INT_transform <- function(x){
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# ---------------------------------------------------------------------------
# Memory-lean chunked row operations.
# apply(mat, 1, FUN) on the FULL ~673k-junction matrix materialises an internal
# row-list of the WHOLE matrix (several full copies -> ~70-100GB spikes). These
# helpers process ROW_CHUNK rows at a time, so the row-list is only ever built
# for a small block. Each row is independent -> results are BIT-IDENTICAL to the
# unchunked apply(); only peak memory changes.
# ---------------------------------------------------------------------------
ROW_CHUNK <- 20000L

# per-row scalar FUN (var/mean/median/min/max/sd/counts...) -> length-nrow(mat) vector
chunked_rowfun_vec <- function(mat, FUN, ...) {
  n <- nrow(mat); res <- numeric(n)
  for (s in seq.int(1L, n, by = ROW_CHUNK)) {
    idx <- s:min(s + ROW_CHUNK - 1L, n)
    res[idx] <- apply(mat[idx, , drop = FALSE], 1, FUN, ...)
  }
  names(res) <- rownames(mat)
  res
}

# per-row FUN returning one value per column (e.g. INT) -> nrow(mat) x ncol(mat) matrix
chunked_row_mat <- function(mat, FUN, ...) {
  n <- nrow(mat)
  out <- matrix(NA_real_, n, ncol(mat), dimnames = dimnames(mat))
  for (s in seq.int(1L, n, by = ROW_CHUNK)) {
    idx <- s:min(s + ROW_CHUNK - 1L, n)
    out[idx, ] <- t(apply(mat[idx, , drop = FALSE], 1, FUN, ...))
  }
  out
}

# ---------------------------------------------------------------------------
# PRE-FILTER helper: read the per-site raw read-count summary (junction & IR) and return, keyed by the
# combined-matrix feature id (junction_<site_id> / IR_<site_id>), the per-feature missingness and mean
# read depth. Used to restrict the TOP-FEATURES PCA basis to well-measured features. Small TSVs (one row
# per analysed site), so plain fread is cheap. Returns NULL if no usable summary is found.
# ---------------------------------------------------------------------------
read_site_quality <- function(summary_dir) {
  files <- list(junction = file.path(summary_dir, "site_readcount_summary_junction.tsv"),
                IR       = file.path(summary_dir, "site_readcount_summary_IR.tsv"))
  miss <- numeric(0); mc <- numeric(0)
  for (mod in names(files)) {
    f <- files[[mod]]
    if (!file.exists(f)) next
    df <- tryCatch(fread(f, data.table = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    if (!all(c("site_id", "missingness", "mean_count") %in% colnames(df))) next
    feat <- paste0(mod, "_", df$site_id)            # reconstruct combined-matrix feature id
    m <- as.numeric(df$missingness); names(m) <- feat
    d <- as.numeric(df$mean_count);  names(d) <- feat
    miss <- c(miss, m); mc <- c(mc, d)
  }
  if (length(miss) == 0) return(NULL)
  list(miss = miss, mc = mc)
}

# Validate inputs
if (is.null(args$psi_matrix)) {
  stop("--psi_matrix is required!")
}
if (is.null(args$gte)) {
  stop("--gte is required!")
}

message("===== Splicing QC with PCA/MDS Outlier Detection =====")
message(paste("PSI matrix:", args$psi_matrix))
message(paste("GTE file:", args$gte))
message(paste("SD threshold:", args$sd))
message(paste("Output directory:", args$output_dir))

# Create output directories
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$output_dir, "splice_data_QCd"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$output_dir, "splice_data_summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$output_dir, "splice_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$output_dir, "splice_PCs"), recursive = TRUE, showWarnings = FALSE)

# Read PSI matrix
message("Reading PSI matrix...")
# Handle both .gz and uncompressed files
if (grepl("\\.gz$", args$psi_matrix)) {
  psi <- fread(cmd = paste("gunzip -c", args$psi_matrix), data.table = FALSE)
} else {
  psi <- fread(args$psi_matrix, data.table = FALSE)
}
rownames(psi) <- psi[, 1]
psi <- psi[, -1, drop = FALSE]

# Drop annotation columns carried by the combined matrix (gene, strand); keep only sample columns.
psi <- psi[, !(colnames(psi) %in% c("gene", "strand")), drop = FALSE]

message(paste("Initial PSI dimensions:", nrow(psi), "junctions ×", ncol(psi), "samples"))

# Read GTE file
message("Reading GTE file...")
gte <- fread(args$gte, header = FALSE, data.table = TRUE)
colnames(gte) <- c("V1", "V2")  # V1 = GenotypeID, V2 = RnaID

# Read genotype filter if provided.
# geno_filter is a plink .fam (gen_data_QCd/*_ToImputation.fam). With --const-fid the FID (col V1)
# is a constant "0" and the genotype sample ID lives in the IID column (V2). Match gte genotype IDs
# (V1) against the .fam IID (V2), exactly like eQTLGen ProcessExpression.R (gte$V1 %in% geno_fam$IID).
if (!is.null(args$geno_filter) && file.exists(args$geno_filter)) {
  message("Reading genotype filter...")
  geno_filter <- fread(args$geno_filter, header = FALSE, data.table = FALSE)
  gte <- gte[gte$V1 %in% geno_filter$V2, ]
  message(paste("Filtered to", nrow(gte), "genotyped samples"))
  if (nrow(gte) == 0) {
    stop("No genotyped samples after geno_filter: gte column 1 (genotype IDs) did not match the ",
         ".fam IID column (V2). Check that gte genotype IDs correspond to the .fam IID, not FID.")
  }
}

# Filter PSI matrix to samples in GTE
psi_samples <- colnames(psi)
gte_rna_samples <- gte$V2
common_samples <- intersect(psi_samples, gte_rna_samples)

message(paste("Samples in PSI matrix:", length(psi_samples)))
message(paste("Samples in GTE file:", length(gte_rna_samples)))
message(paste("Common samples:", length(common_samples)))

psi_filtered <- psi[, common_samples, drop = FALSE]

# Summary table (mirrors ProcessExpression.R funnel)
summary_table <- data.table(
  Stage = "Unprocessed matrix",
  Nr_of_junctions = nrow(psi),
  Nr_of_samples = ncol(psi)
)
summary_table <- rbind(summary_table, data.table(
  Stage = "Samples which overlap with QCd genotype info",
  Nr_of_junctions = nrow(psi_filtered),
  Nr_of_samples = ncol(psi_filtered)
))

# Free the original full matrix (all input samples) now that psi_filtered (genotyped samples) is built;
# for a 673k-junction cohort this reclaims ~25-50GB of data.frame before the memory-heavy steps below.
rm(psi); gc()

# ===================================================================
# Junction QC: Missing rate and variance filtering
# ===================================================================
message("\n===== Junction QC: Missing rate and variance filtering =====")

# Calculate missing rate per junction (chunked: avoids a full 12GB is.na() logical copy)
missing_rate <- chunked_rowfun_vec(psi_filtered, function(x) mean(is.na(x)))
message(paste("Mean missing rate across junctions:", round(mean(missing_rate), 3)))

# Filter by missing rate (essentially no filter: allow up to 99% missing)
# Set to very permissive threshold to effectively disable filtering
missing_threshold <- 0.99
keep_missing <- missing_rate < missing_threshold
psi_filtered <- psi_filtered[keep_missing, , drop = FALSE]

message(paste("After missing rate filter (<", missing_threshold, "):", nrow(psi_filtered), "junctions"))

# Calculate variance per junction (chunked: avoids full-matrix row-list ~70GB spike)
junction_var <- chunked_rowfun_vec(psi_filtered, var, na.rm = TRUE)
message(paste("Mean variance across junctions:", round(mean(junction_var, na.rm = TRUE), 4)))

# Filter low-variance junctions (essentially no filter: keep everything with var > 0)
# Set to minimum threshold to effectively disable filtering
var_threshold <- 0.0
keep_var <- junction_var > var_threshold | is.na(junction_var)
psi_filtered <- psi_filtered[keep_var, , drop = FALSE]

message(paste("After variance filter (> 0):", nrow(psi_filtered), "junctions"))

# ===================================================================
# TOP-FEATURES: pick the N most-variable junctions for the PCA/MDS steps ONLY.
# Rationale: standard prcomp on the full ~673k-junction matrix needs ~hundreds of GB (hence the
# Gram/dual PCA in 1_DataQC_GRAM_PCA). Computing sample-outlier detection + covariate PCs on
# the top-N most-variable junctions (~like expression's ~20k genes) keeps standard prcomp at ~tens
# of GB and needs NO Gram PCA. The FULL junction matrix is still written out as the phenotype;
# only the PCA/MDS use this subset. Selection is on RAW per-junction variance (pre-INT), because the
# later per-junction INT forces every junction to unit variance (making "most variable" meaningless).
#
# PRE-FILTER (--site_summary_dir): raw-variance ranking alone over-selects high-missingness / low-depth /
# IR-artifact features. When a read-count summary is supplied, candidates are first restricted to
# well-measured features (missingness < max_pca_missingness AND mean_count >= min_pca_meancount) and the
# top-N most-variable are taken FROM that set. Only the PCA/MDS basis changes; the full matrix is unchanged.
# ===================================================================
n_top <- args$n_top_features

cand <- rownames(psi_filtered)                 # QC'd candidates (missing<0.99 & var>0)
jvar_cand <- junction_var[keep_var]            # per-candidate variance, aligned to psi_filtered rows
names(jvar_cand) <- cand

# Apply the splicing-PCA pre-filter if a read-count summary was supplied.
pass <- rep(TRUE, length(cand)); names(pass) <- cand
prefilter_applied <- FALSE
if (!is.null(args$site_summary_dir) && dir.exists(args$site_summary_dir)) {
  qc <- read_site_quality(args$site_summary_dir)
  if (!is.null(qc)) {
    matched <- cand %in% names(qc$miss)
    message(sprintf("PRE-FILTER: matched %d/%d (%.1f%%) candidate features to the read-count summary.",
                    sum(matched), length(cand), 100 * mean(matched)))
    if (mean(matched) >= 0.5) {
      m <- qc$miss[cand]; d <- qc$mc[cand]
      pass <- matched & !is.na(m) & !is.na(d) &
              (m < args$max_pca_missingness) & (d >= args$min_pca_meancount)
      names(pass) <- cand
      prefilter_applied <- TRUE
      message(sprintf("PRE-FILTER: %d/%d features pass missingness<%.2f & mean_count>=%g.",
                      sum(pass), length(cand), args$max_pca_missingness, args$min_pca_meancount))
    } else {
      message("PRE-FILTER: summary coverage <50% of candidates -> SKIPPING pre-filter ",
              "(check --site_summary_dir / feature-id formats).")
    }
  } else {
    message("PRE-FILTER: no usable read-count summary in --site_summary_dir -> skipping pre-filter.")
  }
} else {
  message("PRE-FILTER: no --site_summary_dir -> skipping pre-filter (rank by variance over all QC'd features).")
}

cand_pass <- cand[pass]
if (sum(pass) == 0) {                           # fail-safe: never empty the PCA basis
  message("PRE-FILTER: WARNING no feature passed the pre-filter -> falling back to ALL QC'd features.")
  cand_pass <- cand
  prefilter_applied <- FALSE
}

if (!is.null(n_top) && n_top > 0 && n_top < length(cand_pass)) {
  ord <- order(jvar_cand[cand_pass], decreasing = TRUE)
  top_junctions <- cand_pass[ord[seq_len(n_top)]]
  message(sprintf("TOP-FEATURES: PCA/MDS on the %d most-variable %sfeatures (of %d candidates%s); full matrix still written out.",
                  n_top, if (prefilter_applied) "pre-filtered " else "", length(cand_pass),
                  if (prefilter_applied) sprintf(", pre-filtered from %d", length(cand)) else ""))
} else {
  top_junctions <- cand_pass
  message(sprintf("TOP-FEATURES: n_top_features=%s -> using ALL %d %scandidate features for PCA/MDS (full prcomp).",
                  as.character(n_top), length(cand_pass), if (prefilter_applied) "pre-filtered " else ""))
}
psi_pca <- psi_filtered[top_junctions, , drop = FALSE]   # subset fed to PCA/MDS

# Persist the top-feature id list so a downstream finalize (FINALIZE_SPLICING, when expression QC is active)
# can recompute the splicing covariate PCs on the final geno n splice n expr sample set from the same basis.
writeLines(top_junctions, file.path(args$output_dir, "splice_PCs", "top_junctions.txt"))

# ===================================================================
# Iterative PCA-based outlier detection
# ===================================================================
message("\n===== Iterative PCA-based outlier detection =====")

IterativeOutlierDetection <- function(
  psi_mat, sd_threshold = 4
) {

  outliers <- c()
  it_round <- 0
  nr_outliers <- 1
  list_ggplots <- list()
  plot_it_round <- 0

  message(paste("Removing samples which deviate more than", sd_threshold, "SDs from the mean of PC1 or PC2."))
  message("Starting iterative outlier detection...")

  # Prepare data for PCA: transpose so samples are rows
  psi_t <- t(psi_mat)

  while (nr_outliers > 0) {
    it_round <- it_round + 1
    plot_it_round <- plot_it_round + 1

    message(paste("\n--- Iteration", it_round, "---"))

    # Handle missing values: replace NA with row mean
    psi_imputed <- psi_t
    for (i in 1:ncol(psi_imputed)) {
      col_mean <- mean(psi_imputed[, i], na.rm = TRUE)
      psi_imputed[is.na(psi_imputed[, i]), i] <- col_mean
    }

    # Scale the data
    psi_scaled <- scale(psi_imputed)

    # PCA
    message("  Calculating PCA...")
    pca_res <- prcomp(psi_scaled, center = FALSE, scale. = FALSE)

    # Extract all PCs dynamically (following eQTLGen pattern)
    # Number of PCs = min(n_samples - 1, n_junctions)
    num_pcs <- ncol(pca_res$x)
    message(sprintf("  Extracted %d PCs", num_pcs))

    PCs <- as.data.frame(pca_res$x, stringsAsFactors = FALSE)
    PCs$sample <- rownames(psi_t)
    PCs <- PCs[, c("sample", setdiff(colnames(PCs), "sample"))]

    message("  PCA calculation finished!")

    # Calculate means and SDs
    PCs1_mean <- mean(PCs$PC1)
    PCs1_sd <- sd(PCs$PC1)
    PCs2_mean <- mean(PCs$PC2)
    PCs2_sd <- sd(PCs$PC2)

    # Identify outliers based on PC1 and PC2
    PCs$outlier <- "no"

    # PC1 outliers - use vector indexing instead of dataframe subsetting
    pc1_outliers <- (PCs$PC1 > PCs1_mean + sd_threshold * PCs1_sd) |
                    (PCs$PC1 < PCs1_mean - sd_threshold * PCs1_sd)
    PCs$outlier[pc1_outliers] <- "yes"

    # PC2 outliers - use vector indexing instead of dataframe subsetting
    pc2_outliers <- (PCs$PC2 > PCs2_mean + sd_threshold * PCs2_sd) |
                    (PCs$PC2 < PCs2_mean - sd_threshold * PCs2_sd)
    PCs$outlier[pc2_outliers] <- "yes"

    # Create plots for PC1-8 (following eQTLGen DataQC style: 4 plots showing PC1 vs PC2, PC3 vs PC4, PC5 vs PC6, PC7 vs PC8)
    p1 <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    list_ggplots[[plot_it_round]] <- p1

    plot_it_round <- plot_it_round + 1
    p2 <- ggplot(PCs, aes(x = PC3, y = PC4, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    list_ggplots[[plot_it_round]] <- p2

    plot_it_round <- plot_it_round + 1
    p3 <- ggplot(PCs, aes(x = PC5, y = PC6, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    list_ggplots[[plot_it_round]] <- p3

    plot_it_round <- plot_it_round + 1
    p4 <- ggplot(PCs, aes(x = PC7, y = PC8, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    list_ggplots[[plot_it_round]] <- p4

    # Get non-outliers
    non_outliers <- PCs[PCs$outlier == "no", ]$sample

    # Remove outlier samples
    psi_t <- psi_t[non_outliers, , drop = FALSE]

    nr_outliers <- sum(PCs$outlier == "yes")

    if (nr_outliers > 0) {
      message(paste0("  Iteration ", it_round, ": Removed ", nr_outliers,
                     " outlier(s). Re-running PCA..."))
      outliers <- c(outliers, PCs[PCs$outlier == "yes", ]$sample)
    } else {
      message(paste0("  Iteration ", it_round, ": No outliers detected. Finalizing."))
    }

    # Safety: max 10 iterations
    if (it_round >= 10) {
      message("  Reached maximum 10 iterations. Stopping.")
      break
    }
  }

  # Final PCA on clean data
  message("\nFinal PCA on outlier-free data...")
  psi_imputed <- psi_t
  for (i in 1:ncol(psi_imputed)) {
    col_mean <- mean(psi_imputed[, i], na.rm = TRUE)
    psi_imputed[is.na(psi_imputed[, i]), i] <- col_mean
  }
  psi_scaled <- scale(psi_imputed)
  pca_final <- prcomp(psi_scaled, center = FALSE, scale. = FALSE)

  # Extract all PCs dynamically
  num_pcs_final <- ncol(pca_final$x)
  message(sprintf("Final PCA: extracted %d PCs", num_pcs_final))

  PCs_final <- as.data.frame(pca_final$x, stringsAsFactors = FALSE)
  PCs_final$sample <- rownames(psi_t)
  PCs_final$outlier <- "no"
  PCs_final <- PCs_final[, c("sample", setdiff(colnames(PCs_final), c("sample", "outlier")), "outlier")]

  # Variance explained (up to 10 PCs or all if fewer)
  num_var_pcs <- min(10, num_pcs_final)
  variance_explained <- summary(pca_final)$importance[2, 1:num_var_pcs]

  return(list(
    psi_mat_clean = t(psi_t),
    outliers = outliers,
    PCs_final = PCs_final,
    pca_object = pca_final,
    variance_explained = variance_explained,
    plots = list_ggplots,
    total_iterations = it_round
  ))
}

# Run iterative outlier detection
# PCA/MDS run on the top-feature subset (psi_pca); only SAMPLES are removed, never junctions.
iterative_result <- IterativeOutlierDetection(psi_pca, sd_threshold = args$sd)

message(paste("\nTotal iterations:", iterative_result$total_iterations))
message(paste("Total outliers removed:", length(iterative_result$outliers)))
message(paste("Remaining samples:", ncol(iterative_result$psi_mat_clean)))

summary_table <- rbind(summary_table, data.table(
  Stage = "After removal of all splicing outliers (PCA)",
  Nr_of_junctions = nrow(psi_filtered),                       # full junction count (outlier detection removes samples, not junctions)
  Nr_of_samples = ncol(iterative_result$psi_mat_clean)
))

# Save PCA iteration plots
message("\nSaving PCA iteration plots...")

# Following eQTLGen pattern: combine 4 plots (PC1-8) into 2x2 grid for each iteration
# Each iteration now has 4 plots (PC1vs2, PC3vs4, PC5vs6, PC7vs8)
num_iterations <- length(iterative_result$plots) / 4

message(sprintf("  Creating combined 2x2 grid plots for %d iteration(s)...", num_iterations))

# Save first iteration as "PCA_before.png" (following eQTLGen naming)
if (num_iterations >= 1) {
  # Get first 4 plots (iteration 1)
  p_combined_before <- (iterative_result$plots[[1]] | iterative_result$plots[[2]]) /
                       (iterative_result$plots[[3]] | iterative_result$plots[[4]])

  ggsave(file.path(args$output_dir, "splice_plots", "PCA_before.png"),
         p_combined_before, height = 9, width = 10.5, units = "in", dpi = 300)
  ggsave(file.path(args$output_dir, "splice_plots", "PCA_before.pdf"),
         p_combined_before, height = 9, width = 10.5, units = "in", dpi = 300)
}

# Save subsequent iterations as iteration-numbered combined plots (if any)
if (num_iterations > 1) {
  for (iter in 2:num_iterations) {
    plot_idx_start <- (iter - 1) * 4 + 1
    p_combined <- (iterative_result$plots[[plot_idx_start]] | iterative_result$plots[[plot_idx_start + 1]]) /
                  (iterative_result$plots[[plot_idx_start + 2]] | iterative_result$plots[[plot_idx_start + 3]])

    png_filename <- sprintf("PCA_iteration_%02d.png", iter)
    pdf_filename <- sprintf("PCA_iteration_%02d.pdf", iter)

    ggsave(file.path(args$output_dir, "splice_plots", png_filename),
           p_combined, height = 9, width = 10.5, units = "in", dpi = 300)
    ggsave(file.path(args$output_dir, "splice_plots", pdf_filename),
           p_combined, height = 9, width = 10.5, units = "in", dpi = 300)
  }
}

# ===================================================================
# MDS-based outlier detection
# ===================================================================
message("\n===== MDS-based outlier detection =====")

# Prepare data for MDS
psi_clean <- iterative_result$psi_mat_clean
psi_clean_t <- t(psi_clean)

# Impute missing values
psi_imputed <- psi_clean_t
for (i in 1:ncol(psi_imputed)) {
  col_mean <- mean(psi_imputed[, i], na.rm = TRUE)
  psi_imputed[is.na(psi_imputed[, i]), i] <- col_mean
}

# Calculate correlation distance
message("  Calculating correlation distance...")
cor_mat <- cor(t(psi_imputed))
dist_mat <- 1 - cor_mat

# Run MDS
message("  Running isoMDS...")
mds_res <- isoMDS(dist_mat, k = 2)

MDS <- data.frame(
  sample = rownames(psi_clean_t),
  MDS1 = mds_res$points[, 1],
  MDS2 = mds_res$points[, 2],
  stringsAsFactors = FALSE
)

# Identify outliers based on MDS
MDS1_mean <- mean(MDS$MDS1)
MDS1_sd <- sd(MDS$MDS1)
MDS2_mean <- mean(MDS$MDS2)
MDS2_sd <- sd(MDS$MDS2)

MDS$outlier <- "no"

# MDS1 outliers - use vector indexing instead of dataframe subsetting
mds1_outliers <- (MDS$MDS1 > MDS1_mean + args$sd * MDS1_sd) |
                 (MDS$MDS1 < MDS1_mean - args$sd * MDS1_sd)
MDS$outlier[mds1_outliers] <- "yes"

# MDS2 outliers - use vector indexing instead of dataframe subsetting
mds2_outliers <- (MDS$MDS2 > MDS2_mean + args$sd * MDS2_sd) |
                 (MDS$MDS2 < MDS2_mean - args$sd * MDS2_sd)
MDS$outlier[mds2_outliers] <- "yes"

mds_outliers <- MDS[MDS$outlier == "yes", ]$sample
message(paste("  MDS outliers detected:", length(mds_outliers)))

# Remove MDS outliers. NOTE (top-features): psi_clean is the top-N PCA subset, used only to
# determine the clean SAMPLE set; psi_final is rebuilt from the FULL matrix (psi_filtered) so the
# written phenotype keeps all junctions.
if (length(mds_outliers) > 0) {
  keep_samples <- setdiff(colnames(psi_clean), mds_outliers)
  psi_final <- psi_filtered[, keep_samples, drop = FALSE]
} else {
  psi_final <- psi_filtered[, colnames(psi_clean), drop = FALSE]
  message("  No MDS outliers detected.")
}

# Always report the MDS step (mirrors ProcessExpression.R), even when no outliers were removed
summary_table <- rbind(summary_table, data.table(
  Stage = "After removal of all splicing outliers (MDS)",
  Nr_of_junctions = nrow(psi_final),
  Nr_of_samples = ncol(psi_final)
))

# Plot MDS (following eQTLGen DataQC style)
p_mds <- ggplot(MDS, aes(x = MDS1, y = MDS2, colour = outlier)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  scale_colour_manual(values = c("no" = "black", "yes" = "red"))

# Save both PNG and PDF (following eQTLGen pattern)
ggsave(file.path(args$output_dir, "splice_plots", "MDS_outliers.png"),
       p_mds, width = 8, height = 6)
ggsave(file.path(args$output_dir, "splice_plots", "MDS_outliers.pdf"),
       p_mds, width = 8, height = 6)

# ===================================================================
# Per-phenotype inverse normal transform (eQTLGen ProcessExpression.R-style)
# Applied to the final QC'd matrix BEFORE the final PCA, so the splicing PCs, the written
# splice_data_preprocessed.txt, and the per-junction summary are all on the INT'd phenotypes.
# Forces each junction to N(0,1) across samples (the per-gene step that ProcessExpression.R does
# but the upstream LeafCutter-style per-sample qqnorm does not). NOTE: this Gaussianizes every
# junction including near-constant ones, so combine with a usage/variance filter (var_threshold)
# to avoid fabricating signal from degenerate phenotypes.
# ===================================================================
if (isTRUE(args$int_transform)) {
  message("\n===== Per-phenotype inverse normal transform (eQTLGen ProcessExpression-style) =====")
  message(paste("  Applying per-junction rank-INT to", nrow(psi_final), "junctions x",
                ncol(psi_final), "samples"))
  # chunked per-junction INT (identical to t(apply(psi_final,1,INT_transform)); avoids the
  # full-matrix row-list + t() temporaries that spike to ~100GB). Pre-allocates one output matrix.
  psi_final <- chunked_row_mat(psi_final, INT_transform)
} else {
  message("\n  --int_transform FALSE: skipping per-phenotype INT (per-sample qqnorm only).")
}

# ===================================================================
# Final PCA on fully cleaned data
# ===================================================================
message("\n===== Final PCA on fully cleaned data =====")

# TOP-FEATURES: compute the final covariate PCs on the same top-N most-variable junctions
# (now INT'd), so standard prcomp stays at ~tens of GB. The written phenotype matrix stays full.
psi_final_t <- t(psi_final[top_junctions, , drop = FALSE])
psi_imputed_final <- psi_final_t
for (i in 1:ncol(psi_imputed_final)) {
  col_mean <- mean(psi_imputed_final[, i], na.rm = TRUE)
  psi_imputed_final[is.na(psi_imputed_final[, i]), i] <- col_mean
}
psi_scaled_final <- scale(psi_imputed_final)
pca_final <- prcomp(psi_scaled_final, center = FALSE, scale. = FALSE)

# Extract all PCs dynamically
num_pcs_final <- ncol(pca_final$x)
message(sprintf("Final PCA: extracted %d PCs", num_pcs_final))

PCs_final <- as.data.frame(pca_final$x, stringsAsFactors = FALSE)
PCs_final$sample <- rownames(psi_final_t)
PCs_final <- PCs_final[, c("sample", setdiff(colnames(PCs_final), "sample"))]

# Plot final PCA (only if we have enough PCs) - following eQTLGen DataQC style
p_final_12 <- ggplot(PCs_final, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.3) +
  theme_bw()

# Only create additional plots if we have those PCs
if (num_pcs_final >= 4) {
  p_final_34 <- ggplot(PCs_final, aes(x = PC3, y = PC4)) +
    geom_point(alpha = 0.3) +
    theme_bw()
} else {
  p_final_34 <- NULL
}

if (num_pcs_final >= 6) {
  p_final_56 <- ggplot(PCs_final, aes(x = PC5, y = PC6)) +
    geom_point(alpha = 0.3) +
    theme_bw()
} else {
  p_final_56 <- NULL
}

if (num_pcs_final >= 8) {
  p_final_78 <- ggplot(PCs_final, aes(x = PC7, y = PC8)) +
    geom_point(alpha = 0.3) +
    theme_bw()
} else {
  p_final_78 <- NULL
}

# Scree plot - show up to 50 PCs or all if fewer (matching the number of PCs saved)
num_scree_pcs <- min(50, num_pcs_final)
variance_explained <- summary(pca_final)$importance[2, 1:num_scree_pcs]
cumulative_variance <- cumsum(variance_explained)

scree_data <- data.frame(
  PC = paste0("PC", 1:num_scree_pcs),
  Variance = variance_explained,
  PC_numeric = 1:num_scree_pcs,
  Cumulative = cumulative_variance
)
scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)

# Individual variance explained (scree plot) - following eQTLGen DataQC style
p_scree <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Proportion of Variance Explained") +
  xlab("Principal Component") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

# Cumulative variance explained (following eQTLGen DataQC style - bar chart, not line)
p_cumulative <- ggplot(scree_data, aes(x = PC, y = Cumulative)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))

# Save final PCA plots - following eQTLGen pattern: combined 2x2 grid as "PCA_after.png"
if (!is.null(p_final_34) && !is.null(p_final_56) && !is.null(p_final_78)) {
  # All 8 PCs available - create 2x2 grid
  p_final_combined <- (p_final_12 | p_final_34) / (p_final_56 | p_final_78)

  ggsave(file.path(args$output_dir, "splice_plots", "PCA_after.png"),
         p_final_combined, height = 7.5, width = 9, units = "in", dpi = 300)
  ggsave(file.path(args$output_dir, "splice_plots", "PCA_after.pdf"),
         p_final_combined, height = 7.5, width = 9, units = "in", dpi = 300)
} else if (!is.null(p_final_34)) {
  # Only first 4 PCs available - create 1x2 layout
  p_final_combined <- (p_final_12 | p_final_34)

  ggsave(file.path(args$output_dir, "splice_plots", "PCA_after.png"),
         p_final_combined, height = 5, width = 9, units = "in", dpi = 300)
  ggsave(file.path(args$output_dir, "splice_plots", "PCA_after.pdf"),
         p_final_combined, height = 5, width = 9, units = "in", dpi = 300)
} else {
  # Only PC1 vs PC2 available
  ggsave(file.path(args$output_dir, "splice_plots", "PCA_after.png"),
         p_final_12, height = 5, width = 5, units = "in", dpi = 300)
  ggsave(file.path(args$output_dir, "splice_plots", "PCA_after.pdf"),
         p_final_12, height = 5, width = 5, units = "in", dpi = 300)
}

ggsave(file.path(args$output_dir, "splice_plots", "PCA_final_scree.png"),
       p_scree, width = 8, height = 6)
ggsave(file.path(args$output_dir, "splice_plots", "PCA_final_scree.pdf"),
       p_scree, width = 8, height = 6)

# Save cumulative variance plot (following eQTLGen pattern)
ggsave(file.path(args$output_dir, "splice_plots", "PCA_cumulative_variance.png"),
       p_cumulative, width = 10, height = 6)
ggsave(file.path(args$output_dir, "splice_plots", "PCA_cumulative_variance.pdf"),
       p_cumulative, width = 10, height = 6)

# ===================================================================
# Save outputs
# ===================================================================
message("\n===== Saving outputs =====")

# Save cleaned PSI matrix
psi_final_df <- as.data.frame(psi_final)
psi_final_df <- cbind(junction = rownames(psi_final_df), psi_final_df)

# Rename sample IDs from RNA IDs to genotype IDs using GTE mapping,
# so splice_data_preprocessed.txt uses the same ID space as CovariatePCs.txt
rna_to_geno <- setNames(as.character(gte$V1), as.character(gte$V2))
sample_cols <- colnames(psi_final_df)[-1]
colnames(psi_final_df)[-1] <- rna_to_geno[sample_cols]
message(paste("  Renamed sample IDs from RNA IDs to genotype IDs using GTE mapping.",
              sum(!is.na(rna_to_geno[sample_cols])), "of", length(sample_cols), "samples mapped."))

preprocessed_file <- file.path(args$output_dir, "splice_data_QCd", "splice_data_preprocessed.txt")
fwrite(psi_final_df, preprocessed_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Generate md5sum checksum (following eQTLGen pattern)
message("  Generating md5sum checksum...")
md5_cmd <- sprintf("md5sum %s > %s",
                   preprocessed_file,
                   file.path(args$output_dir, "splice_data_QCd", "md5sums.txt"))
system(md5_cmd)

# Align PC sample IDs to genotype-ID space, exactly like the QC'd matrix above and like
# eQTLGen ProcessExpression.R (which renames expression columns to genotype IDs BEFORE the
# final PCA, so exp_PCs.txt is written in genotype-ID space). This lets the covariate step
# merge genotype PCs + splicing PCs on SampleID directly, with no GTE remap (MakeCovariateFile.R).
PCs_final$sample <- rna_to_geno[as.character(PCs_final$sample)]

# Save PCA coordinates (for covariate file)
# Following eQTLGen pattern: save up to 50 PCs, or all if fewer
pc_cols <- grep("^PC[0-9]+$", colnames(PCs_final), value = TRUE)
num_pcs_to_save <- min(50, length(pc_cols))
cols_to_save <- c("sample", pc_cols[1:num_pcs_to_save])
message(sprintf("Saving %d PCs to splicing_PCs.txt", num_pcs_to_save))

fwrite(PCs_final[, cols_to_save],
       file.path(args$output_dir, "splice_PCs", "splicing_PCs.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)

# Save variance explained summary (following eQTLGen pattern)
# Create summary_pcs.txt with variance explained for saved PCs
variance_all <- summary(pca_final)$importance[2, 1:num_pcs_to_save]
cumulative_var <- cumsum(variance_all)
pcs_summary <- data.frame(
  PC = paste0("PC", 1:num_pcs_to_save),
  Variance_Explained = variance_all,
  Cumulative_Variance = cumulative_var
)
fwrite(pcs_summary,
       file.path(args$output_dir, "splice_data_summary", "summary_pcs.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)

# Save summary table
fwrite(summary_table,
       file.path(args$output_dir, "splice_data_summary", "summary_table.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)

# Save outlier lists
all_outliers <- c(iterative_result$outliers, mds_outliers)
all_outliers_unique <- unique(all_outliers)
fwrite(data.frame(sample = all_outliers_unique),
       file.path(args$output_dir, "splice_data_summary", "outlier_samples.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Per-junction QC statistics (equivalent to eQTL processed_gene_summary.txt)
message("  Computing per-junction QC statistics...")
# chunked per-row stats (identical to apply(psi_final,1,FUN); each walk is memory-lean)
junction_summary <- data.table(
  junction         = rownames(psi_final),
  mean             = chunked_rowfun_vec(psi_final, mean,   na.rm = TRUE),
  median           = chunked_rowfun_vec(psi_final, median, na.rm = TRUE),
  min              = chunked_rowfun_vec(psi_final, min,    na.rm = TRUE),
  max              = chunked_rowfun_vec(psi_final, max,    na.rm = TRUE),
  sd               = chunked_rowfun_vec(psi_final, sd,     na.rm = TRUE),
  nr_values        = chunked_rowfun_vec(psi_final, function(x) sum(!is.na(x))),
  nr_unique_values = chunked_rowfun_vec(psi_final, function(x) length(unique(x[!is.na(x)])))
)
fwrite(junction_summary,
       file.path(args$output_dir, "splice_data_summary", "processed_junction_summary.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)

# NOTE: HTML report is rendered by the RenderReport process in DataQC.nf
# (bin/Report_DataQc.Rmd), which combines genotype + splicing QC. This QC
# script intentionally produces data + plots only.

message("\n===== Splicing QC Complete! =====")
message(paste("Output directory:", args$output_dir))
message(paste("Final PSI matrix:", nrow(psi_final), "junctions ×", ncol(psi_final), "samples"))
message(paste("Total outliers removed:", length(all_outliers_unique)))
message(paste("PCA outliers:", length(iterative_result$outliers)))
message(paste("MDS outliers:", length(mds_outliers)))
