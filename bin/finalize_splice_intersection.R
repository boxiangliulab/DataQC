#!/usr/bin/env Rscript

# finalize_splice_intersection.R
# ACTIVE expression QC (DataQC mod 1): re-align the splicing phenotype + covariate PCs to the final
# analysed sample set E = genotype n splice n expression-QC.
#
# ProcessSplicing.R wrote the phenotype (per-feature rank-INT) + covariate PCs on the splice-QC'd set S.
# ProcessExpression.R (--remove_outliers true) then removed expression outliers, giving E (subset of S) in
# exp_data_summary/final_samples.txt. Here we:
#   (1) subset the phenotype SAMPLE columns to E,
#   (2) RE-APPLY the per-feature rank-INT across E. This is BIT-EXACT INT-on-E: rank-INT is a monotone
#       transform, so re-ranking the already-INT'd values equals the INT of the originals restricted to E
#       (no need to re-read the raw Stage-0 matrix),
#   (3) recompute the splicing covariate PCs on E from the saved top_junctions basis, using the SAME
#       prcomp recipe as ProcessSplicing (mean-impute NA per column -> scale() -> prcomp(center=FALSE)),
#   (4) append a summary stage.
#
# This operates on the TALL phenotype (features x samples): data.table::fread handles tall matrices fine,
# whereas it segfaults on the WIDE (samples x 270k-feature) layout. The caller (FINALIZE_SPLICING) therefore
# transposes wide->tall before this script and tall->wide after. If E == S (no expression outliers removed)
# this is an idempotent no-op.

suppressPackageStartupMessages({ library(optparse); library(data.table) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--tall", type = "character", help = "TALL phenotype (feature id col + sample cols), rewritten in place"),
  make_option("--splice_dir", type = "character", help = "outputfolder_splice (for top_junctions / PCs / summary)"),
  make_option("--final_samples", type = "character", help = "final_samples.txt (E, one genotype ID per line)"),
  make_option("--num_splicing_pcs", type = "integer", default = 50L, help = "max PCs to save [default %default]")
)))

INT_transform <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
ROW_CHUNK <- 20000L

topj_file <- file.path(opt$splice_dir, "splice_PCs", "top_junctions.txt")
pcs_file  <- file.path(opt$splice_dir, "splice_PCs", "splicing_PCs.txt")
summ_file <- file.path(opt$splice_dir, "splice_data_summary", "summary_table.txt")
if (!file.exists(opt$tall)) stop("tall phenotype not found: ", opt$tall)
if (!file.exists(opt$final_samples)) stop("final_samples not found: ", opt$final_samples)

E <- readLines(opt$final_samples); E <- E[E != ""]
message(sprintf("FINALIZE_SPLICING: final set E = %d samples.", length(E)))

# ---- read TALL phenotype (feature x sample) -----------------------------------------------------------
tall <- fread(opt$tall, data.table = FALSE)          # col 1 = feature id; remaining = samples
feat <- as.character(tall[[1]])
samp <- colnames(tall)[-1]
keep <- samp %in% E
message(sprintf("FINALIZE_SPLICING: phenotype has %d samples; keeping %d in E (dropping %d).",
                length(samp), sum(keep), length(samp) - sum(keep)))
if (sum(keep) < 10) stop("Fewer than 10 phenotype samples fall in E; aborting (check final_samples / ID space).")

M <- as.matrix(tall[, c(FALSE, keep), drop = FALSE]); storage.mode(M) <- "double"
rownames(M) <- feat
kept_samp <- samp[keep]

# ---- (2) re-apply per-feature rank-INT across E (exact INT-on-E; chunked per row) ---------------------
for (s in seq.int(1L, nrow(M), by = ROW_CHUNK)) {
  idx <- s:min(s + ROW_CHUNK - 1L, nrow(M))
  M[idx, ] <- t(apply(M[idx, , drop = FALSE], 1, INT_transform))
}

# write finalized TALL back (feature id col named like the input's first column)
out_tall <- data.frame(M, check.names = FALSE)
out_tall <- cbind(feat, out_tall); colnames(out_tall)[1] <- colnames(tall)[1]
fwrite(out_tall, opt$tall, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("FINALIZE_SPLICING: re-INT'd %d features on E and rewrote the tall phenotype.", nrow(M)))

# ---- (3) recompute splicing covariate PCs on E (ProcessSplicing recipe) -------------------------------
if (file.exists(topj_file)) {
  topj <- readLines(topj_file); topj <- topj[topj != ""]
  sel  <- feat %in% topj
  Mt <- t(M[sel, , drop = FALSE])                    # samples_E x top (INT'd on E), rownames = samples
  for (i in seq_len(ncol(Mt))) { cm <- mean(Mt[, i], na.rm = TRUE); Mt[is.na(Mt[, i]), i] <- cm }
  Ms <- scale(Mt)
  ok <- apply(Ms, 2, function(z) all(is.finite(z)))  # drop zero-variance cols so prcomp is well-defined
  Ms <- Ms[, ok, drop = FALSE]
  pca <- prcomp(Ms, center = FALSE, scale. = FALSE)
  PCs <- as.data.frame(pca$x, stringsAsFactors = FALSE)
  pc_cols <- grep("^PC[0-9]+$", colnames(PCs), value = TRUE)
  np <- min(opt$num_splicing_pcs, length(pc_cols))
  out <- data.frame(sample = rownames(Mt), PCs[, pc_cols[1:np], drop = FALSE], check.names = FALSE)
  fwrite(out, pcs_file, sep = "\t", quote = FALSE, row.names = FALSE)
  ve <- summary(pca)$importance[2, 1:np]
  fwrite(data.frame(PC = paste0("PC", 1:np), Variance_Explained = ve, Cumulative_Variance = cumsum(ve)),
         file.path(opt$splice_dir, "splice_data_summary", "summary_pcs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("FINALIZE_SPLICING: recomputed %d splicing PCs on E from %d top features.", np, ncol(Ms)))
} else {
  message("FINALIZE_SPLICING: WARNING top_junctions.txt not found -> covariate PCs NOT recomputed.")
}

# ---- (4) append a summary stage -----------------------------------------------------------------------
if (file.exists(summ_file)) {
  st <- fread(summ_file, data.table = FALSE)
  if (all(c("Stage", "Nr_of_junctions", "Nr_of_samples") %in% colnames(st))) {
    add <- data.frame(Stage = "After expression-QC intersection (geno n splice n expr)",
                      Nr_of_junctions = nrow(M), Nr_of_samples = length(kept_samp), stringsAsFactors = FALSE)
    st <- rbind(st, add[, colnames(st)])
    fwrite(st, summ_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
message("FINALIZE_SPLICING: done.")
