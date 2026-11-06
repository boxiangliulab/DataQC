#!/usr/bin/env Rscript

# ProcessExpression.R  —  FAITHFUL port of the eQTLGen ProcessExpression.R (RNA-seq path) for the
# splicing workflow. The normalization (remove 0-variance -> CPM>0.5 in >=1% of samples -> TMM -> cpm
# -> per-gene INT), the iterative PCA outlier detection, the MDS outlier detection, and the
# expression-based sex / contamination check (XIST vs Y genes, ContaminationArea/Slope geometry), plus
# the output files, are kept identical to eQTLGen. The ONLY intentional differences are:
#   (1) all three QC stages are DIAGNOSTIC by default — outliers / sex-mismatches / contaminated samples
#       are detected, plotted and written out, but NOT removed. Pass --remove_outliers true to restore
#       eQTLGen's original removal behaviour.
#   (2) the analysed sample set is aligned to the SPLICING-QC'd samples (--sample_list, read via cut -f1),
#       so expression and splicing share an identical sample set for the Stage-4 HASE encoding.
# Only the RNA-seq platform is ported (the cohorts here are RNA-seq gene counts).

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(optparse)
  library(patchwork)
  library(MASS)
  library(dplyr)
})
setDTthreads(4)

option_list <- list(
  make_option(c("-e", "--expression_matrix"), type = "character",
    help = "Unprocessed RNA-seq gene count matrix. Genes in rows (first col = gene ID), samples in columns (RNA IDs)."),
  make_option(c("-l", "--gte"), type = "character",
    help = "GTE file, no header: col1 = genotype ID, col2 = RNA ID."),
  make_option(c("-m", "--emp_probe_mapping"), type = "character",
    help = "Ensembl->chromosome map (cols Ensembl, chromosome_name); used to locate XIST + Y genes for the sex check."),
  make_option(c("-s", "--sd"), type = "double", default = 4,
    help = "SD threshold for PC1/PC2 and MDS1/MDS2 outlier flagging (eQTLGen default 4) [default %default]"),
  make_option(c("-c", "--contamination_area"), type = "double", default = 0.3,
    help = "Angular width of the contamination gray-zone around the sex line [default %default]"),
  make_option(c("-a", "--contamination_slope"), type = "double", default = 45,
    help = "Angle of the sex-discrimination line (degrees) [default %default]"),
  make_option(c("-i", "--sex_info"), type = "character",
    help = "Genotype SexCheck.txt (cols incl. IID, SNPSEX 1=male/2=female). Used for the sex check + MDS sex annotation."),
  make_option(c("-f", "--geno_filter"), type = "character",
    help = "Plink .fam of QC'd genotype samples (col2 = genotype ID)."),
  make_option(c("--sample_list"), type = "character",
    help = "Splicing-QC'd matrix (samples x features); only its first column (cut -f1) = QC'd sample IDs (genotype-ID space) to align to."),
  make_option(c("--remove_outliers"), type = "logical", default = FALSE,
    help = "FALSE (default) = diagnostic: detect & report, keep all samples. TRUE = remove outliers like eQTLGen."),
  make_option(c("-o", "--output"), type = "character")
)
args <- parse_args(OptionParser(option_list = option_list))
REMOVE <- isTRUE(args$remove_outliers)

dir.create(args$output, showWarnings = FALSE)
dir.create(paste0(args$output, "/exp_plots"), showWarnings = FALSE)
dir.create(paste0(args$output, "/exp_data_QCd"), showWarnings = FALSE)
dir.create(paste0(args$output, "/exp_PCs"), showWarnings = FALSE)
dir.create(paste0(args$output, "/exp_data_summary"), showWarnings = FALSE)

message(sprintf("Expression QC (%s). SD threshold: %s",
                if (REMOVE) "REMOVING outliers (eQTLGen behaviour)" else "DIAGNOSTIC: flag & keep", args$sd))

#############
# functions #
#############
INT_transform <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
# Map genetic sex (SNPSEX 1=male, 2=female) to a plotting factor; anything else (0 / unknown / NA) ->
# "Unknown" so the diagnostic PCA/MDS points still render for cohorts with no determined sex (e.g. the
# genotype has no usable chrX, so the sex check returns 0 for all samples).
to_sex_factor <- function(x) factor(ifelse(x %in% c(1, "1"), "M", ifelse(x %in% c(2, "2"), "F", "Unknown")), levels = c("F", "M", "Unknown"))
shap_test <- function(x) if (length(unique(x)) > 1) shapiro.test(x)$p.value else NA

# splicing-QC'd sample universe (genotype-ID space), read cheaply via cut -f1 (never loads the wide matrix)
samples_keep <- as.character(fread(cmd = paste("cut -f1", shQuote(args$sample_list)), header = TRUE)[[1]])
message(sprintf("Aligning expression to %d splicing-QC'd samples.", length(samples_keep)))

# gte (global): used by IterativeOutlierDetection to map genotype IDs back to RNA IDs, like eQTLGen
gte_main <- fread(args$gte, header = FALSE, colClasses = "character")

# eQTLGen RNAseq_preprocess, with the sample filter aligned to the splicing-QC'd set (samples_keep)
RNAseq_preprocess <- function(exp, normalize = TRUE, gene_inclusion = NULL) {
  colnames(exp)[1] <- "Probe"
  exp$Probe <- as.character(exp$Probe)
  exp <- as.data.frame(exp)
  rownames(exp) <- exp[, 1]
  exp <- exp[, -1]
  exp <- abs(as.matrix(exp))

  gte <- fread(args$gte, header = FALSE, colClasses = "character")
  geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1, 2)))
  gte <- gte[gte$V1 %in% geno_fam$IID, ]
  gte <- gte[gte$V1 %in% samples_keep, ]                       # <-- align to splicing-QC'd samples
  gte <- gte[gte$V2 %in% as.character(colnames(exp)), ]

  exp <- exp[, colnames(exp) %in% gte$V2]
  exp <- exp[, base::order(colnames(exp))]
  gte <- gte[base::order(gte$V2), ]
  if (!all(colnames(exp) == gte$V2)) stop("Something went wrong matching genotype and expression IDs.")
  colnames(exp) <- gte$V1

  gene_variance <- apply(exp, 1, var)
  exp <- exp[gene_variance > 0, ]

  exp_keep <- DGEList(counts = exp)
  n_pass <- rowSums(cpm(exp_keep, log = FALSE) > 0.5)
  keep <- n_pass >= round(ncol(exp) / 100, 0)
  if (!is.null(gene_inclusion) && length(gene_inclusion) > 0) {
    incl <- intersect(gene_inclusion, names(keep))
    keep[incl] <- TRUE
  }
  exp <- exp[rownames(exp) %in% names(keep[keep == TRUE]), ]

  if (normalize) {
    exp_n <- calcNormFactors(DGEList(counts = exp), method = "TMM")
    exp_n <- cpm(exp_n, log = FALSE)
  } else {
    exp_n <- exp
  }
  return(exp_n)
}

exp_summary <- function(x) {
  data.table(
    gene = rownames(x),
    mean = apply(x, 1, mean), median = apply(x, 1, median),
    min = apply(x, 1, min), max = apply(x, 1, max), sd = apply(x, 1, sd),
    nr_values = apply(x, 1, length),
    nr_unique_values = apply(x, 1, function(z) length(unique(z))),
    nr_of_zero_values = apply(x, 1, function(z) length(z[z == 0])),
    shapiro_P = apply(x, 1, shap_test))
}

rnaseq_cpm_summary <- function(exp) {
  colnames(exp)[1] <- "Probe"
  exp$Probe <- as.character(exp$Probe)
  exp <- as.data.frame(exp); rownames(exp) <- exp[, 1]; exp <- exp[, -1]; exp <- abs(as.matrix(exp))
  gte <- fread(args$gte, header = FALSE, colClasses = "character")
  geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1, 2)))
  gte <- gte[gte$V1 %in% geno_fam$IID, ]; gte <- gte[gte$V1 %in% samples_keep, ]
  gte <- gte[gte$V2 %in% as.character(colnames(exp)), ]
  exp <- exp[, colnames(exp) %in% gte$V2]; exp <- exp[, base::order(colnames(exp))]; gte <- gte[base::order(gte$V2), ]
  colnames(exp) <- gte$V1
  gene_variance <- apply(exp, 1, var); exp <- exp[gene_variance > 0, ]
  exp_keep <- DGEList(counts = exp)
  n1 <- rowSums(cpm(exp_keep, log = FALSE) > 0.5)
  exp_keep <- calcNormFactors(exp_keep)
  n2 <- rowSums(cpm(exp_keep, log = FALSE) > 0.5)
  data.frame(n_samples_with_cpm_threshold_passed = n1,
             n_samples_with_cpm_threshold_passed_normalized = n2,
             n_samples_threshold = round(ncol(exp) / 100, 0))
}

# eQTLGen iterative PCA outlier detection (RNA-seq). Internally converges to find the full outlier set;
# the CALLER decides whether to remove them (diagnostic vs REMOVE).
IterativeOutlierDetection <- function(input_exp, sd_threshold = 4, gene_inclusion = NULL) {
  and <- input_exp
  list_ggplots <- list(); summary_pcs <- list(); outliers <- c()
  it_round <- 0; plot_it_round <- 0; nr_outliers <- 1; PCs_it1 <- NULL
  message("Starting iterative PCA outlier detection...")
  while (nr_outliers > 0) {
    and_p <- RNAseq_preprocess(and, gene_inclusion = gene_inclusion)
    and_p <- log2(and_p + 0.25)
    it_round <- it_round + 1
    pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
    PCs <- as.data.table(pcs$x); PCs$sample <- rownames(pcs$x)
    np <- min(50, ncol(pcs$x))
    importance <- pcs$sdev^2 / sum(pcs$sdev^2)
    summary_pcs[[it_round]] <- data.table(PC = paste0("PC", 1:np), explained_variance = importance[1:np])
    PCs$outlier <- "no"
    m1 <- mean(PCs$PC1); s1 <- sd(PCs$PC1)
    PCs$outlier[PCs$PC1 > m1 + sd_threshold * s1 | PCs$PC1 < m1 - sd_threshold * s1] <- "yes"
    m2 <- mean(PCs$PC2); s2 <- sd(PCs$PC2)
    PCs$outlier[PCs$PC2 > m2 + sd_threshold * s2 | PCs$PC2 < m2 - sd_threshold * s2] <- "yes"
    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier)) +
      geom_point(alpha = 0.3) + theme_bw() + scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    non_outliers <- PCs$sample[PCs$outlier != "yes"]
    non_outliers_rna <- gte_main[gte_main$V1 %in% non_outliers, ]$V2
    and <- and[, colnames(and) %in% c(non_outliers_rna, "Feature"), with = FALSE]
    nr_outliers <- length(PCs$sample[PCs$outlier == "yes"])
    if (nr_outliers > 0) outliers <- c(outliers, PCs$sample[PCs$outlier == "yes"])
    if (it_round == 1) { PCs_it1 <- PCs } else { PCs_it1$outlier[PCs_it1$sample %in% outliers] <- "yes" }
    if (it_round >= 10) { message("Reached 10 iterations; stopping."); break }
  }
  list(exp_mat = and_p, plots = list_ggplots, summary_pcs = summary_pcs,
       PCs_first = PCs_it1, outliers = unique(outliers))
}

# eQTLGen expression-based sex / contamination check (XIST vs Y genes, area/slope geometry).
# Returns the flags; removal is applied by the caller only if REMOVE.
ExpressionBasedSampleSwap <- function(and) {
  emp <- fread(args$emp_probe_mapping, keepLeadingZeros = TRUE)
  sex_specific_genes <- c("ENSG00000229807", emp$Ensembl[emp$chromosome_name == "Y"])
  and_pp <- RNAseq_preprocess(and, gene_inclusion = sex_specific_genes)
  xist <- and_pp[rownames(and_pp) == "ENSG00000229807", ]
  y_ids <- emp$Ensembl[emp$chromosome_name == "Y"]
  y_genes <- and_pp[rownames(and_pp) %in% y_ids, , drop = FALSE]
  nr_of_y_genes <- nrow(y_genes)
  message(paste("Nr. of Y genes:", nr_of_y_genes))
  if (nr_of_y_genes == 0) { message("No Y genes; omitting sex check."); return(NULL) }

  y_mean <- apply(y_genes, 2, mean); min_y_mean <- min(y_mean)
  if (length(xist) == 0) {
    message("XIST not detected; using mean of X-chr genes instead.")
    xg <- and_pp[rownames(and_pp) %in% emp$Ensembl[emp$chromosome_name == "X"], , drop = FALSE]
    xist_mean <- apply(xg, 2, mean); min_xist <- min(xist_mean)
    yd <- data.frame(sample = colnames(y_genes), xist = xist_mean - min_xist, y_genes = y_mean - min_y_mean)
  } else {
    min_xist <- min(xist)
    yd <- data.frame(sample = colnames(y_genes), xist = xist - min_xist, y_genes = y_mean - min_y_mean)
  }

  geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1, 2)))
  gf <- geno_fam[, c(2, 4), with = FALSE]; colnames(gf) <- c("sample", "Sex"); gf$Sex <- as.character(gf$Sex)
  yd <- merge(yd, gf, by = "sample")
  max_exp <- max(yd$y_genes, yd$xist)

  lower_slope  <- tan((args$contamination_slope - args$contamination_area / 2) / 180 * pi)
  upper_slope  <- tan((args$contamination_slope + args$contamination_area / 2) / 180 * pi)
  middle_slope <- tan(args$contamination_slope / 180 * pi)

  yd$expressionSexNaive <- case_when(yd$y_genes > yd$xist * middle_slope ~ 1,
                                     yd$y_genes < yd$xist * middle_slope ~ 2)
  x_expr_med <- if (nrow(yd[yd$Sex == 1 & yd$expressionSexNaive == 1, ]) > 0) median(yd[yd$Sex == 1 & yd$expressionSexNaive == 1, "xist"]) else 0
  y_expr_med <- if (nrow(yd[yd$Sex == 2 & yd$expressionSexNaive == 2, ]) > 0) median(yd[yd$Sex == 2 & yd$expressionSexNaive == 2, "y_genes"]) else 0
  x_expr_min <- if (nrow(yd[yd$Sex == 2 & yd$expressionSexNaive == 2, ]) > 0) min(yd[yd$Sex == 2 & yd$expressionSexNaive == 2, "xist"]) else 0
  y_expr_min <- if (nrow(yd[yd$Sex == 1 & yd$expressionSexNaive == 1, ]) > 0) min(yd[yd$Sex == 1 & yd$expressionSexNaive == 1, "y_genes"]) else 0
  y_max <- max(yd$y_genes); x_max <- max(yd$xist)
  x_anchor <- min(x_expr_med, x_expr_min - 0.01 * x_max)
  y_anchor <- min(y_expr_med, y_expr_min - 0.01 * y_max)
  yd$xist_corrected <- yd$xist - x_anchor

  yd$contaminated <- case_when(
    (yd$y_genes > (yd$xist_corrected * lower_slope + y_anchor) &
     yd$y_genes < (yd$xist_corrected * upper_slope + y_anchor)) ~ "yes", TRUE ~ "no")
  yd$expressionSex <- case_when(
    yd$y_genes < (yd$xist_corrected * middle_slope + y_anchor) ~ 2,
    yd$y_genes > (yd$xist_corrected * middle_slope + y_anchor) ~ 1)
  yd$mismatch <- case_when(yd$Sex == 0 ~ "unknown",
                           yd$expressionSex == yd$Sex ~ "no",
                           yd$expressionSex != yd$Sex ~ "yes")
  yd$status <- case_when(
    yd$contaminated == "yes" & yd$mismatch == "yes" ~ "Contaminated and\nsex mismatch",
    yd$contaminated == "yes" ~ "Likely contaminated",
    yd$mismatch == "yes" ~ "Sex mismatch", TRUE ~ "Passed")

  exclusion_zone <- tibble(x = c(x_expr_med, max_exp)) %>%
    mutate(lower_bound = (x - x_expr_med) * lower_slope + y_expr_med,
           upper_bound = (x - x_expr_med) * upper_slope + y_expr_med,
           middle_line = (x - x_expr_med) * middle_slope + y_expr_med)
  p_sex <- ggplot(exclusion_zone, aes(x = x, ymin = lower_bound, ymax = upper_bound)) +
    geom_ribbon(alpha = 0.2) + geom_line(aes(x = x, y = middle_line), linetype = 2, colour = "blue") +
    geom_point(data = yd, inherit.aes = FALSE, aes(col = status, shape = Sex, x = xist, y = y_genes)) +
    scale_colour_manual(values = alpha(c("Passed" = "black", "Likely contaminated" = "orange",
      "Sex mismatch" = "red", "Contaminated and\nsex mismatch" = "darkmagenta"), 0.5), name = "Status") +
    coord_cartesian(ylim = c(0, max_exp), xlim = c(0, max_exp)) + theme_bw() +
    ylab(paste0("mean Y genes - min (n=", nr_of_y_genes, ")")) + xlab("XIST - min(XIST)") +
    ggtitle("Expression-based sex / contamination check")
  ggsave(paste0(args$output, "/exp_plots/SexSpecificGenesXIST.png"), p_sex, height = 5, width = 7, dpi = 300)
  ggsave(paste0(args$output, "/exp_data_summary/expression_sex_check.png"), p_sex, height = 5, width = 7, dpi = 300)

  comb <- merge(gte_main, yd, by.x = "V1", by.y = "sample")
  colnames(comb)[c(1, 2)] <- c("Genotype_SampleID", "Expression_SampleID")
  fwrite(comb, paste0(args$output, "/exp_data_summary/SexMismatchStatus.txt"), sep = "\t", quote = FALSE)
  writeLines(yd$sample[yd$mismatch == "yes"], paste0(args$output, "/exp_data_summary/sex_discordant_samples.txt"))
  message(sprintf("Sex check: %d mismatch, %d likely-contaminated (diagnostic%s).",
                  sum(yd$mismatch == "yes", na.rm = TRUE), sum(yd$contaminated == "yes"),
                  if (REMOVE) "; REMOVED" else ", not removed"))
  return(yd)
}

############
# Analysis #
############
and <- fread(args$expression_matrix, header = TRUE, keepLeadingZeros = TRUE)
colnames(and) <- as.character(colnames(and)); colnames(and)[1] <- "Feature"; and$Feature <- as.character(and$Feature)
message(paste("Initially:", nrow(and), "genes and", ncol(and) - 1, "samples"))
summary_table <- data.table(Stage = "Unprocessed matrix", Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1)

# Restrict samples to the QC'd genotype set AND the splicing-QC'd set
geno_fam <- fread(args$sex_info, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1, 2)))
geno_fam$IID <- as.character(geno_fam$IID)
# Sex lookup (genotype-ID space) used as the point SHAPE in the diagnostic expression PCA/MDS plots, like eQTLGen.
sex_lu <- geno_fam[, c(2, 4), with = FALSE]; colnames(sex_lu) <- c("Sample", "Sex")
sex_lu$Sample <- as.character(sex_lu$Sample); sex_lu$Sex <- to_sex_factor(sex_lu$Sex)
gen_filter <- fread(args$geno_filter, header = FALSE, keepLeadingZeros = TRUE, colClasses = list(character = c(1, 2)))
gte <- copy(gte_main)
gte <- gte[gte$V1 %in% geno_fam$IID & gte$V1 %in% as.character(gen_filter$V2) & gte$V1 %in% samples_keep, ]
and <- and[, colnames(and) %in% c("Feature", gte$V2), with = FALSE]
summary_table <- rbind(summary_table, data.table(Stage = "Samples overlapping splicing-QC'd + genotype-QC'd set",
                       Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1))
if (ncol(and) - 1 < 10) stop("Fewer than 10 samples overlap splicing-QC'd + genotype data.")

# ---- (1) PCA outlier detection (diagnostic unless REMOVE) ----
io <- IterativeOutlierDetection(and, sd_threshold = args$sd)
pca_outliers <- io$outliers
# diagnostic plot (expression report): PCA BEFORE outlier removal, outliers from all iterations in red,
# genetic/reported sex as the point shape (eQTLGen Fig 12 analog).
if (!is.null(io$PCs_first)) {
  pcb <- merge(io$PCs_first, sex_lu, by.x = "sample", by.y = "Sample", all.x = TRUE)
  ggplot(pcb, aes(x = PC1, y = PC2, colour = outlier, shape = Sex)) + geom_point(alpha = 0.4) +
    theme_bw() + scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
    ggtitle("Expression PCA (before outlier removal)")
  ggsave(paste0(args$output, "/exp_plots/PCA_before.png"), height = 7.5, width = 9, dpi = 300)
}
if (REMOVE && length(pca_outliers) > 0) {
  keep_rna <- gte[gte$V1 %in% setdiff(samples_keep, pca_outliers), ]$V2
  and <- and[, colnames(and) %in% c("Feature", keep_rna), with = FALSE]
}
summary_table <- rbind(summary_table, data.table(
  Stage = if (REMOVE) "After removal of PCA outliers" else "PCA outliers (flagged, not removed)",
  Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1))

# ---- (2) MDS outlier detection (diagnostic unless REMOVE) ----
and_p <- RNAseq_preprocess(and); and_p <- log2(and_p + 0.25)
dist <- cor(and_p, method = "pearson")
mds <- as.data.frame(isoMDS(1 - dist, k = 2)$points); colnames(mds) <- c("MDS1", "MDS2"); mds$Sample <- colnames(and_p)
mds$outlier <- "no"
m1 <- mean(mds$MDS1); s1 <- sd(mds$MDS1); m2 <- mean(mds$MDS2); s2 <- sd(mds$MDS2)
mds$outlier[mds$MDS1 > m1 + args$sd * s1 | mds$MDS1 < m1 - args$sd * s1 |
            mds$MDS2 > m2 + args$sd * s2 | mds$MDS2 < m2 - args$sd * s2] <- "yes"
sex2 <- geno_fam[, c(2, 4), with = FALSE]; colnames(sex2) <- c("Sample", "Sex"); sex2$Sample <- as.character(sex2$Sample); sex2$Sex <- to_sex_factor(sex2$Sex)
mds <- merge(mds, sex2, by = "Sample")
ggplot(mds, aes(x = MDS1, y = MDS2, colour = outlier, shape = Sex)) + geom_point(alpha = 0.3) + theme_bw() +
  scale_colour_manual(values = c("no" = "black", "yes" = "red")) + ggtitle("MDS outlier detection")
ggsave(paste0(args$output, "/exp_plots/MDS_before.png"), height = 5, width = 6, dpi = 300)
mds_outliers <- mds$Sample[mds$outlier == "yes"]
if (REMOVE && length(mds_outliers) > 0) {
  keep_rna <- gte[gte$V1 %in% setdiff(colnames(and_p), mds_outliers), ]$V2
  and <- and[, colnames(and) %in% c("Feature", keep_rna), with = FALSE]
}
summary_table <- rbind(summary_table, data.table(
  Stage = if (REMOVE) "After removal of MDS outliers" else "MDS outliers (flagged, not removed)",
  Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1))

# ---- (3) expression-based sex / contamination check (diagnostic unless REMOVE) ----
yd <- ExpressionBasedSampleSwap(and)
sex_mismatch <- if (!is.null(yd)) yd$sample[yd$mismatch == "yes"] else character(0)
if (REMOVE && !is.null(yd)) {
  keep_gids <- yd$sample[yd$mismatch != "yes" & yd$contaminated != "yes"]
  keep_rna <- gte[gte$V1 %in% keep_gids, ]$V2
  and <- and[, colnames(and) %in% c("Feature", keep_rna), with = FALSE]
  summary_table <- rbind(summary_table, data.table(Stage = "After removal of sex/contamination errors",
                         Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1))
}

# combined diagnostic outlier list (PCA + MDS), genotype-ID space
writeLines(unique(c(pca_outliers, mds_outliers)), paste0(args$output, "/exp_data_summary/potential_outlier_samples.txt"))

# ---- Final: normalize + per-gene INT + write (eQTLGen output set) ----
gene_cpm_summary <- rnaseq_cpm_summary(and)
fwrite(gene_cpm_summary, paste0(args$output, "/exp_data_summary/cpm_summary.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

and_pp <- RNAseq_preprocess(and)
summary_table <- rbind(summary_table, data.table(Stage = "After filtering genes (CPM>0.5 in >=1% of samples)",
                       Nr_of_features = nrow(and_pp), Nr_of_samples = ncol(and_pp)))

and_p <- t(apply(and_pp, 1, INT_transform))       # per-gene INT -> genes x samples
# order samples to match the splicing-QC'd order (for the shared Stage-4 HASE encoding)
ord <- samples_keep[samples_keep %in% colnames(and_p)]
# HARD GUARD: the expression matrix feeds the Stage-4 permuted / LD-estimation path, whose LD reference
# MUST be estimated on exactly the same individuals as the sQTL (splicing) meta-analysis. The expression
# set is therefore required to cover the FULL splicing-QC'd set. If any splicing-QC'd sample lacks
# expression data, stop (resolve upstream) rather than silently shipping a smaller / mismatched matrix.
missing_expr <- setdiff(samples_keep, colnames(and_p))
if (length(missing_expr) > 0) {
  writeLines(missing_expr, paste0(args$output, "/exp_data_summary/splicing_samples_missing_expression.txt"))
  if (!REMOVE) {
    # DIAGNOSTIC mode: the splicing set is NOT modified, so the expression (LD) set must fully cover it.
    stop(sprintf(paste0("%d splicing-QC'd sample(s) have no expression data, so the expression matrix ",
      "(Stage-4 permuted / LD path) would NOT match the sQTL sample set. ",
      "See exp_data_summary/splicing_samples_missing_expression.txt. ",
      "Resolve upstream (e.g. drop these samples from the splicing set) and re-run."), length(missing_expr)))
  } else {
    # ACTIVE mode: these samples (expression outliers and/or lacking expression) are intentionally dropped;
    # FINALIZE_SPLICING re-aligns the splicing phenotype + PCs to the intersection, so LD set == sQTL set.
    message(sprintf(paste0("Expression QC ACTIVE: %d splicing-QC'd sample(s) are not in the final expression set ",
      "(removed as expression outliers and/or lacking expression data). The splicing phenotype + PCs will be ",
      "re-aligned to this %d-sample intersection (geno n splice n expr) by FINALIZE_SPLICING."),
      length(missing_expr), length(ord)))
  }
}
and_p <- and_p[, ord, drop = FALSE]; and_pp <- and_pp[, ord, drop = FALSE]
# Final analysed sample set E = geno n splice n expr (genotype-ID space, splicing order). FINALIZE_SPLICING
# re-aligns the splicing phenotype + covariate PCs to exactly this set so both Stage-4 paths share one set.
writeLines(colnames(and_p), paste0(args$output, "/exp_data_summary/final_samples.txt"))

pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
PCs <- data.table(Sample = rownames(pcs$x), pcs$x)
# diagnostic plot (expression report): PCA AFTER, on the fully processed (INT-normalised) matrix; samples
# flagged by PCA/MDS are outlined red, sex as shape (eQTLGen Fig 14 analog). Diagnostic: nothing removed.
pca_after <- merge(PCs[, .(Sample, PC1, PC2)], sex_lu, by = "Sample", all.x = TRUE)
pca_after$outlier <- ifelse(pca_after$Sample %in% unique(c(pca_outliers, mds_outliers)), "yes", "no")
ggplot(pca_after, aes(x = PC1, y = PC2, colour = outlier, shape = Sex)) + geom_point(alpha = 0.4) +
  theme_bw() + scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
  ggtitle("Expression PCA (final, INT-normalised)")
ggsave(paste0(args$output, "/exp_plots/PCA_after.png"), height = 7.5, width = 9, dpi = 300)

# INT'd matrix (samples x genes) — the Stage-4 permuted/LD input
and_p2 <- data.table(ID = colnames(and_p), as.data.table(t(and_p)))
fwrite(and_p2, paste0(args$output, "/exp_data_QCd/exp_data_preprocessed.txt"), sep = "\t", quote = FALSE)
# non-INT log2 matrix (eQTLGen also writes this)
and_log2 <- log2(and_pp + 1)
fwrite(data.table(ID = colnames(and_log2), as.data.table(t(and_log2))),
       paste0(args$output, "/exp_data_QCd/exp_data_preprocessed_noINT_log2.txt"), sep = "\t", quote = FALSE)

importance <- pcs$sdev^2 / sum(pcs$sdev^2)
np <- min(100, nrow(PCs), length(importance))
summary_pcs <- data.table(PC = paste0("PC", 1:np), explained_variance = importance[1:np])
fwrite(PCs[, 1:min(101, ncol(PCs)), with = FALSE], paste0(args$output, "/exp_PCs/exp_PCs.txt"), sep = "\t", quote = FALSE)
fwrite(summary_pcs, paste0(args$output, "/exp_data_summary/summary_pcs.txt"), sep = "\t", quote = FALSE)
ggplot(summary_pcs, aes(x = factor(PC, levels = PC), y = explained_variance)) + geom_bar(stat = "identity") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab("PC")
ggsave(paste0(args$output, "/exp_plots/PCA_final_scree_plot.png"), height = 6, width = 17, dpi = 300)

# per-gene summaries (raw + processed), eQTLGen
fwrite(exp_summary(RNAseq_preprocess(and, normalize = FALSE)),
       paste0(args$output, "/exp_data_summary/raw_gene_summary.txt"), sep = "\t", quote = FALSE)
fwrite(exp_summary(log2(RNAseq_preprocess(and, normalize = TRUE) + 0.25)),
       paste0(args$output, "/exp_data_summary/processed_gene_summary.txt"), sep = "\t", quote = FALSE)

fwrite(summary_table, paste0(args$output, "/exp_data_summary/summary_table.txt"), sep = "\t", quote = FALSE)
message("Expression QC done -> ", args$output)
