#!/usr/bin/env Rscript
# Filter blanks, QC samples and normalize feature table.

library(optparse)
library(matrixStats)

option_list <- list(
    make_option("--incsv", type="character", help="Annotated feature CSV"),
    make_option("--samples", type="character", help="Sample sheet CSV"),
    make_option("--outfiltered", type="character", help="Filtered feature CSV"),
    make_option("--outnorm", type="character", help="Normalized feature CSV"),
    make_option("--blank_cutoff", type="double", default=0.25, help="Blank cutoff ratio"),
    make_option("--qc_cv", type="double", default=0.3, help="Max CV for QC samples"),
    make_option("--norm", type="character", default="TIC", help="Normalization method [TIC|PQN|none]")
)
opt <- parse_args(OptionParser(option_list=option_list))

feat <- read.csv(opt$incsv)
samples <- read.csv(opt$samples)

# Find sample columns (exclude metadata cols)
meta_cols <- c("mz", "rt_min", "rt_max", "mzmed", "rtmed", "rtmin", "rtmax",
               "isotopes", "adduct", "pcgroup", "neutral_mass", "charge")
sample_cols <- setdiff(colnames(feat), meta_cols)
sample_cols <- sample_cols[sample_cols %in% samples$sample_id]

cat("[FILTER] Samples found:", length(sample_cols), "\n")

# Blank filter: remove features where blank mean >= blank_cutoff * sample mean
blank_cols <- sample_cols[samples$sample_type[samples$sample_id %in% sample_cols] == "blank"]
sample_only_cols <- sample_cols[samples$sample_type[samples$sample_id %in% sample_cols] == "sample"]

if (length(blank_cols) > 0 && length(sample_only_cols) > 0) {
    blank_means <- rowMeans(feat[, blank_cols, drop=FALSE], na.rm=TRUE)
    sample_means <- rowMeans(feat[, sample_only_cols, drop=FALSE], na.rm=TRUE)
    keep <- blank_means < opt$blank_cutoff * sample_means
    feat_filt <- feat[keep, ]
    cat("[FILTER] Blank filter: removed", sum(!keep), "features, kept", sum(keep), "\n")
} else {
    feat_filt <- feat
}

# QC CV filter
qc_cols <- sample_cols[samples$sample_type[samples$sample_id %in% sample_cols] == "QC"]
if (length(qc_cols) > 0) {
    qc_means <- rowMeans(feat_filt[, qc_cols, drop=FALSE], na.rm=TRUE)
    qc_sds <- rowSds(as.matrix(feat_filt[, qc_cols, drop=FALSE]), na.rm=TRUE)
    qc_cv <- qc_sds / (qc_means + 1e-10)
    keep_cv <- qc_cv < opt$qc_cv
    feat_filt <- feat_filt[keep_cv, ]
    cat("[FILTER] QC CV filter: kept", sum(keep_cv), "features\n")
}

write.csv(feat_filt, opt$outfiltered, row.names=FALSE)

# Normalization
feat_norm <- feat_filt
if (opt$norm == "TIC") {
    sample_only <- setdiff(colnames(feat_norm), meta_cols)
    tic <- colSums(feat_norm[, sample_only, drop=FALSE], na.rm=TRUE)
    norm_factor <- tic / mean(tic)
    for (s in sample_only) {
        feat_norm[, s] <- feat_norm[, s] / norm_factor[s]
    }
    cat("[NORM] TIC normalization applied\n")
} else if (opt$norm == "PQN") {
    cat("[NORM] PQN normalization not yet implemented, using none\n")
}

write.csv(feat_norm, opt$outnorm, row.names=FALSE)
cat("[FILTER] Done.\n")
