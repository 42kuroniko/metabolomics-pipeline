#!/usr/bin/env Rscript
# Statistical analysis: PCA, PLS-DA, volcano plot.

library(optparse)
library(ggplot2)
library(mixOmics)

option_list <- list(
    make_option("--incsv", type="character", help="Normalized feature table CSV"),
    make_option("--samples", type="character", help="Sample sheet CSV"),
    make_option("--outpca", type="character", help="PCA plot PDF"),
    make_option("--outplsda", type="character", help="PLS-DA plot PDF"),
    make_option("--outvolcano", type="character", help="Volcano plot PDF"),
    make_option("--outdiff", type="character", help="Differential features CSV"),
    make_option("--group_col", type="character", default="group", help="Group column name"),
    make_option("--pval", type="double", default=0.05, help="P-value threshold"),
    make_option("--fc", type="double", default=2.0, help="Fold change threshold"),
    make_option("--ncomp", type="integer", default=2, help="Number of components")
)
opt <- parse_args(OptionParser(option_list=option_list))

feat <- read.csv(opt$incsv)
samples <- read.csv(opt$samples)

# Prepare matrix
meta_cols <- c("mz", "rt_min", "rt_max", "mzmed", "rtmed", "rtmin", "rtmax",
               "isotopes", "adduct", "pcgroup", "neutral_mass", "charge")
sample_cols <- intersect(colnames(feat), samples$sample_id)
mat <- as.matrix(feat[, sample_cols, drop=FALSE])
mat[is.na(mat)] <- 0

# Match sample order
samples <- samples[match(sample_cols, samples$sample_id), ]
groups <- samples[, opt$group_col]

# --- PCA ---
cat("[STATS] Running PCA...\n")
pca_res <- pca(mat, ncomp=opt$ncomp, center=TRUE, scale=TRUE)

df_pca <- data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], Group=groups)
p_pca <- ggplot(df_pca, aes(x=PC1, y=PC2, color=Group)) +
    geom_point(size=3) + stat_ellipse() +
    labs(x=paste0("PC1 (", round(pca_res$vip[1]*100), "%)"),
         y=paste0("PC2 (", round(pca_res$vip[2]*100), "%)"),
         title="PCA") + theme_bw()

pdf(opt$outpca, width=6, height=5)
print(p_pca)
dev.off()

# --- PLS-DA ---
cat("[STATS] Running PLS-DA...\n")
plsda_res <- plsda(mat, groups, ncomp=opt$ncomp)

pdf(opt$outplsda, width=6, height=5)
plotIndiv(plsda_res, comp=c(1,2), group=groups, legend=TRUE,
          title="PLS-DA", ellipse=TRUE)
dev.off()

# --- Volcano plot ---
cat("[STATS] Running volcano plot...\n")
groups_u <- unique(groups)
if (length(groups_u) >= 2) {
    g1 <- groups_u[1]; g2 <- groups_u[2]
    idx1 <- which(groups == g1)
    idx2 <- which(groups == g2)

    mean1 <- rowMeans(mat[, idx1, drop=FALSE])
    mean2 <- rowMeans(mat[, idx2, drop=FALSE])
    fc <- mean2 / (mean1 + 1e-10)
    log2fc <- log2(fc)

    # Welch's t-test
    t_stat <- sapply(1:nrow(mat), function(i) {
        x <- mat[i, idx1]; y <- mat[i, idx2]
        t.test(x, y, var.equal=FALSE)$statistic
    })
    pvals <- pt(abs(t_stat), df=length(c(idx1,idx2))-2, lower.tail=FALSE) * 2
    pvals[pvals < 1e-100] <- 1e-100
    logp <- -log10(pvals)

    volcano_df <- data.frame(mz=feat$mz[1:nrow(mat)], rt=feat$rt_min[1:nrow(mat)],
                             log2FC=log2fc, logP=logp,
                             significance=ifelse(pvals < opt$pval & abs(log2fc) > log2(opt$fc),
                                                 ifelse(log2fc > 0, "Up", "Down"), "NS"))

    p_vol <- ggplot(volcano_df, aes(x=log2FC, y=logP, color=significance)) +
        geom_point(size=1.5) +
        scale_color_manual(values=c("Up"="red", "Down"="blue", "NS"="grey")) +
        geom_hline(yintercept=-log10(opt$pval), lty=2) +
        geom_vline(xintercept=c(-log2(opt$fc), log2(opt$fc)), lty=2) +
        labs(title="Volcano Plot", x="log2 Fold Change", y="-log10 P-value") + theme_bw()

    pdf(opt$outvolcano, width=6, height=5)
    print(p_vol)
    dev.off()

    # Save differential features
    diff_idx <- which(pvals < opt$pval & abs(log2fc) > log2(opt$fc))
    diff_feat <- cbind(feat[diff_idx, c("mz", "rt_min", "rt_max")],
                       log2FC=log2fc[diff_idx], Pvalue=pvals[diff_idx])
    write.csv(diff_feat, opt$outdiff, row.names=FALSE)
    cat("[STATS] Differential features:", length(diff_idx), "\n")
} else {
    cat("[STATS] Less than 2 groups, skipping volcano\n")
    file.create(opt$outvolcano)
    file.create(opt$outdiff)
}

cat("[STATS] Done.\n")
