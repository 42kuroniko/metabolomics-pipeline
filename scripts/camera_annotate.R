#!/usr/bin/env Rscript
# Author: Rong Wu
# CAMERA annotation: adducts, isotopes, fragments.

library(CAMERA)
library(optparse)

option_list <- list(
    make_option("--xset", type="character", help="XCMS XCMSnExp object (.RData)"),
    make_option("--incsv", type="character", help="Feature table from xcms_process"),
    make_option("--outcsv", type="character", help="Annotated feature CSV"),
    make_option("--outrdata", type="character", help="Annotated CAMERA object (.RData)"),
    make_option("--polarity", type="character", default="positive", help="Ionization mode")
)
opt <- parse_args(OptionParser(option_list=option_list))

load(opt$xset)

cat("[CAMERA] Annotating adducts and isotopes...\n")

# Convert to xsAnnotate
an <- annotate(xs=xdata, sample=NA)

# Group by retention time windows
an <- groupFWHM(an, sigma=6)
an <- findIsotopes(an, maxiso=3)
an <- findAdducts(an, polarity=opt$polarity)

# Get annotated results
feat_annot <- getAnnotatedFeatureTable(an)

if (nrow(feat_annot) > 0) {
    write.csv(feat_annot, opt$outcsv, row.names=FALSE)
} else {
    feat_df <- read.csv(opt$incsv)
    write.csv(feat_df, opt$outcsv, row.names=FALSE)
    cat("[CAMERA] No CAMERA annotations, keeping original features.\n")
}

save(an, file=opt$outrdata)
cat("[CAMERA] Done. Annotated features:", nrow(feat_annot), "\n")
