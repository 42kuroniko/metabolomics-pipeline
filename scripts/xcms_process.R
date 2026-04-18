#!/usr/bin/env Rscript
# XCMS peak picking, alignment, and grouping for LC-MS metabolomics data.

library(xcms)
library(optparse)

option_list <- list(
    make_option("--indir", type="character", help="Directory with mzML files"),
    make_option("--outcsv", type="character", help="Output feature table CSV"),
    make_option("--outrdata", type="character", help="Output xcms XCMSnExp object (.RData)"),
    make_option("--polarity", type="character", default="positive", help="Ionization mode [positive|negative]"),
    make_option("--ppm", type="integer", default=15, help="ppm tolerance"),
    make_option("--peakwidth", type="character", default="5,20", help="Peak width range (sec)"),
    make_option("--snthresh", type="integer", default=6, help="Signal-to-noise threshold")
)
opt <- parse_args(OptionParser(option_list=option_list))

pw <- as.numeric(strsplit(opt$peakwidth, ",")[[1]])

cat("[XCMS] Loading mzML files from:", opt$indir, "\n")
mzml_files <- list.files(opt$indir, pattern="\\.mzML$", full.names=TRUE)
if (length(mzml_files) == 0) {
    cat("[ERROR] No mzML files found in", opt$indir, "\n")
    quit(status=1)
}

cat("[XCMS] Found", length(mzml_files), "mzML files\n")

# OnDiskMSnExp
data_raw <- readMSData(mzml_files, mode="onDisk")

# Peak picking
cat("[XCMS] CentWave peak picking...\n")
cwp <- CentWaveParam(ppm=opt$ppm, peakwidth=pw, snthresh=opt$snthresh,
                     polarity=opt$polarity, integrate=1)

xdata <- findChromPeaks(data_raw, param=cwp)

# Alignment / grouping
cat("[XCAMS] Grouping features across samples...\n")
pd <- data.frame(sample=basename(mzml_files))
sampleNames(xdata) <- pd$sample

# Peak density grouping
pdp <- PeakDensityParam(sampleGroups=rep(1, length(mzml_files)),
                        minFrac=0.5, bw=0.5, binSize=0.05)

xdata <- groupChromPeaks(xdata, param=pdp)

# Fill peaks (missing data)
cat("[XCMS] Filling missing peaks...\n")
xdata <- adjustRtime(xdata, param=PeakGroupsParam(minDiff=1))

# Get feature table
cat("[XCMS] Extracting feature table...\n")
ft <- featureDefinitions(xdata)
ftab <- featureValues(xdata, value="into")

# Combine
out_df <- cbind(as.data.frame(ft[, c("mzmed", "rtmed", "rtmin", "rtmax")]),
                as.data.frame(ftab))
rownames(out_df) <- NULL

# Add m/z and RT in min
out_df$mz <- out_df$mzmed
out_df$rt_min <- out_df$rtmin / 60
out_df$rt_max <- out_df$rtmax / 60

write.csv(out_df, opt$outcsv, row.names=FALSE)
save(xdata, file=opt$outrdata)

cat("[XCMS] Done. Wrote:", opt$outcsv, "\n")
cat("[XCMS] Features detected:", nrow(out_df), "\n")
