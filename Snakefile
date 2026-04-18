# =============================================================================
# metabolomics-pipeline: LC-MS Non-targeted Metabolomics
# Study: MTBLS2824 - Urine Metabolomics of Genitourinary Urothelial Cancer
# Author: Rong Wu
# =============================================================================

import os
import pandas as pd

configfile: "config.yaml"

PROJECT_DIR = os.path.abspath(workflow.basedir)
DATA_DIR = os.path.join(PROJECT_DIR, config.get("data_dir", "data"))
RESULTS = os.path.join(PROJECT_DIR, config.get("results_dir", "results"))
SAMPLE_SHEET = os.path.join(DATA_DIR, config.get("sample_sheet", "sample_sheet.csv"))

os.makedirs(RESULTS, exist_ok=True)

# =============================================================================
# Target outputs
# =============================================================================
rule all:
    input:
        f"{RESULTS}/05_stats/pca_plot.pdf",
        f"{RESULTS}/05_stats/plsda_plot.pdf",
        f"{RESULTS}/05_stats/volcano_plot.pdf",
        f"{RESULTS}/06_pathway/enrichment_table.csv",
        f"{RESULTS}/export/feature_table.tsv",
        f"{RESULTS}/export/annotations.tsv",

# -----------------------------------------------------------------------------
# 1. Convert raw Thermo/Bruker files to mzML (via msconvert)
# -----------------------------------------------------------------------------
rule convert_raw:
    output:
        directory(f"{RESULTS}/01_converted")
    log:
        f"{RESULTS}/logs/01_convert.log"
    params:
        raw_dir = config["raw_data_dir"],
        software = config.get("ms_convert", {}).get("software", "thermo"),
    shell:
        """
        mkdir -p {output}
        if [ "{params.software}" = "thermo" ]; then
            find {params.raw_dir} -name "*.raw" -o -name "*.d" | while read f; do
                wine msconvert "$f" --mzml --outdir {output} 2>> {log} || echo "msconvert failed: $f" >> {log}
            done
        fi
        echo "Convert done" >> {log}
        """

# -----------------------------------------------------------------------------
# 2. XCMS peak picking, alignment, grouping
# -----------------------------------------------------------------------------
rule xcms_process:
    input:
        f"{RESULTS}/01_converted"
    output:
        feature_table = f"{RESULTS}/02_xcms/features.csv",
        object = f"{RESULTS}/02_xcms/xset.RData"
    log:
        f"{RESULTS}/logs/02_xcms.log"
    params:
        polarity = config.get("xcms", {}).get("polarity", "positive"),
        ppm = config.get("xcms", {}).get("ppm", 15),
        peakwidth = config.get("xcms", {}).get("peakwidth", "5,20"),
        snthresh = config.get("xcms", {}).get("snthresh", 6),
    shell:
        """
        Rscript {PROJECT_DIR}/scripts/xcms_process.R \
            --indir {input} \
            --outcsv {output.feature_table} \
            --outrdata {output.object} \
            --polarity {params.polarity} \
            --ppm {params.ppm} \
            --peakwidth {params.peakwidth} \
            --snthresh {params.snthresh} \
            2>> {log}
        """

# -----------------------------------------------------------------------------
# 3. CAMERA annotation (adducts, isotopes, fragments)
# -----------------------------------------------------------------------------
rule camera_annotate:
    input:
        xset = f"{RESULTS}/02_xcms/xset.RData",
        feature_csv = f"{RESULTS}/02_xcms/features.csv"
    output:
        annotated_csv = f"{RESULTS}/03_annotated/annotated_features.csv",
        annotated_rdata = f"{RESULTS}/03_annotated/xset_annot.RData"
    log:
        f"{RESULTS}/logs/03_camera.log"
    params:
        polarity = config.get("xcms", {}).get("polarity", "positive"),
    shell:
        """
        Rscript {PROJECT_DIR}/scripts/camera_annotate.R \
            --xset {input.xset} \
            --incsv {input.feature_csv} \
            --outcsv {output.annotated_csv} \
            --outrdata {output.annotated_rdata} \
            --polarity {params.polarity} \
            2>> {log}
        """

# -----------------------------------------------------------------------------
# 4. Filter blanks, QC samples, normalization
# -----------------------------------------------------------------------------
rule filter_normalize:
    input:
        annotated_csv = f"{RESULTS}/03_annotated/annotated_features.csv",
        sample_sheet = SAMPLE_SHEET
    output:
        filtered_csv = f"{RESULTS}/04_filtered/filtered_features.csv",
        normalized_csv = f"{RESULTS}/04_filtered/normalized_features.csv"
    log:
        f"{RESULTS}/logs/04_filter.log"
    params:
        blank_cutoff = config.get("filter", {}).get("blank_cutoff", 0.25),
        qc_cv_cutoff = config.get("filter", {}).get("qc_cv_cutoff", 0.3),
        norm_method = config.get("filter", {}).get("norm_method", "TIC"),
    shell:
        """
        Rscript {PROJECT_DIR}/scripts/normalize_filter.R \
            --incsv {input.annotated_csv} \
            --samples {input.sample_sheet} \
            --outfiltered {output.filtered_csv} \
            --outnorm {output.normalized_csv} \
            --blank_cutoff {params.blank_cutoff} \
            --qc_cv {params.qc_cv_cutoff} \
            --norm {params.norm_method} \
            2>> {log}
        """

# -----------------------------------------------------------------------------
# 5. Statistical analysis (PCA, PLS-DA, volcano)
# -----------------------------------------------------------------------------
rule stats_analysis:
    input:
        normalized_csv = f"{RESULTS}/04_filtered/normalized_features.csv",
        sample_sheet = SAMPLE_SHEET
    output:
        pca_plot = f"{RESULTS}/05_stats/pca_plot.pdf",
        plsda_plot = f"{RESULTS}/05_stats/plsda_plot.pdf",
        volcano_plot = f"{RESULTS}/05_stats/volcano_plot.pdf",
        diff_features = f"{RESULTS}/05_stats/differential_features.csv"
    log:
        f"{RESULTS}/logs/05_stats.log"
    params:
        pvalue_thresh = config.get("stats", {}).get("pvalue_thresh", 0.05),
        fc_thresh = config.get("stats", {}).get("fc_thresh", 2.0),
        group_col = config.get("stats", {}).get("group_col", "group"),
        ncomp = config.get("stats", {}).get("ncomp", 2),
    shell:
        """
        Rscript {PROJECT_DIR}/scripts/stats_plot.R \
            --incsv {input.normalized_csv} \
            --samples {input.sample_sheet} \
            --outpca {output.pca_plot} \
            --outplsda {output.plsda_plot} \
            --outvolcano {output.volcano_plot} \
            --outdiff {output.diff_features} \
            --group_col {params.group_col} \
            --pval {params.pvalue_thresh} \
            --fc {params.fc_thresh} \
            --ncomp {params.ncomp} \
            2>> {log}
        """

# -----------------------------------------------------------------------------
# 6. Pathway enrichment (MetaboAnalystR / KEGG)
# -----------------------------------------------------------------------------
rule pathway_enrichment:
    input:
        diff_features = f"{RESULTS}/05_stats/differential_features.csv",
        normalized_csv = f"{RESULTS}/04_filtered/normalized_features.csv"
    output:
        enrichment_csv = f"{RESULTS}/06_pathway/enrichment_table.csv",
        pathway_pdf = f"{RESULTS}/06_pathway/pathway_plots.pdf"
    log:
        f"{RESULTS}/logs/06_pathway.log"
    params:
        pathway_db = config.get("pathway", {}).get("database", "hsa"),
    shell:
        """
        Rscript {PROJECT_DIR}/scripts/pathway_enrichment.R \
            --diff {input.diff_features} \
            --full {input.normalized_csv} \
            --outenrich {output.enrichment_csv} \
            --outpdf {output.pathway_pdf} \
            --db {params.pathway_db} \
            2>> {log}
        """

# -----------------------------------------------------------------------------
# 7. Export R-ready tables
# -----------------------------------------------------------------------------
rule export_r:
    input:
        normalized_csv = f"{RESULTS}/04_filtered/normalized_features.csv",
        annotated_csv = f"{RESULTS}/03_annotated/annotated_features.csv",
        diff_features = f"{RESULTS}/05_stats/differential_features.csv"
    output:
        feature_tsv = f"{RESULTS}/export/feature_table.tsv",
        annotations_tsv = f"{RESULTS}/export/annotations.tsv",
        diff_tsv = f"{RESULTS}/export/differential_features.tsv"
    shell:
        """
        mkdir -p {RESULTS}/export
        cp {input.normalized_csv} {output.feature_tsv}
        cp {input.annotated_csv} {output.annotations_tsv}
        cp {input.diff_features} {output.diff_tsv}
        """
