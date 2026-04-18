# metabolomics-pipeline

LC-MS non-targeted metabolomics analysis using XCMS + CAMERA + MetaboAnalystR.

## Study

**MTBLS2824** - Urine Metabolomics of Genitourinary Urothelial Cancer
National Taiwan University Hospital

Human urine samples. LC-QTOF, non-targeted metabolomics. Comparison of urothelial carcinoma patients vs healthy controls.

## Pipeline

```
raw (.d/.raw) -> msconvert (Thermo) -> mzML
mzML -> xcms (peak picking, alignment, grouping)
     -> CAMERA (adduct/isotope annotation)
     -> feature table
     -> blank filter / QC filter / normalization
     -> HMDB/KEGG annotation
     -> PCA / PLS-DA / volcano plot
     -> pathway enrichment
     -> R-ready export
```

## Setup

```bash
bash scripts/setup_env.sh
snakemake -n
snakemake --use-conda -j 4
```

## Output

```
results/
  01_converted/       mzML files
  02_xcms/            peak table, aligned features
  03_annotated/       CAMERA annotated features
  04_filtered/        filtered and normalized table
  05_stats/           PCA, PLS-DA, volcano plots
  06_pathway/         pathway enrichment
  export/             R/phyloseq-ready tables
```

## Citation

The Clinical Experiences of Urine Metabolomics of Genitourinary Urothelial Cancer in a Tertiary Hospital in Taiwan. MetaboLights MTBLS2824.
