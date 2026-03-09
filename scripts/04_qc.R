#!/usr/bin/env Rscript
# =============================================================================
# 04_qc.R  |  Pipeline Step 4: QC Plots and Threshold Configuration
# =============================================================================
# Generates QC diagnostic plots for each sample and writes a qc_config.yml
# with default thresholds for the user to review and edit before filtering.
#
# No normalization or clustering is performed here — objects are raw counts
# with QC metrics attached (from steps 2 and 3).
#
# Merge strategy: DietSeurat strips all assays except RNA before merging.
# All QC metrics live in @meta.data so no count matrices are needed.
# This makes the merge fast even for large multiome objects.
#
# Plots generated (one page per metric, saved individually and in combined PDF):
#   violins/   — VlnPlot per metric with red cutoff lines + mean summary table
#   density/   — DensityScatter per sample: nCount_peaks vs TSS.enrichment
#   doublets/  — scDblFinder.score violin across all samples
#
# Output tree:
#   output/plots/{prefix}-qc/
#   ├── {prefix}-QC.pdf          <- all plots combined
#   ├── violins/QC-{metric}.pdf
#   ├── density/QC-density-{sample}.pdf
#   └── doublets/QC-doublets.pdf
#
# After reviewing the plots, edit configs/qc_config.yml to set your
# per-sample thresholds, then run 05_filter.R.
#
# Usage:
#   Rscript scripts/04_qc.R configs/samplesheet.csv \
#       --project_prefix myproject \
#       --RDS_file_in output/RDS-files/myproject-03-callpeaks-obj-list.RDS
#
# Input:  output of 03_callpeaks.R  ({prefix}-03-callpeaks-obj-list.RDS)
# Output: output/plots/{prefix}-qc/{prefix}-QC.pdf
#         output/plots/{prefix}-qc/violins/
#         output/plots/{prefix}-qc/density/
#         output/plots/{prefix}-qc/doublets/
#         output/tables/{prefix}-doublet-summary.csv
#         configs/qc_config.yml
# =============================================================================

library(argparser,  quietly = TRUE)
library(yaml,       quietly = TRUE)
library(Seurat,     quietly = TRUE)
library(Signac,     quietly = TRUE)
library(dplyr,      quietly = TRUE)
library(ggplot2,    quietly = TRUE)
library(gridExtra,  quietly = TRUE)
library(patchwork,  quietly = TRUE)
library(glue,       quietly = TRUE)

# 1. Parse arguments
p <- arg_parser("Step 4: Generate QC plots and write qc_config.yml")
p <- add_argument(p, "samplesheet",
                  help = "Path to samplesheet CSV", type = "character")
p <- add_argument(p, "--pipeline_config",
                  help    = "Path to pipeline_config.yml",
                  default = "configs/pipeline_config.yml", type = "character")
p <- add_argument(p, "--project_prefix",
                  help    = "Output file prefix (no .RDS extension)",
                  default = "multiome", type = "character")
p <- add_argument(p, "--RDS_file_in",
                  help = "Input RDS file — output of 03_callpeaks.R", type = "character")
argv <- parse_args(p)

if (is.na(argv$RDS_file_in) || argv$RDS_file_in == "NA") {
    stop("--RDS_file_in is required. Provide the path to output of 03_callpeaks.R.")
}

# 2. Load config and helper functions
pipeline_config <- yaml::read_yaml(argv$pipeline_config)
samplesheet     <- read.csv(argv$samplesheet)

source("scripts/utils.R")

# 3. Load input objects from previous step
message("Starting Step 4: qc")
message("Loading: ", argv$RDS_file_in)
seu_obj_list <- readRDS(argv$RDS_file_in)

# 4. DietSeurat (RNA only) then merge for multi-sample VlnPlot
#    Strips peaks assay matrices before merge — all QC metrics live in @meta.data
#    so count matrices are never needed. Keeps merge fast on large multiome objects.
message("  Slimming objects with DietSeurat before merge...")
seu_obj_list_diet <- lapply(seu_obj_list, function(seu_obj) {
    DietSeurat(seu_obj, assays = "RNA")
})
seu_obj_merged <- merge(seu_obj_list_diet[[1]], seu_obj_list_diet[-1])
message("  Merge complete.")

# 5. Create output directory tree
qc_dir      <- file.path("output/plots", paste0(argv$project_prefix, "-qc"))
violin_dir  <- file.path(qc_dir, "violins")
density_dir <- file.path(qc_dir, "density")
doublet_dir <- file.path(qc_dir, "doublets")

for (d in c(violin_dir, density_dir, doublet_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# 6. Define default QC cutoffs
#    These mirror the defaults written to qc_config.yml by create_qc_yaml_file()
#    and are shown as red dotted lines on the violin plots.
#    min = lower bound (keep cells above), max = upper bound (keep cells below)
qc_cutoffs <- list(
    nFeature_RNA      = c(min = 250,   max = 10000),
    nFeature_peaks    = c(min = 200,   max = 25000),
    percent.mt        = c(max = 20),
    nucleosome_signal = c(max = 2),
    TSS.enrichment    = c(min = 2),
    nCount_peaks      = c(max = 50000)
)

# 7. Open combined PDF — all plots appended here as well as saved individually
combined_pdf_path <- file.path(qc_dir, paste0(argv$project_prefix, "-QC.pdf"))
pdf(combined_pdf_path, height = 8, width = 16)

# 8. Violin plots per QC metric
#    - One page per metric in combined PDF + individual PDF saved to violins/
#    - VlnPlot grouped by orig.ident with red dotted lines at default cutoffs
#    - tableGrob mean-per-sample summary table alongside violin (3:1 width ratio)
for (metric in names(qc_cutoffs)) {
    cut     <- qc_cutoffs[[metric]]
    has_min <- !is.na(cut["min"])
    has_max <- !is.na(cut["max"])

    # ---- build cutoff label for plot title
    cutoff_label <- dplyr::case_when(
        has_min & has_max ~ glue("QC: {cut['min']} – {cut['max']}"),
        has_min           ~ glue("QC ≥ {cut['min']}"),
        has_max           ~ glue("QC ≤ {cut['max']}"),
        TRUE              ~ ""
    )

    # ---- mean per sample summary table
    summary_df <- seu_obj_merged@meta.data %>%
        group_by(orig.ident) %>%
        summarise(mean = round(mean(.data[[metric]], na.rm = TRUE), 2), .groups = "drop")

    table_grob <- tableGrob(summary_df, rows = NULL)

    # ---- violin plot with cutoff lines
    p <- VlnPlot(seu_obj_merged,
                 features = metric,
                 group.by = "orig.ident",
                 pt.size  = 0) +
        NoLegend() +
        ggtitle(paste0(metric, "  (", cutoff_label, ")"))

    if (has_min) {
        p <- p + geom_hline(yintercept = cut["min"], color = "red",
                            linetype = "dotted", linewidth = 0.8)
    }
    if (has_max) {
        p <- p + geom_hline(yintercept = cut["max"], color = "red",
                            linetype = "dotted", linewidth = 0.8)
    }

    # ---- compose final layout: violin (3) + table (1)
    final_plot <- p + wrap_elements(table_grob) +
        plot_layout(ncol = 2, widths = c(3, 1))

    # ---- print to combined PDF
    print(final_plot)

    # ---- save individual PDF to violins/
    pdf(file.path(violin_dir, glue("QC-{metric}.pdf")), height = 8, width = 16)
    print(final_plot)
    dev.off()
}

# 9. Density scatter per sample: nCount_peaks vs TSS.enrichment
#    - One page per sample in combined PDF + individual PDFs saved to density/
#    - Quantile lines help identify natural cutoffs for each sample independently
lapply(names(seu_obj_list), function(sample) {
    p <- DensityScatter(seu_obj_list[[sample]],
                        x = "nCount_peaks", y = "TSS.enrichment",
                        log_x = TRUE, quantiles = TRUE) +
        ggtitle(sample)

    print(p)

    pdf(file.path(density_dir, glue("QC-density-{sample}.pdf")), height = 8, width = 10)
    print(p)
    dev.off()
})

# 10. Doublet score violin across all samples
#     - scDblFinder.score: 0 = confident singlet, 1 = confident doublet
doublet_plot <- VlnPlot(seu_obj_merged,
                        features = "scDblFinder.score",
                        group.by = "orig.ident",
                        pt.size  = 0) +
    NoLegend() +
    ggtitle("Doublet Score by Sample (scDblFinder)")

print(doublet_plot)

pdf(file.path(doublet_dir, "QC-doublets.pdf"), height = 8, width = 12)
print(doublet_plot)
dev.off()

# Close combined PDF
dev.off()
message("  QC plots saved: ", combined_pdf_path)

# 11. Save doublet summary CSV: total cells, doublet count, and % per sample
doublet_summary <- bind_rows(lapply(seu_obj_list, function(seu_obj) {
    n_total    <- ncol(seu_obj)
    n_doublets <- sum(seu_obj$scDblFinder.class == "doublet")
    data.frame(
        sample       = seu_obj@project.name,
        total_cells  = n_total,
        n_doublets   = n_doublets,
        pct_doublets = round(n_doublets / n_total * 100, 2)
    )
}))

doublet_csv_path <- paste0("output/tables/", argv$project_prefix, "-doublet-summary.csv")
write.csv(doublet_summary, file = doublet_csv_path, row.names = FALSE)
message("  Doublet summary saved: ", doublet_csv_path)

# 12. Write qc_config.yml with default thresholds for each sample
#     Review the plots above, then edit configs/qc_config.yml before running 05_filter.R
#     Note: remove_doublets defaults to FALSE — set to TRUE to remove predicted doublets
create_qc_yaml_file(seu_obj_list)
message("  Default thresholds written to: configs/qc_config.yml")
message("  --> Review QC plots and edit configs/qc_config.yml before running 05_filter.R")

message("Step 4 complete.")
