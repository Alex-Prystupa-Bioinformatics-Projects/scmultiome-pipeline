# =============================================================================
# utils.R  |  QC Utility Functions
# =============================================================================
# Helper functions for generating and applying QC thresholds.
# Used by 04_qc.R (write thresholds) and 05_filter.R (apply thresholds).
# =============================================================================

# Write a qc_config.yml with default per-sample thresholds
# The user should review QC plots from 04_qc.R and edit this file
# before running 05_filter.R
create_qc_yaml_file <- function(seu_obj_list) {

    # 1. Build per-sample filter thresholds with sensible defaults
    qc_config <- lapply(names(seu_obj_list), function(sample) {
        list(
            filters = list(
                nFeature_RNA_min   = 250,
                nFeature_RNA_max   = 10000,
                nFeature_peaks_min = 200,
                nFeature_peaks_max = 25000,
                percent.mt         = 20,
                nucleosome_signal  = 2,
                TSS.enrichment_min = 2,
                nCount_peaks       = 50000
            )
        )
    })
    names(qc_config) <- names(seu_obj_list)

    # 2. Add global doublet filtering flag
    #    Default is FALSE — doublets are flagged but not removed
    #    Set to TRUE in qc_config.yml to remove predicted doublets in 05_filter.R
    qc_config[["remove_doublets"]] <- FALSE

    # 3. Write to YAML file
    write_yaml(qc_config, "configs/qc_config.yml")
}

# Filter each sample's Seurat object using thresholds from qc_config.yml
# Optionally removes predicted doublets if remove_doublets is set to TRUE
filter_seu_list_by_qc <- function(seu_obj_list, qc_configs) {

    # 1. Read global doublet filtering flag
    remove_doublets <- isTRUE(qc_configs$remove_doublets)

    # 2. Apply per-sample QC filters to each object
    setNames(lapply(names(seu_obj_list), function(sample) {
        seu_obj  <- seu_obj_list[[sample]]
        filters  <- qc_configs[[sample]]$filters
        n_before <- ncol(seu_obj)

        # Apply cell-level QC thresholds
        seu_obj <- subset(seu_obj, subset =
            nFeature_RNA      >= filters$nFeature_RNA_min    &
            nFeature_RNA      <  filters$nFeature_RNA_max    &
            nFeature_peaks    >= filters$nFeature_peaks_min  &
            nFeature_peaks    <  filters$nFeature_peaks_max  &
            percent.mt        <  filters$percent.mt           &
            nucleosome_signal <  filters$nucleosome_signal    &
            nCount_peaks      <  filters$nCount_peaks         &
            TSS.enrichment    >= filters$TSS.enrichment_min
        )

        # Optionally remove predicted doublets (scDblFinder.class == "doublet")
        if (remove_doublets) {
            n_before_doublet <- ncol(seu_obj)
            seu_obj          <- subset(seu_obj, scDblFinder.class == "singlet")
            message("  ", sample, ": removed ", n_before_doublet - ncol(seu_obj), " doublets")
        }

        message("  ", sample, ": ", n_before, " -> ", ncol(seu_obj), " cells retained")
        return(seu_obj)
    }), names(seu_obj_list))
}
