# combine atac_fragments files with the GEX folder
CombineDirectories <- function(data.dir, new.dir, sample){
   file.copy(from = paste0(data.dir, "/atac_fragments.tsv.gz"),
             to = paste0("data/raw/", sample, "/fragments.tsv.gz"))
   
   file.copy(from = paste0(data.dir, "/atac_fragments.tsv.gz.tbi"),
             to = paste0("data/raw/", sample, "/fragments.tsv.gz.tbi"))
   
   gex.files <- list.files(path = paste0(data.dir, "/filtered_feature_bc_matrix"))
   lapply(gex.files, 
          function(file.to.copy){
             file.copy(from = paste0(data.dir, "/filtered_feature_bc_matrix/", file.to.copy),
                       to = paste0("data/raw/", sample, "/", file.to.copy))
          })
}

# load in data and create assays for both Gene Expression and ATAC Peaks
CreateMultiomeSeurat <- function(data.dir, my.annotation = annotation, frag.path=NULL){
    require(Seurat)
    require(Signac)
    mat.list <- Read10X(data.dir)
    seu <- CreateSeuratObject(counts = mat.list[["Gene Expression"]])

    if(is.null(frag.path)){
        frag.path <- paste0(data.dir, "/fragments.tsv.gz")
    }

    # create ATAC assay and add it to the object
    seu[["ATAC"]] <- CreateChromatinAssay(
        counts = mat.list$Peaks,
        sep = c(":", "-"),
        fragments = frag.path,
        annotation = my.annotation
    )
    seu <- subset(seu, nFeature_RNA > 250)
    # if(seqlevelsStyle(seu[["ATAC"]]@ranges) != "UCSC"){
    #   seqlevelsStyle(seu[["ATAC"]]@ranges) <- "UCSC"
    #   rownames(seu[["ATAC"]]@counts) <- paste0("chr", rownames(seu[["ATAC"]]@counts))
    #   rownames(seu[["ATAC"]]@data) <- paste0("chr", rownames(seu[["ATAC"]]@data))
    # }
    return(seu)
}

# Call peaks using MACS2 and create a new "peaks" assay on the Seurat object.
# my.blacklist is required — pass genome$blacklist from genome_utils.R
CallMyPeaks <- function(seu, fragpath=NULL, grouping.var=NULL, my.macs2.path=NULL,
                        my.annotation=NULL, my.blacklist){
    require(Seurat)
    require(Signac)

    atac.assay <- ifelse("ATAC" %in% names(seu@assays), "ATAC", "peaks")
    DefaultAssay(seu) <- atac.assay

    # Set blacklist seqlevels style to match peak calls (NCBI format)
    seqlevelsStyle(my.blacklist) <- "NCBI"

    # New Code Below, was failing when grouping.var was NA
    if (!is.null(grouping.var)) {
      seu@meta.data[[grouping.var]] <- stringr::str_replace_all(seu@meta.data[[grouping.var]], " ", "_")
      seu@meta.data[[grouping.var]] <- stringr::str_replace_all(seu@meta.data[[grouping.var]], "[(]", "")
      seu@meta.data[[grouping.var]] <- stringr::str_replace_all(seu@meta.data[[grouping.var]], "[)]", "")
    }

    peaks <- CallPeaks(seu, 
      group.by = grouping.var,
      outdir = "./data/raw/macs-peaks",
      cleanup=FALSE,  
      macs2.path = my.macs2.path,
      verbose = TRUE
    )
    
    # dont want to lose the older peaks we had when peaks were called all together
    if(!is.null(grouping.var)){
      atac.assay <- ifelse("ATAC" %in% names(seu@assays), "ATAC", "peaks")
      old.peaks <- seu[[atac.assay]]@ranges
      frag.path <- seu@assays[[atac.assay]]@fragments[[1]]@path

      # combine the peaks into one range
      peaks <- GenomicRanges::reduce(c(peaks, old.peaks))
    }

    # Remove peaks on nonstandard chromosomes and in species-specific blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x = peaks, ranges = my.blacklist, invert = TRUE)

    # quantify counts in each peak
    macs2_counts <- FeatureMatrix(
      fragments = Fragments(seu),
      features = peaks,
      cells = colnames(seu)
    )

    # create a new assay using the MACS2 peak set and add it to the Seurat object
    seu[["peaks"]] <- 
    CreateChromatinAssay(
      counts = macs2_counts,
      fragments = Fragments(seu),
      annotation = my.annotation
    )
    return(seu)
}

# performing classical normalization and dimensionality reduction
Preprocess.and.Reduce.Dims <- function(seu, harmony=FALSE, 
                                        rna.pcs=30, atac.pcs=30, harmony.vars = NULL, 
                                        rna.theta = 0.5, atac.theta = 0.5,
                                        vars.to.regress=NULL, residual.features = NULL, 
                                        residual.peaks=NULL){
    require(Seurat)
    require(Signac)

    # RNA analysis
    DefaultAssay(seu) <- "RNA"
    seu <- SCTransform(seu, verbose = FALSE, residual.features = residual.features,
       vars.to.regress = vars.to.regress) %>% RunPCA(npc = rna.pcs, verbose = FALSE) 
    
    dims.to.use <- 1:rna.pcs
    if(harmony){
      require(harmony)
      message("Running harmony integration on RNA: ", harmony.vars)
      seu <- RunHarmony(seu, group.by.vars = harmony.vars, reduction = "pca", assay.use = "SCT", dims.use = 1:rna.pcs, reduction.save = "rna.harmony", project.dim = FALSE, theta = rep(rna.theta, length(harmony.vars)), max_iter = 25)
      dims.to.use <- 1:rna.pcs
    }
    reduction.to.use <- ifelse(harmony, "rna.harmony", "pca")
    seu <- seu %>% RunUMAP(dims = dims.to.use, reduction = reduction.to.use, reduction.name = "umap.rna", reduction.key = 'rnaUMAP_', verbose = FALSE)

    message("RNA analysis of ",seu@project.name," complete. Moving on to ATAC analysis")
    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    if("peaks" %in% names(seu@assays)){
      DefaultAssay(seu) <- "peaks"
    } else {
      DefaultAssay(seu) <- "ATAC"
    }
    seu <- RunTFIDF(seu, verbose = FALSE)
    seu <- FindTopFeatures(seu, min.cutoff = 'q50', verbose = FALSE)
    seu <- RunSVD(seu, n = atac.pcs, features = residual.peaks, verbose = FALSE)

    dims.to.use <- 1:atac.pcs
    if(harmony){
    message("Running harmony integration on ATAC: ", harmony.vars)
    # seu <- ScaleData(seu, assay = "ATAC")
    seu <- RunHarmony(seu, group.by.vars = harmony.vars, reduction = "lsi",assay.use = "peaks", dims.use = 2:atac.pcs, reduction.save = "atac.harmony", project.dim = FALSE, theta = rep(atac.theta, length(harmony.vars)), max_iter = 25)
    dims.to.use <- 1:(atac.pcs - 1)
    }
    reduction.to.use <- ifelse(harmony, "atac.harmony", "lsi")
    seu <- RunUMAP(seu, reduction = reduction.to.use, dims = dims.to.use, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = FALSE)

    message("ATAC analysis of ",seu@project.name," complete.")
    return(seu)
}

# constructing the WNN grapoh that incorporates both ATAC and RNA reduction
ConstructWNNGraph <- function(seu, harmony = FALSE, resolution = 0.8, rna.pcs=30, atac.pcs=30){
    require(Seurat)
    require(Signac)
  if(harmony){
    reduction_list <- list("rna.harmony", "atac.harmony")
    dims.to.use <- list(c(1:rna.pcs), c(1:(atac.pcs - 1)))
  } else {
    reduction_list <-  list("pca", "lsi")
    dims.to.use <- list(c(1:rna.pcs), c(2:atac.pcs))
  }
  message("Constructing WNN graph.")
  seu <- FindMultiModalNeighbors(seu, reduction.list = reduction_list, dims.list = dims.to.use, verbose = FALSE)
  seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", verbose = FALSE)
  seu <- FindClusters(seu, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)
  message("WNN graph complete.")
  return(seu)
}

# links peaks to genes to determine if there is a link between accessibility and expression
LinkMyPeaks <- function(seu, genes = NULL, peak.genome,distance.to.use=1000000){
    require(Seurat)
    require(Signac)
  DefaultAssay(seu) <- "peaks"
  seu <- RegionStats(seu, genome = peak.genome)
  # seu <- LinkPeaks(seu, peak.assay = "peaks", expression.assay = "RNA", distance = 5000, score_cutoff = 0.1)
  # seu <- FindVariableFeatures(seu, assay = "RNA")
  # seu <- LinkPeaks(seu, peak.assay = "peaks", expression.assay = "RNA", genes.use = seu@assays[["RNA"]]@var.features, distance = 10000)
  seu <- NormalizeData(seu, assay = "RNA")
  seu <- JoinLayers(seu, assay = "RNA") # Alex Added Code 
  seu <- 
    LinkPeaks(seu, 
      peak.assay = "peaks", 
      expression.assay = "RNA", 
      genes.use = genes,
      distance = distance.to.use)
  return(seu)
}

# EnrichR function
RunEnrichR <- function(genes.to.test, dbs=NULL, plot=T, title =NULL, top.n = 30){
    require(enrichR)
    require(ggplot2)
  if(is.null(dbs)){
    dbs <- c("GO_Biological_Process_2021" ,
             # "GO_Cellular_Component_2021",  "GO_Molecular_Function_2021", "MSigDB_Hallmark_2020",
             "KEGG_2021_Human",
             "Reactome_2022"
             )
  }
  res <- enrichr(genes = genes.to.test, databases = dbs) 
  Sys.sleep(1) # enrichR gets confused if too quick between different gene lists, must be an API thing on the internet
  res <- res[which(unlist(lapply(res, nrow)) > 0)] %>% dplyr::bind_rows(.id = "db")

  res <-  res %>%
        mutate(nGenes = sapply(strsplit(Overlap, split = "/"), getElement, 1)) %>%
               # db = ifelse(grepl("Biological", db), "BP",
               #             ifelse(grepl("Cellular", db), "CC", "MF"))) %>%
        dplyr::filter(nGenes > 1,
                      Adjusted.P.value < 0.05,
               Combined.Score > 0)
  
  if(nrow(res) == 0){
    warning("nrow == 0, skipping")
    return(NULL)
  }
  
  if(plot){
    p <- 
      res %>%
        group_by(db) %>%
        top_n(n = top.n, wt = -Adjusted.P.value) %>%
        arrange(-Adjusted.P.value) %>%
        mutate(sig = ifelse(Adjusted.P.value < 0.05, "*", "ns"),
               Term = factor(Term, levels = Term)) %>%
        ggplot(aes(x = -log(P.value), y = Term, fill = sig)) +
        geom_bar(stat = "identity") +
      geom_text(aes(label = Overlap)) +
        ggtitle(title) +
        theme_classic() +
        facet_grid(db~., scale = "free") +
        scale_fill_manual(values = c("*" = "#990000","ns" = "#999999")) +
        theme(axis.text = element_text(color = "black"),
              strip.text.y = element_text(color = "black", face = "bold"),
              plot.title = element_text(hjust = 1))
    return(p)
  } 
  return(res)
}


# keep standard chromosomes of a BSgenome object
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

# taken from https://stackoverflow.com/questions/25149520/append-filename-with-date 
date.time.append <- function(str, sep = '-', date.format ="%Y_%m_%d_%H_%M_%S") {
  stopifnot(is.character(str))
  return(paste(str, format(Sys.time(), date.format), sep = sep))  
}

# footprint motifs and specific peaks
FootprintMyPeaks <- function(obj, peaks.to.test, motifs.to.test=NULL, peak.genome=NULL, motif.set=NULL){
  require(motifmatchr)
  require(TFBSTools)
  require(JASPAR2020)

  if(is.null(motif.set)){
    motif.set <- getMatrixSet(
      x = JASPAR2020,
      opts = list(species = 9606, all_versions = TRUE)
    )
    # tf.names <- unlist(lapply(motif.set@listData, function(mini.list){return(mini.list@name)}), use.names = F)
    # tf.in.object <- tf.names[tf.names %in% rownames(obj[["RNA"]])]
    # motif.set[tf.names %in% tf.in.object]
  }

  if(is.null(peak.genome)){
    #peak.genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    peak.genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    seqlevelsStyle(peak.genome) <- "UCSC"
  }

  # # subset to retain only peaks in the standard chromosomes and remove non 'peaks' assays
  # # lose information but footprinting functionality is very memory intensive and this should help it
  # main.chroms <- standardChromosomes(peak.genome)
  # keep.peaks <- as.logical(seqnames(granges(obj)) %in% main.chroms)
  # obj <- obj[keep.peaks, ]

  # add motif information
  DefaultAssay(obj) <- "peaks"
  obj <- AddMotifs(obj, genome = peak.genome, pfm = motif.set)

  if(is.null(motifs.to.test)){
    enriched.motifs <- FindMotifs(
      object = obj,
      features = peaks.to.test
    )
    motifs.to.test <-
      enriched.motifs %>%
      dplyr::filter(p.adjust < 0.001) %>%
      top_n(n = 250, wt = -log(p.adjust)) %>%
      pull(motif) %>% 
      unique()
  }

  obj <- Footprint(
    object = obj,
    motif.name = motifs.to.test,
    genome = peak.genome, 
    in.peaks = TRUE, 
    compute.expected = FALSE,
    upstream = 1000,
    downstream = 1000
  )
  return(obj)
}

# Build a consensus peak set across all samples, recount each sample, and return
# an updated seu_obj_list ready to merge. Each object is stripped to RNA + new
# consensus peaks assay so all samples share the same peak feature space.
# ── DE Markers ────────────────────────────────────────────────────────────────

# Runs FindAllMarkers (Wilcoxon, positive markers only) and formats output in a
# Libra-like column structure. Returns a named list with three elements:
#   Raw-Markers     — all markers pre-filtering (12 cols)
#   Filtered-Markers — list of per-cluster data frames, columns renamed dynamically
#   Top-Markers     — pivot table of top N genes per cluster
custom_all_markers_function <- function(
  seu_obj,
  ident_col,
  assay = "RNA",
  p_adj_max = 0.05,
  lfc_min = 0.5,
  pct_min = 0.1,
  top_n = 100,
  background_name = "background"
) {

  # 1. Set identities
  Seurat::Idents(seu_obj) <- ident_col

  # 2. If identities are purely numeric-like, rename to cluster_<n> with proper numeric ordering
  ident_levels <- levels(Seurat::Idents(seu_obj))

  if (all(grepl("^[0-9]+$", ident_levels))) {
    numeric_levels <- sort(as.numeric(ident_levels))
    new_levels <- paste0("cluster_", numeric_levels)
    names(new_levels) <- as.character(numeric_levels)

    Seurat::Idents(seu_obj) <- factor(
      x = paste0("cluster_", as.character(Seurat::Idents(seu_obj))),
      levels = paste0("cluster_", numeric_levels)
    )
  }

  # 3. Pull normalized data layer for mean-expression columns
  expr_mat <- SeuratObject::LayerData(
    object = seu_obj,
    assay = assay,
    layer = "data"
  )

  # 4. Run FindAllMarkers (positive markers only)
  raw_markers <- Seurat::FindAllMarkers(
    object = seu_obj,
    assay = assay,
    only.pos = TRUE,
    logfc.threshold = 0,
    min.pct = 0,
    test.use = "wilcox"
  )

  # 5. Ensure gene column exists
  if (!"gene" %in% colnames(raw_markers)) {
    raw_markers <- tibble::rownames_to_column(raw_markers, var = "gene")
  }

  # 6. Rename columns to Libra-like structure
  raw_markers <- raw_markers %>%
    dplyr::rename(
      cell_type = cluster,
      avg_logFC = avg_log2FC,
      celltype_pct = pct.1,
      background_pct = pct.2,
      p_val_adj_BH = p_val_adj
    )

  # 7. Preserve proper numeric cluster ordering if applicable
  raw_cell_types <- unique(as.character(raw_markers$cell_type))

  if (all(grepl("^cluster_[0-9]+$", raw_cell_types))) {
    ordered_cell_types <- raw_cell_types[order(as.numeric(sub("^cluster_", "", raw_cell_types)))]
    raw_markers$cell_type <- factor(raw_markers$cell_type, levels = ordered_cell_types)
  }

  # 8. Compute mean normalized expression for cell type vs all other cells
  ident_vector <- as.character(Seurat::Idents(seu_obj))
  all_cells <- colnames(seu_obj)
  cell_types <- levels(Seurat::Idents(seu_obj))

  mean_exp_list <- lapply(cell_types, function(ct) {
    ct_cells <- all_cells[ident_vector == ct]
    bg_cells <- all_cells[ident_vector != ct]

    tibble::tibble(
      cell_type = ct,
      gene = rownames(expr_mat),
      celltype_exp = Matrix::rowMeans(expr_mat[, ct_cells, drop = FALSE]),
      background_exp = Matrix::rowMeans(expr_mat[, bg_cells, drop = FALSE])
    )
  })

  mean_exp_df <- dplyr::bind_rows(mean_exp_list)

  if (all(grepl("^cluster_[0-9]+$", unique(mean_exp_df$cell_type)))) {
    ordered_cell_types <- unique(as.character(mean_exp_df$cell_type))
    ordered_cell_types <- ordered_cell_types[order(as.numeric(sub("^cluster_", "", ordered_cell_types)))]
    mean_exp_df$cell_type <- factor(mean_exp_df$cell_type, levels = ordered_cell_types)
  }

  # 9. Join mean expression and add metadata / Bonferroni correction
  raw_markers <- raw_markers %>%
    dplyr::left_join(mean_exp_df, by = c("cell_type", "gene")) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::mutate(
      p_val_adj_Bonferroni = stats::p.adjust(p_val, method = "bonferroni")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      de_family = "singlecell",
      de_method = "wilcox"
    ) %>%
    dplyr::select(
      cell_type,
      gene,
      avg_logFC,
      celltype_pct,
      background_pct,
      celltype_exp,
      background_exp,
      p_val,
      p_val_adj_BH,
      p_val_adj_Bonferroni,
      de_family,
      de_method
    )

  # 10. Filter markers
  markers_filtered_up <- raw_markers %>%
    dplyr::filter(
      p_val_adj_Bonferroni < p_adj_max,
      avg_logFC > lfc_min,
      celltype_pct > pct_min
    )

  # 11. Split filtered markers by cell type and rename pct/exp columns dynamically
  filtered_list <- split(markers_filtered_up, markers_filtered_up$cell_type)

  filtered_list <- lapply(names(filtered_list), function(ct) {
    df <- filtered_list[[ct]]

    colnames(df)[colnames(df) == "celltype_pct"]    <- paste0(ct, ".pct")
    colnames(df)[colnames(df) == "background_pct"]  <- paste0(background_name, ".pct")
    colnames(df)[colnames(df) == "celltype_exp"]    <- paste0(ct, ".exp")
    colnames(df)[colnames(df) == "background_exp"]  <- paste0(background_name, ".exp")

    df
  })

  names(filtered_list) <- names(split(markers_filtered_up, markers_filtered_up$cell_type))

  # 12. Top markers table (pivot: one column per cluster, ranked by p_val_adj_Bonferroni then logFC)
  top_up_markers <- markers_filtered_up %>%
    dplyr::group_by(cell_type) %>%
    dplyr::arrange(p_val_adj_Bonferroni, dplyr::desc(avg_logFC), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::select(cell_type, gene, rank) %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = gene) %>%
    head(top_n)

  # 13. Return list
  list(
    "Raw-Markers"      = raw_markers,
    "Filtered-Markers" = filtered_list,
    "Top-Markers"      = top_up_markers
  )
}

# ── Consensus Peaks ───────────────────────────────────────────────────────────

create_consensus_peaks <- function(seu_obj_list, genome) {
    require(Seurat)
    require(Signac)
    require(GenomicRanges)

    # 1. Extract peak GRanges from each sample
    message("  Extracting peak ranges per sample...")
    peak_list <- lapply(seu_obj_list, function(seu_obj) {
        DefaultAssay(seu_obj) <- "peaks"
        granges(seu_obj[["peaks"]])
    })

    # 2. Build consensus peak set: union, width filter, standard chroms only
    message("  Building consensus peak set...")
    combined_peaks <- reduce(unlist(GRangesList(peak_list)))
    peak_widths    <- width(combined_peaks)
    combined_peaks <- combined_peaks[peak_widths > 20 & peak_widths < 10000]
    combined_peaks <- keepStandardChromosomes(combined_peaks, pruning.mode = "coarse")
    message("  Consensus peaks after filtering: ", length(combined_peaks))

    # 3. Per sample: recount with consensus peaks and rebuild peaks assay
    setNames(lapply(names(seu_obj_list), function(sample) {
        seu_obj <- seu_obj_list[[sample]]
        message("  Recounting consensus peaks for: ", sample)

        # Get fragment path and create fragment object scoped to filtered cells
        DefaultAssay(seu_obj) <- "peaks"
        frag_path <- seu_obj[["peaks"]]@fragments[[1]]@path
        frag_obj  <- CreateFragmentObject(path = frag_path, cells = colnames(seu_obj))

        # Recount cells against consensus peaks
        feat_matrix <- FeatureMatrix(
            fragments = frag_obj,
            features  = combined_peaks,
            cells     = colnames(seu_obj)
        )

        # Strip all assays except RNA, then add new consensus peaks assay
        DefaultAssay(seu_obj) <- "RNA"
        seu_obj <- DietSeurat(seu_obj, assays = "RNA")
        seu_obj[["peaks"]] <- CreateChromatinAssay(
            counts     = feat_matrix,
            fragments  = frag_obj,
            annotation = genome$annotation
        )

        return(seu_obj)
    }), names(seu_obj_list))
}
