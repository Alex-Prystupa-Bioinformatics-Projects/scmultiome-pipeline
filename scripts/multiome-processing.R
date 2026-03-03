#!/usr/bin/env Rscript

# This is the master R script that will mimick the R markdown tester file 'sc-mulit-organoid-psor-001.Rmd'
# it will only factor in a few arguments, depending on stage of pipline
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Run single-cell RNA + ATAC Multiomic Analysis from 10X Genomics platform")

# Add required command line arguments
p <- add_argument(p, "pipeline", 
                help="Comma delimted combinations of: init, create, callpeaks, qc, cluster, merge, linkpeaks", 
                type="character")
p <- add_argument(p, "samplesheet", help="samplesheet in csv format", type="character")

# Add optional command line flags
p <- add_argument(p, "--project_prefix", help="outfile name, no .RDS!", type="character", default="multiome")
p <- add_argument(p, "--grouping.var", help="grouping variable for peak calling algorithm", type = "character", default="NA")
p <- add_argument(p, "--RDS.file.in", help="RDS in-file for the pipeline desired", default="NA")
p <- add_argument(p, "--RunHarmony",flag=TRUE, help="Run Harmony batch correction")
p <- add_argument(p, "--footprint-peaks", help="atac peaks to run footprinting analysis on. txt file", type="character")
p <- add_argument(p, "--my.macs.path", help="path to macs environment. only used if callpaks pipeline is run", type="character")
#p <- add_argument(p, "--qc.sheet", help="path to csv file containing clusters to remove per sample", type="character")
p <- add_argument(p, "--SoupOrCellDF", help="path to csv file containing barcodes and their assigned samples", type="character", default="NA")
p <- add_argument(p, "--qc.split", flag=TRUE, help="whether or not to perform qc clustering on objects individually or combined")


# Parse the command line arguments
argv <- parse_args(p)

pipelines.to.run <- unlist(strsplit(argv$pipeline, split = ","))
samplesheet <- argv$samplesheet

library(future, quietly=TRUE)
library(future.apply, quietly=TRUE)

options(future.globals.maxSize = Inf)
options(future.rng.onMisuse = "ignore")

# parrallelize the processing
plan("multicore", workers = as.numeric(future::availableCores()))
plan()

typeof(argv$RunHarmony)
# .libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.2/", .libPaths()))
#library(EnsDb.Hsapiens.v86, quietly=TRUE)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Mmusculus.v79, quietly=TRUE) # mouse
library(BSgenome.Mmusculus.UCSC.mm10, quietly=TRUE) # mouse
library(Seurat, quietly=TRUE)
library(Signac, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(yaml, quietly=TRUE)
# library(enrichR, quietly=TRUE)

# source in wrapper functions
source("scripts/functions.R")
source('scripts/utils.R')

# being processing
# always read in the samplesheet for referemnce
samplesheet <- read.csv(samplesheet)
samples <- samplesheet$sampleName

# load in annotation files for peaks and ranges
#annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79) # mouse
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
peak.genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#peak.genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 # mouse
#seqlevels(peak.genome) <- paste0('chr', seqlevels(peak.genome)) ALEX COMMENTED THIS OUT INTENTIONALLY PEAK GENOME ALREADY HAD CHR attached to it

# load in data if desired
if(!(argv$RDS.file.in == "NA")){
    message("Loading in: ", argv$RDS.file.in)
    obj.list <- readRDS(argv$RDS.file.in)
    if(class(obj.list) != "list"){
        obj.list <- list(obj.list)
    }
}

# intitial run combines directories into local directory structure
if("init" %in% pipelines.to.run){
    message("Begining scMulitome processing")
    lapply(samples, 
            function(sample){
            data.dir <- samplesheet$path[samplesheet$sampleName == sample]
            print(data.dir)
            new.dir <- paste0("data/raw/", sample)
            dir.create(new.dir)
            CombineDirectories(data.dir, new.dir, sample)
    })
}

# pipeline to create the seurat objects for each instance in the samplesheet
if("create" %in% pipelines.to.run){
    obj.list <- 
        lapply(samples, 
            function(sample){
                # create objects
                seu <- CreateMultiomeSeurat(data.dir = paste0("data/raw/",sample))
                
                seu@project.name <- sample
                seu$orig.ident <- sample    
                # adding meta data columns to seurat object
                path.col <- grep("path", colnames(samplesheet))
                if(path.col < ncol(samplesheet)){
                    for(md.col.to.add in colnames(samplesheet)[c((path.col+1):ncol(samplesheet))]){
                        seu[[paste0(md.col.to.add)]] <- 
                            samplesheet[,md.col.to.add][samplesheet$sampleName == sample]
                    }
                }
                # assay-specific metrics
                #seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
                seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-") # mouse
                seu <- NucleosomeSignal(seu, assay = "ATAC")
                seu <- TSSEnrichment(seu, assay = "ATAC")
                return(seu)
            })
    names(obj.list) <- samples
    saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix, "-create-obj-list.RDS"))
    print("Done Creating")
}

# Splitting up object(s) by SouporCell called assignment
if(!(argv$SoupOrCellDF == "NA")){
    souporcelldf <- read.csv(argv$SoupOrCellDF)
    obj.list <-
        lapply(obj.list, 
        function(obj){
            souporcelldf <- souporcelldf[souporcelldf$barcode %in% colnames(obj), ] # This should be 100% for one object, less if one data frame is for multiple samplesheet instances
            cells.not.in.SoC <- colnames(obj)[!(colnames(obj) %in% souporcelldf$barcode)]
            residual.df <- data.frame(barcode = cells.not.in.SoC, assignment = "NA")
            souporcelldf <- rbind(souporcelldf, residual.df)
            rownames(souporcelldf) <- souporcelldf$barcode
            obj <- AddMetaData(obj, metadata = souporcelldf)
            obj <- subset(obj, assignment == "NA", invert = TRUE)

            # split object into multiple objects based on assignment
            if(argv$qc.split){
                message("Splitting object into ", length(unique(souporcelldf$assignment))-1, " objects")
                mini.obj.list <- SplitObject(obj, split = "assignment")
                mini.obj.list <- lapply(mini.obj.list,
                    function(mini.obj){
                            mini.obj@project.name <- paste0("sample ", mini.obj$assignment[1])
                            mini.obj$orig.ident <- paste0("control.skin.", mini.obj$assignment[1])
                            return(mini.obj)})
                return(mini.obj.list)
            }

            return(obj)
        })
    obj.list <- unlist(obj.list, recursive = FALSE)
    print("Done Soup or Cell")
}

# pipeline to call our own peaks using MACS2
if("callpeaks" %in% pipelines.to.run){
    # dir.create("./data/raw/macs-peaks")
    obj.list <- 
        lapply(obj.list, 
            function(obj){
                if(argv$grouping.var == "NA"){
                    # obj <- CallMyPeaks(obj)
                    obj <- CallMyPeaks(obj, my.macs2.path=argv$my.macs.path, my.annotation = annotation)
                } else {
                    # obj <- CallMyPeaks(obj, grouping.var=argv$grouping.var)
                    obj <- CallMyPeaks(obj, grouping.var=argv$grouping.var, my.macs2.path=argv$my.macs.path, my.annotation = annotation)
                } 
                return(obj)
            })
    saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-callpeaks-obj-list.RDS"))
    print("Done Call Peaks")
}

# quality control plots and clustering for each object individually
# no subsetting is done in this pipeline, must be done manually (for now)
if("qc" %in% pipelines.to.run){

    pdf(file = paste0("output/plots/", argv$project_prefix, "-qc-plots.pdf"),
        height = 8, width = 12)
    list.of.vars <- list("1" = c("nCount_RNA",  "nCount_peaks"),
                     "2" = c("nFeature_RNA","nFeature_peaks"),
                     "3" = "percent.mt",
                     "4" = c("nucleosome_signal" ,"TSS.enrichment"))
    p.list <-
        lapply(list.of.vars, function(vars.to.plot){
            p <- 
                bind_rows(lapply(obj.list, function(seu){return(seu@meta.data)})) %>%
                    dplyr::select(orig.ident, all_of(vars.to.plot)) %>%
                    reshape2::melt() %>%
                    ggplot(aes(x = orig.ident, y = value, fill = orig.ident)) +
                    geom_violin() +
                    facet_grid(~variable) +
                    # scale_fill_manual(values = c("#F7BE9F", "#EA7580")) +
                    theme_classic() +
                    theme(axis.text = element_text(color = "black"))
                    if(!("percent.mt" %in% vars.to.plot | "nucleosome_signal" %in% vars.to.plot)){
                        p <- p +scale_y_log10() 
                    }
        return(p)
    })
    print(p.list)

    density.scatter.list <- 
        lapply(obj.list, 
        function(seu){
            #p <- FeatureScatter(seu, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
	    p <- DensityScatter(seu, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
            return(p)
        })
    print(density.scatter.list)

    # preprocessing
    obj.list <- lapply(obj.list,
    function(seu){
        DefaultAssay(seu) <- "RNA"
        seu <- SCTransform(seu)
        seu <-   
            seu %>%
            RunPCA() %>%
            FindNeighbors() %>%
            FindClusters()
        return(seu)
    })

    p.list.2 <- 
        lapply(obj.list, 
            function(seu){
            p <- VlnPlot(seu, 
                pt.size = 0, 
                features = unlist(list.of.vars), 
                group.by = "seurat_clusters",
                stack = T, 
                fill.by = "ident", 
                log = T) + 
                NoLegend() +
                ggtitle(seu@project.name)
            return(p)
            })
    print(p.list.2)
    dev.off()

    saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-qc-obj-list.RDS"))

    # Create QC Yaml File
    create_qc_yaml_file(obj.list)
}

# filter objects based on qc output and user input
if ("filter" %in% pipelines.to.run) {

    # Load qc_configs
    qc_configs <- yaml::read_yaml("configs/qc_config.yml")

    # Filter by yaml file
    seu_list_filtered <- filter_seu_list_by_qc(obj.list, qc_configs)
    names(seu_list_filtered) <- names(obj.list)

    # Save Filtered List (used as input to merge step)
    saveRDS(seu_list_filtered, file = paste0("output/RDS-files/", argv$project_prefix,"-cluster-obj-list-filtered.RDS"))
}

# construct WNN graphs 
if("cluster" %in% pipelines.to.run){
    library(clustree, quietly=TRUE)
    message("Running Clustering Pipeline.")
    obj.list <- lapply(obj.list, function(obj){
        obj@project.name <- obj$cell.lineage[1]
        return(obj)
    })

    obj.list <- 
        lapply(obj.list, 
            function(obj){
                obj <- Preprocess.and.Reduce.Dims(obj, 
                                                harmony = argv$RunHarmony,
                                                harmony.vars = "orig.ident")
                obj <- ConstructWNNGraph(obj, 
                                        harmony = argv$RunHarmony,
                                        resolution = seq(0,1,0.1))
                return(obj)
            })
    saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-cluster-obj-list.RDS"))
    
    # plots to help decide resolution to use
    lapply(obj.list, function(obj){        
        pdf(file = paste0("output/plots/", argv$project_prefix, "-",obj@project.name, "-cluster-plots.pdf"),
            height = 8, width = 12)
            
        # print(clustree(obj@meta.data, prefix = "wsnn_res."))
        print(DimPlot(obj, reduction = "wnn.umap", group.by = paste0("wsnn_res.", seq(0.1, 1, 0.1)), label = T) & NoLegend())
        p1 <- DimPlot(obj, reduction = "pca", group.by = "wsnn_res.0.5", label = TRUE, label.size = 5, repel = FALSE) + ggtitle("PCA")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p2 <- DimPlot(obj, reduction = "umap.rna", group.by = "wsnn_res.0.5", label = TRUE, label.size = 5, repel = FALSE) + ggtitle("RNA")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p3 <- DimPlot(obj, reduction = "umap.atac", group.by = "wsnn_res.0.5", label = TRUE, label.size = 5, repel = FALSE) + ggtitle("ATAC") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p4 <- DimPlot(obj, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 5, repel = FALSE) + ggtitle("WNN") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))

        if(!(argv$SoupOrCellDF == "NA")){
            p5 <- DimPlot(obj, reduction = "pca", group.by = "assignment", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("PCA") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
            p6 <- DimPlot(obj, reduction = "umap.rna", group.by = "assignment", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("RNA") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
            p7 <- DimPlot(obj, reduction = "umap.atac", group.by = "assignment", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("ATAC") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
            p8 <- DimPlot(obj, reduction = "wnn.umap", group.by = "assignment", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("WNN")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
            print(
                ggpubr::ggarrange(
                    ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1),
                    ggpubr::ggarrange(p5, p6, p7, p8, ncol = 4, nrow = 1),
                    ncol = 1, nrow = 2)
            )
        } else {
            print(ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1))
        }
        dev.off()
    })

    # calculate marker genes at each resolution
    lapply(obj.list, 
        function(obj){
            M.list <- 
            lapply(seq(0.1, 1, 0.1), function(res){
                Idents(obj) <- paste0("wsnn_res.", res)
                obj <- PrepSCTFindMarkers(obj)
                M <- FindAllMarkers(obj, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25)
                M <- M %>% mutate(resolution = res)
                write.csv(M, file = paste0("output/tables/cluster-markers-res.", res, ".csv"))
                return(M) 
            })
            M <- bind_rows(M.list)
            write.csv(M, file = paste0("output/tables/cluster-markers-bound.csv"))
        })

    saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-cluster-obj-list.RDS"))
}

# perform sub-clustering analysis on a list of objects
if("subcluster" %in% pipelines.to.run){
    message("Running sub-clustering Pipeline.")
    obj.list <- SplitObject(obj.list[[1]], split.by = "cell.lineage")
    obj.list <- lapply(obj.list, function(obj){
        message("Subclustering: ", obj$cell.lineage[1])
        obj@project.name <- obj$cell.lineage[1]
        mini.obj.list <- SplitObject(obj, split.by = "orig.ident")
        for(i in 1:length(mini.obj.list)){
            message("Subclustering: ", mini.obj.list[[i]]$cell.lineage[1], "-", mini.obj.list[[i]]$orig.ident[1])
            mini.obj.list[[i]] <- subset(mini.obj.list[[i]], percent.mt < 30 & nFeature_peaks >= 1000 & nFeature_SCT >= 500)
            mini.obj.list[[i]][["peaks"]]@fragments <- mini.obj.list[[i]][["peaks"]]@fragments[i]
            mini.obj.list[[i]] <- RunTFIDF(mini.obj.list[[i]], assay = "peaks")
            mini.obj.list[[i]] <- FindTopFeatures(mini.obj.list[[i]], min.cutoff = 'q50', assay = "peaks", verbose = FALSE)
        }
        var.peaks <- lapply(
            mini.obj.list, function(mini.obj){
            which.variable <- rownames(mini.obj[["peaks"]]) %in% VariableFeatures(mini.obj)
            return(mini.obj[["peaks"]]@ranges[which.variable, ])
        })
        residual.peaks <- lapply(var.peaks, paste0)
        residual.peaks.to.use <- names(table(unlist(residual.peaks)))[table(unlist(residual.peaks)) == 6]

        # residual.peaks <- GenomicRanges::reduce(unlist(GenomicRanges::GRangesList(var.peaks)))
        # residual.peaks.string <- paste0(residual.peaks)
        residual.peaks.string <- stringr::str_replace(residual.peaks.to.use, pattern = ":", replacement = "-")

        residual.features <- lapply(mini.obj.list,
                                    function(mini.obj){
                                    # mini.obj <- subset(mini.obj, percent.mt < 30 & nFeature_peaks >= 1000 & nFeature_SCT >= 500)
                                    DefaultAssay(mini.obj) <- "RNA"
                                    mini.obj <- SCTransform(mini.obj) 
                                    res.feats <- VariableFeatures(mini.obj)
                                    return(res.feats)}
                                    )
        residual.features.to.use <- names(table(unlist(residual.features)))[table(unlist(residual.features)) >= 5]

        merged.obj <- Preprocess.and.Reduce.Dims(obj, 
                harmony = TRUE,
                harmony.vars = c("insitution", "orig.ident"),
                residual.features = residual.features.to.use,
                rna.pcs = 30, 
                atac.pcs = 30,
                rna.theta = 10,
                atac.theta = 10,
                residual.peaks =  residual.peaks.string,
                vars.to.regress = c("insitution", "orig.ident"))

        merged.obj <- ConstructWNNGraph(merged.obj, 
                                    harmony = TRUE, 
                                    rna.pcs = 30, 
                                    atac.pcs = 30,
                                    resolution = 0.4)

        message("Completed subclustering of ", merged.obj$cell.lineage[1])
        return(merged.obj)
    })
    saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-subcluster-obj-list.RDS"))
}

# merge objects and create conserved peaks across objects
if("merge" %in% pipelines.to.run){
    message("Running Merging Pipeline")

    # Load Filtered RDS File if it exists
    file_path <- paste0("output/RDS-files/", argv$project_prefix, "-cluster-obj-list-filtered.RDS")

    if (file.exists(file_path)) {
        obj.list <- readRDS(file_path)
        message("Loaded existing filtered object list from file.")
        message(file_path)
    } else {
        message("Using unfiltered object for merge")
    }

    # Run Downstream Merge Pipeline
    peak.list <- lapply(obj.list, function(seu){
        assay.to.use <- "ATAC"
        if(!(assay.to.use %in% names(seu@assays))){
            assay.to.use <- "peaks"
        }
        peak.granges <- seu[[assay.to.use]]@ranges
        return(peak.granges)
    })

    # intersecting the peak list
    combined.peaks <- reduce(unlist(GRangesList(peak.list)))
    
    # Filter out bad peaks based on length
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
    combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
    
    frag.paths <- lapply(obj.list, function(seu){
        assay.to.use <- "ATAC"
        if(!(assay.to.use %in% names(seu@assays))){
            assay.to.use <- "peaks"
        }
        frag.path <- seu@assays[[assay.to.use]]@fragments[[1]]@path
        return(frag.path)
    })

    # done post-qc so we know which cells to keep already
    frag.list <- lapply(names(frag.paths),function(f.path.name){
        f.path <- frag.paths[[f.path.name]]
        cells.to.keep <- colnames(obj.list[[f.path.name]])
        frag.obj <- 
            CreateFragmentObject(path = f.path, cells = cells.to.keep)
        return(frag.obj)
    })
    names(frag.list) <- names(frag.paths)

    # feature-matrix generation for peaks assay
    feat.mat.list <- lapply(names(frag.paths),function(f.path.name){
        f.path <- frag.paths[[f.path.name]]
        cells.to.keep <- colnames(obj.list[[f.path.name]])
        frag.obj <- frag.list[[f.path.name]]
        frag.matrix <- 
            FeatureMatrix(fragments = frag.obj,
                        features = combined.peaks,
                        cells = colnames(obj.list[[f.path.name]]))
        return(frag.matrix)
    })
    names(feat.mat.list) <- names(frag.paths)
    
    # re-create the peaks assay using the new counts
    obj.list <- lapply(names(frag.paths),function(f.path.name){
        seu <- obj.list[[f.path.name]]
        DefaultAssay(seu) <- "RNA"
        seu <- DietSeurat(seu, assays = c("RNA", "SCT"))
        seu[["peaks"]] <- 
            CreateChromatinAssay(
                    counts = feat.mat.list[[f.path.name]], 
                    fragments = frag.list[[f.path.name]], 
                    annotation = annotation)
        return(seu)
    })

    # merge the objects now with a combined peaks set
    # if(length(obj.list) == 2){
    #     merged.obj <- merge(obj.list[[1]], obj.list[[2]])
    # } else {
    #     merged.obj <- merge(obj.list[[1]], obj.list[c(2:length(obj.list))])
    # }
    message("Length of obj.list: ", length(obj.list))
    merged.obj <- merge(obj.list[[1]], obj.list[-1])
    
    # create union of variable genes
    residual.features.to.use <- Reduce(intersect, lapply(obj.list, function(obj){return(obj[["SCT"]]@var.features)}))

    # # create union of variable peaks
    # top.feats.list <- lapply(obj.list, function(seu){
    #     DefaultAssay(seu) <- "peaks"
    #     seu <- RunTFIDF(seu)
    #     seu <- FindTopFeatures(seu, min.cutoff = "q50")
    #     var.peaks <- seu[["peaks"]]@var.features
    #     return(var.peaks)
    # })
    # residual.peaks.to.use <- Reduce(intersect, top.feats.list)
    # message("There are: ", length(residual.peaks.to.use), " residual peaks found.")
    residual.peaks.to.use <- NULL
    
    # rerun the clustering on the merged object
    if(argv$RunHarmony){
        message("RunHarmony flagged as TRUE")
        my.harmony.vars <- "orig.ident"
    } else {
        message("RunHarmony flagged as not TRUE")
        my.harmony.vars <- NULL
    }
    merged.obj <- Preprocess.and.Reduce.Dims(merged.obj, 
        harmony = argv$RunHarmony,
        harmony.vars = my.harmony.vars,
        residual.features = residual.features.to.use,
        residual.peaks =  residual.peaks.to.use,
        vars.to.regress = NULL)
    merged.obj <- ConstructWNNGraph(merged.obj, 
                                    harmony = argv$RunHarmony, 
                                    resolution = seq(0,1,0.1))
    merged.obj <- list(merged.obj)
    saveRDS(merged.obj, file = paste0("output/RDS-files/", argv$project_prefix,"-merged-obj-list-improved.RDS"))
}

# link peaks to genes
if("linkpeaks" %in% pipelines.to.run){
    future::plan("sequential")
    obj.list <- 
        lapply(obj.list,
            function(obj){
                    #genes.to.link <- c("B2m", "Fst", "H2-K1", "Sst", "Apoe", "Calca", "Tac1", "Apod", "C1qa", "Adgre1", "Ccr2", "Ccr5")
                    genes.to.link <- NULL
                    obj <- LinkMyPeaks(obj,
                                        genes = genes.to.link,
                                        distance.to.use = 250001,
                                        peak.genome = peak.genome)
                    gc()
                    return(obj)
            })
    saveRDS(obj.list, file = paste0("output/RDS-files/",argv$project_prefix, "-peak-gene-linked-obj-list.RDS"))
}

# footprinting TF activity by motifs
if("footprint" %in% pipelines.to.run){
    library(motifmatchr, quietly=TRUE)
    library(TFBSTools, quietly=TRUE)
    library(JASPAR2020, quietly=TRUE)

    peaks.to.footprint <- readr::read_delim(argv$footprint-peaks, rownames = F)
    footprint.list <- 
        lapply(obj.list,
            function(obj){
                footprint.obj <- 
                    FootprintMyPeaks(obj, 
                        peaks.to.test = peaks.to.footprint,
                        peak.genome = peak.genome)
            })
    saveRDS(obfootprintj.list, file = paste0(argv$outfilename, "-footprinted-obj-list.RDS"))
}

message("disco complete.")
