# genome_utils.R
# Loads the correct genome libraries and returns annotation, peak.genome, and blacklist
# based on the species field in pipeline_config.yml.
# Source this file in any step that needs genome/annotation objects.

load_genome <- function(pipeline_config) {
    species <- pipeline_config$species

    if (species == "mouse") {
        library(EnsDb.Mmusculus.v79, quietly = TRUE)
        library(BSgenome.Mmusculus.UCSC.mm10, quietly = TRUE)

        annotation  <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
        seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

        peak.genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
        blacklist   <- Signac::blacklist_mm10
        seqlevelsStyle(blacklist) <- "NCBI"

        mt_pattern  <- "^mt-"

    } else if (species == "human") {
        library(EnsDb.Hsapiens.v86, quietly = TRUE)
        library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)

        annotation  <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
        seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

        peak.genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        blacklist   <- Signac::blacklist_hg38_unified

        mt_pattern  <- "^MT-"

    } else {
        stop("Unsupported species: '", species, "'. Set species to 'mouse' or 'human' in pipeline_config.yml.")
    }

    message("Loaded genome resources for species: ", species)
    return(list(
        annotation  = annotation,
        peak.genome = peak.genome,
        blacklist   = blacklist,
        mt_pattern  = mt_pattern
    ))
}
