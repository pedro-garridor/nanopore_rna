# Nanopore RNA-Seq workflow
# Copyright (C) 2024, Pedro Garrido Rodríguez

# Libraries
suppressPackageStartupMessages(library(bambu))

# Arguments
args <- commandArgs(trailingOnly = TRUE)
ncores <- args[1]
out.dir <- 'results/bambu/'
ref.genome <- as.character(args[2])
ref.tx <- as.character(args[3])
bams <- unlist(lapply(args[-c(1:3)], as.character))
bais <- gsub('$', '\\.bai', bams)

# Import libraries and data
hg38.annot <- prepareAnnotations(ref.tx)
BAM <- Rsamtools::BamFileList(bams,
                              index = bais)


# Run bambu
livers <- bambu(reads = BAM,
                genome = ref.genome,
                annotations = hg38.annot,
                ncore = ncores,
                verbose = F)

# Output bambu results
writeBambuOutput(livers, path = out.dir)


# Session info
save.image(paste(out.dir, 'bambu.RData', sep = '/'))
sessionInfo()
