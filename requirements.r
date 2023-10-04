
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17")

bio_pkgs <- c("DESeq2", "BiocParallel", "limma", 
                "WGCNA", "clusterProfiler", "org.Hs.eg.db", 
                "enrichplot", "biomaRt", "gcrma", "oligo")

R_pkgs <- c("data.table", "dplyr", "devtools")

# install:
BiocManager::install(bio_pkgs)

# load all at once
invisible(lapply(bio_pkgs, function(x) library(x, character.only=TRUE)))
install.packages(R_pkgs)
invisible(lapply(R_pkgs, function(x) library(x, character.only=TRUE)))