if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
## These packages should be downloaded from Bioconductor
## Uncomment and run these lines
## BiocManager::install("IRanges")
## BiocManager::install("BiocParallel")
## Do the same with any other missing package, e.g.,
## BiocManager::install("haven")

suppressPackageStartupMessages({
    library(haven)
    library(readr)
    library(dplyr)
    library(IRanges)
    library(stringdist)
    library(purrr)
    library(tidyr)
    library(tibble)
    library(msm)
    library(abind)
    library(BiocParallel)
})
