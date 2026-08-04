#-------------------------------------------------------------------------
# 0. PACKAGES
#-------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Increase timeout for large Bioconductor annotation packages
# (org.Hs.eg.db, hgu133a.db, illuminaHumanv4.db can exceed the 60s default)
options(timeout = 600)
BiocManager::install(
  c(
    "GEOquery", "limma", "illuminaHumanv4.db", "hgu133a.db",
    "clusterProfiler", "org.Hs.eg.db", "enrichplot"
  ),
  ask = FALSE, update = FALSE
)

install.packages(c("pheatmap", "ggplot2", "dplyr"))

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(illuminaHumanv4.db)
library(hgu133a.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(pheatmap)