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

#-------------------------------------------------------------------------
# 1. REUSABLE DEG FUNCTION
#
# Wraps: GEO download -> log2 check -> group factor -> design matrix ->
#        limma fit -> contrast -> topTable -> probe annotation
#-------------------------------------------------------------------------

run_deg_analysis <- function(gse_id,
                             group_meta_column,
                             annotation_db,
                             contrast_groups = NULL,
                             keytype = "PROBEID",
                             p_cutoff = 0.05,
                             dataset_label = gse_id) {
  
  message("[", dataset_label, "] Downloading ", gse_id, " from GEO...")
  gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
  ex <- exprs(gset)
  
  qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  needs_log <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (needs_log) {
    ex[ex <= 0] <- NA
    ex <- log2(ex)
  }
  
  group_info <- pData(gset)[[group_meta_column]]
  groups <- make.names(group_info)
  gset$group <- factor(groups)
  group_names <- levels(gset$group)
  
  if (is.null(contrast_groups)) {
    if (length(group_names) < 2) {
      stop("[", dataset_label, "] Fewer than 2 groups found in '",
           group_meta_column, "'. Check the metadata column name.")
    }
    contrast_groups <- group_names[1:2]
    message("[", dataset_label, "] No contrast_groups given. Using '",
            contrast_groups[1], "' vs '", contrast_groups[2],
            "'. Confirm this is the direction you intend (order matters for logFC sign).")
  }
  
  design <- model.matrix(~0 + gset$group)
  colnames(design) <- group_names
  contrast_formula <- paste(contrast_groups[1], "-", contrast_groups[2])
  
  fit <- lmFit(ex, design)
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  deg_table <- topTable(
    fit2, adjust = "fdr", sort.by = "B", number = Inf, p.value = p_cutoff
  )
  
  probe_ids <- rownames(deg_table)
  gene_annotation <- AnnotationDbi::select(
    annotation_db, keys = probe_ids,
    columns = c("SYMBOL", "GENENAME"), keytype = keytype
  )
  
  deg_table$PROBEID <- rownames(deg_table)
  deg_table <- merge(deg_table, gene_annotation,
                     by.x = "PROBEID", by.y = keytype, all.x = TRUE)
  deg_table$PROBEID <- rownames(deg_table)
  deg_table <- merge(deg_table, gene_annotation,
                     by.x = "PROBEID", by.y = keytype, all.x = TRUE)
  deg_table <- deg_table[order(deg_table$adj.P.Val), ]
  
  message("[", dataset_label, "] Done. ", nrow(deg_table),
          " significant probes at adj.P.Val < ", p_cutoff, ".")
  
  list(
    dataset_label = dataset_label,
    gse_id = gse_id,
    gset = gset,
    ex = ex,
    design = design,
    contrast_formula = contrast_formula,
    fit2 = fit2,
    deg_table = deg_table
  )
}