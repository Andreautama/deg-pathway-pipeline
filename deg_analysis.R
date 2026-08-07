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
                             dataset_label = gse_id,
                             local_file = NULL) {
  
  if (!is.null(local_file)) {
    message("[", dataset_label, "] Loading ", gse_id,
            " from local file (GEO server download was unreliable): ", local_file)
    gset <- getGEO(filename = local_file, getGPL = FALSE)
  } else {
    message("[", dataset_label, "] Downloading ", gse_id, " from GEO...")
    gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
  }
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

#-------------------------------------------------------------------------
# Volcano plot: logFC vs. significance, highlighting up/down regulated genes
#-------------------------------------------------------------------------

volcano_data <- lung_result$deg_table
volcano_data$direction <- "Not significant"
volcano_data$direction[volcano_data$adj.P.Val < 0.01 & volcano_data$logFC > 1] <- "Up in tumor"
volcano_data$direction[volcano_data$adj.P.Val < 0.01 & volcano_data$logFC < -1] <- "Down in tumor"

top_up <- volcano_data[volcano_data$direction == "Up in tumor", ]
top_up <- top_up[order(top_up$adj.P.Val), ][1:8, ]

top_down <- volcano_data[volcano_data$direction == "Down in tumor", ]
top_down <- top_down[order(top_down$adj.P.Val), ][1:8, ]

top_labels <- rbind(top_up, top_down)

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up in tumor" = "firebrick",
                                "Down in tumor" = "steelblue",
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_labels, aes(label = SYMBOL),
                  size = 3, color = "black", max.overlaps = 20) +
  labs(title = "Volcano Plot: Lung Adenocarcinoma vs Normal (GSE10072)",
       x = "log2 Fold Change", y = "-log10(adj.P.Val)", color = NULL) +
  theme_minimal()

#-------------------------------------------------------------------------
# Heatmap: top 30 most significant genes across all samples
#-------------------------------------------------------------------------

top_genes <- lung_result$deg_table[1:30, ]
heatmap_matrix <- lung_result$ex[top_genes$PROBEID, ]
rownames(heatmap_matrix) <- ifelse(is.na(top_genes$SYMBOL) | top_genes$SYMBOL == "",
                                   top_genes$PROBEID, top_genes$SYMBOL)

sample_groups <- data.frame(
  Group = lung_result$gset$group,
  row.names = colnames(heatmap_matrix)
)

pheatmap(heatmap_matrix,
         scale = "row",
         annotation_col = sample_groups,
         show_colnames = FALSE,
         main = "Top 30 DE Genes: Lung Adenocarcinoma vs Normal (GSE10072)")

#-------------------------------------------------------------------------
# 2. PATHWAY ENRICHMENT FUNCTION
#
# Takes the DEG table from run_deg_analysis() and asks: which biological
# processes (GO) and pathways (KEGG) are over-represented among the
# significant genes?
#-------------------------------------------------------------------------

run_enrichment <- function(deg_table, pvalue_cutoff = 0.05) {
  
  gene_symbols <- unique(na.omit(deg_table$SYMBOL))
  gene_symbols <- gene_symbols[gene_symbols != ""]
  
  if (length(gene_symbols) < 5) {
    warning("Fewer than 5 annotated gene symbols in DEG list; ",
            "enrichment results may be unreliable or empty.")
  }
  
  id_map <- bitr(
    gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  go_result <- enrichGO(
    gene = id_map$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = pvalue_cutoff,
    readable = TRUE
  )
  
  kegg_result <- enrichKEGG(
    gene = id_map$ENTREZID,
    organism = "hsa",
    pvalueCutoff = pvalue_cutoff
  )
  
  list(go = go_result, kegg = kegg_result, id_map = id_map)
}

#-------------------------------------------------------------------------
#LUNG ADENOCARCINOMA DATASET (GSE10072)
#
# source_name_ch1 metadata contains two labels: "Adenocarcinoma of the Lung"
# and "Normal Lung Tissue" (confirmed by inspecting the GEO series matrix
# directly). After make.names() conversion these become
# "Adenocarcinoma.of.the.Lung" and "Normal.Lung.Tissue". Contrast is set
# explicitly (tumor vs normal) so positive logFC means higher expression
# in tumor tissue, rather than relying on factor level order.

# NOTE: GEO's automatic download for this series was unreliable (repeated
# connection failures), so this run uses a manually downloaded local copy
# of the series matrix via the local_file argument, with getGPL = FALSE
# (platform annotation is not needed since hgu133a.db supplies gene
# symbols downstream).
#-------------------------------------------------------------------------

lung_result <- run_deg_analysis(
  gse_id            = "GSE10072",
  group_meta_column = "source_name_ch1",
  annotation_db     = hgu133a.db,
  contrast_groups   = c("Adenocarcinoma.of.the.Lung", "Normal.Lung.Tissue"),
  p_cutoff          = 0.01,
  dataset_label     = "Lung adenocarcinoma (GSE10072)",
  local_file        = "DATASET_GSE10072.txt"
)

head(lung_result$deg_table[, c("PROBEID", "SYMBOL", "logFC", "adj.P.Val")])

lung_enrichment <- run_enrichment(lung_result$deg_table, pvalue_cutoff = 0.05)

if (!is.null(lung_enrichment$go) && nrow(lung_enrichment$go@result) > 0) {
  dotplot(lung_enrichment$go, showCategory = 15) +
    ggtitle("GO Biological Process Enrichment: Lung Adenocarcinoma (GSE10072)")
}

if (!is.null(lung_enrichment$kegg) && nrow(lung_enrichment$kegg@result) > 0) {
  dotplot(lung_enrichment$kegg, showCategory = 15) +
    ggtitle("KEGG Pathway Enrichment: Lung Adenocarcinoma (GSE10072)")
}