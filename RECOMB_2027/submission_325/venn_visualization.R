suppressMessages (library (ggplot2))
suppressMessages (library (ggVennDiagram))

# Rscript venn_visualization.R rawLog2FC_RESULT normLog2FC_RESULT DESeq2_RESULT CLUSTER REGULATION OUTPUT_PREFIX

args <- commandArgs (trailingOnly = TRUE)
pathRFC <- args[[1]]; pathNFC <- args[[2]]; pathDESeq2 <- args[[3]]
cluster <- args[[4]]; regulation <- args[[5]]; output <- args[[6]]

RFC_results <- read.table (pathRFC, header = TRUE, row.names = NULL, sep = "\t", check.names = FALSE)
NFC_results <- read.table (pathNFC, header = TRUE, row.names = NULL, sep = "\t", check.names = FALSE)
DESeq2_results <- read.table (pathDESeq2, header = TRUE, row.names = NULL, sep = "\t", check.names = FALSE)
allClusters <- sort (unique (c(RFC_results$cluster, NFC_results$cluster, DESeq2_results$cluster)))

color <- ifelse (regulation == "DOWN", "steelblue", "firebrick")
set <- list (raw = sort (unique (RFC_results[RFC_results$cluster == cluster & RFC_results$regulation == regulation, "feature"])),
             norm = sort (unique (NFC_results[NFC_results$cluster == cluster & NFC_results$regulation == regulation, "feature"])),
             DESeq2 = sort (unique (DESeq2_results[DESeq2_results$cluster == cluster & DESeq2_results$regulation == regulation, "feature"])))
p <- ggVennDiagram (set, category.names = c("raw log2FC", "norm. log2FC", "DESeq2"), set_size = 6, label_size = 7.5) +
    scale_fill_gradient (low = "grey95", high = color) +
    scale_x_continuous (expand = expansion (mult = 0.25))
ggsave (paste0 (output, "_venn.png"), plot = p)
p <- ggVennDiagram (set, category.names = c("raw log2FC", "norm. log2FC", "DESeq2"),
                    force_upset = TRUE, order.set.by = "name")
ggsave (paste0 (output, "_upset.png"), plot = p)


