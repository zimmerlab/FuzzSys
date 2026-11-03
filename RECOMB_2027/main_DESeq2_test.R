suppressMessages (library (DESeq2))

# Rscript main_DESeq2_test.R COUNT_MATRIX METADATA OUTPUT_DIRECTORY

args <- commandArgs (trailingOnly = TRUE)
mtxPath <- args[[1]]; metadataPath <- args[[2]]; outputDir <- args[[3]]

mtx <- read.table (mtxPath, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata <- read.table (metadataPath, header = TRUE, row.names = NULL, sep = "\t", check.names = FALSE)
allClusters <- unique (metadata$context); allClusters <- allClusters[order (allClusters)]

colData <- data.frame (sample = c(metadata$tumor, metadata$normal),
                       cluster = c(metadata$context, metadata$context),
                       type = c(rep ("tumor", nrow (metadata)), rep ("normal", nrow (metadata))),
                       stringsAsFactors = TRUE)
colData <- colData[!duplicated (colData), ]; rownames (colData) <- colData$sample; colData <- colData[colnames (mtx), ]
dds <- DESeqDataSetFromMatrix (round (mtx), colData = colData, design = ~ type); dds <- estimateSizeFactors (dds)

effect_size = data.frame (row.names = rownames (mtx)); significance = data.frame (row.names = rownames (mtx))
for (cluster in allClusters)
{
    print (cluster)
    cluster_dds <- DESeq (dds[, dds$cluster == cluster])
    res <- results (cluster_dds, contrast = c("type", "tumor", "normal"))
    effect_size[, cluster] <- res$log2FoldChange; significance[, cluster] <- res$padj
}
significance[significance < 1e-10] <- 1e-10; significance <- -log10 (significance)

if (!dir.exists (outputDir))
{
    dir.create (outputDir, recursive = TRUE)
}
write.table (counts (dds, normalized = TRUE), file.path (outputDir, "normalized_counts.tsv"), sep = "\t", quote = FALSE)
write.table (effect_size, file.path (outputDir, "DESeq2_log2FC.tsv"), sep = "\t", quote = FALSE)
write.table (significance, file.path (outputDir, "DESeq2_padj.tsv"), sep = "\t", quote = FALSE)


