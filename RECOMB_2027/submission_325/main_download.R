suppressMessages (library (dplyr))
suppressMessages (library (Matrix))
suppressMessages (library (TCGAbiolinks))
suppressMessages (library (GenomicRanges))

# Rscript main_download.R BASE_DIRECTORY

args <- commandArgs (trailingOnly = TRUE)
setwd (args[[1]])
if (!dir.exists ("miRNA"))
    dir.create ("miRNA")
if (!dir.exists ("./miRNA/data/"))
    dir.create ("./miRNA/data/")
if (!dir.exists ("RNA"))
    dir.create ("RNA")
if (!dir.exists ("./RNA/data/"))
    dir.create ("./RNA/data/")


projects <- c(); barcodes <- data.frame ()
for (proj in sort (filter (getGDCprojects (), grepl ("^TCGA-", id))$id))
{
    print(proj)
    miRNA_data <- GDCquery (project = proj, data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification",
                            experimental.strategy = "miRNA-Seq")$results[[1]]
    tumor <- miRNA_data[miRNA_data$sample_type == "Primary Tumor",  c("cases", "cases.submitter_id", "sample.submitter_id")]
    normal <- miRNA_data[miRNA_data$sample_type == "Solid Tissue Normal",  c("cases", "cases.submitter_id", "sample.submitter_id")]
    miRNA_data <- merge (tumor, normal, by = "cases.submitter_id", all = FALSE, suffixes = c("_tumor", "_normal"))
    if (nrow (miRNA_data) == 0)
        next
    RNA_data <- GDCquery (project = proj, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts")$results[[1]]
    tumor <- RNA_data[RNA_data$sample_type == "Primary Tumor",  c("cases", "cases.submitter_id", "sample.submitter_id")]
    normal <- RNA_data[RNA_data$sample_type == "Solid Tissue Normal",  c("cases", "cases.submitter_id", "sample.submitter_id")]
    RNA_data <- merge (tumor, normal, by = "cases.submitter_id", all = FALSE, suffixes = c("_tumor", "_normal"))
    if (nrow (RNA_data) == 0)
        next
    
    metadata <- merge (miRNA_data, RNA_data, by = c("cases.submitter_id", "sample.submitter_id_tumor", "sample.submitter_id_normal"),
                       all = FALSE, suffixes = c("_miRNA", "_RNA"))
    metadata <- metadata[, c("cases.submitter_id", "cases_tumor_miRNA", "cases_normal_miRNA", "cases_tumor_RNA", "cases_normal_RNA")]
    miRNA_data <- metadata[, c("cases.submitter_id", "cases_tumor_miRNA", "cases_normal_miRNA")]
    RNA_data <- metadata[, c("cases.submitter_id", "cases_tumor_RNA", "cases_normal_RNA")]
    miRNA_data <- miRNA_data[!duplicated (miRNA_data), ]; print (paste ("miRNA:", nrow (miRNA_data), "sample pairs"))
    RNA_data <- RNA_data[!duplicated (RNA_data), ]; print (paste ("RNA:", nrow (RNA_data), "sample pairs"))
    metadata <- metadata[!duplicated (metadata), ]; print (paste ("merged:", nrow (metadata), "sample pairs"))
    if (nrow (metadata) < 30)
        next
    
    colnames (metadata) <- c("sample", "miRNA_tumor", "miRNA_normal", "RNA_tumor", "RNA_normal")
    metadata$context <- strsplit (proj, "-")[[1]][2]
    projects <- c(projects, proj); barcodes <- rbind (barcodes, metadata)
}
metadata <- data.frame ()
for (sample in unique (barcodes$sample))
{
    tmp <- barcodes[barcodes$sample == sample, ]
    tmp$sample <- paste (tmp$sample, 0:(nrow (tmp) - 1), sep = "_")
    metadata <- rbind (metadata, tmp)
}
miRNA_barcodes <- metadata[, c("sample", "miRNA_tumor", "miRNA_normal", "context")]
colnames (miRNA_barcodes) <- c("sample", "tumor", "normal", "context")
miRNA_barcodes <- miRNA_barcodes[!duplicated (miRNA_barcodes[, c("tumor", "normal", "context")]), ]
write.table (miRNA_barcodes, "./miRNA/data/metadata.tsv", row.names = FALSE, sep = "\t", quote = FALSE)
RNA_barcodes <- metadata[, c("sample", "RNA_tumor", "RNA_normal", "context")]
colnames (RNA_barcodes) <- c("sample", "tumor", "normal", "context")
RNA_barcodes <- RNA_barcodes[!duplicated (RNA_barcodes[, c("tumor", "normal", "context")])]
write.table (RNA_barcodes, "./RNA/data/metadata.tsv", row.names = FALSE, sep = "\t", quote = FALSE)
rm (tmp, tumor, normal, sample, metadata, RNA_data, miRNA_data)

miRNA_df <- data.frame (); RNA_df <- data.frame ()
for (proj in projects)
{
    print (proj)
    barcode <- miRNA_barcodes[miRNA_barcodes$context == strsplit (proj, "-")[[1]][2], ]
    barcode <- unique (c(barcode$tumor, barcode$normal))
    miRNA_data <- GDCquery (project = proj, data.category = "Transcriptome Profiling", data.type = "miRNA Expression Quantification",
                            experimental.strategy = "miRNA-Seq", barcode = barcode)
    GDCdownload (miRNA_data, files.per.chunk = 50)
    miRNA_data <- GDCprepare (miRNA_data); miRNA <- miRNA_data[, grep ("read_count_", colnames (miRNA_data))]
    rownames (miRNA) <- miRNA_data$miRNA_ID; colnames (miRNA) <- gsub ("read_count_", "", colnames (miRNA))
    miRNA_df <- cbind (miRNA_df, miRNA); rm (miRNA, miRNA_data)
    
    barcode <- RNA_barcodes[RNA_barcodes$context == strsplit (proj, "-")[[1]][2], ]
    barcode <- unique (c(barcode$tumor, barcode$normal))
    RNA_data <- GDCquery (project = proj, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts", barcode = barcode)
    GDCdownload (RNA_data, files.per.chunk = 50)
    RNA_data <- GDCprepare (RNA_data)
    RNA <- data.frame (RNA_data@assays@data@listData$unstranded); rowRanges <- data.frame (RNA_data@rowRanges)
    RNA <- cbind (data.frame (gene = rowRanges$gene_name), RNA)
    RNA <- aggregate (. ~ gene, RNA, sum); rownames (RNA) <- RNA$gene; colnames (RNA) <- barcode; RNA <- RNA[, -1]
    RNA <- RNA[unique (rowRanges[(!rowRanges$seqnames %in% c("chrX", "chrY")) & (rowRanges$gene_type == "protein_coding"), "gene_name"]), ]
    RNA_df <- cbind (RNA_df, RNA); rm (RNA, RNA_data, rowRanges)
}
miRNA_df <- cbind (data.frame (feature = rownames (miRNA_df)), miRNA_df); colnames (miRNA_df)[1] <- ""
write.table (miRNA_df, "./miRNA/data/raw_counts.tsv", row.names = FALSE, sep = "\t", quote = FALSE)
RNA_df <- cbind (data.frame (gene = rownames (RNA_df)), RNA_df); colnames (RNA_df)[1] <- ""
write.table (RNA_df, "./RNA/data/raw_counts.tsv", row.names = FALSE, sep = "\t", quote = FALSE)


