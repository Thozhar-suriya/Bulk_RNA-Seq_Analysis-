# ------ Global Knitr Options ------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  error = FALSE
)

# ------ R Package Installation ------
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "Rsubread", "DESeq2", "biomaRt",  "pheatmap", "enrichplot", "dplyr",
    "ggplot2", "ggrepel", "forcats"
))

# ------ Software Versions ------
library(Rsubread)
packageVersion("Rsubread")
library(DESeq2)
packageVersion("DESeq2")

# ------ 4. Build Rsubread Index ------
setwd("~/Suriya/Align")
library(Rsubread)

buildindex(
  basename = "Ref_Genome_GRCh38",
  reference = "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
  memory = 16000,
  indexSplit = TRUE
)

# ------ 5. Alignment in Rsubread ------
path <- "/home/vamouda/Suriya/Align/Seq"

reads1 <- list.files(path = path, pattern = "*_R1_val_1.fq.gz$", full.names = TRUE)
reads2 <- list.files(path = path, pattern = "*_R2_val_2.fq.gz$", full.names = TRUE)

if (length(reads1) != length(reads2)) stop("Mismatch between R1 and R2 files!")

sample_names <- gsub("_R1_val_1.fq.gz$", "", basename(reads1))

align(
  index = "Ref_Genome_GRCh38",
  readfile1 = reads1,
  readfile2 = reads2,
  input_format = "gzFASTQ",
  output_file = paste0(sample_names, ".bam"),
  output_format = "BAM",
  nthreads = 32,
  unique = TRUE
)

# ------ 6. featureCounts ------
bam_path <- "/home/vamouda/Suriya/Align/Bam"
bam_files <- list.files(path = bam_path, pattern = "\\.bam$", full.names = TRUE)
sample_names <- gsub("\\.bam$", "", basename(bam_files))

feature_Count <- featureCounts(
  files = bam_files,
  annot.ext = "/home/vamouda/Suriya/Align/GTF/Homo_sapiens.GRCh38.115.gtf.gz",
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  isPairedEnd = TRUE,
  requireBothEndsMapped = TRUE,
  minMQS = 10,
  strandSpecific = 0,
  countMultiMappingReads = FALSE,
  allowMultiOverlap = TRUE,
  largestOverlap = TRUE,
  nthreads = 32
)

write.csv(feature_Count$counts, "gene_expression_counts.csv")

# ------ 7. DESeq2 Analysis ------
library(DESeq2)
library(dplyr)

data <- read.csv("Count_data.csv", header = TRUE, row.names = 1)
meta <- read.csv("Meta_data.csv", row.names = 1)

head(data, 6)
head(meta, 6)

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = meta,
  design = ~ Tissue
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$Tissue <- factor(dds$Tissue, levels = c("Tumor", "Normal"))
dds$Tissue <- relevel(dds$Tissue, ref = "Normal")
dds$Tissue <- droplevels(dds$Tissue)

dds <- DESeq(dds)

res <- results(dds)
res05 <- results(dds, alpha = 0.05)

# Export DE results
de_genes <- subset(res, res$pvalue <= 0.05 &
                     (res$log2FoldChange > 1 | res$log2FoldChange < -1))
write.csv(de_genes, "Data_UpandDown.csv")

up_genes <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange > 1)
write.csv(up_genes, "Data_Upregulated.csv")

down_genes <- subset(res, res$pvalue <= 0.05 & res$log2FoldChange < -1)
write.csv(down_genes, "Data_Downregulated.csv")

# ------ Annotation Using biomaRt ------
library(biomaRt)
df <- read.csv("Data_UpandDown.csv")
genes <- df$ensembl_gene_id

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = genes,
  mart = mart
)

annotated_results <- merge(df, G_list, by = "ensembl_gene_id")
write.csv(annotated_results, "Data_UpandDown_withsymbol.csv")

# ------ 8. Normalization ------
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

write.csv(assay(vsd), "Data_vst_norm.csv")
write.csv(assay(rld), "Data_rlog_norm.csv")

# ------ Visualization Data Prep ------
DGE.results.sorted <- read.csv("Data_UpandDown.csv", row.names = 1)
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))

normalizeddata <- read.csv("Data_vst_norm.csv", row.names = 1)
Visualization <- normalizeddata[DGEgenes, ]
write.csv(Visualization, "Data_UpandDown_Visualization.csv")
