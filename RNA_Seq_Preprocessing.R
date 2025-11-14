raw_path <- "1_RawData/"
output_path <- "3_Outputs/"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages(c("readr", "dplyr", "tibble", "ggplot2", "ggrepel"))

library(readr)
library(dplyr)
library(tibble)

# IMPORT METADATA
colData_raw <- readr::read_csv(paste0(raw_path, "Master_Sample_List.csv"), show_col_types = FALSE)

colData_clean <- colData_raw |>
  dplyr::filter(Dataset == "GSE277906") |>
  dplyr::select(Sample_ID, Condition = Phenotype, Sample_title) |>
  tibble::column_to_rownames(var = "Sample_title")

colData_clean$Condition <- factor(colData_clean$Condition)
colData_clean$Condition <- relevel(colData_clean$Condition, ref = "Normal")
cat(paste(nrow(colData_clean), "samples are found."))

# IMPORT COUNT DATA
countData <- read.delim(
  file = paste0(raw_path, "GSE277906_counts_anno.txt.gz"),
  sep = "\t",
  header = TRUE,
  row.names = 1
)
countData <- as.matrix(countData)
cat(paste(nrow(countData), "genes are found."))

# ALIGNMENT
library(DESeq2)

valid_samples <- intersect(colnames(countData), rownames(colData_clean))

countData <- countData[, valid_samples]
colData_clean <- colData_clean[valid_samples, ]
countData <- countData[, rownames(colData_clean)]

alignment_check <- all(colnames(countData) == rownames(colData_clean))
cat(paste("Sample Alignment Check:", alignment_check))

countData[is.na(countData)] <- 0
storage.mode(countData) <- "integer"

# DESeq2 OBJECT
if (alignment_check) {
  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData_clean,
    design = ~ Condition
  )
  print(dds)
  saveRDS(dds, paste0(output_path, "dds_raw.rds"))
}

# FILTER LOW-EXPRESSION GENES
library(ggplot2)
library(ggrepel)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(paste(nrow(dds), "genes are filtered for analysis."))

# VARIANCE-STABILIZING TRANSFORMATION (VST)
vsd <- vst(dds, blind = FALSE)

# PCA QC Plot 
pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Condition, shape = Condition)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = name), size = 3, max.overlaps = 50) +
  xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
  labs(
    title = "QC: Principal Component Analysis (PCA)",
    subtitle = "PCOS vs. Normal Cumulus Cells"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
ggsave(paste0(output_path, "PCA_plot.png"), plot = pca_plot, width = 8, height = 7)
cat(paste0("PCA plot saved to: ", output_path))