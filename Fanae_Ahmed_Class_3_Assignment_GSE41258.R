# Clear workspace
rm(list = ls())

setwd("D:/R_Projects/AI_Omics_Internship_2025/Module_3/GSE41258")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.rstudio.com/")
}

# Install Bioconductor packages
BiocManager::install(c("GEOquery", "affy", "arrayQualityMetrics", 
                       "hgu133acdf", "hgu133a.db"), 
                     update = FALSE, ask = FALSE)

# Install CRAN packages
install.packages(c("dplyr", "ggplot2", "RColorBrewer"), 
                 repos = "http://cran.rstudio.com/")

# Load libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(affy)
  library(arrayQualityMetrics)
  library(hgu133acdf)
  library(hgu133a.db)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
})

cat("✓ All packages loaded successfully\n\n")

cat("=== Downloading Metadata ===\n")
gse_data <- getGEO("GSE41258", GSEMatrix = TRUE)
phenotype_data <- pData(gse_data[[1]])

# Extract and clean tissue types
tissue_clean <- gsub("tissue: ", "", phenotype_data$characteristics_ch1)

# Define sample groups
normal_tissues <- c("Normal Colon", "Normal Liver", "Normal Lung")
cancer_tissues <- c("Primary Tumor", "Liver Metastasis", "Lung Metastasis", 
                    "Polyp", "Polyp, high grade", "Microadenoma")

# Create grouping variable
sample_groups <- case_when(
  tissue_clean %in% normal_tissues ~ "Normal",
  tissue_clean %in% cancer_tissues ~ "Cancer",
  grepl("cell line:", tissue_clean) ~ "CellLine"
)
sample_groups <- factor(sample_groups, levels = c("Normal", "Cancer", "CellLine"))

# Filter to exclude cell lines
samples_to_keep <- sample_groups %in% c("Normal", "Cancer")
filtered_groups <- droplevels(sample_groups[samples_to_keep])

cat("\nSample Distribution:\n")
print(table(sample_groups))
cat("\nAnalysis Groups (excluding cell lines):\n")
print(table(filtered_groups))

cat("\n=== Loading CEL Files ===\n")
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")

# Get probe count before filtering
n_probes_before <- length(featureNames(raw_data))
cat("Number of samples:", ncol(raw_data), "\n")
cat("Number of probes (before filtering):", n_probes_before, "\n")

# Match sample order with phenotype data
sample_names <- sampleNames(raw_data)
cel_gsm <- sapply(strsplit(sample_names, "_"), function(x) x[1])
match_indices <- match(cel_gsm, phenotype_data$geo_accession)
groups_ordered <- sample_groups[match_indices]

cat("\n=== Performing RMA Normalization ===\n")
normalized_data <- rma(raw_data)
expr_matrix <- exprs(normalized_data)

cat("Normalization complete!\n")
cat("Expression matrix dimensions:", nrow(expr_matrix), "×", ncol(expr_matrix), "\n")

cat("\n=== Running Quality Control ===\n")
cat("This may take 10-20 minutes...\n\n")

arrayQualityMetrics(
  expressionset = normalized_data,
  outdir = "Results/QC_Normalized_Data",
  force = TRUE,
  do.logtransform = FALSE
)

cat("\n QC Report generated: Results/QC_Normalized_Data/index.html\n")

cat("\n=== Creating Diagnostic Plots ===\n")

# Create Plots directory
if (!dir.exists("Plots")) dir.create("Plots")

# Custom color palette
sample_colors <- colorRampPalette(brewer.pal(12, "Set3"))(ncol(expr_matrix))

# 1. Boxplot
png("Plots/Boxplot_Normalized.png", width = 1800, height = 900, res = 150)
par(mar = c(8, 5, 3, 2))
boxplot(expr_matrix,
        main = "RMA Normalized Expression Distribution (GSE41258)",
        col = sample_colors,
        las = 2,
        cex.axis = 0.5,
        xlab = "",
        ylab = "Log2 Expression",
        outline = FALSE)
dev.off()

# 2. Density Plot
png("Plots/Density_Normalized.png", width = 1600, height = 900, res = 150)
plot(density(expr_matrix[,1]), 
     main = "Density Plot - Normalized Expression",
     xlab = "Log2 Expression",
     ylab = "Density",
     col = sample_colors[1],
     lwd = 2,
     ylim = c(0, max(apply(expr_matrix, 2, function(x) max(density(x)$y)))))
for(i in 2:ncol(expr_matrix)) {
  lines(density(expr_matrix[,i]), col = sample_colors[i], lwd = 2)
}
dev.off()

# 3. PCA Plot (Beautiful version)
pca_result <- prcomp(t(expr_matrix), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Group = groups_ordered
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = c("Normal" = "#2E86AB", 
                                "Cancer" = "#E63946", 
                                "CellLine" = "#F77F00")) +
  labs(title = "PCA Plot - GSE41258 (Normalized Data)",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid.major = element_line(color = "gray90"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))

ggsave("Plots/PCA_Normalized.png", pca_plot, 
       width = 10, height = 8, dpi = 300)

cat("Plots saved in Plots/ directory\n")

cat("\n=== Filtering Low-Intensity Probes ===\n")

# Calculate median intensity
row_medians <- rowMedians(expr_matrix)

# Determine threshold
threshold <- 5.0

# Apply filter
expr_filtered <- expr_matrix[row_medians > threshold, ]

n_probes_after <- nrow(expr_filtered)

cat("Probes before filtering:", n_probes_before, "\n")
cat("Probes after filtering:", n_probes_after, "\n")
cat("Probes removed:", n_probes_before - n_probes_after, "\n")

cat("\n=== Saving Results ===\n")

# Save expression data
write.csv(expr_filtered, "Results/GSE41258_filtered_expression.csv")
saveRDS(expr_filtered, "Results/GSE41258_filtered_expression.rds")

# Save sample groups
sample_info <- data.frame(
  Sample_ID = colnames(expr_filtered),
  GSM_ID = cel_gsm,
  Group = groups_ordered,
  Tissue_Type = tissue_clean[match_indices]
)
write.csv(sample_info, "Results/GSE41258_sample_info.csv", row.names = FALSE)

cat("All results saved\n")

# Open QC report
browseURL("Results/QC_Normalized_Data/index.html")