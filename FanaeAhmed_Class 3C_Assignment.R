rm(list = ls())
gc()

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.rstudio.com/")
}

BiocManager::install(c("limma", "AnnotationDbi", "hgu133a.db"), update = FALSE, ask = FALSE)
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"), repos = "http://cran.rstudio.com/", quiet = TRUE)

suppressPackageStartupMessages({
  library(limma)           
  library(AnnotationDbi)   
  library(hgu133a.db)      
  library(dplyr)           
  library(tibble)         
  library(ggplot2)         
  library(pheatmap)       
})

cat("All required packages loaded successfully.\n")

getwd()

setwd("D:/R_Projects/AI_Omics_Internship_2025/Module_2/GSE41258") 

data <- readRDS("GSE41258_expression_filtered.RDS") 
data <- as.matrix(data) 
sample_info <- read.csv("Results/sample_groups_detailed.csv")

# Load Processed Data and Prepare Groups 
data <- readRDS("GSE41258_expression_filtered.RDS")
data <- as.matrix(data) 

sample_info <- read.csv("Results/sample_groups_detailed.csv")

groups <- factor(sample_info$Group, levels = c("Normal", "Cancer", "CellLine"))

# Filter both data and groups to keep only Normal and Cancer samples
samples_to_keep <- groups %in% c("Normal", "Cancer")
data_filtered <- data[, samples_to_keep]
groups_final <- droplevels(groups[samples_to_keep])

cat("\n--- Data Ready for Analysis ---\n")
cat("Expression Matrix Dimensions:", dim(data_filtered), "\n")
cat("Sample Groups:\n")
print(table(groups_final))


# Probe IDs to Gene Mapping and Duplicate Handling 

# Map Probe IDs to Gene Symbols
probe_ids <- rownames(data_filtered)

# Map probe IDs using the hgu133a.db package
gene_symbols <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"  
)

# Convert mapping to a data frame and remove probes without a symbol (NA)
gene_map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2) %>%
  dplyr::filter(!is.na(SYMBOL))

# Prepare data for averaging: only keep probes that successfully mapped
data_annotated <- data_filtered[gene_map_df$PROBEID, ]
symbols_to_use <- gene_map_df$SYMBOL


# Check and Handle Duplicates
# Count how many probes map to the same gene
duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  filter(probes_per_gene > 1)
total_probes_to_collapse <- sum(duplicate_summary$probes_per_gene)

cat("\n--- Probe Mapping Summary ---\n")
cat("Total probes involved in duplicate mapping:", total_probes_to_collapse, "\n")

# Collapse multiple probes per gene using the average expression
averaged_data <- limma::avereps(data_annotated, ID = symbols_to_use)
data_final <- as.matrix(averaged_data)

cat("Final number of unique genes for DEG analysis:", nrow(data_final), "\n")


# Differential Gene Expression Analysis (Limma) 
# Create the design matrix for linear modeling
design <- model.matrix(~0 + groups_final)
colnames(design) <- levels(groups_final) # Rename columns to 'Normal' and 'Cancer'

# Fit linear model
fit_1 <- lmFit(data_final, design)

# Define the contrast: Cancer vs Normal
contrast_matrix <- makeContrasts(cancer_vs_normal = Cancer - Normal,
                                 levels = design)

# Apply contrasts and compute moderated t-statistics (eBayes)
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)
fit_2 <- eBayes(fit_contrast)

# Extract all DEG results
deg_results <- topTable(fit_2,
                        coef = "cancer_vs_normal", # The comparison we defined
                        number = Inf,               # Get all genes
                        adjust.method = "BH")       # Benjamini-Hochberg correction


# Classify and Save DEGs
# Define significance cut-offs: adj. P-value < 0.05 AND |logFC| > 1
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "No")
))

# Subset genes based on regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")
deg_updown <- rbind(upregulated, downregulated) 

cat("\n--- Final DEG Count ---\n")
num_upregulated <- nrow(upregulated)
num_downregulated <- nrow(downregulated)
cat("Upregulated Genes:", num_upregulated, "\n")
cat("Downregulated Genes:", num_downregulated, "\n")


# Save results to CSV files
if (!dir.exists("Results")) dir.create("Results")
write.csv(deg_results, file = "Results/GSE41258_DEGs_Complete.csv", row.names = TRUE)
write.csv(upregulated, file = "Results/GSE41258_Upregulated_DEGs.csv", row.names = TRUE)
write.csv(downregulated, file = "Results/GSE41258_Downregulated_DEGs.csv", row.names = TRUE)
cat("\nDEG results saved in 'Results/' folder.\n")


# Data Visualization: Volcano Plot 
# Start plotting to a PNG file
png("Results/Volcano_Plot_GSE41258.png", width = 2000, height = 1500, res = 300)

volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  # Set colors for the groups
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: Cancer vs Normal (GSE41258)",
       x = expression(log[2]~"Fold Change"), 
       y = expression(-log[10]~"Adjusted P-value"),
       color = "Regulation") +
  # Add cut-off lines for logFC and p-value
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

print(volcano_plot)
dev.off() 
cat("Volcano Plot saved: Results/Volcano_Plot_GSE41258.png\n")


# Data Visualization: Heatmap of Top 25 DEGs
# Select top 25 genes based on the smallest adjusted p-values from the DEG list
top_genes_25 <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 25)

# Subset the final gene-level expression matrix
heatmap_data <- data_final[top_genes_25, ]

# Create annotation data frame for samples (columns)
sample_annotation <- data.frame(Group = groups_final)
rownames(sample_annotation) <- colnames(heatmap_data)

# Define colors for the sample groups in the annotation bar
ann_colors <- list(
  Group = c(Normal = "forestgreen", Cancer = "firebrick")
)

# Start plotting to a PNG file
png("Results/Heatmap_Top25_DEGs_GSE41258.png", width = 2000, height = 2500, res = 300)

pheatmap(
  heatmap_data,
  scale = "row",                     
  cluster_rows = TRUE,               
  cluster_cols = TRUE,               
  show_rownames = TRUE,              
  show_colnames = FALSE,             
  annotation_col = sample_annotation,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  main = "Heatmap of Top 25 Differentially Expressed Genes"
)

dev.off() # Close the PNG device
cat("Heatmap saved: Results/Heatmap_Top25_DEGs_GSE41258.png\n")


# Results Summary
cat("i. A total of", total_probes_to_collapse, "probes were involved in mapping to the same gene (multiple probes per gene symbol). These duplicates were handled by using the **limma::avereps()** function, which computes the average expression value for all probes associated with a single gene symbol.\n")
cat("ii. The contrast performed for differential gene expression analysis was: cancer_vs_normal.\n")
cat("iii. Based on the cut-offs (adj. P-value < 0.05 and |logFC| > 1), we found", num_upregulated, "genes to be upregulated and", num_downregulated, "genes to be downregulated in the cancer samples compared to normal samples.\n")