setwd("D:/R_Projects/AI_Omics_Internship_2025/Module_2")
cat("Working directory set to:\n", getwd(), "\n\n")

results_dir <- "Results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
  cat("Created folder:", results_dir, "\n\n")
} else {
  cat("Results folder already exists:", results_dir, "\n\n")
}

# detect actual column name 
detect_col <- function(df, candidates) {
  # returns the actual column name present in df (matching candidate case-insensitively)
  for (cand in candidates) {
    matches <- names(df)[tolower(names(df)) == tolower(cand)]
    if (length(matches) == 1) return(matches)
  }
  return(NULL)
}

# Input: single values logFC and padj
classify_gene <- function(logFC, padj) {
  logFC <- as.numeric(logFC)
  padj  <- as.numeric(padj)
  if (!is.na(padj) && padj < 0.05 && !is.na(logFC) && logFC > 1) return("Upregulated")
  if (!is.na(padj) && padj < 0.05 && !is.na(logFC) && logFC < -1) return("Downregulated")
  return("Not_Significant")
}

expected_files <- c("DEGs_data_1.csv", "DEGs_data_2.csv")

files_found <- expected_files[file.exists(expected_files)]

if (length(files_found) == 0) {
  cat("Expected filenames not found. Searching for files matching 'DEGs' pattern...\n")
  files_found <- list.files(pattern = "^DEGs.*\\.csv$", ignore.case = TRUE)
}

if (length(files_found) == 0) {
  stop("No DEGs CSV files found in the working directory. Please place DEGs_data_1.csv and DEGs_data_2.csv here, or change setwd().")
}

cat("Files to be processed:\n")
print(files_found)
cat("\n")

#Column name candidates
padj_candidates <- c("padj", "adj.P.Val", "adj_p", "FDR", "qvalue", "p.adjust")
logfc_candidates <- c("logFC", "log2FoldChange", "log2.FC", "log2_fold_change", "log2FoldChange:logFC")

#Process each file in a for-loop
total_summary <- list(Upregulated = 0, Downregulated = 0, Not_Significant = 0)

for (file_name in files_found) {
  cat("--------------------------------------------------\n")
  cat("Processing file:", file_name, "\n")
  df <- tryCatch({
    read.csv(file_name, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    cat("Error reading file:", file_name, " - ", e$message, "\n")
    next
  })
  
  # Detect columns
  padj_col <- detect_col(df, padj_candidates)
  logfc_col <- detect_col(df, logfc_candidates)
  geneid_col <- detect_col(df, c("Gene_Id", "GeneID", "gene", "Gene", "Gene_Id"))
  
  if (is.null(geneid_col)) {
    # If there is no gene id column, add row numbers as gene id (safe fallback)
    df$Gene_Id <- paste0("gene_", seq_len(nrow(df)))
    geneid_col <- "Gene_Id"
    cat("No Gene_Id column detected; created 'Gene_Id' with row identifiers.\n")
  }
  
  if (is.null(logfc_col)) {
    cat("ERROR: No logFC column found in", file_name, "\n")
    cat("Checked candidate names:", paste(logfc_candidates, collapse = ", "), "\n")
    cat("Skipping this file.\n\n")
    next
  }
  
  if (is.null(padj_col)) {
    cat("WARNING: No padj/adjusted p-value column found in", file_name, "\n")
    cat("Checked candidate names:", paste(padj_candidates, collapse = ", "), "\n")
    cat("If only raw p-values are present, please rename that column to 'padj' or one of the candidates.\n")
    cat("Skipping this file.\n\n")
    next
  }
  
  # Convert padj and logFC to numeric 
  df[[padj_col]] <- as.numeric(as.character(df[[padj_col]]))
  df[[logfc_col]] <- as.numeric(as.character(df[[logfc_col]]))
  
  # Replace missing padj values with 1 
  na_padj_count <- sum(is.na(df[[padj_col]]))
  if (na_padj_count > 0) {
    cat("Found", na_padj_count, "NA(s) in", padj_col, "- replacing with 1.\n")
    df[[padj_col]][is.na(df[[padj_col]])] <- 1
  }
  
  # Apply classify_gene across rows 
  df$status <- mapply(classify_gene, df[[logfc_col]], df[[padj_col]])
  
  # Add a boolean 'significant' column 
  df$significant <- df$status != "Not_Significant"
  
  # Print summary counts using table()
  cat("\nSummary for file:", file_name, "\n")
  print(table(df$status))
  up_count <- sum(df$status == "Upregulated", na.rm = TRUE)
  down_count <- sum(df$status == "Downregulated", na.rm = TRUE)
  nonsig_count <- sum(df$status == "Not_Significant", na.rm = TRUE)
  total_sig <- up_count + down_count
  cat("\nCounts --> Upregulated:", up_count,
      "| Downregulated:", down_count,
      "| Not_Significant:", nonsig_count,
      "| Total significant:", total_sig, "\n\n")
  
 
  total_summary$Upregulated <- total_summary$Upregulated + up_count
  total_summary$Downregulated <- total_summary$Downregulated + down_count
  total_summary$Not_Significant <- total_summary$Not_Significant + nonsig_count
  
  # Save processed file to Results folder
  out_name <- paste0(tools::file_path_sans_ext(basename(file_name)), "_processed.csv")
  out_path <- file.path(results_dir, out_name)
  write.csv(df, file = out_path, row.names = FALSE)
  cat("Saved processed file to:", out_path, "\n\n")
}



