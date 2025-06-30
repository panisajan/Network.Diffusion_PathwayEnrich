
load("g.Rdata")

library("org.Hs.eg.db")
library(dplyr)
library(igraph)
library(KEGGREST)
library(diffuStats)
library(openxlsx)

go_bp <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "GO"), columns = c("GO", "SYMBOL", "ONTOLOGY"), keytype = "GO")

ErbB = c("GO:0038127", "GO:0007173", "GO:0038134", "GO:0038128", "GO:0038129", "GO:0038130",
         "GO:1901185", "GO:1901186", "GO:1901184")
go_bp_ErbB <- go_bp %>%
  filter(GO %in% ErbB & ONTOLOGY == "BP") %>%  # Filter for Biological Process (BP)
  dplyr::select(GO, SYMBOL)
gene_ErbB_go = unique(go_bp_ErbB$SYMBOL[go_bp_ErbB$SYMBOL %in% V(g)$name]) ## 56 genes

GSK3 = c("GO:0016055", "GO:0008286")
go_bp_GSK3 <- go_bp %>%
  filter(GO %in% GSK3 & ONTOLOGY == "BP") %>%  # Filter for Biological Process (BP)
  dplyr::select(GO, SYMBOL)
gene_GSK3_go = unique(go_bp_GSK3$SYMBOL[go_bp_GSK3$SYMBOL %in% V(g)$name]) ## 194 genes

p53 = c("GO:0072332", "GO:0042771", "GO:1902253")#,"GO:0006915")
go_bp_p53 <- go_bp %>%
  filter(GO %in% p53 & ONTOLOGY == "BP") %>%  # Filter for Biological Process (BP)
  dplyr::select(GO, SYMBOL)
gene_p53_go = unique(go_bp_p53$SYMBOL[go_bp_p53$SYMBOL %in% V(g)$name]) ## 46 genes


kegg_data <- read.csv("kegg_full.csv")

ErbB_kegg = c("hsa04012","hsa01521")
GSK3_kegg = c('hsa04310','hsa04910')
p53_kegg = c('hsa04115')#,'hsa04210')

kegg_ErbB_dat = subset(kegg_data, kegg_data$pathway_id %in% ErbB_kegg)
kegg_GSK3_dat = subset(kegg_data, kegg_data$pathway_id %in% GSK3_kegg)
kegg_p53_dat = subset(kegg_data, kegg_data$pathway_id %in% p53_kegg)

kegg_ErbB_dat = na.omit(kegg_ErbB_dat)
gene_ErbB_kegg = unique(strsplit(paste(kegg_ErbB_dat$gene_symbols, collapse = ","), ", ")[[1]])
gene_GSK3_kegg = unique(strsplit(paste(kegg_GSK3_dat$gene_symbols, collapse = ","), ", ")[[1]])
gene_p53_kegg = unique(strsplit(paste(kegg_p53_dat$gene_symbols, collapse = ","), ", ")[[1]])

gene_ErbB_kegg = gene_ErbB_kegg[gene_ErbB_kegg %in% V(g)$name] ## 81 genes
gene_GSK3_kegg = gene_GSK3_kegg[gene_GSK3_kegg %in% V(g)$name] ## 292 genes
gene_p53_kegg = gene_p53_kegg[gene_p53_kegg %in% V(g)$name] ## 71 genes

X_vals <- c("ErbB", "GSK3", "p53")
Y_vals <- c("kegg", "go")

## read file
for (x in X_vals) {
  for (y in Y_vals) {
    file_name <- paste0("gene_", x, "_", y, "5")
    assign(paste0("data_", file_name), read.table(paste0(file_name,".csv"), header = TRUE, sep = ","))
  }
}

wb <- createWorkbook()

# Loop through X and Y
for (X in X_vals) {
  for (Y in Y_vals) {
    # Construct the dynamic variable names for data and gene lists
    data_gene_name <- paste0("data_gene_", X, "_", Y, "5") # Assuming the '5' suffix is consistent
    gene_name <- paste0("gene_", X, "_", Y)
    
    # Check if the data exist 
    if (!exists(data_gene_name) || !exists(gene_name)) {
      cat(paste0("Skipping combination: ", X, "_", Y, " - Missing '", data_gene_name, "' or '", gene_name, "' in environment.\n"))
      next # Skip to the next iteration if data is missing
    }
    
    cat(paste0("\n--- Processing combination: ", X, "_", Y, " ---\n"))
    
    # Get the actual data frames
    data_gene_current <- get(data_gene_name)
    gene_current <- get(gene_name)
    
    # Initialize list to store 100 sets of genes from setdiff
    d_list_current <- vector("list", 100)
    
    # Loop through each NLC_i
    for (i in 1:100) {
      col_name <- paste0("NLC", i)
      # Ensure the column exists in the data_gene_current dataframe
      if (!col_name %in% colnames(data_gene_current)) {
        warning(paste0("Column '", col_name, "' not found in '", data_gene_name, "'. Skipping for this combination and NLC iteration.\n"))
        d_list_current[[i]] <- character(0) # Assign empty character vector to avoid errors
        next
      }
      values <- data_gene_current[[col_name]]
      values <- values[!is.na(values)]  # Remove NA
      d_list_current[[i]] <- setdiff(values, gene_current)
    }
    
    # name the list
    names(d_list_current) <- paste0("NLC", 1:100)
    
    # Combine all genes from the list into a single vector
    all_genes <- unlist(d_list_current)
    
    # Prepare an empty data frame
    if (length(all_genes) == 0) {
      cat("No genes found after setdiff for this combination. Creating an empty data frame.\n")
      gene_counts_df <- data.frame(Gene = character(0), Appearance_Count = integer(0))
    } else {
      # Count the frequency of each gene
      gene_counts <- table(all_genes)
      
      # Convert to a data frame
      gene_counts_df <- as.data.frame(gene_counts)
      colnames(gene_counts_df) <- c("Gene", "Appearance_Count")
      
      # Order the results
      gene_counts_df <- gene_counts_df[order(gene_counts_df$Appearance_Count, decreasing = TRUE), ]
      
      # top 10 most frequent genes
      cat("Top 10 most frequent genes:\n")
      print(head(gene_counts_df, 10))
    }
    
    # Save to a data frame 
    df_name <- paste0("hit_df_", X, "_", Y)
    assign(df_name, gene_counts_df, envir = .GlobalEnv) # Assign to global environment
    cat(paste0("Results saved to data frame: ", df_name, "\n"))
    
    # Add the gene_counts_df as a new sheet
    sheet_name <- paste0(X, "_", Y) # e.g., "ErbB_kegg"
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = gene_counts_df)
    cat(paste0("Results added to Excel sheet: '", sheet_name, "'\n"))
  }
}

# Save workbook
excel_file_name <- "all_hit_counts.xlsx"
saveWorkbook(wb, excel_file_name, overwrite = TRUE)
cat(paste0("\nAll results saved to '", excel_file_name, "' with multiple sheets.\n"))

########### for check ONLY ######################

for (x in X_vals) {
  for (y in Y_vals) {
    df_name <- paste0("hit_df_", x, "_", y)
    gene_set <- gene_list[[x]][[y]]
    
    if (exists(df_name)) {
      hit_df <- get(df_name)
      cat("Checking:", df_name, "\n")
      print(table(hit_df$Gene %in% gene_set))
      cat("\n")
    } else {
      cat("Missing:", df_name, "\n\n")
    }
  }
}

