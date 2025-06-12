# Load required libraries
library(diffuStats)
library(igraph)
library(pROC)
library(PRROC)

# Load data
load("K_laplacian.Rdata")
load("K_normlap.Rdata")
load("g.Rdata")
load("transition_matrix.Rdata")

# Load function
source("rwr.R")
kegg <- read.csv("kegg_pathway.csv", header = TRUE)

# Parameters
n_iterations <- 100

# Main loop over pathways
for (k in 1:170) {
  print(k)
  
  # Extract pathway data
  pathway_id <- kegg$pathway_id[k]
  gene_symbols <- strsplit(kegg$gene_symbols[k], ",\\s*")[[1]]
  namestrue <- intersect(gene_symbols, V(g)$name)
  test <- as.numeric(V(g)$name %in% namestrue)
  names(test) <- V(g)$name
  num_test_data <- sum(test)
  
  if (num_test_data == 0) next  # Skip if no matching genes
  
  # Set iteration limits
  max_iterations <- if (num_test_data > 100) 100 else num_test_data
  percentages <- round(seq(1, 100, length.out = max_iterations) * num_test_data / 100)
  
  # Initialize results dataframe
  results <- data.frame(
    index = integer(max_iterations),
    mean_auc_lc = numeric(max_iterations),
    var_auc_lc = numeric(max_iterations),
    mean_aupr_lc = numeric(max_iterations),
    var_aupr_lc = numeric(max_iterations),
    mean_auc_ht = numeric(max_iterations),
    var_auc_ht = numeric(max_iterations),
    mean_aupr_ht = numeric(max_iterations),
    var_aupr_ht = numeric(max_iterations),
    mean_auc_rwr = numeric(max_iterations),
    var_auc_rwr = numeric(max_iterations),
    mean_aupr_rwr = numeric(max_iterations),
    var_aupr_rwr = numeric(max_iterations),
    mean_auc_neigh = numeric(max_iterations),
    var_auc_neigh = numeric(max_iterations),
    mean_aupr_neigh = numeric(max_iterations),
    var_aupr_neigh = numeric(max_iterations),
    mean_auc_neigh2 = numeric(max_iterations),
    var_auc_neigh2 = numeric(max_iterations),
    mean_aupr_neigh2 = numeric(max_iterations),
    var_aupr_neigh2 = numeric(max_iterations),
    percentage_of_seed = numeric(max_iterations),
    num_test_data = integer(max_iterations)
  )
  
  # Iteration over seed percentages
  for (j in seq_len(max_iterations)) {
    num_seeds <- if (num_test_data > 100) percentages[j] else j
    auc_lc_values <- numeric(n_iterations)
    auc_ht_values <- numeric(n_iterations)
    auc_rwr_values <- numeric(n_iterations)
    auc_neigh_values <- numeric(n_iterations)
    auc_neigh2_values <- numeric(n_iterations)
    aupr_lc_values <- numeric(n_iterations)
    aupr_ht_values <- numeric(n_iterations)
    aupr_rwr_values <- numeric(n_iterations)
    aupr_neigh_values <- numeric(n_iterations)
    aupr_neigh2_values <- numeric(n_iterations)
    
    
    # Randomized iterations
    for (i in seq_len(n_iterations)) {
      seed_gene <- sample(namestrue, num_seeds)
      y <- as.numeric(V(g)$name %in% seed_gene)
      
      # Diffusion scores
      F_diffusion <- K %*% y
      F_dif <- F_diffusion[, 1]
      
      F_heat <- K2 %*% y
      F_ht <- F_heat[, 1]
      
      rwr_scores <- rwr(transition_matrix, start_genes = seed_gene, restart_prob = 0.3)
      
      # Neighborhood lists
      neighbors_list <- unique(unlist(neighbors(g, seed_gene, mode = "all")))
      neighbors_list2 <- unique(unlist(neighbors(g, neighbors_list, mode = "all")))
      
      # Scores for different methods
      neigh_scores <- as.numeric(V(g)$name %in% c(seed_gene, as_ids(neighbors_list)))
      neigh_scores2 <- as.numeric(V(g)$name %in% c(seed_gene, as_ids(neighbors_list), as_ids(neighbors_list2)))
      
      # Calculate AUC and AUPR
      roc_lc <- roc(test, F_dif)
      auc_lc_values[i] <- auc(roc_lc)
      aupr_lc_values[i] <- pr.curve(scores.class0 = F_dif[test == 1], scores.class1 = F_dif[test == 0], curve = TRUE)$auc.integral
      
      roc_ht <- roc(test, F_ht)
      auc_ht_values[i] <- auc(roc_ht)
      aupr_ht_values[i] <- pr.curve(scores.class0 = F_ht[test == 1], scores.class1 = F_ht[test == 0], curve = TRUE)$auc.integral
      
      roc_rwr <- roc(test, rwr_scores)
      auc_rwr_values[i] <- auc(roc_rwr)
      aupr_rwr_values[i] <- pr.curve(scores.class0 = rwr_scores[test == 1], scores.class1 = rwr_scores[test == 0], curve = TRUE)$auc.integral
      
      roc_neigh <- roc(test, neigh_scores)
      auc_neigh_values[i] <- auc(roc_neigh)
      aupr_neigh_values[i] <- pr.curve(scores.class0 = neigh_scores[test == 1], scores.class1 = neigh_scores[test == 0], curve = TRUE)$auc.integral
      
      roc_neigh2 <- roc(test, neigh_scores2)
      auc_neigh2_values[i] <- auc(roc_neigh2)
      aupr_neigh2_values[i] <- pr.curve(scores.class0 = neigh_scores2[test == 1], scores.class1 = neigh_scores2[test == 0], curve = TRUE)$auc.integral
    }
    
    # Store results
    results[j, ] <- data.frame(
      index = j,
      mean_auc_lc = mean(auc_lc_values),
      var_auc_lc = var(auc_lc_values),
      mean_aupr_lc = mean(aupr_lc_values),
      var_aupr_lc = var(aupr_lc_values),
      mean_auc_ht = mean(auc_ht_values),
      var_auc_ht = var(auc_ht_values),
      mean_aupr_ht = mean(aupr_ht_values),
      var_aupr_ht = var(aupr_ht_values),
      mean_auc_rwr = mean(auc_rwr_values),
      var_auc_rwr = var(auc_rwr_values),
      mean_aupr_rwr = mean(aupr_rwr_values),
      var_aupr_rwr = var(aupr_rwr_values),
      mean_auc_neigh = mean(auc_neigh_values),
      var_auc_neigh = var(auc_neigh_values),
      mean_aupr_neigh = mean(aupr_neigh_values),
      var_aupr_neigh = var(aupr_neigh_values),
      mean_auc_neigh2 = mean(auc_neigh2_values),
      var_auc_neigh2 = var(auc_neigh2_values),
      mean_aupr_neigh2 = mean(aupr_neigh2_values),
      var_aupr_neigh2 = var(aupr_neigh2_values),
      percentage_of_seed = (num_seeds / num_test_data) * 100,
      num_test_data = num_test_data
    )
  }
  
  # Save results
  pathway_id_no_colon <- gsub(":", "_", pathway_id)
  #write.csv(results, paste0(pathway_id_no_colon, "_kegg_new.csv"), row.names = FALSE)
  write.csv(results, paste0(pathway_id_no_colon, "_kegg_5method.csv"), row.names = FALSE)
}

