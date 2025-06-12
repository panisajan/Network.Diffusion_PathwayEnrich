setwd("C:/Users/Asus/Desktop/AomamxAsus/diffusion_neighbor/")

# Load required libraries
library(diffuStats)
library(igraph)
library(pROC)
library(PRROC)
library(STRINGdb)
library(gprofiler2)

# Load data
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 900)
all_genes <- string_db$get_proteins()

# kegg <- read.csv("kegg_pathway.csv", header = TRUE)

data <- data.frame(gene = all_genes$preferred_name)
mapped_data <- string_db$map(data, "gene", removeUnmappedRows = FALSE)

interaction_network <- string_db$get_interactions(mapped_data$STRING_id)
interaction_network <- merge(interaction_network, all_genes, 
                             by.x = "from", by.y = "protein_external_id", all.x = TRUE)
interaction_network <- merge(interaction_network, all_genes, 
                             by.x = "to", by.y = "protein_external_id", all.x = TRUE, 
                             suffixes = c("_A", "_B"))
interaction_network <- interaction_network[, c("preferred_name_A", "preferred_name_B", "combined_score")]
colnames(interaction_network) <- c("Gene_A", "Gene_B", "Combined_Score")
head(interaction_network)


g <- graph_from_data_frame(d = interaction_network, directed = FALSE, vertices = NULL)
components <- igraph::clusters(g, mode="weak")
largest_cluster <- which.max(components$csize)
large_nodes <- V(g)[components$membership == largest_cluster]
g <- induced_subgraph(g, large_nodes)
is_connected(g)

g <- simplify(g)

K <- regularisedLaplacianKernel(g,sigma2=1, add_diag=1,normalized=FALSE)
# save(K,file = "K_laplacian.Rdata")

K2 <- regularisedLaplacianKernel(g,sigma2=1,normalized=TRUE)
# save(K2,file = "K_normlap.Rdata")
# save(g,file = "g.Rdata")

# load("g.Rdata")

adj_matrix = as_adjacency_matrix(g)
# Normalize the adjacency matrix to create a transition matrix
transition_matrix <- t(apply(adj_matrix, 1, function(row) row / sum(row)))

# save(transition_matrix,file = "transition_matrix.Rdata")
