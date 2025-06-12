Network diffusion and Neighborhood method for pathway enrichment
analysis
================
PJ
2025-06-12

## 1. PPI network construction from STRINGdb

##### Loading Required Libraries

``` r
library(STRINGdb)
library(igraph)
library(diffuStats)
library(gprofiler2)
```

##### Configuration

- Set the STRINGdb version to 12.0 (latest version)
- Specify the species as Homo sapiens (9606)
- Set the interaction confidence score threshold to 900 to include only
  high-confidence interactions
- Load all gene data to ensure a comprehensive analysis

``` r
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 900)
all_genes <- string_db$get_proteins()
```

##### Retrieve Interactions

- Extract interaction data from STRINGdb
- Merge the interaction data with gene symbol annotations

``` r
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
```

| Gene_A | Gene_B    | Combined_Score |
|:-------|:----------|---------------:|
| XYLT2  | B4GALT7   |            979 |
| XYLT2  | B4GALT7   |            979 |
| RB1CC1 | GABARAPL2 |            953 |
| RB1CC1 | GABARAPL2 |            953 |
| MRPS10 | MRPS35    |            999 |
| MRPS10 | MRPS35    |            999 |

##### Graph Construction

- Convert the interaction data into a graph structure
- Use the combined score as edge weights

``` r
g <- graph_from_data_frame(d = interaction_network, directed = FALSE, vertices = NULL)
components <- igraph::clusters(g, mode="weak")
largest_cluster <- which.max(components$csize)
large_nodes <- V(g)[components$membership == largest_cluster]
g <- induced_subgraph(g, large_nodes)
g <- simplify(g)
is_connected(g)
```

    ## [1] TRUE

## 2. Seed nodes Assignment and Expansion

##### Selecting Seed Nodes

- Seed nodes (Genes of interest) can come from any prior analysis in the
  gene symbol format
- For example,

``` r
seed_genes <- c("TP53", "BRCA1", "EGFR", "AKT1", "PIK3CA", 
                "MTOR", "PTEN", "CDK2", "KRAS", "BRAF")
```

#### **Neighborhood expansion**

Neighborhood method considers seed genes and seed genes’ neighbors.

``` r
neighbors_list <- unique(unlist(lapply(seed_genes, function(gene) {
  V(g)$name[neighbors(g, gene, mode = "all")]
})))
neighbors_list = unique(c(neighbors_list, seed_genes))
length(neighbors_list)
```

    ## [1] 892

#### **Network Diffusion**

Network diffusion spreads values from the initial seed nodes across the
network. The matrix $D$ represents the sum of each row in the adjacency
matrix $A$, and the graph Laplacian is given by $𝐿=𝐷−𝑊$, where $W$ is
the adjacency matrix. The Laplacian kernel and normalized Laplacian are
calculated as $𝐾_L=(𝐼+𝛼𝐿)^{−1}$ and
$K_{NL}=(I+α(I-D^{-1/2} AD^{-1/2}))^{-1}$ , where $I$ is the identity
matrix and $α$ is the diffusion rate, set to 1 in this study. The
diffusion score is calculated as:

$$ Diffusion=K⋅y $$ where $𝑦$ is a vector with a value of 1 for nodes in
seed nodes and 0 for others. Nodes that are more connected to the seed
nodes will have higher diffusion scores, meaning they are more
influenced by the seed group.

First, we use diffuStats package to compute normalized regularize
laplacian kernel ($K_{NL}$).

``` r
K_nl <- regularisedLaplacianKernel(g,sigma2=1, add_diag=1,normalized=TRUE)
```

Next, the vector $y$ is computed to obtain diffusion score.

``` r
y <- ifelse(V(g)$name %in% seed_genes, 1, 0)
```

Note that the network diffusion method assigns a score to all genes in
the network. We need to select the size of the top scores, for example,
by choosing the 99th percentile. This means we consider the genes whose
diffusion scores are in the top 1% of all scores, focusing on the most
strongly influenced genes in the network.

``` r
size = length(y)*0.01
F_NL <- K_nl %*% y
F_NL <- data.frame(F_NL)

F_NL$name = rownames(F_NL)
F_NL = F_NL[order(F_NL$F_NL, decreasing = TRUE), ]
F_NL_select = head(F_NL$name, size)
F_NL_select
```

    ##   [1] "PTEN"         "TP53"         "KRAS"         "PIK3CA"       "AKT1"        
    ##   [6] "BRAF"         "EGFR"         "MTOR"         "CDK2"         "BRCA1"       
    ##  [11] "KIAA1549"     "NPAT"         "K7ERJ3_HUMAN" "MAST2"        "PREX2"       
    ##  [16] "NF1"          "MAPKAP1"      "SPRY2"        "HRAS"         "PIK3CD"      
    ##  [21] "CDKN2A"       "PIK3R3"       "PIK3CB"       "PIK3R2"       "PIK3R1"      
    ##  [26] "NRAS"         "PI3"          "SOS1"         "AKT3"         "AKT2"        
    ##  [31] "BRAT1"        "CCNO"         "ZPR1"         "ERRFI1"       "RAF1"        
    ##  [36] "FKBP3"        "EML4"         "PDGFRB"       "ERBB2"        "RPS6KB2"     
    ##  [41] "MRAS"         "BRAP"         "RRAS2"        "PIK3R6"       "PIK3R5"      
    ##  [46] "ARAF"         "IRS1"         "SPDYC"        "SPDYE4"       "SPDYE3"      
    ##  [51] "SPDYE2"       "SPDYE2B"      "SPDYE6"       "SPDYE16"      "PRKCZ"       
    ##  [56] "SPDYE1"       "SPDYE5"       "PIK3CG"       "RASGRP4"      "RRAS"        
    ##  [61] "ERBB3"        "HSP90AB1"     "PDGFRA"       "RICTOR"       "HSP90AA1"    
    ##  [66] "YWHAZ"        "CDC37"        "CCND1"        "PRKCA"        "FRK"         
    ##  [71] "ARID3A"       "E4F1"         "ZNF385A"      "PCBP4"        "TP53AIP1"    
    ##  [76] "AIFM2"        "S100A2"       "GPNMB"        "PLCH1"        "PLCH2"       
    ##  [81] "RPS6KB1"      "MLST8"        "FKBP5"        "MAP2K2"       "CDK4"        
    ##  [86] "MYC"          "LRIG1"        "PRKCG"        "PIP5KL1"      "AKT1S1"      
    ##  [91] "MET"          "INPPL1"       "KSR1"         "UBASH3B"      "SPDYE11"     
    ##  [96] "CDKN1A"       "CDC25C"       "TNS3"         "GRB2"         "PRKCB"       
    ## [101] "RPTOR"        "FOXO4"        "RHEB"         "YWHAG"        "EIF4EBP1"    
    ## [106] "CHUK"         "MAPK1"        "GSK3B"        "INPP5D"       "ESR1"        
    ## [111] "MAP2K1"       "MAPK3"        "SPDYA"        "FKBP4"        "IGF1R"       
    ## [116] "EIF4EBP2"

## 3. Pathway enrichment analysis

Pathway enrichment analysis is performed using the gprofiler2 package.
Here, the analysis uses the results from the network diffusion scores,
focusing on pathways from “<GO:BP>” (Gene Ontology Biological
Processes), “KEGG” (Kyoto Encyclopedia of Genes and Genomes), “REAC”
(Reactome), and “WP” (WikiPathways).

``` r
result <- gost(query = F_NL_select,
                 organism = "hsapiens",  
                 sources = c("GO:BP", "KEGG", "REAC", "WP"))
table <- result$result
head(table)
```

| query | significant | p_value | term_size | query_size | intersection_size | precision | recall | term_id | source | term_name | effective_domain_size | source_order | parents |
|:---|:---|---:|---:|---:|---:|---:|---:|:---|:---|:---|---:|---:|:---|
| query_1 | TRUE | 0 | 1341 | 112 | 72 | 0.6428571 | 0.0536913 | <GO:0006468> | <GO:BP> | protein phosphorylation | 21017 | 2207 | <GO:0016310>, <GO:0036211> |
| query_1 | TRUE | 0 | 1575 | 112 | 74 | 0.6607143 | 0.0469841 | <GO:0016310> | <GO:BP> | phosphorylation | 21017 | 5198 | <GO:0006796> |
| query_1 | TRUE | 0 | 2550 | 112 | 79 | 0.7053571 | 0.0309804 | <GO:0006796> | <GO:BP> | phosphate-containing compound metabolic process | 21017 | 2482 | <GO:0006793> |
| query_1 | TRUE | 0 | 2553 | 112 | 79 | 0.7053571 | 0.0309440 | <GO:0006793> | <GO:BP> | phosphorus metabolic process | 21017 | 2479 | <GO:0008152> |
| query_1 | TRUE | 0 | 901 | 112 | 52 | 0.4642857 | 0.0577137 | <GO:0001932> | <GO:BP> | regulation of protein phosphorylation | 21017 | 546 | <GO:0006468>, <GO:0031399>, <GO:0042325> |
| query_1 | TRUE | 0 | 2946 | 112 | 77 | 0.6875000 | 0.0261371 | <GO:0035556> | <GO:BP> | intracellular signal transduction | 21017 | 9090 | <GO:0007165> |
