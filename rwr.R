
rwr <- function(transition_matrix, start_genes, restart_prob = 0.3, max_iter = 100, tol = 1e-6) {
  # Initialize restart vector (starting probabilities)
  start_indices <- match(start_genes, rownames(transition_matrix))
  if (any(is.na(start_indices))) stop("Invalid gene names in start_genes!")
  
  restart_vec <- numeric(nrow(transition_matrix))
  restart_vec[start_indices] <- 1 / length(start_indices)  # Equal probability for multiple start nodes
  
  # Initialize scores (current state of the walk)
  scores <- restart_vec
  
  for (i in seq_len(max_iter)) {
    # Update scores: (1 - r) * Transition + r * Restart
    new_scores <- (1 - restart_prob) * (transition_matrix %*% scores) + restart_prob * restart_vec
    
    # Check convergence
    if (sum(abs(new_scores - scores)) < tol) {
      break
    }
    scores <- new_scores
  }
  
  return(scores)
}



