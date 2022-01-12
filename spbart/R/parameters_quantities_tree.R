# -------------------------------------------------------------------------#
# Description: this script contains functions that are used to generate    #
#              the predictions, update variance and compute the tree prior #
#              and the marginalised likelihood                             #
# -------------------------------------------------------------------------#

# 1. simulate_mu: generate the predicted values (mu's)
# 2. update_sigma2: updates the parameters sigma2
# 3. update_z: updates the latent variables z. This is required for MOTR-BART for classification.
# 4. get_tree_prior: returns the tree log prior score
# 5. tree_full_conditional: computes the marginalised likelihood for all nodes for a given tree
# 6. get_number_distinct_cov: counts the number of distinct covariates that are used in a tree to create the splitting rules
# Compute the full conditionals -------------------------------------------------

tree_full_conditional = function(tree, R, sigma2, sigma2_mu, common_vars) {

  # Function to compute log full conditional distirbution for an individual tree
  # R is a vector of partial residuals

  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  # Identify those terminals that don't have a double split or are stumps
  if(nrow(tree$tree_matrix) != 1) {
    terminal_ancestors = get_ancestors(tree) # get the ancestor for all terminal nodes
    aux_table = table(terminal_ancestors[,1], terminal_ancestors[,2]) # create a table
    which_terminals_one_split = which(apply(aux_table,1,sum) == 1) # Find terminals which have only one covariate as ancestor
    get_index = apply(aux_table[which_terminals_one_split, ],1, function(x) which(x == 1)) # get the index of the column associated to the splitting rule
    check_common_vars = get_index %in% which(colnames(aux_table) %in% common_vars) # check whether the covariate is common to X1 and X2
    which_terminal_no_double_split = names(get_index)[check_common_vars] # terminals with only one ancestor where the ancestor is common to X1 and X2

  } else{
    which_terminal_no_double_split = 1 # stump
  }

  # set up sigma2_mu = 0 for all which_terminal_no_double_split
  sigma2_mu_aux = rep(sigma2_mu, length(which_terminal))
  sigma2_mu_aux[which(which_terminal %in% which_terminal_no_double_split)] = 0

  # Get node sizes for each terminal node
  nj = as.numeric(tree$tree_matrix[which_terminal,'node_size'])

  # Get sum of residuals and sum of residuals squared within each terminal node
  sumRsq_j = aggregate(R, by = list(tree$node_indices), function(x) sum(x^2))[,2]
  S_j = aggregate(R, by = list(tree$node_indices), sum)[,2]

  # Now calculate the log posterior
  log_post = 0.5 * ( sum(log( sigma2 / (nj*sigma2_mu_aux + sigma2))) +
                       sum( (sigma2_mu_aux* S_j^2) / (sigma2 * (nj*sigma2_mu_aux + sigma2))))
  return(log_post)
}

simulate_mu = function(tree, R, sigma2, sigma2_mu, common_vars) {

  # Simulate mu values for a given tree

  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  # Identify those terminals that don't have a double split on g and e
  if(nrow(tree$tree_matrix) != 1) {
    terminal_ancestors = get_ancestors(tree) # get the ancestor for all terminal nodes
    aux_table = table(terminal_ancestors[,1], terminal_ancestors[,2]) # create a table
    which_terminals_one_split = which(apply(aux_table,1,sum) == 1) # Find terminals which have only one covariate as ancestor
    get_index = apply(aux_table[which_terminals_one_split, ],1, function(x) which(x == 1)) # get the index of the column associated to the splitting rule
    check_common_vars = get_index %in% which(colnames(aux_table) %in% common_vars) # check whether the covariate is common to X1 and X2
    which_terminal_no_double_split = as.numeric(names(get_index)[check_common_vars]) # terminals with only one ancestor where the ancestor is common to X1 and X2

  } else{
    which_terminal_no_double_split = 1 # stump
  }
  # set up sigma2_mu = 0 for all which_terminal_no_double_split
  sigma2_mu_aux = rep(sigma2_mu, length(which_terminal))
  sigma2_mu_aux[which(which_terminal %in% which_terminal_no_double_split)] = 0

  # Get node sizes for each terminal node
  nj = as.numeric(tree$tree_matrix[which_terminal,'node_size'])

  # Get sum of residuals in each terminal node
  sumR = aggregate(R, by = list(tree$node_indices), sum)[,2]

  # Now calculate mu values
  mu = rnorm(length(nj) ,
             mean = (sumR / sigma2) / (nj/sigma2 + sigma2_mu_aux),
             sd = sqrt(1/(nj/sigma2 + sigma2_mu_aux)))

  # Wipe all the old mus out for other nodes
  tree$tree_matrix[,'mu'] = NA

  # Put in just the ones that are useful
  tree$tree_matrix[which_terminal,'mu'] = mu
  tree$tree_matrix[which_terminal_no_double_split, 'mu'] = 0 # set to zero the terminal node with no interaction which descends straight from the root node

  return(tree)
}


# Update sigma2 -------------------------------------------------------------

update_sigma2 <- function(S, n, nu, lambda){
  u = 1/rgamma(1, shape = (n + nu)/2, rate = (S + nu*lambda)/2)
  return(u)
}

# Update the latent variable z ---------------

update_z = function(y, prediction){

  ny0 = sum(y==0)
  ny1 = sum(y==1)
  z = rep(NA, length(y))

  z[y==0] = rtruncnorm(ny0, a = -Inf, b=0,   mean = prediction[y==0], 1)
  z[y==1] = rtruncnorm(ny1, a = 0   , b=Inf, mean = prediction[y==1], 1)

  return(z)
}

# Get tree priors ---------------------------------------------------------

get_tree_prior = function(tree, alpha, beta, common_vars) {

  # First find the level of each node, then the depth is the maximum of the level
  level = rep(NA, nrow(tree$tree_matrix))
  level[1] = 0 # First row always level 0

  # Escpae quickly if tree is just a stump
  if(nrow(tree$tree_matrix) == 1) {
    return(log(1 - alpha)) # Tree depth is 0
  }

  for(i in 2:nrow(tree$tree_matrix)) {
    # Find the current parent
    curr_parent = as.numeric(tree$tree_matrix[i,'parent'])
    # This child must have a level one greater than it's current parent
    level[i] = level[curr_parent] + 1
  }

  # Only compute for the internal nodes
  internal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 0)
  log_prior = 0
  for(i in 1:length(internal_nodes)) {
    log_prior = log_prior + log(alpha) - beta * log(1 + level[internal_nodes[i]])
  }
  # Now add on terminal nodes
  terminal_nodes = which(as.numeric(tree$tree_matrix[,'terminal']) == 1)
  for(i in 1:length(terminal_nodes)) {
    log_prior = log_prior + log(1 - alpha * ((1 + level[terminal_nodes[i]])^(-beta)))
  }

  return(log_prior)

  }

update_beta <- function(y, x, sigma2, omega_inv){

  sigma_beta = solve(t(x)%*%x / sigma2 + omega_inv)
  mu_beta = sigma_beta%*%(t(x)%*%y/sigma2)

  beta_hat = rmvnorm(1,
                     mean = mu_beta,
                     sigma = sigma_beta)

  return(t(beta_hat))
}

update_omega <- function(beta_hat, b, V, v){
  out = riwish(v + 1, (beta_hat - b)%*%t(beta_hat - b) + V)
  return(out)
}