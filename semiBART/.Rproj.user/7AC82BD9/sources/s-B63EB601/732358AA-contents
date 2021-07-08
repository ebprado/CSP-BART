#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'as.formula' 'terms'
#' @importFrom MCMCpack 'rdirichlet' 'riwish'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom lme4 'lFormula'
#' @importFrom dbarts 'makeModelMatrixFromDataFrame'
#'
# x1 = x[,4:5]
# x2 = x
# y = y
# sparse = FALSE
# ntrees = 50
# node_min_size = 5
# alpha = 0.95
# beta = 2
# nu = 3
# lambda = 0.1
# mu_mu = 0
# sigma2 = 1
# sigma2_mu = 1
# nburn = 100
# npost = 100
# nthin = 1


semibart = function(formula,
                   x1, # it needs to contain the response
                   x2, # it doesn't need to contain the response
                   sparse = FALSE,
                   ntrees = 10,
                   node_min_size = 5,
                   alpha = 0.95,
                   beta = 2,
                   nu = 3,
                   lambda = 0.1,
                   mu_mu = 0,
                   sigma2 = 1,
                   sigma2_mu = 1,
                   nburn = 1000,
                   npost = 1000,
                   nthin = 1) {

  if (class(x1) != 'data.frame' || class(x2) != 'data.frame') {stop('X1 and X2 need to be data frames.')}

  # formula = as.formula(formula)
  data = MakeDesignMatrix(formula, x1)
  y = data$y
  x1 = as.matrix(data$X) # matrix to be used in the linear predictor
  x2 = makeModelMatrixFromDataFrame(x2, drop = FALSE) # matrix to be used in the BART component

  colnames_x1 = colnames(x1)
  colnames_x2 = colnames(x2)

  common_variables = which(colnames_x2%in%colnames_x1)

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  bart_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x2))
  var_count_store = matrix(0, ncol = ncol(x2), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x2), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))
  beta_store = matrix(NA, ncol = ncol(x1), nrow=store_size)
  colnames(beta_store) = colnames_x1

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p1 = ncol(x1)
  p2 = ncol(x2)
  p = p1 + p2
  s = rep(1/p2, p2)
  Omega = diag(p1)
  Omega_inv = solve(Omega)
  b = rep(0, p1)
  V = diag(p1)
  v = p1
  beta_hat = rep(0, p1)
  current_partial_residuals = y_scale

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = x2)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  yhat_bart = get_predictions(curr_trees, x2, single_tree = ntrees == 1)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = y_hat
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      beta_store[curr,] = beta_hat
      bart_store[curr,] = yhat_bart
    }

    # Update linear predictor -------
    beta_hat = update_beta(y_scale - yhat_bart, x1, sigma2, Omega_inv)
    yhat_linear = x1%*%beta_hat

    # Update covariance matrix of the linear predictor
    Omega = update_omega(beta_hat, b, V, v)
    Omega_inv = solve(Omega)

      # Start looping through trees
      for (j in 1:ntrees) {

        current_partial_residuals = y_scale - yhat_bart - yhat_linear + tree_fits_store[,j]

        # Propose a move (grow, prune, change, or swap)
        type = sample_move(curr_trees[[j]], i, nburn, common_variables)

        # Generate a new tree based on the current
        new_trees[[j]] = update_tree(y = current_partial_residuals,
                                     X = x2,
                                     type = type,
                                     curr_tree = curr_trees[[j]],
                                     node_min_size = node_min_size,
                                     s = s,
                                     common_vars = common_variables)

        # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_old = tree_full_conditional(curr_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu) +
          get_tree_prior(curr_trees[[j]], alpha, beta, common_variables)

        # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
        l_new = tree_full_conditional(new_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu) +
          get_tree_prior(new_trees[[j]], alpha, beta, common_variables)

        # Exponentiate the results above
        a = exp(l_new - l_old)

        if(a > runif(1)) {
          curr_trees[[j]] = new_trees[[j]]

          if (type =='change'){
            var_count[curr_trees[[j]]$var[1]] = var_count[curr_trees[[j]]$var[1]] + 1
            var_count[curr_trees[[j]]$var[2]] = var_count[curr_trees[[j]]$var[2]] - 1
          }

          if (type=='grow'){
            var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] + 1 }

          if (type=='prune'){
            var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] - 1 }
        }

        # Update mu whether tree accepted or not
        curr_trees[[j]] = simulate_mu(curr_trees[[j]],
                                      current_partial_residuals,
                                      sigma2,
                                      sigma2_mu)
      # Updating BART predictions
      current_fit = get_predictions(curr_trees[j], x2, single_tree = TRUE)
      yhat_bart = yhat_bart - tree_fits_store[,j] # subtract the old fit
      yhat_bart = yhat_bart + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

      } # End loop through trees

    # Updating the final predictions
    y_hat = yhat_linear + yhat_bart

    sum_of_squares = sum((y_scale - y_hat)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

    cat('\n') # Make sure progress bar ends on a new line

    if (colnames_x1[1]!="(Intercept)") {beta_hat = beta_store*y_sd}
    if (colnames_x1[1]=="(Intercept)") {
      beta_hat = beta_store*y_sd
       beta_hat[,1] = beta_hat[,1] + y_mean
    }

  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat = y_hat_store*y_sd + y_mean,
              beta_hat = beta_hat,
              bart_hat = bart_store*y_sd,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              var_count_store = var_count_store,
              s = s_prob_store,
              formula = formula,
              colnames.x1 = colnames_x1,
              colnames.x2 = colnames_x2))

} # End main function


#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'as.formula' 'model.matrix'
#' @importFrom MCMCpack 'rdirichlet' 'riwish'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom lme4 'lFormula'
#' @importFrom dbarts 'makeModelMatrixFromDataFrame'
#'

cl_semibart = function(formula,
                    x1, # it needs to contain the response
                    x2, # it doesn't need to contain the response
                    sparse = FALSE,
                    ntrees = 10,
                    node_min_size = 5,
                    alpha = 0.95,
                    beta = 2,
                    nu = 3,
                    lambda = 0.1,
                    mu_mu = 0,
                    sigma2 = 1,
                    sigma2_mu = 1,
                    nburn = 1000,
                    npost = 1000,
                    nthin = 1) {

  if (class(x1) != 'data.frame' || class(x2) != 'data.frame') {stop('X1 and X2 need to be data frames.')}

  # formula = as.formula(formula)
  data = MakeDesignMatrix(formula, x1)

  y = data$y
  x1 = as.matrix(data$X) # matrix to be used in the linear predictor
  x2 = makeModelMatrixFromDataFrame(x2, drop = FALSE) # matrix to be used in the BART component

  colnames_x1 = colnames(x1)
  colnames_x2 = colnames(x2)

  common_variables = which(colnames_x2%in%colnames_x1)

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  bart_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x2))
  var_count_store = matrix(0, ncol = ncol(x2), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x2), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))
  beta_store = matrix(NA, ncol = ncol(x1), nrow=store_size)
  colnames(beta_store) = colnames_x1

  # Scale the response target variable

  n = length(y)
  p1 = ncol(x1)
  p2 = ncol(x2)
  p = p1 + p2
  s = rep(1/p2, p2)
  Omega = diag(p1)
  Omega_inv = solve(Omega)
  b = rep(0, p1)
  V = diag(p1)
  v = p1
  beta_hat = rep(0, p1)
  z = ifelse(y == 0, -3, 3)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = z,
                            X = x2)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  yhat_bart = get_predictions(curr_trees, x2, single_tree = ntrees == 1)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      y_hat_store[curr,] = pnorm(y_hat)
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      beta_store[curr,] = beta_hat
      bart_store[curr,] = yhat_bart
    }

    # Update linear predictor -------
    beta_hat = update_beta(z - yhat_bart, x1, sigma2, Omega_inv)
    yhat_linear = x1%*%beta_hat

    # Update covariance matrix of the linear predictor
    Omega = update_omega(beta_hat, b, V, v)
    Omega_inv = solve(Omega)

    # Start looping through trees
    for (j in 1:ntrees) {

      current_partial_residuals = z - yhat_bart - yhat_linear + tree_fits_store[,j]

      # Propose a move (grow, prune, change, or swap)
      type = sample_move(curr_trees[[j]], i, nburn, common_variables)

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = z,
                                   X = x2,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s,
                                   common_vars = common_variables)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu) +
        get_tree_prior(curr_trees[[j]], alpha, beta, common_variables)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu) +
        get_tree_prior(new_trees[[j]], alpha, beta, common_variables)

      # Exponentiate the results above
      if (type == 'stump') {a=1} else {a = exp(l_new - l_old)}

      if(a > runif(1)) {
        curr_trees[[j]] = new_trees[[j]]

        if (type =='change'){
          var_count[curr_trees[[j]]$var[1]] = var_count[curr_trees[[j]]$var[1]] + 1
          var_count[curr_trees[[j]]$var[2]] = var_count[curr_trees[[j]]$var[2]] - 1
        }

        if (type=='grow'){
          var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] + 1 }

        if (type=='prune'){
          var_count[curr_trees[[j]]$var] = var_count[curr_trees[[j]]$var] - 1 }
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_mu(curr_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu)
      # Updating BART predictions
      current_fit = get_predictions(curr_trees[j], x2, single_tree = TRUE)
      yhat_bart = yhat_bart - tree_fits_store[,j] # subtract the old fit
      yhat_bart = yhat_bart + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Updating the final predictions
    y_hat = yhat_linear + yhat_bart

    # Update z (latent variable)
    z = update_z(y, y_hat)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              y_hat = y_hat_store,
              beta_hat = beta_store,
              bart_hat = bart_store,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              var_count_store = var_count_store,
              s = s_prob_store,
              formula = formula,
              colnames.x1 = colnames_x1,
              colnames.x2 = colnames_x2))

} # End main function
