#' Combined Semi-Parametric BART
#'
#' @param y x
#' @param x x
#' @param sparse x
#' @param ntrees x
#' @param node_min_size x
#' @param alpha x
#' @param beta x
#' @param nu x
#' @param lambda x
#' @param mu_mu x
#' @param sigma2 x
#' @param sigma2_mu x
#' @param nburn x
#' @param npost x
#' @param nthin x
#'
#' @return x
#' @importFrom stats 'rgamma' 'rexp' 'dnorm' 'sd' 'rchisq' 'rnorm' 'pnorm' 'as.formula' 'terms' 'xtabs'
#' @importFrom truncnorm 'rtruncnorm'
#' @importFrom lme4 'lFormula'
#' @importFrom dbarts 'makeModelMatrixFromDataFrame'
#' @export
#'
#' @examples
#' #
#'
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
bart = function(   y,
                   x, # it needs to contain the response
                   sparse = FALSE,
                   ntrees = 50,
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


  aux_identify_factor_variables = NULL
  common_variables = NULL

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
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)
  current_partial_residuals = y_scale

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = x)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  yhat_bart = get_predictions(curr_trees, x, single_tree = ntrees == 1)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in seq_len(TotIter)) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = y_hat
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      bart_store[curr,] = yhat_bart
    }

    # Start looping through trees
    for (j in seq_len(ntrees)) {

      current_partial_residuals = y_scale - yhat_bart + tree_fits_store[,j]

      # Propose a move (grow, prune, change, or swap)
      type = sample_move(curr_trees[[j]], i, nburn)

      # Generate a new tree based on the current
      new_trees[[j]] = update_tree(y = current_partial_residuals,
                                   X = x,
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s,
                                   common_vars = common_variables,
                                   aux_factor_var = aux_identify_factor_variables)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu,
                                    common_variables,
                                    aux_identify_factor_variables) +
        get_tree_prior(curr_trees[[j]], alpha, beta, common_variables)

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu,
                                    common_variables,
                                    aux_identify_factor_variables) +
        get_tree_prior(new_trees[[j]], alpha, beta, common_variables)

      # Exponentiate the results above
      if(isTRUE(new_trees[[j]]$ForceStump)) {a=1} else {a = l_new - l_old}

      if(a > 0 || a > -rexp(1)) {
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
                                    sigma2_mu,
                                    common_variables,
                                    aux_identify_factor_variables)
      # Updating BART predictions
      current_fit = get_predictions(curr_trees[j], x, single_tree = TRUE)
      yhat_bart = yhat_bart - tree_fits_store[,j] # subtract the old fit
      yhat_bart = yhat_bart + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Updating the final predictions
    y_hat = yhat_bart

    sum_of_squares = sum((y_scale - y_hat)^2)

    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if(isTRUE(sparse) & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  results <- list(trees = tree_store,
                  sigma2 = sigma2_store*y_sd^2,
                  y_hat = y_hat_store*y_sd + y_mean,
                  bart_hat = bart_store*y_sd,
                  npost = npost,
                  nburn = nburn,
                  nthin = nthin,
                  ntrees = ntrees,
                  y_mean = y_mean,
                  y_sd = y_sd,
                  var_count_store = var_count_store,
                  s = s_prob_store)

  class(results) <- "cspbart"
  return(results)
} # End main function
