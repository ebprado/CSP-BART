#' @method predict "cspbart"
#' @rdname cspbart
#' @param object x
#' @param newdata_x1 x
#' @param newdata_x2 x
#' @param type x
#' @param ... Catches unused arguments.
#' @usage
#' \method{predict}{cspbart}(object,
#'         newdata_x1,
#'         newdata_x2,
#'         type = c('all', 'median', 'mean'),
#'         ...)
#' @export
predict.cspbart = function(object, newdata_x1, newdata_x2,
                         type = c('all', 'median', 'mean'), ...) {

  data_new = MakeDesignMatrixPredict(object$formula, newdata_x1)
  X1 = as.matrix(data_new$X) # matrix to be used in the linear predictor
  X2 = makeModelMatrixFromDataFrame(newdata_x2, drop = FALSE)
  beta_hat = colMeans(object$beta_hat)

  # Create holder for predicted values
  n_newX = dim(X2)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(X2))

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    X2,
                                    single_tree = length(curr_trees) == 1)
  }

  if (colnames(X1)[1] == '(Intercept)') {
    # Sort out what to return
    out = switch(type,
                 all = sweep(object$y_sd * y_hat_mat, 2, as.numeric(X1%*%beta_hat), "+"),
                 mean = X1%*%beta_hat + object$y_sd * colMeans(y_hat_mat),
                 median = X1%*%beta_hat + object$y_sd * apply(y_hat_mat,2,'median'))
  } else {

    out = switch(type,
           all = sweep(object$y_mean + object$y_sd * y_hat_mat, 2, as.numeric(X1%*%beta_hat), "+"),
           mean = X1%*%beta_hat + object$y_mean + object$y_sd * colMeans(y_hat_mat),
           median = X1%*%beta_hat + object$y_mean + object$y_sd * apply(y_hat_mat,2,'median'))
  }

  return(out)

} # end of predict function

########################################################################################################
# Predictions for classification
########################################################################################################

#' @method predict "cl_cspbart"
#' @rdname cl_cspbart
#' @param object x
#' @param newdata_x1 x
#' @param newdata_x2 x
#' @param type x
#' @param ... Catches unused arguments.
#' @usage
#' \method{predict}{cl_cspbart}(object,
#'         newdata_x1,
#'         newdata_x2,
#'         type = c('all', 'median', 'mean'),
#'         ...)
#' @export
predict.cl_cspbart = function(object, newdata_x1, newdata_x2,
                              type = c('all', 'median', 'mean'),
                              ...) {

  data_new = MakeDesignMatrixPredict(object$formula, newdata_x1)
  X1 = as.matrix(data_new$X) # matrix to be used in the linear predictor
  X2 = makeModelMatrixFromDataFrame(newdata_x2, drop = FALSE)
  beta_hat = colMeans(object$beta_hat)

  # Create holder for predicted values
  n_newX = dim(X2)[1]
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(X2))

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    X2,
                                    single_tree = length(curr_trees) == 1)
  }

    # Sort out what to return
    out = switch(type,
           all = sweep(y_hat_mat, 2, as.numeric(X1%*%beta_hat), "+"),
           mean = pnorm(X1%*%beta_hat + colMeans(y_hat_mat)),
           median = pnorm(X1%*%beta_hat + apply(y_hat_mat,2,'median')))

  return(out)

} # end of predict function
