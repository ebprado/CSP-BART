# -------------------------------------------------------------------------#
# Description: this script contains auxiliar functions needed to update    #
# the trees with details and to map the predicted values to each obs       #
# -------------------------------------------------------------------------#

# 1. fill_tree_details: takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in each terminal node
# 2. get_predictions: gets the predicted values from a current set of trees
# 3. get_children: it's a function that takes a node and, if the node is terminal, returns the node. If not, returns the children and calls the function again on the children
# 4. resample: an auxiliar function
# 5. get_ancestors: get the ancestors of all terminal nodes in a tree
# 6. update_s: full conditional of the vector of splitting probability.
# 7. get_number_distinct_cov: given a tree, it returns the number of distinct covariates used to create its structure
# 8. MakeDesignMatrix: it's a function that creates the design matrix for the linear based on the formula
# 9. var_used_trees: create a data frame with the frequency of the covariates used in trees.
# Fill_tree_details -------------------------------------------------------

fill_tree_details = function(curr_tree, X) {

  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix

  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix

  # Start with dummy node indices
  node_indices = rep(1, nrow(X))

  # For all but the top row, find the number of observations falling into each one
  for(i in 2:nrow(tree_matrix)) {

    # Get the parent
    curr_parent = as.numeric(tree_matrix[i,'parent'])

    # Find the split variable and value of the parent
    split_var = as.numeric(tree_matrix[curr_parent,'split_variable'])
    split_val = as.numeric(tree_matrix[curr_parent, 'split_value'])

    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')
    if(left_or_right == 'left') {
      # If left use less than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i
    } else {
      # If right use greater than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i
    }
  } # End of loop through table

  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))

} # End of function

# Get predictions ---------------------------------------------------------

get_predictions = function(trees, X, single_tree = FALSE) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]

  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] =
          trees$tree_matrix[unique_node_indices[i], 'mu']
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    # Do a recursive call to the function
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  +
      get_predictions(partial_trees, X,
                      single_tree = length(partial_trees) == 1)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }

  return(predictions)
}

# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}

# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]

update_s  <- function(var_count, p, alpha_s) {
  shape   <- alpha_s/p + var_count
  temp    <- rgamma(length(shape), shape, rate=1)
    temp/sum(temp)
}

get_number_distinct_cov <- function(tree){

  # Select the rows that correspond to internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 0)
  # Get the covariates that are used to define the splitting rules
  num_distinct_cov = length(unique(tree$tree_matrix[which_terminal,'split_variable']))

  return(num_distinct_cov)
}

sample_move = function(curr_tree, i, nburn){

  if (nrow(curr_tree$tree_matrix) == 1 || i < max(floor(0.1*nburn), 10)) {
    type = 'grow'
  } else {
    type = sample(c('grow', 'prune', 'change'), 1)
  }
  return(type)
}

## ---------------------------------------------------

ChangeColnames <- function(x) {
  colnames(x) <- paste0(".", colnames(x))
  x
}

MakeDesignMatrix <- function(formula, data){

  IsThereRandomEffects = try(silent = TRUE,
                             # When there is at least one random effect term in the formula
                             {parsedFormula = lFormula(formula = formula, data = data)
                             y_name = gsub('\\().*$', '', parsedFormula$formula[2]) # get the response variable name
                             y = data[,y_name]
                             Z = t(as.matrix(parsedFormula$reTrms$Zt)) # design matrix for random effects

                             # Number of random effect terms
                             number_random_effect_terms = length(parsedFormula$reTrms$cnms)

                             # When there are more than one random effect term
                             if (number_random_effect_terms > 1){
                               aux_indx_ini = 1
                               aux_indx_end = 0
                               for (i in 1:number_random_effect_terms){
                                 term = parsedFormula$reTrms$cnms[i]
                                 term_name = names(term)
                                 unique_values_cov = length(unique(data[,term_name]))
                                 aux_indx_end = aux_indx_ini + unique_values_cov * length(term[[1]]) - 1
                                 colnames(Z)[aux_indx_ini:aux_indx_end] = paste(colnames(Z)[aux_indx_ini:aux_indx_end], rep(parsedFormula$reTrms$cnms[[i]],unique_values_cov), sep='')
                                 colnames(Z)[aux_indx_ini:aux_indx_end] = gsub('\\.*\\(Intercept\\)',term_name , colnames(Z)[aux_indx_ini:aux_indx_end])
                                 aux_indx_ini = aux_indx_end + 1
                               }
                               X_Z = as.data.frame(as.matrix(cbind(parsedFormula$X, Z)))
                             } else {
                               colnames(Z) = paste(colnames(Z), parsedFormula$reTrms$cnms[[1]], sep='')
                               colnames(Z) = gsub('\\.*\\(Intercept\\)', names(parsedFormula$reTrms$cnms), colnames(Z))
                               X_Z = as.data.frame(as.matrix(cbind(parsedFormula$X, Z)))
                             }
                             return(list(y = y,
                                         X = X_Z))})

  # When there is no random effect terms (only fixed effects)
  if (is.na(IsThereRandomEffects[2])) {
    aux_termsFormula = formula(terms(formula, data=data))
    termsFormula = terms(aux_termsFormula)
    getIntercept = attr(termsFormula, 'intercept')
    getCovariates = attr(termsFormula, 'term.labels')

    backup_options <- options()
    options(contrasts = rep ("contr.sum", 2))
    X <- model.matrix(formula, data = data)
    options(backup_options)

    # I did this because it seems that the model.matrix function has a bug in the renaming when there isn't an intercept.
    # Specifically, model.matrix doesn't keep column names the same for factors - and that's why this HORRIBLE manipulation.
    vars_factor_char = lapply(data[, sapply(data, is.factor) | sapply(data, is.character)],
                              function(x) levels(x)
    )

    which_factors_have_different_name = which(getCovariates %in% names(vars_factor_char))
    var_different_name = getCovariates[which_factors_have_different_name]
    corrected_colnames = NULL

    if (length(var_different_name) > 0){

      for (i in 1:length(var_different_name)){
        var_different_name_levels = vars_factor_char[[var_different_name[i]]]
        corrected_colnames = append(corrected_colnames, paste(var_different_name[i], var_different_name_levels[-length(var_different_name_levels)], sep='')) # remove the reference level
      }

      which_colnames_to_change = which(attr(X, 'assign') %in% which_factors_have_different_name)
      colnames(X)[which_colnames_to_change] = corrected_colnames
    }

    # If I don't require intercept, model.matrix messes up wit the contrasts

    if (getIntercept == 1){
      X = X[,-1]
    } else {
      print('Add the intercept to the formula!')
      stop()
    }

    if (length(vars_factor_char) > 0){
      orig_names = NULL
      orig_names_dot = NULL

      for (i in 1:length(vars_factor_char)){
        sub_list_name = names(vars_factor_char[i])
        values = vars_factor_char[i][[1]]
        orig_names = append(orig_names, paste(sub_list_name, values, sep=''))
        orig_names_dot = append(orig_names_dot, paste(sub_list_name, values, sep='.'))
      }
      column_names = colnames(X)
      for (j in 1:ncol(X)){
        is_factor_var = which(orig_names %in% column_names[j])
        if(length(is_factor_var) > 0)
          colnames(X)[j] = orig_names_dot[orig_names %in% column_names[j]]
      }
    }

    y_name = gsub('\\().*$', '', formula[2]) # get the response variable name
    y = data[,y_name]
    return(list(y = y,
                X = X))
  }
}

MakeDesignMatrixPredict <- function(formula, data){

  aux1 = gsub("[^[:alnum:]_.]", ' ', formula[3]) # take the terms in the predictor and remove numbers and special characters (expect _ and .)
  aux2 = sort(unlist(strsplit(aux1, ' ')), decreasing = TRUE)[1] # get the first non-empty character
  new_formula = as.formula(paste(aux2, '~' ,formula[3])) # remove the response, as the new dataset don't have it

  IsThereRandomEffects = try(silent = TRUE,
    # When there is at least one random effect term in the formula
    {parsedFormula = lFormula(formula = new_formula, data = data)
    Z = t(as.matrix(parsedFormula$reTrms$Zt)) # design matrix for random effects

    # Number of random effect terms
    number_random_effect_terms = length(parsedFormula$reTrms$cnms)

    # When there are more than one random effect term
    if (number_random_effect_terms > 1){
      aux_indx_ini = 1
      aux_indx_end = 0
      for (i in 1:number_random_effect_terms){
        term = parsedFormula$reTrms$cnms[i]
        term_name = names(term)
        unique_values_cov = length(unique(data[,term_name]))
        aux_indx_end = aux_indx_ini + unique_values_cov * length(term[[1]]) - 1
        colnames(Z)[aux_indx_ini:aux_indx_end] = paste(colnames(Z)[aux_indx_ini:aux_indx_end], rep(parsedFormula$reTrms$cnms[[i]],unique_values_cov), sep='')
        colnames(Z)[aux_indx_ini:aux_indx_end] = gsub('\\.*\\(Intercept\\)',term_name , colnames(Z)[aux_indx_ini:aux_indx_end])
        aux_indx_ini = aux_indx_end + 1
      }
      X_Z = as.data.frame(as.matrix(cbind(parsedFormula$X, Z)))
    } else {
      colnames(Z) = paste(colnames(Z), parsedFormula$reTrms$cnms[[1]], sep='')
      colnames(Z) = gsub('\\.*\\(Intercept\\)', names(parsedFormula$reTrms$cnms), colnames(Z))
      X_Z = as.data.frame(as.matrix(cbind(parsedFormula$X, Z)))
    }
    return(list(X = X_Z))})

  # When there is no random effect terms (only fixed effects)
  if (is.na(IsThereRandomEffects[2])) {
    aux_termsFormula = formula(terms(formula, data=data))
    termsFormula = terms(aux_termsFormula)
    getIntercept = attr(termsFormula, 'intercept')
    getCovariates = attr(termsFormula, 'term.labels')

    backup_options <- options()
    options(contrasts = rep ("contr.sum", 2))
    X <- model.matrix(formula, data = data)
    options(backup_options)

    # I did this because it seems that the model.matrix function has a bug in the renaming when there isn't an intercept.
    # Specifically, model.matrix doesn't keep column names the same for factors - and that's why this HORRIBLE manipulation.
    vars_factor_char = lapply(data[, sapply(data, is.factor) | sapply(data, is.character)],
                              function(x) levels(x)
    )

    which_factors_have_different_name = which(getCovariates %in% names(vars_factor_char))
    var_different_name = getCovariates[which_factors_have_different_name]
    corrected_colnames = NULL

    if (length(var_different_name) > 0){

      for (i in 1:length(var_different_name)){
        var_different_name_levels = vars_factor_char[[var_different_name[i]]]
        corrected_colnames = append(corrected_colnames, paste(var_different_name[i], var_different_name_levels[-length(var_different_name_levels)], sep='')) # remove the reference level
      }

      which_colnames_to_change = which(attr(X, 'assign') %in% which_factors_have_different_name)
      colnames(X)[which_colnames_to_change] = corrected_colnames
    }

    # If I don't require intercept, model.matrix messes up wit the contrasts

    if (getIntercept == 1){
      X = X[,-1]
    } else {
      print('Add the intercept to the formula!')
      stop()
    }

    # I did this because it seems that the model.matrix function has a bug in the renaming when there isn't an intercept
    vars_factor_char = lapply(data[, sapply(data, is.factor) | sapply(data, is.character)],
                              function(x) levels(x)
    )

    if (length(vars_factor_char) > 0){
      orig_names = NULL
      orig_names_dot = NULL

      for (i in 1:length(vars_factor_char)){
        sub_list_name = names(vars_factor_char[i])
        values = vars_factor_char[i][[1]]
        orig_names = append(orig_names, paste(sub_list_name, values, sep=''))
        orig_names_dot = append(orig_names_dot, paste(sub_list_name, values, sep='.'))
      }
      column_names = colnames(X)
      for (j in 1:ncol(X)){
        colnames(X)[j] = orig_names_dot[orig_names %in% column_names[j]]
      }
    }
    return(list(X = X))
  }
}

#' Get covariates used as splits per iteration
#' @export
#' @param object x
#' @param raw x
#'
#' @return x
#' @examples
#'
var_used_trees = function(object, raw = FALSE) {

  # Create holder for predicted values
  n_its = object$npost
  ntrees = object$ntrees
  colnames_x2 = object$colnames.x2

  # Get which covariates are used by each tree in each MCMC iteration
  vars_trees = matrix(NA, nrow=n_its, ncol=ntrees)
  for (i in seq_len(n_its)) {
    for (j in seq_len(ntrees)) {
      new_tree = object$trees[[i]][[j]]$tree_matrix
      internal_nodes = as.numeric(which(new_tree[,'terminal'] == 0))
      aux = new_tree[internal_nodes, 'split_variable']
      if(length(aux) > 1){
        vars_trees[i,j] = paste(colnames_x2[unique(sort(aux[!is.na(aux)]))], collapse = ',')
      }
    }
  }
if(raw == TRUE) {return(data.frame(vars_trees))}
if(raw == FALSE) {
  vars_trees = as.data.frame(table(vars_trees))
  aux_count_var = strsplit(as.character(vars_trees[,'vars_trees']),',')
  vars_trees$count = sapply(aux_count_var, function(x) length(x))
  return(vars_trees)
  }
}

get_ancestors = function(tree){

  save_ancestor = NULL
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  if(nrow(tree$tree_matrix) == 1) {
    save_ancestor = cbind(terminal = NULL,
                          ancestor = NULL)
  } else {
    for (k in seq_len(length(which_terminal))) {
      get_parent = as.numeric(as.character(tree$tree_matrix[which_terminal[k], 'parent'])) # get the 1st parent
      get_split_var = as.character(tree$tree_matrix[get_parent, 'split_variable']) # then, get the covariate associated to the row of the parent

      save_ancestor = rbind(save_ancestor,
                            cbind(terminal = which_terminal[k],
                                  # parent   = get_parent,
                                  ancestor = get_split_var))
      while (get_parent > 1){
        get_parent = as.numeric(as.character(tree$tree_matrix[get_parent,'parent'])) # then, get the subsequent parent
        get_split_var = as.character(tree$tree_matrix[get_parent, 'split_variable']) # then, get the covariate associated to the row of the new parent
        save_ancestor = rbind(save_ancestor,
                              cbind(terminal = which_terminal[k],
                                    # parent   = get_parent,
                                    ancestor = get_split_var))
      }
    }
    save_ancestor = unique(save_ancestor) # remove duplicates
    save_ancestor = save_ancestor[order(save_ancestor[,1], save_ancestor[,2]),] # sort by terminal and ancestor
  }

  return(save_ancestor)
}

get_ancestors_internal = function(tree){
  save_ancestor = NULL
  tree = tree$tree_matrix
  which_internal = which(tree[,'terminal'] == 0)


  if(nrow(tree) == 1) {
    save_ancestor = cbind(internal = NULL,
                          ancestor = NULL)
  } else {
    for (k in length(which_internal):1) {
      internal_node = which_internal[k]
      parent = tree[internal_node, 'parent']
      get_split_var = tree[internal_node, 'split_variable']

      save_ancestor = rbind(save_ancestor,
                            cbind(internal = internal_node,
                                  parent   = parent,
                                  split_var = get_split_var))
      while (is.na(parent) == FALSE && parent > 0) {
        get_split_var = tree[parent, 'split_variable']
        parent = tree[parent, 'parent']
        save_ancestor = rbind(save_ancestor,
                              cbind(internal = internal_node,
                                    parent   = parent,
                                    split_var = get_split_var))
      }
    }
  }
  return(save_ancestor[,,drop=FALSE])

}
