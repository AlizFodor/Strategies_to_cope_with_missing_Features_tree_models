# FUNCTIONS SCRIPT

# The following script contains all functions necessary for the workflow. 


# create a stratified test/train partition
split_data <- function(X, y, train_frac = 0.7, seed = 455) {
  set.seed(seed)
  train_i <- caret::createDataPartition(y, p = train_frac, list = FALSE)
  list(
    X_train = X[train_i, , drop = FALSE],
    y_train = y[train_i],
    X_test  = X[-train_i, , drop = FALSE],
    y_test  = y[-train_i]
  )
}


# compute pre-determined evaluation metrics by comparing true response to the 
# estimated predictions
compute_metrics_from_preds <- function(y_true, prob, threshold) {
  
  roc_obj <- suppressWarnings(pROC::roc(y_true, prob))
  
  auc <- as.numeric(roc_obj$auc)
  pred_class <- ifelse(prob > threshold, 1, 0)
  
  acc <- mean(pred_class == y_true)
  tp <- sum(pred_class == 1 & y_true == 1)
  fp <- sum(pred_class == 1 & y_true == 0)
  fn <- sum(pred_class == 0 & y_true == 1)
  tn <- sum(pred_class == 0 & y_true == 0)
  
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  spec <- ifelse(tn + fp == 0, 0, tn / (tn + fp))
  
  bal_acc <- (recall + spec) / 2 # class balanced accurcay 
  
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  list(accuracy = acc, balanced_accuracy = bal_acc, AUC = auc, precision = precision,
       recall = recall, sensitivity = recall, specificity = spec, F1 = f1)
}



# return a dataset with certain amount of removed features, two options of masking features:
# zero masking: for logistic regr and zero sum regr set all values of one features to zero
# NA masking: for tree based models set values of one feature to NA 



### --- LOGISTIC REGRESSION

# fit a generalized linear model (logistic regression)
train_logistic <- function(X, y, weights=NULL) {
  D <- data.frame(y = y, X)
  
  model <- glm(y ~ ., 
               data = D, 
               family = binomial(link = "logit"),
               weights = weights) # weigths to be used in fitting process
  return(model)
}

# run fitted logistic regression model on an unseen testing dataset
predict_logistic <- function(model, X, y_true) {
  Xdf <- as.data.frame(X)
  prob <- predict(model, newdata = Xdf, type = "response")
  
  return(list(probabilities = prob, y_test = y_true))
}



### --- ZERO SUM REGRESSION 

# under partial feature availability during prediction, zero sum regression model is retrained on the available features
# to preserve the zero sum constraint, same amount of negative and positive coefficients are kept for the retraining process.

# this function determines which coefficients to keep for refitting a model (based on the zero sum constraint, sum of coefficients needs to be zero)
betas_keep <- function(zs_beta, missing_cols) {
  
  # split into pos and neg coefficients
  beta_pos <- zs_beta[zs_beta > 0]
  beta_neg <- zs_beta[zs_beta < 0]
  
  # sort by magnitude
  beta_pos <- sort(beta_pos, decreasing=TRUE)
  beta_neg <- sort(beta_neg, decreasing=FALSE)
  
  # how many features missing per group (neg, pos)
  m_pos <- mean(names(beta_pos) %in% missing_cols)
  m_neg <- mean(names(beta_neg) %in% missing_cols)
  
  # How many features to keep 
  k_pos <- floor((1 - m_pos) * length(beta_pos))
  k_neg <- floor((1 - m_neg) * length(beta_neg))
  
  # keep at least one from each sign
  k_pos <- max(1, k_pos)
  k_neg <- max(1, k_neg)
  
  # which features to keep 
  keep_pos <- names(beta_pos)[1:k_pos]
  keep_neg <- names(beta_neg)[1:k_neg]
  
  return(unique(c(keep_pos, keep_neg)))
}


### --- RANDOM FOREST 

# Train a random forest (ranger package)
train_rf <- function(X, y, num_trees = 150, mtry = NULL, max_depth = NULL, weights=NULL) {
  
  if (is.null(mtry)) mtry <- floor(sqrt(ncol(X))) # default parameter in CART
  dat <- data.frame(y = as.factor(y), X)
  
  model <- ranger(
    y ~ .,
    data = dat,
    case.weights = weights,
    num.trees = num_trees,
    mtry = mtry,
    probability = TRUE,
    importance = "impurity",
    max.depth = max_depth
  )
  return(model)
}


# run random forest on an unseen test data
predict_rf <- function(model, X, y_test = NULL) {
  pred <- predict(model, data = data.frame(X))
  prob <- pred$predictions[, "1"]# prob for class 1 
  
  return(list(probabilities = prob, y_test = y_test))
}


### --- RANDOM FOREST WITH SURROGATE SPLITS

# the following code is designed to extend an existing random forest model with internal surrogate
# splitting nodes. The initial input is a random forest object (of the form of the ranger package). 
# The resulting hybrid_surrogate_forest is an extension of this object that contains alternative splitting rules
# for each internal node. Additionally, the output object contains options of majority routing (direction of largest impurity reduction). 
# The functions are implemented according to the logic proposed in CART. 


# ===== 1.) Best surrogate split search
# orienation of the routings: 
# - 0: normal orientation (L -> Lower, U -> Upper)
# - 1: flipped orientation (L -> Upper, U -> Lower)

best_surrogate_split <- function(x, mainrouting, majority_score) {
  # x is vector of values of the surrogate feature
  # mainrouting is the routing of those samples encoded in {-1, 1}
  # majority score: agreement by always routing samples to majority side of primary split (used as baseline for comparison)
  
  no_na <- !is.na(x) & !is.na(mainrouting)
  x <- x[no_na]
  y <- mainrouting[no_na]
  m <- length(x)
  if (m < 2) return(NULL) # at least two samples present
  
  ord <- order(x)
  x <- x[ord] # ordered by increasing surrogate feature values
  y <- y[ord] # rearranged according to surrogate vector
  
  if (x[1] == x[m]) return(NULL)  # no valid split possible
  
  TU <- sum(y == 1L) # number of true upper in primary split
  TL <- m - TU # number of true lower
  
  # S[t] = S[t-1] + y_sorted[t], with S0 = 0, sum of +/- labels in y until t
  S <- c(0L, cumsum(as.integer(y)))
  
  best_score <- -Inf
  best_t <- NA_integer_
  best_orient <- NA_integer_
  
  for (t in 1:(m - 1)) {
    if (x[t] == x[t + 1]) next
    
    # S index shift because of first 0
    correct_norm <- TU - S[t + 1L]  # orient 0
    correct_flip <- TL + S[t + 1L]  # orient 1
    
    if (correct_norm > best_score) {
      best_score <- correct_norm
      best_t <- t
      best_orient <- 0L
    }
    if (correct_flip > best_score) {
      best_score <- correct_flip
      best_t <- t
      best_orient <- 1L
    }
  }
  
  if (is.na(best_t)) return(NULL)
  
  agree <- best_score / m
  if (agree <= majority_score) return(NULL)
  
  thr <- (x[best_t] + x[best_t + 1L]) / 2
  list(threshold = thr, agree = agree, orient = best_orient)
}


# ==== 2.) Compute Surrogates for one node:
# - for one node store all surrogate information in a matrix or NULL
# - Each row one surrogate routing
# - cols: surr_index in X, split threshold, agreement, orient
# - top K surrogates kept (matrix sorted by decreasing)
# - Additionally only consider surrogates with abs(Pearson correlation > 0.2)

compute_surrogates_for_node_fast <- function(X_node, primary_feature_id, primary_split,
                                             rho_min = 0.2, topK = 5) {
  # X_node: subset of rows (samples) reaching the node with same cols as X
  # primary_feature_id: integer indexing primary feature in X
  # primary_split: threshold
  # rho_min: threshold on |p| (absolute pearson correlation) to determine which surroagtes to consider
  # top K: Max number of surrogates kept
  
  X_node <- as.matrix(X_node)
  m <- nrow(X_node)
  if (m < 2) return(NULL)
  
  x_primary <- X_node[, primary_feature_id]
  # primary routing labels: -1 = Left/Lower, +1 = Right/Upper
  mainlabel <- ifelse(x_primary <= primary_split, -1L, +1L)
  
  TU <- sum(mainlabel == 1L)
  TL <- m - TU
  baseline <- max(TU, TL) / m
  
  p <- ncol(X_node)
  cand <- setdiff(seq_len(p), primary_feature_id)
  if (length(cand) == 0) return(NULL)
  
  # Pearson correlation prefilter (vectorized)
  rho <- suppressWarnings(cor(X_node[, cand, drop = FALSE], x_primary, use = "complete.obs"))
  keep <- cand[!is.na(rho) & (abs(rho) >= rho_min)]
  q <- length(keep)
  if (q == 0) return(NULL)
  
  out_mat <- matrix(NA_real_, nrow = q, ncol = 4)
  colnames(out_mat) <- c("feature", "split", "agree", "orient")
  k_out <- 0L
  
  for (k in keep) {
    xk <- X_node[, k]
    if (all(is.na(xk))) next
    if (length(unique(xk[!is.na(xk)])) < 2) next
    
    out <- best_surrogate_split(xk, mainlabel, baseline)
    if (is.null(out)) next
    
    k_out <- k_out + 1L
    out_mat[k_out, 1] <- k
    out_mat[k_out, 2] <- out$threshold
    out_mat[k_out, 3] <- out$agree
    out_mat[k_out, 4] <- out$orient
  }
  
  if (k_out == 0L) return(NULL)
  
  out_mat <- out_mat[seq_len(k_out), , drop = FALSE]
  # sort decreasing agree
  out_mat <- out_mat[order(-out_mat[, "agree"]), , drop = FALSE]
  
  if (!is.null(topK) && nrow(out_mat) > topK) {
    out_mat <- out_mat[seq_len(topK), , drop = FALSE]
  }
  
  out_mat[, "feature"] <- as.integer(out_mat[, "feature"])
  out_mat[, "orient"]  <- as.integer(out_mat[, "orient"])
  
  out_mat
}

# ==== 3) Compute node membership for samples:
# - Given a tree, computes set of sample indices that reach each node
# - top down approach

build_node_samples_topdown <- function(tree_df, X_train) {
  node_ids <- tree_df$nodeID
  node_samples <- vector("list", length(node_ids))
  names(node_samples) <- as.character(node_ids)
  
  # all samples in the root:
  root <- 0L
  node_samples[[as.character(root)]] <- seq_len(nrow(X_train))
  
  # process nodes by increasing depth
  ord_nodes <- seq_len(nrow(tree_df))
  
  for (rr in ord_nodes) {
    row <- tree_df[rr, ]
    nid <- row$nodeID
    idx <- node_samples[[as.character(nid)]]
    if (length(idx) == 0L) next
    if (isTRUE(row$terminal)) next
    
    if (is.na(row$splitvarID)) next
    j <- row$splitvarID + 1L
    s <- row$splitval
    
    left  <- row$leftChild
    right <- row$rightChild
    
    xj <- X_train[idx, j]
    left_idx  <- idx[xj <= s]
    right_idx <- idx[xj >  s]
    
    if (!is.na(left) && length(left_idx)) {
      keyL <- as.character(left)
      node_samples[[keyL]] <- c(node_samples[[keyL]], left_idx)
    }
    if (!is.na(right) && length(right_idx)) {
      keyR <- as.character(right)
      node_samples[[keyR]] <- c(node_samples[[keyR]], right_idx)
    }
  }
  node_samples
}


# ==== 4) Reconstruct random forest with additional stored surrogate routing rules
# - each tree stored as list of nodes 
# - each node stores surrogate routing under missingness + primary split as matrices
# - this is the main function of building the forest

build_hybrid_forest_from_ranger <- function(rf_model, X_train, y_train,
                                            rho_min = 0.2, topK = 5) {
  X_train <- as.matrix(X_train)
  y_train <- as.integer(y_train)
  
  n_trees <- rf_model$num.trees
  forest <- vector("list", n_trees)
  
  for (t in seq_len(n_trees)) {
    # extract the tree inforatmion from the ranger obj stored as lists: 
    tree_df <- ranger::treeInfo(rf_model, tree = t)
    
    node_samples <- build_node_samples_topdown(tree_df, X_train)
    
    nodes <- vector("list", nrow(tree_df))
    node_ids <- tree_df$nodeID
    names(nodes) <- as.character(node_ids)
    
    for (r in seq_len(nrow(tree_df))) {
      row <- tree_df[r, ]
      nid <- row$nodeID
      
      splitvarID <- if (!is.na(row$splitvarID)) row$splitvarID + 1L else NA_integer_
      splitval   <- row$splitval
      left       <- row$leftChild
      right      <- row$rightChild
      is_terminal <- isTRUE(row$terminal)
      
      idx <- node_samples[[as.character(nid)]]
      
      node <- list(
        nodeID = nid,
        splitvarID = splitvarID,
        splitval = splitval,
        leftChild = left,
        rightChild = right,
        isTerminal = is_terminal,
        surrogates = NULL,
        classProb = if(is_terminal) row$pred.1 else NA_real_
      )
      
      if (!is_terminal && length(idx) > 0) {
        node$classProb <- mean(y_train[idx] == 1L)
      }
      
      if (!is_terminal) {
        n_left  <- if (!is.na(left))  length(node_samples[[as.character(left)]])  else 0L
        n_right <- if (!is.na(right)) length(node_samples[[as.character(right)]]) else 0L
        
        if (!is.na(left) && !is.na(right)) {
          node$majorityChild <- if (n_left >= n_right) left else right
        } else if (!is.na(left)) {
          node$majorityChild <- left
        } else if (!is.na(right)) {
          node$majorityChild <- right
        }
      }
      
      if (!is_terminal && !is.na(splitvarID) && length(idx) > 1) {
        X_node <- X_train[idx, , drop = FALSE]
        node$surrogates <- compute_surrogates_for_node_fast(
          X_node = X_node,
          primary_feature_id = splitvarID,
          primary_split = splitval,
          rho_min = rho_min,
          topK = topK
        )
      }
      nodes[[as.character(nid)]] <- node
    }
    forest[[t]] <- list(nodes = nodes, rootID = 0L)
  }
  list(forest = forest, n_trees = n_trees)
}

#


# 5) Prediction with itnegrated surrogate routing rules
# - orient 0=normal:(x <= split)-> left, else right
# - orient 1 = flipped:(x<= split)->right, else left
# to fallbacks possible if no surrogate found: majority routing or tree truncation

predict_tree_with_surrogates <- function(tree, x_row, fallback = c("majority", "prob")) {
  fallback <- match.arg(fallback)
  nodeID <- tree$rootID
  
  route_used <- "primary"

  repeat {
    node <- tree$nodes[[as.character(nodeID)]]
    
    if (isTRUE(node$isTerminal)) {return(list(prob = node$classProb, route = route_used))}
    
    j <- node$splitvarID
    s <- node$splitval
    v <- x_row[j]
    
    if (!is.na(v)) {
      nodeID <- if (v <= s) node$leftChild else node$rightChild
      next
    }
    
    # primary missing -> try surrogates
    used <- FALSE
    sur <- node$surrogates
    used_surrogate <- FALSE

    if (!is.null(sur) && nrow(sur) > 0) {
      for (r in seq_len(nrow(sur))) {
        f_id <- sur[r, "feature"]
        sp   <- sur[r, "split"]
        ori  <- sur[r, "orient"]
        vs   <- x_row[f_id]
        
        if (!is.na(vs)) {
          if (ori == 0L) {
            nodeID <- if (vs <= sp) node$leftChild else node$rightChild
          } else {
            nodeID <- if (vs <= sp) node$rightChild else node$leftChild
          }
          route_used <- "surrogate"
          used_surrogate <- TRUE
          break
        }
      }
    }
    
    if (used_surrogate) next
    
    # fallback
    if (fallback == "majority" && !is.na(node$majorityChild)) {
      nodeID <- node$majorityChild
      route_used <- "majority"
      next
    } else {
      return(list(prob = node$classProb, route = "prob"))
    }
  }
}


predict_hybrid_forest<- function(hybrid_forest, X_new, fallback = c("majority", "prob")) {
  fallback <- match.arg(fallback)
  X_new <- as.matrix(X_new)
  
  n <- nrow(X_new)
  T_total <- hybrid_forest$n_trees
  prob_mat <- matrix(NA_real_, nrow = n, ncol = T_total)
  
  route_counts <- c(primary=0L, surrogate=0L, majority=0L, prob=0L)
  
  for (t in seq_len(T_total)) {
    tree <- hybrid_forest$forest[[t]]
    for (i in seq_len(n)) {
      res <- predict_tree_with_surrogates(tree, X_new[i, ], fallback = fallback)
      prob_mat[i, t] <-res$prob
      route_counts[res$route] <- route_counts[res$route] +1L
    }
  }
  list(probabilities = rowMeans(prob_mat, na.rm = TRUE),
       routing_counts = route_counts)
}



### --- RANDOM FOREST WITH TREE DROPPING

# drop the entire tree in a random forest model that contains a missing feature
# return the prediction using only the remaining forsests

# helper for the following function:
predict_tree_drop <- function(tree, x_row) {
  
  nodeID <- tree$rootID
  
  repeat {
    node <- tree$nodes[[as.character(nodeID)]]
    
    if (isTRUE(node$isTerminal)) {return(list(prob = node$classProb, used = TRUE))}
    
    j <- node$splitvarID
    s <- node$splitval
    v <- x_row[j]
    
    if (is.na(v)) {return(list(prob = NA_real_, used = FALSE))}
    
    nodeID <- ifelse(v <= s, node$leftChild, node$rightChild)
  }
}

# input: the hybrid forest and the reduced dataset
predict_rf_drop_trees <- function(hybrid_forest, X_new) {
  
  X_new <- as.matrix(X_new)
  n <- nrow(X_new)
  T_total <- hybrid_forest$n_trees
  
  prob_mat <- matrix(NA_real_, n, T_total)
  used_mat <- matrix(FALSE, n, T_total)
  
  for (t in seq_len(T_total)) {
    tree <- hybrid_forest$forest[[t]]
    
    for (i in seq_len(n)) {
      res <- predict_tree_drop(tree, X_new[i, ])
      prob_mat[i, t] <- res$prob
      used_mat[i, t] <- res$used
    }
  }
  list(probabilities = rowMeans(prob_mat, na.rm = TRUE),n_used_trees  = rowSums(used_mat))
}


### --- TREE TRUNCTATION

# Truncate trees with missing features, return probabilites at nodes

predict_tree_early_stop <- function(tree, x_row) {
  
  nodeID <- tree$rootID
  
  repeat {
    node <- tree$nodes[[as.character(nodeID)]]
    
    if (isTRUE(node$isTerminal)) {return(node$classProb)}
    
    j <- node$splitvarID
    s <- node$splitval
    v <- x_row[j]
    
    # return the probability:
    if (is.na(v)) {
      return(node$classProb)
    }
    
    nodeID <- ifelse(v <= s, node$leftChild, node$rightChild)
  }
}

# takes as input the random forest with surrogates
predict_rf_early_stop <- function(hybrid_forest, X_new) {
  X_new <- as.matrix(X_new)
  n <- nrow(X_new)
  T_total <- hybrid_forest$n_trees
  
  prob_mat <- matrix(NA_real_, n, T_total)
  
  for (t in seq_len(T_total)) {
    tree <- hybrid_forest$forest[[t]]
    for (i in seq_len(n)) {prob_mat[i, t] <- predict_tree_early_stop(tree, X_new[i, ])}
  }
  rowMeans(prob_mat, na.rm = TRUE)
}


### --- XG BOOST 

# train one XG Boost model, tune parameters with grid search over cv runs
train_xgb <- function(X, y,
                      folds = 5,
                      nrounds_max = 1000, 
                      early_stopping_rounds=50,
                      seed=468) {
  set.seed(seed)
  
  n0 <- sum(y == 0)
  n1 <- sum(y == 1)
  
  scale_pos_weight <- n0 / n1
  
  dtrain <- xgb.DMatrix(data=as.matrix(X), label=y, missing=NA)
  
  # parameter tuning
  grid <- expand.grid(
    eta=c(0.01,0.05,0.1,0.15,0.02),
    max_depth=c(2,4,6),
    min_child_weight=c(1,5),
    subsample=c(0.7,0.9),
    colsample_bytree=c(0.7,0.9),
    gamma=c(0,1)
  )
  
  best_auc <- -Inf
  best_iter <- NA_integer_
  best_params <- NULL

  for (i in seq_len(nrow(grid))){
    g <- grid[i,]
    
    params <- list(
      booster="gbtree",
      objective="binary:logistic",
      eval_metric="auc",
      eta=g$eta,
      max_depth=g$max_depth,
      min_child_weight=g$min_child_weight,
      subsample=g$subsample,
      colsample_bytree=g$colsample_bytree,
      gamma=g$gamma,
      scale_pos_weight=scale_pos_weight
    )
    cv<-xgb.cv(
      params=params,
      data=dtrain,
      nrounds=nrounds_max,
      nfold=folds,
      stratified=TRUE,
      early_stopping_rounds = early_stopping_rounds,
      maximize = TRUE,
      verbose=0
    )
    auc_mean <- max(cv$evaluation_log$test_auc_mean)

    if (auc_mean > best_auc){
      best_auc <- auc_mean
      best_iter <- cv$evaluation_log$iter[which.max(cv$evaluation_log$test_auc_mean)]
      best_params <- params
    }
    
  }
  model<-xgb.train(params=best_params,data=dtrain,nrounds=best_iter,verbose=0)
  
  return(list(model=model, 
              best_params=best_params,
              best_iteration=best_iter, 
              cv_test_auc=best_auc))
  }


# run xg boost model on new unseen data 
predict_xgb <- function(model, X, y_test = NULL) {
  dtest <- xgb.DMatrix(data = as.matrix(X), missing=NA)
  prob  <- predict(model, dtest)
  return(list(probabilities = prob, y_test = y_test))
}





