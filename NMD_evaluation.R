# Evaluation on the NMD datasets
# Script for evaluation of all models on the NMD dataset

# The following script contains following evaluation steps: 
# 1. Training of each model on the dataset 
# 2. Visualization of results
# 3. Gradual removal of features of dataset
# 4. Running each model on the partial dataset 


library(dplyr)
library(tidyr)
library(ggplot2)
library(zeroSum)  
library(caret)
library(pROC)
library(ranger)
library(rpart)
library(patchwork)
library(xgboost)
library(tibble)
library(stats)


#--------------------------------------------------------------------------------
# DATASETS 
#--------------------------------------------------------------------------------

smartSeq   <- as.data.frame(readRDS("Datasets/SMARTSEQ_FINAL_128genes_noscale.rds"))
perturbSeq <- as.data.frame(readRDS("Datasets/PERTURB_FINAL_128genes_noscale.rds"))

sparsity <- function(X){
  mean(X == 0)
}

sparsity(smartSeq)
sparsity(perturbSeq)


colnames(smartSeq)      <- make.names(colnames(smartSeq))
colnames(perturbSeq)    <- make.names(colnames(perturbSeq))

# load response vectors (binary vectors, 0=NMD_inactive, 1=NMD_active)
y_smart    <- readRDS("./Datasets/y_smart.rds")
y_pert     <- readRDS("./Datasets/y_pert.rds")


# ---- Downsample the Perturbseq Dataset 

set.seed(093)
table(y_pert)

i_0 <- which(y_pert==0)
i_1 <- sample(which(y_pert==1), length(i_0)*2)

keep <- c(i_0, i_1)

perturbSeq <- perturbSeq[keep,]
y_pert     <- y_pert[keep]

table(y_pert)
table(y_smart)

roc_runs <- vector("list", B)


# FUNCTIONS: See script "Functions"

##### --- ITERATION --- #####
# In the following part all models are trained over B iterations on the complete data (no feature removal)

# number of iterations: 
B <- 1

set.seed(1)
seeds <- sample.int(1e6, B)

all_runs <- vector("list", B)

for (b in seq_len(B)) { 
  
  set.seed(seeds[b])
  
  ##### ---- Z-scaling 
  
  # this is coloumn wise standardization: 
  normalize <- function(x, mean, sd) { 
    return((x - mean) / sd) 
  }
  
  # Scaling of SmartSeq:
  smart_mean <- colMeans(smartSeq)
  smart_sd   <- apply(smartSeq, 2, sd)
  smart_scaled   <- normalize(smartSeq, smart_mean, smart_sd)
  
  sp <- split_data(smart_scaled, y_smart, train_frac = 0.7, seed = seeds[b])
  
  X_train_scaled <- sp$X_train
  X_test_scaled <- sp$X_test
  
  # Scaling of PerturbSeq:
  pert_mean <- colMeans(perturbSeq)
  pert_sd   <- apply(perturbSeq, 2, sd)
  pert_scaled <- normalize(perturbSeq, pert_mean, pert_sd)
  
  
  #### ---- Class weights
  
  n0 <- sum(sp$y_train == 0)
  n1 <- sum(sp$y_train == 1)
  
  w0 <- 1 / n0
  w1 <- 1 / n1
  
  class_weights <- ifelse(sp$y_train == 1, w1, w0)
  
  #--------------------------------------------------------------------------------
  # MODEL TRAINING
  #--------------------------------------------------------------------------------
  
  ### --- LOGISTIC REGRESSION
  
  logist_model <- train_logistic(X_train_scaled, sp$y_train, weights=class_weights)
  
  pred_logi_smart <- predict_logistic(logist_model, X_test_scaled, sp$y_test)
  pred_logi_pert  <- predict_logistic(logist_model, pert_scaled, y_pert)
  
  metr_logi_smart <- compute_metrics_from_preds(y_true=pred_logi_smart$y_test, prob = pred_logi_smart$probabilities, threshold=0.5)
  metr_logi_pert  <- compute_metrics_from_preds(y_true=pred_logi_pert$y_test, prob = pred_logi_pert$probabilities, threshold=0.5)
  
  logi_results <- bind_rows(tibble::tibble(dataset = "SmartSeq", model   = "Logistic regression", !!!metr_logi_smart),
                            tibble::tibble(dataset = "PerturbSeq",model   = "Logistic regression",!!!metr_logi_pert))
  
  ### --- ZERO SUM REGRESSION 
  
  fit_zerosum <- zeroSum::zeroSum(x = as.matrix(X_train_scaled), y = sp$y_train, family  = "binomial", weights = class_weights, alpha   = 0)
  
  prob_zs_smart <- zeroSum:::predict.zeroSum(fit_zerosum, newx = as.matrix(X_test_scaled), type = "response")
  prob_zs_pert <- zeroSum:::predict.zeroSum(fit_zerosum, newx = as.matrix(pert_scaled), type = "response")
  
  metr_zs_smart <- compute_metrics_from_preds(y_true=sp$y_test,prob=prob_zs_smart,threshold=0.5)
  metr_zs_pert <- compute_metrics_from_preds(y_true=y_pert,prob=prob_zs_pert,threshold=0.5)
  
  zs_results <- bind_rows(
    tibble::tibble(dataset= "SmartSeq",model="ZeroSum",!!!metr_zs_smart),
    tibble::tibble(dataset= "PerturbSeq",model= "ZeroSum",!!!metr_zs_pert))
  
  # save to reduce time: 
  #saveRDS(fit_zerosum, "fitzerosum.rds")
  #fit_zerosum <- readRDS("fitzerosum.rds")
  
  
  ### --- RANDOM FOREST 
  
  rf_model <- train_rf(X_train_scaled, sp$y_train, num_trees = 100, max_depth = NULL, weights=class_weights)
  
  pred_rf_smart <- predict_rf(rf_model, X_test_scaled, y_test = sp$y_test)
  pred_rf_pert  <- predict_rf(rf_model, pert_scaled, y_test = y_pert)
  
  metr_rf_smart <- compute_metrics_from_preds(y_true=pred_rf_smart$y_test,prob= pred_rf_smart$probabilities,threshold = 0.5)
  metr_rf_pert <- compute_metrics_from_preds(y_true=pred_rf_pert$y_test,prob= pred_rf_pert$probabilities,threshold = 0.5)
  
  rf_results <- bind_rows(
    tibble::tibble(dataset="SmartSeq", model="RandomForest", !!!metr_rf_smart),
    tibble::tibble(dataset = "PerturbSeq", model="RandomForest", !!!metr_rf_pert))
  
  ### --- RANDOM FOREST WITH SURROGATE ROUTINGS (=HYBRID FOREST)
  
  hybrid_forest_train <- build_hybrid_forest_from_ranger(
    rf_model = rf_model,
    X_train  = X_train_scaled,
    y_train  = sp$y_train
  )
  
  pred_hybrid_smart <- predict_hybrid_forest(hybrid_forest_train, X_test_scaled)
  pred_hybrid_pert <- predict_hybrid_forest(hybrid_forest_train, pert_scaled)
  
  metr_hybrid_smart <- compute_metrics_from_preds(y_true= sp$y_test,prob= pred_hybrid_smart$probabilities,threshold = 0.5)
  metr_hybrid_pert <- compute_metrics_from_preds(y_true= y_pert,prob= pred_hybrid_pert$probabilities,threshold = 0.5)
  
  hybrd_results <- bind_rows(
    tibble::tibble(dataset = "SmartSeq", model = "SurrogateHybrid", !!!metr_hybrid_smart),
    tibble::tibble(dataset = "PerturbSeq", model = "SurrogateHybrid", !!!metr_hybrid_pert))
  
  ### --- XG BOOST 
  
  xgb_model <- train_xgb(X = X_train_scaled,y = sp$y_train)
  
  pred_xgb_smart <- predict_xgb(xgb_model$model, X_test_scaled, y_test= sp$y_test)
  pred_xgb_pert  <- predict_xgb(xgb_model$model, pert_scaled, y_test= y_pert)
  
  metr_xgb_smart <- compute_metrics_from_preds(y_true= pred_xgb_smart$y_test,prob= pred_xgb_smart$probabilities,threshold = 0.5)
  metr_xgb_pert <- compute_metrics_from_preds(y_true= pred_xgb_pert$y_test,prob= pred_xgb_pert$probabilities,threshold = 0.5)
  
  xgb_results <- bind_rows(
    tibble::tibble(dataset = "SmartSeq",  model = "XGBoost", !!!metr_xgb_smart),
    tibble::tibble(dataset = "PerturbSeq", model = "XGBoost", !!!metr_xgb_pert))
  
  #saveRDS(xgb_model, "xgb_model.rds")
  #xgb_model <- readRDS("xgb_model.rds")
  
  all_runs[[b]] <- bind_rows(logi_results, zs_results, rf_results, xgb_results) %>%
    mutate(run = b, seed = seeds[b])
  
  roc_runs[[b]] <- bind_rows(
    tibble(dataset="SmartSeq",  model="Logistic regr.", roc_obj=list(roc(sp$y_test, pred_logi_smart$probabilities, quiet=TRUE))),
    tibble(dataset="SmartSeq",  model="Zero Sum regr.", roc_obj=list(roc(sp$y_test, prob_zs_smart, quiet=TRUE))),
    tibble(dataset="SmartSeq",  model="Random Forest",  roc_obj=list(roc(pred_rf_smart$y_test, pred_rf_smart$probabilities, quiet=TRUE))),
    tibble(dataset="SmartSeq",  model="XGBoost",        roc_obj=list(roc(pred_xgb_smart$y_test, pred_xgb_smart$probabilities, quiet=TRUE))),
    
    tibble(dataset="PerturbSeq", model="Logistic regr.", roc_obj=list(roc(y_pert, pred_logi_pert$probabilities, quiet=TRUE))),
    tibble(dataset="PerturbSeq", model="Zero Sum regr.", roc_obj=list(roc(y_pert, prob_zs_pert, quiet=TRUE))),
    tibble(dataset="PerturbSeq", model="Random Forest",  roc_obj=list(roc(pred_rf_pert$y_test, pred_rf_pert$probabilities, quiet=TRUE))),
    tibble(dataset="PerturbSeq", model="XGBoost",        roc_obj=list(roc(pred_xgb_pert$y_test, pred_xgb_pert$probabilities, quiet=TRUE)))
  ) %>% mutate(run = b)
  
  }

results_repeated <- bind_rows(all_runs)
roc_df_all <- bind_rows(roc_runs)

#saveRDS(results_repeated, "R_objects/results_repeated.rds")
results_repeated <- readRDS("R_objects/results_repeated.rds")

results_repeated %>%
  group_by(dataset, model) %>%
  summarise(
    mean_auc = mean(AUC),
    sd_auc   = sd(AUC),
    mean_balacc = mean(balanced_accuracy),
    sd_balacc   = sd(balanced_accuracy),
    mean_f1 = mean(F1),
    mean_sens = mean(sensitivity),
    mean_spec = mean(specificity),
    .groups = "drop"
  )

delta_auc <- results_repeated %>%
  select(run, model, dataset, AUC) %>%
  pivot_wider(names_from = dataset, values_from = AUC) %>%
  mutate(delta_auc = SmartSeq - PerturbSeq)


###### ======== PLOT ======== ###### 


plot_sum <- results_repeated %>%
  select(dataset, model, balanced_accuracy, AUC, F1, sensitivity, specificity) %>%
  pivot_longer(cols = c(balanced_accuracy, AUC, F1, sensitivity, specificity),
    names_to = "metric", values_to = "value") %>%
  group_by(dataset, model, metric) %>% summarise(mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE), .groups = "drop")

plot_sum$dataset <- factor(plot_sum$dataset, levels = c("SmartSeq", "PerturbSeq"))

plot_sum$metric <- factor(plot_sum$metric,
                         levels = c("balanced_accuracy", "AUC", "F1", "sensitivity", "specificity"),
                         labels = c("Balanced acc.", "AUC", "F1", "Sens.", "Spec."))

plot_sum$model <- factor(plot_sum$model,
                         levels = c("Logistic regression", "ZeroSum", "RandomForest", "XGBoost"),
                         labels = c("Logistic regr.", "Zero Sum regr.", "Random Forest", "XGBoost"))

pal_print <- c(
  "Logistic regr." = "#4D4D4D",
  "Zero Sum regr." = "#9C9C9C",
  "Random Forest"  = "#E69F00",
  "XGBoost"        = "#CC4C02"
)


ggplot(plot_sum, aes(x = metric, y = mean, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
    position = position_dodge(width = 0.7), width = 0.5, linewidth = 0.6) +
  facet_wrap(~ dataset, ncol = 2,
             labeller = labeller(
               dataset = c(
                 "SmartSeq" = "(A)   NMD Smart-Seq2 data \n(source domain)",
                 "PerturbSeq" = "(B)   NMD Perturb-Seq data \n(external validation)"
               )
               )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "", y = "Score", fill = NULL) +
  scale_fill_manual(values = pal_print) +
  theme_minimal(base_size = 20) +
  theme(
    aspect.ratio = 0.8,
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    # (1) facet header box
    strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.3),
    strip.text       = element_text(face = "bold"),
    
    # (3) space between facet panels
    panel.spacing    = grid::unit(1.2, "lines")
  )


# ROC plots 

library(pROC)
library(tibble)

roc_df <- bind_rows(
  tibble(
    dataset = "SmartSeq",
    model   = "Logistic regr.",
    roc_obj = list(roc(sp$y_test, pred_logi_smart$probabilities,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
  ),
  tibble(
    dataset = "SmartSeq",
    model   = "Zero Sum regr.",
    roc_obj = list(roc(sp$y_test, prob_zs_smart,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
  ),
  tibble(
    dataset = "SmartSeq",
    model   = "Random Forest",
    roc_obj = list(roc(pred_rf_smart$y_test, pred_rf_smart$probabilities,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
  ),
  tibble(
    dataset = "SmartSeq",
    model   = "XGBoost",
    roc_obj = list(roc(pred_xgb_smart$y_test, pred_xgb_smart$probabilities,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
  ),
  tibble(
    dataset = "PerturbSeq",
    model   = "Logistic regr.",
    roc_obj = list(roc(y_pert, pred_logi_pert$probabilities,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
  ),
  tibble(
    dataset = "PerturbSeq",
    model   = "Zero Sum regr.",
    roc_obj = list(roc(y_pert, prob_zs_pert,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
  ),
  tibble(
    dataset = "PerturbSeq",
    model   = "Random Forest",
    roc_obj = list(roc(pred_rf_pert$y_test, pred_rf_pert$probabilities,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
  ),
  tibble(
    dataset = "PerturbSeq",
    model   = "XGBoost",
    roc_obj = list(roc(pred_xgb_pert$y_test, pred_xgb_pert$probabilities,
                       levels = c(0, 1), direction = "<", quiet = TRUE))
  )
)

roc_plot_df <- bind_rows(lapply(seq_len(nrow(roc_df)), function(i) {
  g <- pROC::ggroc(roc_df$roc_obj[[i]], legacy.axes = TRUE)
  d <- g$data
  d$dataset <- roc_df$dataset[i]
  d$model   <- roc_df$model[i]
  d$auc     <- as.numeric(pROC::auc(roc_df$roc_obj[[i]]))
  d
}))

roc_plot_df$model <- factor(
  roc_plot_df$model,
  levels = c("Logistic regr.", "Zero Sum regr.", "Random Forest", "XGBoost")
)

roc_plot_df$dataset <- factor(
  roc_plot_df$dataset,
  levels = c("SmartSeq", "PerturbSeq")
)

roc <- ggplot(roc_plot_df, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  facet_wrap(~ dataset, ncol = 2,
             labeller = labeller(
               dataset = c(
                 "SmartSeq"   = "NMD Smart-Seq2 (source domain)",
                 "PerturbSeq" = "NMD Perturb-Seq (external validation)"
               )
             )) +
  scale_color_manual(values = pal_print) +
  coord_equal() +
  labs(
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    panel.spacing = grid::unit(2, "lines")
  )


ggsave(
  filename = "roc_nmd1.png",
  plot = roc,
  path = "Plots/",
  width = 3200,
  height = 2000,
  units = "px",
  dpi = 300
)


#--------------------------------------------------------------------------------
# SIMULATING INCREASING FEATURE MISSINGNESS
#--------------------------------------------------------------------------------

missing_fracs <- seq(0, 0.9, by = 0.1)

# Number of iterations: 
iter <- 10 

# Progressive feature loss: 
set.seed(856)

p <- ncol(X_test_scaled)

all_masked_data <- vector("list", iter)

for (i in 1:iter){
  
  feat_order_smart <- sample(seq_len(p))
  feat_order_pert  <- sample(seq_len(p))
  
  masked_data <- lapply(missing_fracs, function(frac){
    k <- floor(frac * p)
    cs <- feat_order_smart[seq_len(k)]
    cp <- feat_order_pert[seq_len(k)]
    list(
      frac = frac,
      smart_zero = { X <- X_test_scaled; if (k > 0) X[, cs] <- 0; X },
      pert_zero  = { X <- pert_scaled;   if (k > 0) X[, cp] <- 0; X },
      smart_na   = { X <- X_test_scaled; if (k > 0) X[, cs] <- NA_real_; X },
      pert_na    = { X <- pert_scaled;   if (k > 0) X[, cp] <- NA_real_; X }
    )
  })
  
  names(masked_data) <- paste0("frac_", missing_fracs)
  
  all_masked_data[[i]] <- masked_data
}

names(all_masked_data) <- paste0("iter_", 1:iter)




# --- collecting results:
logi_miss <- list()
zs_miss   <- list()
rf_miss   <- list()
hyb_miss  <- list()
rf_drop_miss <- list()
rf_early_miss <- list()
rf_retrain <- list()
xgb_miss  <- list()



# --- LOGISTIC REGRESSION 

for (i in seq_along(all_masked_data)){
  
  masked_data <- all_masked_data[[i]]
  
  for (obj in masked_data) {
    frac <- obj$frac
    
    pred_l_smart <- predict_logistic(logist_model, obj$smart_zero, y_true = sp$y_test)
    pred_l_pert  <- predict_logistic(logist_model, obj$pert_zero, y_true = y_pert)
    
    metr_l_smart <- compute_metrics_from_preds(sp$y_test, prob = pred_l_smart$probabilities, 0.5)
    metr_l_pert  <- compute_metrics_from_preds(y_pert, prob = pred_l_pert$probabilities, 0.5)
    
    logi_miss[[length(logi_miss)+1]] <- bind_rows(
      tibble(dataset="SmartSeq",  missing_frac=frac, model="Logistic",!!!metr_l_smart),
      tibble(dataset="PerturbSeq", missing_frac=frac, model="Logistic",!!!metr_l_pert))
  }
}

logi_miss <- bind_rows(logi_miss)
saveRDS(logi_miss, paste0("R_objects/logi_miss{iter}.rds"))


# ---- ZERO SUM RETRAINING

i_z <- fit_zerosum$lambdaMinIndex  
beta_sparse <- fit_zerosum$coef[[i_z]]
zs_betas <- as.numeric(beta_sparse)
names(zs_betas) <- fit_zerosum$variables.names


for (i in seq_along(all_masked_data)){
  
  masked_data <- all_masked_data[[i]]
  
  for (obj in masked_data){
    
    frac <- obj$frac
    S_masked <- obj$smart_zero
    S_missing_cols <- colnames(S_masked)[colSums(S_masked==0)==nrow(S_masked)]
    
    zs_betas_use <- zs_betas[names(zs_betas) %in% colnames(S_masked)]
    keep_feats <- betas_keep(zs_betas_use, S_missing_cols)
    
    # Retraining
    X_train_reduced <- X_train_scaled[, keep_feats, drop=FALSE]
    X_test_smart    <- X_test_scaled[,  keep_feats, drop=FALSE]
    X_test_pert     <- pert_scaled[, keep_feats, drop=FALSE]
    
    fit_zs_reduced <- zeroSum::zeroSum(
      x = as.matrix(X_train_reduced),
      y = sp$y_train,
      family  = "binomial",
      weights = class_weights,
      alpha   = 0
    )
    saveRDS(fit_zs_reduced, paste0("R_objects/fit_zs_iter", i, "_frac_", frac, ".rds"))

    prob_smart <- zeroSum:::predict.zeroSum(fit_zs_reduced, newx = as.matrix(X_test_smart), type="response")
    prob_pert  <- zeroSum:::predict.zeroSum(fit_zs_reduced, newx = as.matrix(X_test_pert),  type="response")
    
    metr_smart <- compute_metrics_from_preds(sp$y_test, prob_smart, 0.5)
    metr_pert  <- compute_metrics_from_preds(y_pert,   prob_pert,  0.5)
    
    zs_miss[[length(zs_miss)+1]] <- bind_rows(
      tibble(iteration=i, dataset="SmartSeq",  missing_frac=frac, model="ZeroSum_retrained",  !!!metr_smart),
      tibble(iteration=i, dataset="PerturbSeq", missing_frac=frac, model="ZeroSum_retrained", !!!metr_pert)
    )
  }
}
zs_miss <- bind_rows(zs_miss)
saveRDS(zs_miss, "R_objects/zs_miss.rds")



# ----- Random Forest Impurity reduction

for (i in seq_along(all_masked_data)){
  
  masked_data <- all_masked_data[[i]]
  
  for (obj in masked_data){
    
    frac <- obj$frac
  
    prob_smart <- predict(rf_model, data = data.frame(obj$smart_na))$predictions[, "1"]
    prob_pert  <- predict(rf_model, data = data.frame(obj$pert_na))$predictions[, "1"]
    
    metr_smart <- compute_metrics_from_preds(sp$y_test, prob_smart, 0.5)
    metr_pert  <- compute_metrics_from_preds(y_pert, prob_pert,  0.5)
    
    rf_miss[[length(rf_miss)+1]] <- bind_rows(
      tibble(iteration=i, dataset="SmartSeq",  missing_frac=frac, model="RandomForest",    !!!metr_smart),
      tibble(iteration=i, dataset="PerturbSeq", missing_frac=frac, model="RandomForest",   !!!metr_pert)
    )
  }
}

rf_miss <- bind_rows(rf_miss)
saveRDS(rf_miss, "R_objects/rf_miss.rds")


# ---- RANDOM FOREST WITH SURROGATES

for (i in seq_along(all_masked_data)){
  masked_data <- all_masked_data[[i]]
  for (obj in masked_data){

    frac <- obj$frac
    
    res_smart <- predict_hybrid_forest(hybrid_forest_train, obj$smart_na)
    res_pert  <- predict_hybrid_forest(hybrid_forest_train, obj$pert_na)
    
    metr_smart <- compute_metrics_from_preds(sp$y_test, res_smart$probabilities, 0.5)
    metr_pert  <- compute_metrics_from_preds(y_pert, res_pert$probabilities, 0.5)
    
    hyb_miss[[length(hyb_miss) + 1]] <- bind_rows(
      bind_cols(
        tibble(iteration=i, dataset = "SmartSeq",model="SurrogateHybrid", missing_frac = frac, !!!metr_smart),
        as_tibble_row(res_smart$routing_counts)),
      bind_cols(tibble(iteration=i, dataset = "PerturbSeq", model= "SurrogateHybrid",missing_frac = frac, !!!metr_pert),
        as_tibble_row(res_pert$routing_counts)))
  }
}

hyb_miss <- bind_rows(hyb_miss)
saveRDS(hyb_miss, "R_objects/hyb_miss.rds")


# ----- Random Forest Tree dropping 

for (i in seq_along(all_masked_data)){
  masked_data <- all_masked_data[[i]]
  for (obj in masked_data){

    frac <- obj$frac
    
    res_smart <- predict_rf_drop_trees(hybrid_forest_train, obj$smart_na)
    res_pert  <- predict_rf_drop_trees(hybrid_forest_train, obj$pert_na)
    
    # ---- SMARTSEQ
    prob_smart <- res_smart$probabilities
    valid_smart <- !is.na(prob_smart)
    
    metr_smart <- if (sum(valid_smart) == 0 || length(unique(sp$y_test[valid_smart])) < 2) {
      list(accuracy = NA_real_,balanced_accuracy = NA_real_,AUC = NA_real_,precision = NA_real_,recall = NA_real_,sensitivity = NA_real_,specificity = NA_real_,F1 = NA_real_)
    } else {
      compute_metrics_from_preds(sp$y_test[valid_smart],prob_smart[valid_smart],0.5)}
    
    # ---- PERTURBSEQ
    prob_pert <- res_pert$probabilities
    valid_pert <- !is.na(prob_pert)
    
    metr_pert <- if (sum(valid_pert) == 0 || length(unique(y_pert[valid_pert])) < 2) {
      list(accuracy = NA_real_,balanced_accuracy = NA_real_,AUC = NA_real_,precision = NA_real_,recall = NA_real_,sensitivity = NA_real_,specificity = NA_real_,F1 = NA_real_)
    } else {compute_metrics_from_preds(y_pert[valid_pert],prob_pert[valid_pert],0.5)
    }
    
    rf_drop_miss[[length(rf_drop_miss) + 1]] <- bind_rows(
      tibble(iteration=i, dataset = "SmartSeq",model = "RF_DropTrees",missing_frac = frac,n_used_trees = mean(res_smart$n_used_trees),
             frac_no_prediction = mean(!valid_smart),!!!metr_smart),
      tibble(iteration=i, dataset = "PerturbSeq",model = "RF_DropTrees",missing_frac = frac,n_used_trees = mean(res_pert$n_used_trees),
             frac_no_prediction = mean(!valid_pert),!!!metr_pert))
  }
}

rf_drop_miss <- bind_rows(rf_drop_miss)
saveRDS(rf_drop_miss, "R_objects/rf_drop_miss.rds")



# ----- Random Forest Early Stopping

for (i in seq_along(all_masked_data)){
  masked_data <- all_masked_data[[i]]
  for (obj in masked_data){
    frac <- obj$frac
    
    prob_smart <- predict_rf_early_stop(hybrid_forest_train, obj$smart_na)
    prob_pert <- predict_rf_early_stop(hybrid_forest_train, obj$pert_na)
    
    rf_early_miss[[length(rf_early_miss) + 1]] <- bind_rows(
      tibble(iteration=i, dataset = "SmartSeq", model = "RF_EarlyStop", missing_frac = frac,
        !!!compute_metrics_from_preds(sp$y_test, prob_smart, 0.5)), 
      tibble(iteration=i, dataset = "PerturbSeq", model = "RF_EarlyStop", missing_frac = frac,
        !!!compute_metrics_from_preds(y_pert, prob_pert, 0.5)))
  }
}

rf_early_miss <- bind_rows(rf_early_miss)
saveRDS(rf_early_miss, "R_objects/rf_early_miss.rds")


# ------ Retraining 

for (i in seq_along(all_masked_data)){
  masked_data <- all_masked_data[[i]]
  for (obj in masked_data){
    frac <- obj$frac
    
    miss_s <- colnames(obj$smart_na)[colSums(is.na(obj$smart_na)) == nrow(obj$smart_na)]
    keep_s <- setdiff(colnames(X_train_scaled), miss_s)
    
    # train on available features only, predict on masked test
    m_s <- train_rf(X_train_scaled[, keep_s, drop=FALSE], sp$y_train, weights=class_weights, num_trees=100)
    p_s <- predict(m_s, data.frame(obj$smart_na[, keep_s, drop=FALSE]))$predictions[, "1"]
    met_s <- compute_metrics_from_preds(sp$y_test, p_s, 0.5)
    
    # PERT: same idea, but train on full SmartSeq and predict on PerturbSeq
    miss_p <- colnames(obj$pert_na)[colSums(is.na(obj$pert_na)) == nrow(obj$pert_na)]

    p_p <- predict(m_s, data.frame(obj$pert_na[, keep_s, drop=FALSE]))$predictions[, "1"]
    met_p <- compute_metrics_from_preds(y_pert, p_p, 0.5)
    
    rf_retrain[[length(rf_retrain)+1]] <- bind_rows(
      tibble(iteration=i, dataset="SmartSeq",  model="RF_Retrain", missing_frac=frac, !!!met_s),
      tibble(iteration=i, dataset="PerturbSeq", model="RF_Retrain", missing_frac=frac, !!!met_p)
    )
  }
}

rf_retrain <- bind_rows(rf_retrain)
saveRDS(rf_retrain, "R_objects/rf_retrain.rds")


# ----- XG Boost
train_cols <- colnames(X_train_scaled)

for (i in seq_along(all_masked_data)){
  masked_data <- all_masked_data[[i]]
  for (obj in masked_data){

    frac <- obj$frac
    
    Xs <- as.matrix(obj$smart_na[, train_cols, drop = FALSE])
    Xp <- as.matrix(obj$pert_na[, train_cols, drop = FALSE])
    
    pred_x_smart <- predict_xgb(xgb_model$model, Xs, y_test = sp$y_test)
    pred_x_pert  <- predict_xgb(xgb_model$model, Xp, y_test = y_pert)
    
    metr_smart <- compute_metrics_from_preds(sp$y_test, pred_x_smart$probabilities, 0.5)
    metr_pert  <- compute_metrics_from_preds(y_pert, pred_x_pert$probabilities,  0.5)
    
    xgb_miss[[length(xgb_miss) + 1]] <- bind_rows(
      tibble(iteration=i, dataset = "SmartSeq",   missing_frac = frac, model = "XGBoost", !!!metr_smart),
      tibble(iteration=i, dataset = "PerturbSeq", missing_frac = frac, model = "XGBoost", !!!metr_pert)
    )
  }
}

xgb_miss <- bind_rows(xgb_miss)

#saveRDS(xgb_miss, "R_objects/xgb_miss.rds")

logi_miss <- readRDS("R_objects/logi_miss{iter}.rds")
zs_miss <- readRDS("R_objects/zs_miss.rds")
rf_miss <- readRDS("R_objects/rf_miss.rds")
hyb_miss <- readRDS("R_objects/hyb_miss.rds")
xgb_miss <- readRDS("R_objects/xgb_miss.rds")
rf_drop_miss <- readRDS("R_objects/rf_drop_miss.rds")
rf_early_miss <- readRDS("R_objects/rf_early_miss.rds")
rf_retrain <- readRDS("R_objects/rf_retrain.rds")


# ========================================================
# CHOOSE WHAT TO SHOW IN THE PLOT
# ========================================================
models_keep <- c(
  #"Logistic",
  "ZeroSum_retrained",
  "RF_Retrain",
  "RandomForest",
  "SurrogateHybrid",
  #"RF_DropTrees",
  #"RF_EarlyStop",
  "XGBoost"
)

metrics_keep <- c(
  "AUC",
  "balanced_accuracy",
  "F1",
  "sensitivity",
  "specificity"
)
# ======================================================


all_miss <- bind_rows(
  zs_miss,
  logi_miss,
  rf_miss,
  rf_drop_miss,
  rf_early_miss,
  hyb_miss,
  xgb_miss, 
  rf_retrain,
) %>%
  mutate(
    dataset = factor(dataset, levels = c("SmartSeq", "PerturbSeq")),
    model   = factor(model)
  )


metrics_long <- all_miss %>%
  pivot_longer(
    cols = c(accuracy, balanced_accuracy, AUC, F1, sensitivity, specificity),
    names_to = "metric", values_to = "value") %>% 
  filter(model %in% models_keep, metric %in% metrics_keep)


metrics_sum <- metrics_long %>%
  group_by(dataset, model, missing_frac, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = sum(!is.na(value)),
    se   = sd / sqrt(n),
    .groups = "drop")


model_cols <- c(
  "Logistic"          = "#7F7F7F",
  
  "ZeroSum_retrained" = "#FF7F0E",
  "RF_Retrain"        = "#FDE724FF", 
  
  "RandomForest"      = "#D62728",
  "SurrogateHybrid"   = "#9467BD",
  "RF_DropTrees"      = "#7AD151FF", 
  "RF_EarlyStop"      = "#56B4E9",  
  
  "XGBoost"           = "#000004"  
)

model_labels <- c(
  Logistic            = "Logistic regr.",
  ZeroSum_retrained   = "Zero-sum regr.",
  RF_Retrain          = "RF retrained",
  RandomForest        = "RF majority",
  SurrogateHybrid     = "RF surrogates",
  RF_DropTrees        = "RF drop trees",
  RF_EarlyStop        = "RF early stop",
  XGBoost             = "XGBoost"
)


plot_row <- function(metric_name, y_label = NULL) {
  
  df <- metrics_sum %>% filter(metric == metric_name)
  
  df_smart <- df %>% filter(dataset == "SmartSeq")
  df_pert  <- df %>% filter(dataset == "PerturbSeq")
  
  p1 <- ggplot(df_smart, aes(x = missing_frac, 
                             y = mean, color = model, 
                             group = model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_x_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) +
    scale_y_continuous(
      limits = if (metric_name == "AUC") c(0.5, 1) else c(0, 1),
      breaks = if (metric_name == "AUC") {c(0.5, 0.7, 0.9, 1.0)} else {c(0, 0.2, 0.4, 0.6, 0.8, 1.0)}) +
    labs(x = NULL, y = y_label) +
    scale_color_manual(
      values = model_cols,
      breaks = models_keep,
      labels = model_labels[models_keep])+
    guides(color = guide_legend(override.aes = list(
      shape = 16,
      size  = 3,
      linewidth = 0
    ))) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank()
    )

  p2 <- ggplot(df_pert, aes(x = missing_frac, y = mean, color = model, 
                            group = model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_x_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) +
    scale_y_continuous(
      limits = if (metric_name == "AUC") c(0.5, 1) else c(0, 1),
      breaks = if (metric_name == "AUC") {c(0.5, 0.7, 0.9, 1.0)} else {c(0, 0.2, 0.4, 0.6, 0.8, 1.0)}) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(
      values = model_cols,
      breaks = models_keep,
      labels = model_labels[models_keep])+
    guides(color = guide_legend(override.aes = list(
      shape = 16,
      size  = 3,
      linewidth = 0
    ))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank())
  p1 + p2
}


p_auc   <- plot_row("AUC",               "AUC")
p_bacc  <- plot_row("balanced_accuracy", "Balanced accuracy")
p_f1    <- plot_row("F1",                "F1-score")
p_sens  <- plot_row("sensitivity",       "Sensitivity")
p_spec  <- plot_row("specificity",       "Specificity")


header_smart <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, 
           label = "NMD Smart-Seq2 data (source domain)", size = 4) +
  theme_void()

header_pert <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, 
           label = "NMD Perturb-Seq data (external validation)", size = 4) +
  theme_void()

final_plot <-
  (header_smart + header_pert) / 
  p_auc /
  p_bacc /
  p_f1 /
  p_sens /
  p_spec +
  plot_layout(guides = "collect", heights = c(0.8, 1.7, 1.7, 1.7, 1.7, 1.7)) &
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

final_plot



ggsave(
  filename = "p_11.png",
  plot = final_plot,
  path = "Plots/",
  width = 3000,
  height = 4000,
  units = "px",
  dpi = 300
)


# ========================================================
# CHOOSE WHAT TO SHOW IN THE PLOT
# ========================================================
models_keep <- c(
  #"Logistic",
  "ZeroSum_retrained",
  "RF_Retrain",
  "RandomForest",
  "SurrogateHybrid"
  #"RF_DropTrees",
  #"RF_EarlyStop",
  #"XGBoost"
)

metrics_keep <- c(
  "AUC",
  "balanced_accuracy"
  #"F1",
  #"sensitivity",
  #"specificity"
)
# ======================================================


all_miss <- bind_rows(
  zs_miss,
  logi_miss,
  rf_miss,
  rf_drop_miss,
  rf_early_miss,
  hyb_miss,
  xgb_miss, 
  rf_retrain,
) %>%
  mutate(
    dataset = factor(dataset, levels = c("SmartSeq", "PerturbSeq")),
    model   = factor(model)
  )


metrics_long <- all_miss %>%
  pivot_longer(
    cols = c(accuracy, balanced_accuracy, AUC, F1, sensitivity, specificity),
    names_to = "metric", values_to = "value") %>% 
  filter(model %in% models_keep, metric %in% metrics_keep)


metrics_sum <- metrics_long %>%
  group_by(dataset, model, missing_frac, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = sum(!is.na(value)),
    se   = sd / sqrt(n),
    .groups = "drop")


model_cols <- c(
  "Logistic"          = "#7F7F7F",
  
  "ZeroSum_retrained" = "#FF7F0E",
  "RF_Retrain"        = "#FDE724FF", 
  
  "RandomForest"      = "#D62728",
  "SurrogateHybrid"   = "#9467BD",
  "RF_DropTrees"      = "#7AD151FF", 
  "RF_EarlyStop"      = "#56B4E9",  
  
  "XGBoost"           = "#000004"  
)

model_labels <- c(
  Logistic            = "Logistic regr.",
  ZeroSum_retrained   = "Zero-sum regr.",
  RF_Retrain          = "RF retrained",
  RandomForest        = "RF majority",
  SurrogateHybrid     = "RF surrogates",
  RF_DropTrees        = "RF drop trees",
  RF_EarlyStop        = "RF early stop",
  XGBoost             = "XGBoost"
)


plot_row <- function(metric_name, y_label = NULL) {
  
  df <- metrics_sum %>% filter(metric == metric_name)
  
  df_smart <- df %>% filter(dataset == "SmartSeq")
  df_pert  <- df %>% filter(dataset == "PerturbSeq")
  
  p1 <- ggplot(df_smart, aes(x = missing_frac, 
                             y = mean, color = model, 
                             group = model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_x_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) +
    scale_y_continuous(
      limits = if (metric_name == "AUC") c(0.5, 1) else c(0, 1),
      breaks = if (metric_name == "AUC") {c(0.5, 0.7, 0.9, 1.0)} else {c(0, 0.2, 0.4, 0.6, 0.8, 1.0)}) +
    labs(x = NULL, y = y_label) +
    scale_color_manual(
      values = model_cols,
      breaks = models_keep,
      labels = model_labels[models_keep])+
    guides(color = guide_legend(override.aes = list(
      shape = 16,
      size  = 3,
      linewidth = 0
    ))) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank()
    )
  
  p2 <- ggplot(df_pert, aes(x = missing_frac, y = mean, color = model, 
                            group = model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_x_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) +
    scale_y_continuous(
      limits = if (metric_name == "AUC") c(0.5, 1) else c(0, 1),
      breaks = if (metric_name == "AUC") {c(0.5, 0.7, 0.9, 1.0)} else {c(0, 0.2, 0.4, 0.6, 0.8, 1.0)}) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(
      values = model_cols,
      breaks = models_keep,
      labels = model_labels[models_keep])+
    guides(color = guide_legend(override.aes = list(
      shape = 16,
      size  = 3,
      linewidth = 0
    ))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank())
  p1 + p2
}


p_auc   <- plot_row("AUC",               "AUC")
p_bacc  <- plot_row("balanced_accuracy", "Balanced accuracy")
p_f1    <- plot_row("F1",                "F1-score")
p_sens  <- plot_row("sensitivity",       "Sensitivity")
p_spec  <- plot_row("specificity",       "Specificity")


header_smart <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, 
           label = "NMD Smart-Seq2 data (source domain)", size = 4) +
  theme_void()

header_pert <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, 
           label = "NMD Perturb-Seq data (external validation)", size = 4) +
  theme_void()

final_plot <-
  (header_smart + header_pert) / 
  p_auc /
  p_bacc +
  #p_f1 /
  #p_sens /
  #p_spec +
  plot_layout(guides = "collect", 
              #heights = c(0.8, 1.7, 1.7, 1.7, 1.7, 1.7)) &
              heights = c(0.8, 1.7, 1.7)) &
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10))

final_plot



ggsave(
  filename = "p_12.png",
  plot = final_plot,
  path = "Plots/",
  width = 2500,
  height = 1600,
  units = "px",
  dpi = 300
)


#################### NOW AGAIN BUT WITH STANDARD DEVIATION INCLUDED: 

# ========================================================
# CHOOSE WHAT TO SHOW IN THE PLOT
# ========================================================
models_keep <- c(
  #"Logistic",
  "ZeroSum_retrained",
  "RF_Retrain",
  "RandomForest",
  "SurrogateHybrid",
  #"RF_DropTrees",
  #"RF_EarlyStop",
  "XGBoost"
)

metrics_keep <- c(
  "AUC",
  "balanced_accuracy",
  "F1",
  "sensitivity",
  "specificity"
)
# ======================================================


all_miss <- bind_rows(
  zs_miss,
  logi_miss,
  rf_miss,
  rf_drop_miss,
  rf_early_miss,
  hyb_miss,
  xgb_miss, 
  rf_retrain
) %>%
  mutate(
    dataset = factor(dataset, levels = c("SmartSeq", "PerturbSeq")),
    model   = factor(model)
  )


metrics_long <- all_miss %>%
  pivot_longer(
    cols = c(accuracy, balanced_accuracy, AUC, F1, sensitivity, specificity),
    names_to = "metric", values_to = "value"
  ) %>%
  filter(model %in% models_keep, metric %in% metrics_keep)


metrics_sum <- metrics_long %>%
  group_by(dataset, model, missing_frac, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = sum(!is.na(value)),
    se   = sd / sqrt(n),
    .groups = "drop"
  )


model_cols <- c(
  "Logistic"          = "#7F7F7F",
  "ZeroSum_retrained" = "#FF7F0E",
  "RF_Retrain"        = "#FDE724FF",
  "RandomForest"      = "#D62728",
  "SurrogateHybrid"   = "#9467BD",
  "RF_DropTrees"      = "#7AD151FF",
  "RF_EarlyStop"      = "#56B4E9",
  "XGBoost"           = "#000004"
)

model_labels <- c(
  Logistic            = "Logistic regr.",
  ZeroSum_retrained   = "Zero-sum regr.",
  RF_Retrain          = "RF retrained",
  RandomForest        = "RF majority",
  SurrogateHybrid     = "RF surrogates",
  RF_DropTrees        = "RF drop trees",
  RF_EarlyStop        = "RF early stop",
  XGBoost             = "XGBoost"
)


plot_row <- function(metric_name, y_label = NULL) {
  
  df <- metrics_sum %>% filter(metric == metric_name)
  
  df_smart <- df %>% filter(dataset == "SmartSeq")
  df_pert  <- df %>% filter(dataset == "PerturbSeq")
  
  p1 <- ggplot(df_smart, aes(x = missing_frac, y = mean, color = model, group = model)) +
    geom_ribbon(
      aes(ymin = if (metric_name == "AUC") pmax(0.5, mean - sd) else pmax(0, mean - sd),
          ymax = pmin(1, mean + sd),
          fill = model), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_x_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) +
    scale_y_continuous(
      limits = if (metric_name == "AUC") c(0.5, 1) else c(0, 1),
      breaks = if (metric_name == "AUC") c(0.5, 0.7, 0.9, 1.0) else c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
    ) +
    labs(x = NULL, y = y_label) +
    scale_color_manual(
      values = model_cols,
      breaks = models_keep,
      labels = model_labels[models_keep]
    ) +
    scale_fill_manual(
      values = model_cols,
      breaks = models_keep,
      labels = model_labels[models_keep]
    ) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3, linewidth = 0))) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank()
    )
  
  p2 <- ggplot(df_pert, aes(x = missing_frac, y = mean, color = model, group = model)) +
    geom_ribbon(
      aes(ymin = if (metric_name == "AUC") pmax(0.5, mean - sd) else pmax(0, mean - sd),
          ymax = pmin(1, mean + sd),
          fill = model), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_x_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) +
    scale_y_continuous(
      limits = if (metric_name == "AUC") c(0.5, 1) else c(0, 1),
      breaks = if (metric_name == "AUC") c(0.5, 0.7, 0.9, 1.0) else c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
    ) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(
      values = model_cols,
      breaks = models_keep,
      labels = model_labels[models_keep]
    ) +
    scale_fill_manual(
      values = model_cols,
      breaks = models_keep,
      labels = model_labels[models_keep]
    ) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3, linewidth = 0))) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank()
    )
  
  p1 + p2
}


p_auc   <- plot_row("AUC",               "AUC")
p_bacc  <- plot_row("balanced_accuracy", "Balanced accuracy")
p_f1    <- plot_row("F1",                "F1-score")
p_sens  <- plot_row("sensitivity",       "Sensitivity")
p_spec  <- plot_row("specificity",       "Specificity")


header_smart <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = "NMD Smart-Seq2 data (source domain)", size = 4) +
  theme_void()

header_pert <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = "NMD Perturb-Seq data (external validation)", size = 4) +
  theme_void()


final_plot <-
  (header_smart + header_pert) /
  p_auc /
  p_bacc /
  p_f1 /
  p_sens /
  p_spec +
  plot_layout(guides = "collect", heights = c(0.8, 1.7, 1.7, 1.7, 1.7, 1.7)) &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

final_plot


ggsave(
  filename = "p_12.png",
  plot = final_plot,
  path = "Plots/",
  width = 3000,
  height = 4000,
  units = "px",
  dpi = 300
)




#######################################################################################
##### ======== SURROGATE ANALYSIS ============== 
#######################################################################################

hyb_miss <- readRDS("R_objects/hyb_miss.rds")


hyb_summary <- hyb_miss %>%
  mutate(missing_frac = as.numeric(missing_frac), dataset = as.character(dataset)) %>%
  group_by(dataset, missing_frac) %>%
  summarise(
    accuracy = mean(accuracy, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    AUC = mean(AUC, na.rm = TRUE),
    F1 = mean(F1, na.rm = TRUE),
    
    # mean routing counts across iterations
    primary   = round(mean(primary,   na.rm = TRUE)),
    surrogate = round(mean(surrogate, na.rm = TRUE)),
    majority  = round(mean(majority,  na.rm = TRUE)),
    prob      = round(mean(prob,      na.rm = TRUE)),
    
    total_routes = primary + surrogate + majority + prob,
    .groups = "drop"
  ) %>%
  arrange(dataset, missing_frac)

hyb_summary


hyb_summary2 <- hyb_summary %>%
  mutate(
    missing_primary = surrogate + majority + prob,
    fallback = majority + prob,
    surrogate_rate = surrogate / total_routes,
    fallback_rate  = fallback  / total_routes,
    primary_rate   = primary   / total_routes
  ) %>%
  select(dataset, missing_frac, accuracy, balanced_accuracy, AUC, F1,
         primary, missing_primary, surrogate, fallback,
         primary_rate, surrogate_rate, fallback_rate)

hyb_summary2



routing_plot_df <- hyb_summary %>%
  mutate(
    fallback = majority + prob
  ) %>%
  select(dataset, missing_frac, primary, surrogate, fallback) %>%
  pivot_longer(
    cols = c(primary, surrogate, fallback),
    names_to = "routing",
    values_to = "count"
  )

ggplot(routing_plot_df,
       aes(x = missing_frac, y = count, fill = routing)) +
  geom_col(position = "fill") +
  facet_wrap(~dataset, ncol = 2) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Fraction of missing features",
    y = "Routing proportion", fill = "Routing type") +
  scale_fill_manual(values = c(
    primary   = "#4CAF50",
    surrogate = "#9467BD",
    fallback  = "#E64A19"
  )) + theme_minimal(base_size = 18)


routing_plot_df <- hyb_summary %>% mutate(fallback = majority + prob) %>%
  select(dataset, missing_frac, primary, surrogate, fallback) %>%
  pivot_longer(cols = c(primary, surrogate, fallback), names_to = "routing", values_to = "count")

routing_plot_df$dataset <- factor(routing_plot_df$dataset, levels = c("SmartSeq", "PerturbSeq"))

hyb_summary$dataset <- factor(hyb_summary$dataset, levels = c("SmartSeq", "PerturbSeq"))

surr_plot <- ggplot(routing_plot_df,
       aes(x = missing_frac, y = count, fill = routing)) +
  geom_col(position = "fill") +
  # --- AUC LINE
  geom_line(data = hyb_summary, aes(x = missing_frac, y = AUC, group = 1),
    inherit.aes = FALSE, color = "black", linewidth = 1.3) +
  geom_point(data = hyb_summary, aes(x = missing_frac, y = AUC), inherit.aes = FALSE,
    color = "black", size = 2) +
  facet_wrap(~dataset, ncol = 2,
             labeller = labeller(dataset = c(
                 "SmartSeq" = "NMD Smart-Seq2 data \n(source domain)",
                 "PerturbSeq" = "NMD Perturb-Seq data \n(external validation)"))) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Fraction of missing features", y = "Routing proportion / AUC",
    fill = "Routing type") +
  scale_fill_manual(values = c(
    primary   = "#4CAF50",
    surrogate = "#9467BD",
    fallback  = "#E64A19")) +
  theme_minimal(base_size = 27) + 
  theme (legend.position = "bottom", legend.direction = "horizontal")

surr_plot

ggsave(
  filename = "surrogate_3.png",
  plot = surr_plot,
  path = "Plots/",
  width = 5000,
  height = 3000,
  units = "px",
  dpi = 300
  )


##############################################################################
#### RUN MODEL ON PCs 
####################################################################
set.seed(859)
##### ---- Z-scaling 

# this is coloumn wise standardization: 
normalize <- function(x, mean, sd) { 
  return((x - mean) / sd) 
}

B <- 20 
seeds <- sample.int(1e6, B)
all_runs_pca <- vector("list", B)


for (b in seq_len(B)){
  set.seed(seeds[b])
  
  # Scaling of SmartSeq:
  smart_mean <- colMeans(smartSeq)
  smart_sd   <- apply(smartSeq, 2, sd)
  smart_scaled   <- normalize(smartSeq, smart_mean, smart_sd)
  
  sp <- split_data(smart_scaled, y_smart, train_frac = 0.7, seed = seeds[b])
  
  X_train_scaled <- sp$X_train
  X_test_scaled <- sp$X_test
  
  # Scaling of PerturbSeq:
  pert_mean <- colMeans(perturbSeq)
  pert_sd   <- apply(perturbSeq, 2, sd)
  pert_scaled <- normalize(perturbSeq, pert_mean, pert_sd)
  
  
  #### ---- Class weights
  
  n0 <- sum(sp$y_train == 0)
  n1 <- sum(sp$y_train == 1)
  
  w0 <- 1 / n0
  w1 <- 1 / n1
  
  class_weights <- ifelse(sp$y_train == 1, w1, w0)
  
  
  ### PCA 
  pca <- prcomp(X_train_scaled, center = FALSE, scale. = FALSE)
  PC_train <- pca$x
  
  PC_test <- predict(pca, X_test_scaled)
  PC_pert <- predict(pca, pert_scaled)
  
  ## Train RF on genes 
  
  rf_genes <- train_rf(X_train_scaled, sp$y_train, num_trees = 100, max_depth = NULL, weights=class_weights)
  hyb_genes <- build_hybrid_forest_from_ranger(rf_model=rf_genes, X_train  = X_train_scaled, y_train  = sp$y_train)
  
  # train RF on PCs
  rf_pcs <- train_rf(PC_train, sp$y_train, num_trees = 100, max_depth=NULL, weights = class_weights)
  hyb_pcs <- build_hybrid_forest_from_ranger(rf_model = rf_pcs, X_train  = PC_train, y_train  = sp$y_train)
  
  # Prediction genes
  pred_genes_smart <- predict_hybrid_forest(hyb_genes, X_test_scaled)
  pred_genes_pert  <- predict_hybrid_forest(hyb_genes, pert_scaled)
  
  
  metr_genes_smart <- compute_metrics_from_preds(sp$y_test, pred_genes_smart$probabilities, 0.5)
  metr_genes_pert <- compute_metrics_from_preds(y_pert, pred_genes_pert$probabilities, 0.5)
  
  # Predict PCs
  pred_pcs_smart <- predict_hybrid_forest(hyb_pcs, PC_test)
  pred_pcs_pert  <- predict_hybrid_forest(hyb_pcs, PC_pert)
  
  metr_pcs_smart <- compute_metrics_from_preds(sp$y_test, pred_pcs_smart$probabilities, 0.5)
  metr_pcs_pert <- compute_metrics_from_preds(y_pert, pred_pcs_pert$probabilities, 0.5)
  

  all_runs_pca[[b]] <- bind_rows(
    tibble(run = b, seed = seeds[b],
      model = "Genes", dataset = "SmartSeq", !!!metr_genes_smart),
    tibble(run = b, seed = seeds[b], model = "Genes",
      dataset = "PerturbSeq", !!!metr_genes_pert),
    tibble(run = b, seed = seeds[b], model = "PCs",
      dataset = "SmartSeq", !!!metr_pcs_smart),
    tibble(run = b, seed = seeds[b], model = "PCs",
      dataset = "PerturbSeq", !!!metr_pcs_pert))
}

results_compare_all <- bind_rows(all_runs_pca)

summary_pca_nice <- results_compare_all %>%
  group_by(model, dataset) %>%
  summarise(AUC = sprintf("%.3f \u00B1 %.3f",
                  mean(AUC, na.rm = TRUE),
                  sd(AUC, na.rm = TRUE)),
    balanced_accuracy = sprintf("%.3f \u00B1 %.3f",
                                mean(balanced_accuracy, na.rm = TRUE),
                                sd(balanced_accuracy, na.rm = TRUE)), 
    .groups = "drop")

summary_pca_nice


pred_genes_smart$routing_counts
pred_pcs_smart$routing_counts


# genes
cor_genes <- cor(X_train_scaled)
mean(abs(cor_genes[upper.tri(cor_genes)]))

# PCs
cor_pcs <- cor(PC_train)
mean(abs(cor_pcs[upper.tri(cor_pcs)]))


##### MASKING OF FEATURES ON PCS
iter <- 10


# Missingness settings
results_perf <- list()
results_routing <- list()

missing_fracs <- seq(0, 0.9, by = 0.1)

p_genes <- ncol(X_test_scaled)
p_pcs   <- ncol(PC_test)

# same number of masked columns for genes and PCs
# generate masking separately for Smartseq and Perturbseq

for (i in 1:iter) {
  
  set.seed(7777 + i)
  
  for (frac in missing_fracs) {
    
    k_genes <- floor(frac * p_genes)
    k_pcs   <- floor(frac * p_pcs)
    
    mask_smart_genes <- if (k_genes > 0) sample(seq_len(p_genes), k_genes) else integer(0)
    mask_pert_genes  <- if (k_genes > 0) sample(seq_len(p_genes), k_genes) else integer(0)
    
    mask_smart_pcs <- if (k_pcs > 0) sample(seq_len(p_pcs), k_pcs) else integer(0)
    mask_pert_pcs  <- if (k_pcs > 0) sample(seq_len(p_pcs), k_pcs) else integer(0)
    
    # masked sets for genes data :
    
    X_test_genes_miss <- X_test_scaled
    pert_genes_miss   <- pert_scaled
    
    if (k_genes > 0) {
      X_test_genes_miss[, mask_smart_genes] <- NA_real_
      pert_genes_miss[, mask_pert_genes]    <- NA_real_
    }
    
    # masked sets for PCs: 
    
    PC_test_miss <- PC_test
    PC_pert_miss <- PC_pert
    
    if (k_pcs > 0) {
      PC_test_miss[, mask_smart_pcs] <- NA_real_
      PC_pert_miss[, mask_pert_pcs]  <- NA_real_
    }
    
    pred_genes_smart <- predict_hybrid_forest(hyb_genes, X_test_genes_miss)
    pred_genes_pert  <- predict_hybrid_forest(hyb_genes, pert_genes_miss)
    
    metr_genes_smart <- compute_metrics_from_preds(y_true = sp$y_test, prob = pred_genes_smart$probabilities, threshold = 0.5)
    metr_genes_pert <- compute_metrics_from_preds(y_true = y_pert, prob = pred_genes_pert$probabilities, threshold = 0.5)
    
    
    pred_pcs_smart <- predict_hybrid_forest(hyb_pcs, PC_test_miss)
    pred_pcs_pert  <- predict_hybrid_forest(hyb_pcs, PC_pert_miss)
    
    metr_pcs_smart <- compute_metrics_from_preds(y_true = sp$y_test, prob = pred_pcs_smart$probabilities, threshold = 0.5)
    metr_pcs_pert <- compute_metrics_from_preds(y_true = y_pert, prob = pred_pcs_pert$probabilities, threshold = 0.5)
    
    results_perf[[length(results_perf) + 1]] <- bind_rows(
      tibble(iteration=i, representation="Genes", dataset="SmartSeq",
             missing_frac=frac, !!!metr_genes_smart),
      tibble(iteration=i, representation="Genes", dataset="PerturbSeq",
             missing_frac=frac, !!!metr_genes_pert),
      tibble(iteration=i, representation="PCs", dataset="SmartSeq",
             missing_frac=frac, !!!metr_pcs_smart),
      tibble(iteration=i, representation="PCs", dataset="PerturbSeq",
             missing_frac=frac, !!!metr_pcs_pert)
    )
    
    results_routing[[length(results_routing) + 1]] <- bind_rows(
      
      tibble(iteration=i, representation="Genes", dataset="SmartSeq",
             missing_frac=frac,
             primary=pred_genes_smart$routing_counts["primary"],
             surrogate=pred_genes_smart$routing_counts["surrogate"],
             majority=pred_genes_smart$routing_counts["majority"],
             prob=pred_genes_smart$routing_counts["prob"]),
      
      tibble(iteration=i, representation="Genes", dataset="PerturbSeq",
             missing_frac=frac,
             primary=pred_genes_pert$routing_counts["primary"],
             surrogate=pred_genes_pert$routing_counts["surrogate"],
             majority=pred_genes_pert$routing_counts["majority"],
             prob=pred_genes_pert$routing_counts["prob"]),
      
      tibble(iteration=i, representation="PCs", dataset="SmartSeq",
             missing_frac=frac,
             primary=pred_pcs_smart$routing_counts["primary"],
             surrogate=pred_pcs_smart$routing_counts["surrogate"],
             majority=pred_pcs_smart$routing_counts["majority"],
             prob=pred_pcs_smart$routing_counts["prob"]),
  
      tibble(iteration=i, representation="PCs", dataset="PerturbSeq",
             missing_frac=frac,
             primary=pred_pcs_pert$routing_counts["primary"],
             surrogate=pred_pcs_pert$routing_counts["surrogate"],
             majority=pred_pcs_pert$routing_counts["majority"],
             prob=pred_pcs_pert$routing_counts["prob"])
    )
  }
  
  
}


results_perf_tbl <- bind_rows(results_perf)
results_routing_tbl <- bind_rows(results_routing)

results_all_tbl <- results_perf_tbl %>%
  left_join(results_routing_tbl, by = c("iteration", "representation", "dataset", "missing_frac"))

plot_sum <- results_all_tbl %>%
  group_by(representation, dataset, missing_frac) %>%
  summarise(mean_AUC = mean(AUC, na.rm = TRUE),
            sd_AUC = sd(AUC, na.rm = TRUE),
            mean_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
            sd_balanced_accuracy = sd(balanced_accuracy, na.rm = TRUE),
            mean_F1 = mean(F1, na.rm = TRUE),
            sd_F1 = sd(F1, na.rm = TRUE),
            mean_sensitivity = mean(sensitivity, na.rm = TRUE),
            sd_sensitivity = sd(sensitivity, na.rm = TRUE),
            mean_specificity = mean(specificity, na.rm = TRUE),
            sd_specificity = sd(specificity, na.rm = TRUE),
            .groups = "drop"
            )

plot_df <- bind_rows(
  plot_sum %>%
    transmute(
      representation, dataset, missing_frac,
      metric = "AUC", mean = mean_AUC, sd = sd_AUC),
  plot_sum %>%
    transmute(representation, dataset, missing_frac,
      metric = "Balanced accuracy", mean = mean_balanced_accuracy,
      sd = sd_balanced_accuracy),
  plot_sum %>%
    transmute(representation, dataset, missing_frac,
      metric = "F1-score", mean = mean_F1, sd = sd_F1),
  plot_sum %>%
    transmute(
      representation, dataset, missing_frac,
      metric = "Sensitivity", mean = mean_sensitivity, sd = sd_sensitivity),
  plot_sum %>%
    transmute(
      representation, dataset, missing_frac, metric = "Specificity",
      mean = mean_specificity, sd = sd_specificity)
) %>% mutate(dataset = factor(dataset, levels = c("SmartSeq", "PerturbSeq")),
    representation = factor(representation, levels = c("Genes", "PCs")),
    metric = factor(metric, levels = c("AUC", "Balanced accuracy", "F1-score", "Sensitivity", "Specificity")))

model_cols <- c("Genes" = "#D62728", "PCs"   = "#1F77B4")

model_labels <- c("Genes" = "Surrogate RF (genes)", "PCs"   = "Surrogate RF (PCs)")

plot_row <- function(metric_name, y_label = NULL) {
  df <- plot_df %>% filter(metric == metric_name)
  df_smart <- df %>% filter(dataset == "SmartSeq")
  df_pert  <- df %>% filter(dataset == "PerturbSeq")
  p1 <- ggplot(df_smart, aes(x = missing_frac, y = mean, group = representation)) +
    geom_ribbon(
      aes(ymin = if (metric_name == "AUC") pmax(0.5, mean - sd) else pmax(0, mean - sd),
        ymax = pmin(1, mean + sd), fill = representation),
      alpha = 0.18, color = NA, show.legend = FALSE) +
    geom_line(aes(color = representation), linewidth = 0.9, show.legend = FALSE) +
    geom_point(aes(fill = representation), shape = 21, size = 2, show.legend = TRUE) +
    scale_x_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) +
    scale_y_continuous(
      limits = if (metric_name == "AUC") c(0.5, 1) else c(0, 1),
      breaks = if (metric_name == "AUC") c(0.5, 0.7, 0.9, 1.0) else c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    labs(x = NULL, y = y_label) +
    scale_color_manual(values = model_cols, labels = model_labels, guide = "none") +
    scale_fill_manual(values = model_cols, labels = model_labels, guide = guide_legend(
        override.aes = list(shape = 21, size = 4, alpha = 1))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_blank(), strip.text = element_blank(), strip.background = element_blank())
  
  p2 <- ggplot(df_pert, aes(x = missing_frac, y = mean, group = representation)) +
    geom_ribbon(aes( ymin = if (metric_name == "AUC") pmax(0.5, mean - sd) else pmax(0, mean - sd),
        ymax = pmin(1, mean + sd), fill = representation), alpha = 0.18,
      color = NA, show.legend = FALSE) +
    geom_line(aes(color = representation), linewidth = 0.9, show.legend = FALSE) +
    geom_point(aes(fill = representation), shape = 21, size = 2, show.legend = TRUE) +
    scale_x_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9)) +
    scale_y_continuous(limits = if (metric_name == "AUC") c(0.5, 1) else c(0, 1),
      breaks = if (metric_name == "AUC") c(0.5, 0.7, 0.9, 1.0) else c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(values = model_cols, labels = model_labels, guide = "none") +
    scale_fill_manual(values = model_cols, labels = model_labels, guide = guide_legend(
        override.aes = list(shape = 21, size = 4, alpha = 1))) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", plot.title = element_blank(),
      strip.text = element_blank(), strip.background = element_blank())
  p1 + p2
}

p_auc   <- plot_row("AUC", "AUC")
p_bacc  <- plot_row("Balanced accuracy", "Balanced accuracy")
p_f1    <- plot_row("F1-score", "F1-score")
p_sens  <- plot_row("Sensitivity", "Sensitivity")
p_spec  <- plot_row("Specificity", "Specificity")

header_smart <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "NMD SmartSeq2 \n(source domain)", size = 7) +
  theme_void()

header_pert <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "NMD PerturbSeq \n(external validation)", size = 7) +
  theme_void()

final_plot <-
  (header_smart + header_pert) /
  p_auc /
  p_bacc /
  p_f1 /
  p_sens /
  p_spec + plot_layout(
    guides = "collect",
    heights = c(0.6, 1.5, 1.5, 1.5, 1.5, 1.5)
  ) & theme(legend.position = "bottom",
    legend.title = element_blank(), legend.text = element_text(size = 22))

final_plot

ggsave(
  filename = "pcas_features4.png",
  plot = final_plot,
  path = "Plots/",
  width = 3000,
  height = 3500,
  units = "px",
  dpi = 300
)


###################################################################################

routing_overview <- results_routing_tbl %>%
  mutate(
    total_routes = primary + surrogate + majority + prob,
    surrogate_pct = 100 * surrogate / total_routes
  ) %>%
  group_by(representation, missing_frac) %>%
  summarise(
    mean_surrogate_pct = mean(surrogate_pct, na.rm = TRUE),
    sd_surrogate_pct   = sd(surrogate_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    remaining_features = case_when(
      representation == "Genes" ~ p_genes - floor(missing_frac * p_genes),
      representation == "PCs"   ~ p_pcs   - floor(missing_frac * p_pcs)
    )
  ) %>%
  select(missing_frac, representation, remaining_features, mean_surrogate_pct) %>%
  pivot_wider(
    names_from = representation,
    values_from = c(remaining_features, mean_surrogate_pct)
  ) %>%
  rename(
    remaining_genes = remaining_features_Genes,
    remaining_pcs = remaining_features_PCs,
    surrogate_genes = mean_surrogate_pct_Genes,
    surrogate_pcs = mean_surrogate_pct_PCs
  ) %>%
  mutate(
    surrogate_genes = round(surrogate_genes, 2),
    surrogate_pcs   = round(surrogate_pcs, 2)
  ) %>%
  arrange(missing_frac)

routing_overview


colnames(routing_overview) <- c(
  "Missing fraction",
  "Remaining genes",
  "Remaining PCs",
  "Surrogate routing (%) genes",
  "Surrogate routing (%) PCs"
)

print(routing_overview)
