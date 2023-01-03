# __________________________________________________________________________________________________
# Boosted Regression trees
# Resources for mlr3 with xgboost:
# https://github.com/KI-Research-Institute/LearningWithExternalStats/blob/f6ead01200c2a85ac0d8e5c67f0d99a7061eabb2/extras/mlWrappers/wxgboost.R
# https://mlr3book.mlr-org.com/performance.html

# TODO, read on cross validation
# https://machinelearningmastery.com/k-fold-cross-validation/
# And on boosting parameters (e.g. understand nrounds better)
# Variable importance in boosting
# __________________________________________________________________________________________________

# CWM values & TU data
data_cwm <- readRDS(file.path(path_cache, "data_cwm.rds"))
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))

# Combine max tu and cwm for each region
data_cwm <- lapply(data_cwm, function(y)
  y[max_tu, max_log_tu := i.max_log_tu, on = c(site = "TSITE_NO_WQ")])

# TODO: Site T12073525 does not exist in chemical data
# max_tu[TSITE_NO_WQ %like% "T12073525", ]
data_cwm$PN <- data_cwm$PN[site != "T12073525",]

# Trait names
trait_names <-
  grep("feed.*|resp.*|locom.*|size.*|volt.*|sensitivity.*",
       names(dat),
       value = TRUE)
# dat <- data_cwm$California
# dat <- dcast(dat, site + max_log_tu ~ trait, value.var = "cwm_val")

## Boosting ----
perfor_ls <- list()
plot_ls <- list()
imp_ls <- list()
model_ls <- list()

for(region in names(data_cwm)[c(1,3:5)]) {
  set.seed(1234)
  
  # Data
  dat <- data_cwm[[region]]
  dat <-
    dcast(dat, site + max_log_tu ~ trait, value.var = "cwm_val")
  
  # Task
  task <- TaskRegr$new(id = region,
                       backend = dat[, .SD, .SDcols = c(trait_names,
                                                        "max_log_tu")],
                       target = "max_log_tu")
  # Learner
  xgboost_learner <- lrn(
    "regr.xgboost",
    objective = "reg:squarederror",
    # for the learning task
    eval_metric = 'rmse',
    # evaluation for test data set
    gamma = 10,
    # Minimum loss reduction required to make a further partition on a leaf node.
    # The larger, the more conservative the algorithm will be
    booster = "gbtree",
    eta = to_tune(p_dbl(lower = 0.001, upper = 0.3)),
    # controlling the learning rate (scaling of the contrib. of each tree by a factor 0-1 when added)
    # to prevent overfitting, low values means more nrounds (robust to overfitting but also slower to compute)
    max_depth = 6L,
    nrounds = to_tune(p_int(lower = 1, upper = 50)),
    # number of boosting rounds
    predict_type = "response"
  )
  
  # Tuning
  # use the auto_tuner function, no need to extract the best set of hyperparameters at the end
  instance <- tune(
    method = tnr("grid_search", resolution = 5),
    task = task,
    learner = xgboost_learner,
    resampling = rsmp("cv", folds = 3),
    measures = msr("regr.mse"),
    terminator = trm("evals", n_evals = 5)
  )
  # instance$archive
  # instance$result
  
  # Model performance with nested resampling
  # We could use nested resampling to evaluate the model performance
  # Nested resampling is an additional step after fitting the final model and should not
  # be used to find optimal hyperparameters.
  # Idea: If the same data for model selection and evaluation of the model is used, we bias the
  # performance estimate (eval. on test data could leak information about the test data's structure into
  # the model)
  # Steps:
  # - outer resampling: cv to get different testing an training datasets
  # - inner resampling: Within the training data use cv to get different inner testing and training data sets
  # - Tuning the hyperparameters with the inner data splits
  # - Fits learner on the outer training data set with the tuned hyperparameters from the inner resampling
  # - Evaluates the performance of the learner on the outer testing data
  # - Repeat inner resampling + fitting + evaluation for all folds of the outer resampling
  # - Aggregate perfromance values to get an unbiased perfromance estimate
  # Could also use a auto tuner with a different number of cv
  at <- auto_tuner(
    method = tnr("grid_search", resolution = 5),
    learner = xgboost_learner,
    resampling = rsmp("cv", folds = 4),
    # for inner resampling
    measure = msr("regr.mse"),
    terminator = trm("evals", n_evals = 5),
  )
  
  outer_resampling <- rsmp("cv", folds = 3)
  rr <- resample(task,
                 at,
                 outer_resampling,
                 store_models = TRUE) # investigate inner tuning with store models set to TRUE
  # rr$predictions()
  
  perfor_ls[[region]] <- list(
    "perform_inner_rsmp" = extract_inner_tuning_results(rr),
    # Performance results on the inner resampling
    # These values should not be used to fit a final model
    "perform_outer_resamp" = rr$score(),
    # Outer resampling performance results
    # Only problematic if there's a huge difference
    "agg_perform_outer_rsmp" = rr$aggregate()
    # Aggregated performance of the outer resampling should be reported (unbiased with optimal
    # hyperparameters)
  )
  
  ## Final Model ----
  xgboost_learner$param_set$values <-
    instance$result_learner_param_vals
  xgboost_learner$train(task)
  model_ls[[region]] <- xgboost_learner$model
  
  ## Interpret model ----
  xgboost_exp <- explain_mlr3(
    xgboost_learner,
    data     = dat[, .SD, .SDcols = trait_names],
    y        = dat$max_log_tu,
    label    = "XGBOOST",
    colorize = FALSE
  )
  
  ### Variable importance ----
  # Uses mean dropout loss for xgboost
  imp_ls[[region]] <-
    list(
      "permutation_imp" = model_parts(xgboost_exp),
      "model_imp" = xgboost_learner$importance()
    )
  
  # Similar but not equal to the importance reported from the model
  plot_ls[[region]] <-
    list("plot_permutation_imp" = plot(model_parts(xgboost_exp))) 
  
  ### PDPs ----
  # TODO:Might extend to PDPs and even further
  # pdp <- model_profile(
  #   xgboost_exp,
  #   variables = c(
  #     "locom_swim",
  #     "feed_predator",
  #     "size_medium",
  #     "sensitivity_organic",
  #     "size_large"
  #   )
  # )$agr_profiles
  #
  # plot(pdp) +
  #   scale_y_continuous("Max Log TU") +
  #   theme_bw()
  
}
# TODO: Check performance values in detail
# perfor_ls
# model_ls

# Plots permutation importance
plot_ls$California
ggsave(
  filename = file.path(path_out,
                       "Graphs",
                       "Most_imp_traits_California.png"),
  width = 50,
  height = 30,
  units = "cm"
)

# Overview 5 most important traits  
permutation_imp <- lapply(imp_ls, function(y) y$permutation_imp)
most_imp_permutation_imp <- lapply(permutation_imp, function(y)
  as.data.table(y) %>%
    .[, .(mean_dropout_loss = mean(dropout_loss)), by = "variable"] %>%
    .[!variable %in% c("_baseline_", "_full_model_"), ] %>%
    .[order(-mean_dropout_loss), .SD[1:5, ]])
most_imp_permutation_imp <- rbindlist(most_imp_permutation_imp, idcol = "region")
most_imp_permutation_imp[, region_QA := fcase(
  region == "California",
  "CSQA",
  region == "Northeast",
  "NESQA",
  region == "PN",
  "PNSQA",
  region == "Southeast",
  "SESQA"
)] # TODO: this needs to be standardized in all figures and tables
most_imp_permutation_imp[, mean_dropout_loss := round(mean_dropout_loss, digits = 2)]
setcolorder(most_imp_permutation_imp,
            c("region", "region_QA"))
fwrite(most_imp_permutation_imp,
       file = file.path(path_paper,
                        "Tables",
                        "Most_imp_traits.csv"))
