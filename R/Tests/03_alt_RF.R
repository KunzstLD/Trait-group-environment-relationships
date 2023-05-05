# CWM values & TU data
data_cwm <- readRDS(file.path(path_cache, "data_cwm.rds"))
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

# Combine max tu and cwm for each region
# make an exception for Midwest (merge via STAID)
data_cwm <- lapply(data_cwm, function(y) {
  on <- if("STAID" %in% names(y)) c("STAID" = "site") else "site"
  y[max_tu, max_log_tu := i.max_log_tu, on = on]
})

# There are sites that don't have chemical data
# lapply(data_cwm, function(y)
#   y[is.na(max_log_tu), ])  
# Few sites from Midwest:
# unique(data_cwm$Midwest[is.na(max_log_tu), STAID])
# T03611200, T03318800, T393247089260701
# and PN: T12073525
# unique(data_cwm$PN[is.na(max_log_tu), site])

# Remove sites with no chemical information 
# Note to Ian that there were some duplicate STAID (but different sites)
# for Midwest and that for three sites max_log_tu could not be matched
# max_tu[site %like% "T12073525", ]
# data_cwm$PN[site == "T12073525",]
data_cwm$Midwest <- data_cwm$Midwest[!is.na(max_log_tu), ]
data_cwm$PN <- data_cwm$PN[!is.na(max_log_tu),]

# Trait names
trait_names <- unique(data_cwm$California$trait)

## Boosting ----
perfor_ls <- list()
model_ls <- list()
archive_ls <- list()
imp_ls <- list()
pdp_ls <- list()

for(region in names(data_cwm)[1]) {
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
  rf_learner <- lrn("regr.ranger",
                    predict_type = "prob",
                    importance = "permutation",
                    mtry = p_int(lower = 1, upper = length(dat) -1), 
                    min.node.size = p_int(lower = 1, upper = 10))
  
  # Tuning
  # use the auto_tuner function, no need to extract the best set of hyperparameters at the end
  instance <- tune(
    method = tnr("grid_search", resolution = 20),
    task = task,
    learner = rf_learner,
    resampling = rsmp("cv", folds = 3),
    measures = msr("regr.rmse"), # evaluation performance of training data
    terminator = trm("none")
  )
  archive_ls[[region]] <- instance$archive
  
  # Final Model
  rf_tuned <- lrn("regr.xgboost", id = "Xgboost tuned")
  rf_tuned$param_set$values <- instance$result_learner_param_vals
  rf_tuned$train(task)
  model_ls[[region]] <- rf_tuned$model
  
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
  # - Aggregate perfromance values to get an unbiased performance estimate
  # Could also use a auto tuner with a different number of cv
  at <- auto_tuner(
    method = tnr("grid_search", resolution = 20),
    learner = rf_learner,
    resampling = rsmp("cv", folds = 3),
    # for inner resampling
    measure = msr("regr.rmse"),
    terminator = trm("none"),
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
    "perform_full_archive" = extract_inner_tuning_archives(rr),
    # Full tuning archive, check if hyperparameters converge? 
    "perform_outer_resamp" = rr$score(),
    # Outer resampling performance results
    # Only problematic if there's a huge difference between outer and inner resampling
    "agg_perform_outer_rsmp" = rr$aggregate()
    # Aggregated performance of the outer resampling should be reported (unbiased with optimal
    # hyperparameters)
  )
  
  # Interpret model
  RF_exp <- explain_mlr3(
    rf_tuned,
    data     = dat[, .SD, .SDcols = trait_names],
    y        = dat$max_log_tu,
    label    = "RF",
    colorize = FALSE
  )
  
  # Variable importance
  # Uses mean dropout loss for xgboost?
  imp_ls[[region]] <-
    list(
      "permutation_imp" = model_parts(RF_exp,
                                      B = 50,
                                      loss_function = loss_root_mean_square),
      "model_imp" = rf_tuned$importance()
    )
  
  # PDPs
  pdp_ls[[region]] <- model_profile(RF_exp)$agr_profiles
  
}
