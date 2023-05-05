# __________________________________
# Testing boosted regression trees
# __________________________________
data_cwm <- readRDS(file.path(path_cache, "data_cwm.rds"))
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

# traits
trait_names <- unique(data_cwm$California$trait)

# Combine max tu and cwm for each region
# make an exception for Midwest (merge via STAID)
data_cwm <- lapply(data_cwm, function(y) {
  on <- if ("STAID" %in% names(y)) c("STAID" = "site") else "site"
  y[max_tu, max_log_tu := i.max_log_tu, on = on]
})
data <- dcast(data_cwm$California, site + max_log_tu ~ trait,
  value.var = "cwm_val"
)
# View(data)

# XGBoost basic R implementation
ind <- sample(1:nrow(data), size = round(nrow(data) * 0.8))

# Training set
dtrain <- data[ind, .SD, .SDcols = trait_names]
tu_train <- data[ind, .SD, .SDcols = "max_log_tu"]
dtrain <- as.matrix(dtrain)
tu_train <- as.matrix(tu_train)
dtrain <- xgb.DMatrix(dtrain,
  label = tu_train
)

# Test set
dtest <- data[-ind, .SD, .SDcols = trait_names]
tu_test <- data[-ind, .SD, .SDcols = "max_log_tu"]
dtest <- as.matrix(dtest)
tu_test <- as.matrix(tu_test)
dtest <- xgb.DMatrix(dtest,
                     label = tu_test)
watchlist <- list(train = dtrain,
                  eval = dtest)

## A simple xgb.train example:
param <- list(
  max_depth = 15,
  eta = 0.1,
  verbose = 0,
  nthread = 2,
  objective = "reg:squarederror",
  eval_metric = "rmse"
)
# param <- list(
#   max_depth = 15,
#   alpha = 6.743465,
#   lambda = 5.793523,
#   eta = 0.1949129,
#   subsample = 0.5, # 0.3164783
#   verbose = 0,
#   nthread = 2,
#   objective = "reg:squarederror",
#   eval_metric = "rmse"
# )
bst <- xgb.train(param, dtrain,
  nrounds = 100,
  watchlist
)
# gives probably prediction on test data

# Using conditional inference trees
# library(party)
# library(mboost)

# MLR3 approach with tuning ----
task <- TaskRegr$new(id = "Cal",
                     backend = data[, .SD, .SDcols = c(trait_names,
                                                           "max_log_tu")],
                     target = "max_log_tu")
# Learner
xgboost_learner <- lrn(
  "regr.xgboost",
  objective = "reg:squarederror",
  gamma = to_tune(p_dbl(lower = 0, upper = 10)),
  # Minimum loss reduction required to make further
  # prediction on a leaf node of a tree (the larger
  # the more conservative)
  # alpha = to_tune(p_dbl(lower = 0, upper = 10)),
  # removes features in the end, L1 regularization
  # lambda = to_tune(p_dbl(lower = 0, upper = 10)),
  # shrinks features without removing them, L2 regularization
  booster = "gbtree",
  eta = to_tune(p_dbl(lower = 0.001, upper = 1)),
  # Learning rate, scaling the contrib. of 
  # each treee
  max_depth = to_tune(p_int(lower = 2L, upper = 15L)),
  subsample = to_tune(p_dbl(lower = 0.1, upper = 1)),
  nrounds = 100,
  # to_tune(p_int(lower = 10, upper = 50))
  predict_type = "response"
)

# Tuning
instance <- tune(
  method = tnr("irace"),
  task = task,
  learner = xgboost_learner,
  resampling = rsmp("cv", folds = 3),
  measures = msr("regr.rmse"), # evaluation performance of training data
  terminator = trm("evals", n_evals = 1000)
)

# Performance on validation data (holdout from cv)
summary_funs <- list(
  "min" = min,
  "max" = max,
  "mean" = mean,
  "median" = median
)
as.data.table(instance$archive)[, lapply(summary_funs, function(f)
  f(regr.rmse))]
as.data.table(instance$archive)[regr.rmse == min(regr.rmse), ]
# autoplot(instance, type = "surface")

# Check predictions to see if there are any problems
xgboost_learner$param_set$values <- instance$result_learner_param_vals
xgboost_learner$train(task)
predictions <- xgboost_learner$predict(task)
# predictions$score()
# Performance on the training data: regr.rmse: 0.54
autoplot(predictions)
predictions$score()

as.data.table(predictions)[, .(row_ids, sqdiff = (truth - response) ^ 2, truth, response)] %>%
  .[order(-sqdiff),]
