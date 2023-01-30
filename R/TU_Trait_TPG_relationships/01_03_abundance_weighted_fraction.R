# ________________________________________________
# Abundance weighted fraction of TPGs ----
# Same families as in the CWM approach should be covered
# TODO: 
# - 3 families cannot be classified to a certain TPG
# (California)
## _______________________________________________

# Abundance weighted fraction
trait_family <- readRDS(file.path(path_cache, "trait_family_tpg.rds"))
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))
abund <- abund[Region == "California", .(
    site,
    species,
    genus,
    family,
    order,
    taxon,
    taxonomic_level,
    abundance
)]

# Merge TPG
# No 11 and 14
abund[trait_family, group := i.group, on = "family"]
# abund[, .N, by = "group"] %>% 
# .[order(group), ]

# 3 families which are not classified to a certain TPG
# (California)
# throw out for now
# unique(abund[is.na(group), .(family, order)])
# abund[is.na(group), abundance]
# noa_trait_conus[family %in% c(
#   "Ptiliidae",
#   "Limoniidae",
#   "Dolichopodidae"
# ), ]
# noa_zuelig[family %in% c(
#   "Ptiliidae",
#   "Limoniidae",
#   "Dolichopodidae"
# ), 
abund <- abund[!is.na(group),]

# Calc weighted fraction
abund[, tot_abund_site := sum(abundance), by = "site"]
abund[, `:=`(
    weighted_fraction = sum(abundance),
    weighted_fraction_rel = sum(abundance) / tot_abund_site
),
by = c(
    "site",
    "group"
)
]

# TPGs as columns for BRT analyses
abund[, group := paste0("T", group)]
abund <- abund[, .(site, group, weighted_fraction)] %>%
    unique()
trait_groups <- dcast(abund[, .(site, group, weighted_fraction)],
    site ~ group,
    value.var = "weighted_fraction"
)

# Add max TU
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

trait_groups[max_tu, max_log_tu := i.max_log_tu, on = "site"]
var_names <- names(trait_groups)[names(trait_groups) %like% "^T[0-9]"]


############################ Test/Experimental ########################################
# XGBoost
ind <- sample(1:nrow(trait_groups), size = round(nrow(trait_groups) * 0.8))

# Training set
dtrain <- trait_groups[ind, .SD, .SDcols = var_names]
tu_train <- trait_groups[ind, .SD, .SDcols = "max_log_tu"]
dtrain <- as.matrix(dtrain)
tu_train <- as.matrix(tu_train)
dtrain <- xgb.DMatrix(dtrain,
  label = tu_train
)

# Test set
dtest <- trait_groups[-ind, .SD, .SDcols = var_names]
tu_test <- trait_groups[-ind, .SD, .SDcols = "max_log_tu"]
dtest <- as.matrix(dtest)
tu_test <- as.matrix(tu_test)
dtest <- xgb.DMatrix(dtest,
                     label = tu_test)
watchlist <- list(train = dtrain,
                  eval = dtest)

## A simple xgb.train example:
set.seed(1234)
param <- list(
  max_depth = 15,
  eta = 0.95,
  verbose = 0,
  nthread = 2,
  objective = "reg:squarederror",
  eval_metric = "rmse"
)
bst <- xgb.train(param, dtrain,
  nrounds = 200,
  watchlist
)


# Hyperparam. optimization with MLR3
task <- TaskRegr$new(id = "Cal",
                     backend = trait_groups[, .SD, .SDcols = c(var_names,
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
  max_depth = 6L, # to_tune(p_int(lower = 2L, upper = 30L)),
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
# Performance on the training data: regr.rmse: 2.73
autoplot(predictions)
predictions$score()
as.data.table(predictions)[, .(row_ids, sqdiff = (truth - response) ^ 2, truth, response)] %>%
  .[order(-sqdiff),] %>% 
  .[sqdiff >= 1, ]


# Use the previous optimal config and prevent overfitting with 
# regularization
xgboost_learner_new <- lrn(
  "regr.xgboost",
  objective = "reg:squarederror",
  gamma = 9.000126,
  # Minimum loss reduction required to make further
  # prediction on a leaf node of a tree (the larger
  # the more conservative)
  alpha = to_tune(p_dbl(lower = 0, upper = 10)),
  # removes features in the end, L1 regularization
  lambda = to_tune(p_dbl(lower = 0, upper = 10)),
  # shrinks features without removing them, L2 regularization
  booster = "gbtree",
  eta = 0.03408398,
  # Learning rate, scaling the contrib. of 
  # each treee
  max_depth = 6L, # to_tune(p_int(lower = 2L, upper = 30L)),
  subsample = 0.3199171,
  nrounds = 100,
  # to_tune(p_int(lower = 10, upper = 50))
  predict_type = "response"
)
instance_new <- tune(
  method = tnr("irace"),
  task = task,
  learner = xgboost_learner_new,
  resampling = rsmp("cv", folds = 3),
  measures = msr("regr.rmse"), # evaluation performance of training data
  terminator = trm("evals", n_evals = 1000)
)
as.data.table(instance_new$archive)[, lapply(summary_funs, function(f)
  f(regr.rmse))]
as.data.table(instance_new$archive)[regr.rmse == min(regr.rmse), ]
