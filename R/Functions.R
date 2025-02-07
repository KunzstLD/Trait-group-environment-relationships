# ______________________________________________________________________
# Collection of functions used during data processing & analysis 
# ______________________________________________________________________

# Check for completeness of trait dataset ----
create_pattern_ind <- function(x, non_trait_cols) {
  if (missing(non_trait_cols)) {
    trait_names_pattern <- sub("\\_.*|\\..*", "", names(x)) %>%
      unique() %>%
      paste0("^", .)
  } else{
    pat <- paste0(non_trait_cols, collapse = "|")
    # get trait names & create pattern for subset
    trait_names_pattern <-
      grep(pat, names(x), value = TRUE, invert = TRUE) %>%
      sub("\\_.*|\\..*", "", .) %>%
      unique() %>%
      paste0("^", .)
  }
  trait_names_pattern
}

completeness_trait_data <- function(x, non_trait_cols) {
  trait_names_pattern <- create_pattern_ind(
    x = x,
    non_trait_cols = non_trait_cols
  )
  
  # test how complete trait sets are
  output <- matrix(ncol = 2, nrow = length(trait_names_pattern))
  for (i in seq_along(trait_names_pattern)) {
    # vector containing either 0 (no NAs) or a number (> 0) meaning that all
    # entries for this trait contained NA
    vec <-
      x[, apply(.SD, 1, function(y) {
        base::sum(is.na(y))
      }),
      .SDcols = names(x) %like% trait_names_pattern[[i]]
      ]
    
    # How complete is the dataset for each individual trait?
    output[i, ] <-
      c(
        (length(vec[vec == 0]) / nrow(x)) %>% `*`(100) %>% round(),
        trait_names_pattern[[i]]
      )
  }
  output <- as.data.table(output)
  return(output)
}

# Normalization of trait scores ----
# All trait states of one trait are divided by their row sum
# Hence, trait affinities are represented as "%" or ratios
# this works for traits named: "groupingfeature_trait"
# Use for data in wide format!
normalize_by_rowSum <- function(x,
                                non_trait_cols,
                                na.rm = TRUE) {
  # get trait names & create pattern for subset
  trait_names_pattern <- create_pattern_ind(
    x = x,
    non_trait_cols = non_trait_cols
  )

  # loop for normalization (trait categories for each trait sum up to 1)
  for (cols in trait_names_pattern) {
    # get row sum for a specific trait
    x[, rowSum := apply(.SD, 1, sum, na.rm = na.rm),
      .SDcols = names(x) %like% cols
    ]

    # get column names for assignment
    col_name <- names(x)[names(x) %like% cols]

    # divide values for each trait state by
    # the sum of trait state values
    x[, (col_name) := lapply(.SD, function(y) {
      y / rowSum
    }),
    .SDcols = names(x) %like% cols
    ]
  }
  # del rowSum column
  x[, rowSum := NULL]
  return(x)
}


# Trait aggregation ----
# Direct aggregation to family-level
# Works for data.tables
direct_agg <- function(trait_data,
                       non_trait_cols,
                       method,
                       na.rm = TRUE,
                       on = "family") {
  if (!on %in% c("genus", "family")) {
    warning(paste("Aggregation to", on, "is not implemented yet"))
  }
  # get names of trait columns
  pat <- paste0(non_trait_cols, collapse = "|")
  trait_col <-
    grep(pat, names(trait_data), value = TRUE, invert = TRUE)
  
  # aggregate to family-level
  # & merge information on order back
  # subset so that no NA values occur in data
  # (otherwise all NA entries are viewed as a group &
  # aggregated as well)
  if (on == "genus") {
    agg_data <- trait_data[!is.na(genus),
                           lapply(.SD, method, na.rm = na.rm),
                           .SDcols = trait_col,
                           by = on]
    agg_data[trait_data, `:=`(family = i.family,
                              order = i.order), on = on]
    return(agg_data)
  }
  if (on == "family") {
    agg_data <- trait_data[!is.na(family),
                           lapply(.SD, method, na.rm = na.rm),
                           .SDcols = trait_col,
                           by = on]
    agg_data[trait_data, order := i.order, on = on]
    return(agg_data)
  }
}

# Statistical (helper) functions ----

# List of stat. summary functions 
summary_funs <- list(
  "min" = min,
  "max" = max,
  "mean" = mean,
  "median" = median
)

## Geometric mean ---- 
# test <- rnorm(n = 10000, mean = 2, sd = 0.1)
# Problem of overflow, will result in Inf for large lists
# gmean <- function(x){
#   (Reduce("*", x))^(1/length(x))
# }
# TODO: Understand the relation to the logarithm
gmean <- function(x) {
  exp(mean(log(x)))
}

## Simple regression ----
# For inside use in data.tables
# Old version: 
# summary(lm(form, data = x)) %>%
  #     coef(.) %>%
  #     as.data.table(., keep.rownames = "id")
sregr_dt <- function(x, form) {
  lm_summary <- summary(lm(form, data = x))
  cbind(
    coef(lm_summary),
    "adj_r_squared" = lm_summary["adj.r.squared"]
  ) %>%
    as.data.table(., keep.rownames = "id")  
}

# Create a data.table out of an lm summary 
lm_summary_to_dt <- function(lm_obj) {
  lm_summary <- summary(lm_obj)
  data.table(
    coef(lm_summary),
    "r_squared" = lm_summary[["r.squared"]],
    keep.rownames = "id"
  )
}
  
## Interactions ----
fun_interactions <- function(x, formula) {
    lapply(formula, function(fo) {
        lm_obj <- lm(fo, data = x)
        cbind(
            broom::tidy(lm_obj),
            broom::glance(lm_obj)[c("r.squared", "adj.r.squared")],
            fo
        )
    }) %>%
        rbindlist(.)
}

## Cluster analysis ----
mycluster_hc <- function(x, k) {
  list(cluster = cutree(hclust(as.dist(x),
                               method = "ward.D2"),
                        k = k))
}

clustering_traits <- function(x, # data.frame or data.table
                              fc_traits, # fuzzy coded traits
                              q_traits, # quantiative traits
                              taxa,
                              k_max = 20) {
  setDF(x)
  # Add row.names
  row.names(x) <- taxa

  # Convert to ktab object
  vec <- sub("\\_.*", "\\1", fc_traits)
  blocks <- rle(vec)$lengths
  prep_fuzzy_data <- prep.fuzzy(x[, fc_traits], blocks)
  ktab <- ktab.list.df(list(
    prep_fuzzy_data,
    x[, q_traits, drop = FALSE]
  ))
  dist_mat <-
    dist.ktab(ktab, type = c("F", "Q")) # check function description

  # Estimate optimal number of groups
  gap <- clusGap(
    x = as.matrix(dist_mat),
    FUN = mycluster_hc,
    K.max = k_max,
    B = 100
  )
  optimal_nog_gap <- maxSE(gap$Tab[, "gap"],
    gap$Tab[, "SE.sim"],
    method = "Tibs2001SEmax"
  )
  optimal_nog_silh <- NbClust(
    diss = dist_mat,
    distance = NULL,
    min.nc = 2,
    max.nc = k_max,
    method = "ward.D2",
    index = "silhouette"
  )

  # Actual clustering
  # Creation of distance matrix & hierarchical cluster analysis
  hc_taxa <- hclust(dist_mat, method = "ward.D2")
  # Get labels of dendrogram
  dend_label <- hc_taxa %>%
    as.dendrogram() %>%
    labels()

  # save to list
  results_cl <- list(
    "distance_matrix" = dist_mat,
    "hc_wardD2" = hc_taxa,
    "labels_dendrogram" = dend_label,
    "gap_statistic" = optimal_nog_gap,
    "silhouette" = optimal_nog_silh$Best.nc
  )
}

## Calculate mean trait profiles ----
calc_mean_tps <- function(x,
                          taxa_id) {
  x <- melt(
    x,
    id.vars = c(taxa_id,
                "group",
                "n_taxa_group"),
    value.name = "affinity",
    variable.name = "trait"
  )
  mean_tps <- x[, .(mean_affinity = mean(affinity),
                    n_taxa_group),
                by = c("group", "trait")]
  mean_tps
  mean_tps <- unique(mean_tps)
  mean_tps[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]
}

## Community weighted mean trait ----
# TODO: transform abundances?
calc_cwm <- function(abund,
                     trait,
                     trait_names) {
  taxa_names <- names(abund)[names(abund) != c("site")]
  cwm <- list()
  for (i in trait_names) {
    trait_sub <- trait[, .SD, .SDcols = c("taxon", i)]
    trait_sub <- trait_sub[match(taxa_names, taxon), ]
    
    cwm[[i]] <-
      abund[, apply(.SD, 1, function(x)
        weighted.mean(trait_sub[[i]],
                      w = x)),
        .SDcols = !c("site"),
        by = "site"]
  }
  cwm
}
# Test Code:
# taxa_names <- names(abund_wf$California)[names(abund_wf$California) != "site"]
# test_herbivore <-
#   trait_matrix[Region == "California" &
#                  taxon %in% names(abund_wf$California), .(taxon, feed_herbivore)]
# test_herbivore <- test_herbivore[match(taxa_names, taxon) ]
# abund_wf$California[, apply(.SD,
#                             1,
#                             function(x)
#                               weighted.mean(test_herbivore$feed_herbivore, w = x)),
#                     .SDcols = !c("site"),
#                     by = "site"]


## Community weighted sum trait ----
# TODO: transform abundances?
cws <- function(abund,
         trait,
         trait_names) {
  taxa_names <- names(abund)[names(abund) != c("site")]
  cws <- list()
  for (i in trait_names) {
    trait_sub <- trait[, .SD, .SDcols = c("taxon", i)]
    trait_sub <- trait_sub[match(taxa_names, taxon), ]

    cws[[i]] <-
      abund[, apply(.SD, 1, function(x) {
        sum(x * trait_sub[[i]])
      }),
      .SDcols = !c("site"),
      by = "site"
      ]
  }
  cws
}
# Test code
# taxa_names <- names(abund_wf$California)[names(abund_wf$California) != c("site")]
# trait_sub <- trait_matrix_ls$California[match(taxa_names, taxon), ]
# all(names(abund_wf$California)[2:201] == trait_sub$taxon)

# abund_wf$California[1:2, apply(
#   .SD, 1,
#   function(x) sum(x * trait_sub[["feed_filter"]])
# ), .SDcols = !c("site")]


## XGBOOST ----

# XGBOOST with train/val/test framework
# split ratio for train/test split
# features are the variables used
# threads is the number of CPUs used for parallelization
# id for naming each task
# ?xgboost()
perform_xgboost <- function(x,
                            split_ratio = 0.8,
                            features,
                            id,
                            threads = 4) {
  set.seed(1234)

  # Split in train and test data
  ind <- sample(1:nrow(x), size = round(nrow(x) * split_ratio))
  train <- x[ind, .SD, .SDcols = c(
    features,
    "max_log_tu"
  )]
  test <- x[-ind, .SD, .SDcols = c(
    features,
    "max_log_tu"
  )]

  # Create tasks
  task_train <- TaskRegr$new(
    id = paste0(id, "_train"),
    backend = train,
    target = "max_log_tu"
  )

  # Learner
  xgboost_learner <- lrn(
    "regr.xgboost",
    objective = "reg:squarederror",
    # eval_metric = "rmse",
    # evaluation for validation data set
    gamma = to_tune(p_dbl(lower = 0, upper = 10)),
    # Minimum loss reduction required to make a
    # further partition on a leaf node.
    # The larger, the more conservative the algorithm will be
    # lambda = to_tune(p_dbl(lower = 0, upper = 10)),
    # shrinks features without removing them, L2 regularization
    # the larger the more conservative
    booster = "gbtree",
    eta = to_tune(p_dbl(lower = 0.01, upper = 1)),
    # controlling the learning rate to prevent overfitting
    # (scaling of the contrib. of each tree by a factor 0-1 when added)
    # low values means more nrounds
    # (robust to overfitting but also slower to compute)
    max_depth = to_tune(p_int(lower = 2L, upper = 10L)),
    subsample = to_tune(p_dbl(lower = 0.1, upper = 1)),
    nrounds = to_tune(p_int(lower = 10, upper = 100)),
    # number of boosting rounds
    predict_type = "response"
  )

  # Parallelization
  set_threads(xgboost_learner, n = threads)

  # Tuning
  instance <- tune(
    method = tnr("irace"),
    task = task_train,
    learner = xgboost_learner,
    resampling = rsmp("cv", folds = 3),
    measures = msr("regr.rmse"), # evaluation performance of training data
    terminator = trm("evals", n_evals = 1000)
  )

  # Train xgboost on train dataset with optimized parameters
  # Check model as well on test dataset
  xgboost_tuned <- lrn("regr.xgboost", id = "Xgboost tuned")
  xgboost_tuned$param_set$values <- instance$result_learner_param_vals
  xgboost_tuned$train(task_train)
  pred_train <- xgboost_tuned$predict_newdata(newdata = train)$score(mlr3::msr("regr.rmse"))
  pred_test <- xgboost_tuned$predict_newdata(newdata = test)$score(mlr3::msr("regr.rmse"))

  # Retrieve most imp. variable
  most_imp_vars <- xgboost_tuned$importance()

  # Return
  list(
    "pred_train" = pred_train,
    "pred_test" = pred_test,
    "importance" = most_imp_vars,
    "instance" = instance,
    "final_model" = xgboost_tuned
  )
}

# XGBOOST with nested resampling to evaluate model performance
# TODO: 
# - Generalize for different input data
# - Don't use a loop
# perform_xgboost_nestedresamp <- function(x) {
#   for (region in names(x)) {
#     set.seed(1234)

#     # Data
#     dat <- x[[region]]
#     dat <-
#       dcast(dat, site + max_log_tu ~ trait, value.var = "cwm_val")

#     # Task
#     task <- TaskRegr$new(
#       id = region,
#       backend = dat[, .SD, .SDcols = c(
#         trait_names,
#         "max_log_tu"
#       )],
#       target = "max_log_tu"
#     )
#     # Learner
#     xgboost_learner <- lrn(
#       "regr.xgboost",
#       objective = "reg:squarederror",
#       # for the learning task
#       # eval_metric = "rmse",
#       # evaluation for test data set
#       gamma = to_tune(p_dbl(lower = 0, upper = 10)),
#       # Minimum loss reduction required to make a further partition on a leaf node.
#       # The larger, the more conservative the algorithm will be
#       booster = "gbtree",
#       eta = to_tune(p_dbl(lower = 0.001, upper = 0.3)),
#       # controlling the learning rate (scaling of the contrib. of each tree by a factor 0-1 when added)
#       # to prevent overfitting, low values means more nrounds (robust to overfitting but also slower to compute)
#       max_depth = to_tune(p_int(lower = 2L, upper = 15L)),
#       subsample = to_tune(p_dbl(lower = 0.1, upper = 1)),
#       nrounds = to_tune(p_int(lower = 10, upper = 100)),
#       # number of boosting rounds
#       predict_type = "response"
#     )

#     # Tuning
#     instance <- tune(
#       method = tnr("irace"),
#       task = task,
#       learner = xgboost_learner,
#       resampling = rsmp("cv", folds = 3),
#       measures = msr("regr.rmse"), # evaluation performance of training data
#       terminator = trm("evals", n_evals = 1000)
#     )
#     archive_ls[[region]] <- instance$archive

#     # Final Model
#     xgboost_tuned <- lrn("regr.xgboost", id = "Xgboost tuned")
#     xgboost_tuned$param_set$values <- instance$result_learner_param_vals
#     xgboost_tuned$train(task)
#     model_ls[[region]] <- xgboost_tuned$model

#     # TODO: Model performance with nested resampling
#     # We could use nested resampling to evaluate the model performance
#     # Nested resampling is an additional step after fitting the final model and should not
#     # be used to find optimal hyperparameters.
#     # Idea: If the same data for model selection and evaluation of the model is used, we bias the
#     # performance estimate (eval. on test data could leak information about the test data's structure into
#     # the model)
#     # Steps:
#     # - outer resampling: cv to get different testing an training datasets
#     # - inner resampling: Within the training data use cv to get different inner testing and training data sets
#     # - Tuning the hyperparameters with the inner data splits
#     # - Fits learner on the outer training data set with the tuned hyperparameters from the inner resampling
#     # - Evaluates the performance of the learner on the outer testing data
#     # - Repeat inner resampling + fitting + evaluation for all folds of the outer resampling
#     # - Aggregate perfromance values to get an unbiased performance estimate
#     # Could also use a auto tuner with a different number of cv
#     # at <- auto_tuner(
#     #   method = tnr("grid_search", resolution = 20),
#     #   learner = xgboost_learner,
#     #   resampling = rsmp("cv", folds = 3),
#     #   # for inner resampling
#     #   measure = msr("regr.rmse"),
#     #   terminator = trm("none"),
#     # )

#     # outer_resampling <- rsmp("cv", folds = 3)
#     # rr <- resample(task,
#     #                at,
#     #                outer_resampling,
#     #                store_models = TRUE) # investigate inner tuning with store models set to TRUE
#     # rr$predictions()

#     # perfor_ls[[region]] <- list(
#     #   "perform_inner_rsmp" = extract_inner_tuning_results(rr),
#     #   # Performance results on the inner resampling
#     #   # These values should not be used to fit a final model
#     #   "perform_full_archive" = extract_inner_tuning_archives(rr),
#     #   # Full tuning archive, check if hyperparameters converge?
#     #   "perform_outer_resamp" = rr$score(),
#     #   # Outer resampling performance results
#     #   # Only problematic if there's a huge difference between outer and inner resampling
#     #   "agg_perform_outer_rsmp" = rr$aggregate()
#     #   # Aggregated performance of the outer resampling should be reported (unbiased with optimal
#     #   # hyperparameters)
#     # )

#     # # Interpret model
#     # xgboost_exp <- explain_mlr3(
#     #   xgboost_tuned,
#     #   data     = dat[, .SD, .SDcols = trait_names],
#     #   y        = dat$max_log_tu,
#     #   label    = "XGBOOST",
#     #   colorize = FALSE
#     # )

#     # # Variable importance
#     # # Uses mean dropout loss for xgboost?
#     # imp_ls[[region]] <-
#     #   list(
#     #     "permutation_imp" = model_parts(xgboost_exp,
#     #                                     B = 50,
#     #                                     loss_function = loss_root_mean_square),
#     #     "model_imp" = xgboost_tuned$importance()
#     #   )

#     # PDPs
#     # pdp_ls[[region]] <- model_profile(xgboost_exp)$agr_profiles
#   }
# }

# Bootstrapping to estimate confidence intervals for the 
# prediction and test errors
bootstrap_ci_rmse <- function(
    x,
    model,
    split_ratio = 0.8,
    features,
    n_bootstraps = 1000) {
  # split in training and test data
  # same seed that was used during model training
  # to obtain the same training and test data
  set.seed(1234)
  ind <- sample(1:nrow(x), size = round(nrow(x) * split_ratio))
  train <- x[ind, .SD, .SDcols = names(x) %in% c(
    features,
    "max_log_tu"
  )]
  test <- x[-ind, .SD, .SDcols = names(x) %in% c(
    features,
    "max_log_tu"
  )]

  # bootestrapped CIs
  train_rmse_samples <- numeric(n_bootstraps)
  test_rmse_samples <- numeric(n_bootstraps)
  for (i in 1:n_bootstraps) {
    train_sample <- train[sample(.N, replace = TRUE), ]
    test_sample <- test[sample(.N, replace = TRUE), ]
    train_rmse_samples[[i]] <- model$predict_newdata(newdata = train_sample)$score(mlr3::msr("regr.rmse"))
    test_rmse_samples[[i]] <- model$predict_newdata(newdata = test_sample)$score(mlr3::msr("regr.rmse"))
  }
  list(
    "CI_rmse_train" = quantile(train_rmse_samples, c(0.025, 0.975)),
    "CI_rmse_test" = quantile(test_rmse_samples, c(0.025, 0.975))
  )
}

## GAMs ----
# Function to extract EDF and P value of smooth terms
extractEDF <- function(gam_object, idcol) {
  edf_ls <- list()
  for (i in names(gam_object)) {
    edf_ <- broom::tidy(gam_object[[i]], parametric = FALSE)
    edf_ls[[i]] <- data.table(
      edf = edf_$edf,
      p_smooth_terms = edf_$p.value,
      gam_or_lm = ifelse(edf_$edf > 1.01 &
                           edf_$p.value <= 0.05, "gam", "lm")
    )
  }
  edf_dt <- rbindlist(edf_ls, idcol = idcol)
  edf_dt
}


# General helper functions ----

## Google search ----
# For google search from Rstudio
query_google <- function(x) {
  invisible(lapply(x,
                   function(y) {
                     utils::browseURL(url = paste0("https://google.com/search?q=", y))
                   }))
}

## Load batches of similarly named datasets ----
# So far, only RDS files
load_data <- function(path,
                      pattern,
                      format = "rds",
                      name_rm_pattern = NULL) {
  files <- list.files(path = path, pattern = pattern)

  if (format == "rds") {
    data <- lapply(files, function(y) readRDS(file = file.path(path, y)))
    name <- sub("\\.rds", "", files)
  }

  if (format == "csv") {
    data <- lapply(files, function(y) fread(file = file.path(path, y)))
    name <- sub("\\.csv", "", files)
  }

  if (!is.null(name_rm_pattern)) {
    name <- sub(name_rm_pattern, "", name)
  }
  data <- setNames(data, name)
  data
}

## Removing strange notations ----
# E.g., in taxonomic nomenclature
rm_strange_notations <- function(x){
  sub("\\?| group", "", x)
}

# Plotting ----
## Dendrograms ----
fun_dendrog_pl <- function(hc,
                           optimal_nog,
                           labels,
                           hang_height = 0.001) {
  hc %>%
    as.dendrogram() %>%
    color_branches(k = optimal_nog) %>%
    hang.dendrogram(hang_height = hang_height) %>%
    set("labels_cex", 0.8) %>%
    dendextend::ladderize() %>%
    set("labels", labels)
}