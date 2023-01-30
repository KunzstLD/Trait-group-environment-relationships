# ________________________________________________________
# Boosted Regression trees
# Resources for mlr3 with xgboost:
# https://github.com/KI-Research-Institute/LearningWithExternalStats/blob/f6ead01200c2a85ac0d8e5c67f0d99a7061eabb2/extras/mlWrappers/wxgboost.R
# https://mlr3book.mlr-org.com/performance.html
# cross validation:
# https://machinelearningmastery.com/k-fold-cross-validation/
# And on boosting parameters (e.g. understand nrounds better)
# Variable importance in boosting
# ________________________________________________________

# Data preprocessing ----
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

# Boosting ----

## Train/validation/test split approach ----
# TODO: Check why importance values are not returned!
trait_names <- unique(data_cwm$California$trait)
res_xgboost <- list()

for (region in names(data_cwm)) {
  x <- dcast(data_cwm[[region]], site + max_log_tu ~ trait,
    value.var = "cwm_val"
  )
  res_xgboost[[region]] <- perform_xgboost(
    x = x,
    features = trait_names,
    id = region
  )
}
lapply(res_xgboost, function(y) c(y$pred_train, y$pred_test)) %>%
  saveRDS(., file.path(path_cache, "performance_xgboost_train_val_test.rds"))

summary_funs <- list(
  "min" = min,
  "max" = max,
  "mean" = mean,
  "median" = median
)

## Nested resampling ----
# perfor_ls$California$perform_full_archive[, lapply(summary_funs, function(f)
#   f(regr.rmse))]
# as.data.table(archive_ls$California)
# saveRDS(perfor_ls, file.path(path_cache, "performance_brt_combined.rds"))
# saveRDS(model_ls, file.path(path_cache, "model_summary_brt_combined.rds"))

# TODO: regr. mean error still rel. high, why?
# Seems to be a problem with the size of our data and the number of variables
# -> If we do CV then get rmse between 2.5 and 3
# Is there another way? If we train on the full dataset then we get
# TODO Nested resampling
# https://arxiv.org/pdf/2107.05847.pdf
# https://www.tidymodels.org/learn/work/nested-resampling/

# Plots permutation importance
# TODO
# - How to deal with correlation among variables
# - Rather use Impurity?
# https://christophm.github.io/interpretable-ml-book/feature-importance.html
# saveRDS(imp_ls, file.path(path_cache, "variable_importance_brt_combined.rds"))


# Plot indicates loss of original data, loss with permuted data, and results for 100 permutations
# plot(imp_ls$California$permutation_imp)
# # ggsave(
# #   filename = file.path(path_out,
# #                        "Graphs",
# #                        "Most_imp_traits_California.png"),
# #   width = 50,
# #   height = 30,
# #   units = "cm"
# # )

# # Overview 5 most important traits  
# permutation_imp <- lapply(imp_ls, function(y) y$permutation_imp)
# most_imp_permutation_imp <- lapply(permutation_imp, function(y)
#   as.data.table(y) %>%
#     .[, .(mean_dropout_loss = mean(dropout_loss)), by = "variable"] %>%
#     .[!variable %in% c("_baseline_", "_full_model_"), ] %>%
#     .[order(-mean_dropout_loss), .SD[1:5, ]])
# most_imp_permutation_imp <- rbindlist(most_imp_permutation_imp, idcol = "region")
# most_imp_permutation_imp[, region_QA := fcase(
#   region == "California",
#   "CSQA",
#   region == "Northeast",
#   "NESQA",
#   region == "PN",
#   "PNSQA",
#   region == "Southeast",
#   "SESQA"
# )] # TODO: this needs to be standardized in all figures and tables
# most_imp_permutation_imp[, mean_dropout_loss := round(mean_dropout_loss, digits = 2)]
# setcolorder(most_imp_permutation_imp,
#             c("region", "region_QA"))
# fwrite(most_imp_permutation_imp,
#        file = file.path(path_paper,
#                         "Tables",
#                         "Most_imp_traits.csv"))

# # PDPs
# saveRDS(pdp_ls, file.path(path_cache, "pdp_brt_combined.rds"))
# plot(pdp_ls$California) +
#   scale_y_continuous("Max Log TU") +
#   theme_bw() +
#   ggtitle("Partial Dependence profile California") +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text.x = element_text(family = "Roboto Mono",
#                                size = 14),
#     axis.text.y = element_text(family = "Roboto Mono",
#                                size = 14),
#     strip.text = element_text(family = "Roboto Mono",
#                               size = 14),
#     legend.position = "none",
#   )
# ggsave(
#   filename = file.path(path_out,
#                        "Graphs",
#                        "PDP_California.png"),
#   width = 50,
#   height = 35,
#   units = "cm"
# )
