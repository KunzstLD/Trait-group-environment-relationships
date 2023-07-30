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

## CWM, CWS & TU data ----

# Community weighted and community sum traits
data_cwm <- readRDS(file.path(path_cache, "data_cwm.rds"))
data_cws <- readRDS(file.path(path_cache, "data_cws.rds"))
lapply(data_cwm, function(x) x[trait == "sensitivity_organic", range(cwm_val)])

# Toxicitiy
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")
# max_tu[max_log_tu == -5, .N, by = "Region"]

# Few sites have two IDs (for eco and wq)
# because chemical and eco information were taken from slightly
# different positions (mainly because of accessability issues)
# eco site T03611200, corresponding wq site: T03611100
# eco site T12073525, corresponding wq site: T12073425
data_cwm$Midwest[STAID == "T03611200", STAID := "T03611100"]
data_cwm$PN[site == "T12073525", site := "T12073425"]
data_cws$Midwest[STAID == "T03611200", STAID := "T03611100"]
data_cws$PN[site == "T12073525", site := "T12073425"]

# Combine max tu and cwm for each region
# make an exception for Midwest (merge via STAID)
data_cwm <- lapply(data_cwm, function(x) {
  on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
  x[max_tu, max_log_tu := i.max_log_tu, on = on]
})
data_cws <- lapply(data_cws, function(x) {
  on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
  x[max_tu, max_log_tu := i.max_log_tu, on = on]
})

# Few sites that don't have chemical data
# lapply(data_cwm, function(y)
#   y[is.na(max_log_tu), ])
# Midwest:
# unique(data_cwm$Midwest[is.na(max_log_tu), STAID])
# T03318800, T393247089260701
data_cwm$Midwest <- data_cwm$Midwest[!is.na(max_log_tu), ]
data_cws$Midwest <- data_cws$Midwest[!is.na(max_log_tu), ]
# saveRDS(data_cwm, file.path(path_cache, "data_cwm_final.rds"))

# 432 sites are used in the end (chemical information)
# Stream insect abundance data for 435 sites, 
# but three from Midwest have the same chemical information 
# lapply(data_cwm, function(x) {
#   if ("STAID" %in% names(x)) x[, uniqueN(STAID)] else x[, uniqueN(site)]
# }) %>% unlist %>% sum

## TPG data & TU data ----
trait_groups_rel <- readRDS(file.path(path_cache, "trait_groups_rel.rds"))
# trait_groups <- readRDS(file.path(path_cache, "trait_groups.rds")) 

# Final preprocessing and add max logTU
trait_groups_rel <- lapply(trait_groups_rel, function(x) {

  # RM sites with no chemical information
  x$Midwest <- x$Midwest[!STAID %in% c(
    "T03318800",
    "T393247089260701"
  ), ]

  # Change the two site ids to match the chemical information
  x$Midwest[STAID == "T03611200", STAID := "T03611100"]
  x$PN[site == "T12073525", site := "T12073425"]

  # Add max TU
  x <- lapply(x, function(y) {
    on <- if ("STAID" %in% names(y)) c("STAID" = "site") else "site"
    y[max_tu, max_log_tu := i.max_log_tu, on = on]
  })
})
c(trait_groups_rel_family, trait_groups_rel_genus) %<-% trait_groups_rel[c(
  "family_lvl",
  "genus_lvl"
)]
# saveRDS(trait_groups_rel, file.path(path_cache, "trait_groups_rel_final.rds"))

# Boosting ----

## Train/validation/test split approach ----

## CWM ----
trait_names <- unique(data_cwm$California$trait)
res_xgboost_cwm <- list()

for (region in names(data_cwm)) {
  x <- dcast(data_cwm[[region]], site + max_log_tu ~ trait,
    value.var = "cwm_val"
  )
  res_xgboost_cwm[[region]] <- perform_xgboost(
    x = x,
    features = trait_names,
    id = region
  )
}
saveRDS(res_xgboost_cwm, file.path(path_cache, "res_xgboost_cwm.rds"))
# should also save the tuned models!
lapply(res_xgboost_cwm, function(x) c(x$pred_train^2, x$pred_test^2))
lapply(res_xgboost_cwm, function(x) x$importance)
lapply(res_xgboost_cwm, function(x) x$instance)

# Check predictions for test data
# set.seed(1234)
# midwest <- dcast(data_cwm$Midwest, site + max_log_tu ~ trait,
#   value.var = "cwm_val"
# )
# ind <- sample(1:nrow(midwest), size = round(nrow(midwest) * 0.8))
# train_midwest <- midwest[ind, .SD, .SDcols = c(
#   trait_names,
#   "max_log_tu"
# )]
# test <- midwest[-ind, .SD, .SDcols = c(
#   trait_names,
#   "max_log_tu"
# )]

# pred_train_midwest <- res_xgboost_cwm$Midwest$final_model$predict_newdata(newdata = train_midwest)
# as.data.table(pred_train_midwest) %>% 
# .[, abs_dev := abs(truth - response)] %>% 
# ggplot(., aes(x = truth, y = response)) +
# geom_point()


## CWS ----
res_xgboost_cws <- list()
for (region in names(data_cws)) {
  x <- dcast(data_cws[[region]], site + max_log_tu ~ trait,
    value.var = "cws_val"
  )
  res_xgboost_cws[[region]] <- perform_xgboost(
    x = x,
    features = trait_names,
    id = region
  )
}
# saveRDS(res_xgboost_cws, file.path(path_cache, "res_xgboost_cws.rds"))

## TPGs relative fraction ----
res_xgboost_tpgs_rel_family <- list()
for (region in names(trait_groups_rel_family)) {
  var_names <- names(trait_groups_rel_family[[region]])[names(trait_groups_rel_family[[region]]) %like% "^T[0-9]"]
  res_xgboost_tpgs_rel_family[[region]] <- perform_xgboost(
    x = trait_groups_rel_family[[region]],
    features = var_names,
    id = region
  )
}
lapply(res_xgboost_tpgs_rel_family, function(x) c(x$pred_train^2, x$pred_test^2))
lapply(res_xgboost_tpgs_rel_family, function(x) x$importance)
lapply(res_xgboost_tpgs_rel_family, function(x) x$instance)
saveRDS(res_xgboost_tpgs_rel_family, file.path(path_cache, "res_xgboost_tpgs_rel_family.rds"))

res_xgboost_tpgs_rel_genus <- list()
for (region in names(trait_groups_rel_genus)) {
  var_names <- names(trait_groups_rel_genus[[region]])[names(trait_groups_rel_genus[[region]]) %like% "^T[0-9]"]
  res_xgboost_tpgs_rel_genus[[region]] <- perform_xgboost(
    x = trait_groups_rel_genus[[region]],
    features = var_names,
    id = region
  )
}
saveRDS(res_xgboost_tpgs_rel_genus, file.path(path_cache, "res_xgboost_tpgs_rel_genus.rds"))
lapply(res_xgboost_tpgs_rel_genus, function(x) c(x$pred_train^2, x$pred_test^2))
lapply(res_xgboost_tpgs_rel_genus, function(x) x$importance)
lapply(res_xgboost_tpgs_rel_genus, function(x) x$instance)

# ## TPGs absolute fraction ----
# res_xgboost_tpgs <- list()
# for (region in names(trait_groups)) {
#   var_names <- names(trait_groups[[region]])[names(trait_groups[[region]]) %like% "^T[0-9]"]
#   res_xgboost_tpgs[[region]] <- perform_xgboost(
#     x = trait_groups[[region]],
#     features = var_names,
#     id = region
#   )
# }
# saveRDS(res_xgboost_tpgs, file.path(path_cache, "res_xgboost_tpgs.rds"))


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

