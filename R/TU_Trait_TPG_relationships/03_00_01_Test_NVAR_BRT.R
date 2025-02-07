# __________________________________________________________________________________________________
# Test BRTs with 15 most important CWM traits
# __________________________________________________________________________________________________

## Data preproc ----
# Similar to 03_00
# Community weighted and community sum traits
data_cwm <- readRDS(file.path(path_cache, "data_cwm.rds"))
lapply(data_cwm, function(x) x[trait == "sensitivity_organic", range(cwm_val)])

# Toxicitiy
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

# Few sites have two IDs (for eco and wq)
# because chemical and eco information were taken from slightly
# different positions (mainly because of accessability issues)
data_cwm$Midwest[STAID == "T03611200", STAID := "T03611100"]
data_cwm$PN[site == "T12073525", site := "T12073425"]

# Combine max tu and cwm for each region, but make an exception for Midwest (merge via STAID)
data_cwm <- lapply(data_cwm, function(x) {
  on <- if ("STAID" %in% names(x)) c("STAID" = "site") else "site"
  x[max_tu, max_log_tu := i.max_log_tu, on = on]
})

# Few sites that don't have chemical data
data_cwm$Midwest <- data_cwm$Midwest[!is.na(max_log_tu), ]

# List of importance scores for CWM traits, in decreasing order
traits_expanded <- readRDS(file.path(path_cache, "traits_expanded.rds"))

## BRT models with 15 CWM Traits ----
res_xgboost_cwm_subset_test <- list()
data_cwm_subset_ls <- list()

for (region_name in names(data_cwm)) {
  
  # Test with 15 most important traits (remove the last 5 most unimportant)
  selected_traits <- traits_expanded[region==region_name, .SD[1:15]][["traits"]]
  
  # Create subset and save for later
  data_cwm_subset <- data_cwm[[region_name]][trait %in% selected_traits, ]
  data_cwm_subset_ls[[region_name]] <- data_cwm_subset
  
  # trait_names <- unique(data_cwm_subset$trait)
  # x <- dcast(data_cwm_subset, site + max_log_tu ~ trait,
  #            value.var = "cwm_val"
  # )
  # res_xgboost_cwm_subset_test[[region_name]] <- perform_xgboost(
  #   x = x,
  #   features = trait_names,
  #   id = region_name
  # )
}
saveRDS(res_xgboost_cwm_subset_test, file.path(path_cache, "res_xgboost_cwm_subset_test_subset_test.rds"))
# lapply(res_xgboost_cwm_subset_test, function(x) c(x$pred_train^2, x$pred_test^2))
# lapply(res_xgboost_cwm_subset_test, function(x) x$importance)
# lapply(res_xgboost_cwm_subset_test, function(x) x$instance)

# Get prediction on test and training data
# lapply(res_xgboost_cwm_subset_test, function(x) data.table("MSE_train"=x$pred_train^2,
#                                                "MSE_test"=x$pred_test^2)) |> 
#   rbindlist(id="Region") |> 
#   saveRDS(xgboost_perform, file.path(path_cache, "xgboost_perform.rds"))

## Bootstrap for CI ----

# Extract final models
res_xgboost_cwm_subset_test <- readRDS(file.path(path_cache, "res_xgboost_cwm_subset_test_subset_test.rds"))
final_model_cwm <- lapply(res_xgboost_cwm_subset_test, `[[`, "final_model")

# Convert the datasets to wide format
data_cwm_subset_ls <- lapply(data_cwm_subset_ls, function(x) {
  dcast(x, site + max_log_tu ~ trait,
        value.var = "cwm_val"
  )
})

# Loop over datasets and models to calculate CIs
names(data_cwm_subset_ls) == names(final_model_cwm)

cwm_ci_rmse <- Map(
  function(x, y) {
    bootstrap_ci_rmse(
      x = x,
      model = y,
      features = names(x)[!names(x) %in% c("site", "max_log_tu")]
    )
  },
  x = data_cwm_subset_ls,
  y = final_model_cwm
)
# saveRDS(cwm_ci_rmse, file.path(path_cache, "cwm_ci_rmse_subset_test.rds"))
cwm_ci_rmse <- readRDS(file.path(path_cache, "cwm_ci_rmse_subset_test.rds"))

## Create output table with mean squared error and range ----
# Add CI 
perform_tbl <- lapply(res_xgboost_cwm_subset_test, function(x)
  data.table(
    "train_rmse" = x$pred_train,
    "test_rmse" = x$pred_test
  )) |>
  rbindlist(id = "region") |>
  _[, .(
    region,
    train_mse = train_rmse ^ 2,
    test_mse = test_rmse ^ 2
  )] 

ci_table <- lapply(cwm_ci_rmse, function(x)
  as.data.table(unlist(x), keep.rownames = TRUE)) |>
  rbindlist(id = "region") |>
  _[, .(region, V1, V2 = V2 ^ 2)] |> # Squre to obtain the MSE
  dcast(region ~ V1, value.var = "V2")

perform_tbl <- merge.data.table(perform_tbl, ci_table, by="region")
setnames(
  perform_tbl,
  names(perform_tbl),
  c("region",
    "pred_train",
    "pred_test",
    "CI_mse_test_2.5",
    "CI_mse_test_97.5",
    "CI_mse_train_2.5",
    "CI_mse_train_97.5")
)
fwrite(
  perform_tbl,
  file.path(path_paper, "Tables", "xgboost_subset_test_perfom.csv")
)

## Plotting ----
# Create plot of CWM with 20 (pred_error) and 15 traits (perform_tbl) next to each other!
pred_error <- readRDS(file.path(path_cache, "pred_error.rds"))
pred_error[, c("CI_mse_train_2.5",
               "CI_mse_train_97.5",
               "CI_mse_test_2.5",
               "CI_mse_test_97.5",
               "pred_train",
               "pred_test") := lapply(.SD, function(x)
                 x ^ 2), .SDcols = c(
                   "CI_rmse_train_2.5",
                   "CI_rmse_train_97.5",
                   "CI_rmse_test_2.5",
                   "CI_rmse_test_97.5",
                   "pred_train",
                   "pred_test"
                 )]
pred_error_cwm <- pred_error[approach=="cwm", .SD, .SDcols = c("region",
                                                               "approach",
                                                               "pred_train",
                                                               "pred_test",
                                                               "CI_mse_train_2.5",
                                                               "CI_mse_train_97.5",
                                                               "CI_mse_test_2.5",
                                                               "CI_mse_test_97.5")]

perform_tbl <- fread(file.path(path_paper, "Tables", "xgboost_subset_test_perfom.csv"))
perform_tbl[, approach := "cwm_15"]
perform_tbl[region=="PN", region := "Northwest"]

pred_error_cwm <- rbind(pred_error_cwm, perform_tbl)


ggplot(pred_error_cwm) +
  geom_point(aes(x = approach, y = pred_test)) +
  geom_point(aes(x = approach, y = pred_train),
             alpha = 0.3,
             position = position_nudge(x = 0.2)
  ) +
  geom_pointrange(
    aes(
      x = approach,
      y = pred_test,
      ymin = CI_mse_test_2.5,
      ymax = CI_mse_test_97.5
    )
  ) +
  geom_pointrange(
    aes(
      x = approach,
      y = pred_train,
      ymin = CI_mse_train_2.5,
      ymax = CI_mse_train_97.5
    ),
    position = position_nudge(x = 0.2),
    alpha = 0.35
  ) +
  facet_wrap(. ~ region) +
  labs(x = "", y = "MSE") +
  scale_x_discrete(labels = c("20 CWM \n traits", "15 most important \n CWM traits")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16, face="bold"),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 16
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 16
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 16
    )
  )
ggsave(
  filename = file.path(
    path_paper,
    "Graphs",
    "perform_xgboost_20_15_CWM.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)
