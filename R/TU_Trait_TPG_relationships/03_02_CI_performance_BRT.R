# _______________________________________________________
# Confidence intervals for MSE estimates 
# _______________________________________________________

# CWM ----
# Load & transform CWM data
data_cwm <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
trait_names <- unique(data_cwm$California$trait)
data_cwm <- lapply(data_cwm, function(x) {
  dcast(x, site + max_log_tu ~ trait,
    value.var = "cwm_val"
  )
})

# Load BRT output
res_xgboost_cwm <- readRDS(file.path(path_cache, "res_xgboost_cwm.rds"))

# Extract final models
final_model_cwm <- lapply(res_xgboost_cwm, `[[`, "final_model")

# Loop over datasets and models to calculate CIs
names(data_cwm) == names(final_model_cwm)
cwm_ci_rmse <- Map(
  function(x, y) {
    bootstrap_ci_rmse(
      x = x,
      model = y,
      features = trait_names
    )
  },
  x = data_cwm,
  y = final_model_cwm
)
saveRDS(cwm_ci_rmse, file.path(path_cache, "cwm_ci_rmse.rds"))

# TPGs ----
## Family level ----
trait_groups_rel <- readRDS(file.path(path_cache, "trait_groups_rel_final.rds"))
trait_groups_rel_family <- trait_groups_rel$family_lvl
tpg_names <- names(trait_groups_rel_family$Southeast)[
  !names(trait_groups_rel_family$Southeast) %in% c("max_log_tu", "site")
]

# Load BRT output
res_xgboost_tpgs_rel_family <- readRDS(
  file.path(path_cache, "res_xgboost_tpgs_rel_family.rds")
)

# Extract final models
final_model_tpg_family <- lapply(
  res_xgboost_tpgs_rel_family,
  `[[`, "final_model"
)

# Loop over datasets and models to calculate CIs
names(trait_groups_rel_family) == names(final_model_tpg_family)
tpg_family_ci_rmse <- Map(
  function(x, y) {
    bootstrap_ci_rmse(
      x = x,
      model = y,
      features = tpg_names
    )
  },
  x = trait_groups_rel_family,
  y = final_model_tpg_family
)
saveRDS(tpg_family_ci_rmse, file.path(path_cache, "tpg_family_ci_rmse.rds"))

# Genus level ----
trait_groups_rel_genus <- trait_groups_rel$genus_lvl

# Load BRT output
res_xgboost_tpgs_rel_genus <- readRDS(
  file.path(path_cache, "res_xgboost_tpgs_rel_genus.rds")
)

# Extract final models
final_model_tpg_genus <- lapply(
  res_xgboost_tpgs_rel_genus,
  `[[`, "final_model"
)

# Loop over datasets and models to calculate CIs
names(trait_groups_rel_genus) == names(final_model_tpg_genus)
tpg_genus_ci_rmse <- Map(
  function(x, y) {
    bootstrap_ci_rmse(
      x = x,
      model = y,
      features = tpg_names
    )
  },
  x = trait_groups_rel_genus,
  y = final_model_tpg_genus
)
saveRDS(tpg_genus_ci_rmse, file.path(path_cache, "tpg_genus_ci_rmse.rds"))


# Plotting ----
# Prediction errors
pred_error_cwm <- lapply(res_xgboost_cwm, `[`, c("pred_train", "pred_test")) %>%
  rbindlist(., idcol = "region") %>%
  .[, approach := "cwm"]
pred_error_tpgs_family <- lapply(res_xgboost_tpgs_rel_family, `[`, c("pred_train", "pred_test")) %>%
  rbindlist(., idcol = "region") %>%
  .[, approach := "tpg_family"]
pred_error_tpgs_genus <- lapply(res_xgboost_tpgs_rel_genus, `[`, c("pred_train", "pred_test")) %>%
  rbindlist(., idcol = "region") %>%
  .[, approach := "tpg_genus"]
pred_error <- rbind(
  pred_error_cwm,
  pred_error_tpgs_family,
  pred_error_tpgs_genus
)
pred_error[region == "PN", region := "Northwest"]

# CIs
cwm_ci_rmse <- readRDS(file.path(path_cache, "cwm_ci_rmse.rds"))
tpg_genus_ci_rmse <- readRDS(file.path(path_cache, "tpg_genus_ci_rmse.rds"))
tpg_family_ci_rmse <- readRDS(file.path(path_cache, "tpg_family_ci_rmse.rds"))

# Helper function to transform CI results into a list of data.tables
ci_result_to_dt <- function(...) {
  result_list <- list()

  for (lt in list(...)) {
    output <- lapply(lt, as.data.frame) %>%
      lapply(., function(x) as.data.table(x, keep.rownames = TRUE)) %>%
      rbindlist(., idcol = "region")
    result_list <- c(result_list, list(output))
  }

  return(result_list)
}
ci_pred_error <- ci_result_to_dt(cwm_ci_rmse, tpg_family_ci_rmse, tpg_genus_ci_rmse)
names(ci_pred_error) <- c("cwm", "tpg_family", "tpg_genus")
ci_pred_error <- rbindlist(ci_pred_error, idcol = "approach")
setnames(ci_pred_error, "rn", "quantile")
ci_pred_error[region == "PN", region := "Northwest"]
ci_pred_error <- dcast(ci_pred_error, ... ~ quantile,
  value.var = c("CI_rmse_train", "CI_rmse_test")
)

# Create output table
# fwrite(ci_pred_error[, .(
#   approach,
#   region,
#   `CI_mse_train_2.5%` = round(`CI_rmse_train_2.5%` ^
#                                  2, digits = 2),
#   `CI_mse_train_97.5%` = round(`CI_rmse_train_97.5%` ^
#                                   2, digits = 2),
#   `CI_mse_test_2.5%` = round(`CI_rmse_test_2.5%` ^ 2, digits =
#                                 2),
#   `CI_mse_test_97.5` = round(`CI_rmse_test_97.5%` ^ 2, digits =
#                                 2)
# )],
# file.path(path_paper, "Tables", "ci_pred_error.csv"))

# Create paper plot
pred_error[ci_pred_error,
  `:=`(
    CI_rmse_train_2.5 = `i.CI_rmse_train_2.5%`,
    CI_rmse_train_97.5 = `i.CI_rmse_train_97.5%`,
    CI_rmse_test_2.5 = `i.CI_rmse_test_2.5%`,
    CI_rmse_test_97.5 = `i.CI_rmse_test_97.5%`
  ),
  on = c("region", "approach")
]
saveRDS(pred_error, file.path(path_cache, "pred_error.rds"))

ggplot(pred_error) +
  geom_point(aes(x = approach, y = pred_test^2)) +
  geom_point(aes(x = approach, y = pred_train^2),
    alpha = 0.3,
    position = position_nudge(x = 0.2)
  ) +
  geom_pointrange(
    aes(
      x = approach,
      y = pred_test^2,
      ymin = CI_rmse_test_2.5^2,
      ymax = CI_rmse_test_97.5^2
    )
  ) +
  geom_pointrange(
    aes(
      x = approach,
      y = pred_train^2,
      ymin = CI_rmse_train_2.5^2,
      ymax = CI_rmse_train_97.5^2
    ),
    position = position_nudge(x = 0.2),
    alpha = 0.35
  ) +
  facet_wrap(. ~ region) +
  labs(x = "", y = "MSE") +
  scale_x_discrete(labels = c("Traits", "TPG \n (family)", "TPG \n (genus)")) +
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
    "perform_xgboost.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)
