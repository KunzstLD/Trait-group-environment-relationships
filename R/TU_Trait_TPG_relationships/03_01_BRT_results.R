## Performance results  ----
# Train/val/test split approach with regularization
res_xgboost_comb <- load_data(
  path = path_cache,
  pattern = "res\\_xgboost.*",
  name_rm_pattern = "res\\_xgboost\\_"
)

# Get prediction on test and training data
xgboost_perform <- purrr::map(res_xgboost_comb, ~ sapply(.x, function(x) {
  c(
    "train" = x$pred_train,
    "test" = x$pred_test
  )
}))
xgboost_perform <- lapply(xgboost_perform, function(x) as.data.table(x, keep.rownames = TRUE)) %>%
    rbindlist(., id = "method") %>%
    melt(., id.vars = c("rn", "method")) %>%
    dcast(., ... ~ rn, value.var = "value")
setnames(xgboost_perform, "variable", "region")
# saveRDS(xgboost_perform, file.path(path_cache, "xgboost_perform.rds"))

# Plotting test and training RMSE
xgboost_perform <- readRDS(file.path(path_cache, "xgboost_perform.rds"))
xgboost_perform %>%
  melt(.,
    id.vars = c("region", "method"),
    value.name = "rmse"
  ) %>%
  ggplot(., aes(x = method, y = rmse)) +
  geom_col(aes(fill = as.factor(variable)),
    position = "dodge",
    width = 0.6
  ) +
  scale_fill_d3(labels = c("Test", "Train")) +
  coord_flip() +
  labs(
    x = "",
    y = "RMSE",
    fill = ""
  ) +
  facet_wrap(. ~ region) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 14
    )
  )
ggsave(
  filename = file.path(
    path_out,
    "Graphs",
    "perform_xgboost.png"
  ),
  width = 35,
  height = 20,
  units = "cm"
)

# Most important traits ----
# Use impurity importance
xgboost_imp <- purrr::map(res_xgboost_comb, ~ sapply(.x, function(x) {
  "train" <- x$importance
}))

# For some datasets not all variables have importance scores
# There importance score is 0
xgboost_imp <- lapply(xgboost_imp, function(x) {
    lapply(x, function(y) as.data.table(y, keep.rownames = TRUE))
}) %>%
    lapply(., function(x) rbindlist(x, id = "region")) %>%
    rbindlist(., id = "method")
setnames(
    xgboost_imp,
    c("rn", "y"),
    c("trait_TPG", "score")
)

# Top 5 per region
# Add grouping features for cwm and cws
xgboost_imp[
  method %in% c("cwm", "cws"),
  grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait_TPG)
]

# cwm:
# predator 4, gatherer 3 , sensitivity organic 3
xgboost_imp[method == "cwm", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = "trait_TPG"] %>% 
  .[order(-N), ]

# Dominance of feeding mode traits (at least one in each region)
# then locomotion 
xgboost_imp[method == "cwm", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = c("grouping_feature", "region")] %>% 
  .[order(grouping_feature, -N), ]

# Plot most important traits
xgboost_imp[method == "cwm", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>%
  ggplot(
    .,
    aes(
      x = as.factor(region),
      y = score,
      group = as.factor(trait_TPG)
    )
  ) +
  geom_col(
    position = "dodge"
  ) +
  geom_text(
    mapping = aes(
      label = as.factor(trait_TPG),
      group = as.factor(trait_TPG)
    ),
    position = position_dodge(width = .9),
    size = 4.1,
    hjust = -0.1
  ) +
  lims(y = c(0, 0.5)) +
  labs(
    x = "Region",
    y = "Impurity importance"
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    legend.title = element_text(
      family = "Roboto Mono",
      size = 16
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    legend.position = "none",
    panel.grid = element_blank()
  )


# tpgs_rel
# T12, T5 in 4 regions
# T2, T10, T1, T8 in 3 regions 
xgboost_imp[method == "tpgs_rel", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = "trait_TPG"] %>% 
  .[order(-N), ] 

# Load defining trait combinations
# TODO: Maybe look deeper into trait profiles for these groups?
defining_traits <- readRDS(file.path(path_cache, "defining_traits.rds"))
defining_traits[, group := paste0("T", group)]
defining_traits[group %in% c("T5", "T12"), ]
defining_traits[group %in% c("T1", "T2", "T8", "T10"), ]

# Taxonomic composition of TPGs
# First glance: groups look like they are mainly composed of
# taxa from one to two orders 
tpg_taxonomic_comp <- readRDS(file.path(path_cache, "tpg_taxonomic_composition.rds"))
tpg_taxonomic_comp[group %in% c("T5", "T12"), ]
