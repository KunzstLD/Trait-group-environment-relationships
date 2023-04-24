## Performance results  ----
# Train/val/test split approach with regularization
res_xgboost_comb <- load_data(
  path = path_cache,
  pattern = "res\\_xgboost.*",
  name_rm_pattern = "res\\_xgboost\\_"
)
# names(res_xgboost_comb)

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

# Test and training RMSE
xgboost_perform <- readRDS(file.path(path_cache, "xgboost_perform.rds"))

# Publication table
perform_publ_tbl <- xgboost_perform[
  method %in% c("cwm", "tpgs_rel", "tpgs_rel_genus"),
  .(
    method,
    region,
    train_mse = round(train.regr.rmse^2, digits = 2),
    test_mse = round(test.regr.rmse^2, digits = 2)
  )
] %>%
  .[order(method, region), ] %>%
  melt(., measure.vars = c("train_mse", "test_mse")) %>% 
  dcast(., ...~ method+variable)
setcolorder(
  perform_publ_tbl,
  c(
    "region",
    "cwm_train_mse",
    "tpgs_rel_train_mse",
    "tpgs_rel_genus_train_mse",
    "cwm_test_mse",
    "tpgs_rel_test_mse",
    "tpgs_rel_genus_test_mse"
  )
)
fwrite(
  perform_publ_tbl,
  file.path(path_paper, "Tables", "xgboost_perfom.csv")
)

# Plot
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
traits <- unique(xgboost_imp$trait_TPG)[1:20]
lookup_traits <- data.table(
  trait = traits,
  trait_label = c(
    "predator",
    "burrowing",
    "large size",
    "gatherer",
    "sens. organic",
    "shredder",
    "bi/multivolt.",
    "swimming",
    "plast. & spi.",
    "gills",
    "medium size",
    "parasite",
    "crawling",
    "univolt.",
    "herbivore",
    "semivolt.",
    "filterer",
    "tegument",
    "small size",
    "sessil"
  )
)
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
xgboost_imp[lookup_traits,
  trait_label := i.trait_label, on = c("trait_TPG" = "trait")
]

# CWM:
# predator 4, gatherer 3 , sensitivity organic 3
xgboost_imp[method == "cwm", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = "trait_TPG"] %>% 
  .[order(-N), ]
# xgboost_imp[method == "cwm", ] %>%
#   .[order(-score), .SD[1:5, ], by = "region"]  %>% 
#   saveRDS(., file.path(path_cache, "most_important_traits_cwm.rds"))

# Table for publication
xgboost_imp[method == "cwm", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>%
  .[!is.na(score), .(region, trait_label, score = round(score, digits = 2))] %>%
  fwrite(., file.path(path_paper, "Tables", "Most_imp_traits.csv"))

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

# Abundance weighted fraction (TPGS_REL)
# Family-level approach
# T12, T5 in 4 regions
# T2, T10, T1, T8 in 3 regions 
xgboost_imp[method == "tpgs_rel", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = "trait_TPG"] %>% 
  .[order(-N), ] 

# Table for publication
xgboost_imp[method == "tpgs_rel", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>%
  .[!is.na(score), .(region, trait_TPG, score = round(score, digits = 2))] %>% 
  .[, trait_TPG := paste0(trait_TPG, "_genus")]  %>% 
  fwrite(., file.path(path_paper, "Tables", "Most_imp_TPGs.csv"))

# Genus-level approach
xgboost_imp[method == "tpgs_rel_genus", trait_TPG := paste0(trait_TPG, "_genus")]
xgboost_imp[method == "tpgs_rel_genus", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = "trait_TPG"] %>% 
  .[order(-N), ] 

# Table for publication
xgboost_imp[method == "tpgs_rel_genus", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>%
  .[!is.na(score), .(region, trait_TPG, score = round(score, digits = 2))] %>% 
  fwrite(., file.path(path_paper, "Tables", "Most_imp_TPGs_genus.csv"))

# Load defining trait combinations
defining_traits_family <- readRDS(file.path(path_cache, "defining_traits_family.rds"))
defining_traits_family[, group := paste0("T", group)]
defining_traits_family[group %in% c("T5", "T12"), ]
defining_traits_family[group %in% c("T1", "T2", "T8", "T10"), ]

defining_traits_genus <- readRDS(file.path(path_cache, "defining_traits_genus.rds"))
defining_traits_genus[, group := paste0("T", group)]

# Taxonomic composition of TPGs
# First glance: groups look like they are mainly composed of
# taxa from one to two orders 
tpg_taxonomic_comp_family <- readRDS(file.path(path_cache, "tpg_taxonomic_composition_family.rds"))
tpg_taxonomic_comp_family[, group := paste0(group, "_fam")]

# Full overview table taxonomic composition
# TODO: check warning in dcast!
tpg_taxonomic_comp_family %>%
  .[, prop_ord := round(prop_ord, digits = 4) * 100] %>%
  .[, group := factor(group,
    levels = c(
      "T1_fam",
      "T2_fam",
      "T3_fam",
      "T4_fam",
      "T5_fam",
      "T6_fam",
      "T7_fam",
      "T8_fam",
      "T9_fam",
      "T10_fam",
      "T11_fam",
      "T12_fam",
      "T13_fam",
      "T14_fam",
      "T15_fam"
    ),
    ordered = TRUE
  )] %>%
  dcast(., ... ~ order, value.var = "prop_ord") %>%
  fwrite(., file.path(path_paper, "Tables", "taxonomic_composition_tpgs_family.csv"))

# Taxonomic composition
tpg_taxonomic_comp[group %in% c("T5_fam", "T12_fam"), ] %>%
  dcast(., ... ~ order, value.var = "prop_ord")%>%
  fwrite(., file.path(path_paper, "Tables", "taxonomic_composition_tpgs_family_t5_t12.csv"))

# Genus lvl
tpg_taxonomic_comp_genus <- readRDS(file.path(path_cache, "tpg_taxonomic_composition_genus.rds"))
tpg_taxonomic_comp_genus[, group := paste0(group, "_genus")]

# Full overview table taxonomic composition
tpg_taxonomic_comp_genus %>%
  .[, prop_ord := round(prop_ord, digits = 4) * 100] %>%
  .[, group := factor(group,
    levels = c(
      "T1_genus",
      "T2_genus",
      "T3_genus",
      "T4_genus",
      "T5_genus",
      "T6_genus",
      "T7_genus",
      "T8_genus",
      "T9_genus",
      "T10_genus",
      "T11_genus",
      "T12_genus",
      "T13_genus",
      "T14_genus",
      "T15_genus"
    ),
    ordered = TRUE
  )] %>%
  dcast(., ... ~ order, value.var = "prop_ord") %>%
  fwrite(., file.path(path_paper, "Tables", "taxonomic_composition_tpgs_genus.csv"))

# PDPs ----
# For now only for Midwest

# Load data used to fit BRT models
data_cwm_final <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
data_cwm_midwest <- dcast(data_cwm_final$Midwest, site + max_log_tu ~ trait,
  value.var = "cwm_val"
)
trait_names <- unique(data_cwm$Midwest$trait)

trait_groups_rel_final <- readRDS(file.path(path_cache, "trait_groups_rel_final.rds"))
trait_groups_rel_final_midwest <- trait_groups_rel_final$Midwest
var_names <- names(trait_groups_rel_final_midwest)[names(trait_groups_rel_final_midwest) %like% "^T[0-9]"]

# Interpret model and create PDPs
# CWM approach
xgboost_exp_midwest_cwm <- explain_mlr3(
  model = res_xgboost_comb$cwm$Midwest$final_model,
  data = data_cwm_midwest[, .SD, .SDcols = trait_names],
  y = data_cwm_midwest$max_log_tu,
  label = "XGBOOST_cwm_midwest"
)
pdp_midwest_cwm <- model_profile(xgboost_exp_midwest_cwm)$agr_profiles
plot(pdp_midwest_cwm) +
  theme_bw() +
  theme(
    legend.position = "none", 
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
    )
  )
ggsave(
  filename = file.path(path_out, "Graphs", "PDP_Midwest_cwm.png"),
  device = "png",
  width = 50,
  height = 30,
  units = "cm"
)

# TPG_rel approach
xgboost_exp_midwest_tpg <- explain_mlr3(
  model = res_xgboost_comb$tpgs_rel$Midwest$final_model,
  data = trait_groups_rel_final_midwest[, .SD, .SDcols = var_names],
  y = trait_groups_rel_final_midwest$max_log_tu,
  label = "XGBOOST_tpg_rel_midwest"
)
pdp_midwest_tpg <- model_profile(xgboost_exp_midwest_tpg)$agr_profiles
plot(pdp_midwest_tpg) +
  theme_bw() +
  theme(
    legend.position = "none",
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
    )
  )
ggsave(
  filename = file.path(path_out, "Graphs", "PDP_Midwest_tpgs_rel.png"),
  device = "png",
  width = 50,
  height = 30,
  units = "cm"
)

# Plotting CWM values/TPG values
data_midwest[, .SD, .SDcols = trait_names] %>%
  melt(., id.vars = NULL, variable.name = "trait") %>% 
  ggplot(., aes(x = as.factor(trait), y = value)) +
  geom_boxplot() +
  facet_wrap(.~as.factor(trait), scales = "free") +
  theme_bw()

trait_groups_rel_final_midwest[, .SD, .SDcols = var_names] %>%
  melt(.,
    variable.name = "tpg"
  ) %>%
  ggplot(., aes(x = as.factor(tpg), y = value)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw()
