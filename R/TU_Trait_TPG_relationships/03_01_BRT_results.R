## Performance results  ----
# Train/val/test split approach with regularization
res_xgboost_comb <- load_data(
  path = path_cache,
  pattern = "res\\_xgboost.*",
  name_rm_pattern = "res\\_xgboost\\_"
)

# # Get prediction on test and training data
# xgboost_perform <- purrr::map(res_xgboost_comb, ~ sapply(.x, function(x) {
#   c(
#     "train" = x$pred_train,
#     "test" = x$pred_test
#   )
# }))
# xgboost_perform <- lapply(xgboost_perform, function(x) as.data.table(x, keep.rownames = TRUE)) %>%
#     rbindlist(., id = "method") %>%
#     melt(., id.vars = c("rn", "method")) %>%
#     dcast(., ... ~ rn, value.var = "value")
# setnames(xgboost_perform, "variable", "region")
# saveRDS(xgboost_perform, file.path(path_cache, "xgboost_perform.rds"))

# Test and training RMSE
xgboost_perform <- readRDS(file.path(path_cache, "xgboost_perform.rds"))

# Publication table
perform_publ_tbl <- xgboost_perform[
  method %in% c("cwm", "tpgs_rel_family", "tpgs_rel_genus"),
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
    "cwm_test_mse",
    "tpgs_rel_family_train_mse",
    "tpgs_rel_family_test_mse",
    "tpgs_rel_genus_train_mse",
    "tpgs_rel_genus_test_mse"
  )
)
setnames(
  perform_publ_tbl,
  names(perform_publ_tbl),
  c("Region",
  "MSE Train CWM",
  "MSE Test CWM",
  "MSE Train TPG (family)",
  "MSE Test TPG (family)",
  "MSE Train TPGS (genus)",
  "MSE Test TPG (genus)")
)
fwrite(
  perform_publ_tbl,
  file.path(path_paper, "Tables", "xgboost_perfom.csv")
)

# Most important traits ----

# Use impurity importance
xgboost_imp <- purrr::map(res_xgboost_comb, ~ sapply(.x, function(x) {
  "train" <- x$importance
}))
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
traits <- unique(xgboost_imp$trait_TPG)[1:20]
lookup_traits <- data.table(
  trait = traits,
  trait_label = c(
    "predator",
    "large size",
    "gills",
    "swimming",
    "gatherer",
    "burrowing",
    "bi/multivolt.",
    "medium size",
    "herbivore",
    "crawling",
    "filterer",
    "univolt.",
    "semivolt.",
    "plastron & spi.",
    "parasite",
    "small size",
    "S_org",
    "shredder",
    "tegument",
    "sessil"
  )
)

# Add grouping features for cwm and cws
xgboost_imp[
  method %in% c("cwm", "cws"),
  grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait_TPG)
]
xgboost_imp[lookup_traits,
  trait_label := i.trait_label, on = c("trait_TPG" = "trait")
]

#### Correlation CWM and TPG importance ranks ####

# Give correct names for TPGs
xgboost_imp[method == "tpgs_rel_family", trait_TPG := paste0("TPG", sub("T", "", trait_TPG), "_fam")]
xgboost_imp[method == "tpgs_rel_genus", trait_TPG := paste0("TPG", sub("T", "", trait_TPG), "_genus")]

# Split into three datasets
xgboost_cwm <- xgboost_imp[method=="cwm", ]
xgboost_tpg_fam <- xgboost_imp[method=="tpgs_rel_family", ]
xgboost_tpg_gen <- xgboost_imp[method=="tpgs_rel_genus", ]

# Convert to factor
xgboost_cwm[, trait_TPG := factor(trait_TPG)]
levels(xgboost_cwm$trait_TPG)

##### CWM ----

# Add a zero imp. score for each trait that did not occur in the ranking
regions <- unique(xgboost_imp$region)
traits_expanded <- list()

for (i in regions) {
    data <- xgboost_cwm[region == i, ]

    traits_ <- as.character(data[, trait_TPG])
    traits_unranked_ <- lookup_traits[!traits %in% traits_, trait]

    traits_expanded[[i]] <- data.frame(
      "traits" = c(traits_, traits_unranked_),
      "score" = c(data$score, rep(0, length(traits_unranked_)))
    )
}
traits_expanded <- rbindlist(traits_expanded, idcol="region")

# Calc the correlation
# Region grid for comparison
comparison_grid  <- as.data.frame(t(combn(unique(xgboost_imp$region), 2)))
setDT(comparison_grid)

corr_result <- list()
for (j in 1:nrow(comparison_grid)) {
  region1 <- comparison_grid[j, V1]
  region2 <- comparison_grid[j, V2]
  
  df_1 <- traits_expanded[region == region1, ]
  df_2 <- traits_expanded[region == region2, ]
  
  # Match according to traits, then compare by score
  test_stat <- cor.test(df_1[, score], 
                        df_2[order(match(traits, df_1$traits)), score], 
                          method ="spearman",
                          exact = FALSE, 
                          alternative = "two.sided"
                        )
  
  corr_result[[j]] <- data.frame("rho" = test_stat$estimate,
                                 "p_val" = test_stat$p.value)
}
corr_result <- rbindlist(corr_result)
comparison_grid[, `:=`(cwm_corr = corr_result$rho, 
                       p_val = corr_result$p_val)]


##### TPGs ----
tpgs_fam <- unique(xgboost_tpg_fam$trait_TPG)
tpgs_gen <- unique(xgboost_tpg_gen$trait_TPG)

tpgs_fam_expanded <- list()
tpgs_gen_expanded <- list()

# Expand with zero score
for (i in regions) {
  data_fam_ <- xgboost_tpg_fam[region == i, ]
  data_gen_ <- xgboost_tpg_gen[region == i, ]
  
  tpgs_fam_expanded_ <- as.character(data_fam_[, trait_TPG])
  tpgs_gen_expanded_ <- as.character(data_gen_[, trait_TPG])
  
  tpgs_fam_unranked_ <- tpgs_fam[!tpgs_fam %in% tpgs_fam_expanded_]
  tpgs_gen_unranked_ <- tpgs_gen[!tpgs_gen %in% tpgs_gen_expanded_]
  
  tpgs_fam_expanded[[i]] <- data.frame(
    "tpgs_fam" = c(tpgs_fam_expanded_, tpgs_fam_unranked_),
    "score" = c(data_fam_$score, rep(0, length(
      tpgs_fam_unranked_
    )))
  )
  
  tpgs_gen_expanded[[i]] <- data.frame(
    "tpgs_gen" = c(tpgs_gen_expanded_, tpgs_gen_unranked_),
    "score" = c(data_gen_$score, rep(0, length(
      tpgs_gen_unranked_
    )))
  )
}

# Calc the correlation
tpgs_fam_expanded <- rbindlist(tpgs_fam_expanded, idcol="region")
tpgs_gen_expanded <- rbindlist(tpgs_gen_expanded, idcol="region")

corr_result_tpgs_fam <- list()
corr_result_tpgs_gen <- list()

for (j in 1:nrow(comparison_grid)) {
  region1_ <- comparison_grid[j, V1]
  region2_ <- comparison_grid[j, V2]
  
  df_fam_1 <- tpgs_fam_expanded[region==region1_,]
  df_fam_2 <- tpgs_fam_expanded[region==region2_,]
  
  df_gen_1 <- tpgs_gen_expanded[region==region1_,]
  df_gen_2 <- tpgs_gen_expanded[region==region2_,]
  
  test_stat_fam <- cor.test(
    df_fam_1[, score],
    df_fam_2[order(match(tpgs_fam, df_fam_1$tpgs_fam)), score],
     method="spearman",
     exact=FALSE,
     alternative = "two.sided"
  )
  
  corr_result_tpgs_fam[[j]] <- data.frame("rho" = test_stat_fam$estimate,
                                          "p_val" = test_stat_fam$p.value)
  
  test_stat_gen <- cor.test(
    df_gen_1[, score],
    df_gen_2[order(match(tpgs_gen, df_gen_1$tpgs_gen)), score],
      method="spearman",
      exact=FALSE,
      alternative = "two.sided"
  )
  
  corr_result_tpgs_gen[[j]] <- data.frame("rho" = test_stat_gen$estimate,
                                          "p_val" = test_stat_gen$p.value)
}
corr_result_tpgs_fam <- rbindlist(corr_result_tpgs_fam)
corr_result_tpgs_gen <- rbindlist(corr_result_tpgs_gen)

comparison_grid[, `:=`(
  tpgs_fam_corr = corr_result_tpgs_fam$rho,
  tpgs_fam_p_val = corr_result_tpgs_fam$p_val,
  tpgs_gen_corr = corr_result_tpgs_gen$rho,
  tpgs_gen_p_val = corr_result_tpgs_gen$p_val
)]

# Save table
comparison_grid[, "Compared region" := paste(V1, "-", V2)]
setnames(comparison_grid, 
         old = c("cwm_corr", 
                 "p_val",
                 "tpgs_fam_corr",
                 "tpgs_fam_p_val",
                 "tpgs_gen_corr", 
                 "tpgs_gen_p_val"),
         new = c("Corr. CWM traits", 
                 "P val Corr CWM traits", 
                 "TPGs correlation (family level)", 
                 "TPGs P val (family level)", 
                 "TPGs correlation (genus level)", 
                 "TPGs P val (genus level)"))
comparison_grid[, `Compared region` := sub("PN", "Southwest", `Compared region`)]

num_cols <- names(Filter(is.numeric, comparison_grid))
comparison_grid[, (num_cols) := lapply(.SD, function(x)
  round(x, 2)), .SDcols = num_cols]

fwrite(comparison_grid[, .SD, .SDcols = c("Compared region", num_cols)], 
       file.path(path_paper, "Tables", "Correlations_rankings.csv"))

#### Top 5 important/responsive traits/TPGs per region ####

# CWM:
# size large 3 times,
# many traits two times: small, gills, gatherer, volt_bi_multi, S_org,
# volt_semi, predator, swimming, filterer
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
  .[region == "PN", region := "Northwest"] %>%
  .[order(region), ] %>%
  setnames(
    .,
    c("region", "trait_label", "score"),
    c("Region", "Trait", "Impurity score")
  ) %>%
  fwrite(., file.path(path_paper, "Tables", "Most_imp_traits.csv"))

# Dominance of feeding mode traits
# (4 regions)
xgboost_imp[method == "cwm", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = c("grouping_feature", "region")] %>% 
  .[order(grouping_feature, -N), ]

# Abundance weighted fraction (TPGS_REL) ----
# Family-level approach
# T1 in 4 regions
# T2, T5, T8, T9, T12, T10 in 3 regions 
xgboost_imp[method == "tpgs_rel_family", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = "trait_TPG"] %>% 
  .[order(-N), ] 

# Table for publication
xgboost_imp[method == "tpgs_rel_family", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>%
  .[!is.na(score), .(region, trait_TPG, score = round(score, digits = 2))] %>% 
  .[, trait_TPG := paste0(trait_TPG, "_fam")]  %>% 
  .[region == "PN", region := "Northwest"] %>% 
  .[order(region), ] %>%
  setnames(
    .,
    c("region", "trait_TPG", "score"),
    c("Region", "TPG", "Impurity score")
  ) %>% 
  fwrite(., file.path(path_paper, "Tables", "Most_imp_TPGs.csv"))

# Genus-level approach
# T4, T10, T12 in three regions
xgboost_imp[method == "tpgs_rel_genus", trait_TPG := paste0(trait_TPG, "_genus")]
xgboost_imp[method == "tpgs_rel_genus", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>% 
  .[, .N, by = "trait_TPG"] %>% 
  .[order(-N), ] 

xgboost_imp[method == "tpgs_rel_genus", ] %>%
  .[trait_TPG %in% c("T4_genus", "T10_genus", "T12_genus"), ]

# Table for publication
xgboost_imp[method == "tpgs_rel_genus", ] %>%
  .[order(-score), .SD[1:5, ], by = "region"] %>%
  .[!is.na(score), .(region, trait_TPG, score = round(score, digits = 2))] %>% 
  .[, trait_TPG := paste0("TPG", sub("T", "", trait_TPG), "_genus")]  %>% 
  .[region == "PN", region := "Northwest"] %>% 
  .[order(region), ] %>%
  setnames(
    .,
    c("region", "trait_TPG", "score"),
    c("Region", "TPG", "Impurity score")
  ) %>% 
  fwrite(., file.path(path_paper, "Tables", "Most_imp_TPGs_genus.csv"))

# Load defining trait combinations
# defining_traits_family <- readRDS(file.path(path_cache, "defining_traits_family.rds"))
# defining_traits_genus <- readRDS(file.path(path_cache, "defining_traits_genus.rds"))

# Summary table for all methods 
xgboost_imp[method %in% c("cwm", "tpgs_rel_family", "tpgs_rel_genus")] %>% 
  .[, score := round(score, digits = 2)] %>% 
  .[order(-score), .SD[1:5, ], by = c("method", "region")]  %>% 
  dcast(., method+trait_TPG~region, value.var = "score") %>% 
  .[method == "tpgs_rel_family", trait_TPG := paste0("TPG", sub("T", "", trait_TPG), "_fam")]  %>% 
  .[method == "tpgs_rel_genus", trait_TPG := paste0("TPG", sub("T", "", trait_TPG), "_genus")] %>% 
  setnames(
    .,
    c("trait_TPG", "method", "PN"),
    c("TPG", "Method", "Northwest")
  ) %>%
  print()
  fwrite(., file.path(path_paper, "Tables", "Most_imp_TPGs_summary.csv"))


# Taxonomic composition of TPGs ----
## Family level ----
tpg_taxonomic_comp_family <- readRDS(file.path(path_cache, "tpg_taxonomic_composition_family.rds"))
tpg_taxonomic_comp_family[, group := paste0(group, "_fam")]
tpg_taxonomic_comp_family[, group := paste0("TPG", group)]

# Full overview table taxonomic composition
tpg_taxonomic_comp_family %>%
  .[, prop_ord := round(prop_ord, digits = 4) * 100] %>%
  .[, group := factor(group,
    levels = c(
      "TPG1_fam",
      "TPG2_fam",
      "TPG3_fam",
      "TPG4_fam",
      "TPG5_fam",
      "TPG6_fam",
      "TPG7_fam",
      "TPG8_fam",
      "TPG9_fam",
      "TPG10_fam",
      "TPG11_fam",
      "TPG12_fam",
      "TPG13_fam",
      "TPG14_fam",
      "TPG15_fam"
    ),
    ordered = TRUE
  )] %>%
  .[Region == "PN", Region := "Northwest"] %>% 
  setnames(., "group", "TPG") %>% 
  dcast(., ... ~ order, value.var = "prop_ord") %>%
  fwrite(., file.path(path_paper, "Tables", "taxonomic_composition_tpgs_family.csv"))

# Check taxonomic composition for most consistent TPG
tpg_taxonomic_comp_family[TPG == "TPG1_fam", ] %>% 
.[order(Region, -prop_ord), ] 

tpg_taxonomic_comp_family[TPG %in% c(
  "TPG2_fam", 
  "TPG5_fam",
  "TPG8_fam",
  "TPG9_fam",
  "TPG10_fam",
  "TPG12_fam"), ] %>%
.[order(Region, TPG, -prop_ord), ] 

# Composition for those TPGs on family level that did not occur in every region 
tpg_taxonomic_comp_family[TPG %in% c(
  "TPG11_fam",
  "TPG14_fam"
), ]


## Genus lvl ----
tpg_taxonomic_comp_genus <- readRDS(file.path(path_cache, "tpg_taxonomic_composition_genus.rds"))
tpg_taxonomic_comp_genus[, group := paste0("TPG", group, "_genus")]

# Full overview table taxonomic composition
tpg_taxonomic_comp_genus %>%
  .[, prop_ord := round(prop_ord, digits = 4) * 100] %>%
  .[, group := factor(group,
    levels = c(
      "TPG1_genus",
      "TPG2_genus",
      "TPG3_genus",
      "TPG4_genus",
      "TPG5_genus",
      "TPG6_genus",
      "TPG7_genus",
      "TPG8_genus",
      "TPG9_genus",
      "TPG10_genus",
      "TPG11_genus",
      "TPG12_genus",
      "TPG13_genus",
      "TPG14_genus",
      "TPG15_genus"
    ),
    ordered = TRUE
  )] %>%
  .[Region == "PN", Region := "Northwest"] %>% 
  setnames(., "group", "TPG") %>% 
  dcast(., ... ~ order, value.var = "prop_ord") %>%
  fwrite(., file.path(path_paper, "Tables", "taxonomic_composition_tpgs_genus.csv"))

# Check taxonomic composition for most consistent TPG
tpg_taxonomic_comp_genus[TPG == "TPG4_genus", ] %>% 
.[order(Region, -prop_ord), ] 
tpg_taxonomic_comp_genus[TPG == "TPG10_genus" , ] %>% 
.[order(Region, -prop_ord), ] 
tpg_taxonomic_comp_genus[TPG == "TPG12_genus" , ] %>% 
.[order(Region, -prop_ord), ]

# TPG7_genus not found in all regions
tpg_taxonomic_comp_genus[TPG == "TPG7_genus", ] %>% 
 .[order(Region, -prop_ord), ] 

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
