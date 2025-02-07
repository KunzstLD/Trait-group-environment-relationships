# ______________________________________________________________________________
# Comparing Performance , Consistency and Specificity for
# - selected TPGs (responsive towards Pesticdes)
# - selected Traits (responsive towards Pesticdes)
# - SPEAR 
# - EPT Taxa
# towards estimated pesticide toxicity based on max logTU 
# ______________________________________________________________________________

# ______________________________________________________________________________

# Preprocessing ----
# Load responsive TPGs, Traits, max logTU & further env. factors
data_cwm_env <- readRDS(file.path(path_cache, "data_cwm_env.rds"))
data_tpg_env_fam <- readRDS(file.path(path_cache, "data_tpg_env_fam.rds"))
names(data_cwm_env)[4] <- "Northwest"
names(data_tpg_env_fam)[4] <- "Northwest"

# Traits
# feed parasite has been additionally included, but we exclude it here, since 
# it was only responsive in one region
mr_traits <- c(
  "size_large",
  "size_small",
  "resp_gil",
  "feed_gatherer",
  "feed_filter",
  "feed_predator",
  "volt_bi_multi",
  "volt_semi",
  "locom_swim",
  "sensitivity_organic"
)
# TPGs
mr_tpgs <- data_tpg_env_fam$California$tpg |> unique()

# ______________________________________________________________________________
# Individual regression results ----
# Select those traits & TPGs positively and negatively reacting to pesticides 
cwm_table <- readRDS(file.path(path_cache, "cwm_tabl_publ.rds"))
cwm_table <- cwm_table[cwm_trait %in% mr_traits, ]
tpg_table_fam <- readRDS(file.path(path_cache, "tpg_tabl_publ.rds"))

# Calc. direction of effects (for LM & GAM)
cwm_table[term != "(Intercept)", direction_max_log_tu := coalesce(sign(estimate), sign(median_derivative))]
tpg_table_fam[term != "(Intercept)", direction_max_log_tu := coalesce(sign(estimate), sign(median_derivative))]

# Check for direction of effect and stat. Significance
# First check if effect is stat significant. If this is the case then check the direction

# Traits
# The effect of size large is in different directions (depending on the region)
# and when taking statistical significance into account, indeed 
# size large reacts not consistently 
# cwm_table[term != "(Intercept)" &
#             p_p_gam <= 0.05, ] |> 
#   _[, sum(direction_max_log_tu, na.rm = TRUE), by = "cwm_trait"] |> 
#   _[order(V1), ]
trait_direction <- cwm_table[term != "(Intercept)" &
                               p_p_gam_BH <= 0.05, ] |>
  _[,  .(sum_direction= sum(direction_max_log_tu, na.rm = TRUE)), by = "cwm_trait"] |>
  _[order(sum_direction), ]

# TPGs (no difference without BH correction)
tpg_direction <- tpg_table_fam[term != "(Intercept)" &
                                 p_p_gam_BH <= 0.05, ] |>
  _[, .(sum_direction = sum(direction_max_log_tu, na.rm = TRUE)), by = "TPG"] |>
  _[order(sum_direction), ]

# TPG9_fam seems not to be "responsive" in the individual regressions
# tpg_table_fam[TPG=="TPG9_fam", ]

# ______________________________________________________________________________
# Create three datasets ----

# Traits:
# - All traits
# - positive reacting
# - negative reacting
# Put size large into negative, as it has more negative directions when not using the BH correction)
pos_traits <- trait_direction[sum_direction > 0, cwm_trait]
neg_traits <- trait_direction[sum_direction <= 0, cwm_trait]
  
data_cwm_env_all <- lapply(data_cwm_env, function(x) x[trait %in% mr_traits, ])
data_cwm_env_pos <- lapply(data_cwm_env, function(x) x[trait %in% pos_traits,])
data_cwm_env_neg <- lapply(data_cwm_env, function(x) x[trait %in% neg_traits, ])

# TPGs:
pos_tpgs <- tpg_direction[sum_direction > 0, TPG]
neg_tpgs <- tpg_direction[sum_direction <= 0, TPG]

# Dataset for all TPGs doesn't need to be subsetted 
data_tpg_env_fam_pos <- lapply(data_tpg_env_fam, function(x) x[tpg %in% pos_tpgs,])
data_tpg_env_fam_neg <- lapply(data_tpg_env_fam, function(x) x[tpg %in% neg_tpgs, ])

# Transform to wf
transform_to_wf <- function(data_list, var_col, id_col, val_col, env_cols) {
  lapply(data_list, function(x) {
    cols <- c(var_col, id_col, val_col, env_cols)
    
    # Subset data to relevant env factors
    x <- x[, .SD, .SDcols = cols]
    
    # Dynamically construct dcast formula
    # Should be site + max_log_tu + remaining env. factors
    formula_str <- paste(c(id_col, env_cols), collapse = " + ")
    dcast_formula <- as.formula(paste(formula_str, "~", var_col))
    
    # Apply dcast
    x_wf <- dcast(x, dcast_formula, value.var = val_col)
    return(x_wf)
  })
}
# Traits:
data_cwm_env_all_wf <- transform_to_wf(
  data_cwm_env_all,
  var_col = "trait",
  id_col = "site",
  val_col = "cwm_val",
  env_cols = c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median"
  )
)
data_cwm_env_pos_wf <- transform_to_wf(
  data_cwm_env_pos,
  var_col = "trait",
  id_col = "site",
  val_col = "cwm_val",
  env_cols = c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median"
  )
)
data_cwm_env_neg_wf <- transform_to_wf(
  data_cwm_env_neg,
  var_col = "trait",
  id_col = "site",
  val_col = "cwm_val",
  env_cols = c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median"
  )
)

# TPGs
data_tpg_env_fam_wf <- transform_to_wf(
  data_tpg_env_fam,
  var_col = "tpg",
  id_col = "site",
  val_col = "value",
  env_cols = c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median"
  )
)
data_tpg_env_fam_pos_wf <- transform_to_wf(
  data_tpg_env_fam_pos,
  var_col = "tpg",
  id_col = "site",
  val_col = "value",
  env_cols = c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median"
  )
)
data_tpg_env_fam_neg_wf <- transform_to_wf(
  data_tpg_env_fam_neg,
  var_col = "tpg",
  id_col = "site",
  val_col = "value",
  env_cols = c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median"
  )
)

# Create matrices for RDA
create_matrices <- function(data_list, var_col) {
  # Y matrices with traits
  Y <- lapply(data_list, function(x)
    x[, .SD, .SDcols = c("site", var_col)])
  
  # X matrices
  # Max log TU + env. factors
  X_env <- lapply(data_list, function(x)
    x[, .SD, .SDcols = c("site",
                         "max_log_tu",
                         "Riffle.FRC",
                         "Temp.median",
                         "orthoP_4wk.median")])
  
  return (list("Y" = Y, "X" = X_env))
}

# Traits
cwm_data_lists <- list(
  all = list(data = data_cwm_env_all_wf, var_col = mr_traits),
  pos = list(data = data_cwm_env_pos_wf, var_col = pos_traits),
  neg = list(data = data_cwm_env_neg_wf, var_col = neg_traits)  #
)
cwm_matrices <- lapply(cwm_data_lists, function(x)
  create_matrices(x$data, x$var_col))

# TPGs
tpg_data_lists <- list(
  all = list(data = data_tpg_env_fam_wf, var_col = mr_tpgs),
  pos = list(data = data_tpg_env_fam_pos_wf, var_col = pos_tpgs),
  neg = list(data = data_tpg_env_fam_neg_wf, var_col = neg_tpgs)  #
)
tpg_matrices <- lapply(tpg_data_lists, function(x)
  create_matrices(x$data, x$var_col))

# ______________________________________________________________________________
# Check correlations among traits & TPGs 
calc_correlations <- function(matrices) {
  corr_pl_y_list <- list()
  corr_matrix_y_ls <- list()
  for (region in names(matrices$all$Y)) {
    corr_matrix <- cor(matrices$all$Y[[region]][, -"site"], use = "complete.obs")
    corr_matrix_y_ls[[region]] <- as.data.table(corr_matrix, keep.rownames = TRUE)
    # Visualize the correlation matrix
    corrplot(corr_matrix,
             method = "number",
             type = "upper",
             tl.cex = 0.8)
    
    corr_pl_y_list[[region]] <- recordPlot()
  }
  corr_matrix_y <- rbindlist(corr_matrix_y_ls, idcol = "region")
  return(corr_matrix_y)
}

# Traits
corr_matrix_cwm <- calc_correlations(cwm_matrices)
setnames(corr_matrix_cwm, "rn", "traits")
corr_matrix_cwm <- melt(
  corr_matrix_cwm,
  id.vars = c("region", "traits"),
  variable.name = "trait_compared",
  value.name = "corr"
)
corr_matrix_cwm[corr >=0.6 & corr != 1, ]
corr_matrix_cwm[corr <=-0.6 & corr != 1, ]

# TPGs
corr_matrix_tpg <- calc_correlations(tpg_matrices)
setnames(corr_matrix_tpg, "rn", "tpg")
corr_matrix_tpg <- melt(
  corr_matrix_tpg,
  id.vars = c("region", "tpg"),
  variable.name = "tpg_compared",
  value.name = "corr"
)
corr_matrix_tpg[corr >=0.6 & corr != 1, ]
corr_matrix_tpg[corr <=-0.6 & corr != 1, ]

# ______________________________________________________________________________
# RDA ----
# RDA with max logTU + P + T + Riffles
calc_rda <- function(X_env, Y_mat) {
  Y_mat <- lapply(Y_mat, function(x)
    x[, !"site"])
  rda <- Map(
    function(X, Y)
      rda(
        as.matrix(Y) ~ max_log_tu + Riffle.FRC + Temp.median + orthoP_4wk.median,
        data = X,
        na.action = na.omit
      ),
    X_env,
    Y_mat
  )
  
  return(rda)
}
cwm_rda <- lapply(cwm_matrices, function(cwm) calc_rda(cwm$X, cwm$Y))
tpg_rda <- lapply(tpg_matrices, function(tpg) calc_rda(tpg$X, tpg$Y))

# Amount of variance explained by the first
output_variance_rda <- function(rda_results) {
  rda_var_ls <- list()
  for (region in names(rda_results)) {
    tot_var <- summary(rda_results[[region]])$tot.chi
    constr_var <- summary(rda_results[[region]])$constr.chi
    unconstr_var <- summary(rda_results[[region]])$unconst.chi
    
    cum_var <- summary(rda_results[[region]])$concont$importance["Cumulative Proportion", 2]
    
    rda_table <- data.table(
      "constr_var" = constr_var / tot_var,
      "unconstr_var" = unconstr_var / tot_var,
      "cum_var_rda1_rda2" = cum_var
    )
    rda_var_ls[[region]] <- rda_table
    
  }
  return(rda_var_ls)
}
cwm_rda_variance <- lapply(cwm_rda, function(cwm) output_variance_rda(cwm))
cwm_rda_variance <- lapply(cwm_rda_variance, function(x) rbindlist(x, idcol = "region"))
cwm_rda_variance <- rbindlist(cwm_rda_variance, idcol = "direction")

tpg_rda_variance <- lapply(tpg_rda, function(tpg) output_variance_rda(tpg))
tpg_rda_variance <- lapply(tpg_rda_variance, function(x) rbindlist(x, idcol = "region"))
tpg_rda_variance <- rbindlist(tpg_rda_variance, idcol = "direction")


## Permutationstest for marginal effects ----

# Check collinearity of predictors, as marginal effects can be misleading when predictors
# are strongly correlated
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13800
# Correlations are all relatively low
corr_pl_x_list <- list()
for (region in names(cwm_matrices$all$X)) {
  corr_matrix <- cor(cwm_matrices$all$X[[region]][, c("max_log_tu",
                                        "Riffle.FRC",
                                        "Temp.median",
                                        "orthoP_4wk.median")], use = "complete.obs")
  
  # Visualize the correlation matrix
  corrplot(corr_matrix,
           method = "number",
           type = "upper",
           tl.cex = 0.8)
  
  corr_pl_x_list[[region]] <- recordPlot()
}
# corr_pl_x_list$California

# Perform marginal test for max_log_tu, i.e. test how much additional variance is
# explained after accounting for the other env. factors
# Anova uses a permutation test on the Y matrix, permuting its rows. (999 times)
# producing p values when the association between predictor and Y is broken
marginal_effects_rda <- function(rda_results) {
  anova_marginal_trait <- lapply(rda_results, function(x)
    anova.cca(x, by = "margin", step = TRUE))
  anova_marginal_trait <- lapply(anova_marginal_trait, function(x)
    setDT(broom::tidy(x)))
  anova_marginal_trait <- lapply(anova_marginal_trait, function(x)
    x[, total_variance := sum(Variance)] |>
      _[, prop_variance := (Variance / total_variance) * 100])
  anova_marginal_trait <- rbindlist(anova_marginal_trait, idcol = "region")
  setnames(anova_marginal_trait, "term", "env_factor")
  return(anova_marginal_trait)
}
cwm_rda_me <- lapply(cwm_rda, function(cwm) marginal_effects_rda(cwm))
tpg_rda_me <- lapply(tpg_rda, function(tpg) marginal_effects_rda(tpg))

# Check VIF
# Variance inflation values are all close to 1
vif_rda <- function(rda_results){
  vif <- lapply(rda_results, function(x) vif.cca(x)) 
  return(vif)
}
cwm_rda_vif <- lapply(cwm_rda, function(cwm) vif_rda(cwm))
tpg_rda_vif <- lapply(tpg_rda, function(tpg) vif_rda(tpg))

## Direction of each effect ---

# TODO: Just take the sum or the largest absolute value to define the direction for the first two axes!
# Extract biplot scores for environmental variables (first two axes) and check for sign consistency
# biplot_scores <- function(rda_results) {
#   env_scores_trait <- lapply(rda_results, function(x)
#     as.data.table(scores(x, display = "bp"), keep.rownames = TRUE))
#   env_scores_trait <- rbindlist(env_scores_trait, idcol = "region")
#   setnames(env_scores_trait, "rn", "env_factor")
#   env_scores_trait[, direction_env_factor := sign(RDA1)]
#   return(env_scores_trait)
# }

# Extract biplot scores for env. factors  (first two axes) and check the angle of the
# arrow:
# 1. Quadrant (0-90): positive association with RDA1 & RDA2
# 2. QUadrant (90-180): positive association with RDA1, negative with RDA2
# 3. Quadrant (180 -270): negative association with RDA1 & RDA2
# 4. Quadrant (270-360): negative association with RDA1, positive with RDA2
# biplot_scores <- function(rda_results) {
#   env_scores_trait <- lapply(rda_results, function(x)
#     as.data.table(scores(x, display = "bp"), keep.rownames = TRUE))
#   
#   env_scores_trait <- rbindlist(env_scores_trait, idcol = "region")
#   setnames(env_scores_trait, "rn", "env_factor")
#   
#   # Compute angle theta in radians and convert to degrees
#   env_scores_trait[, angle_env_factor := atan2(RDA2, RDA1) * (180 / pi)]
#   
#   return(env_scores_trait)
# }
cwm_rda_scores <- lapply(cwm_rda, function(cwm) biplot_scores(cwm))
tpg_rda_scores <- lapply(tpg_rda, function(tpg) biplot_scores(tpg))


### Consistency & Specificity ----
cwm_rda_me <- rbindlist(cwm_rda_me, idcol="direction")
tpg_rda_me <- rbindlist(tpg_rda_me, idcol="direction")
rda_me <- rbindlist(list("cwm" = cwm_rda_me, "tpg" = tpg_rda_me), idcol="type")
rda_me$env_factor <- factor(
  rda_me$env_factor,
  levels = c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median",
    "Residual"
  )
)

# Is the marginal effect of max lgTU stat. & the env. factors significant across the different regions?
rda_me[p.value <= 0.05, .N, by = c("type", "direction", "env_factor")] |>
  _[order(type, direction, -N), ]
rda_me[, significant := fifelse(p.value <= 0.05, "Y", "N")]

# Direction of coefficient of max log tu & env. factors
# Are the other env. factors consistently responding across the multiple regions?
# TODO: Check only direction when marginal effect significant?
cwm_rda_scores <- rbindlist(cwm_rda_scores, idcol="direction")
tpg_rda_scores <- rbindlist(tpg_rda_scores, idcol="direction")
rda_scores <- rbindlist(list("cwm" = cwm_rda_scores, "tpg" = tpg_rda_scores), idcol="type")
rda_scores[rda_me, significant := i.significant, on = c("type", "direction", "region", "env_factor")]

# All effects
rda_scores[, .(sum_direction = sum(direction_env_factor)), by = c("type", "direction", "env_factor")] |> 
  _[order(type, direction, -sum_direction), ]

# Only significant effects
rda_scores[significant == "Y", .(sum_direction = sum(direction_env_factor)), by = c("type", "direction", "env_factor")] |>
  _[order(type, direction, -sum_direction), ]

rda_scores[type == "tpg" & direction == "neg" & env_factor == "max_log_tu", ]


### Performance & Specificity----
# Marginal effect/variance of max_log_tu & the other env. factors

# Output plot of all marginal effects
ggplot(rda_me[env_factor != "Residual", ], aes(x = region, y = prop_variance, fill = env_factor)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(width = 0.8, preserve = "single"),
    width = 1
  ) +
  geom_text(
    data = rda_me[env_factor != "Residual" & p.value <= 0.05, ],
    aes(label = "*", x = region, y = prop_variance),
    # Position slightly above bar
    position = position_dodge2(width = 0.8, preserve = "single"),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  facet_grid(type ~ direction) +
  scale_fill_manual(
    values = c("yellow3", "steelblue", "tomato1", "mistyrose4"),
    labels = c(
      "max_log_tu" = "Max logTU",
      "Riffle.FRC" = "Riffle Fraction",
      "Temp.median" = "Temperature",
      "orthoP_4wk.median" = "Phosphate"
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(
      size = 16,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 16)
  ) +
  labs(x = "Region", y = "Proportion of Variance", fill = "Environmental Factor")

# Output table 
rda_me_publ <- rda_me[, .(type, direction, region, env_factor, Variance, p.value, prop_variance)]






# Modify EPT script -> include P+T+Riffles in the regressions
# Modify SPEAR -> include P+T+Riffles in the regressions