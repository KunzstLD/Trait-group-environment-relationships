# Create three datasets ----

# Traits:
# - All traits
# - positive reacting
# - negative reacting
# Put size large into negative, as it has more negative directions when not using the BH correction)
pos_traits <- trait_direction[sum_direction > 0, cwm_trait]
neg_traits <- trait_direction[sum_direction <= 0, cwm_trait]

# Standardize all env. varibales and cwm traits
cols_to_stand = c(
  "max_log_tu",
  "SubstrateD84.M",
  "WetWidthIQR.STD",
  "DMax_42day.M",
  "Riffle.FRC",
  "Temp.median",
  "NH3_4wk.median",
  "orthoP_4wk.median",
  "NO3NO2_4wk.median"
)
data_cwm_env <- lapply(data_cwm_env, function(x)
  x[, (cols_to_stand) := lapply(.SD, scale), .SDcols = cols_to_stand])
data_cwm_env <- lapply(data_cwm_env, function(x) x[, cwm_val := scale(cwm_val), by = "trait"])

data_cwm_env_all <- lapply(data_cwm_env, function(x) x[trait %in% mr_traits, ])
data_cwm_env_pos <- lapply(data_cwm_env, function(x) x[trait %in% pos_traits,])
data_cwm_env_neg <- lapply(data_cwm_env, function(x) x[trait %in% neg_traits, ])

# TPGs:
pos_tpgs <- tpg_direction[sum_direction > 0, TPG]
neg_tpgs <- tpg_direction[sum_direction <= 0, TPG]

# Standardization env. and tpg data in TPG dataset
data_tpg_env_fam <- lapply(data_tpg_env_fam, function(x)
  x[, (cols_to_stand) := lapply(.SD, scale), .SDcols = cols_to_stand])
data_tpg_env_fam <- lapply(data_tpg_env_fam, function(x) x[, value := scale(value), by = "tpg"])

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
# saveRDS(cwm_rda, file.path(path_cache, "cwm_rda.rds"))
# saveRDS(tpg_rda, file.path(path_cache, "tpg_rda.rds"))

# Read RDA results
cwm_rda <- readRDS(file.path(path_cache, "cwm_rda.rds"))
tpg_rda <- readRDS(file.path(path_cache, "tpg_rda.rds"))

# Amount of variance explained by the first two axes
output_variance_rda <- function(rda_results) {
  rda_var_ls <- list()
  for (region in names(rda_results)) {
    tot_var <- summary(rda_results[[region]])$tot.chi
    constr_var <- summary(rda_results[[region]])$constr.chi
    adj_constr_var <- RsquareAdj(rda_results[[region]])$adj.r.squared
    unconstr_var <- summary(rda_results[[region]])$unconst.chi
    
    cum_var <- summary(rda_results[[region]])$concont$importance["Cumulative Proportion", 2]
    
    rda_table <- data.table(
      "constr_var" = constr_var / tot_var,
      "adj_constr_var"= adj_constr_var,
      "unconstr_var" = unconstr_var / tot_var,
      "cum_var_rda1_rda2_constrained" = cum_var
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

cwm_rda_variance[tpg_rda_variance, `:=`(
  constr_var_tpgs = i.constr_var,
  adj_constr_var_tpgs = i.adj_constr_var,
  unconstr_var_tpgs = i.unconstr_var,
  cum_var_rda1_rda2_constrained_tpgs = i.cum_var_rda1_rda2_constrained
), on = c("direction", "region")]
cwm_rda_variance[, .(
  direction,
  region,
  constr_var_traits = round(constr_var, digits = 2),
  adj_constr_var_traits = round(adj_constr_var, digits = 2), 
  constr_var_tpgs = round(constr_var_tpgs, digits = 2),
  adj_constr_var_tpgs = round(adj_constr_var_tpgs, digits = 2),
  cum_var_rda1_rda2_traits = round(cum_var_rda1_rda2_constrained, digits = 2),
  cum_var_rda1_rda2_tpgs = round(cum_var_rda1_rda2_constrained_tpgs, digits = 2)
)] |> 
  fwrite(file.path(path_out, "rda_variance.csv"))

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

### Consistency & Specificity ----
# RDA relationship between X and Y
# For scaling=2
# Consistency represented as: stat sign of explanatory variables + Correlation between response and explanatory variables, represented by 
# the angle between species scores and biplot scores 
# "The angles in the biplot between response and explanatory variables, 
# and between response variables themselves or explanatory variables themselves, 
# reflect their correlations."

# Function to calculate cosine similarity (angle-based correlation)
cosine_similarity <- function(rda_results) {
  
  # Apply function to each RDA object
  similarity_list <- lapply(rda_results, function(x) {
    
    # Extract species scores (response variables)
    X <- scores(x, choices = 1:2, scaling = 2, display = "sp")
    
    # Extract biplot scores (explanatory variables)
    Y <- scores(x, choices = 1:2, scaling = 2, display = "bp")
    
    # Normalize each row (vector) to get unit vectors
    X_norm <- X / sqrt(rowSums(X^2))
    Y_norm <- Y / sqrt(rowSums(Y^2))
    
    # Compute cosine similarity matrix
    similarity_matrix <- X_norm %*% t(Y_norm)
    similarity_matrix <- as.data.table(similarity_matrix, keep.rownames = TRUE)
    
    return(similarity_matrix)
  })
  
  return(similarity_list)
}
cwm_rda_cosine <- lapply(cwm_rda, function(x)
  cosine_similarity(x)) |>
  lapply(function(x)
    rbindlist(x, idcol = "region")) |>
  rbindlist(idcol = "direction")

tpg_rda_cosine <- lapply(tpg_rda, function(x)
  cosine_similarity(x)) |> 
  lapply(function(x)
    rbindlist(x, idcol = "region")) |>
  rbindlist(idcol = "direction")

rda_cosine <- rbindlist(list("traits" = cwm_rda_cosine, "tpg" = tpg_rda_cosine), idcol =
                          "type")
saveRDS(rda_cosine, file.path(path_cache, "rda_cosine.rds"))

# Example code:
# test_rda <- cwm_rda$all$California
# spe.sc <- scores(test_rda, choices = 1:2, scaling=2, display="sp")
# biplot.sc <- scores(test_rda, choices = 1:2, scaling=2, display="bp")
# plot(test_rda, scaling=2, type="n") |> 
#    points("sites", pch=16, col="grey") |> # weighted sums of species scores - coord. of the sites as expressed in the space of the response Y
#    points("biplot", pch=16, col="red") |>  # Biplot scores for constraining variables: coord. of tips of the vectors representing the expl. variables
#    text("species", arrows = TRUE, length=0.05, col="blue") |>  # Species scores: coord. of the tips of the vectors representing the response variables
#    text("biplot", arrows = TRUE, length=0.05, col="red")

#_______________________________________________________________________________
# Direction of each effect ---
# Extract biplot scores for environmental variables (first two axes) and check for sign consistency
# Use the sum or the largest absolute value to define the direction for the first two axes!
# RDA space: constrained ordination space -> Variation in Traits/TPGs is explained by env. variables (X)
# biplot_scores <- function(rda_results) {
#   env_scores_trait <- lapply(rda_results, function(x)
#     as.data.table(scores(x, display = "bp"), keep.rownames = TRUE))
#   env_scores_trait <- rbindlist(env_scores_trait, idcol = "region")
#   setnames(env_scores_trait, "rn", "env_factor")
#   env_scores_trait[, direction_env_factor := sign(RDA1 + RDA2)]
#   return(env_scores_trait)
# }
# biplot_scores <- function(rda_results) {
#   env_scores_trait <- lapply(rda_results, function(x)
#     as.data.table(scores(x, display = "bp"), keep.rownames = TRUE))
#   
#   env_scores_trait <- rbindlist(env_scores_trait, idcol = "region")
#   setnames(env_scores_trait, "rn", "env_factor")
#   
#   # Compute angle in degrees
#   env_scores_trait[, angle_env_factor := atan2(RDA2, RDA1) * (180 / pi)]
#   
#   # Compute vector length (magnitude)
#   env_scores_trait[, magnitude := sqrt(RDA1^2 + RDA2^2)]
#   
#   # Compute (pos or neg.) direction as sign of both RDA axes
#   env_scores_trait[, direction_env_factor := sign(RDA1 + RDA2)]
#   
#   return(env_scores_trait)
# }
# cwm_rda_scores <- lapply(cwm_rda, function(cwm) biplot_scores(cwm))
# tpg_rda_scores <- lapply(tpg_rda, function(tpg) biplot_scores(tpg))
#_______________________________________________________________________________

cwm_rda_me <- rbindlist(cwm_rda_me, idcol="direction")
tpg_rda_me <- rbindlist(tpg_rda_me, idcol="direction")
rda_me <- rbindlist(list("traits" = cwm_rda_me, "tpg" = tpg_rda_me), idcol="type")
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
  _[order(type, direction, -N), ] |> 
  dcast(...~type, value.var = "N") |> 
  fwrite(file.path(path_out, "rda_n_sign_marginal_effect.csv"))

rda_me[, significant := fifelse(p.value <= 0.05, "Y", "N")]
# saveRDS(rda_me, file.path(path_cache, "rda_me.rds"))

# Direction of coefficient of max log tu & env. factors
# Are the other env. factors consistently responding across the multiple regions?
cwm_rda_scores <- rbindlist(cwm_rda_scores, idcol="direction")
tpg_rda_scores <- rbindlist(tpg_rda_scores, idcol="direction")
rda_scores <- rbindlist(list("traits" = cwm_rda_scores, "tpg" = tpg_rda_scores), idcol="type")
rda_scores[rda_me, significant := i.significant, on = c("type", "direction", "region", "env_factor")]
saveRDS(rda_scores, file.path(path_cache, "rda_scores.rds"))

# All effects
# rda_scores[, .(sum_direction = sum(direction_env_factor)), by = c("type", "direction", "env_factor")] |> 
#   _[order(type, direction, -sum_direction), ]
# 
# # Only significant effects
# rda_scores[significant == "Y", .(sum_direction = sum(direction_env_factor)), by = c("type", "direction", "env_factor")] |>
#   _[order(type, direction, -sum_direction), ]
# rda_scores[type == "tpg" & direction == "neg" & env_factor == "max_log_tu", ]


### Performance & Specificity----
# Marginal effect/variance of max_log_tu & the other env. factors

# Output plot of all marginal effects
ggplot(rda_me[env_factor != "Residual", ],
       aes(
         x = region,
         y = prop_variance,
         fill = env_factor,
         alpha = significant
       )) +
  geom_bar(stat = "identity",
           position = "dodge",
           width = 0.9) +
  facet_grid(type ~ direction, 
             labeller = labeller(
               type = c("tpg" = "TPGs", 
                        "traits" = "Traits"))) +
  scale_fill_manual(
    values = c("yellow3", "steelblue", "tomato1", "mistyrose4"),
    labels = c(
      "max_log_tu" = "Max logTU",
      "Riffle.FRC" = "Riffle Fraction",
      "Temp.median" = "Temperature",
      "orthoP_4wk.median" = "Phosphate"
    )
  ) +
  scale_alpha_manual(values = c("Y" = 1, "N" = 0.35),
                     labels = c("Y" = "Yes", "N" = "No")) +
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
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  labs(x = "Region",
       y = "Proportion of Variance",
       fill = "Environmental Factor",
       alpha = "Marginal effect stat. significant?")
ggsave(file.path(path_out, "prop_var_me.png"))

# Output table 
rda_me_publ <- rda_me[, .(type, direction, region, env_factor, Variance, p.value, prop_variance)]