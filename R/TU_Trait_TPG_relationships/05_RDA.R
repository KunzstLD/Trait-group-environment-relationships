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
# TPGs (already subsetted to the most responsive)
mr_tpgs <- data_tpg_env_fam$California$tpg |> unique()

# ______________________________________________________________________________
# Standardize all env. variables and cwm traits
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
# data_cwm_env <- lapply(data_cwm_env, function(x)
#   x[, (cols_to_stand) := lapply(.SD, scale), .SDcols = cols_to_stand])
# data_cwm_env <- lapply(data_cwm_env, function(x) x[, cwm_val := scale(cwm_val), by = "trait"])

# Subset to most responsive traits
data_cwm_env_all <- lapply(data_cwm_env, function(x) x[trait %in% mr_traits, ])

# Standardize env. and tpg data in TPG dataset
# data_tpg_env_fam <- lapply(data_tpg_env_fam, function(x)
#   x[, (cols_to_stand) := lapply(.SD, scale), .SDcols = cols_to_stand])
# data_tpg_env_fam <- lapply(data_tpg_env_fam, function(x) x[, value := scale(value), by = "tpg"])
data_tpg_env_fam <- lapply(data_tpg_env_fam, function(x) {
  x[, tpg := factor(
    tpg,
    levels = c(
      "TPG1_fam",
      "TPG2_fam",
      "TPG5_fam",
      "TPG8_fam",
      "TPG9_fam",
      "TPG10_fam",
      "TPG12_fam",
      "TPG15_fam"
    )
  )]
  x
})

# Univariate plots (individual responses vs predictor variables) 
predictor_vector <- c(
  maxTU = "max_log_tu"# ,
  # Riffles = "Riffle.FRC",
  # Temperatue = "Temp.median",
  # Phosphate = "orthoP_4wk.median"
)
regions <- c("California", "Midwest", "Northeast", "Northwest", "Southeast")

# Loop for plotting 
for (region in regions) {
  
  subset_traits <- data_cwm_env[[region]]
  subset_tpgs   <- data_tpg_env_fam[[region]]
  
  # Plot traits
  for (label in names(predictor_vector)) {
    
    predictor_col <- predictor_vector[[label]]
    
    p1 <- ggplot(
      subset_traits[subset_traits$trait %in% mr_traits, ],
      aes(x = .data[[predictor_col]], y = cwm_val)
    ) +
      geom_point() +
      facet_wrap(~ trait) +
      # xlim(-5, 1) +
      ylim(-1, 1) +
      theme_bw() +
      theme(
        legend.position = "right",
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
      ) +
      ggtitle(label = paste(region)) +
      labs(
        y = "CWM Trait Value",
        x = paste0(label)
      )
    
    output_name <- paste0("trait_", label, "_", region, ".png")
    
    ggsave(
      file.path(path_paper, "Figures", output_name),
      p1,
      width = 12,
      height = 6.5
    )
    
    
    # Plot TPGs
    p2 <- ggplot(
      subset_tpgs,
      aes(x = .data[[predictor_col]], y = value)
    ) +
      geom_point() +
      facet_wrap(~ tpg) +
      # xlim(-5, 1) +
      ylim(0, 1) +
      theme_bw() +
      theme(
        legend.position = "right",
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
      ) +
      ggtitle(label = paste(region)) +
      labs(
        y = "CWM TPG Value",
        x = paste0(label)
      )
    
    output_name <- paste0("tpgs_", label, "_", region, ".png")
    
    ggsave(
      file.path(path_paper, "Figures", output_name),
      p2,
      width = 12,
      height = 6.5
    )
  }
}


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
cwm_data_list <- list(data = data_cwm_env_all_wf, var_col = mr_traits)
cwm_matrices <- create_matrices(cwm_data_list$data, cwm_data_list$var_col)

# TPGs
tpg_data_list <- list(data = data_tpg_env_fam_wf, var_col = mr_tpgs)
tpg_matrices <- create_matrices(tpg_data_list$data, tpg_data_list$var_col)

# ______________________________________________________________________________
# Check correlations among traits & TPGs 
calc_correlations <- function(matrices) {
  corr_pl_y_list <- list()
  corr_matrix_y_ls <- list()
  for (region in names(matrices$Y)) {
    corr_matrix <- cor(matrices$Y[[region]][, -"site"], use = "complete.obs")
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
cwm_rda <- calc_rda(cwm_matrices$X, cwm_matrices$Y) 
tpg_rda <- calc_rda(tpg_matrices$X, tpg_matrices$Y) 
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
cwm_rda_variance <- output_variance_rda(cwm_rda)
cwm_rda_variance <- rbindlist(cwm_rda_variance, idcol = "region")
tpg_rda_variance <- output_variance_rda(tpg_rda)
tpg_rda_variance <- rbindlist(tpg_rda_variance, idcol = "region")

cwm_rda_variance[tpg_rda_variance, `:=`(
  constr_var_tpgs = i.constr_var,
  adj_constr_var_tpgs = i.adj_constr_var,
  unconstr_var_tpgs = i.unconstr_var,
  cum_var_rda1_rda2_constrained_tpgs = i.cum_var_rda1_rda2_constrained
), on = "region"]
cwm_rda_variance[, .(
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
for (region in names(cwm_matrices$X)) {
  corr_matrix <- cor(cwm_matrices$X[[region]][, c("max_log_tu",
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
cwm_rda_me <- marginal_effects_rda(cwm_rda)
tpg_rda_me <- marginal_effects_rda(tpg_rda)

# Check VIF
# Variance inflation values are all close to 1
vif_rda <- function(rda_results){
  vif <- lapply(rda_results, function(x) vif.cca(x)) 
  return(vif)
}
vif_rda(cwm_rda)
vif_rda(tpg_rda)

### Consistency & Specificity ----
# RDA relationship between X and Y
# For scaling=2
# Consistency represented as: stat sign of explanatory variables AND
# Correlation between response and explanatory variables, represented by 
# the angle between species scores and biplot scores 
# "The angles in the biplot between response and explanatory variables, 
# and between response variables themselves or explanatory variables themselves, 
# reflect their correlations.(Legrendre)"

# Function to calculate cosine similarity (angle-based correlation)
# Here we assess the directional consistency
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
cwm_rda_cosine <- cosine_similarity(cwm_rda) |>
  rbindlist(idcol = "region")
tpg_rda_cosine <- cosine_similarity(tpg_rda) |>
  rbindlist(idcol = "region")
rda_cosine <- rbindlist(list("traits" = cwm_rda_cosine, "tpg" = tpg_rda_cosine), idcol =
                          "type")
saveRDS(rda_cosine, file.path(path_cache, "rda_cosine.rds"))

# Consistency of significance
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
rda_me[p.value <= 0.05, .N, by = c("type", "env_factor")] |>
  _[order(type, -N), ] |> 
  dcast(...~type, value.var = "N") |> 
  fwrite(file.path(path_out, "rda_n_sign_marginal_effect.csv"))

rda_me[, significant := fifelse(p.value <= 0.05, "Y", "N")]
# saveRDS(rda_me, file.path(path_cache, "rda_me.rds"))


### Performance & Specificity----
# Marginal effect/variance of max_log_tu & the other env. factors

# Output plot of all marginal effects
# rda_me <- readRDS(file.path(path_cache, "rda_me.rds"))
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
  facet_grid(type~., 
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
  scale_alpha_manual(values = c("Y" = 1, "N" = 0.35))+ 
  guides(alpha = "none") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  labs(x = "Region",
       y = "Marginal effect",
       fill = "Environmental Factor")
ggsave(file.path(path_paper, "Figures", "prop_var_marginal_effect.png"))


# Output table 
rda_me_publ <- rda_me[, .(type, region, env_factor, Variance, p.value, prop_variance)]

# Modify EPT script -> include P+T+Riffles in the regressions
# Modify SPEAR -> include P+T+Riffles in the regressions

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