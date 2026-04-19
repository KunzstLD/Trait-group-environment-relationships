# __________________________________________________________________________________________________
# Relationship between fraction of EPT taxa, estimated toxicity and other environmental factors ----
# __________________________________________________________________________________________________

## Data preparation ----
# We load abundances and toxicity data and calculate the fraction of EPT taxa
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))
abund <- abund[, .(
    Region,
    site, 
    STAID,
    taxon,
    species,
    genus,
    family,
    order,
    taxonomic_level,
    abundance)]

# Calc. abundance per site, but merge chemical information per STAID
# Midwest has more sites than STAID, e.g. T391114085205801
abund[Region == "Midwest", STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
abund[Region == "Midwest", STAID := paste0("T", STAID)]

# Few sites have two IDs (for eco and wq)
# because chemical and eco information were taken from slightly
# different positions (mainly because of accessability issues)
abund[STAID == "T03611200", STAID := "T03611100"]
abund[site == "T12073525", site := "T12073425"]
# save this version of the abundance data with improved STAID labels 
# and the site IDs changed for later use
# saveRDS(abund, file.path(path_cache, "total_abund_CEOPT_corrected.rds"))

# Toxicity
max_tu <- readRDS(file.path(path_cache, "max_tu.rds"))
setnames(max_tu, "TSITE_NO_WQ", "site")

# Calculate Fraction EPT of the whole assemblage per site
# apply sqrt transformation
abund[, abundance := sqrt(abundance)]
abund[, abund_order_site := sum(abundance), by = c("site", "order")]
abund[, total_abund := sum(abundance), by = "site"]
abund[order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera"),
    abund_EPT := sum(abundance),
    by = "site"
]
abund_subs <- unique(abund[
    !is.na(abund_EPT),
    .(Region, site, STAID, abund_EPT, total_abund)
])
abund_subs[, frac_EPT := abund_EPT/total_abund]
# Hmisc::describe(abund_subs$frac_EPT)

# Merge toxicity
abund_subs[max_tu[Region != "Midwest", ],
    max_log_tu := i.max_log_tu,
    on = "site"
]
abund_subs[max_tu[Region == "Midwest"],
    max_log_tu := i.max_log_tu,
    on = c("STAID" = "site")
]
abund_subs[Region == "PN", Region := "Northwest"]

# Few sites don't have chemical information
abund_subs <- abund_subs[!is.na(max_log_tu), ]
# Few STAIDs are duplicates (same chemical info, but different ecology)
# abund_subs[!is.na(STAID) & duplicated(STAID),]

# Merge environmental data
data_cwm_env <- readRDS(file.path(path_cache, "data_cwm_env.rds"))
names(data_cwm_env)[4] <- "Northwest"
env_data <- rbindlist(data_cwm_env, idcol = "Region", fill = TRUE)[, .(Region, site, Riffle.FRC, Temp.median, orthoP_4wk.median)] |>
  unique()

abund_subs[env_data, `:=`(
  Riffle.FRC = i.Riffle.FRC,
  Temp.median = i.Temp.median,
  orthoP_4wk.median = i.orthoP_4wk.median
), on = c("Region", "site")]

# Final EPT data
# saveRDS(abund_subs, file.path(path_cache, "ept.rds"))

# Scatterplot with maxTU
ggplot(
  abund_subs, aes(x = max_log_tu, y = frac_EPT)) +
  geom_point() +
  facet_wrap(~ Region) +
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
  labs(
    y = "%EPT taxa",
    x = "maxTU"
    )
ggsave(
  file.path(path_paper, "Figures", "ept_maxTU.png"),
  width = 12,
  height = 6.5
)

# ----------------------------------
# Multiple linear regression ----
# ----------------------------------
lm_ept <- list()
sim_res_ept <- list()
for (i in unique(abund_subs$Region)) {
  lm_ept[[i]] <- lm(frac_EPT ~ max_log_tu + Riffle.FRC + Temp.median + orthoP_4wk.median,
                    data = abund_subs[Region == i, ])
  
  sim_res_ept[[i]] <- simulateResiduals(fittedModel = lm_ept[[i]])
}

# Check Diagnostics
for(i in names(sim_res_ept)){
  plot(sim_res_ept[[i]], main = i)
}

## Summary table ----
ept_table <-lapply(
  lm_ept,
  function(x) {
    broom::tidy(x)
  }
) |> 
  rbindlist(idcol = "Region")

# Add R2 and deviance 
ept_gof <- lapply(lm_ept, broom::glance) |> rbindlist(idcol = "Region")
ept_table[ept_gof, `:=`(r.squared = i.r.squared, adj.r.squared = i.adj.r.squared), on = "Region"]

## Type II Sum of squares ----
# Test each main effect after accounting for the other main effects (assuming no interactions)
ept_type2ss <- lapply(lm_ept, function(x) car::Anova(x, type=2, test="F"))
ept_type2ss <- lapply(ept_type2ss, broom::tidy) |> rbindlist(idcol="Region")

# Calc. prop. explained variance -> Performance measure
ept_type2ss[, total_var := sum(sumsq), by = "Region"]
ept_type2ss[, prop_var := (sumsq/total_var)*100]
saveRDS(ept_type2ss, file.path(path_cache, "ept_type2ss.rds"))

# Calc. N of significant marginal effects -> consistency of significance
ept_type2ss[p.value <= 0.05, .N, by="term"]

# Average across regions
ept_type2ss[, mean(prop_var), by = "term"]

## Check the direction of the effect (parameter level) ----
ept_table[, dir_estimate := sign(estimate)]

# Positive: only Riffle Frc. in two regions 
ept_table[p.value <=0.05 & dir_estimate==1, sum(dir_estimate), by="term"]
# Negative: max_log_tu in 4 regions and temperature in 2 regions
ept_table[p.value <=0.05 & dir_estimate==-1, sum(dir_estimate), by="term"]

## Output for SI ----
ept_table[, `:=`(
  estimate = round(estimate, digits = 2),
  std.error = round(std.error, digits = 2),
  statistic = round(statistic, digits = 2),
  p.value = fifelse(
    round(p.value, digits = 3) == 0,
    "<0.01",
    paste(round(p.value, digits = 3))
  ),
  r.squared = round(r.squared, digits = 2),
  adj.r.squared = round(adj.r.squared, digits = 2)
)]
fwrite(ept_table, file.path(path_paper, "Tables", "EPT_regr_results.csv"))

# Create Publication table ----
spear_type2ss <- readRDS(file.path(path_cache, "spear_type2ss.rds"))
rda_me <- readRDS(file.path(path_cache, "rda_me.rds"))

rda_me_subset <- rda_me[env_factor =="max_log_tu", .(type, region, p.value, prop_variance)]
setnames(rda_me_subset, "prop_variance", "prop_var")

ept_subset <- ept_type2ss[term == "max_log_tu", .(Region, p.value, prop_var)]
setnames(ept_subset, "Region", "region")
ept_subset[, type := "%EPT"]

spear_subset <- spear_type2ss[term == "max_log_tu", .(Region, p.value, prop_var)]
setnames(spear_subset, "Region", "region")
spear_subset[, type := "SPEAR"]

rbind(rda_me_subset, ept_subset, spear_subset) |> 
  _[p.value > 0.05, ]

rbind(rda_me_subset, ept_subset, spear_subset) |>
  dcast(region ~ type, value.var = "prop_var") |>
  _[, .(
    region,
    TPG = round(tpg, 2),
    Traits = round(traits, 2),
    SPEAR = round(SPEAR, 2),
    `%EPT` = round(`%EPT`, 2)
  )] |> 
  fwrite(file.path(path_paper, "Tables", "Comp_tpg_trait_indicators.csv"))

# ______________________________________________________________________________
# GAMs & LMs with pesticide toxicity ----
# lm_ept <- list()
# sim_res_ept <- list()
# gam_ept <- list()
# gam_assump_ept <- list()
# for (i in unique(abund_subs$Region)) {
#     lm_ept[[i]] <- lm(frac_EPT ~ max_log_tu,
#         data = abund_subs[Region == i, ]
#     )
# 
#     sim_res_ept[[i]] <- simulateResiduals(fittedModel = lm_ept[[i]])
# 
#     gam_ept[[i]] <- gam(
#         frac_EPT ~ s(max_log_tu),
#         data = abund_subs[Region == i, ],
#         method = "REML"
#     )
# 
#     gam_assump_ept[[i]] <- gam.check(gam_ept[[i]])
# }
# 
# ## Diagnostics ----
# # Extract EDF & significance of smooth term
# edf_ept <- extractEDF(gam_object=gam_ept, idcol="Region")
# 
# # GAM & LM plots
# par(mfrow = c(2,3))
# for(i in names(gam_ept)){
#     plot(gam_ept[[i]], main = i, residuals = TRUE, pch = 1, ylab = "frac_EPT")
# }
# 
# for(i in names(sim_res_ept)){
#     plot(sim_res_ept[[i]], main = i)
# }

#_______________________________________________________________________________
# Plotting ----
# gam_region_ept <- edf_ept[gam_or_lm == "gam", ][["Region"]]
# lm_region_ept <- edf_ept[gam_or_lm == "lm", ][["Region"]]
# 
# ggplot(abund_subs, aes(x = max_log_tu, y = frac_EPT)) +
#     facet_wrap(as.factor(Region) ~.) +
#     geom_point(alpha = 0.5, shape = 1, size = 2) +
#     geom_smooth(
#         data = ~ .x[Region %in% lm_region_ept, ],
#         method = "lm",
#         formula = y ~ x,
#         se = TRUE,
#         color = "steelblue"
#     ) +
#     labs(
#         x = "Max logTU",
#         y = "Fraction EPT taxa per site"
#     ) + 
#     geom_smooth(
#         data = ~ .x[Region %in% gam_region_ept, ],
#         method = "gam",
#         formula = y ~ s(x),
#         se = TRUE,
#         color = "steelblue"
#         ) +
#     theme_bw() +
#     theme(
#         legend.position = "none",
#         axis.title = element_text(size = 16, face = "bold"),
#         axis.text.x = element_text(
#             family = "Roboto Mono",
#             size = 16
#         ),
#         axis.text.y = element_text(
#             family = "Roboto Mono",
#             size = 16
#         ),
#         strip.text = element_text(
#             family = "Roboto Mono",
#             size = 16
#         )
#     )
# ggsave(file.path(path_paper, "Graphs", "EPT_toxicity.png"),
#     width = 35,
#     height = 20,
#     units = "cm"
# )
# 
# # Summary Table EPT ----
# ept_table <- rbind(
#     lapply(
#         lm_ept[lm_region_ept],
#         function(x) {
#             broom::tidy(x)
#         }
#     ) %>%
#         rbindlist(., idcol = "Region") %>%
#         setnames(., "statistic", "t"),
#     rbind(
#         lapply(
#             gam_ept[gam_region_ept],
#             function(x) {
#                 broom::tidy(x, parametric = TRUE)
#             }
#         ) %>%
#             rbindlist(., idcol = "Region") %>%
#             setnames(., "statistic", "t"),
#         lapply(
#             gam_ept[gam_region_ept],
#             function(x) {
#                 broom::tidy(x, parametric = FALSE)
#             }
#         ) %>%
#             rbindlist(., idcol = "Region") %>%
#             setnames(
#                 .,
#                 c("statistic", "p.value"),
#                 c("F", "p.value_gam")
#             ),
#         fill = TRUE
#     ),
#     fill = TRUE
# )
# ept_table[, `:=`(
#   estimate = round(estimate, digits = 2),
#   std.error = round(std.error, digits = 2),
#   t = round(t, digits = 2),
#   F = round(F, digits = 2),
#   p.value = fifelse(
#     round(p.value, digits = 3) == 0,
#     "<0.01",
#     paste(round(p.value, digits = 3))
#   ),
#   p.value_gam = fifelse(
#     round(p.value_gam, digits = 3) == 0,
#     "<0.01",
#     paste(round(p.value_gam, digits = 3))
#   ), 
#   edf = round(edf, digits = 2),
#   ref.df = round(ref.df, digits = 2)
# )]
# 
# # Add R2 and deviance 
# ept_gof <- lapply(lm_ept[lm_region_ept], broom::glance) %>%
#     rbindlist(., idcol = "Region")
# ept_dev <- lapply(gam_ept[gam_region_ept], function(x) {
#   data.table("dev_explained" = summary(x)$dev.expl)
# }) |>
#   rbindlist(, idcol = "Region") |> 
#   _[, dev_explained := round(dev_explained, digits=3)]
# 
# # Merge R2 and deviance
# ept_table[ept_gof,
#           r.squared := i.r.squared,
#           on = "Region"
# ]
# ept_table[ept_dev, dev_explained := i.dev_explained, on="Region"]
# ept_table[,r_squared_dev_expl := coalesce(r.squared, dev_explained)]
# ept_table[, r_squared_dev_expl := round(r_squared_dev_expl, digits = 3)]
# ept_table[,p_p_gam := coalesce(p.value, p.value_gam)]
# ept_table[term != "(Intercept)", ]
# 
# setcolorder(
#     ept_table,
#     c("Region", "term", "estimate", "std.error", "t", "p_p_gam", "r_squared_dev_expl")
# )
# fwrite(ept_table, file.path(path_paper, "Tables", "EPT_log_tu_results.csv"))
