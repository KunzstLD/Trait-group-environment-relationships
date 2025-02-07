# ______________________________________________________________________________
# Relationship of the most important traits and TPGs based on the 
# importance ranking in the BRTs with max logTU
# Aim:
# - Testing the relationship of traits and TPGs with pesticide toxicity
# - Selecting a subset of most responsive traits: 
#     - Consistently most important, i.e. across multiple regions (at least two)
#     - Stat. significant Relationship in LM/GAM with pesticide toxicity across multiple regions
# ______________________________________________________________________________

# CWM ----
data_cwm_final <- readRDS(file.path(path_cache, "data_cwm_final.rds"))
most_imp_cwm <- readRDS(file.path(path_cache, "most_important_traits_cwm.rds"))
most_imp_cwm[, .N, by = "trait_TPG"] %>% 
.[order(-N)]

# select size large and the other consistently most important traits per region
cwm_combined <- rbindlist(data_cwm_final, idcol = "region", fill = TRUE) %>%
  .[trait %in% c(
    "size_large",
    "size_small",
    "sensitivity_organic",
    "feed_predator",
    "feed_parasite", # not responsive across multiple regions
    "feed_gatherer",
    "feed_filter",
    "resp_gil",
    "volt_bi_multi",
    "volt_semi",
    "locom_swim"
  ), ] %>%
  .[, region := fifelse(region == "PN", "Northwest", region)]

## LM or GAM modelling ----

# bivar. relationship
biv_plot <- list()
for(i in unique(cwm_combined$region)){
  for(j in unique(cwm_combined$trait)){
    biv_plot[[paste(i, "_", j)]] <- ggpairs(cwm_combined[region == i & trait==j, 
                                                         .(max_log_tu, cwm_val)],
                                            title = paste(i, "_", j)) 
  }
}
biv_plot

# Decide on a model for the most consistent responding traits
# Use GAM when GAM suggested (EDF > 1) and significant turning points
lm_cwm <- list()
sim_res_lm <- list()
gam_cwm <- list()
gam_assump <- list()
for (i in unique(cwm_combined$region)) {
  for(j in unique(cwm_combined$trait)){
    lm_cwm[[paste0(i, "_", j)]] <- lm(cwm_val ~ max_log_tu,
                                      data = cwm_combined[region == i & trait==j, ]
    )
    
    sim_res_lm[[paste0(i, "_", j)]] <- simulateResiduals(fittedModel = lm_cwm[[paste0(i, "_", j)]])
    
    gam_cwm[[paste0(i, "_", j)]] <- gam(
      cwm_val ~ s(max_log_tu),
      data = cwm_combined[region == i & trait==j, ],
      method = "REML"
    )
    
    gam_assump[[paste0(i, "_", j)]] <- gam.check(gam_cwm[[paste0(i, "_", j)]]) 
  }
}

# Diagnostics GAMs, visual inspection
par(mfrow = c(5,5))
for(i in names(gam_cwm)){
    plot(gam_cwm[[i]], main = i, residuals = TRUE, pch = 1, ylab = "CWM val")
}

# Extract EDF & significance of smooth term
edf_cwm <- extractEDF(gam_object=gam_cwm, idcol="region_cwm")
# saveRDS(edf_cwm, file.path(path_cache, "edf_cwm.rds"))

# Visual diagnostics for those CWM ~ TU traits relationships where LMs were suggested
lm_region_cwm_all <- edf_cwm[gam_or_lm=="lm", ][["region_cwm"]]
gam_region_cwm_all <- edf_cwm[gam_or_lm=="gam", ][["region_cwm"]]

for(i in names(sim_res_lm)[names(sim_res_lm) %in% lm_region_cwm_all]){
    plot(sim_res_lm[[i]], main = i)
}

# Direction of nonlinear effects (via derivatives)
gam_deriv_ls <- list()
for (i in names(gam_cwm)[names(gam_cwm) %in% gam_region_cwm_all]) {
  gam_deriv_ls[[i]] <- setDT(derivatives(gam_cwm[[i]]))
}
gam_deriv_trait <- rbindlist(gam_deriv_ls, idcol="region_trait")

# Assess based on median
gam_deriv_trait[, median_derivative := median(.derivative), by = "region_trait"]
# saveRDS(gam_deriv_trait, file.path(path_cache, "gam_deriv_trait.rds"))

# Save cwm combined with max_log_tu
# cwm_combined <- cwm_combined[, .(region, trait, cwm_val, max_log_tu)]
# saveRDS(cwm_combined, file.path(path_cache, "cwm_combined.rds"))
   
#_______________________________________________________________________________

# TPGs ----
## Family level ----
# Select consistently most responsive TPGs 
# TPG1_fam,TPG2_fam, TPG5_fam, TPG8_fam, TPG9_fam, TPG10_fam, and TPG12_fam
trait_groups_rel_final <- readRDS(file.path(path_cache, "trait_groups_rel_final.rds"))
tpgs_names <- grep("T[0-9]{1,}" ,names(trait_groups_rel_final$family_lvl$Midwest), value = TRUE)
trait_groups_rel_final$family_lvl <- lapply(
    trait_groups_rel_final$family_lvl,
    function(x) {
        setnames(
            x,
            tpgs_names,
            paste0("TPG", sub("T", "", tpgs_names), "_fam"),
            skip_absent = TRUE
        )
    }
)
tpg_comb_family <- rbindlist(trait_groups_rel_final$family_lvl, idcol = "region", fill = TRUE) %>%
    .[, .(
        region,
        TPG1_fam,
        TPG2_fam,
        TPG5_fam,
        TPG8_fam,
        TPG9_fam,
        TPG10_fam,
        TPG12_fam,
        TPG15_fam, # Added afterwards, not consistently most important, but fifth most important in Midwest
        max_log_tu
    )] %>%
    melt(.,
        id.vars = c("max_log_tu", "region"),
        variable.name = "tpg",
        value.name = "tpg_val"
    ) %>%
    .[, region := fifelse(region == "PN", "Northwest", region)]

## LM or GAM modelling ----
biv_plot_tpg_fam <- list()
for (i in unique(tpg_comb_family$region)) {
    biv_plot_tpg_fam[[i]] <- ggpairs(
        tpg_comb_family[
            region == i & tpg == "TPG1_fam",
            .(max_log_tu, tpg_val)
        ],
        title = i
    )
}
biv_plot_tpg_fam

# Choose a regression model
lm_tpg_fam <- list()
sim_res_lm_tpg_fam <- list()
gam_tpg_fam <- list()
gam_assump_tpg_fam <- list()
for (i in unique(tpg_comb_family$region)) {
    for (j in unique(tpg_comb_family$tpg)) {
        lm_tpg_fam[[paste0(i, "_", j)]] <- lm(tpg_val ~ max_log_tu,
            data = tpg_comb_family[region == i & tpg == j, ]
        )

        sim_res_lm_tpg_fam[[paste0(i, "_", j)]] <- simulateResiduals(fittedModel = lm_tpg_fam[[paste0(i, "_", j)]])

        gam_tpg_fam[[paste0(i, "_", j)]] <- gam(
            tpg_val ~ s(max_log_tu),
            data = tpg_comb_family[region == i & tpg == j, ],
            method = "REML"
        )

        gam_assump_tpg_fam[[paste0(i, "_", j)]] <- gam.check(gam_tpg_fam[[paste0(i, "_", j)]])
    }
}

# GAM Plots
par(mfrow = c(6,7))
for(i in names(gam_tpg_fam)){
  plot(gam_tpg_fam[[i]], main = i, residuals = TRUE, pch = 1, ylab = "TPG val")
}

# EDF & significance of smooth terms for decision between LM or GAM 
edf_tpg_fam <- extractEDF(gam_object=gam_tpg_fam, idcol="region_tpg_fam")
# saveRDS(edf_tpg_fam, file.path(path_cache, "edf_tpg_fam.rds"))

# LM Diagnostics
lm_region_tpg_fam_all <- edf_tpg_fam[gam_or_lm=="lm", ][["region_tpg_fam"]]
gam_region_tpg_fam_all <- edf_tpg_fam[gam_or_lm=="gam", ][["region_tpg_fam"]]
for(i in names(sim_res_lm_tpg_fam)[names(sim_res_lm_tpg_fam) %in% lm_region_tpg_fam_all]) {
  plot(sim_res_lm_tpg_fam[[i]], main = i)
}

# Direction of nonlinear effects (via derivatives)
gam_deriv_ls <- list()
for (i in names(gam_tpg_fam)[names(gam_tpg_fam) %in% gam_region_tpg_fam_all]) {
  gam_deriv_ls[[i]] <- setDT(derivatives(gam_tpg_fam[[i]]))
}
gam_deriv_tpg_fam <- rbindlist(gam_deriv_ls, idcol="region_tpg_fam")
gam_deriv_tpg_fam[, median_derivative := median(.derivative), by = "region_tpg_fam"]
# saveRDS(gam_deriv_tpg_fam, file.path(path_cache, "gam_deriv_tpg_fam.rds"))


## Genus level ----
# Most consistent TPGs on genus level: TPG4_genus, TPG10_genus, and TPG12_genus
# Also Added T9_genus and T15_genus for comparison, as they were the most important TPGs on genus level in some regions
trait_groups_rel_final$genus_lvl <- lapply(
    trait_groups_rel_final$genus_lvl,
    function(x) {
        setnames(
            x,
            tpgs_names,
            paste0(tpgs_names, "_genus"),
            skip_absent = TRUE
        )
    }
)
tpg_comb_genus <- rbindlist(trait_groups_rel_final$genus_lvl, idcol = "region", fill = TRUE) %>%
    .[, .(region, T4_genus, T10_genus, T9_genus, T12_genus, T15_genus, max_log_tu)] %>%
    melt(.,
        id.vars = c("max_log_tu", "region"),
        variable.name = "tpg",
        value.name = "tpg_val"
    ) %>%
    .[region == "PN", region := "Northwest"]

# Choose a model
lm_tpg_gen <- list()
sim_res_lm_tpg_gen <- list()
gam_tpg_gen <- list()
gam_assump_tpg_gen <- list()
for (i in unique(tpg_comb_genus$region)) {
    for (j in unique(tpg_comb_genus$tpg)) {
        lm_tpg_gen[[paste0(i, "_", j)]] <- lm(tpg_val ~ max_log_tu,
            data = tpg_comb_genus[region == i & tpg == j, ]
        )

        sim_res_lm_tpg_gen[[paste0(i, "_", j)]] <-
            simulateResiduals(fittedModel = lm_tpg_gen[[paste0(i, "_", j)]])

        tryCatch(
            {
                gam_tpg_gen[[paste0(i, "_", j)]] <- gam(
                    tpg_val ~ s(max_log_tu),
                    data = tpg_comb_genus[region == i &
                        tpg == j, ],
                    method = "REML"
                )
            },
            error = function(e) {
                cat("ERROR :", conditionMessage(e), "\n")
            }
        )

        tryCatch(
            {
                gam_assump_tpg_gen[[paste0(i, "_", j)]] <- gam.check(gam_tpg_gen[[paste0(i, "_", j)]])
            },
            error = function(e) {
                cat("ERROR :", conditionMessage(e), "\n")
            }
        )
    }
}

# Diagnostics
# GAM
par(mfrow = c(5,3))
for(i in names(gam_tpg_gen)){
  plot(gam_tpg_gen[[i]], main = i, residuals = TRUE, pch = 1, ylab = "TPG val (genus level)")
}

# EDF & significance of smooth terms
edf_tpg_gen <- extractEDF(gam_object=gam_tpg_gen, idcol="region_tpg_gen")
# saveRDS(edf_tpg_gen, file.path(path_cache, "edf_tpg_gen.rds"))

# Diagnostics LM
lm_region_tpg_gen_all <- edf_tpg_gen[gam_or_lm=="lm", ][["region_tpg_gen"]]
gam_region_tpg_gen_all <- edf_tpg_gen[gam_or_lm=="gam", ][["region_tpg_gen"]]
for(i in names(sim_res_lm_tpg_gen)[names(sim_res_lm_tpg_gen) %in% lm_region_tpg_gen_all]) {
  plot(sim_res_lm_tpg_gen[[i]], main = i)
}

# Direction of nonlinear effects (via derivatives)
gam_deriv_ls <- list()
for (i in names(gam_tpg_gen)[names(gam_tpg_gen) %in% gam_region_tpg_gen_all]) {
  gam_deriv_ls[[i]] <- setDT(derivatives(gam_tpg_gen[[i]]))
}
gam_deriv_tpg_gen <- rbindlist(gam_deriv_ls, idcol="region_tpg_gen")
gam_deriv_tpg_gen[, median_derivative := median(.derivative), by = "region_tpg_gen"]
# saveRDS(gam_deriv_tpg_gen, file.path(path_cache, "gam_deriv_tpg_gen.rds"))

# Save TPG combined (family and genus level) with max_log_tu
# tpg_comb_family <- tpg_comb_family[, .(region, tpg, tpg_val, max_log_tu)]
# tpg_comb <- rbindlist(
#     list(
#         "tpg_fam" = tpg_comb_family,
#         "tpg_gen" = tpg_comb_genus 
#     ),
#     use.names = TRUE,
#     idcol = "approach"
# )
# tpg_comb[, id := paste0(region, "_", tpg)]
# saveRDS(tpg_comb, file.path(path_cache, "tpg_combined.rds"))

#_______________________________________________________________________________

# Summary plots ----
# Load CWM & TPG data with max log TU
cwm_combined <- readRDS(file.path(path_cache, "cwm_combined.rds"))
edf_cwm <- readRDS(file.path(path_cache, "edf_cwm.rds"))
tpg_comb <- readRDS(file.path(path_cache, "tpg_combined.rds"))
tpg_comb[tpg %like% "_genus", tpg := sub("T", "TPG", tpg)]
edf_tpg_fam <- readRDS(file.path(path_cache, "edf_tpg_fam.rds"))
edf_tpg_gen <- readRDS(file.path(path_cache, "edf_tpg_gen.rds"))

## CWM plot ----
cwm_combined[, region_cwm := paste0(region, "_", trait)]

# Plot only for the most consistent CWM trait: size large
lm_region_cwm <- edf_cwm[region_cwm %like% ".+size_large" & gam_or_lm=="lm", ][["region_cwm"]]
gam_region_cwm <- edf_cwm[region_cwm %like% ".+size_large" & gam_or_lm=="gam", ][["region_cwm"]]

cwm_plot <- ggplot(cwm_combined[trait=="size_large", ], aes(x = max_log_tu, y = cwm_val)) +
    facet_wrap(.~as.factor(region)) +
    geom_point(alpha = 0.5, shape = 1, size = 2) +
    geom_smooth(
        data = ~ .x[region_cwm %in% lm_region_cwm, ],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    geom_smooth(
        data = ~ .x[region_cwm %in% gam_region_cwm, ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "steelblue"
    ) +
    labs(x = "Max logTU", 
         y = "CWM size large") +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title = element_text(size = 16, face = "bold"),
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
ggsave(file.path(path_paper, "Graphs", "gam_lm_cwm_toxicity.png"),
    width = 35,
    height = 20,
    units = "cm"
)

## TPG plots ----
lm_region_tpg_fam <- edf_tpg_fam[gam_or_lm=="lm",][["region_tpg_fam"]]
gam_region_tpg_fam <- edf_tpg_fam[gam_or_lm=="gam",][["region_tpg_fam"]]

### TPG fam plot ----
tpg_fam_plot <- ggplot(tpg_comb[approach=="tpg_fam", ], 
                   aes(x = max_log_tu, y = tpg_val)) +
    facet_grid(as.factor(tpg) ~ as.factor(region), scales="free_y") +  
    geom_point(alpha = 0.5, shape = 1, size = 2) +
    geom_smooth(
        data = ~ .x[id %in% lm_region_tpg_fam],
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        color = "steelblue"
    ) +
    geom_smooth(
        data = ~ .x[id %in% gam_region_tpg_fam, ],
        method = "gam",
        formula = y ~ s(x),
        se = TRUE,
        color = "steelblue"
    ) +
    labs(x = "Max logTU", 
         y = "CWM TPG") +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title = element_text(size = 16, face = "bold"),
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
        ),
        panel.spacing.x = unit(6, "mm")
    )
ggsave(file.path(path_paper, "Graphs", "gam_lm_tpg_fam_toxicity.png"),
    width = 38,
    height = 45,
    units = "cm"
)

### TPG genus plot ----
# Removed T9 and T15 again, as they were only responsive TPGs in one region
lm_region_tpg_gen <- edf_tpg_gen[gam_or_lm=="lm",][["region_tpg_gen"]]
gam_region_tpg_gen <- edf_tpg_gen[gam_or_lm=="gam",][["region_tpg_gen"]]

ggplot(tpg_comb[approach == "tpg_gen" & !tpg %in% c("TPG9_genus", "TPG15_genus", "TPG15_fam")], 
       aes(x = max_log_tu, y = tpg_val)) +
  facet_grid(as.factor(tpg) ~ as.factor(region), scales="free_y") +
  geom_point(alpha = 0.5, shape = 1, size = 2) +
  geom_smooth(
    data = ~ .x[id %in% c(lm_tpg_gen_ls), ],
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "steelblue"
  ) +
  geom_smooth(
    data = ~ .x[id %in% gam_tpg_gen_ls, ],
    method = "gam",
    formula = y ~ s(x),
    se = TRUE,
    color = "steelblue"
  ) +
  # scale_y_continuous(limits = c(0, 0.5),
  #                    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  labs(x = "Max logTU", 
       y = "CWM TPG") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16, face = "bold"),
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
    ),
    panel.spacing.x = unit(6, "mm")
  )
ggsave(file.path(path_paper, "Graphs", "gam_lm_tpg_genus_toxicity.png"),
       width = 35,
       height = 21,
       units = "cm"
)

# Summary table ----

## CWM table ----
lm_region_cwm_all <- edf_cwm[gam_or_lm=="lm",][["region_cwm"]]
gam_region_cwm_all <- edf_cwm[gam_or_lm=="gam",][["region_cwm"]]

cwm_table <- rbind(
    lapply(
        lm_cwm[lm_region_cwm_all],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., idcol = "Region") %>% 
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_cwm[gam_region_cwm_all],
            function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., idcol = "Region")%>% 
            setnames(., "statistic", "t"),
        lapply(
            gam_cwm[gam_region_cwm_all],
            function(x) {
                broom::tidy(x, parametric = FALSE)
            }
        ) %>%
            rbindlist(., idcol = "Region") %>%
            setnames(., 
            c("statistic", "p.value"), 
            c("F", "p.value_gam")),
        fill = TRUE
    ),
    fill = TRUE
)

# Merge derivatives from GAM
cwm_table[unique(gam_deriv_trait[, .(region_trait, median_derivative)]), median_derivative := i.median_derivative, on =
            c(Region = "region_trait")]

# Apply the Benjamini Hochberg correction across all p-values (GAM & LM)
# - Ranks p-values and adjusts them according to their rank (
# - min(m/rank * p.value, 1), whereby m denotes the total number of hypothesis tested
cwm_table[, p_p_gam := coalesce(p.value, p.value_gam)]
cwm_table[, p_p_gam_BH := p.adjust(p_p_gam, method="BH")]

# R2 & deviance
cwm_r2 <- lapply(lm_cwm[lm_region_cwm_all], function(x) {
  broom::glance(x)
}) |>
  rbindlist(, idcol = 'Region') |>
  _[, c('Region', 'r.squared')] |> 
  _[, r.squared := round(r.squared, digits=3)]

cwm_dev <- lapply(gam_cwm[gam_region_cwm_all], function(x) {
  data.table('dev_explained' = summary(x)$dev.expl)
}) |>
  rbindlist(, idcol = 'Region') |> 
  _[, dev_explained := round(dev_explained, digits=3)]

# Merge R2 and deviance
cwm_table[cwm_r2, r_squared := r.squared, on = 'Region']
cwm_table[cwm_dev, dev_explained := dev_explained, on = 'Region']
cwm_table[,r_squared_dev_expl := coalesce(r_squared, dev_explained)]
cwm_table[,p_p_gam := coalesce(p.value, p.value_gam)]
cwm_table[, cwm_trait := sub("(California|Midwest|Northeast|Northwest|Southeast)_(.+)", "\\2", Region)]
cwm_table[, Region := sub("(California|Midwest|Northeast|Northwest|Southeast)_(.+)", "\\1", Region)]
saveRDS(cwm_table, file.path(path_cache, "cwm_tabl_publ.rds"))

# Polish for SI
# cwm_table[, `:=`(
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
#   )
# )]


## TPG table ----

### Family level ----
tpg_table_fam <- rbind(
    lapply(
        lm_tpg_fam[
            names(lm_tpg_fam) %in% lm_region_tpg_fam],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., idcol = "Region_TPG") %>% 
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_tpg_fam[
                names(gam_tpg_fam) %in% gam_region_tpg_fam
            ], function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., idcol = "Region_TPG") %>% 
            setnames(., "statistic", "t"),
        lapply(
            gam_tpg_fam[
                names(gam_tpg_fam) %in% gam_region_tpg_fam], function(x) {
                broom::tidy(x, parametric = FALSE)
            }
        ) %>%
            rbindlist(., idcol = "Region_TPG") %>% 
            setnames(., 
            c("statistic", "p.value"), 
            c("F", "p.value_gam")),
        fill = TRUE
    ) %>% 
    .[order(Region_TPG), ],
    fill = TRUE
) 

# Merge derivatives from GAM
tpg_table_fam[unique(gam_deriv_tpg_fam[, .(region_tpg_fam, median_derivative)]), median_derivative := i.median_derivative, on =
                c(Region_TPG = "region_tpg_fam")]

# Apply the Benjamini Hochberg correction across all p-values (GAM & LM)
tpg_table_fam[, p_p_gam := coalesce(p.value, p.value_gam)]
tpg_table_fam[, p_p_gam_BH := p.adjust(p_p_gam, method="BH")]

# Explained variance/deviance
tpg_fam_r2 <- lapply(lm_tpg_fam[names(lm_tpg_fam) %in% lm_region_tpg_fam], function(x) {
  broom::glance(x)
}) |>
  rbindlist(, idcol = 'Region_TPG') |>
  _[, c('Region_TPG', 'r.squared')] |>
  _[, r.squared := round(r.squared, digits = 3)]

tpg_fam_dev <- lapply(gam_tpg_fam[names(gam_tpg_fam) %in% gam_region_tpg_fam], function(x) {
  data.table('dev_explained' = summary(x)$dev.expl)
}) |>
  rbindlist(, idcol = 'Region_TPG') |>
  _[, dev_explained := round(dev_explained, digits = 3)]

tpg_table_fam[tpg_fam_r2, r_squared := r.squared, on = 'Region_TPG']
tpg_table_fam[tpg_fam_dev, dev_explained := dev_explained, on = 'Region_TPG']
tpg_table_fam[,r_squared_dev_expl := coalesce(r_squared, dev_explained)]
tpg_table_fam[,p_p_gam := coalesce(p.value, p.value_gam)]
tpg_table_fam[, TPG := sub("(.+)_(.+_.+)", "\\2", Region_TPG)]
tpg_table_fam[, Region := sub("(.+)_(.+_.+)", "\\1", Region_TPG)]
saveRDS(tpg_table_fam, file.path(path_cache, "tpg_tabl_publ.rds"))

# Prepare for SI
# tpg_table_fam[, `:=`(
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
#   )
# )]


### Genus level ----
tpg_table_gen <- rbind(
    lapply(
        lm_tpg_gen[
            names(lm_tpg_gen) %in% lm_region_tpg_gen
        ],
        function(x) {
            broom::tidy(x)
        }
    ) %>%
        rbindlist(., "Region_TPG",
            use.names = TRUE,
            fill = TRUE
        ) %>%
        setnames(., "statistic", "t"),
    rbind(
        lapply(
            gam_tpg_gen[
                names(gam_tpg_gen) %in% gam_region_tpg_gen
            ],
            function(x) {
                broom::tidy(x, parametric = TRUE)
            }
        ) %>%
            rbindlist(., "Region_TPG",
                use.names = TRUE,
                fill = TRUE
            ) %>%
            setnames(
                .,
                c("statistic", "p.value"),
                c("F", "p.value_gam")
            ),
        lapply(
            gam_tpg_gen[
                names(gam_tpg_gen) %in% gam_region_tpg_gen
            ],
            function(x) {
                broom::tidy(x, parametric = FALSE)
            }
        ) %>%
            rbindlist(., "Region_TPG",
                use.names = TRUE,
                fill = TRUE
            ) %>%
            setnames(
                .,
                c("statistic", "p.value"),
                c("F", "p.value_gam")
            ),
        fill = TRUE
    ) %>%  
    .[order(Region_TPG), ],
    fill = TRUE
)

# Merge derivatives from GAM
tpg_table_gen[unique(gam_deriv_tpg_gen[, .(region_tpg_gen, median_derivative)]), median_derivative := i.median_derivative, on =
                c(Region_TPG = "region_tpg_gen")]

# Apply the Benjamini Hochberg correction across all p-values (GAM & LM)
tpg_table_gen[, p_p_gam := coalesce(p.value, p.value_gam)]
tpg_table_gen[, p_p_gam_BH := p.adjust(p_p_gam, method="BH")]

# Prepare for SI
tpg_table_gen[, `:=`(
    estimate = round(estimate, digits = 2),
    std.error = round(std.error, digits = 2),
    t = round(t, digits = 2),
    F = round(F, digits = 2),
    p.value = fifelse(
        round(p.value, digits = 3) == 0,
        "<0.01",
        paste(round(p.value, digits = 3))
    ),
    p.value_gam = fifelse(
        round(p.value_gam, digits = 3) == 0,
        "<0.01",
        paste(round(p.value_gam, digits = 3))
    )
)]

# Explained Variance/Deviance
tpg_gen_r2 <- lapply(lm_tpg_gen[names(lm_tpg_gen) %in% lm_region_tpg_gen], function(x) {
  broom::glance(x)
}) |>
  rbindlist(, idcol = 'Region_TPG') |>
  _[, c('Region_TPG', 'r.squared')] |>
  _[, r.squared := round(r.squared, digits = 3)]

tpg_gen_dev <- lapply(gam_tpg_gen[names(gam_tpg_gen) %in% gam_region_tpg_gen], function(x) {
  data.table('dev_explained' = summary(x)$dev.expl)
}) |>
  rbindlist(, idcol = 'Region_TPG') |>
  _[, dev_explained := round(dev_explained, digits = 3)]

# Merge R2 and dev
tpg_table_gen[tpg_gen_r2, r_squared := r.squared, on = "Region_TPG"]
tpg_table_gen[tpg_gen_dev, dev_explained := dev_explained, on = "Region_TPG"]
tpg_table_gen[,r_squared_dev_expl := coalesce(r_squared, dev_explained)]
tpg_table_gen[,p_p_gam := coalesce(p.value, p.value_gam)]
tpg_table_gen[, TPG := sub("(.+)_(.+_.+)", "\\2", Region_TPG)]
tpg_table_gen[, Region := sub("(.+)_(.+_.+)", "\\1", Region_TPG)]
tpg_table_gen[, TPG := sub("^T", "TPG", TPG)]
tpg_table_gen[Region=="Southeast" & TPG=="TPG12_genus" & term == "(Intercept)", ]
saveRDS(tpg_table_gen, file.path(path_cache, "tpg_tabl_genus_publ.rds"))



# Tables for publication ----
## CWM ----
cwm_table <- readRDS(file.path(path_cache, "cwm_tabl_publ.rds"))
cwm_table[, .(
  Region,
  cwm_trait,
  term,
  estimate,
  std.error,
  t, 
  p_p_gam,
  edf = round(edf, digits=2),
  F,
  r_squared_dev_expl
)] |>
  _[cwm_trait %in%  c(
    "size_large",
    "size_small",
    "feed_predator",
    "feed_gatherer",
    "feed_filter",
    "resp_gil",
    "volt_bi_multi",
    "volt_semi",
    "bi_multi",
    "locom_swim",
    "sensitivity_organic"
  ), ] |> 
  _[order(Region, cwm_trait), ] |> 
  fwrite(file.path(path_paper, "Tables", "CWM_log_tu_results.csv"))

## TPGs ----
tpg_table_fam <- readRDS(file.path(path_cache, "tpg_tabl_publ.rds"))
tpg_table_genus <- readRDS(file.path(path_cache, "tpg_tabl_genus_publ.rds"))

tpg_table_final <- rbind(
    tpg_table_fam[, c(
        "Region",
        "TPG",
        "term",
        "estimate",
        "std.error",
        "t",
        "p_p_gam",
        "edf",
        "F",
        "r_squared_dev_expl"
    )],
    tpg_table_genus[, c(
        "Region",
        "TPG",
        "term",
        "estimate",
        "std.error",
        "t",
        "p_p_gam",
        "edf",
        "F",
        "r_squared_dev_expl"
    )]
) 
tpg_table_final[, Taxonomic_resolution := ifelse(grepl("_fam", TPG), "family_level", "genus_level")]
tpg_table_final = tpg_table_final[order(Taxonomic_resolution, Region), ]
tpg_table_final[, edf := round(edf, digits = 2)] %>%
  .[, .SD, .SDcols = !names(tpg_table_final) %in% "Taxonomic_resolution"] %>%
    fwrite(., file.path(path_paper, "Tables", "TPG_log_tu_results.csv"))