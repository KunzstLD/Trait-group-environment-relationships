# ______________________________________________________________________________
# Plot Dominant TPGs per region
# ______________________________________________________________________________

# Create Data.frame out of list
trait_groups_rel_final <- readRDS(file.path(path_cache, "trait_groups_rel_final.rds"))
tpgs_names <- grep("T[0-9]{1,}" ,
                   names(trait_groups_rel_final$family_lvl$Midwest),
                   value = TRUE)
trait_groups_rel_final$family_lvl <- lapply(trait_groups_rel_final$family_lvl, function(x) {
  setnames(x, tpgs_names, paste0("TPG", sub("T", "", tpgs_names), "_fam"), skip_absent = TRUE)
})
tpg_comb_family <- rbindlist(trait_groups_rel_final$family_lvl,
                             idcol = "region",
                             fill = TRUE)
tpg_comb_family[, STAID := NULL]
tpg_comb_family <- melt(
  tpg_comb_family,
  id.vars = c("max_log_tu", "region", "site"),
  variable.name = "tpg",
  value.name = "tpg_val"
) |>
  _[, region := fifelse(region == "PN", "Northwest", region)] |>
  _[, .(region, site, max_log_tu, tpg, tpg_val)]

# Create the ordered levels
tpg_levels <- paste0("TPG", 1:15, "_fam")
tpg_comb_family[, tpg := factor(tpg, levels = tpg_levels)]

# Create the ordered levels
tpg_levels <- paste0("TPG", 1:15, "_fam")
tpg_comb_family[, tpg := factor(tpg, levels = tpg_levels)]

# Create label mapping (1–15)
tpg_labels <- setNames(as.character(1:15), tpg_levels)

# Dominant TPGs per region based on abundance
tpg_comb_family[, .(mean_rel_abund = mean(tpg_val)), by = .(region, tpg)] |>
  _[, .(region, tpg, mean_rel_abund = mean_rel_abund * 100)] |>
  _[order(region, -mean_rel_abund), ] 

##  Plot ----
tpg_comb_family[, .(mean_rel_abund = mean(tpg_val)), by = .(region, tpg)] |>
  ggplot(aes(x = tpg, y = mean_rel_abund * 100)) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    width = 0.9,
    fill = "steelblue"
  ) +
  geom_text(aes(label = tpg_labels[tpg], y = 0),
            vjust = -0.3,
            size = 5) +
  facet_wrap( ~ region) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  labs(y = "Mean Rel. abundance (%)", x = "TPG (fam)")
ggsave(
  file.path(path_paper, "Figures", "tpgs_per_region.png"),
  width = 12,
  height = 6.5
)



tpg_comb_family[, .(mean_rel_abund = mean(tpg_val)), by = .(region, tpg)] |>
  _[order(region, -mean_rel_abund), ] |>
  _[tpg %in% c("TPG8_fam", "TPG10_fam", "TPG12_fam")] |>
  _[, .(min(mean_rel_abund), max(mean_rel_abund)), by = "tpg"]

