# ____________________________________________
# Mean trait profiles
# ____________________________________________

# Load tpgs
tpgs <- readRDS(file.path(path_cache, "tpgs_family_genus.rds"))

# Load trait data
trait_genera <- readRDS(file.path(path_cache, "trait_genera.rds"))
trait_family <- readRDS(file.path(path_cache, "trait_family.rds"))
trait_names <-
  grep("feed.*|resp.*|locom.*|size.*|volt.*|sensitivity.*",
       names(trait_genera),
       value = TRUE)

# Calculation mean trait profiles ----
# Genus level
trait_genera[tpgs, group := i.group,
             on = c("genus" = "taxon")]
trait_genera[, n_taxa_group := .N, by = "group"]
# saveRDS(trait_genera, file.path(path_cache, "trait_genera_tpg.rds"))
mtps_genera <- calc_mean_tps(
  x = trait_genera[, .SD, .SDcols = !c("family", "order")],
  taxa_id = "genus"
)

# Overview on orders 
# trait_genera[, .N, by = "order"] %>% 
# .[order(-N), ] %>% 
# fwrite(., file.path(path_out, "overview_order_tpg_genus.csv"))

# Family level
trait_family[tpgs, group := i.group,
             on = c("family" = "taxon")]
trait_family[, n_taxa_group := .N, by = "group"]
# saveRDS(trait_family, file.path(path_cache, "trait_family_tpg.rds"))
mtps_family <- calc_mean_tps(
  x = trait_family[, .SD, .SDcols = !c("order")],
  taxa_id = "family"
)

# Overview over orders
trait_family[, .N, by = "order"] %>%
  .[order(-N), ] %>%
  fwrite(., file.path(path_out, "overview_order_tpg_family.csv"))

# Normalisation
# Most mean trait profiles are normalised
# Still, for some there seems to be problems.
# Hence, normalisation is done again
# mtps_genera[, sum(mean_affinity), by = c("group", "grouping_feature")] %>% 
#    .[V1 != 1, ]
mtps_genera[trait %like% "feed.*|resp.*|locom.*|size.*|volt.*",
  mean_affinity := mean_affinity / sum(mean_affinity),
  by = c("group", "grouping_feature")
]
mtps_family[trait %like% "feed.*|resp.*|locom.*|size.*|volt.*",
  mean_affinity := mean_affinity / sum(mean_affinity),
  by = c("group", "grouping_feature")
]

# Comparison between mean trait profiles of TPGs ----
mtps_genera %>% 
  ggplot(., aes(x = trait, y = mean_affinity)) +
  geom_boxplot() +
  coord_flip()
saveRDS(mtps_genera,
        file.path(path_cache, "mtps_genera.rds"))

mtps_family %>% 
  ggplot(., aes(x = trait, y = mean_affinity)) +
  geom_boxplot() +
  coord_flip()
saveRDS(mtps_family,
        file.path(path_cache, "mtps_family.rds"))

# Create an overview table
# TODO: how to call Qmd documents?
# Load Tables.Qmd
# source("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/R/TU_Trait_TPG_relationships/Tables.Qmd")

# Closer look at groups ----
mtps_family[mean_affinity >= 0.5 & grouping_feature == "feed", .N, by = "trait"]

mtps_family[mean_affinity >= 0.5, ] %>%
  .[order(group), ] %>% 
  .[, .(group, trait, grouping_feature)] %>% 
  .[, group := as.character(paste0("TPG_", group))] %>% 
  dcast(., ...~ trait, fun.aggregate = length)  %>% 
  .[feed_filter == 0 & feed_herbivore == 0 & feed_predator == 0 & feed_gatherer == 0 & feed_shredder == 0, ]
