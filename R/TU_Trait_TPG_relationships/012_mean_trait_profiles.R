# __________________________________________________________________________________________________
# Mean trait profiles
# __________________________________________________________________________________________________

# Load tpgs & full trait data
trait_profile_groups <- readRDS(file.path(path_cache, "trait_profile_groups.rds"))
trait_matrix <- readRDS(file.path(path_cache, "trait_matrix_final.rds"))
trait_names <-
  grep("feed.*|resp.*|locom.*|size.*|volt.*",
       names(trait_matrix),
       value = TRUE)

# Calculation mean trait profiles ----
trait_matrix[trait_profile_groups, 
             group := group, 
             on = c("Region", "taxon")]
trait_matrix[, n_taxa_group := .N, by = c("Region", "group")]
trait_matrix_sub <- trait_matrix[, .SD, .SDcols = c("taxon",
                                                    trait_names,
                                                    "sensitivity_organic",
                                                    "group",
                                                    "n_taxa_group",
                                                    "Region")]
trait_matrix_sub <- melt(trait_matrix_sub,
                         id.vars = c("taxon",
                                     "group",
                                     "n_taxa_group",
                                     "Region"), 
                         value.name = "affinity",
                         variable.name = "trait")
mean_tps <- trait_matrix_sub[, .(mean_affinity = mean(affinity), 
                                 n_taxa_group), 
                 by = c("Region", "group", "trait")]
mean_tps <- unique(mean_tps)
mean_tps[, grouping_feature := sub("([a-z]{1,})(\\_)(.+)", "\\1", trait)]

# Mean trait profiles are already normalised!
# mean_tps[, sum(mean_affinity), by = c("Region", "group", "grouping_feature")] %>% 
#   .[grouping_feature %in% c("feed", "resp", "locom", "size", "volt"), ] %>% 
#   .[V1 > 1, ]

# Comparison between mean trait profiles of TPGs ----
# Check if there are any redundancies (i.e. if groups can be put together)
# mean_tps[Region == "California",] %>% 
#   ggplot(., aes(x = trait, y = mean_affinity)) +
#   geom_boxplot() +
#   coord_flip()
mean_tps[trait %in% trait_names & mean_affinity >= 0.5 , ]

# Does not include organic sensitivity
# Create an overview of defining traits (affinity of 0.5 or greater of mean trait profile)
mean_tps[mean_affinity >= 0.5, 
         defining_tc := paste(trait, collapse = ", "),
         by = c("Region", "group")]
saveRDS(mean_tps,
        file.path(path_cache, "mean_tps.rds"))

# Create an overview table
# TODO: how to call Qmd documents?
# Load Tables.RMD
# source(
#   "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/R/TU_Trait_TPG_relationships/Tables.Qmd"
# )

# TODO, create on overview on sensitivity_organic values per TPG
ggplot(trait_matrix_sub[trait == "sensitivity_organic",], 
       aes(x = as.factor(group), y = affinity)) +
  geom_violin() +
  stat_summary(fun = "median", color = "red") +
  stat_summary(fun = "mean", color = "forestgreen") +
  facet_wrap(.~Region,
             scales = "free") +
  theme_bw() +
  labs(x = "TPGs", 
       y = "S_org values")

