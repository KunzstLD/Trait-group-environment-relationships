# _______________________________________________________________
# Functional groups ----
# Clustering the genera and family-level taxa in the 
# North American trait databases
## ______________________________________________________________
trait_genera <- readRDS(file.path(path_cache, "trait_genera.rds"))
trait_family <- readRDS(file.path(path_cache, "trait_family.rds"))

fc_traits <-
  grep("feed.*|resp.*|locom.*|size.*|volt.*",
       names(trait_genera),
       value = TRUE)
q_traits <- "sensitivity_organic"

# Cluster analysis ----
results_genera_cl <- clustering_traits(
  x = trait_genera[, .SD, .SDcols = c(fc_traits, q_traits)],
  fc_traits = fc_traits,
  q_traits = q_traits,
  taxa = trait_genera$genus
)
results_family_cl <- clustering_traits(
  x = trait_family[, .SD, .SDcols = c(fc_traits, q_traits)],
  fc_traits = fc_traits,
  q_traits = q_traits,
  taxa = trait_family$family
)
results_cl <- list(
  "family" = results_family_cl,
  "genus" = results_genera_cl
)
saveRDS(results_cl, file.path(path_cache, "results_cl"))

# Cophenetic distance ----
# Relatively low for both (0.55, 0.56)
cor(results_genera_cl$distance_matrix,
    cophenetic(results_genera_cl$hc_wardD2))
cor(results_family_cl$distance_matrix,
    cophenetic(results_family_cl$hc_wardD2))

# Dendrograms ----
# Genera-lvl: 
#  - Gap statistic and Silhouette width suggest 20 groups for genera
#  - Would probably more if kmax would be higher
#  - TODO: How to decide -> based on Dendrogram?
#  - Could try different numbers: 10, 12, ...?
# Family-lvl:
#  - Gap statistic: 15, silhouette: 20
#  - TODO Check trait profiles
dendrograms <- list()
for (i in names(results_cl)) {
  plot <- fun_dendrog_pl(
    hc = results_cl[[i]]$hc_wardD2,
    optimal_nog = results_cl[[i]]$gap_statistic,
    labels = results_cl[[i]]$labels_dendrogram
  )
  dendrograms[[i]] <- plot
  png(
    file = file.path(path_out,
                     "Graphs",
                     paste0(
                       "Dendrogram_", i, ".png"
                     )),
    width = 1100,
    height = 1400,
    res = 100
  )
  plot(dendrograms[[i]], horiz = TRUE)
  dev.off()
}
plot(dendrograms$family, horiz = TRUE)

# TPGs ----
tpgs <- list()
for(i in names(results_cl)) {
  tpgs[[i]] <-
    data.table(
      taxon = names(
        cutree(
          dendrograms[[i]],
          k = results_cl[[i]]$gap_statistic,
          order_clusters_as_data = FALSE
        )
      ),
      group = cutree(
        dendrograms[[i]],
        k = results_cl[[i]]$gap_statistic,
        order_clusters_as_data = FALSE
      )
    )
}
tpgs <- rbindlist(tpgs, idcol = "dataset")
saveRDS(tpgs, file.path(path_cache, "tpgs_family_genus.rds"))