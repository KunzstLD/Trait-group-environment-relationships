# __________________________________________________________________________________________________
# Dendrograms & trait profile groups ----
# TODO: New strategy: use the optimal number of groups from
# Gap statistic + look on the dendrogram for subclusters!
# __________________________________________________________________________________________________

# Results clustering
results_cl <- readRDS(file.path(path_cache, "Results_clustering.rds"))

# Optimal number of groups according to the gap statistic
# Manual assignment (can be changed later), for now lowest number of groups 
gap <- readRDS(file.path(path_cache, "Gap_statistic_results.rds"))

results_cl$California$optimal_nog <- min(unique(gap$California)) #12
results_cl$Midwest$optimal_nog <- unique(gap$Midwest) # 15
results_cl$Northeast$optimal_nog <- min(unique(gap$Northeast)) # 13
results_cl$PN$optimal_nog <- min(unique(gap$PN)) # 10
results_cl$Southeast$optimal_nog <- unique(gap$Southeast) #15

# Dendrograms ----
# Plot and save dendrograms
dendrograms <- list()
for (i in names(results_cl)) {
  plot <- fun_dendrog_pl(
    hc = results_cl[[i]]$hc_wardD2,
    optimal_nog = results_cl[[i]]$optimal_nog,
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

# TPGs ----
trait_profile_groups <- list()
for(i in names(results_cl)) {
  trait_profile_groups[[i]] <-
    data.table(
      taxon = names(
        cutree(
          dendrograms[[i]],
          k = results_cl[[i]]$optimal_nog,
          order_clusters_as_data = FALSE
        )
      ),
      group = cutree(
        dendrograms[[i]],
        k = results_cl[[i]]$optimal_nog,
        order_clusters_as_data = FALSE
      )
    )
}
trait_profile_groups <- rbindlist(trait_profile_groups, idcol = "Region")
saveRDS(trait_profile_groups, file.path(path_cache, "trait_profile_groups.rds"))
