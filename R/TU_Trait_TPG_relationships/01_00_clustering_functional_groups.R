# __________________________________________________________________________________________________
# Functional groups ----
# - Use cluster analysis to delineate groups with similar trait profiles
# - Then use the abundance per TPG
# TODO: How to incorporate organic sensitivity?
# __________________________________________________________________________________________________

# Trait data ----
# TODO: Method description, comparison overview all detected taxa & av. trait information
trait_matrix <-
  readRDS(file.path(path_cache, "trait_matrix_final.rds"))
# str(trait_matrix)
trait_names <-
  grep("feed.*|resp.*|locom.*|size.*|volt.*",
       names(trait_matrix),
       value = TRUE)

# Total abundance data
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))

# Distance matrices & Clustering ----

## Testing optimal number of groups ----
# Loop for deciding which measure for the optimal number of groups to use
results_total <- vector("list", 10)
results_nog <- list()
for (j in 1:10) {
  for (i in unique(trait_matrix$Region)) {
    dat <- trait_matrix[Region == i, .SD, .SDcols = c("taxon",
                                                      trait_names,
                                                      "sensitivity_organic")]
    
    # convert to data.frame -> data table does not support row.names
    setDF(dat)
    
    # add row.names
    row.names(dat) <- dat$taxon
    dat$taxon <- NULL
    
    # Convert to ktab object
    # TODO this section needs improvement but works atm
    vec <- sub("\\_.*", "\\1", trait_names)
    blocks <- rle(vec)$lengths
    prep_fuzzy_data <- prep.fuzzy(dat[, trait_names], blocks)
    ktab <- ktab.list.df(list(prep_fuzzy_data,
                              dat[, "sensitivity_organic", drop = FALSE]))
    dist_mat <-
      dist.ktab(ktab, type = c("F", "Q")) # check function description
    
    # Estimate optimal number of groups
    gap <- clusGap(
      x = as.matrix(dist_mat),
      FUN = mycluster_hc,
      K.max = 30,
      B = 100
    )
    optimal_nog_gap <- maxSE(gap$Tab[, "gap"],
                             gap$Tab[, "SE.sim"],
                             method = "Tibs2001SEmax")
    optimal_nog_silh <- NbClust(
      diss = dist_mat,
      distance = NULL,
      min.nc = 2,
      max.nc = 30,
      method = "ward.D2",
      index = "silhouette"
    )
    
    results_nog[[i]] <- list("gap_statistic" = optimal_nog_gap,
                             "silhouette" = optimal_nog_silh$Best.nc)
  }
  results_total[[j]] <- results_nog
}
# lapply(results_cl, function(y) y$optimal_nog)

# Gap statistic
lapply(results_total, function(y)
  lapply(y, function(z)
    z$gap_statistic)) %>%
  rbindlist() %>%
  saveRDS(., file.path(path_cache, "Gap_statistic_results.rds"))

# Silhouette width gives always almost the maximum value
# seems not to be a good option!
lapply(results_total, function(y)
  lapply(y, function(z)
    z$silhouette)) %>%
  rbindlist() %>%
  saveRDS(., file.path(path_cache, "Silhouette_results.rds"))

# Gap statistic (after 10 iterations):
# - California: 20, three times also 12
# - Midwest: 15
# - Northeast: 18 (5 times), 13 (4 times), 16 (once)
# - PN: 16 (6 times), 10 (4 times)
# - Southeast 10


## Actual clustering ----
# Use convergence across 
results_cl <- vector("list", 5)
names(results_cl) <- unique(trait_matrix$Region)

for (i in unique(trait_matrix$Region)) {
  dat <- trait_matrix[Region == i, .SD, .SDcols = c("taxon",
                                                    trait_names,
                                                    "sensitivity_organic")]
  
  # convert to data.frame -> data table does not support row.names
  setDF(dat)
  
  # add row.names
  row.names(dat) <- dat$taxon
  dat$taxon <- NULL
  
  # Convert to ktab object
  vec <- sub("\\_.*", "\\1", trait_names)
  blocks <- rle(vec)$lengths
  prep_fuzzy_data <- prep.fuzzy(dat[, trait_names], blocks)
  ktab <- ktab.list.df(list(prep_fuzzy_data,
                            dat[, "sensitivity_organic", drop = FALSE]))
  
  # Organic sensitivity is on a different scale then the other traits,
  # conversion is handled by ade4 package (Gower function should correct with the
  # max range to obtain a value between 0 and 1)
  dist_mat <-
    dist.ktab(ktab, type = c("F", "Q")) # check function description
  
  # Creation of distance matrix & hierarchical cluster analysis
  hc_taxa <- hclust(dist_mat, method = "ward.D2")
  # Get labels of dendrogram
  dend_label <- hc_taxa %>%
    as.dendrogram() %>%
    labels()
  
  # save to list
  results_cl[[i]] <- list(
    "distance_matrix" = dist_mat,
    "hc_wardD2" = hc_taxa,
    "labels_dendrogram" = dend_label
  )
}
saveRDS(results_cl, 
        file.path(path_cache, "Results_clustering.rds"))
