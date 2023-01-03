# __________________________________________________________________________________________________
# Create CWM per trait ----
# TODO: organic sensitivty?
# __________________________________________________________________________________________________

# Trait data (with TPGs)
trait_matrix <- readRDS(file.path(path_cache, "trait_matrix_tpgs.rds"))
trait_matrix_ls <- split(trait_matrix, f = trait_matrix$Region)
trait_matrix_lf <- melt(
  trait_matrix,
  id.vars = c(
    "Region",
    "taxon",
    "group",
    "n_taxa_group"
  ),
  variable.name = "trait"
)
trait_names <-
  grep("feed.*|resp.*|locom.*|size.*|volt.*|sensitivity.*",
       names(trait_matrix),
       value = TRUE)

## Abundance data ----
# TODO Check Midwest, site IDs fixed but somewhere are duplicates or some strange values
# because dcast does fall back to length
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))

# Fix STAID of Midwest data for later match with ecotox data 
# abund[!is.na(STAID), STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
# abund[!is.na(STAID), site := paste0("T", STAID)]
abund_ls <- split(abund, f = abund$Region)

# Few taxa have to be excluded because no trait information is available
abund_ls <- Map(function(x,y) x[taxon %in% y$taxon, ], abund_ls, trait_matrix_ls)
abund_wf <- lapply(abund_ls, function(y)
  y[, .(Region,
        site,
        taxon,
        species,
        genus,
        family,
        order,
        taxonomic_level,
        abundance)] %>%
    dcast(., site ~ taxon, value.var = "abundance"))

# library(ggbeeswarm)
# ggplot(abund, aes(x = Region, 
#                   y = log(abundance+1))) +
#   geom_quasirandom()



## Calc cwm for each trait ----
data_cwm <- Map(
  function(x, y)
    calc_cwm(
      abund = x,
      trait = y,
      trait_names = trait_names
    ),
  abund_wf,
  trait_matrix_ls
)
# Test Code:
# taxa_names <- names(abund_wf$California)[names(abund_wf$California) != "site"]
# test_herbivore <-
#   trait_matrix[Region == "California" &
#                  taxon %in% names(abund_wf$California), .(taxon, feed_herbivore)]
# test_herbivore <- test_herbivore[match(taxa_names, taxon) ]
# abund_wf$California[, apply(.SD,
#                             1,
#                             function(x)
#                               weighted.mean(test_herbivore$feed_herbivore, w = x)),
#                     .SDcols = !c("site"),
#                     by = "site"]
data_cwm <- lapply(data_cwm, function(y) rbindlist(y, idcol = "trait"))
data_cwm <- lapply(data_cwm, function(y) setnames(y, "V1", "cwm_val"))
saveRDS(data_cwm, file.path(path_cache, "data_cwm.rds"))
