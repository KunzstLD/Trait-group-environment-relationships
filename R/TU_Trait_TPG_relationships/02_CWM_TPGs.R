# ___________________________
# Create CWM per trait ----
# ___________________________

# Trait data (with TPGs)
trait_matrix <- readRDS(file.path(path_cache, "trait_matrix_final.rds"))
trait_names <-
  grep("feed.*|resp.*|locom.*|size.*|volt.*|sensitivity.*",
       names(trait_matrix),
       value = TRUE)
trait_matrix <- trait_matrix[, .SD, .SDcols = c(
  "Region",
  "taxon",
  "species",
  "genus",
  "family",
  "order",
  trait_names
)]

trait_matrix_ls <- split(trait_matrix, f = trait_matrix$Region)
trait_matrix_lf <- melt(
  trait_matrix,
  id.vars = c(
    "Region",
    "taxon",
    "species",
    "genus",
    "family",
    "order"
  ),
  variable.name = "trait"
)

## Abundance data ----
# Match trait data
abund <- readRDS(file.path(path_cache, "total_abund_CEOPT.rds"))

# Few taxa cannot be considered (see also trait matching section)
# unique(abund[!family %in% trait_matrix$family, family])

# Some STAID of Midwest are duplicates, though the sites are different
# Taken care of in merge with ecotox data 
# abund[!is.na(STAID), .(.N, site), by = "STAID"] %>% 
#    .[N > 222, ] %>% 
#    unique(.)
abund_ls <- split(abund, f = abund$Region)
abund_ls <- Map(
  function(x, y) x[taxon %in% y$taxon, ],
  abund_ls,
  trait_matrix_ls
)

abund_wf <- lapply(abund_ls, function(y)
  y[, .(Region,
        site,
        STAID,
        taxon,
        species,
        genus,
        family,
        order,
        taxonomic_level,
        abundance)] %>%
    dcast(., site ~ taxon, value.var = "abundance"))


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
data_cwm$Midwest[abund_ls$Midwest, STAID := i.STAID, on = "site"]
data_cwm$Midwest[, STAID := sub("\"([0-9]{1,})\"", "\\1", STAID)]
data_cwm$Midwest[, STAID := paste0("T", STAID)]
data_cwm <- lapply(data_cwm, function(y) setnames(y, "V1", "cwm_val"))
saveRDS(data_cwm, file.path(path_cache, "data_cwm.rds"))
