# Prepare Trait DB Bob Zuelig
# traits are all binary coded and mutual exclusive
noa_trait_enhanced <- fread(file.path(path_in, "NoA", "TraitsEnhancedExpandedV2.csv"),
                            na.strings = c("NA", ""))
# str(noa_trait_enhanced)
# rowSums(noa_trait_enhanced[, .SD, .SDcols = patterns("locom.*")])
# View(noa_trait_enhanced)

# rename traits
setnames(
  noa_trait_enhanced,
  old = grep(
    "Volt.*|Resp.*|Size.*|Habi.*|Trop.*",
    names(noa_trait_enhanced),
    value = TRUE
  ),
  new = c(
    "volt_semi",
    "volt_uni",
    "volt_bi_multi",
    "resp_teg",
    "resp_gil",
    "resp_pls_spi",
    "size_small",
    "size_medium",
    "size_large",
    "locom_burrow",
    "locom_climb",
    "locom_sprawl",
    "locom_cling",
    "locom_swim",
    "locom_skate",
    "feed_gatherer",
    "feed_filter",
    "feed_herbivore",
    "feed_predator",
    "feed_shredder"
  )
)

# Harmonise traits ----
# volt ok
# resp ok
# size ok

## Locomotion ----
# crawl (crawler, walkers, climber, clinger)
# swim (swimmer, skater)
# burrow
# sessil is missing -> could be a problem
noa_trait_enhanced[, locom_swim := apply(.SD, 1, sum),
                   .SDcols = c("locom_swim",
                               "locom_skate")]
noa_trait_enhanced[, locom_crawl := apply(.SD, 1, sum),
                   .SDcols = c("locom_cling",
                               "locom_climb",
                               "locom_sprawl")]
noa_trait_enhanced[, c("locom_cling",
                       "locom_climb",
                       "locom_sprawl",
                       "locom_skate") := NULL]

## Feeding mode ----
# shredder (chewer, miners, xylophagus, herbivore piercers)
# gatherer (collector gatherer, detrivore)
# filter: collector filterer (active filterers, passive filterers, absorbers)
# herbivore: herbivore piercer & scraper (grazer)
# predator: predator
# parasite: parasite -> is missing

# Postprocessing ----
# Create subset and remove taxa with missing  values
trait_cols <- grep("volt.*|resp.*|size.*|locom.*|feed.*",
                   names(noa_trait_enhanced),
                   value = TRUE)
noa_trait_enhanced_subset <-
  noa_trait_enhanced[, .SD, .SDcols = c("Phylum",
                                        "Class",
                                        "Order",
                                        "Family",
                                        "Genus",
                                        "OTU",
                                        trait_cols)]

# Either complete or just NAs
# View(noa_trait_enhanced_subset[is.na(volt_uni), ])
# ind_na <- noa_trait_enhanced_subset[, apply(.SD, 1, function(y) {
#   sum(!is.na(y))
# }),
# .SDcols = trait_cols]
# all(!data.table::between(ind_na, 0, 17, incbounds = FALSE))
noa_trait_enhanced_subset <-
  na.omit(noa_trait_enhanced_subset[, .SD, .SDcols = c("OTU", trait_cols)])

# merge back taxonomy
noa_trait_enhanced_subset[noa_trait_enhanced,
                          `:=`(
                            Phylum = i.Phylum,
                            Class = i.Class,
                            Order = i.Order,
                            Family = i.Family,
                            Genus = i.Genus
                          ),
                          on = "OTU"]
# Change colorder so that traits of a grouping feature are next to each other
setcolorder(noa_trait_enhanced_subset,
            neworder = c(
              grep("feed.*", names(noa_trait_enhanced_subset), value = TRUE),
              grep("resp.*", names(noa_trait_enhanced_subset), value = TRUE),
              grep("locom.*", names(noa_trait_enhanced_subset), value = TRUE),
              grep("size.*", names(noa_trait_enhanced_subset), value = TRUE),
              grep("volt.*", names(noa_trait_enhanced_subset), value = TRUE)
            ))

# Handle ambiguous names in OTU
noa_trait_enhanced_subset[, OTU_new := sub("(.+)(\\/)(.+)", "\\1",OTU)]
noa_trait_enhanced_subset[, OTU_new := sub("(.+)(genus nr\\. )(.+)", "\\3", OTU_new)]
noa_trait_enhanced_subset[, OTU_new := sub(" genus.+", "", OTU_new)]
noa_trait_enhanced_subset[, `:=`(
  OTU_new = sub("(.+)( cf\\.)(.+)", "\\3", OTU_new),
  Genus = sub("(.+)( cf\\.)(.+)", "\\3", OTU_new)
)]

# By cleaning some of the OTUs, duplicates are created. Check if these duplicates contain
# the similar information
duplicates <- noa_trait_enhanced_subset[duplicated(OTU_new), OTU_new]
duplicates <- duplicates[duplicates != "Orthocladiinae"]
# Orthocladiinae multiple times, exclude and check manually

res_similarity <- list()
for(i in duplicates) {
  comp_tp <-
    noa_trait_enhanced_subset[OTU_new == i, .SD, .SDcols = c(trait_cols, "OTU_new")] %>%
    .[, OTU_new := paste0(OTU_new, "_", 1:.N)] %>%
    melt(., id.vars = "OTU_new") %>%
    dcast(., ... ~ OTU_new, value.var = "value") %>%
    .[, .(all(.SD[, 1] == .SD[, 2])), .SDcols = !"variable"]
  res_similarity[[i]] <- comp_tp
}
res_similarity <- unlist(res_similarity)

# Fix some of these manually
# Agnetina
res_similarity[res_similarity == 0]

# Manual fixes
noa_trait_enhanced_subset[OTU_new == "Agnetina" & grepl("\\/", OTU), OTU_new := "Paragnetina"]
# have similar trait profiles
noa_trait_enhanced_subset[OTU_new == "Orthocladiinae", OTU_new := OTU]
noa_trait_enhanced_subset[OTU_new == "Orthocladiinae", `:=`(Genus = OTU, OTU_new = OTU)]

# All duplicates that are now present have the same trait profiles
# Can thus be deleted
noa_trait_enhanced_subset <- noa_trait_enhanced_subset[!duplicated(OTU_new), ]
noa_trait_enhanced_subset[, OTU := NULL]

# save
saveRDS(
  noa_trait_enhanced_subset,
  file.path(path_in, "NoA", "trait_db_enhanced_preproc.rds")
)
