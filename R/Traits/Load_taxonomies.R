# Load taxonomies
if (params$region == "California") {
  ca_taxonomy <-
    readRDS(
      "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Cache/ca_taxonomy.rds"
    )
  taxonomy <- copy(ca_taxonomy)
}

if (params$region == "Midwest") {
  midwest_taxonomy <- readRDS("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Cache/midwest_taxonomy.rds")
  taxonomy <- copy(midwest_taxonomy)
}

if (params$region == "Northeast") {
  northeast_taxonomy <- readRDS("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Cache/northeast_taxonomy.rds")
  taxonomy <- copy(northeast_taxonomy)
}

if (params$region == "PN") {
  pn_taxonomy <- readRDS("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Cache/pn_taxonomy.rds")
  taxonomy <- copy(pn_taxonomy)
}

if (params$region == "Southeast") {
  southeast_taxonomy <- readRDS("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Cache/southeast_taxonomy.rds")
  taxonomy <- copy(southeast_taxonomy)
}