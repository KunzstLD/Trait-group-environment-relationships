if (params$region == "California") {
  # Add taxonomy to abundance data
  # Check for which taxon trait information is there, how abundant are these?
  ca_abund <- fread(
    "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Data/NoA/CA_NoAmbig_Species Matrix n82.csv"
  )

  # Conversion to long format
  # Some integer columns which will be coerced to double
  ca_abund_lf <- melt(
    ca_abund,
    id.vars = c(
      "TSTAID",
      "CollectionDate",
      "SMCOD",
      "BU_ID"
    ),
    variable.name = "taxon",
    value.name = "abundance"
  )
  ca_abund_lf[taxon %like% "Concha", taxon := "Conchapelopia"]
  abund <- copy(ca_abund_lf)
}

if (params$region == "Midwest") {
  midwest_abund <- read_xlsx("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Data/NoA/MSQA_Opt4_Comb_MV.xlsx")
  setDT(midwest_abund)

  # What about these sites, different SCODES but same STAID? 
  # midwest_abund[STAID %like% c("06808495|391114085205801|412911089540101"), ]
  
  midwest_abund_lf <- melt(
    midwest_abund,
    id.vars = c(
      "SCODE",
      "SUID",
      "STAID",
      "REACH",
      "DATE",
      "SampleID",
      "SMCOD"
    ),
    variable.name = "taxon",
    value.name = "abundance"
  )

  # For now, sum up duplicate taxon
  midwest_abund_lf[taxon %like% "Centroptilum", abundance := sum(abundance), by = "SCODE"]
  midwest_abund_lf <- midwest_abund_lf[taxon != "Centroptilum/Procloeon", ]

  midwest_abund_lf[taxon %like% "Leucrocuta", abundance := sum(abundance), by = "SCODE"]
  midwest_abund_lf <- midwest_abund_lf[taxon != "Leucrocuta/Nixe", ]

  midwest_abund_lf[taxon %like% "Corduliidae", abundance := sum(abundance), by = "SCODE"]
  midwest_abund_lf <- midwest_abund_lf[taxon != "Corduliidae/Libellulidae", ]

  midwest_abund_lf[taxon %like% "Phaenopsectra", abundance := sum(abundance), by = "SCODE"]
  midwest_abund_lf <- midwest_abund_lf[taxon != "Phaenopsectra/Tribelos", ]

  midwest_abund_lf[taxon %like% "Micropsectra", abundance := sum(abundance), by = "SCODE"]
  midwest_abund_lf <- midwest_abund_lf[taxon != "Micropsectra/Tanytarsus", ]

  midwest_abund_lf[taxon %like% "Cricotopus", abundance := sum(abundance), by = "SCODE"]
  midwest_abund_lf <- midwest_abund_lf[taxon != "Cricotopus/Orthocladius", ]

  midwest_abund_lf[, taxon := sub(" \\(.+\\)", "", taxon)]
  midwest_abund_lf[, taxon := sub("\\/.+", "", taxon)]
  midwest_abund_lf[, taxon := sub("(.+)(\\s)(.+)(\\s)(.+)", "\\1\\2\\3", taxon)]

  abund <- copy(midwest_abund_lf)
}

if (params$region == "Northeast") {
  northeast_abund <- fread(
    "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Data/NoA/NE_NoAmbig_Species_Matrix n92_Primer.csv"
  )
  # Conchapelopia (Helopelopia) & Cricotopus (Isocladius)` -> duplicates with zero Abundance, hence removed
  northeast_abund$`Conchapelopia (Helopelopia)` <- NULL
  northeast_abund$`Cricotopus (Isocladius)` <- NULL
  northeast_abund_lf <- melt(
    northeast_abund,
    id.vars = "TSTAID",
    variable.name = "taxon",
    value.name = "abundance"
  )
  northeast_abund_lf[taxon == "Conchapelopia (Helopelopia)", taxon := "Conchapelopia"]
  northeast_abund_lf[taxon == "Cricotopus (Isocladius)", taxon := "Cricotopus"]
  abund <- copy(northeast_abund_lf)
}

if (params$region == "PN") {
  pn_abund <- fread(
    "/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Data/NoA/PNSQA Species Matrix_Primer.csv"
  )
  pn_abund_lf <- melt(
    pn_abund,
    id.vars = "TSTAID",
    variable.name = "taxon",
    value.name = "abundance"
  )
  pn_abund_lf[, unique_id := paste0("id_", 1:nrow(pn_abund_lf))]
  
  # After name cleaning, a few duplicate taxa are "created"
  # names(pn_abund)[names(pn_abund) %like% c("Sperchon|Phaenopsectra|Micropsectra|Cricotopus|Eukiefferiella")]
  # Add up abundances of those (for now)
  pn_abund_lf[, taxon := sub("\\/.+", "", taxon)]
  # pn_abund_lf[, .N, by = taxon] %>% 
  #   .[N > 87, ]
  # pn_abund_lf[taxon == "Sperchon", ] %>% 
  #   .[order(TSTAID), ]
  pn_abund_lf[taxon %in% c("Sperchon",
                           "Phaenopsectra",
                           "Micropsectra",
                           "Cricotopus",
                           "Eukiefferiella"), abundance := sum(abundance), by = TSTAID]
  for(taxa in c("Sperchon",
                "Phaenopsectra",
                "Micropsectra",
                "Cricotopus",
                "Eukiefferiella")) {
    id <- pn_abund_lf[taxon == taxa,] %>%
      .[duplicated(TSTAID), unique_id]
    pn_abund_lf <- pn_abund_lf[!unique_id %in% id,]
  }
  pn_abund_lf[, taxon := sub(" group", "", taxon)]
  pn_abund_lf[, taxon := sub(" complex", "", taxon)]
  abund <- copy(pn_abund_lf)
}

if (params$region == "Southeast") {
  southeast_abund <- fread("/home/kunzst/Dokumente/Projects/Trait_DB/Trait-group-environment-relationships/Data/NoA/SESQA_Inverts_transp_TSTAID_n75_clean.csv")
  southeast_abund_lf <- melt(
    southeast_abund,
    id.vars = c(
      "TSTAID",
      "SHORT_NAME",
      "SiteType2USE",
      "UrbanCenter_tier"
    ),
    variable.name = "taxon",
    value.name = "abundance"
  )
  abund <- southeast_abund_lf
}
