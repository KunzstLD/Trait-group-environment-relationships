# ______________________________________________________________________
# Data processing
# TODO SumTU or some other Toxicity measure?
# ______________________________________________________________________

# Quality assessment California ---- 

# Abundance data: 270 taxa
ca_qa <- fread(file.path(path_in, "NoA", "CA_NoAmbig_Species Matrix n82.csv"))
# names(ca_qa)
# str(ca_qa)

# ____________________________________________________________________
## Establishing taxonomy ----
# ____________________________________________________________________

# Load trait data 
noa_traits <- readRDS(file.path(path_in, "NoA", "Traits_US_LauraT_pp_harmonized.rds"))
noa_traits[, taxon := coalesce(species, genus, family, order)]
noa_traits <- noa_traits[, .SD, 
                         .SDcols = patterns("species|genus|family|order|taxon|feed|resp|size|locom|volt")]

# Conchapelopia (Helopelopia) - Chirononmidae genus, which one?
# Assigned Conchapelopia 
taxa_names <-
  names(ca_qa)[!names(ca_qa) %in% c("TSTAID", "CollectionDate", "SMCOD", "BU_ID")]
taxa_names[taxa_names == "Conchapelopia (Helopelopia)"] <- "Conchapelopia"
ca_taxa <- as.data.table(taxa_names)
setnames(ca_taxa, "taxa_names", "taxon")

ca_taxa[noa_traits,
        `:=`(
          species = i.species,
          genus = i.genus,
          family = i.family,
          order = i.order
        ),
        on = "taxon"]

# Using taxize
ID_gbif <-
  get_gbifid(ca_taxa[is.na(species) &
                       is.na(genus) & is.na(family) & is.na(order),
                     taxon])
ID_merge <-
  data.table("taxon" = ca_taxa[is.na(species) &
                                 is.na(genus) &
                                 is.na(family) &
                                 is.na(order), taxon],
             "query" = ID_gbif)
tax_classf <- classification(
  id = ID_gbif,
  db = "gbif"
)
tax_classf <- rbind(tax_classf)
setDT(tax_classf)
# merge back original names of the queried taxa
tax_classf[ID_merge,
           merge_taxa := i.taxon,
           on = "query"]
tax_classf <- dcast(tax_classf, query + merge_taxa ~ rank, value.var = "name")

# merge to CA taxonomy
ca_taxa[tax_classf,
        `:=`(
          species = i.species,
          genus = i.genus,
          family = i.family,
          order = i.order,
          class = i.class,
          phylum = i.phylum
        ),
        on = c(taxon = "merge_taxa")]

# Search for remaining taxa
ca_taxa[taxon == "Prosobranchia", class := taxon]
ca_taxa[taxon == "Forcipomyiinae", `:=`(
  subfamily = taxon,
  family = "Ceratopogonidae",
  order = "Diptera"
)]
ca_taxa[taxon == "Hydroporinae", `:=`(
  subfamily = taxon,
  family = "Dytiscidae",
  order = "Coleoptera"
)]
ca_taxa[taxon == "Macropelopiini", `:=`(
  subfamily = taxon,
  family = "Chironomidae",
  order = "Diptera"
)]
ca_taxa[taxon == "Hemerodromiinae", `:=`(
  subfamily = taxon,
  family = "Empididae",
  order = "Diptera"
)]
ca_taxa[taxon == "Pentaneurini", `:=`(
  tribe = taxon,
  subfamily = "Tanypodinae",
  family = "Chironomidae",
  order = "Diptera"
)]
ca_taxa[taxon == "Macropelopiini", `:=`(
  tribe = taxon,
  subfamily = "Tanypodinae",
  family = "Chironomidae",
  order = "Diptera"
)]
ca_taxa[taxon == "Acari", class := taxon]
ca_taxa[taxon == "Ceratopogoninae", `:=`(
  subfamily = taxon,
  family = "Ceratopogonidae",
  order = "Diptera"
)]
ca_taxa[taxon %like% ".*inae", subfamily := taxon]

# Postprocessing taxonomy
ca_taxa[grepl("sp\\.", species), species := NA_character_]
ca_taxa[!grepl("\\s", species) &
          !is.na(species), species := paste(genus, species)]

ca_taxa[, taxonomic_level := fcase(
  !is.na(species),
  "species",
  is.na(species) & !is.na(genus),
  "genus",
  is.na(species) & is.na(genus) & !is.na(tribe),
  "tribe",
  is.na(species) & is.na(genus) & is.na(tribe) & !is.na(subfamily),
  "subfamily",
  is.na(species) & is.na(genus) & is.na(tribe) & is.na(subfamily) & !is.na(family),
  "family",
  is.na(species) &
    is.na(genus) & is.na(tribe) & 
    is.na(subfamily) & is.na(family) & !is.na(order),
  "order",
  is.na(species) &
    is.na(genus) & is.na(tribe) & is.na(subfamily) & is.na(family) & is.na(order) & !is.na(class),
  "class",
  is.na(species) &
    is.na(genus) & is.na(tribe) & is.na(subfamily) &
    is.na(family) & is.na(order) & is.na(class) & !is.na(phylum),
  "phylum"
)]
setcolorder(
  ca_taxa,
  c(
    "taxon",
    "species",
    "genus",
    "subfamily",
    "family",
    "order",
    "class",
    "phylum",
    "taxonomic_level"
  )
)

# Save taxonomy
saveRDS(ca_taxa, file = file.path(path_cache, "ca_taxonomy.rds"))
