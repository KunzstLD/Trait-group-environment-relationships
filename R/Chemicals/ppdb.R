# __________________________________________________________________________________________________
# PPDB database
# https://sitem.herts.ac.uk/aeru/ppdb/en/search.htm
# Endpoint EC50 (non-lethal according to the ppdb folks)
# __________________________________________________________________________________________________
ppdb <- readRDS(file.path(path_in, "Chemicals", "ppdb_full_20180427.rds"))

# DB is ordered according to chemicals
# extract ecotox information
# names(ppdb[[1]])
ppdb_ecotox <-
  rbindlist(lapply(ppdb, `[[`, "etox"), idcol = "pesticide")

# CAS RN 
ppdb_cas <- lapply(ppdb, function(y) y$general[y$general$variable == "CAS RN", ])
ppdb_cas <- rbindlist(ppdb_cas, idcol = "pesticide")
setnames(ppdb_cas, "value", "cas")
ppdb_cas[, variable := NULL]

# PPDB data for sediment dwelling organisms and aquatic invertebrates
ppdb_ecotox <- ppdb_ecotox[grepl("Aquatic invertebrates", Property), ]

# Merge CAS RN
ppdb_ecotox[ppdb_cas, cas := i.cas, on = "pesticide"]

# Now cross check ppdb with
ecotox_missing <- readRDS(file.path(data_cache, "ecotox_missing.rds"))
ppdb_ecotox_subet <- ppdb_ecotox[cas %in% na.omit(ecotox_missing$CASRN), ]

# 48 to 96 hours, EC/LC50
ppdb_ecotox_subet <- ppdb_ecotox_subet[grepl("48.+EC50", Property), ]

# Value column needs to be numeric
# Values that have >?
# Units are all in mg/l, convert to ug/l
ppdb_ecotox_subet[, Value := as.numeric(sub("> ", "", Value))] 
ppdb_ecotox_subet <- ppdb_ecotox_subet[, .(pesticide,
                                           cas,
                                           conc_ug_l = Value * 1000,
                                           Property,
                                           `Source/Quality`)]
fwrite(ppdb_ecotox_subet,
       file.path(path_in, "Chemicals", "PPDB_queried.csv"))

# Some degradates of Fentin are listed, but Fentin itself not
# ecotox_missing[Chemname %like% "Fentin", ]
# meta_rsqa_cmax[Chemname %like% "Fentin", ]
# ppdb_ecotox[pesticide %like% "Fentin" & grepl("Aquatic.*", Property), ]
# stx_query(
#     cas = c("962-58-3")
#     , # CAS RN from chem_lup
#     concentration_type = "active ingredient",
#     endpoint = "XX50",
#     ecotox_grp = "invertebrate",
#     habitat = "freshwater",
#     duration = c(48, 96),
#     taxa = NULL,
#     effect = "Mortality"
#   )
