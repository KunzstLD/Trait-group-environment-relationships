# Taxa that could not be merged on genus or family level by the Inidcate tool

# California
abund_ls$California[Taxa_new %in% c(
    "Stratiomyidae",
    "Brundiniella",
    "Psephenus",
    "Eubrianax",
    "Zaitzevia",
    "Ordobrevia",
    "Optioservus",
    "Cleptelmis",
    "Ampumixis",
    "Heteroplectron",
    "Amiocentrus",
    "Pteronarcyidae",
    "Rickera",
    "Kogotus",
    "Calineuria",
    "Peltoperlidae",
    "Brechmorhoga",
    "Octogomphus",
    "Diphetor",
    "Leptohyphidae",
    "Matriella"
), Taxa_new := fifelse(Taxa_new %like% ".*dae$", order, family)]


# Midwest
abund_ls$Midwest[Taxa_new %in% c(
    "Stelechomyia",
    "Psephenus",
    "Ancyronyx",
    "Potamyia",
    "Pteronarcyidae",
    "Pseudiron",
    "Stenonema",
    "Paracloeodes",
    "Fallceon",
    "Leptohyphidae",
    "Amercaenis",
    "Baetiscidae"
), Taxa_new := fifelse(Taxa_new %like% ".*dae$", order, family)]
abund_ls$Midwest[Taxa_new == "Pseudironidae", Taxa_new := order]

# Northeast
abund_ls$Northeast[Taxa_new %in% c(
    "Psephenus",
    "Optioservus",
    "Microcylloepus",
    "Ancyronyx",
    "Dineutus",
    "Mayatrichia",
    "Leucotrichia",
    "Psilotreta",
    "Pteronarcyidae",
    "Peltoperlidae",
    "Stylogomphus",
    "Hagenius",
    "Stenonema",
    "Diphetor",
    "Acerpenna",
    "Leptohyphidae",
    "Teloganopsis"
), Taxa_new := fifelse(Taxa_new %like% ".*dae$", order, family)]

# Northwest
abund_ls$Northwest[Taxa_new %in% c(
    "Pelecorhynchidae",
    "Brundiniella",
    "Zaitzevia",
    "Ordobrevia",
    "Narpus",
    "Heterlimnius",
    "Cleptelmis",
    "Ampumixis",
    "Leucotrichia",
    "Onocosmoecus",
    "Pteronarcys",
    "Pteronarcyidae",
    "Perlinodes",
    "Hesperoperla",
    "Calineuria",
    "Peltoperlidae",
    "Zapada",
    "Visoka",
    "Despaxia",
    "Diphetor",
    "Timpanoga",
    "Matriella",
    "Attenella"
), Taxa_new := fifelse(Taxa_new %like% ".*dae$", order, family)]
abund_ls$Northwest[Taxa_new == "Pteronarcyidae", Taxa_new := order]

# Southeast
abund_ls$Southeast[Taxa_new %in% c(
    "Dineutus",
    "Argia",
    "Stenonema",
    "Diphetor",
    "Sperchopsis",
    "Leucotrichia",
    "Pteronarcyidae",
    "Hagenius",
    "Acerpenna",
    "Dannella",
    "Stelechomyia",
    "Anchytarsus",
    "Peltoperlidae",
    "Stylogomphus",
    "Stenacron",
    "Optioservus",
    "Attenella",
    "Helichus",
    "Teloganopsis",
    "Microcylloepus",
    "Psephenus",
    "Xylotopus"
), Taxa_new := fifelse(Taxa_new %like% ".*dae$", order, family)]
abund_ls$Southeast[Taxa_new %in%
    c(
        "Ptilodactylidae",
        "Leptohyphidae",
        "Ancyronyx"
    ), Taxa_new := order]