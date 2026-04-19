# ______________________________________________________________________________
# Create radar chart
# Here we compare Performance, Consistency and Specificity of individual traits and TPGs 
# to detect effects from pesticides
# TODO: 
# - 2-3 Panels in one figure
# - Labels: Avg. marginal effect?
# - Phosphate to PO4
# - ADD EPT & SPEAR?
# ______________________________________________________________________________

# ---- Performance ----

# Marginal effect/variance of max_log_tu & the other env. factors for traits and TPGs
rda_me <- readRDS(file.path(path_cache, "rda_me.rds"))
rda_me[type=="tpg", type := "TPGs"]
rda_me[type=="traits", type := "Traits"]

# For SPEAR
spear_type2ss <- readRDS(file.path(path_cache, "spear_type2ss.rds"))
spear_me <- spear_type2ss[, .(mean_marginal_effect = mean(prop_var)), by = c("term")]|> 
  dcast(.~term, value.var="mean_marginal_effect")|>
  setnames(
    c(
      "max_log_tu",
      "Riffle.FRC",
      "Temp.median",
      "orthoP_4wk.median"
    ),
    c("max logTU", "Riffles", "Temperature", "PO4")
  ) |> 
  _[, c("Residuals", ".") :=NULL] |> 
  _[, type := "SPEAR"]

# For EPT
ept_type2ss <- readRDS(file.path(path_cache, "ept_type2ss.rds"))
ept_me <- ept_type2ss[, .(mean_marginal_effect = mean(prop_var)), by = c("term")] |> 
  dcast(.~term, value.var="mean_marginal_effect")|>
  setnames(
    c(
      "max_log_tu",
      "Riffle.FRC",
      "Temp.median",
      "orthoP_4wk.median"
    ),
    c("max logTU", "Riffles", "Temperature", "PO4")
  ) |> 
  _[, c("Residuals", ".") := NULL] |>
  _[, type := "EPT"]


# Calc mean marginal effect across regions and prepare for radarchart
# Respect stat. significance when calculating mean marginal variance?
rda_me_wf <- rda_me[, .(type, region, env_factor, prop_variance)] |>
  _[, .(mean_marginal_effect = mean(prop_variance)), by = c("env_factor", "type")] |>
  dcast(type ~ env_factor, value.var = "mean_marginal_effect") |>
  setnames(
    c(
      "max_log_tu",
      "Riffle.FRC",
      "Temp.median",
      "orthoP_4wk.median"
    ),
    c("max logTU", "Riffles", "Temperature", "PO4")
  ) |> 
  _[, "Residual" := NULL]

# Add also explained variance of X in trait or TPG composition 
# Remove Adj. R2 + Residuen from radarchart 
rda_variance <- fread(file.path(path_out, "rda_variance.csv"))
rda_variance <- rda_variance[, lapply(.SD, mean), .SDcols = c("adj_constr_var_traits", "adj_constr_var_tpgs")]
rda_variance <- melt(rda_variance, variable.name = "type", value.name = "Adj. R2")
rda_variance[, type := sub("adj_constr_var_", "", type)]
rda_variance[type=="tpgs", type := sub("tpgs", "tpg", type)]

# Create radar chart
summary_me <- rbind(rda_me_wf, spear_me, ept_me)
radarchart_performance <-  summary_me[, .SD, .SDcols = c("type",
                                                       "max logTU",
                                                       "Riffles",
                                                       "Temperature",
                                                       "PO4")] 

create_ggplot_radarchart <- function(data,
                                     color = c("#E15759", "#4A90E2","#00AFBB", "#E7B800"),
                                     title = NULL,
                                     axis.labels = colnames(data)[-1],
                                     caxislabels = NULL, 
                                     values.radar = NULL,
                                     legend.position = "bottom",
                                     plot.legend=TRUE) {
  # Assumes 'data' has a 'group' column followed by variables to be plotted
  ggradar(
    data,
    axis.labels = axis.labels,
    grid.min = min(caxislabels),
    grid.mid = median(caxislabels),
    grid.max = max(caxislabels),
    values.radar = values.radar,
    group.line.width = 1.5,
    group.point.size = 3.5,
    group.colours = color,
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    gridline.min.colour = "grey",
    gridline.max.colour = "grey",
    axis.label.size = 7.5,
    legend.position = legend.position,
    plot.title = title,
    legend.text.size=20,
    plot.legend = plot.legend
  )
}
plot_performance <- create_ggplot_radarchart(radarchart_performance, 
                         title = "(a) - Average marginal effect",
                         caxislabels = c(0, 5, 10, 15, 20), 
                         values.radar = c("0", "10", "20"), 
                         plot.legend = FALSE) +
  coord_equal(clip = "off") + 
  theme(plot.margin = margin(10, 10, 20, 10))

# ggsave(file.path(path_paper, "Graphs", "radarchart_perf_specif.png"),
#     width = 27,
#     height = 17,
#     units = "cm",
#     dpi=600
# )

# ______________________________________________________________________________
# Consistency & Specificity ----
# rda_cosine <- readRDS(file.path(path_cache, "rda_cosine.rds"))
# ______________________________________________________________________________

## Consistency of stat. significance ----

# Is the marginal effect of max lgTU stat. & the env. factors significant across the different regions?
# Traits & TPGs
rda_consistency <- rda_me[p.value <= 0.05, .N, by = c("type", "env_factor")] |>
  dcast(type ~ env_factor, value.var="N") |> 
  _[, .(type, max_log_tu, Riffle.FRC, Temp.median)] |> 
  setnames(
    c("max_log_tu", "Riffle.FRC", "Temp.median"),
    c(
      "max logTU",
      "Riffles",
      "Temperature"
    )
  )

# Calc. consistency of sig. for EPT & SPEAR 
ept_consistency <- ept_type2ss[p.value <= 0.05, .N, by="term"] |>
  dcast(.~ term, value.var = "N") |> 
  _[, "." := NULL] |> 
  _[, type := "EPT"] |> 
  setnames(
    c("max_log_tu", "Riffle.FRC", "Temp.median"),
    c(
      "max logTU",
      "Riffles",
      "Temperature"
    )
  )
spear_consistency <- spear_type2ss[p.value <= 0.05, .N, by="term"] |>
  dcast(.~ term, value.var = "N") |> 
  _[, "." := NULL] |> 
  _[, type := "SPEAR"] |> 
  setnames(
    c("max_log_tu", "orthoP_4wk.median"),
    c(
      "max logTU",
      "PO4"
    )
  )
summary_consistency <- rbind(rda_consistency, ept_consistency, spear_consistency, fill=TRUE)

# Replace all NANs with zeros
summary_consistency[, c("max logTU", "Riffles", "Temperature", "PO4") := lapply(.SD, function(x)
  fifelse(is.na(x), 0, x)), .SDcols = c("max logTU", "Riffles", "Temperature", "PO4")]
                             
# Plot
plot_consistency <- create_ggplot_radarchart(summary_consistency,
                         title = "(b) - Statistical significance",
                         caxislabels = c(0, 3, 5),
                         values.radar = c("0", "3", "5"),
                         legend.position = "bottom") +
  coord_equal(clip = "off") +
  theme(plot.margin = margin(10, 10, 20, 10))

# ggsave(
#   file.path(
#     path_paper,
#     "Graphs",
#     "radarchart_cons_specif_significance.png"
#   ),
#   width = 27,
#   height = 17,
#   units = "cm",
#   dpi=600
# )

## Directional consistency ----
# Is the correlation between traits and TPGs with the env. factors mainly positive, 
# mainly negative, or rather no clear pattern? 
# Based on the cosine similarity of traits and TPGs with env. factors 
rda_cosine <- readRDS(file.path(path_cache, "rda_cosine.rds"))
setnames(rda_cosine, "rn", "response")
rda_cosine[type=="tpg", type := "TPGs"]
rda_cosine[type=="traits", type := "Traits"]
rda_cosine_agg <- rda_cosine[, lapply(.SD, mean), .SDcols = c("max_log_tu",
                                                              "Riffle.FRC",
                                                              "Temp.median",
                                                              "orthoP_4wk.median"), by = "response"]
rda_cosine_agg[, type := fifelse(response %like% "TPG.*", "tpg", "trait")]
rda_cosine_agg[type=='tpg' & max_log_tu <0, .N]
rda_cosine_agg[type=='tpg' & max_log_tu >0, .N]
rda_cosine_agg[type=='trait' & max_log_tu <0, .N]
rda_cosine_agg[type=='trait' & max_log_tu >0, .N]

rda_cosine_agg[, c("max_log_tu",
                   "Riffle.FRC",
                   "Temp.median",
                   "orthoP_4wk.median") := lapply(.SD, function(x)
                     round(x, digits = 2)), .SDcols = c("max_log_tu",
                                                        "Riffle.FRC",
                                                        "Temp.median",
                                                        "orthoP_4wk.median")] |>
  _[, .(response,
        max_log_tu,
        Riffle.FRC,
        Temp.median,
        orthoP_4wk.median)] |>
  fwrite(file.path(path_paper, 'Tables', 'rda_cosine_avg_region.csv'))

# Agg. positive responding
rda_cosine_agg_pos <- melt(rda_cosine_agg, id.vars = c("response", "type")) |>
  _[value >= 0, mean(value), by = c("type", "variable")] |>
  dcast(type ~ variable, value.var = "V1")
setnames(
  rda_cosine_agg_pos,
  c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median"
  ),
  c(
    "max logTU (Pos.)",
    "Riffles (Pos.)",
    "Temperature (Pos.)",
    "PO4 (Pos.)"
  )
)

# Agg. negative responding
rda_cosine_agg_neg <- melt(rda_cosine_agg, id.vars = c("response", "type")) |>
  _[value < 0, mean(value), by = c("type", "variable")] |>
  dcast(type ~ variable, value.var = "V1")
setnames(
  rda_cosine_agg_neg,
  c(
    "max_log_tu",
    "Riffle.FRC",
    "Temp.median",
    "orthoP_4wk.median"
  ),
  c(
    "max logTU (Neg.)",
    "Riffles (Neg.)",
    "Temperature (Neg.)",
    "PO4 (Neg.)"
  )
)

# Prepare radarchart
radarchart_consistency_similarity <- merge(rda_cosine_agg_pos, rda_cosine_agg_neg, 
                                           by = "type")

# Plot
radarchart_consistency_similarity[type=="tpg", type := "TPGs"]
radarchart_consistency_similarity[type=="traits", type := "Traits"]

plot_cosine <- radarchart_consistency_similarity[, .(
  type,
  `max logTU (Pos.)`,
  `Riffles (Pos.)`,
  `Temperature (Pos.)`,
  `max logTU (Neg.)` = abs(`max logTU (Neg.)`),
  `Riffles (Neg.)` = abs(`Riffles (Neg.)`),
  `Temperature (Neg.)` = abs(`Temperature (Neg.)`)
)] |>
  create_ggplot_radarchart(
    title = "(c) - Cosine similarity",
    caxislabels = c(0, 0.25, 0.50, 0.75, 1),
    values.radar = c("0", "50", "100"),
    color = c("#00AFBB", "#E7B800"),
    plot.legend = FALSE
  ) +
  coord_equal(clip = "off") + 
  theme(plot.margin = margin(10, 10, 20, 10))
# ggsave(
#   file.path(
#     path_paper,
#     "Graphs",
#     "radarchart_cons_specif_similarity.png"
#   ),
#   width = 27,
#   height = 17,
#   units = "cm",
#   dpi = 600
# )

# --- Multipanel plot ----
# combined_plot <- (plot_performance + plot_consistency) / plot_cosine 
combined_plot <- wrap_plots(plot_performance, plot_consistency, plot_cosine, nrow=2)

ggsave(
  plot = combined_plot,
  file.path(
    path_paper,
    "Graphs",
    "multipanel_comp_figure.png"
  ),
  width = 45,
  height = 31,
  units = "cm",
  dpi = 600
)


# # Exclude Phosphate, as it did  not show stat. significant results in the marginal effects 
# radarchart_consistency_similarity <- radarchart_consistency_similarity[, !c("Pos. cos. sim. Phosphate", "Neg. cos. sim. Phosphate")]
# max_min_consistency_similarity <- data.table(
#   "Pos. cos. sim. max logTU" = c(1, 0),
#   "Pos. cos. sim. Riffles" = c(1,0),
#   "Pos. cos. sim. Temperature" = c(1,0),
#   #  "Pos. cos. sim. Phosphate" = c(1,0), 
#   "Neg. cos. sim. max logTU" = c(-1, 0),
#   "Neg. cos. sim. Riffles" = c(-1, 0),
#   "Neg. cos. sim. Temperature" = c(-1, 0)
#   # "Neg. cos. sim. Phosphate" = c(0, -1)
# )
#   
# png(file.path(path_paper, "Graphs", "radarchart_cons_specif_similarity.png"), width = 800, height = 600)
# radarchart_consistency_similarity[, -"type"] |>
#   (\(x) rbind(max_min_consistency_similarity, x, fill = TRUE))() |>
#   create_beautiful_radarchart(
#     caxislabels = c(0, 25, 50, 75, 100),
#     color = c("#00AFBB", "#E7B800"),
#     title = "TPGs - Trait Consistency & Specificity"
#   )
# # Add an horizontal legend
# legend(
#   x = -2,
#   y = 1.2,
#   legend = c("TPGs", "Traits"),
#   horiz = FALSE,
#   bty = "n",
#   pch = 20 ,
#   col = c("#00AFBB", "#E7B800"),
#   text.col = "black",
#   cex = 1.3,
#   pt.cex = 1.2
# )
# dev.off()


#_______________________________________________________________________________

## Direction Marginal effect in ordination space
# rda_scores <- readRDS(file.path(path_cache, "rda_scores.rds"))
# rda_scores[, quadrant := fcase(
#   angle_env_factor >= 0 & angle_env_factor < 90,
#   "q1",
#   angle_env_factor >= 90 &
#     angle_env_factor < 180,
#   "q2",
#   angle_env_factor >= -180 &
#     angle_env_factor < -90,
#   "q3",
#   angle_env_factor >= -90 &
#     angle_env_factor < 0,
#   "q4"
# )]
# In which quadrant is the env. factor or max logTU pointing?
# rda_scores[, .N, by = c("direction", "type", "env_factor", "quadrant")] |> 
#   _[env_factor=="max_log_tu", ]
# 
# # Identify the most common arrow direction in the quadrant and then compute the fraction. Max is 5
# rda_scores_all <- rda_scores[, .N, by = c("direction", "type", "env_factor", "quadrant")] |>
#   _[direction == "neg", ] |>
#   _[, .(fraction_quadrants = max(N) / 5), by = c("type", "env_factor")] |>
#   dcast(type ~ env_factor, value.var = "fraction_quadrants")
# rda_scores_all <- rda_scores_all[, .(type, Riffle.FRC, Temp.median, max_log_tu)]
# 
# setnames(
#   rda_scores_all,
#   c(
#     "max_log_tu",
#     "Riffle.FRC",
#     "Temp.median"
#   ),
#   c(
#     "Fr. similar quadrant max logTU",
#     "Fr. similar quadrant Riffles",
#     "Fr. similar quadrant Temperature"
#   )
# )

