# Taxa richness plots
pl_taxa_richn_aq <- ggplot(taxonomy[order %in% aq_insects, ]) +
  geom_pointrange(aes(
    x = as.factor(order),
    y = taxon_richn,
    ymin = 0,
    ymax = taxon_richn
  )) +
  geom_text(aes(
    x = 8.5,
    y = 75,
    label = paste0("Aq insect taxa = ", aq_ins_total_N)
  )) +
  geom_text(aes(
    x = 8.5,
    y = 65,
    label = paste0("Taxa in total = ", total_N)
  )) +
  labs(
    x = "Order",
    y = "Taxa richness",
    title = paste("Taxonomy overview aquatic insects", params$region)
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 12,
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 12
    ),
    legend.title = element_text(
      family = "Roboto Mono",
      size = 16
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    legend.position = "none"
  )

# Remaining taxa
pl_taxa_richn_non_aq <- ggplot(taxonomy[!order %in% aq_insects, ]) +
  geom_pointrange(aes(
    x = as.factor(order),
    y = taxon_richn,
    ymin = 0,
    ymax = taxon_richn
  )) +
  labs(
    x = "Order",
    y = "Taxa richness",
    title = paste("Taxonomy overview non insects", params$region)
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 12,
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(
      family = "Roboto Mono",
      size = 12
    ),
    legend.title = element_text(
      family = "Roboto Mono",
      size = 16
    ),
    legend.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    strip.text = element_text(
      family = "Roboto Mono",
      size = 14
    ),
    legend.position = "none"
  )
