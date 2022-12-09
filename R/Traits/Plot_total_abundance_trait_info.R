# Violinplots on total abundance and available trait information
# Add % of complete information
abund[, richness_order := uniqueN(taxon), by = order]
abund[trait_info == "complete",
  completeness_order := round(uniqueN(taxon) / richness_order, digits = 2),
  by = order
]
completeness_dat <- unique(abund[
  !is.na(completeness_order),
  .(order, completeness_order, type)
])

abund_insects <- abund[type == "insect",
  .(family,
    order,
    total_abund = round(sum(abundance)),
    trait_info,
    completeness_order
  ),
  by = taxon
] %>%
  .[order(-total_abund), ] %>%
  unique(.) %>%
  ggplot(., ) +
  geom_violin(aes(x = as.factor(order), y = log(total_abund))) +
  geom_jitter(aes(
    x = as.factor(order), y = log(total_abund),
    col = trait_info
  ),
  width = 0.1,
  size = 3
  ) +
  geom_text(
    data = completeness_dat[type == "insect", ],
    aes(
      x = as.factor(order), y = 17,
      label = paste("complete:", completeness_order * 100, "%")
    )
  ) +
  coord_flip() +
  ylim(0, 20) +
  labs(x = "Order", y = "log(Total abundance)", color = "Trait information") +
  scale_color_manual(values = c("steelblue", "forestgreen", "#5d5d5d")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 12
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
    )
  )

abund_non_insects <- abund[type == "non_insect",
  .(family, order, total_abund = round(sum(abundance)), trait_info),
  by = taxon
] %>%
  .[order(-total_abund), ] %>%
  unique(.) %>%
  ggplot(., ) +
  geom_violin(aes(x = as.factor(order), y = log(total_abund))) +
  geom_jitter(aes(
    x = as.factor(order), y = log(total_abund),
    col = trait_info
  ),
  width = 0.1,
  size = 3
  ) +
  geom_text(
    data = completeness_dat[type == "non_insect", ],
    aes(
      x = as.factor(order), y = 17,
      label = paste("complete:", completeness_order * 100, "%")
    )
  ) +
  coord_flip() +
  ylim(0, 20) +
  labs(x = "Order", y = "log(Total abundance)", color = "Trait information") +
  scale_color_manual(values = c("steelblue", "forestgreen", "#5d5d5d")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      family = "Roboto Mono",
      size = 12
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
    )
  )