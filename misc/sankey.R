library(ggsankey)

df <- mtcars %>%
  make_long(cyl, vs, am, gear, carb)

ggplot(
  df,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node)
  )
) +
  geom_sankey() +
  scale_fill_discrete(drop = FALSE)


ggplot(
  df,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = node
  )
) +
  geom_sankey(
    flow.alpha = .6,
    node.color = "gray30"
  ) +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = .5)
  ) +
  ggtitle("Car features")


ggplot(
  df,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = node
  )
) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = .5)
  ) +
  ggtitle("Car features")


df <- gapminder %>%
  dplyr::group_by(continent, year) %>%
  summarise(
    gdp = (sum(pop * gdpPercap) / 1e9) %>% round(0),
    .groups = "keep"
  ) %>%
  ungroup()

ggplot(
  df,
  aes(
    x = year,
    node = continent,
    fill = continent,
    value = gdp
  )
) +
  geom_sankey_bump(
    space = 0,
    type = "alluvial",
    color = "transparent",
    smooth = 6
  ) +
  scale_fill_viridis_d(option = "A", alpha = .8) +
  theme_sankey_bump(base_size = 16) +
  labs(
    x = NULL,
    y = "GDP ($ bn)",
    fill = NULL,
    color = NULL
  ) +
  theme(legend.position = "bottom") +
  labs(title = "GDP development per continent")
