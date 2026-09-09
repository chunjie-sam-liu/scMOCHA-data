# Prismatic Usage

Use `prismatic` to inspect, transform, and validate color vectors created with `paletteer` or stored in project color files.

## Preview

```r
plot(prismatic::color(celltype_colors))
```

This is useful before committing color assignments to a final figure.

## Derived Variants

Keep canonical identity colors unchanged. Create explicit derived objects:

```r
celltype_colors_light <- prismatic::clr_lighten(celltype_colors, 0.25)
celltype_colors_dark <- prismatic::clr_darken(celltype_colors, 0.20)
celltype_colors_alpha <- prismatic::clr_alpha(celltype_colors, 0.75)
```

Use derived colors for labels, borders, backgrounds, low-emphasis annotations, or secondary panels. Do not replace the base named vector unless the user asks to change the canonical mapping.

## Color Vision Deficiency Checks

Preview common simulations:

```r
plot(prismatic::clr_deutan(celltype_colors))
plot(prismatic::clr_protan(celltype_colors))
plot(prismatic::clr_tritan(celltype_colors))
```

If important categories become visually close under simulation, recommend changing specific category colors or choosing a different palette.

## Contrast Checks

Use contrast helpers for text or label colors:

```r
prismatic::contrast_ratio("#000000", "#FFFFFF")
prismatic::best_contrast("#4E79A7", c("#000000", "#FFFFFF"))
```

Use contrast checks when colors appear behind text, labels, or legends.

