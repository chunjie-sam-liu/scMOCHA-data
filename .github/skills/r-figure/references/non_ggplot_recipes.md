# Non-ggplot2 R figure recipes

Recipes for R graphics systems beyond plain ggplot2. Each block is
self-contained: it builds its own synthetic input, so it runs as written and
can be executed before being adapted. Every recipe here was executed before
being written down.

Each recipe names its family, because that decides the save path:

- **A** returns an object -> `saveplot(file, p, width, height)`
- **B** draws on a device -> `pdf(...)`; draw; `dev.off()`
- **C** returns a list containing a ggplot -> pull the element, then family A

Colors are written as literal hexes below only so each block runs standalone.
In real code every one of them comes from the track `color.R`. -> `palette`

Shared preamble:

```r
suppressMessages({
  library(jutils)
  load_pkg(ggsankey, ggupset, ggvenn, ComplexHeatmap, circlize,
           survminer, survival, ggnewscale, scales)
})

fn_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5, color = "grey40", size = 10
      )
    )
}
```

---

## ggsankey - flow between categorical stages (family A)

`make_long()` reshapes one row per subject into the node / next_node long form.
Counting nodes before plotting is what lets each box carry its own `n`.

```r
set.seed(1)
cohort <- tibble::tibble(
  sex = factor(sample(c("Female", "Male"), 400, TRUE),
               levels = c("Female", "Male")),
  dose = factor(sample(c("0", "<250", "\u2265250", "Missing"), 400, TRUE),
                levels = c("0", "<250", "\u2265250", "Missing")),
  outcome = factor(sample(c("Case", "Control"), 400, TRUE),
                   levels = c("Case", "Control"))
)

fn_sankey_plot <- function(data, cols, title, palette = NULL, label_size = 3.5) {
  df_long <- data |> make_long(!!!syms(cols))

  node_counts <- df_long |>
    filter(!is.na(node)) |>
    group_by(x, node) |>
    tally(name = "n") |>
    ungroup() |>
    mutate(label = paste0(node, "\n(n=", comma(n), ")"))

  df_long <- df_long |> left_join(node_counts, by = c("x", "node"))

  p <- ggplot(
    df_long,
    aes(x = x, next_x = next_x, node = node, next_node = next_node,
        fill = factor(node), label = label)
  ) +
    geom_sankey(flow.alpha = 0.4, node.color = "grey30", show.legend = FALSE) +
    geom_sankey_label(size = label_size, color = "white",
                      fill = "grey40", hjust = 0.5) +
    labs(title = title, x = NULL) +
    theme_sankey(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text.x = element_text(size = 13, face = "bold"),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )

  if (!is.null(palette)) {
    p + scale_fill_manual(values = palette, na.value = "grey80")
  } else {
    p + scale_fill_viridis_d(option = "H", alpha = 0.85, na.value = "grey80")
  }
}

p <- fn_sankey_plot(cohort, c("sex", "dose", "outcome"), "Cohort flow")

# Node labels contain a non-ASCII glyph, so the default pdf device is not enough
saveplot("sankey.pdf", p, width = 16, height = 10, device = cairo_pdf)
```

Notes:

- `theme_sankey()` replaces `fn_theme()` here; keep it and add only the title
  and axis-text styling on top.
- Every stage column must be a `factor` with explicit levels, or the node order
  changes between runs. Bin continuous variables first and put `"Missing"`
  last.
- Sankeys need width: `16 x 10` for three to five stages, up to `26` wide for
  more.
- `scales::comma()` for the counts, so `n=1234` reads `n=1,234`.

---

## ggupset - set intersections (family A)

The input is a data frame with a **list-column** of set names, one row per
item. `scale_x_upset()` turns that column into the intersection matrix axis.

```r
lookup <- data.table::data.table(
  locus = rep(paste0("L", 1:60), each = 3),
  trait = sample(c("EF", "GLS", "LAVI", "CM"), 180, TRUE),
  P = runif(180, 1e-9, 1e-3)
)

mem <- unique(lookup[P < 1e-5, .(locus, trait)])
sets <- mem[, .(traits = list(sort(unique(as.character(trait))))), by = locus]
sets[, degree := lengths(traits)]

p <- ggplot(sets, aes(x = traits)) +
  geom_bar(fill = "#E64B35") +
  scale_x_upset(n_intersections = 25) +
  labs(
    x = "Trait combination",
    y = "Shared loci",
    title = "Locus sharing across traits (P < 1e-5)"
  ) +
  theme_bw()

saveplot("upset.pdf", p, width = 9, height = 5)
```

Notes:

- `sort()` inside the list-column matters: without it `c("EF","CM")` and
  `c("CM","EF")` become two different intersections.
- Cap with `n_intersections` or the axis becomes unreadable.
- Use UpSet, not Venn, above three sets.
- `theme_bw()` rather than `fn_theme()`: `scale_x_upset()` draws the
  combination matrix in the axis area, and blanking the panel grid makes it
  harder to read.

---

## ggvenn - two or three set overlap (family A)

`ggvenn()` returns a plain ggplot, so several can go to one multi-page PDF.

```r
a <- paste0("S", 1:600)
b <- paste0("S", 400:900)
cc <- paste0("S", 700:1200)

p_ab <- ggvenn::ggvenn(
  list("All" = a, "Group A" = b),
  fill_color = c("blue", "red"),
  stroke_size = 0.5, set_name_size = 4
)
p_abc <- ggvenn::ggvenn(
  list("All" = a, "Group A" = b, "Group B" = cc),
  fill_color = c("blue", "red", "#374E55"),
  stroke_size = 0.5, set_name_size = 4
)

saveplot("venn.pdf", list(p_ab, p_abc), device = "pdf", width = 6, height = 4)
```

`fill_color` is positional, matching the order of the input list.

---

## ComplexHeatmap - annotated matrix (family B)

This is the one system `saveplot()` cannot handle. It must be wrapped in an
explicit device block.

```r
set.seed(3)
nr <- 40
traits <- c("EF", "GLS", "LAVI", "CM")
M <- matrix(runif(nr * 4, 0, 9), nrow = nr,
            dimnames = list(paste0("L", seq_len(nr)), traits))
M[sample(length(M), 30)] <- NA

meta <- data.table::data.table(
  status = sample(c("replicated", "not tested"), nr, TRUE),
  gw = runif(nr) < 0.15,
  label = paste0("GENE", seq_len(nr))
)

SUGG_LOGP <- 5
COL_NA <- "#F2F2F2"
COL_GRAY <- "#BDBDBD"
STATUS_COLORS <- c(replicated = "#00A087", `not tested` = "#CCCCCC")
TRAIT_COLORS <- c(EF = "#3C5488", GLS = "#4DBBD5",
                  LAVI = "#00A087", CM = "#E64B35")

col_fun <- circlize::colorRamp2(
  c(SUGG_LOGP, max(M, na.rm = TRUE)), c("white", "#CB181D")
)

# cell_fun + rect_gp = gpar(type = "none") is how a cell gets a color rule that
# a single continuous scale cannot express.
cell_fun <- function(j, i, x, y, ww, hh, fill) {
  v <- M[i, j]
  col <- if (is.na(v)) COL_NA else if (v < SUGG_LOGP) COL_GRAY else col_fun(v)
  grid::grid.rect(x, y, ww, hh, gp = grid::gpar(fill = col, col = "grey90"))
}

ra <- rowAnnotation(
  Replication = meta$status,
  `Genome-wide` = ifelse(meta$gw %in% TRUE, "yes", "no"),
  col = list(
    Replication = STATUS_COLORS[intersect(names(STATUS_COLORS),
                                          unique(meta$status))],
    `Genome-wide` = c(yes = "#E64B35", no = "#EEEEEE")
  ),
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size = unit(3, "mm")
)

# anno_mark labels only the rows worth naming; row names stay hidden.
rma <- NULL
at <- which(meta$gw %in% TRUE)
if (length(at) > 0) {
  rma <- rowAnnotation(
    gene = anno_mark(at = at, labels = meta$label[at], side = "right",
                     labels_gp = gpar(fontsize = 6))
  )
}

M_clu <- M
M_clu[is.na(M_clu)] <- 0        # clustering cannot see NA; the color still can

ht <- Heatmap(
  M_clu,
  name = "-log10(P)",
  col = col_fun,
  rect_gp = gpar(type = "none"),
  cell_fun = cell_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9,
                         col = unname(TRAIT_COLORS[traits]),
                         fontface = "bold"),
  column_names_rot = 45,
  left_annotation = ra,
  right_annotation = rma,
  column_title = glue("Landscape ({nr} loci x {length(traits)} traits)"),
  use_raster = nrow(M_clu) > 200,
  raster_quality = 2,
  heatmap_legend_param = list(title = "-log10(P)")
)

state_leg <- Legend(
  labels = c("P > 1e-5 (tested)", "absent / NA"),
  legend_gp = gpar(fill = c(COL_GRAY, COL_NA)),
  title = "state",
  border = "black"
)

pdf("landscape.pdf", width = 9, height = max(7, 0.12 * nr))
draw(ht, annotation_legend_list = list(state_leg), merge_legend = TRUE)
dev.off()
log_info("Wrote landscape.pdf ({nrow(M)} loci x {ncol(M)} traits)")
```

Notes:

- `column_title` is the only title slot; there is no subtitle. Put the N there
  or in the log line.
- A discrete state the continuous scale cannot express (one grey for "tested
  but null", another for "absent") needs its own `Legend()` passed through
  `annotation_legend_list` with `merge_legend = TRUE`.
- Set `use_raster = TRUE` above ~200 rows or the PDF becomes huge.
- Height scales with rows: `height = max(7, 0.12 * nrow(M))`.
- `pdf()` and `dev.off()` stay adjacent. If anything between them can error,
  put `on.exit(dev.off())` right after `pdf()`.

---

## survminer - KM and CIF curves (family C)

`ggsurvplot()` returns a list, not a ggplot. The plot is `p$plot`.

```r
set.seed(4)
n <- 400
df <- data.frame(
  entry_age = runif(n, 5, 20),
  sex = factor(sample(c("Female", "Male"), n, TRUE),
               levels = c("Female", "Male"))
)
df$exit_age <- df$entry_age + rexp(n, 0.03)
df$event <- as.integer(runif(n) < 0.25)

color_map <- c(Female = "#F39B7F", Male = "#8491B4")

fit <- survfit(Surv(entry_age, exit_age, event) ~ sex, data = df)

# Recover the level order survfit actually used; do not assume alphabetical.
strata_vals <- sub("^.*=", "", names(fit$strata))
color_map <- color_map[strata_vals]

p <- ggsurvplot(
  fit,
  data = df,
  pval = FALSE,           # see the counting-process note below
  pval.method = FALSE,
  conf.int = TRUE,
  risk.table = FALSE,     # left truncation: build it manually below
  surv.median.line = "hv",
  palette = unname(color_map),
  break.time.by = 10,
  title = "Event-free survival by sex",
  xlab = "Age (years)",
  ylab = "Event-free survival",
  legend.title = NULL,
  legend.labs = strata_vals,
  ggtheme = fn_theme() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          legend.position = "top")
)

cox_p <- summary(
  coxph(Surv(entry_age, exit_age, event) ~ sex, data = df)
)$sctest[["pvalue"]]
p$plot <- p$plot + labs(subtitle = glue("score-test P = {signif(cox_p, 3)}"))

# n at risk under delayed entry = entered before a AND still at risk after a
x_breaks <- seq(0, 80, by = 10)
risk_tbl <- tibble::tibble(
  age = x_breaks,
  n_risk = sapply(x_breaks, \(a) sum(df$entry_age <= a & df$exit_age > a))
)

risk_plot <- ggplot(risk_tbl, aes(x = age, y = 1, label = n_risk)) +
  geom_text(size = 5) +
  scale_x_continuous(breaks = x_breaks, limits = c(-5, max(x_breaks)),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  labs(title = "Number at risk", x = "Age (years)", y = "strata") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_km <- p$plot / risk_plot + patchwork::plot_layout(heights = c(3, 1))
saveplot("km.pdf", p_km, width = 9, height = 7)
```

Three traps, all of which this recipe already avoids:

- **`pval = TRUE` errors on counting-process data.** With
  `Surv(entry, exit, event)` (delayed entry / age scale), `ggsurvplot()` calls
  `survdiff()`, which raises
  `survdiff not defined for counting process data` and kills the script. Set
  `pval = FALSE` and take the P from a `coxph()` score test instead.
- **The built-in risk table is wrong under left truncation**, because it does
  not subtract subjects who have not entered yet. Build it by hand.
- **`palette` and `legend.labs` are positional.** Recover the order with
  `sub("^.*=", "", names(fit$strata))` and index the color map with it.

The composed object is a patchwork, so it saves as family A.

---

## ggnewscale - two color scales on one plot (family A)

Used when a background reference cloud and the foreground points need
independent palettes.

```r
ref <- data.frame(PC1 = rnorm(600), PC2 = rnorm(600),
                  pop = sample(c("AFR", "EUR", "EAS"), 600, TRUE))
study <- data.frame(PC1 = rnorm(80, 0.5), PC2 = rnorm(80, -0.3),
                    grp = sample(c("Discovery", "Replication"), 80, TRUE))

p <- ggplot() +
  geom_point(data = ref, aes(PC1, PC2, color = pop), size = 0.4, alpha = 0.5) +
  scale_color_manual(
    values = c(AFR = "#009E73", EUR = "#0072B2", EAS = "#CC79A7"),
    name = "Reference"
  ) +
  ggnewscale::new_scale_color() +
  geom_point(data = study, aes(PC1, PC2, color = grp), size = 1.4) +
  scale_color_manual(
    values = c(Discovery = "#E64B35", Replication = "#F39B7F"),
    name = "Study"
  ) +
  labs(title = "Study samples projected onto the reference PCA") +
  fn_theme()

saveplot("pca_newscale.pdf", p, width = 7.5, height = 5.5)
```

`new_scale_color()` goes between the two layer/scale pairs. Order matters: the
background layer and its scale come first, or the reference cloud is drawn on
top of the study points.
