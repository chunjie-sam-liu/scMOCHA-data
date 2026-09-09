---
name: r-figure
description: Write R plotting code in this repository, with ggplot2 or any other R graphics system. Use whenever generating, editing, or reviewing R code that produces a figure - histogram, scatter, forest, bar, heatmap, Manhattan, QQ, PCA, KM/CIF survival curve, Sankey/alluvial, UpSet, Venn, ComplexHeatmap, or a multi-panel assembly - or when choosing a theme, axis label, subtitle, figure size, output path, or save call. Covers the object-returning vs device-drawing families and which save path each needs, the shared fn_theme contract, saveplot instead of ggsave, package loading with load_pkg, colors from the track color file, glue subtitles that carry N and thresholds, log-axis flooring, patchwork and multi-page PDF assembly, and the ASCII-source rule.
---

# R Figure

Figures here are read by collaborators without the script next to them. A plot
is finished when someone can tell what was plotted, on how many observations,
and against which threshold, from the image alone.

ggplot2 is the default, but it is not the only system in use. Decide the family
first, because it determines how the figure is saved.

## Decide the family first

| Family              | Produces                        | Save with                          | Used here for                               |
| ------------------- | ------------------------------- | ---------------------------------- | ------------------------------------------- |
| A. Object-returning | a `ggplot` / `patchwork` object | `saveplot(file, p, width, height)` | ggplot2, ggsankey, ggupset, ggvenn, ggrepel |
| B. Device-drawing   | side effects on an open device  | `pdf(...)`; draw; `dev.off()`      | ComplexHeatmap, grid, base graphics         |
| C. Hybrid           | a list containing a `ggplot`    | pull the element, then family A    | survminer `ggsurvplot()` -> `p$plot`        |

**`saveplot()` only handles family A.** Passing a `Heatmap` object to it does
not produce the figure you expect. A ComplexHeatmap must be wrapped in an
explicit device block. Getting this wrong is the most common failure in figure
code.

This skill is portable and names no project value. The track color file, its
object names, the stage config filename, and the available PDF tooling come
from the repository's `.github/instructions/` bindings or `AGENTS.md`.

## Non-negotiables

1. **Family A saves with `saveplot()`, never `ggsave()`.** `saveplot()` is the
   jutils wrapper; it creates the directory, defaults to 300 dpi and a white
   background, and accepts a list of plots for a multi-page PDF.
2. **Family B always closes its device.** Wrap in `on.exit(dev.off())` or keep
   the `pdf()` / `dev.off()` pair adjacent so an error cannot leave the device
   open and the file truncated.
3. **PDF is the default device.** Use PNG only when the consumer needs a raster
   (slide deck, GitHub comment) or when a Manhattan plot has too many points.
4. **Every ggplot uses the track's theme helper.** The canonical helper is
   `fn_theme()`, defined once in the stage `config.R`. A track that predates
   that contract may bind a different helper name in the repository bindings;
   use whatever the bindings name for the track the script lives in. Never
   call a bare `theme_bw()` where a helper exists. Systems with their own
   theme (`theme_sankey()`) keep it and add the stage title styling on top.
5. **No raw hex in a plotting script.** Every color comes from the track
   `color.R`, including one-off status, highlight, and annotation colors.
   Neutral chrome greys (`"grey70"` for a reference line, `"grey80"` for
   `na.value`) are the only inline exception. Load `palette` before
   introducing a new one.
6. **Source files are ASCII.** Write `\u00b7`, `\u2265`, `\u2013`, not the
   literal characters. If a rendered label contains one of those glyphs, save
   with `device = cairo_pdf` - the default `pdf()` device silently substitutes
   `>=` for `\u2265` and only warns.
7. **`log_info()` after every save**, naming the file and the row count it was
   built from.

## Script skeleton

```r
#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: <Name>
# @CONTACT: <email>
# @DATE: YYYY-MM-DD
# @DESCRIPTION: Figures for stage NN: one line per panel, saying what each
#               panel shows and what it is evidence for.
# @VERSION: v0.1.0

# Reproducibility ----------------------------------------------------------

# Library ------------------------------------------------------------------
suppressMessages({
  library(jutils)
  load_pkg(ggsankey, ggrepel, scales)   # only the extras this script needs
})

# Args ---------------------------------------------------------------------

# Logger -------------------------------------------------------------------
log_layout(layout_glue_colors)
log_threshold(INFO)

# Load data ----------------------------------------------------------------
dotenv()                            # the track's env file
set.seed(9527)
suppressMessages({
  conflicted::conflicts_prefer(dplyr::filter, fs::path)
})

# Source -------------------------------------------------------------------
source("src/NN-stage/config.R")     # stage constants, fn_theme(), path helper
source("src/color.R")               # the track color file, every color used

paths <- stage_paths()
dir_create(paths$figdir)
```

The section order is the one every step script uses (Metainfo /
Reproducibility / Library / Args / Logger / Load data / Source / ...; see
`analysis-pipeline/references/script-templates.md` section E). A figure script
is a step script that happens to draw.

Attaching `jutils` already brings in `ggplot2`, `patchwork`, `paletteer`,
`prismatic`, `glue`, `fs`, `data.table`, and `logger`. Do not add `library()`
calls for any of them. Everything else goes through `load_pkg()` with unquoted
names: `load_pkg(ComplexHeatmap, circlize)`.

Without `jutils`, `ggsave()` is the fallback for family A, but then the output
directory, `dpi`, and `bg` must all be set explicitly, and a list of plots has
to be looped by hand.

## The theme

Define `fn_theme()` once in the stage `config.R` and source it. If the stage
has no config, define it at the top of the figure script and keep it identical
to this. When the repository bindings name a different helper for the track
(a stage that predates `config.R`, for example), match that helper's existing
body instead of introducing a second one beside it:

```r
fn_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        color = "grey40",
        size = 10
      )
    )
}
```

Add `fn_theme() + theme(legend.position = "none")` for per-panel deviations.
Never edit `fn_theme()` to suit one panel.

Base size by output: `12` for a standalone figure, `10-11` for a panel inside a
patchwork grid, `7-9` for a dense multi-facet diagnostic sheet.

## Title and subtitle carry the numbers

The title says what is plotted. The subtitle carries N, thresholds, and the one
caveat a reader needs. Build it with `glue()` from live objects so it can never
go stale.

```r
labs(
  title = "Minor allele frequency at the discovery leads",
  subtitle = glue::glue(
    "{nrow(d)} leads \u00b7 {n_zero} monomorphic, floored for display ",
    "\u00b7 dashed line at MAF = {THRESH}"
  ),
  x = "MAF (log scale; monomorphic floored for display)",
  y = "leads",
  fill = NULL
)
```

Never hard-code a count into a title string. Use `\u00b7` as the separator
between subtitle clauses. Set `fill = NULL` / `color = NULL` when the legend
keys are self-explanatory; a `qc_status` legend header helps nobody.

For a device-drawing figure the same information goes in `column_title` and the
log line, since there is no subtitle slot.

## Axis labels

Statistical notation uses `expression()`, not text:

```r
y = expression(-log[10](italic(P)))
x = expression(Expected ~ -log[10](italic(P)))
x = expression(-log[10](P) ~ "  discovery (QC-filtered genotypes)")
```

Units belong in the axis label, not the tick labels:
`"LV Mass Index (g/m2)"`, `"Age (years)"`, `"Follow-up (years)"`.

## Colors

```r
scale_fill_manual(values = color_sex, na.value = "grey80")
scale_color_manual(values = color_ancestry)
scale_fill_gradient2(low = color_scatter_point, mid = "white",
                     high = color_highlight,
                     midpoint = 0, limits = c(-1, 1), na.value = "grey95")

# Device-drawing systems take a color *function*, not a scale
col_fun <- circlize::colorRamp2(c(0, 5, 10), color_ordered_ramp)
```

Every name above is illustrative. Use the objects the track `color.R` actually
defines, and define any missing one there first — never inline a hex here.

- Discrete categories: a named vector whose names match the factor levels
  exactly. A mismatch silently drops to grey.
- Continuous in ggplot2: `paletteer::paletteer_c()` or an explicit
  `scale_*_gradient2()` for signed values, always with `midpoint = 0`.
- Continuous in ComplexHeatmap: `circlize::colorRamp2()` with explicit breaks.
- Always set `na.value` (ggplot2) or an explicit NA color (ComplexHeatmap). The
  default grey reads as a real category.
- A status palette (`"kept by QC"` / `"removed by QC"`) is a color like any
  other: define it in `color.R` as `color_<concept>_status` and reference the
  name here. Scripts that still carry one locally are frozen history - leave
  them alone unless asked to migrate.

## Log axes and zeros

`scale_x_log10()` drops zeros without warning. Floor them and say so in the
axis label:

```r
# MAF = 0 cannot be drawn on a log axis; floor it at half the smallest non-zero
# value so monomorphic leads stay visible as their own stripe.
nz <- d[maf > 0, min(maf)]
floor_maf <- nz / 2
d[, maf_plot := pmax(maf, floor_maf)]
```

## Reference lines before points

Draw `geom_hline` / `geom_vline` / `geom_abline` before the data layer so the
data sits on top.

```r
geom_abline(slope = 1, intercept = 0, color = "grey70", linewidth = 0.3) +
geom_hline(yintercept = THRESH, linetype = "dashed",
           color = "grey30", linewidth = 0.35) +
geom_point(alpha = 0.6, size = 1.5)
```

Conventions: identity/null line solid `grey60`-`grey70` at `linewidth = 0.3`; a
threshold dashed `grey30`-`grey40` at `linewidth = 0.35`; a genome-wide line in
the highlight color. Every dashed line must be named in the subtitle.

## Bars with value labels

Expand the axis or the labels are clipped:

```r
geom_col(width = 0.7) +
geom_text(aes(label = glue("{n}\n({round(pct, 1)}%)")), vjust = -0.3, size = 3.5) +
scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
```

For a horizontal bar chart use `hjust = -0.05` with `coord_flip()` and
`limits = c(0, max * 1.35)`.

## Point labels

Use `ggrepel::geom_text_repel()` and label only the subset worth naming:

```r
ggrepel::geom_text_repel(
  data = hits[tier == "genome_wide"],
  aes(label = variant_id),
  size = 2.2,
  min.segment.length = 0,
  max.overlaps = 20
)
```

Labelling every point makes the figure unreadable; filter first.

## Beyond ggplot2

Recipes with working code are in
[references/non_ggplot_recipes.md](references/non_ggplot_recipes.md). Summary of
what to reach for:

| Need                             | Package                       | Family | Note                                               |
| -------------------------------- | ----------------------------- | ------ | -------------------------------------------------- |
| Flow between categorical stages  | `ggsankey`                    | A      | `make_long()` + `geom_sankey()` + `theme_sankey()` |
| Set intersections                | `ggupset`                     | A      | list-column of set names + `scale_x_upset()`       |
| 2-3 set overlap                  | `ggvenn`                      | A      | more than 3 sets: use UpSet instead                |
| Annotated matrix / locus x trait | `ComplexHeatmap` + `circlize` | B      | `pdf()` / `draw()` / `dev.off()`                   |
| Survival curve, KM or CIF        | `survminer`                   | C      | `ggsurvplot()$plot`, then patchwork                |
| Second color scale on one plot   | `ggnewscale`                  | A      | `new_scale_color()` between layers                 |

Do not introduce a new plotting package when one of these covers the need. If
one genuinely does not, say so and name the replacement before writing code.

## Assembling panels

`patchwork` for a fixed layout, and scale the saved size with the panel count:

```r
p <- Reduce(`+`, plots) + patchwork::plot_layout(nrow = 1)
saveplot(as.character(paths$figdir / "rg_heatmap.pdf"), p,
         width = 5 * length(plots), height = 5)
```

A list of plots for a multi-page PDF, one page per plot:

```r
saveplot(
  filename = as.character(outdir / "overlap_venn.pdf"),
  plot = list(p_ab, p_ac, p_abc),
  device = "pdf",
  width = 6,
  height = 4
)
```

Related panels that share axes are one faceted plot, not a patchwork of
separate plots. Family B figures cannot be combined with patchwork; give each
its own file, or draw them onto successive pages of one `pdf()` device.

## Saving

```r
# Family A
saveplot(
  as.character(paths$figdir / "maf_by_qc_status.pdf"),
  p_hist,
  width = 8,
  height = 5
)
log_info("Saved maf_by_qc_status.pdf ({nrow(d)} leads)")

# Family B
pdf(as.character(paths$figdir / "discovery_landscape.pdf"),
    width = 9, height = max(7, 0.12 * nrow(M)))
ComplexHeatmap::draw(ht, annotation_legend_list = list(state_leg),
                     merge_legend = TRUE)
dev.off()
log_info("Wrote discovery_landscape.pdf ({nrow(M)} loci x {ncol(M)} traits)")
```

- `saveplot()` and `pdf()` both need a character path; wrap an `fs_path` in
  `as.character()`.
- Always pass `width` and `height` in inches. The default is the last device
  size and is not reproducible.
- Pass `device = cairo_pdf` whenever any rendered label carries a non-ASCII
  glyph. The default device does not fail, it silently substitutes.
- Starting sizes:

  | Figure                           | width x height                     |
  | -------------------------------- | ---------------------------------- |
  | Single scatter / histogram       | `7.5 x 5.5`, `8 x 5`               |
  | Forest, one variant per PDF page | `7.5 x 3.5`                        |
  | Forest, many rows on one page    | `9 x 6`                            |
  | ggplot heatmap panel             | `5 x 5` per panel, `coord_fixed()` |
  | ComplexHeatmap, row-scaled       | `9 x max(7, 0.12 * n_rows)`        |
  | UpSet                            | `9 x 5`                            |
  | Venn, multi-page                 | `6 x 4`                            |
  | Sankey                           | `16 x 10`, up to `26` wide         |
  | Manhattan (PNG)                  | `10 x 4`                           |
  | QQ (PNG)                         | `4.5 x 4.5`                        |
  | Multi-panel patchwork grid       | `12 x 8` to `16 x 10`              |

- When the row count drives the height, scale it with a floor:
  `height = max(7, 0.12 * nrow(mat))`, and `width = 5 * length(plots)` for a
  one-row patchwork.
- Filenames are lowercase, `_`-separated, and describe the content:
  `maf_by_qc_status.pdf`, `upset_sharing.pdf`, `discovery_landscape.pdf`.
- Figures go under the stage figure directory from the stage paths helper - see
  `data-result-layout`.

## Guard empty inputs

A figure script must not fail the stage because one upstream file is missing or
one panel has no rows.

```r
if (nrow(both) > 0) {
  ...
  saveplot(...)
  log_info("Saved p_old_vs_new.pdf ({nrow(both)} rows)")
} else {
  log_warn("No row is tested in both cohorts; skipping p_old_vs_new.pdf")
}
```

```r
read_if <- function(p) {
  p <- as.character(p)
  if (file.exists(p)) as.data.table(import(p, lazy = FALSE)) else NULL
}
```

Warn and skip, never write an empty PDF silently. For family B this matters
more, because an unguarded `pdf()` with nothing drawn still leaves a file on
disk that looks like a success.

## Before reporting success

```bash
ls -l --time-style=+%F_%T <figure dir>/
for f in <figure dir>/*.pdf; do
  printf '%s  %s\n' "$(pdfinfo "$f" | awk '/^Pages/{print $2" pages"}')" "$f"
done
```

`pdfinfo` (poppler) is the check to reach for when it is installed; otherwise
use the R `pdftools` package. Confirm which one the project has before relying
on it. A PDF that neither can parse is a truncated family-B file from a missing
`dev.off()`.

Confirm each expected file exists, is non-empty (a truncated PDF is often only
a few hundred bytes), reports the expected page count, has a mtime from this
run, and that any skipped panel was reported as skipped with its reason.

## References

Both reference files are self-contained: every block builds its own synthetic
input and runs as written, so a recipe can be executed before it is adapted.

- [references/worked_examples.md](references/worked_examples.md) - complete
  ggplot2 figures: log-axis histogram with floored zeros, log-log scatter,
  multi-page forest, labelled bars assembled with patchwork, Manhattan and QQ
  with point thinning.
- [references/non_ggplot_recipes.md](references/non_ggplot_recipes.md) -
  ggsankey, ggupset, ggvenn, ComplexHeatmap, survminer, ggnewscale.

Related skills: `jutils` for `saveplot()`, `load_pkg()`,
`fn_xy_breaks_limits()`, `human_read()`, `human_read_latex_pval()`; `palette`
before adding any color; `data-result-layout` to choose the output path.
