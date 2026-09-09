---
name: palette
description: Use whenever any color enters R code - choosing, deriving, standardizing, validating, or applying it - with paletteer, prismatic, ggplot2 scales, ComplexHeatmap color functions, openxlsx2 fills, and project-level color files such as src/color.R, src_<variant>/color.R, src/colors.R, or src/plot_colors.R. This includes a one-off status or highlight color, an Excel header or block fill, and a color added to a stage config, since every color belongs in the track's one color file rather than inline in the consuming script.
---

# Palette

Use this skill for R scientific plot color work. It is optimized for stable, reusable color choices in manuscripts, multi-panel figures, and collaborator-shared analysis projects.

This skill is portable and names no project value. Which color file each track
owns, which object names already exist, and which stage-local palettes are
frozen come from the repository's `.github/instructions/` bindings or
`AGENTS.md`.

## Default Workflow

1. Identify whether the plot needs discrete, continuous, or diverging colors.
2. Check for an existing project color file, in this order:
   - the track's own file: `<track>/color.R` -- `src/color.R`,
     `src_<variant>/color.R`, `pipeline/color.R`, `<name>_pipeline/color.R`.
     Use the one belonging to the track the plot lives in; the bindings list
     the tracks.
   - `src/colors.R`
   - `src/plot_colors.R`
3. If one of those files already exists, use it. If multiple exist, prefer the one already sourced by plotting scripts.
4. If no project color file exists, ask before creating one. Default new file: `src/color.R`.
5. Any figure, table, or workbook committed to the repository takes its colors
   as named objects from that file.
6. Inline `scale_*_paletteer_*()` calls are acceptable only in a throwaway
   scratch plot that is never committed.
7. Derive and validate the color before adding it; see "Derive, Do Not Invent".

## Project Color File Rules

- One color file per track is the single source of truth. Every color goes into
  that file first, chosen through this skill, and only then gets referenced by
  the consuming script. This covers figure palettes, one-off status and
  highlight colors, Excel header and block fills, and diagram colors. No
  category of color is exempt.
- Never write a raw hex anywhere except the color file: not in a plotting
  script, not in an export script, not in a stage `config.R` or the older
  `00-config.R`. A stage config that needs a palette aliases one
  (`COHORT_COLORS <- color_newgroup`) instead of defining a hex.
- Existing stage-local palettes are frozen history. Leave them exactly as they
  are - their figures were produced with those hexes, so migrating one silently
  changes a deliverable. Migrate only when the user asks for that specific
  migration.
- The one exception is meaningless chrome: neutral greys for `na.value`, a
  reference line, or a panel border. Those stay inline because they encode no
  category. A color that stands for a cohort, trait, ancestry, outcome, or QC
  state is never chrome.
- Preserve existing named color mappings unless the user explicitly asks to change them.
- Do not reorder categories silently.
- Use named vectors for discrete category colors.
- Use functions for continuous or diverging palettes.
- Names in discrete vectors must exactly match factor levels, metadata values, or annotation labels.
- Keep canonical biological identity colors separate from derived variants.

### Naming and loading

Match the convention the file already uses. A common one is a `color_` prefix
on every object - `color_sex`, `color_ancestry`, `color_highlight`. A new
object then follows the same shape: `color_<concept>` for a named vector,
`color_<concept>_<variant>` for a derived one. Read an existing object from the
file before adding one; do not assume a prefix.

Load the file belonging to the track the script lives in:

```r
source("src/color.R")                          # a src/ script
source("src_<variant>/color.R")                # a src_<variant>/ script
source("<track>/color.R")                      # any other track
source(fs::path(repodir, "src", "color.R"))    # same file, built from the env
```

In a multi-track repository the files are independent. They may hold the same
object names and values today, which makes a wrong-track `source()` appear to
work right up until one track's palette changes. Bind to your own track. A
literal path is repository-relative, so the script must be launched from the
repository root; the `fs::path(repodir, ...)` form works from anywhere once
`dotenv()` has run.

Recommended search before editing. Match the filename, not the call form:
scripts build the path with `fs::path()` or `here()` as often as they write a
literal string, and a literal-only regex misses those.

```bash
rg -n '(color|colors|plot_colors)\.R' --glob '*.R' --glob '*.sh' --glob '*.lsf'
```

## Derive, Do Not Invent

A new color is computed from something that already exists, never hand-picked.
In order of preference:

1. Derive from an anchor already in the color file: `clr_lighten()`,
   `clr_darken()`, `clr_alpha()`, `clr_desaturate()`.
2. Pull from a named `paletteer` palette, so the choice is reproducible and
   citable: `paletteer_d("ggsci::nrc_npg")`, `paletteer_c("scico::batlow")`.
3. Only when neither fits, choose a hex by hand, and record in the color file
   why that specific value was needed.

Validate before committing the color:

```r
prismatic::contrast_ratio(new_color, "white")            # >= 3 behind text
plot(prismatic::color(prismatic::clr_deutan(palette)))   # still separable
```

Two colors that must never be confused - two ancestries, a discovery versus a
replication cohort, a reference cloud under study points - are checked against
**each other** under `clr_deutan()`, not only against the background.

## Paletteer Guidance

Use `paletteer` to retrieve palettes and ggplot scales.

- Discrete fixed palette: `paletteer::paletteer_d("package::palette")`
- Dynamic palette: `paletteer::paletteer_dynamic("package::palette", n)`
- Continuous palette: `paletteer::paletteer_c("package::palette", n = 256)`
- ggplot exploration: `scale_color_paletteer_d()`, `scale_fill_paletteer_c()`, etc.
- final shared figures: prefer `scale_color_manual(values = object)` or `scale_fill_manual(values = object)`.

Read `references/paletteer_usage.md` for concrete examples and selection rules.

## Prismatic Guidance

Use `prismatic` as the companion package for color inspection and controlled transformations.

- Preview: `prismatic::color(x)` and `plot(prismatic::color(x))`
- Light/dark variants: `clr_lighten()`, `clr_darken()`
- Alpha variants: `clr_alpha()`
- Saturation variants: `clr_saturate()`, `clr_desaturate()`
- Color vision deficiency simulation: `clr_deutan()`, `clr_protan()`, `clr_tritan()`
- Contrast checks: `contrast_ratio()`, `best_contrast()`

Do not overwrite canonical identity colors with transformed colors. Store transformed variants under explicit names such as `celltype_colors_light` or `celltype_colors_alpha`.

Read `references/prismatic_usage.md` when transformations or validation are requested.

## Color File Pattern

Use this style for final/shared figures:

```r
color_celltype <- c(
  "B cell" = "#4E79A7",
  "T cell" = "#F28E2B",
  "Monocyte" = "#E15759"
)

color_celltype_light <- prismatic::clr_lighten(color_celltype, 0.25)
color_celltype_alpha <- prismatic::clr_alpha(color_celltype, 0.75)

color_expression <- function(n = 256) {
  paletteer::paletteer_c("scico::batlow", n = n)
}

color_signed_score <- function(n = 256) {
  paletteer::paletteer_c("scico::vik", n = n)
}
```

Read `references/color_file_patterns.md` for ggplot2, ComplexHeatmap, pheatmap, and Seurat usage patterns.

## Helper Scripts

Use scripts only when they save time or provide concrete validation.

- `scripts/inspect_paletteer.R`: list/search installed paletteer palettes.
- `scripts/preview_palette.R`: preview a palette and optionally save a PDF swatch.
- `scripts/audit_color_file.R`: source a color file, audit a named vector against expected levels, and optionally save colorblindness previews.

Prefer the repository's R runner. If the project owns `pixi.toml`, run scripts through that environment, for example:

```bash
pixi run Rscript /path/to/palette/scripts/inspect_paletteer.R --type discrete --search ggsci
```
