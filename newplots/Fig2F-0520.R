#!/usr/bin/env Rscript
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-05-20 13:11:48
# @DESCRIPTION: this script is used for ...

# Reproducibility ----------------------------------------------------------
set.seed(1)
# Library -----------------------------------------------------------------

suppressMessages({
  library(jutils)
})

# Args --------------------------------------------------------------------

# s: string, i: integer, f: float, !: boolean, @: array, %: list
VERSION = "v0.0.1"

GetoptLong.options(help_style = "two-column")

# default: default value specified here.

nthread = 8
GetoptLong(
  "nthread=i",
  "Number of threads to use",
  "verbose",
  "Enable verbose logging"
)


# Logger ------------------------------------------------------------------

log_layout(layout_glue_colors)

if (isTRUE(verbose)) {
  log_threshold(TRACE)
  log_info("Verbose mode enabled")
} else {
  log_threshold(INFO)
}


# Load data ---------------------------------------------------------------
library(jutils)
dotenv()
suppressMessages({
  conflicted::conflicts_prefer(dplyr::filter, fs::path)
})


source(path(
  Sys.getenv("HIGHRESDIR"),
  "00-colors.R"
))
outdir <- path(Sys.getenv("OUTDIR"))

METAFULL <- import(outdir / "SAMPLES-METADATA-FULL.xlsx")
ALLVARIANTS <- import(
  outdir / "SAMPLE-VARIANT-CLASSIFICATION-CLUSTER-BULK-AF.xlsx"
) |>
  dplyr::mutate(
    coord = parallel::mclapply(
      X = variant,
      FUN = \(.v) {
        # .v <- gse_data_variant_classification_clusteraf_bulkaf$variant[[1]]
        pos <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\1",
          x = .v
        ) |>
          as.integer()
        ref <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\2",
          x = .v
        )
        alt <- gsub(
          pattern = "(\\d+)([A-Z])>([A-Z])",
          replacement = "\\3",
          x = .v
        )
        data.table(
          seqnames = "MT",
          start = pos,
          end = pos,
          ref = ref,
          alt = alt
        )
      },
      mc.cores = 10
    )
  ) |>
  tidyr::unnest(
    cols = coord
  )

SOMATIC_VARIANTS <- ALLVARIANTS |>
  dplyr::filter(variant_type == "somatic")

somatic_srrids <- SOMATIC_VARIANTS$srrid |> unique()
somatic_variants <- SOMATIC_VARIANTS$variant |> unique()

somatic_variants |> unique() |> length()

SOMATIC_VARIANTS |>
  select(gseid, srrid, variant) |>
  distinct() |>
  mutate(
    gseid_srrid_variant = paste0(gseid, "_", srrid, "_", variant)
  ) -> somatic_gseid_srrid_variant

# Source ---------------------------------------------------------------------

# Conn ---------------------------------------------------------------
conn <- db_conn(
  Sys.getenv("DUCKDB_PATH")
)
# Function ----------------------------------------------------------------

# Main --------------------------------------------------------------------

tbl_ls(conn)

tbl_all_hetero_af_cell <- tbl(conn, "all_hetero_af_cell")
tbl_allvariants_af_cell <- tbl(conn, "allvariants_af_cell")
# tbl_allvariants_af_cluster <- tbl(conn, "allvariants_af_cluster")

tbl_all_hetero_af_cell <- tbl(conn, "allvariants_cell") |>
  filter(srrid %in% somatic_srrids) |>
  filter(variant %in% somatic_variants)

dt_all_hetero_af_cell <- tbl_all_hetero_af_cell |> as.data.table()

dt_all_hetero_af_cell |>
  count(variant, celltype, variant_type) |>
  pivot_wider(
    names_from = variant_type,
    values_from = n
  ) |>
  mutate(
    celltype = gsub("_", " ", celltype)
  ) |>
  mutate(
    celltype = factor(celltype, levels = names(color_celltype))
  ) |>
  select(
    `Cell type` = celltype,
    Heteroplasmy = colorful,
    `Sufficient reads` = black,
    `No sufficient reads` = grey,
    `No reads` = white,
    Variant = variant
  ) |>
  mutate(
    `Exist variant` = if_else(
      !is.na(Heteroplasmy) &
        !is.na(`Sufficient reads`) &
        Heteroplasmy + `Sufficient reads` >= 10 &
        Heteroplasmy >= 2,
      1,
      0
    )
  ) -> table_s5

writexl::write_xlsx(
  table_s5,
  path(Sys.getenv("OUTDIRNOTUSE")) / "StJude-New" / "Table S5 - PBMC.xlsx"
)


table_s5 |>
  select(Variant, `Cell type`, `Exist variant`) |>
  mutate(
    `Exist variant` = if_else(
      `Exist variant` == 1,
      "Present",
      "Absent"
    )
  ) |>
  ggplot(aes(
    x = `Cell type`,
    y = Variant,
    fill = `Exist variant`
  )) +
  geom_tile() +
  scale_fill_manual(
    name = "Status",
    values = c(
      "Present" = "#DE7D64",
      "Absent" = "#AED6E5"
    )
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) -> p_heatmtap

saveplot(
  filename = path(Sys.getenv("OUTDIRNOTUSE")) /
    "StJude-New" /
    "Fig2F - PBMC - HEATMAP.pdf",
  plot = p_heatmtap,
  width = 6,
  height = 26
)


table_s5 |>
  select(Variant, `Cell type`, `Exist variant`) |>
  group_by(
    `Cell type`
  ) |>
  summarize(
    nmutated = sum(`Exist variant`, na.rm = TRUE),
    total = n(),
    prop = nmutated / total
  ) |>
  arrange(`Cell type`) |>
  ggplot(aes(
    x = `Cell type`,
    y = prop,
    fill = `Cell type`
  )) +
  geom_bar(
    stat = "identity",
    # fill = color_celltype,
  ) +
  scale_fill_manual(
    values = color_celltype,
    name = "Cell type"
  ) +
  geom_text(
    aes(label = scales::percent(prop, accuracy = 0.1)),
    vjust = -1.6,
    size = 3,
    fontface = "bold"
  ) +
  # geom_text(
  #   aes(label = paste0("n=", scales::comma(nmutated))),
  #   vjust = -0.3,
  #   size = 3
  # ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(add = c(0.01, 0.05))
  ) +
  labs(
    x = "Cell type",
    y = "Proportion",
    title = "Occurence of somatic variants across cell types (PBMC)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) -> p_bar


saveplot(
  filename = path(Sys.getenv("OUTDIRNOTUSE")) /
    "StJude-New" /
    "Fig2F - PBMC - BAR.pdf",
  plot = p_bar,
  width = 6,
  height = 4
)


tbl_ls(conn)
tbl_barcode <- tbl(conn, "barcode") |>
  select(srrid, celltype = celltype_name, barcode) |>
  distinct()


tbl_barcode |>
  count(celltype) |>
  as.data.table() -> allcells

dt_all_hetero_af_cell <- tbl(conn, "allvariants_cell") |>
  filter(srrid %in% somatic_srrids) |>
  filter(variant %in% somatic_variants) |>
  as.data.table()

# dt_all_hetero_af_cell |> filter(!srrid %in% somatic_srrids)

dt_all_hetero_af_cell |>
  filter(variant_type == "colorful") |>
  mutate(
    gseid_srrid_variant = paste0(gseid, "_", srrid, "_", variant)
  ) |>
  filter(
    gseid_srrid_variant %in% somatic_gseid_srrid_variant$gseid_srrid_variant
  ) |>
  mutate(srrid_barcode = paste0(srrid, "_", barcode)) -> mutated_cells_info


mutated_cells_info |>
  dplyr::distinct(srrid_barcode, srrid, celltype, barcode) |>
  count(celltype, name = "nmutated") |>
  as.data.table() -> mutated_cells


allcells |>
  left_join(mutated_cells, by = "celltype") |>
  mutate(nmutated = if_else(is.na(nmutated), 0L, nmutated)) |>
  mutate(prop = nmutated / n) -> celltype_mutation_summary


celltype_mutation_summary |>
  mutate(
    celltype = gsub("_", " ", celltype)
  ) |>
  mutate(
    celltype = factor(celltype, levels = names(color_celltype))
  ) -> forplot

forplot |>
  ggplot(aes(
    x = celltype,
    y = n,
    fill = celltype
  )) +
  geom_bar(
    stat = "identity",
    # fill = color_celltype,
  ) +
  scale_fill_manual(
    values = color_celltype,
    name = "Cell type"
  ) +
  geom_text(
    aes(label = scales::comma(n, accuracy = 1)),
    vjust = -0.5,
    size = 3,
    fontface = "bold"
  ) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 1),
    # limits = c(0, 1),
    expand = expansion(mult = c(0.0, 0.05))
  ) +
  labs(
    x = "Cell type",
    y = "Number of total cells",
    title = "# of total cells (PBMC)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) -> p_allcells


forplot |>
  ggplot(aes(
    x = celltype,
    y = nmutated,
    fill = celltype
  )) +
  geom_bar(
    stat = "identity",
    # fill = color_celltype,
  ) +
  scale_fill_manual(
    values = color_celltype,
    name = "Cell type"
  ) +
  geom_text(
    aes(label = scales::comma(nmutated, accuracy = 1)),
    vjust = -0.5,
    size = 3,
    fontface = "bold"
  ) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 1),
    limits = c(0, 5000),
    expand = expansion(mult = c(0.0, 0.0))
  ) +
  labs(
    x = "Cell type",
    y = "Number of mutated cells",
    title = "# of mutated cells (PBMC)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) -> p_mutatedcells

forplot |>
  ggplot(aes(
    x = celltype,
    y = prop,
    fill = celltype
  )) +
  geom_bar(
    stat = "identity",
    # fill = color_celltype,
  ) +
  scale_fill_manual(
    values = color_celltype,
    name = "Cell type"
  ) +
  geom_text(
    aes(label = scales::percent(prop, accuracy = 0.01)),
    vjust = -0.5,
    size = 3,
    fontface = "bold"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.01),
    limits = c(0, 0.01),
    expand = expansion(add = c(0.0, 0.0))
  ) +
  labs(
    x = "Cell type",
    y = "Proportion",
    title = "Prone to accumulate mtDNA somatic mutations (PBMC)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) -> p_bar_allcells_zoom

p_bar_allcells_zoom +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(add = c(0.01, 0.0))
  ) +
  labs(
    title = "Prone to accumulate mtDNA somatic mutations (PBMC) Scale to 100%"
  ) -> p_bar_allcells


mutated_cells_info |>
  count(
    celltype,
    srrid_barcode
  ) |>
  # mutated_cells_info_n |>
  mutate(
    celltype = gsub("_", " ", celltype)
  ) |>
  mutate(
    celltype = factor(celltype, levels = names(color_celltype))
  ) |>
  mutate(n = factor(n, levels = sort(unique(n), decreasing = TRUE))) |>
  ggplot(aes(x = celltype, fill = n)) +
  geom_bar() +
  scale_fill_manual(
    values = paletteer_dynamic("cartography::purple.pal", n = 5) |> rev(),
    name = "# of variants in a cell"
  ) +
  geom_text(
    data = \(.d) .d |> count(celltype, name = "total"),
    aes(x = celltype, y = total, label = scales::comma(total)),
    vjust = -0.5,
    size = 3,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 1),
    limits = c(0, 5000),
    expand = expansion(mult = c(0.0, 0.0))
  ) +
  labs(
    x = "Cell type",
    y = "Number of mutated cells",
    title = "# of mutated cells (PBMC)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),

    legend.title.position = "left",
    legend.title = element_text(face = "bold", angle = -90),
  ) -> p_bar_allcells_n_mutations

saveplot(
  filename = path(Sys.getenv("OUTDIRNOTUSE")) /
    "StJude-New" /
    "Fig2F-PBMC-ALL-CELLS-SOMATIC-MUTATED-CELLS-AND-PROPORTION-0527.pdf",
  plot = list(
    p_allcells,
    p_mutatedcells,
    p_bar_allcells_n_mutations,
    p_bar_allcells,
    p_bar_allcells_zoom
  ),
  width = 6,
  height = 4
)


mutated_cells_info |>
  count(
    celltype,
    srrid_barcode
  ) |>
  filter(n > 1) -> a

mutated_cells_info |>
  mutate(hasvariant = 1) |>
  filter(srrid_barcode %in% a$srrid_barcode) |>
  arrange(srrid_barcode) -> for_heatmap

for_heatmap |> count(variant) |> arrange(-n) -> rank_variant

for_heatmap |> count(srrid_barcode) |> arrange(-n) -> rank_cell

# for_heatmap |>
#   filter(variant == "12865A>G")

top5_variants <- rank_variant$variant[1:100]
variant_labels <- setNames(
  ifelse(rank_variant$variant %in% top5_variants, rank_variant$variant, ""),
  rank_variant$variant
)

n_variants <- nrow(rank_variant)
n_cells <- nrow(rank_cell)
plot_width <- max(8, n_variants * 0.25)
plot_height <- max(6, n_cells * 0.08)

for_heatmap |>
  mutate(
    variant = factor(variant, levels = rank_variant$variant),
    srrid_barcode = factor(srrid_barcode, levels = rank_cell$srrid_barcode)
  ) |>
  ggplot(
    aes(
      x = variant,
      y = srrid_barcode,
      fill = hasvariant
    )
  ) +
  geom_tile(color = NA) +
  scale_x_discrete(labels = variant_labels, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 20, r = 10, b = 5, l = 5)
  ) +
  labs(
    x = "Variant",
    y = "Cell"
  ) -> p_multi_variant_heatmap

saveplot(
  filename = path(Sys.getenv("OUTDIRNOTUSE")) /
    "StJude-New" /
    "Fig2F-PBMC-MULTI-VARIANT-HEATMAP-0527.pdf",
  plot = p_multi_variant_heatmap,
  width = plot_width,
  height = plot_height
)

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
