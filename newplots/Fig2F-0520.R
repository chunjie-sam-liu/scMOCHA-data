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

somatic_srrids <- SOMATIC_VARIANTS$srrid
somatic_variants <- SOMATIC_VARIANTS$variant

somatic_variants |> unique() |> length()

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
  geom_text(
    aes(label = paste0("n=", scales::comma(nmutated))),
    vjust = -0.3,
    size = 3
  ) +
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

dt_all_hetero_af_cell |>
  filter(variant_type == "colorful") |>
  dplyr::distinct(srrid, celltype, barcode) |>
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
    # limits = c(0, 1),
    expand = expansion(mult = c(0.0, 0.05))
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
    aes(label = scales::percent(prop, accuracy = 0.1)),
    vjust = -0.5,
    size = 3,
    fontface = "bold"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    # limits = c(0, 1),
    expand = expansion(add = c(0.0, 0.005))
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
    expand = expansion(add = c(0.01, 0.05))
  ) +
  labs(
    title = "Prone to accumulate mtDNA somatic mutations (PBMC) Scale to 100%"
  ) -> p_bar_allcells

saveplot(
  filename = path(Sys.getenv("OUTDIRNOTUSE")) /
    "StJude-New" /
    "Fig2F-PBMC-ALL-CELLS-SOMATIC-MUTATED-CELLS-AND-PROPORTION.pdf",
  plot = list(
    p_allcells,
    p_mutatedcells,
    p_bar_allcells,
    p_bar_allcells_zoom
  ),
  width = 6,
  height = 4
)

# Save  --------------------------------------------------------------

# Session info -------------------------------------------------------------
if (isTRUE(verbose)) {
  sessionInfo()
}
