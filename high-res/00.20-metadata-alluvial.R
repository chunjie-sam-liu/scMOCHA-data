dotenv(".env")
basedir <- path(Sys.getenv("DATADIR"))
# outdir <- "/home/liuc9/github/scMOCHA-data/analysis/zzz/plot-basic"
outdir <- path(Sys.getenv("OUTDIR"))
cleandatadir <- path(Sys.getenv("CLEANDATADIR"))

# sex_pred <- import(
#   "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir_sex.qs"
# ) |>
#   dplyr::select(
#     srrid,
#     sex_pred = sex
#   )

gse_dataset_metadata_full <- import(
  cleandatadir / "gse_dataset_metadata_full.qs"
) |>
  dplyr::mutate(
    Gender = SEXPRED
  )


source(path(
  Sys.getenv("HIGHRESDIR"),
  "00-colors.R"
))
# ! Sankey plot for samples --------------------------------------------------------------------

# meta --------------------------------------------------------------------

gse_dataset_metadata_full |>
  # dplyr::filter(!gseid %in% gseids_tobe_excluded) |>
  # dplyr::filter(gseid != "GSE220189") |>
  dplyr::select(
    gseid,
    srrid,
    Race,
    Ethnicity,
    Gender,
    Age_group,
    disease,
    Chemistry
  ) -> gse_dataset_metadata_full_selected

gse_dataset_metadata_full_selected |>
  dplyr::group_by(Gender) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  dplyr::mutate(Gender_str = glue::glue("{Gender}\n(n={n})")) |>
  dplyr::mutate(
    Gender_str = factor(Gender_str, levels = Gender_str)
  ) -> Gender_str

gse_dataset_metadata_full_selected |>
  dplyr::group_by(Race) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  dplyr::mutate(Race_str = glue::glue("{Race}\n(n={n})")) |>
  dplyr::mutate(Race_str = factor(Race_str, levels = Race_str)) -> Race_str


gse_dataset_metadata_full_selected |>
  dplyr::group_by(Ethnicity) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  dplyr::mutate(Ethnicity_str = glue::glue("{Ethnicity}\n(n={n})")) |>
  dplyr::mutate(
    Ethnicity_str = factor(Ethnicity_str, levels = Ethnicity_str)
  ) -> Ethnicity_str

gse_dataset_metadata_full_selected |>
  dplyr::group_by(disease) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  dplyr::mutate(Disease_str = glue::glue("{disease}\n(n={n})")) |>
  dplyr::mutate(
    Disease_str = factor(Disease_str, levels = Disease_str)
  ) -> Disease_str

gse_dataset_metadata_full_selected |>
  dplyr::group_by(Chemistry) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  dplyr::mutate(Chemistry_str = glue::glue("{Chemistry}\n(n={n})")) |>
  dplyr::mutate(
    Chemistry_str = factor(Chemistry_str, levels = Chemistry_str)
  ) -> Chemistry_str

gse_dataset_metadata_full_selected |>
  dplyr::mutate(
    Age = ifelse(Age_group == "Unknown", "Unknown", "Known")
  ) |>
  dplyr::group_by(Age) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  dplyr::mutate(Age_str = glue::glue("{Age}\n(n={n})")) |>
  dplyr::mutate(Age_str = factor(Age_str, levels = Age_str)) -> Age_str


gse_dataset_metadata_full_selected |>
  dplyr::mutate(
    Age = ifelse(Age_group == "Unknown", "Unknown", "Known")
  ) |>
  dplyr::select(-gseid, -srrid, -Age_group) |>
  dplyr::group_by(Chemistry, Age, Gender, Race, Ethnicity, disease) |>
  dplyr::count() |>
  dplyr::ungroup() |>
  dplyr::left_join(
    Chemistry_str,
    by = c("Chemistry")
  ) |>
  dplyr::left_join(
    Age_str,
    by = c("Age")
  ) |>
  dplyr::left_join(
    Gender_str,
    by = "Gender"
  ) |>
  dplyr::left_join(
    Race_str,
    by = "Race"
  ) |>
  dplyr::left_join(
    Ethnicity_str,
    by = "Ethnicity"
  ) |>
  dplyr::left_join(
    Disease_str,
    by = "disease"
  ) -> for_sankey_plot

library(ggalluvial)

chem_levels <- c("SC3Pv2", "SC3Pv3", "SC5P-R2", "SC5P-PE") |> rev()
ggsci::pal_aaas()(3) |> prismatic::color()
chem_colors <- viridis::viridis_pal(option = "D")(4) |>
  prismatic::color()

for_sankey_plot |>
  dplyr::mutate(
    Chemistry = factor(
      Chemistry,
      levels = chem_levels
    )
  ) |>
  ggplot(aes(
    axis1 = Chemistry_str,
    axis2 = Age_str,
    axis3 = Gender_str,
    axis4 = Race_str,
    axis5 = Ethnicity_str,
    axis6 = Disease_str,
    y = n.x
  )) +
  ggalluvial::geom_alluvium(
    aes(fill = Chemistry),
    width = 1 / 12,
    alpha = 0.8,
  ) +
  ggalluvial::geom_stratum(
    # aes(fill = Chemistry),
    width = 0.3
  ) +
  scale_fill_manual(values = color_chemistry) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum))
  ) +
  scale_x_discrete(
    limits = c(
      "Chemistry_str",
      "Age_str",
      "Gender_str",
      "Race_str",
      "Ethnicity_str",
      "Disease_str"
    ),
    labels = gsub(
      "_str",
      "",
      c(
        "Chemistry_str",
        "Age_str",
        "Sex_str",
        "Race_str",
        "Ethnicity_str",
        "Disease_str"
      )
    ),
    expand = c(0.2, 0.05)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  # ggsci::scale_fill_aaas() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    legend.position = "top"
  ) +
  guides(
    fill = guide_legend(title = "Sequencing", nrow = 1)
  ) -> meta_plot_sankey

meta_plot_sankey

ggsave(
  file.path(outdir, "NUMBER-OF-SAMPLES-WITH-METADATA.pdf"),
  meta_plot_sankey,
  width = 15,
  height = 8
)
