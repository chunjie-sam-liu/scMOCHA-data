conn_all_hetero_af <- DBI::dbConnect(
  duckdb::duckdb(),
  "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
)
DBI::dbListTables(conn_all_hetero_af)


variants_disease_df <- data.frame(
  short = c(
    "3571insC",
    "3664G>A",
    "3916G>A",
    "4063G>A",
    "4142G>A",
    "5910G>A",
    "5949G>A",
    "6124T>C",
    "6253T>C",
    "6340C>T",
    "6567C>T",
    "6663A>G",
    "6924G>T",
    "8932C>T",
    "10398A>G",
    "11778G>A",
    "11872insC",
    "12425delA",
    "13802C>T",
    "13937delAC",
    "14429delG",
    "15342insT",
    "3243A>G",
    "8344A>G",
    "8993T>G"
  ),
  detailed = c(
    "OC m.3571insC",
    "PC, OC m.3664G>A",
    "OC m.3916G>A",
    "OC m.4063G>A",
    "PC m.4142G>A",
    "PC m.5910G>A",
    "PC m.5949G>A",
    "PC m.6124T>C",
    "PC m.6253T>C",
    "PC 6340C>T",
    "OC m.6567C>T",
    "PC m.6663A>G",
    "PC m.6924G>T",
    "PC m.8932C>T",
    "PC m.10398A>G",
    "LHON m.11778G>A",
    "OC m.11872insC",
    "OC m.12425delA",
    "PC m.13802C>T",
    "OC m.13937delAC",
    "OC m.14429delG",
    "OC m.15342insT",
    "MELAS syndrome m.3243A>G",
    "MERRF syndrome m.8344A>G",
    "NARP m.8993T>G"
  ),
  note = c(
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "LHON",
    "",
    "",
    "",
    "",
    "",
    "",
    "MELAS syndrome",
    "MERRF syndrome",
    "NARP"
  ),
  stringsAsFactors = FALSE
)

variants_disease_new_df <- data.frame(
  Gene = c(
    "MT-ND1",
    "MT-ND1",
    "MT-ND1",
    "MT-ND1",
    "MT-ND1",
    "MT-ND1",
    "MT-ND1",
    "MT-ND1",
    "MT-ND3",
    "MT-ND4",
    "MT-ND5",
    "MT-ND6",
    "MT-ND6",
    "MT-ND6",
    "MT-CYB",
    "MT-CO1",
    "MT-CO1",
    "MT-ATP6",
    "MT-ATP6",
    "MT-ATP6"
  ),
  Variant = c(
    "m.3365T>C",
    "m.3376G>A",
    "m.3460G>A",
    "m.3481G>A",
    "m.3796A>G",
    "m.3946G>A",
    "m.3949T>C",
    "m.4175G>A",
    "m.10134C>A",
    "m.11778G>A",
    "m.13513G>A",
    "m.14459G>A",
    "m.14484T>C",
    "m.14487T>C",
    "m.14864T>C",
    "m.6708G>A",
    "m.6930G>A",
    "m.8993T>G",
    "m.8993T>C",
    "m.9176T>C"
  ),
  Complex = c(
    rep("Complex I", 14),
    "Complex III",
    rep("Complex IV", 2),
    rep("Complex V", 3)
  ),
  Findings_Context = c(
    "Found in a single individual (VAF ≥60%) with residual activity <1% to 36%.",
    "Found in a single individual (VAF ≥60%) with residual activity <1% to 36%.",
    "A primary LHON variant included in the 69 Complex I cases reviewed.",
    "Found in a single individual (VAF ≥60%) with residual activity <1% to 36%.",
    "Found in a single individual (VAF ≥60%) with residual activity <1% to 36%.",
    "Co-occurred in a single patient with m.3949 and m.4175.",
    "Co-occurred in a single patient with m.3946 and m.4175.",
    "Specifically noted in a case study (co-occurring with m.3946/m.3949). VAF ≥60% in muscle with severe biochemical defect.",
    "High VAF (near homoplasmy) showed ~17% residual activity; noted for protein loss.",
    "Primary LHON variant. The only variant found at 100% VAF in lymphoblasts with normal activity.",
    "Reviewed for threshold; cases with >90% VAF sometimes retained high residual activity.",
    "Reported in 4 cases; all had VAF >60% with residual muscle activity 36%–72%.",
    "Primary LHON variant included in the review.",
    "Described as one of the most commonly reported MT-ND6 variants (N=4 cases).",
    "The specific variant cited for the Complex III cases (p.Cys40Arg).",
    "Noted to likely alter assembly/stability of Complex IV (along with m.6930).",
    "Introduces a premature stop codon; analyzed for assembly defects.",
    "Most frequent (N=13). Strong correlation (VAF vs Activity) in muscle (τ=-0.58) but not fibroblasts.",
    "Included in Complex V analysis; generally milder phenotype than T>G.",
    "Reviewed in 4 cases; high VAF (>95%) in fibroblasts showed variable activity (74% to >100%)."
  ),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(short = gsub("m\\.", "", Variant))


mt_nd1_ <- c("3460G>A", "3697G>A", "3634A>G", "3380G>A")
disease_variant <- c(
  variants_disease_df$short,
  variants_disease_new_df$short
)

dplyr::tbl(
  conn_all_hetero_af,
  "allvariants"
) |>
  dplyr::filter(variant %in% disease_variant)
