color_disease <- c(
  "Alzheimer's Disease" = "#BC3C29FF",
  "COVID-19" = "#E18727FF",
  "Healthy" = "#0072B5FF",
  "Unknown" = "grey50",
  "Other" = "grey"
)

color_chemistry <- c(
  "SC5P-PE" = "#440154FF",
  "SC5P-R2" = "#31688EFF",
  "SC3Pv3" = "#35B779FF",
  "SC3Pv2" = "#FDE725FF"
)

color_gender <- c(
  "Female" = "#EE0000FF",
  "Male" = "#3B4992FF",
  "Unknown" = "grey50"
)

color_celltype <- c(
  "B" = "#66C2A5FF",
  "CD4 T" = "#FC8D62FF",
  "CD8 T" = "#8DA0CBFF",
  "other T" = "#E5C494FF",
  "NK" = "#FFD92FFF",
  "DC" = "#E78AC3FF",
  "Mono" = "#A6D854FF",
  "other" = "grey50"
)

color_circos_track <- c(
  "phastCons100way" = "#FFD700",
  "gnomad" = "#000000",
  "homo_paf" = "#ae00ff",
  "hete_paf" = "#00b3ff",
  "somatic_paf" = "#FF0000",
  "homo_af" = "#ae00ff",
  "hete_af" = "#00b3ff",
  "somatic_af" = "#FF0000",
  "gene_name_bg" = "#EAF7FF",
  "coverage" = "#3FB5FF"
)

color_mtdna_gene_type <- c(
  "D-Loop" = "#FFFFB3",
  "MT tRNA" = "#FB8072",
  "MT rRNA" = "#8DD3C7",
  "Protein coding" = "#BEBADA",
  "MT OLR" = "#FDB462"
)

color_variantcell <- c(
  "Heteroplasmy" = "red",
  "Sufficient reads" = "darkblue",
  "No sufficient reads" = "gray",
  "No reads" = "white"
)

color_haplogroup <- c(
  "A" = "#C1A72F",
  "B" = "#FAD2D9",
  "C" = "#ED2891",
  "D" = "#F6B667",
  "F" = "#104A7F",
  "H" = "#9EDDF9",
  "HV" = "#3953A4",
  "I" = "#007EB5",
  "J" = "#B2509E",
  "K" = "#97D1A9",
  "L" = "#ED1C24",
  "M" = "#F8AFB3",
  "N" = "#EA7075",
  "R" = "#754C29",
  "T" = "#D49DC7",
  "U" = "#CACCDB",
  "V" = "#D3C3E0",
  "W" = "#A084BD",
  "X" = "#542C88",
  "Y" = "#D97D25",
  "Z" = "#6E7BA2"
)

# Variant-caller comparison (newplots/compare-with-mgatk) --------------------
# Orange/blue anchors taken from the ggsci NEJM palette, the same source as
# color_disease; the pair stays separable under clr_deutan().
color_caller <- c(
  "Original mgatk" = "#E18727",
  "scMOCHA" = "#0072B5"
)

color_variant_set <- c(
  "Both" = "#7F7F7F",
  "Original mgatk only" = "#E18727",
  "scMOCHA only" = "#0072B5"
)

# Three calling arms. The scMOCHA pair is one hue at two lightnesses because
# they are the same caller before and after its downstream AF gate;
# prismatic::clr_darken("#0072B5", 0.35) gives the darker one.
color_arm <- c(
  "Original mgatk" = "#E18727",
  "scMOCHA variant call" = "#0072B5",
  "scMOCHA AF>5%" = "#044A77"
)

# Why a variant is absent from an arm. ggsci NEJM anchors, same source as
# color_disease and color_caller. "Retained" carries the NEJM red because it is
# the one outcome a reader looks for first in 07d; every exclusion reason is
# therefore a non-red hue, and the four that co-occur in 07d (retained, strand
# r only, VMR and strand r, VMR only) are mutually distinguishable.
color_exclusion <- c(
  "Retained" = "#BC3C29",
  "Not proposed as candidate" = "#BEBEBE",
  "< 3 confident cells" = "#7876B1",
  "VMR and strand r" = "#6F99AD",
  "VMR only" = "#EE4C97",
  "strand r only" = "#FFDC91",
  "Blacklisted position" = "#20854E",
  "< 10 cells at AF >= 0.05" = "#7F7F7F"
)

color_cell_inclusion <- c(
  "Kept by both" = "#7F7F7F",
  "Dropped by original mgatk" = "#E18727"
)

color_mgatk_gate <- c(
  "Passes mgatk VMR/strand gate" = "#7F7F7F",
  "Rejected by mgatk VMR/strand gate" = "#BC3C29"
)

# paletteer::paletteer_c("scico::batlow", n = 6); hard-coded so that sourcing
# this file stays free of a runtime palette dependency.
color_af_bin <- c(
  "#001959",
  "#185461",
  "#577646",
  "#B28D2E",
  "#FAA588",
  "#F9CCF9"
)

# The eight samples compared in newplots/compare-with-mgatk. Each takes its
# chemistry's anchor from color_chemistry so a sample's colour stays tied to
# its library type. Five samples share SC3Pv3, so they are separated along a
# lightness ladder off that anchor: prismatic::clr_lighten("#35B779", 0.45) and
# 0.20, the anchor itself, then clr_darken() at 0.28 and 0.52. Under
# clr_deutan() the ladder stays monotonic in lightness, so the five remain
# distinguishable.
color_sample <- c(
  "GSE149689_GSM4509019_3PV3" = "#6CE1A3",
  "GSE155673_GSM4712895_3PV3" = "#50C98B",
  "GSE163314_GSM4976997_3PV2" = "#FDE725",
  "GSE163668_GSM4995445_5PR2" = "#31688E",
  "GSE181279_GSM5494116_5PPE" = "#440154",
  "GSE188632_GSM5687372_3PV3" = "#35B779",
  "GSE220189_GSM6793474_3PV3" = "#058250",
  "GSE271107_GSM8369876_3PV3" = "#025734"
)

color_xlsx_hdr <- "#4D4D4D"
color_xlsx_note <- "#7F7F7F"
color_xlsx_white <- "#FFFFFF"
