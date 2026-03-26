row1 <- import(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE171555/final/GSM5227148/variant_info_from_heatmap.qs"
) |>
  filter(variant == "961T>C")
# GSE206283 - GSM6249274
row2 <- import(
  "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE206283/final/GSM6249274/variant_info_from_heatmap.qs"
) |>
  filter(variant == "961T>C")


row1

row2
