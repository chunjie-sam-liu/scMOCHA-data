library(jutils)
load_pkg(qs)

thevariant <- "4175G>A"
base_dir <- "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/MANUSCRIPTFIGURES-notuse/deg/4175G>A/deg_merge_vaf"
files <- list.files(base_dir, pattern = "\\.qs$", full.names = TRUE)

for (f in files) {
  message("Checking file: ", f)
  d <- tryCatch({
    import(f)
  }, error = function(e) {
    message("Error importing: ", e$message)
    return(NULL)
  })
  
  if (is.null(d)) next
  
  sig <- d[d$p_val_adj < 0.05 & abs(d$avg_log2FC) > 0.25, ]
  message("Total rows: ", nrow(d))
  message("Significant DEGs (FDR < 0.05, |log2FC| > 0.25): ", nrow(sig))
  
  if (nrow(sig) > 0) {
    # Check if they map to ENTREZID
    library(clusterProfiler)
    library(org.Hs.eg.db)
    ids <- bitr(rownames(sig), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    message("Mapped to ENTREZID: ", nrow(ids))
  }
}
