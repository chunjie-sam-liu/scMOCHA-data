library(Rsamtools)
library(GenomicAlignments)
bamfile <- "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-cell_cluster_annotation/execution/MT_cluster.bam"


bam <- BamFile(bamfile, asMates = TRUE)
flagstat <- scanBamFlag(isDuplicate = TRUE)

dup_reads <- scanBam(bam, flag = flagstat)


rg1 <- readr::read_lines(
  "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cell/temp/barcoded_bams/barcodes.1.sort.rg"
)

rg2 <- readr::read_lines(
  "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cell/temp/ready_bam/barcodes.1.qc.rg"
)
rgg <- setdiff(rg1, rg2)
head(rgg)


bcell_bamfile <- "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cluster/temp/temp_bam/B.temp1.barcoded.qc.bam"
