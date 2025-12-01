library(VariantAnnotation)

vcffile = "/mnt/isilon/xing_lab/liuc9/refdata/gtexv8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"


vcf <- readVcf(vcffile, "hg38")
