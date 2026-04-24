#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-19 12:31:19
# @DESCRIPTION:
# @VERSION: v0.0.1



import pysam

vcffile = "/mnt/isilon/xing_lab/liuc9/refdata/gtexv8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"


vcf = pysam.VariantFile(vcffile)

"chrM" in list(vcf.header.contigs)

for rec in vcf.fetch("chrM", 0, 1000):
    print(rec)
