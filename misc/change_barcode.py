#!/usr/bin/env python

import os
import sys

import pysam

bamfile = sys.argv[1]
bamfile = "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cluster/temp/temp_bam/B.temp1.bam"


base = os.path.basename(bamfile)
outfolder = os.path.dirname(bamfile)
basename = os.path.splitext(base)[0]
mtchr = "MT"
barcodeTag = "CB"
umitag = "UB"


def getBarcode(read, tag_get):
    """
    Parse out the barcode per-read
    """
    # for tg in read.tags:
    # 	if(tag_get == tg[0]):
    # 		return(tg[1])
    # return("AA")
    # Using get_tag to get the results
    try:
        return read.get_tag(tag_get)
    except Exception:
        return "AA"


bam = pysam.AlignmentFile(bamfile, "rb")
outname = outfolder + "/" + basename + ".barcoded" + ".bam"
out = pysam.AlignmentFile(outname, "wb", template=bam)

# Make a DNA inspired additional barcode to account for potentially different sample indices
bases = "ACGT"
fauxdon = [
    a + b + c + d for a in bases for b in bases for c in bases for d in bases
]

# Filter for reads that match the set of possible barcodes for this sample
try:
    Itr = bam.fetch(str(mtchr), multiple_iterators=False)
    for read in Itr:
        barcode_id = getBarcode(read, barcodeTag)

        # Now check for true UMI
        if umitag != "XX":
            umi_id = getBarcode(read, umitag)
        else:
            umi_id = ""

        # Make a fake UMI from 1) cell barcode + 2) captured umi + 3) experiment
        # all with just ACGTs so that picard doesn't bark at us.
        # Only do this if the last string element is a number (i.e channel in 10x convention)
        if barcode_id[-1].isnumeric():
            split_two = barcode_id.split("-")
            faux_umi = split_two[0] + umi_id + fauxdon[(int(split_two[1]) - 1)]
        else:
            faux_umi = barcode_id + umi_id
        read.tags = read.tags + [("MU", faux_umi)]
        out.write(read)


except OSError:  # Truncated bam file from previous iteration handle
    print("Finished parsing bam")

bam.close()
out.close()

pysam.index(outname)
