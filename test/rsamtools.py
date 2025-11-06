import pysam

bamfile = "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-cell_cluster_annotation/execution/MT_cluster.bam"

bam = pysam.AlignmentFile(bamfile, "rb")
Itr = bam.fetch(str("MT"), multiple_iterators=False)

bases = "ACGT"
fauxdon = [
    a + b + c + d for a in bases for b in bases for c in bases for d in bases
]
fauxdon


outdir = "~/github/scMOCHA-data/data/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cluster/"
sample = "B"
output_bam = outdir + "/temp/ready_bam/" + sample + ".qc.bam"
rmlog = output_bam.replace(".qc.bam", ".rmdups.log").replace(
    "/temp/ready_bam/", "/logs/rmdupslogs/"
)
filtlog = output_bam.replace(".qc.bam", ".filter.log").replace(
    "/temp/ready_bam/", "/logs/filterlogs/"
)
temp_bam0 = output_bam.replace(".qc.bam", ".temp0.bam").replace(
    "/temp/ready_bam/", "/temp/temp_bam/"
)
temp_bam1 = output_bam.replace(".qc.bam", ".temp1.bam").replace(
    "/temp/ready_bam/", "/temp/temp_bam/"
)
prefixSM = outdir + "/temp/sparse_matrices/" + sample
outputdepth = outdir + "/qc/depth/" + sample + ".depth.txt"

python = "python"
filtclip_py = "/scr1/users/liuc9/tools/mgatk/mgatk/bin/python/filterClipBam.py"
input_bam = "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cluster/temp/barcoded_bams/B.bam"
mito_genome = "MT"
proper_paired = "False"
NHmax = "1"
NMmax = "4"
output_bam
pycall = (
    " ".join(
        [
            python,
            filtclip_py,
            input_bam,
            filtlog,
            mito_genome,
            proper_paired,
            NHmax,
            NMmax,
        ]
    )
    + " > "
    + temp_bam0
)
max_javamem = "8000m"
script_dir = "/scr1/users/liuc9/tools/mgatk/mgatk/"
picardCall = (
    "java"
    + " -Xmx"
    + max_javamem
    + " -jar "
    + script_dir
    + "/bin/picard.jar MarkDuplicates"
)

umi_extra = " BARCODE_TAG=MU"
mdc_long = (
    picardCall
    + " I="
    + temp_bam1
    + " O="
    + output_bam
    + " M="
    + rmlog
    + " REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT QUIET=true VERBOSITY=ERROR USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
    + umi_extra
)

chunk_bam_py = script_dir + "/bin/python/chunk_barcoded_bam.py"
input = "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cluster/temp/temp_bam/B.temp1.bam"
bcbd = "/mnt/isilon/u01_project/large-scale/ting/raw/GSE235050/cromwell-executions/scMOCHABatch/66a50c5e-b7f0-4d1c-9c2f-9dec97c3a7f7/call-scMOCHA/shard-0/sub.scMOCHA/bf8e6c39-e9e6-4a4c-828e-af5935932eaf/call-call_mt_variants/execution/cluster/temp/temp_bam"
barcode_tag = ("CB",)

pycall = " ".join(
    [
        "python",
        chunk_bam_py,
        input,
        bcbd,
        barcode_tag,
        one_barcode_file,
        mito_chr,
        umi_barcode,
    ]
)
