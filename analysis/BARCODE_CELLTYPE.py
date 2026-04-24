#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-24 01:27:47
# @DESCRIPTION:
# @VERSION: v0.0.1


from pathlib import Path

import polars as pl

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)
SRR

dfs = []
for row in SRR.iter_rows(named=True):
    gseid = row["gseid"]
    srrid = row["srrid"]
    srrdir = row["srrdir"]
    print(f"Processing {gseid} {srrid} {srrdir}")
    df = pl.read_csv(
        f"{srrdir}/barcode_cluster.tsv",
        separator="\t",
        has_header=False,
        new_columns=["barcode", "cj", "celltype"],
    )
    # Only keep barcode and celltype columns
    df = df.select(["barcode", "celltype"])
    df = df.with_columns(
        [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
    )
    # Reorder columns to ensure gseid and srrid are first and second
    df = df.select(["gseid", "srrid", "barcode", "celltype"])
    dfs.append(df)


df = pl.concat(dfs)
# Reorder columns
df = df.select(["gseid", "srrid", "barcode", "celltype"])
df.write_csv(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/barcode_celltype.csv"
)
