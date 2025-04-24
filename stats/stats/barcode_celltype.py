#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-24 01:27:47
# @DESCRIPTION:
# @VERSION: v0.0.1


import concurrent.futures
import functools
import logging
import time
from pathlib import Path

import duckdb
import pandas as pd

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pd.read_csv(SRR_FILENAME)
SRR

dfs = []
for _, row in SRR.iterrows():
    gseid = row["gseid"]
    srrid = row["srrid"]
    srrdir = row["srrdir"]
    df = pd.read_csv(
        f"{srrdir}/barcode_cluster.tsv",
        sep="\t",
        header=None,
        names=["barcode", "cj", "celltype"],
        dtype={"barcode": str, "celltype": str},
        usecols=["barcode", "celltype"],
    )
    df["gseid"] = gseid
    df["srrid"] = srrid
    dfs.append(df)


df = pd.concat(dfs, ignore_index=True)
# Reorder columns
df = df[["gseid", "srrid", "barcode", "celltype"]]
df.to_csv(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/barcode_celltype.csv",
    index=False,
)
