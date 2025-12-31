#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-30 18:23:46
# @DESCRIPTION:
# @VERSION: v0.0.1


import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/data/scfoundation/out/gse_dataset_metadata_full.csv"
)
SRR = pl.read_csv(SRR_FILENAME, ignore_errors=True)
SRRHEAD = SRR.select(["gseid"]).unique()
basedir = "/mnt/isilon/u01_project/large-scale/liuc9/raw"


def process_row(row):
    gseid = row["gseid"]
    # gsedir = row["gsedir"]
    cmd_arr = [
        "/scr1/users/liuc9/tools/miniforge3/envs/renv/bin/Rscript",
        "/home/liuc9/github/scMOCHA-data/src/06.1-collect-variants-new.R",
        "-g",
        gseid,
        "-b",
        basedir,
    ]
    # print(f"Processing {gseid} in {basedir} ...")
    print(" ".join(cmd_arr))
    # try:
    #     subprocess.run(cmd_arr)
    # except Exception as e:
    #     print(f"Error processing {gseid}: {e}")


with ProcessPoolExecutor(max_workers=5) as executor:
    executor.map(process_row, SRRHEAD.iter_rows(named=True))
