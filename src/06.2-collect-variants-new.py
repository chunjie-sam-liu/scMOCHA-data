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
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)
SRRHEAD = (
    SRR.with_columns(
        pl.col("srrdir")
        .map_elements(lambda x: str(Path(x).parent.parent.parent))
        .alias("gsedir")
    )
    .select(["gseid", "gsedir"])
    .unique()
)


def process_row(row):
    gseid = row["gseid"]
    gsedir = row["gsedir"]
    subprocess.run(
        [
            "/scr1/users/liuc9/tools/anaconda3/envs/renv/bin/Rscript",
            "/home/liuc9/github/scMOCHA-data/src/06.1-collect-variants-new.R",
            "-g",
            gseid,
            "-b",
            gsedir,
        ]
    )


with ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(process_row, SRRHEAD.iter_rows(named=True))
