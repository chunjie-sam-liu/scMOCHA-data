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
SRRHEAD = SRR


def process_row(row):
    gseid = row["gseid"]
    # srrid = row["srrid"]
    srrdir = Path(row["srrdir"])
    # Change to srrdir and run the R script
    # os.chdir(srrdir)
    basedir = srrdir.parent.parent.parent
    subprocess.run(
        [
            "/scr1/users/liuc9/tools/anaconda3/envs/renv/bin/Rscript",
            "/home/liuc9/github/scMOCHA-data/src/06.1-collect-variants-new.R",
            "-g",
            gseid,
            "-b",
            str(basedir),
        ]
    )


with ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(process_row, SRRHEAD.iter_rows(named=True))
