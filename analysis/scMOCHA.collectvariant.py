#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-10 16:42:01
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
    srrdir = row["srrdir"]
    # Change to srrdir and run the R script
    os.chdir(srrdir)
    subprocess.run(
        [
            "Rscript",
            "/home/liuc9/github/scMOCHA-data/analysis/scMOCHA.collectvariant.R",
        ]
    )


with ProcessPoolExecutor(max_workers=20) as executor:
    executor.map(process_row, SRRHEAD.iter_rows(named=True))
