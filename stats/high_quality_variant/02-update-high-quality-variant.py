#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-28 13:35:56
# @DESCRIPTION:
# @VERSION: v0.0.1


import os
import shutil
from pathlib import Path

import polars as pl
import typer
from rich import print

BASEDIR = Path("/home/liuc9/github/scMOCHA-data/data")
HQVDIR = BASEDIR / "high_quality_variant"
SRRFILENAME = Path(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRRFILENAME)


def check_stats():
    """
    Check the stats of the high quality variant files.
    """
    for row in SRR.iter_rows(named=True):
        gseid = row["gseid"]
        srrid = row["srrid"]
        srrdir = HQVDIR / gseid / "final" / srrid
        cell_level = srrdir / "cell.variant_stats.tsv.gz"
        if not cell_level.exists():
            print(f"[red]File not found: {cell_level}[/red]")
            continue


if __name__ == "__main__":
    check_stats()
