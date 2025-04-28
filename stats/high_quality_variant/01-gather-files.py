#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-28 13:06:24
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

app = typer.Typer()


@app.command()
def create_sh():
    sh_filename = HQVDIR / "all_srrid.sh"
    for row in SRR.iter_rows(named=True):
        gseid = row["gseid"]
        srrdir = Path(row["srrdir"])
        outdir = HQVDIR / gseid / "final"
        outdir.mkdir(parents=True, exist_ok=True)

        cmd = f"rsync -av  --partial --progress {srrdir} {outdir}/"
        # print(cmd)
        with open(sh_filename, "a") as f:
            f.write(cmd + "\n")
    print(f"sh file created: {sh_filename}")


if __name__ == "__main__":
    app()
