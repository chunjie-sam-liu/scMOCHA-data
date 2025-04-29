#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-28 13:35:56
# @DESCRIPTION:
# @VERSION: v0.0.1


import logging
import math
import os
import shutil
from enum import Enum
from pathlib import Path
from typing import Annotated, Literal, Optional

import polars as pl
import typer
from matplotlib import pyplot as plt
from rich import print
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")

BASEDIR = Path("/home/liuc9/github/scMOCHA-data/data")
HQVDIR = BASEDIR / "high_quality_variant"
SRRFILENAME = Path(
    "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRRFILENAME)


def load_hetero_df(
    srrdir: Path,
    cl: str = "cell",
) -> pl.DataFrame:
    """
    Load the heteroplasmic variant files.
    """
    hetero_filename = srrdir / f"{cl}.cell_heteroplasmic_df.tsv.gz"
    df = pl.read_csv(hetero_filename, separator="\t")
    log.info(f"Loaded hetero_df: {hetero_filename}")
    return df


def plot_vmr_strand(
    df: pl.DataFrame, cutoff_vmr: float, cutoff_strand_correlation: float
):
    # Create a scatter plot of vmr_log vs strand_correlation
    plt.figure(figsize=(10, 6))
    plt.scatter(
        df["strand_correlation"].to_numpy(), df["vmr_log"].to_numpy(), alpha=0.5
    )
    plt.axvline(x=cutoff_strand_correlation, color="r", linestyle="--")
    plt.axhline(y=math.log10(cutoff_vmr), color="r", linestyle="--")
    plt.xlabel("Strand Correlation")
    plt.ylabel("VMR (log10)")
    plt.title(
        f"Strand Correlation vs VMR (Cutoffs: VMR={cutoff_vmr}, Strand Corr={cutoff_strand_correlation})"
    )


def load_stats(
    srrdir: Path,
    cl: str = "cell",
    cutoff_vmr: float = 0.01,
    cutoff_strand_correlation: float = 0.3,
) -> pl.DataFrame:
    """
    Load the stats of the high quality variant files.
    """
    cell_level_filename = srrdir / f"{cl}.variant_stats.tsv.gz"
    if not cell_level_filename.exists():
        print(f"[red]File not found: {cell_level_filename}[/red]")
        return None
    df = pl.read_csv(cell_level_filename, separator="\t").with_columns(
        vmr_log=pl.col("vmr").log10(),
    )

    df_f = df.filter(
        (pl.col("vmr_log") > math.log10(cutoff_vmr))
        & (pl.col("strand_correlation") > cutoff_strand_correlation)
    )
    log.info(
        f"Loaded stats_df: {cell_level_filename} with cutoff_vmr: {cutoff_vmr}, cutoff_strand_correlation: {cutoff_strand_correlation}"
    )
    return df_f


def load_variant_anno(
    srrdir: Path,
) -> pl.DataFrame:
    """
    Load the variant annotation files.
    """
    variant_anno_filename = srrdir / "variant_annotation.tsv"
    if not variant_anno_filename.exists():
        print(f"[red]File not found: {variant_anno_filename}[/red]")
        return None
    df = pl.read_csv(variant_anno_filename, separator="\t").with_columns(
        variant=pl.concat_str(
            [pl.col("Position"), pl.col("Ref"), pl.lit(">"), pl.col("Alt")]
        )
    )
    log.info(f"Loaded variant_anno_df: {variant_anno_filename}")
    return df


def hqv(
    srrdir: Path,
    cl: str = "cell",
    cutoff_vmr: float = 0.01,
    cutoff_strand_correlation: float = 0.3,
) -> None:
    """
    Load the high quality variant files.
    """
    # hetero_df = load_hetero_df(srrdir, cl)
    stats_df = load_stats(srrdir, cl, cutoff_vmr, cutoff_strand_correlation)
    anno_df = load_variant_anno(srrdir)

    anno_df_filter = anno_df.filter(
        pl.col("variant").is_in(stats_df["variant"])
    ).select(pl.col("variant"))

    anno_df_filter.write_csv(
        srrdir / "high_quality_variant.tsv",
        separator="\t",
        include_header=True,
    )
    log.info(
        f"Saved high quality variant file: {srrdir / 'high_quality_variant.tsv'}"
    )


app = typer.Typer()


class Cluster(str, Enum):
    cell = "cell"
    cluster = "cluster"


@app.command()
def check_stats(cluster: Cluster = Cluster.cell):
    """
    Check the stats of the high quality variant files.
    """
    cl = cluster.value
    for row in SRR.iter_rows(named=True):
        # row = next(SRR.iter_rows(named=True))
        gseid = row["gseid"]
        srrid = row["srrid"]
        srrdir = HQVDIR / gseid / "final" / srrid
        cell_level = srrdir / f"{cl}.variant_stats.tsv.gz"
        if not cell_level.exists():
            print(f"[red]File not found: {cell_level}[/red]")
            continue


@app.command()
def generate_hqv(
    srrdir: Path,
    cluster: Cluster = Cluster.cell,
    cutoff_vmr: float = 0.01,
    cutoff_strand_correlation: float = 0.3,
):
    """
    Generate high-quality variants based on the specified parameters.
    This function acts as a wrapper for the hqv function, applying appropriate
    cutoffs for variance-mean ratio and strand correlation to identify reliable variants.
    Args:
        srrdir (Path): Path to the directory containing sequencing data.
        cluster (Cluster, optional): Type of cluster to analyze. Defaults to Cluster.cell.
        cutoff_vmr (float, optional): Cutoff threshold for variance-mean ratio. Defaults to 0.01.
        cutoff_strand_correlation (float, optional): Cutoff threshold for strand correlation. Defaults to 0.3.
    Returns:
        None: Results are saved to disk by the underlying hqv function.
    """

    cl = cluster.value
    hqv(
        srrdir,
        cl,
        cutoff_vmr=cutoff_vmr,
        cutoff_strand_correlation=cutoff_strand_correlation,
    )


@app.command()
def create_sh():
    """
    Creates a shell script that generates high-quality variant commands for each SRR in the dataset.
    The function:
    1. Defines the output shell script path at HQVDIR / "update-high_quality_variant.sh"
    2. Removes any existing script file if present
    3. For each row in the SRR table:
        - Extracts GEO series ID and SRR ID
        - Constructs the SRR directory path
        - Creates a command to run the high-quality variant generation script with parameters:
         * Using cell clustering
         * VMR cutoff of 0.01
         * Strand correlation cutoff of 0.3
    4. Writes each command to the shell script
    The generated script can be executed to process all SRR datasets in batch.
    """

    sh_filename = HQVDIR / "update-high_quality_variant.sh"
    sh_filename.unlink(missing_ok=True)
    with open(sh_filename, "a") as f:
        for row in SRR.iter_rows(named=True):
            # row = next(SRR.iter_rows(named=True))
            gseid = row["gseid"]
            srrid = row["srrid"]
            srrdir = HQVDIR / gseid / "final" / srrid
            # srrdir.resolve()

            cmds = [
                "/scr1/users/liuc9/tools/anaconda3/envs/renv/bin/python3.13",
                "/home/liuc9/github/scMOCHA-data/stats/high_quality_variant/02-update-high-quality-variant.py",
                "generate-hqv",
                str(srrdir),
                "--cluster",
                "cell",
                "--cutoff-vmr",
                "0.01",
                "--cutoff-strand-correlation",
                "0.3",
            ]
            cmd = " ".join(cmds)
            f.write(cmd + "\n")


if __name__ == "__main__":
    app()
