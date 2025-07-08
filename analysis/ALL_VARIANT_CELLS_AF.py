#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-24 18:27:13
# @DESCRIPTION:
# @VERSION: v0.0.1

import logging
import multiprocessing as mp
from functools import partial
from pathlib import Path
from typing import Annotated, Literal

import polars as pl
import typer
from rich import print
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)
SRR


ALL_VARIANT_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant.csv"
)
ALL_VARIANT_DF = pl.read_csv(ALL_VARIANT_FILENAME)
HETEROPLASMIC_DF = (
    ALL_VARIANT_DF.filter(pl.col("issomatic") == "heteroplasmic")
    .select(["variant", "Position"])
    .sort("Position")
)

HOMOPLASMIC_DF = (
    ALL_VARIANT_DF.filter(pl.col("issomatic") == "homoplasmic")
    .select(["variant", "Position"])
    .sort("Position")
)

CONFIG = {
    "HETEROPLASMIC": {
        "POSISTIONS": HETEROPLASMIC_DF["Position"].to_list(),
        "VARIANTS": HETEROPLASMIC_DF["variant"].to_list(),
        "ALTS": [
            variant.split(">")[1]
            for variant in HETEROPLASMIC_DF["variant"].to_list()
        ],
        "POSALTS": [
            f"{pos}_{alt}"
            for pos, alt in zip(
                HETEROPLASMIC_DF["Position"].to_list(),
                [
                    variant.split(">")[1]
                    for variant in HETEROPLASMIC_DF["variant"].to_list()
                ],
            )
        ],
        "table_path": Path(
            "/home/liuc9/github/scMOCHA-data/data/zzz/db/HETEROPLASMIC"
        ),
        "suffix": "hetero",
    },
    "HOMOPLASMIC": {
        "POSISTIONS": HOMOPLASMIC_DF["Position"].to_list(),
        "VARIANTS": HOMOPLASMIC_DF["variant"].to_list(),
        "ALTS": [
            variant.split(">")[1]
            for variant in HOMOPLASMIC_DF["variant"].to_list()
        ],
        "POSALTS": [
            f"{pos}_{alt}"
            for pos, alt in zip(
                HOMOPLASMIC_DF["Position"].to_list(),
                [
                    variant.split(">")[1]
                    for variant in HOMOPLASMIC_DF["variant"].to_list()
                ],
            )
        ],
        "table_path": Path(
            "/home/liuc9/github/scMOCHA-data/data/zzz/db/HOMOPLASMIC"
        ),
        "suffix": "homo",
    },
}

# Extract positions from HETEROPLASMIC_DF
# POSISTIONS = HETEROPLASMIC_DF["Position"].to_list()
# VARIANTS = HETEROPLASMIC_DF["variant"].to_list()
# ALTS = [variant.split(">")[1] for variant in VARIANTS]
# POSALTS = [f"{pos}_{alt}" for pos, alt in zip(POSISTIONS, ALTS)]

TABLEDIR = Path("/home/liuc9/github/scMOCHA-data/analysis/zzz/db/TABLES")


def load_table_pl(gseid, srrid, cluster, hh):
    # gseid, srrid, cluster, hh = "GSE147794", "GSM4446059", "cell", "HETEROPLASMIC"
    tablename = (
        f"{cluster}_cov.{gseid}_{srrid}"
        if cluster != "bulk"
        else f"cluster_cov.{gseid}_{srrid}"
    )
    log.info(f"Loading {tablename}")
    table_path = TABLEDIR / f"{tablename}.csv"
    df = pl.read_csv(
        table_path,
        has_header=True,
    ).select(
        ["base", "barcode", *[str(pos) for pos in CONFIG[hh]["POSISTIONS"]]]
    )

    df_alt_cov = (
        df.pivot(
            index="barcode",
            on="base",
            values=[str(pos) for pos in CONFIG[hh]["POSISTIONS"]],
            aggregate_function="sum",
        )
        .select(["barcode", *CONFIG[hh]["POSALTS"]])
        .sort("barcode")
    )

    if cluster == "bulk":
        # For bulk, we need to sum the values across all barcodes
        # For bulk, we need to aggregate across all barcodes
        df_alt_cov = (
            df_alt_cov.select(
                [
                    pl.sum(str(pos_alt)).alias(pos_alt)
                    for pos, pos_alt in zip(
                        CONFIG[hh]["POSISTIONS"], CONFIG[hh]["POSALTS"]
                    )
                ]
            )
            .with_columns(pl.lit("bulk").alias("barcode"))
            .select(["barcode", *CONFIG[hh]["POSALTS"]])
        )

    df_sum_cov = (
        df.group_by("barcode")
        .agg(
            [
                pl.sum(str(pos)).alias(str(pos))
                for pos in CONFIG[hh]["POSISTIONS"]
            ]
        )
        .sort("barcode")
    )

    df_sum_cov.columns = [
        "barcode",
        *[str(pos) for pos in CONFIG[hh]["POSALTS"]],
    ]

    if cluster == "bulk":
        # For bulk, we need to sum the values across all barcodes
        # For bulk, we need to aggregate across all barcodes
        df_sum_cov = (
            df_sum_cov.select(
                [
                    pl.sum(str(pos_alt)).alias(pos_alt)
                    for pos, pos_alt in zip(
                        CONFIG[hh]["POSISTIONS"], CONFIG[hh]["POSALTS"]
                    )
                ]
            )
            .with_columns(pl.lit("bulk").alias("barcode"))
            .select(["barcode", *CONFIG[hh]["POSALTS"]])
        )

    df_alt_af = df_alt_cov.select(CONFIG[hh]["POSALTS"]) / df_sum_cov.select(
        CONFIG[hh]["POSALTS"]
    )
    df_alt_af.columns = CONFIG[hh]["VARIANTS"]

    # Column bind barcode from df_sum_cov with df_alt_af
    result_df = pl.concat(
        [df_sum_cov.select(["barcode"]), df_alt_af], how="horizontal"
    )

    df_af = (
        result_df.with_columns(
            [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
        )
        .select(["gseid", "srrid", "barcode", *CONFIG[hh]["VARIANTS"]])
        .sort("barcode")
    )

    df_sum_depth = (
        df_sum_cov.with_columns(
            [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
        )
        .select(["gseid", "srrid", "barcode", *CONFIG[hh]["POSALTS"]])
        .sort("barcode")
    )
    df_sum_depth.columns = [
        "gseid",
        "srrid",
        "barcode",
        *CONFIG[hh]["VARIANTS"],
    ]
    df_alt_depth = (
        df_alt_cov.with_columns(
            [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
        )
        .select(["gseid", "srrid", "barcode", *CONFIG[hh]["POSALTS"]])
        .sort("barcode")
    )
    df_alt_depth.columns = [
        "gseid",
        "srrid",
        "barcode",
        *CONFIG[hh]["VARIANTS"],
    ]

    return df_af, df_sum_depth, df_alt_depth


def save_table_pl(gseid, srrid, cluster, hh):
    tablename = f"{cluster}_cov.{gseid}_{srrid}"
    # tablename = (
    #     f"{cluster}_cov.{gseid}_{srrid}"
    #     if cluster != "bulk"
    #     else f"cluster_cov.{gseid}_{srrid}"
    # )
    log.info(f"Loading {tablename}")
    df_af, df_sum_depth, df_alt_depth = load_table_pl(gseid, srrid, cluster, hh)
    table_path_af = (
        CONFIG[hh]["table_path"] / f"{tablename}.{CONFIG[hh]['suffix']}.csv"
    )
    df_af.write_csv(table_path_af, include_header=True)
    log.info(f"Saved {tablename} to {table_path_af}")
    table_path_sum_depth = (
        CONFIG[hh]["table_path"]
        / f"{tablename}.sumdepth.{CONFIG[hh]['suffix']}.csv"
    )
    df_sum_depth.write_csv(table_path_sum_depth, include_header=True)
    log.info(f"Saved {tablename} sum depth to {table_path_sum_depth}")
    table_path_alt_depth = (
        CONFIG[hh]["table_path"]
        / f"{tablename}.altdepth.{CONFIG[hh]['suffix']}.csv"
    )
    df_alt_depth.write_csv(table_path_alt_depth, include_header=True)
    log.info(f"Saved {tablename} alt depth to {table_path_alt_depth}")


def load_file(row, cluster, hh, filename):
    """Load a single heteroplasmic file.

    This function must be outside of merge_pl to be picklable for multiprocessing.
    """

    file_path = CONFIG[hh]["table_path"] / filename
    log.info(f"Loading {file_path}")
    try:
        df = pl.read_csv(file_path, has_header=True)

        df = df.with_columns(
            [pl.col(variant).fill_nan(0) for variant in CONFIG[hh]["VARIANTS"]]
        )
        # Create a binary mask where variant value > 0
        mask = df[CONFIG[hh]["VARIANTS"]] > 0

        # Calculate the row sum (number of variants > 0 per cell)
        df = df.with_columns(pl.sum_horizontal(mask).alias("num_variants"))
        log.info(
            f"Loaded {file_path} with {df.shape[0]} rows and {df.shape[1]} columns"
        )

        return df
    except Exception as e:
        log.error(f"Error loading {file_path}: {e}")
        return None


def merge_pl(cluster, hh):
    # Create a pool with multiple workers

    # Convert the dataframe to a list of dictionaries for multiprocessing
    row_list = SRR.to_dicts()

    # Process the rows in parallel
    dfs = []
    dfs_altdepth = []
    dfs_sumdepth = []
    for row in row_list:
        gseid = row["gseid"]
        srrid = row["srrid"]
        filename = f"{cluster}_cov.{gseid}_{srrid}.{CONFIG[hh]['suffix']}.csv"
        df = load_file(row, cluster, hh, filename)
        if df is not None:
            dfs.append(df)

        filename_altdepth = (
            f"{cluster}_cov.{gseid}_{srrid}.altdepth.{CONFIG[hh]['suffix']}.csv"
        )
        df_altdepth = load_file(row, cluster, hh, filename_altdepth)
        if df_altdepth is not None:
            dfs_altdepth.append(df_altdepth)

        filename_sumdepth = (
            f"{cluster}_cov.{gseid}_{srrid}.sumdepth.{CONFIG[hh]['suffix']}.csv"
        )
        df_sumdepth = load_file(row, cluster, hh, filename_sumdepth)
        if df_sumdepth is not None:
            dfs_sumdepth.append(df_sumdepth)

    # Remove any None values that may have occurred due to errors
    dfs = [df for df in dfs if df is not None]
    dfs_altdepth = [df for df in dfs_altdepth if df is not None]
    dfs_sumdepth = [df for df in dfs_sumdepth if df is not None]

    # Concatenate all dataframes
    all_df = pl.concat(dfs, how="vertical")

    # Save the merged result
    output_path = Path(
        f"/home/liuc9/github/scMOCHA-data/data/zzz/clean-data/all_{CONFIG[hh]['suffix']}_af.{cluster}.csv"
    )
    all_df.write_csv(output_path, include_header=True)
    log.info(f"Saved merged {CONFIG[hh]['suffix']} data to {output_path}")

    all_def_altdepth = pl.concat(dfs_altdepth, how="vertical")
    output_path_altdepth = Path(
        f"/home/liuc9/github/scMOCHA-data/data/zzz/clean-data/all_{CONFIG[hh]['suffix']}_altdepth.{cluster}.csv"
    )
    all_def_altdepth.write_csv(output_path_altdepth, include_header=True)
    log.info(
        f"Saved merged {CONFIG[hh]['suffix']} alt depth data to {output_path_altdepth}"
    )

    all_def_sumdepth = pl.concat(dfs_sumdepth, how="vertical")
    output_path_sumdepth = Path(
        f"/home/liuc9/github/scMOCHA-data/data/zzz/clean-data/all_{CONFIG[hh]['suffix']}_sumdepth.{cluster}.csv"
    )
    all_def_sumdepth.write_csv(output_path_sumdepth, include_header=True)
    log.info(
        f"Saved merged {CONFIG[hh]['suffix']} sum depth data to {output_path_sumdepth}"
    )


app = typer.Typer()

ClusterType = Annotated[
    str, typer.Argument(help="Type of cluster: cell or cluster or bulk")
]
HHType = Annotated[
    str, typer.Argument(help="Type of variant: HETEROPLASMIC or HOMOPLASMIC")
]


@app.command()
def create_sh(
    cluster: ClusterType,
    hh: HHType,
):
    """Create shell script for processing."""
    if cluster not in ["cell", "cluster", "bulk"]:
        typer.echo(f"Error: cluster must be 'cell' or 'cluster', got {cluster}")
        raise typer.Abort()

    if hh not in ["HETEROPLASMIC", "HOMOPLASMIC"]:
        typer.echo(
            f"Error: hh must be 'HETEROPLASMIC' or 'HOMOPLASMIC', got {hh}"
        )
        raise typer.Abort()

    typer.echo(f"Creating shell script for {cluster} {hh}")
    cmds = []
    for row in SRR.iter_rows(named=True):
        gseid = row["gseid"]
        srrid = row["srrid"]
        cmd = [
            "/scr1/users/liuc9/tools/anaconda3/envs/renv/bin/python3.13",
            "/home/liuc9/github/scMOCHA-data/analysis/ALL_VARIANT_CELLS_AF.py",
            "heteroplasmic-af",
            gseid,
            srrid,
            cluster,
            hh,
        ]
        cmds.append(" ".join(cmd))
        output_path = Path(
            f"/home/liuc9/github/scMOCHA-data/data/zzz/db/all_variant_cells.{cluster}.{CONFIG[hh]['suffix']}.run_all.sh"
        )
    with open(
        output_path,
        "w",
    ) as f:
        log.info(f"Writing commands to {output_path}")
        for cmd in cmds:
            f.write(f"{cmd}\n")


@app.command()
def create_all_sh():
    CLUSTER = ["cell", "cluster", "bulk"]
    HH = ["HETEROPLASMIC", "HOMOPLASMIC"]
    for cluster in CLUSTER:
        for hh in HH:
            create_sh(cluster, hh)


@app.command()
def heteroplasmic_af(
    gseid: str,
    srrid: str,
    cluster: ClusterType,
    hh: HHType,
):
    """
    Generate heteroplasmic allele frequency table for all cells.
    """
    if cluster not in ["cell", "cluster", "bulk"]:
        typer.echo(f"Error: cluster must be 'cell' or 'cluster', got {cluster}")
        raise typer.Abort()

    if hh not in ["HETEROPLASMIC", "HOMOPLASMIC"]:
        typer.echo(
            f"Error: hh must be 'HETEROPLASMIC' or 'HOMOPLASMIC', got {hh}"
        )
        raise typer.Abort()

    typer.echo(
        f"Generating heteroplasmic allele frequency table for {gseid} {srrid} {cluster} {hh}"
    )
    log.info("Generating heteroplasmic allele frequency table for all cells.")
    log.info(f"Processing {gseid} {srrid} {cluster} {hh}")
    save_table_pl(gseid, srrid, cluster, hh)
    log.info("Finished generating heteroplasmic allele frequency table.")


@app.command()
def merge_csv(
    cluster: ClusterType,
    hh: HHType,
):
    """
    Merge all heteroplasmic allele frequency tables.
    """
    if cluster not in ["cell", "cluster", "bulk"]:
        typer.echo(f"Error: cluster must be 'cell' or 'cluster', got {cluster}")
        raise typer.Abort()

    if hh not in ["HETEROPLASMIC", "HOMOPLASMIC"]:
        typer.echo(
            f"Error: hh must be 'HETEROPLASMIC' or 'HOMOPLASMIC', got {hh}"
        )
        raise typer.Abort()

    typer.echo(
        f"Merging all heteroplasmic allele frequency tables for {cluster} {hh}"
    )
    log.info("Merging all heteroplasmic allele frequency tables.")
    merge_pl(cluster, hh)
    log.info("Finished merging all heteroplasmic allele frequency tables.")


if __name__ == "__main__":
    # Use spawn method for better compatibility across platforms
    app()
