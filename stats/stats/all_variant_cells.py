#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-24 18:27:13
# @DESCRIPTION:
# @VERSION: v0.0.1

import logging
import multiprocessing as mp
from functools import partial
from pathlib import Path

import polars as pl
import typer
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)
SRR


ALL_VARIANT_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/all_variant.csv"
)
ALL_VARIANT_DF = pl.read_csv(ALL_VARIANT_FILENAME)
HETEROPLASMIC_DF = (
    ALL_VARIANT_DF.filter(pl.col("issomatic") == "heteroplasmic")
    .select(["variant", "Position"])
    .sort("Position")
)
# Extract positions from HETEROPLASMIC_DF
POSISTIONS = HETEROPLASMIC_DF["Position"].to_list()
VARIANTS = HETEROPLASMIC_DF["variant"].to_list()
ALTS = [variant.split(">")[1] for variant in VARIANTS]
POSALTS = [f"{pos}_{alt}" for pos, alt in zip(POSISTIONS, ALTS)]

TABLEDIR = Path("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/db/TABLES")


def load_table_pl(gseid, srrid):
    # gseid, srrid = "GSE147794", "GSM4446059"
    tablename = f"cell_cov.{gseid}_{srrid}"
    log.info(f"Loading {tablename}")
    table_path = TABLEDIR / f"{tablename}.csv"
    df = pl.read_csv(
        table_path,
        has_header=True,
    ).select(["base", "barcode", *[str(pos) for pos in POSISTIONS]])

    df_alt_cov = (
        df.pivot(
            index="barcode",
            on="base",
            values=[str(pos) for pos in POSISTIONS],
            aggregate_function="sum",
        )
        .select(["barcode", *POSALTS])
        .sort("barcode")
    )

    df_sum_cov = (
        df.group_by("barcode")
        .agg([pl.sum(str(pos)).alias(str(pos)) for pos in POSISTIONS])
        .sort("barcode")
    )

    df_sum_cov.columns = ["barcode", *[str(pos) for pos in POSALTS]]

    df_alt_af = df_alt_cov.select(POSALTS) / df_sum_cov.select(POSALTS)
    df_alt_af.columns = VARIANTS

    # Column bind barcode from df_sum_cov with df_alt_af
    result_df = pl.concat(
        [df_sum_cov.select(["barcode"]), df_alt_af], how="horizontal"
    )

    df_out = (
        result_df.with_columns(
            [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
        )
        .select(["gseid", "srrid", "barcode", *VARIANTS])
        .sort("barcode")
    )

    return df_out


def save_table_pl(gseid, srrid):
    tablename = f"cell_cov.{gseid}_{srrid}"
    log.info(f"Loading {tablename}")
    table_path = (
        Path("/home/liuc9/github/scMOCHA-data/data/zzz/db/HETEROPLASMIC")
        / f"{tablename}.hetero.csv"
    )
    df = load_table_pl(gseid, srrid)
    df.write_csv(table_path, include_header=True)
    log.info(f"Saved {tablename} to {table_path}")


def load_file(row):
    """Load a single heteroplasmic file.

    This function must be outside of merge_pl to be picklable for multiprocessing.
    """
    gseid = row["gseid"]
    srrid = row["srrid"]
    filename = f"cell_cov.{gseid}_{srrid}.hetero.csv"
    file_path = (
        Path("/home/liuc9/github/scMOCHA-data/data/zzz/db/HETEROPLASMIC")
        / filename
    )
    log.info(f"Loading {file_path}")
    try:
        df = pl.read_csv(file_path, has_header=True)

        df = df.with_columns(
            [pl.col(variant).fill_nan(0) for variant in VARIANTS]
        )
        # Create a binary mask where variant value > 0
        mask = df[VARIANTS] > 0

        # Calculate the row sum (number of variants > 0 per cell)
        df = df.with_columns(pl.sum_horizontal(mask).alias("num_variants"))
        log.info(
            f"Loaded {file_path} with {df.shape[0]} rows and {df.shape[1]} columns"
        )

        return df
    except Exception as e:
        log.error(f"Error loading {file_path}: {e}")
        return None


def merge_pl():
    # Create a pool with multiple workers

    # Convert the dataframe to a list of dictionaries for multiprocessing
    row_list = SRR.to_dicts()

    # Process the rows in parallel
    dfs = []
    for row in row_list:
        df = load_file(row)
        if df is not None:
            dfs.append(df)

    # Remove any None values that may have occurred due to errors
    dfs = [df for df in dfs if df is not None]

    # Concatenate all dataframes
    all_df = pl.concat(dfs, how="vertical")

    # Save the merged result
    output_path = Path(
        "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data/all_heteroplasmic_af.csv"
    )
    all_df.write_csv(output_path, include_header=True)
    log.info(f"Saved merged heteroplasmic data to {output_path}")


def create_sh():
    cmds = []
    for row in SRR.iter_rows(named=True):
        gseid = row["gseid"]
        srrid = row["srrid"]
        cmd = [
            "/scr1/users/liuc9/tools/anaconda3/envs/renv/bin/python3.13",
            "/home/liuc9/github/scMOCHA-data/stats/stats/all_variant_cells.py",
            gseid,
            srrid,
        ]
        cmds.append(" ".join(cmd))
    with open(
        "/home/liuc9/github/scMOCHA-data/data/zzz/db/all_variant_cells.run_all.sh",
        "w",
    ) as f:
        for cmd in cmds:
            f.write(f"{cmd}\n")


def process_row(row):
    gseid = row["gseid"]
    srrid = row["srrid"]
    log.info(f"Processing {gseid} {srrid}")
    return load_table_pl(gseid, srrid)


def parallel_run():
    # Create a pool with 80 CPU workers
    pool = mp.Pool(processes=20)

    # Convert the dataframe to a list of dictionaries for multiprocessing
    row_list = SRR.to_dicts()

    # Process the rows in parallel
    dfs = pool.map(process_row, row_list)

    # Clean up the pool
    pool.close()
    pool.join()

    all_df = pl.concat(dfs, how="vertical")
    all_df.write_csv(
        Path(
            "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/db/all_heteroplasmic_af.csv"
        ),
        has_header=True,
    )


app = typer.Typer()


@app.command()
def heteroplasmic_af(gseid: str, srrid: str):
    """
    Generate heteroplasmic allele frequency table for all cells.
    """
    log.info("Generating heteroplasmic allele frequency table for all cells.")
    save_table_pl(gseid, srrid)
    log.info("Finished generating heteroplasmic allele frequency table.")


@app.command()
def merge():
    """python /home/liuc9/github/scMOCHA-data/stats/stats/all_variant_cells.py
    Merge all heteroplasmic allele frequency tables.
    """
    log.info("Merging all heteroplasmic allele frequency tables.")
    merge_pl()
    log.info("Finished merging all heteroplasmic allele frequency tables.")


if __name__ == "__main__":
    # Use spawn method for better compatibility across platforms
    app()
