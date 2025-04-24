#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-23 23:58:00
# @DESCRIPTION:
# @VERSION: v0.0.1

import asyncio
import concurrent.futures
import functools
import logging
import multiprocessing as mp
import time
from pathlib import Path

import aiofiles
import duckdb
import pandas as pd
import polars as pl  # Add polars import
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")
log.info("Hello, World!")

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)
SRR

FASTA_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta.df.csv"
)
FASTA = pl.read_csv(FASTA_FILENAME)
FASTA
# Create a Series of positions for use in the join operation
POS = pl.DataFrame({"pos": FASTA["pos"]})


BASES = ["A", "C", "G", "T"]
CLUSTERS = ["cell", "cluster"]
TABLEDIR = Path("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/TABLES")


def load_pl_from_base(filename, base):
    log.info(f"Loading {filename} for base {base}")

    # Use polars to read the CSV file
    df = pl.read_csv(
        filename,
        has_header=False,
        new_columns=["pos", "barcode", "fw", "rv"],
        try_parse_dates=False,
    )

    # Left join POS with df to ensure all positions are represented
    full_df = POS.join(df, on="pos", how="left")

    # Fill NaN values with 0 for counts
    full_df = full_df.fill_null(0)

    # Convert to integer type first
    full_df = full_df.with_columns(
        [
            pl.col("fw").cast(pl.Int32),
            pl.col("rv").cast(pl.Int32),
        ]
    )

    # Now calculate coverage after ensuring fw and rv are properly typed
    full_df = full_df.with_columns(
        (pl.col("fw") + pl.col("rv")).alias("cov").cast(pl.Int32)
    )

    # Create fw pivot
    fw_pivot = (
        full_df.pivot(
            index="barcode",
            on="pos",
            values="fw",
            aggregate_function="first",
        )
        .fill_null(0)
        .with_columns([pl.lit(base).alias("base")])
        .filter(pl.col("barcode").is_not_null())
    )

    # Create rv pivot
    rv_pivot = (
        full_df.pivot(
            index="barcode",
            on="pos",
            values="rv",
            aggregate_function="first",
        )
        .fill_null(0)
        .with_columns([pl.lit(base).alias("base")])
        .filter(pl.col("barcode").is_not_null())
    )

    # Create cov pivot
    cov_pivot = (
        full_df.pivot(
            index="barcode",
            on="pos",
            values="cov",
            aggregate_function="first",
        )
        .fill_null(0)
        .with_columns([pl.lit(base).alias("base")])
        .filter(pl.col("barcode").is_not_null())
    )

    # Reorganize columns to match the expected output format
    base_col = pl.col("base")

    fw_pivot = fw_pivot.select(
        [base_col, *[col for col in fw_pivot.columns if col != "base"]]
    )
    rv_pivot = rv_pivot.select(
        [base_col, *[col for col in rv_pivot.columns if col != "base"]]
    )
    cov_pivot = cov_pivot.select(
        [base_col, *[col for col in cov_pivot.columns if col != "base"]]
    )

    log.info(
        f"Loaded {filename} for base {base}: fw_pivot: {fw_pivot.shape}, rv_pivot: {rv_pivot.shape}, cov_pivot: {cov_pivot.shape}"
    )
    # Convert back to pandas for compatibility with downstream code
    return fw_pivot, rv_pivot, cov_pivot


def load_pl(srrdir, cluster="cell"):
    filenames = [f"{srrdir}/{cluster}.{base}.txt.gz" for base in BASES]

    # Load data for each base in parallel using polars
    results = []
    for filename, base in zip(filenames, BASES):
        results.append(load_pl_from_base(filename, base))

    # Unpack results
    df_fw, df_rv, df_cov = zip(*results)
    # Use polars for faster concatenation
    df_fw_pl = pl.concat(df_fw)
    df_rv_pl = pl.concat(df_rv)
    df_cov_pl = pl.concat(df_cov)

    # Convert back to pandas for compatibility
    return df_fw_pl, df_rv_pl, df_cov_pl


def srrid2csv(gseid, srrid, srrdir, cluster="cell"):
    log.info(f"Processing {gseid}_{srrid} in {srrdir} for cluster {cluster}")
    # Load the data
    pl_fw, pl_rv, pl_cov = load_pl(srrdir, cluster)
    tablename = f"{gseid}_{srrid}"
    # f"{TABLEDIR}/{cluster}_fw.{tablename}.csv"
    pl_fw.write_csv(
        f"{TABLEDIR}/{cluster}_fw.{tablename}.csv", include_header=True
    )
    pl_rv.write_csv(
        f"{TABLEDIR}/{cluster}_rv.{tablename}.csv", include_header=True
    )
    pl_cov.write_csv(
        f"{TABLEDIR}/{cluster}_cov.{tablename}.csv", include_header=True
    )
    log.info(
        f"Saved {gseid}_{srrid} in {srrdir} for cluster {cluster} to {TABLEDIR}"
    )
    # Save the data to CSV files


def process_srr_entry(row):
    """Process a single SRR entry for both cell and cluster."""
    gseid = row["gseid"]
    srrid = row["srrid"]
    srrdir = row["srrdir"]
    cluster = row["cluster"]
    log.info(f"Processing {gseid}_{srrid} in {srrdir} for cluster {cluster}")
    srrid2csv(gseid, srrid, srrdir, cluster=cluster)
    return f"Processed {gseid}_{srrid}"


def all_srrid2csv(max_workers=mp.cpu_count()):
    """Process all SRR entries in parallel using ProcessPoolExecutor."""
    start_time = time.time()
    log.info(
        f"Starting parallel processing of {len(SRR)} SRR entries with {max_workers} workers"
    )

    # Convert DataFrame rows to dictionaries for easier processing
    row_dicts = [
        {**row, "cluster": cluster}
        for row in SRR.to_dicts()
        for cluster in CLUSTERS
    ]

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=max_workers
    ) as executor:
        # Map the process_srr_entry function to all rows
        results = list(executor.map(process_srr_entry, row_dicts))

    end_time = time.time()
    log.info(
        f"Completed processing {len(results)} SRR entries in {end_time - start_time:.2f} seconds"
    )
    return results


if __name__ == "__main__":
    # srrid2csv(
    #     "GSE155673",
    #     "SRR11512399",
    #     "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE155673/final/GSM4712885",
    #     cluster="cell",
    # )
    # load_df_from_base(
    #     filename="/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE155673/final/GSM4712885/cell.A.txt.gz",
    #     base="A",
    # )
    # load_df(
    #     srrdir="/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE155673/final/GSM4712885",
    #     cluster="cell",
    # )
    all_srrid2csv()
