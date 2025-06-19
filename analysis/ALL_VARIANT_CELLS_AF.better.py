#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu (Improved by Assistant)
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-19 (Improved version)
# @DESCRIPTION: Improved version with better parallelization and optimization
# @VERSION: v1.0.0

import gc
import logging
import multiprocessing as mp
import os
from concurrent.futures import (
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Annotated, Any, Dict, List, Optional, Tuple

import polars as pl
import typer
from rich import print
from rich.logging import RichHandler
from rich.progress import Progress, TaskID

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

log = logging.getLogger("rich")


# Configuration class for better organization
@dataclass
class ProcessingConfig:
    """Configuration for processing parameters"""

    positions: List[int]
    variants: List[str]
    alts: List[str]
    posalts: List[str]
    table_path: Path
    suffix: str


# Load configuration data
SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)

ALL_VARIANT_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_variant.csv"
)
ALL_VARIANT_DF = pl.read_csv(ALL_VARIANT_FILENAME)

# Optimize dataframe filtering with lazy evaluation
HETEROPLASMIC_DF = (
    ALL_VARIANT_DF.lazy()
    .filter(pl.col("issomatic") == "heteroplasmic")
    .select(["variant", "Position"])
    .sort("Position")
    .collect()
)

HOMOPLASMIC_DF = (
    ALL_VARIANT_DF.lazy()
    .filter(pl.col("issomatic") == "homoplasmic")
    .select(["variant", "Position"])
    .sort("Position")
    .collect()
)


# Create configuration using dataclass for better type safety
def create_config(
    df: pl.DataFrame, suffix: str, table_path: str
) -> ProcessingConfig:
    """Create processing configuration from dataframe"""
    variants = df["variant"].to_list()
    positions = df["Position"].to_list()
    alts = [variant.split(">")[1] for variant in variants]
    posalts = [f"{pos}_{alt}" for pos, alt in zip(positions, alts)]

    return ProcessingConfig(
        positions=positions,
        variants=variants,
        alts=alts,
        posalts=posalts,
        table_path=Path(table_path),
        suffix=suffix,
    )


CONFIG = {
    "HETEROPLASMIC": create_config(
        HETEROPLASMIC_DF,
        "hetero",
        "/home/liuc9/github/scMOCHA-data/data/zzz/db/HETEROPLASMIC",
    ),
    "HOMOPLASMIC": create_config(
        HOMOPLASMIC_DF,
        "homo",
        "/home/liuc9/github/scMOCHA-data/data/zzz/db/HOMOPLASMIC",
    ),
}

TABLEDIR = Path("/home/liuc9/github/scMOCHA-data/analysis/zzz/db/TABLES")

# Optimized CPU count detection
N_CORES = min(
    mp.cpu_count(), 16
)  # Cap at 16 cores to avoid overwhelming system
CHUNK_SIZE = 100  # Process files in chunks to manage memory


def load_table_pl_optimized(
    gseid: str, srrid: str, cluster: str, hh: str
) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    # gseid, srrid, cluster, hh = "GSE147794", "GSM4446059", "cell", "HETEROPLASMIC"
    """Optimized table loading with lazy evaluation and memory management"""
    config = CONFIG[hh]

    tablename = f"{cluster}_cov.{gseid}_{srrid}"
    table_path = TABLEDIR / f"{tablename}.csv"

    if not table_path.exists():
        log.warning(f"File not found: {table_path}")
        return None, None, None

    try:
        # Use lazy loading and streaming for large files
        df = (
            pl.scan_csv(table_path, has_header=True)
            .select(
                ["base", "barcode", *[str(pos) for pos in config.positions]]
            )
            .collect(engine="streaming")
        )

        # Pivot operation (must be done on DataFrame, not LazyFrame)
        df_alt_cov = (
            df.pivot(
                index="barcode",
                on="base",
                values=[str(pos) for pos in config.positions],
                aggregate_function="sum",
            )
            .select(["barcode", *config.posalts])
            .sort("barcode")
        )

        if cluster == "bulk":
            # Optimize bulk processing
            df_alt_cov = (
                df_alt_cov.lazy()
                .select(
                    [
                        pl.sum(str(pos_alt)).alias(pos_alt)
                        for pos, pos_alt in zip(
                            config.positions, config.posalts
                        )
                    ]
                )
                .with_columns(pl.lit("bulk").alias("barcode"))
                .select(["barcode", *config.posalts])
                .collect()
            )

        # Optimize sum coverage calculation
        df_sum_cov = (
            df.lazy()
            .group_by("barcode")
            .agg([pl.sum(str(pos)).alias(str(pos)) for pos in config.positions])
            .sort("barcode")
            .collect()
        )

        df_sum_cov.columns = ["barcode", *config.posalts]

        if cluster == "bulk":
            df_sum_cov = (
                df_sum_cov.lazy()
                .select(
                    [
                        pl.sum(str(pos_alt)).alias(pos_alt)
                        for pos, pos_alt in zip(
                            config.positions, config.posalts
                        )
                    ]
                )
                .with_columns(pl.lit("bulk").alias("barcode"))
                .select(["barcode", *config.posalts])
                .collect()
            )

        # Calculate allele frequency with null handling
        df_alt_af = (
            df_alt_cov.select(config.posalts).fill_null(0)
            / df_sum_cov.select(config.posalts).fill_null(
                1
            )  # Avoid division by zero
        ).fill_null(0)
        df_alt_af.columns = config.variants

        # Combine results efficiently
        result_df = pl.concat(
            [df_sum_cov.select(["barcode"]), df_alt_af], how="horizontal"
        )

        df_af = (
            result_df.lazy()
            .with_columns(
                [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
            )
            .select(["gseid", "srrid", "barcode", *config.variants])
            .sort("barcode")
            .collect()
        )

        # Prepare depth dataframes
        df_sum_depth = (
            df_sum_cov.lazy()
            .with_columns(
                [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
            )
            .select(["gseid", "srrid", "barcode", *config.posalts])
            .sort("barcode")
            .collect()
        )
        df_sum_depth.columns = ["gseid", "srrid", "barcode", *config.variants]

        df_alt_depth = (
            df_alt_cov.lazy()
            .with_columns(
                [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
            )
            .select(["gseid", "srrid", "barcode", *config.posalts])
            .sort("barcode")
            .collect()
        )
        df_alt_depth.columns = ["gseid", "srrid", "barcode", *config.variants]

        # Clean up intermediate dataframes
        del df, df_alt_cov, df_sum_cov, df_alt_af, result_df
        gc.collect()

        return df_af, df_sum_depth, df_alt_depth

    except Exception as e:
        log.error(f"Error processing {tablename}: {e}")
        return None, None, None


def save_table_pl_parallel(args: Tuple[str, str, str, str]) -> bool:
    """Wrapper for parallel processing of save_table_pl"""
    gseid, srrid, cluster, hh = args
    return save_table_pl_optimized(gseid, srrid, cluster, hh)


def save_table_pl_optimized(
    gseid: str, srrid: str, cluster: str, hh: str
) -> bool:
    """Optimized table saving with error handling"""
    config = CONFIG[hh]
    tablename = f"{cluster}_cov.{gseid}_{srrid}"

    try:
        df_af, df_sum_depth, df_alt_depth = load_table_pl_optimized(
            gseid, srrid, cluster, hh
        )

        if df_af is None:
            return False

        # Ensure output directory exists
        config.table_path.mkdir(parents=True, exist_ok=True)

        # Save files with better error handling
        table_path_af = config.table_path / f"{tablename}.{config.suffix}.csv"
        df_af.write_csv(table_path_af, include_header=True)

        table_path_sum_depth = (
            config.table_path / f"{tablename}.sumdepth.{config.suffix}.csv"
        )
        df_sum_depth.write_csv(table_path_sum_depth, include_header=True)

        table_path_alt_depth = (
            config.table_path / f"{tablename}.altdepth.{config.suffix}.csv"
        )
        df_alt_depth.write_csv(table_path_alt_depth, include_header=True)

        log.info(f"Successfully processed {tablename}")
        return True

    except Exception as e:
        log.error(f"Error saving {tablename}: {e}")
        return False


def load_file_optimized(
    args: Tuple[Dict[str, str], str, str, str],
) -> Optional[pl.DataFrame]:
    """Optimized file loading with better error handling and memory management"""
    row, cluster, hh, filename = args
    config = CONFIG[hh]
    file_path = config.table_path / filename

    if not file_path.exists():
        log.warning(f"File not found: {file_path}")
        return None

    try:
        # Use lazy loading for better memory management
        df = (
            pl.scan_csv(file_path, has_header=True)
            .with_columns(
                [
                    pl.col(variant).fill_null(0).fill_nan(0)
                    for variant in config.variants
                ]
            )
            .collect(engine="streaming")
        )

        # Create binary mask for variants > 0 (more efficient)
        mask_cols = [
            (pl.col(variant) > 0).cast(pl.Int32) for variant in config.variants
        ]
        df = df.with_columns(
            pl.sum_horizontal(*mask_cols).alias("num_variants")
        )

        log.debug(
            f"Loaded {file_path} with {df.shape[0]} rows and {df.shape[1]} columns"
        )
        return df

    except Exception as e:
        log.error(f"Error loading {file_path}: {e}")
        return None


def merge_pl_parallel(cluster: str, hh: str) -> None:
    """Highly optimized parallel merging with progress tracking and memory management"""
    config = CONFIG[hh]
    row_list = SRR.to_dicts()

    log.info(
        f"Starting parallel merge for {cluster} {hh} with {len(row_list)} files"
    )

    # Process different file types in parallel
    file_types = [
        (f"{cluster}_cov", f".{config.suffix}.csv"),
        (f"{cluster}_cov", f".altdepth.{config.suffix}.csv"),
        (f"{cluster}_cov", f".sumdepth.{config.suffix}.csv"),
    ]

    all_results = {}

    with Progress() as progress:
        for file_prefix, file_suffix in file_types:
            task = progress.add_task(
                f"Processing {file_suffix} files...", total=len(row_list)
            )

            # Prepare arguments for parallel processing
            args_list = []
            for row in row_list:
                gseid = row["gseid"]
                srrid = row["srrid"]
                filename = f"{file_prefix}.{gseid}_{srrid}{file_suffix}"
                args_list.append((row, cluster, hh, filename))

            # Process files in parallel with progress tracking
            dfs = []
            with ProcessPoolExecutor(max_workers=N_CORES) as executor:
                # Submit all tasks
                future_to_args = {
                    executor.submit(load_file_optimized, args): args
                    for args in args_list
                }

                # Collect results as they complete
                for future in as_completed(future_to_args):
                    try:
                        df = future.result()
                        if df is not None:
                            dfs.append(df)
                        progress.advance(task)
                    except Exception as e:
                        args = future_to_args[future]
                        log.error(f"Error processing {args[3]}: {e}")
                        progress.advance(task)

            # Store results for this file type
            all_results[file_suffix] = dfs
            log.info(f"Collected {len(dfs)} valid dataframes for {file_suffix}")

    # Concatenate and save results for each file type
    output_configs = [
        (
            f".{config.suffix}.csv",
            "af",
            all_results.get(f".{config.suffix}.csv", []),
        ),
        (
            f".altdepth.{config.suffix}.csv",
            "altdepth",
            all_results.get(f".altdepth.{config.suffix}.csv", []),
        ),
        (
            f".sumdepth.{config.suffix}.csv",
            "sumdepth",
            all_results.get(f".sumdepth.{config.suffix}.csv", []),
        ),
    ]

    for suffix, output_type, dfs in output_configs:
        if not dfs:
            log.warning(f"No data found for {suffix}")
            continue

        try:
            # Use lazy concatenation for better memory management
            log.info(f"Concatenating {len(dfs)} dataframes for {output_type}")
            all_df = pl.concat(dfs, how="vertical")

            # Prepare output path
            output_dir = Path(
                "/home/liuc9/github/scMOCHA-data/data/zzz/clean-data"
            )
            output_dir.mkdir(parents=True, exist_ok=True)

            if output_type == "af":
                output_path = (
                    output_dir / f"all_{config.suffix}_af.{cluster}.csv"
                )
            else:
                output_path = (
                    output_dir
                    / f"all_{config.suffix}_{output_type}.{cluster}.csv"
                )

            # Save with streaming for large files
            all_df.write_csv(output_path, include_header=True)
            log.info(
                f"Saved merged {config.suffix} {output_type} data to {output_path}"
            )

            # Clean up memory
            del all_df, dfs
            gc.collect()

        except Exception as e:
            log.error(f"Error saving merged {output_type} data: {e}")


def process_batch_parallel(
    args_batch: List[Tuple[str, str, str, str]],
) -> List[bool]:
    """Process a batch of files in parallel"""
    with ProcessPoolExecutor(
        max_workers=min(N_CORES, len(args_batch))
    ) as executor:
        futures = [
            executor.submit(save_table_pl_parallel, args) for args in args_batch
        ]
        results = []
        for future in as_completed(futures):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                log.error(f"Batch processing error: {e}")
                results.append(False)
    return results


# CLI Application
app = typer.Typer(help="Optimized variant allele frequency analysis tool")

ClusterType = Annotated[
    str, typer.Argument(help="Type of cluster: cell, cluster, or bulk")
]
HHType = Annotated[
    str, typer.Argument(help="Type of variant: HETEROPLASMIC or HOMOPLASMIC")
]


@app.command()
def create_sh(cluster: ClusterType, hh: HHType):
    """Create optimized shell script for processing with better organization"""
    if cluster not in ["cell", "cluster", "bulk"]:
        typer.echo(
            f"Error: cluster must be 'cell', 'cluster', or 'bulk', got {cluster}"
        )
        raise typer.Abort()

    if hh not in ["HETEROPLASMIC", "HOMOPLASMIC"]:
        typer.echo(
            f"Error: hh must be 'HETEROPLASMIC' or 'HOMOPLASMIC', got {hh}"
        )
        raise typer.Abort()

    config = CONFIG[hh]
    output_dir = Path("/home/liuc9/github/scMOCHA-data/data/zzz/db")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = (
        output_dir
        / f"all_variant_cells.{cluster}.{config.suffix}.run_all_optimized.sh"
    )

    with open(output_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Optimized parallel processing script\n")
        f.write(f"# Generated for {cluster} {hh} processing\n\n")

        # Add error handling and logging
        f.write("set -euo pipefail\n")
        f.write('SCRIPT_DIR=$(dirname "$0")\n')
        f.write('LOG_FILE="${SCRIPT_DIR}/processing_${cluster}_${hh}.log"\n\n')

        # Process in batches for better resource management
        batch_size = 10
        rows = SRR.to_dicts()

        for i in range(0, len(rows), batch_size):
            batch = rows[i : i + batch_size]
            f.write(f"# Batch {i // batch_size + 1}\n")

            for row in batch:
                gseid = row["gseid"]
                srrid = row["srrid"]
                cmd = [
                    "/scr1/users/liuc9/tools/anaconda3/envs/renv/bin/python3.13",
                    "/home/liuc9/github/scMOCHA-data/analysis/ALL_VARIANT_CELLS_AF.better.py",
                    "heteroplasmic-af",
                    gseid,
                    srrid,
                    cluster,
                    hh,
                    ">>",
                    "${LOG_FILE}",
                    "2>&1",
                    "&",
                ]
                f.write(" ".join(cmd) + "\n")

            f.write("wait  # Wait for batch to complete\n\n")

    # Make script executable
    os.chmod(output_path, 0o755)
    log.info(f"Created optimized shell script: {output_path}")


@app.command()
def create_all_sh():
    """Create all shell scripts with parallel generation"""
    clusters = ["cell", "cluster", "bulk"]
    hh_types = ["HETEROPLASMIC", "HOMOPLASMIC"]

    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = []
        for cluster in clusters:
            for hh in hh_types:
                future = executor.submit(create_sh, cluster, hh)
                futures.append(future)

        # Wait for all scripts to be created
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error(f"Error creating shell script: {e}")


@app.command()
def heteroplasmic_af(gseid: str, srrid: str, cluster: ClusterType, hh: HHType):
    """Generate optimized allele frequency table for specific sample"""
    if cluster not in ["cell", "cluster", "bulk"]:
        typer.echo(
            f"Error: cluster must be 'cell', 'cluster', or 'bulk', got {cluster}"
        )
        raise typer.Abort()

    if hh not in ["HETEROPLASMIC", "HOMOPLASMIC"]:
        typer.echo(
            f"Error: hh must be 'HETEROPLASMIC' or 'HOMOPLASMIC', got {hh}"
        )
        raise typer.Abort()

    log.info(
        f"Processing {gseid} {srrid} {cluster} {hh} with optimized pipeline"
    )
    success = save_table_pl_optimized(gseid, srrid, cluster, hh)

    if success:
        log.info("Successfully completed processing")
    else:
        log.error("Processing failed")
        raise typer.Exit(1)


@app.command()
def merge_csv(cluster: ClusterType, hh: HHType):
    """Merge all files with highly optimized parallel processing"""
    if cluster not in ["cell", "cluster", "bulk"]:
        typer.echo(
            f"Error: cluster must be 'cell', 'cluster', or 'bulk', got {cluster}"
        )
        raise typer.Abort()

    if hh not in ["HETEROPLASMIC", "HOMOPLASMIC"]:
        typer.echo(
            f"Error: hh must be 'HETEROPLASMIC' or 'HOMOPLASMIC', got {hh}"
        )
        raise typer.Abort()

    log.info(f"Starting optimized merge for {cluster} {hh}")
    merge_pl_parallel(cluster, hh)
    log.info("Merge completed successfully")


@app.command()
def batch_process(cluster: ClusterType, hh: HHType, batch_size: int = 50):
    """Process all samples in optimized batches"""
    if cluster not in ["cell", "cluster", "bulk"]:
        typer.echo(
            f"Error: cluster must be 'cell', 'cluster', or 'bulk', got {cluster}"
        )
        raise typer.Abort()

    if hh not in ["HETEROPLASMIC", "HOMOPLASMIC"]:
        typer.echo(
            f"Error: hh must be 'HETEROPLASMIC' or 'HOMOPLASMIC', got {hh}"
        )
        raise typer.Abort()

    rows = SRR.to_dicts()
    total_samples = len(rows)

    log.info(f"Processing {total_samples} samples in batches of {batch_size}")

    # Process in batches
    successful = 0
    failed = 0

    for i in range(0, total_samples, batch_size):
        batch = rows[i : i + batch_size]
        log.info(
            f"Processing batch {i // batch_size + 1}/{(total_samples + batch_size - 1) // batch_size}"
        )

        # Prepare batch arguments
        args_batch = [
            (row["gseid"], row["srrid"], cluster, hh) for row in batch
        ]

        # Process batch in parallel
        results = process_batch_parallel(args_batch)

        batch_successful = sum(results)
        batch_failed = len(results) - batch_successful

        successful += batch_successful
        failed += batch_failed

        log.info(
            f"Batch completed: {batch_successful} successful, {batch_failed} failed"
        )

    log.info(
        f"Total processing completed: {successful} successful, {failed} failed"
    )


if __name__ == "__main__":
    # Set start method for multiprocessing
    mp.set_start_method("spawn", force=True)
    app()
