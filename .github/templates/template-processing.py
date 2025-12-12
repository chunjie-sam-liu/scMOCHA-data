#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-12 10:00:00
# @DESCRIPTION: Template for Python processing script
# @VERSION: v0.0.1

# Standard library imports
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Annotated, Optional

# Third-party imports
import polars as pl
import typer
from rich import print
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

# Constants and configuration
BASEDIR = Path("/liulab/chunjie/data/scMOCHA")
DATADIR = BASEDIR / "data"
OUTDIR = BASEDIR / "output"
MAX_WORKERS = 20

CONFIG = {
    "param1": "value1",
    "param2": "value2",
}

# Logging setup
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)
logger = logging.getLogger("scMOCHA")


# Class definitions
class DataProcessor:
    """Process scMOCHA data with configurable parameters."""

    def __init__(self, config: dict):
        """Initialize processor with configuration."""
        self.config = config
        logger.info("Initialized DataProcessor")

    def process(self, data: pl.DataFrame) -> pl.DataFrame:
        """
        Process input data.

        Args:
            data: Input polars DataFrame

        Returns:
            Processed polars DataFrame
        """
        # Processing logic
        processed = data.filter(pl.col("value").is_not_null()).with_columns(
            (pl.col("value") * 2).alias("processed_value")
        )
        return processed


# Function definitions
def process_file(filepath: Path) -> Optional[pl.DataFrame]:
    """
    Process a single file.

    Args:
        filepath: Path to input file

    Returns:
        Processed DataFrame or None if processing fails
    """
    try:
        # Read data
        data = pl.read_csv(filepath)

        # Process
        processor = DataProcessor(CONFIG)
        result = processor.process(data)

        logger.info(f"Processed {filepath.name}: {result.shape[0]} rows")
        return result

    except Exception as e:
        logger.error(f"Failed to process {filepath}: {e}")
        return None


def process_sample(row: dict) -> dict:
    """
    Process a single sample (for parallel execution).

    Args:
        row: Dictionary containing sample metadata

    Returns:
        Dictionary with processing results
    """
    gseid = row["gseid"]
    srrid = row["srrid"]

    try:
        # Construct paths
        input_path = DATADIR / gseid / srrid / "input.csv"
        output_path = DATADIR / gseid / srrid / "output.parquet"

        # Check if input exists
        if not input_path.exists():
            logger.warning(f"Input file not found: {input_path}")
            return {
                "srrid": srrid,
                "status": "skipped",
                "reason": "missing_input",
            }

        # Process data
        data = pl.read_csv(input_path)
        processed = data.filter(pl.col("value").is_not_null())

        # Save result
        output_path.parent.mkdir(parents=True, exist_ok=True)
        processed.write_parquet(output_path)

        logger.info(f"Processed {srrid}: {processed.shape[0]} rows")
        return {"srrid": srrid, "status": "success", "nrow": processed.shape[0]}

    except Exception as e:
        logger.error(f"Failed to process {srrid}: {e}")
        return {"srrid": srrid, "status": "failed", "error": str(e)}


def process_batch(
    metadata_file: Path, max_workers: int = MAX_WORKERS
) -> pl.DataFrame:
    """
    Process batch of samples in parallel.

    Args:
        metadata_file: Path to metadata CSV file
        max_workers: Number of parallel workers

    Returns:
        DataFrame with processing results
    """
    # Read metadata
    metadata = pl.read_csv(metadata_file)
    logger.info(f"Loaded {metadata.shape[0]} samples from metadata")

    # Process in parallel with progress bar
    results = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task(
            "[cyan]Processing samples...", total=metadata.shape[0]
        )

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(process_sample, row)
                for row in metadata.iter_rows(named=True)
            ]

            for future in futures:
                result = future.result()
                results.append(result)
                progress.update(task, advance=1)

    # Convert results to DataFrame
    results_df = pl.DataFrame(results)

    # Log summary
    success_count = results_df.filter(pl.col("status") == "success").shape[0]
    failed_count = results_df.filter(pl.col("status") == "failed").shape[0]
    skipped_count = results_df.filter(pl.col("status") == "skipped").shape[0]

    logger.info("Processing summary:")
    logger.info(f"  Success: {success_count}")
    logger.info(f"  Failed: {failed_count}")
    logger.info(f"  Skipped: {skipped_count}")

    return results_df


# Main CLI function
def main(
    input_file: Annotated[str, typer.Argument(help="Input file path")],
    output_file: Annotated[
        str, typer.Option("--output", "-o", help="Output file path")
    ] = "output.parquet",
    gseid: Annotated[
        Optional[str], typer.Option("--gseid", "-g", help="GSE accession ID")
    ] = None,
    batch_mode: Annotated[
        bool,
        typer.Option("--batch", "-b", help="Process batch from metadata file"),
    ] = False,
    max_workers: Annotated[
        int, typer.Option("--workers", "-w", help="Number of parallel workers")
    ] = MAX_WORKERS,
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable verbose logging")
    ] = False,
):
    """
    Process scMOCHA data files.

    Examples:
        # Process single file
        python script.py input.csv -o output.parquet

        # Process batch in parallel
        python script.py metadata.csv --batch --workers 10

        # Verbose mode
        python script.py input.csv -v
    """
    # Set logging level
    if verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose mode enabled")

    # Convert to Path
    input_path = Path(input_file)
    output_path = Path(output_file)

    # Validate input file exists
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        raise typer.Exit(code=1)

    # Process data
    logger.info(f"Processing {input_path}...")

    if batch_mode:
        # Batch processing mode
        results = process_batch(input_path, max_workers=max_workers)

        # Save results
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results.write_parquet(output_path)
        logger.info(f"Saved batch results to {output_path}")

    else:
        # Single file processing mode
        result = process_file(input_path)

        if result is None:
            logger.error("Processing failed")
            raise typer.Exit(code=1)

        # Save result
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result.write_parquet(output_path)
        logger.info(f"Saved output to {output_path}")

    logger.info("✓ Processing complete!")


if __name__ == "__main__":
    typer.run(main)
