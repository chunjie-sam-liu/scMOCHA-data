#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-22 18:31:34
# @DESCRIPTION: Convert CSV files to Parquet format with deduplication and parallel processing
# @VERSION: v1.0.0

import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Literal, Optional, Set, Tuple

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

# Initialize typer app and console
app = typer.Typer(
    help="Convert CSV files to Parquet format with efficient processing"
)
console = Console()

# Default configurations
DEFAULT_SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
DEFAULT_TABLEPATH = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/TABLES"
)
DEFAULT_PARQUETPATH = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/PARQUET"
)
TABLENAMES = {"cov": "covall", "fw": "covfw", "rv": "covrv"}


# Setup logging
def setup_logging(log_level: str = "INFO") -> None:
    """Setup logging with rich handler."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=console)],
    )


class ParquetConverter:
    """Main class for handling CSV to Parquet conversion."""

    def __init__(
        self, table_path: Path, parquet_path: Path, srr_data: pl.DataFrame
    ):
        self.table_path = table_path
        self.parquet_path = parquet_path
        self.srr_data = srr_data

    def get_existing_parquet_files(self, root_path: Path) -> Set[str]:
        """Get set of existing parquet file identifiers to avoid reprocessing."""
        existing_files = set()
        if root_path.exists():
            for parquet_file in root_path.rglob("*.parquet"):
                # Extract gseid_srrid from path structure
                path_parts = parquet_file.parts
                # Find gseid and srrid from partition structure
                gseid = srrid = None
                for part in path_parts:
                    if part.startswith("gseid="):
                        gseid = part.split("=")[1]
                    elif part.startswith("srrid="):
                        srrid = part.split("=")[1]
                if gseid and srrid:
                    existing_files.add(f"{gseid}_{srrid}")
        return existing_files

    def process_one_sample_to_arrow(
        self,
        file_path: Path,
        gseid: str,
        srrid: str,
    ) -> Optional[pa.Table]:
        """Convert a single CSV file to Arrow table format."""
        try:
            if not file_path.exists():
                logging.warning(f"File not found: {file_path}")
                return None

            # Use streaming for better memory efficiency
            df = (
                pl.scan_csv(file_path, try_parse_dates=False)
                .with_columns(
                    [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
                )
                .collect(streaming=True)
            )
            return df.to_arrow()
        except Exception as e:
            logging.error(f"Error processing {file_path}: {e}")
            return None

    def process_one_sample(
        self,
        gseid: str,
        srrid: str,
        cluster: str,
        cov_fw_rv: Literal["cov", "fw", "rv"],
        root_path: Path,
    ) -> bool:
        """Process a single sample and write to parquet."""
        try:
            csv_filename = f"{cluster}_{cov_fw_rv}.{gseid}_{srrid}.csv"
            csv_path = self.table_path / csv_filename

            if not csv_path.exists():
                logging.warning(f"CSV file not found: {csv_path}")
                return False

            logging.debug(f"Processing {csv_path}")

            arrow_table = self.process_one_sample_to_arrow(
                file_path=csv_path,
                gseid=gseid,
                srrid=srrid,
            )

            if arrow_table is None:
                return False

            logging.debug(f"Writing {gseid}_{srrid} to Parquet")

            pq.write_to_dataset(
                arrow_table,
                root_path=root_path,
                partition_cols=["gseid", "srrid", "base"],
                compression="zstd",
                existing_data_behavior="overwrite_or_ignore",
            )
            return True
        except Exception as e:
            logging.error(f"Error processing sample {gseid}_{srrid}: {e}")
            return False

    def process_samples_for_type(
        self,
        cov_fw_rv: Literal["cov", "fw", "rv"],
        cluster_filter: Optional[str] = None,
        num_workers: int = 8,
        skip_existing: bool = True,
    ) -> Tuple[int, int]:
        """Process all samples for a specific coverage type."""
        root_path = self.parquet_path / TABLENAMES[cov_fw_rv]
        root_path.mkdir(parents=True, exist_ok=True)

        # Get existing files if skip_existing is True
        existing_files = (
            self.get_existing_parquet_files(root_path)
            if skip_existing
            else set()
        )

        # Prepare tasks
        tasks = []
        for row in self.srr_data.iter_rows(named=True):
            gseid = row["gseid"]
            srrid = row["srrid"]
            cluster = cluster_filter or "cell"  # Default cluster if not in data

            # Skip if already processed
            if skip_existing and f"{gseid}_{srrid}" in existing_files:
                continue

            tasks.append((gseid, srrid, cluster, cov_fw_rv, root_path))

        if not tasks:
            console.print(
                f"[green]No new files to process for {cov_fw_rv}[/green]"
            )
            return 0, 0

        processed = 0
        failed = 0

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=console,
            refresh_per_second=1,  # Reduce refresh rate for better performance
        ) as progress:
            task_id = progress.add_task(
                f"Processing {cov_fw_rv} files", total=len(tasks)
            )

            with ThreadPoolExecutor(max_workers=num_workers) as executor:
                future_to_task = {
                    executor.submit(
                        self.process_one_sample,
                        gseid,
                        srrid,
                        cluster,
                        cov_fw_rv,
                        root_path,
                    ): (gseid, srrid)
                    for gseid, srrid, cluster, cov_fw_rv, root_path in tasks
                }

                for future in as_completed(future_to_task):
                    gseid, srrid = future_to_task[future]
                    try:
                        success = future.result()
                        if success:
                            processed += 1
                            logging.info(f"✓ Processed {gseid}_{srrid}")
                        else:
                            failed += 1
                            logging.warning(f"✗ Failed {gseid}_{srrid}")
                    except Exception as e:
                        logging.error(f"Task failed for {gseid}_{srrid}: {e}")
                        failed += 1

                    progress.advance(task_id)

        return processed, failed

    def process_samples_for_type_optimized(
        self,
        cov_fw_rv: Literal["cov", "fw", "rv"],
        cluster_filter: Optional[str] = None,
        num_workers: int = 8,
        skip_existing: bool = True,
        batch_size: int = 100,
    ) -> Tuple[int, int]:
        """Process all samples for a specific coverage type with optimizations."""
        root_path = self.parquet_path / TABLENAMES[cov_fw_rv]
        root_path.mkdir(parents=True, exist_ok=True)

        # Get existing files if skip_existing is True
        existing_files = (
            self.get_existing_parquet_files(root_path)
            if skip_existing
            else set()
        )

        # Prepare tasks
        tasks = []
        for row in self.srr_data.iter_rows(named=True):
            gseid = row["gseid"]
            srrid = row["srrid"]
            cluster = cluster_filter or "cell"  # Default cluster if not in data

            # Skip if already processed
            if skip_existing and f"{gseid}_{srrid}" in existing_files:
                continue

            tasks.append((gseid, srrid, cluster, cov_fw_rv, root_path))

        if not tasks:
            console.print(
                f"[green]No new files to process for {cov_fw_rv}[/green]"
            )
            return 0, 0

        processed = 0
        failed = 0
        total_tasks = len(tasks)

        console.print(
            f"[blue]Processing {total_tasks} tasks for {cov_fw_rv}[/blue]"
        )

        # Process in batches for better progress tracking
        for batch_start in range(0, total_tasks, batch_size):
            batch_end = min(batch_start + batch_size, total_tasks)
            batch_tasks = tasks[batch_start:batch_end]
            batch_num = batch_start // batch_size + 1
            total_batches = (total_tasks + batch_size - 1) // batch_size

            console.print(
                f"[cyan]Batch {batch_num}/{total_batches}: Processing {len(batch_tasks)} files[/cyan]"
            )

            with ThreadPoolExecutor(max_workers=num_workers) as executor:
                future_to_task = {
                    executor.submit(
                        self.process_one_sample,
                        gseid,
                        srrid,
                        cluster,
                        cov_fw_rv,
                        root_path,
                    ): (gseid, srrid)
                    for gseid, srrid, cluster, cov_fw_rv, root_path in batch_tasks
                }

                batch_processed = 0
                batch_failed = 0

                for future in as_completed(future_to_task):
                    gseid, srrid = future_to_task[future]
                    try:
                        success = future.result()
                        if success:
                            batch_processed += 1
                            processed += 1
                        else:
                            batch_failed += 1
                            failed += 1
                    except Exception as e:
                        logging.error(f"Task failed for {gseid}_{srrid}: {e}")
                        batch_failed += 1
                        failed += 1

                console.print(
                    f"[green]Batch {batch_num} complete: {batch_processed} processed, {batch_failed} failed[/green]"
                )

        return processed, failed


@app.command()
def convert(
    cov_types: List[str] = typer.Option(
        ["cov", "fw", "rv"],
        "--type",
        "-t",
        help="Coverage types to process (cov, fw, rv)",
    ),
    cluster: Optional[str] = typer.Option(
        "cell", "--cluster", "-c", help="Cluster type (cell, cluster, bulk)"
    ),
    num_workers: int = typer.Option(
        8, "--workers", "-w", help="Number of worker threads"
    ),
    srr_file: Path = typer.Option(
        DEFAULT_SRR_FILENAME,
        "--srr-file",
        "-s",
        help="Path to SRR CSV file",
        exists=True,
    ),
    table_path: Path = typer.Option(
        DEFAULT_TABLEPATH,
        "--table-path",
        help="Path to input CSV tables directory",
        exists=True,
    ),
    parquet_path: Path = typer.Option(
        DEFAULT_PARQUETPATH,
        "--parquet-path",
        "-o",
        help="Path to output Parquet directory",
    ),
    skip_existing: bool = typer.Option(
        True,
        "--skip-existing/--no-skip-existing",
        help="Skip files that already exist in parquet format",
    ),
    log_level: str = typer.Option(
        "INFO",
        "--log-level",
        "-l",
        help="Logging level (DEBUG, INFO, WARNING, ERROR)",
    ),
) -> None:
    """Convert CSV files to Parquet format with parallel processing."""
    setup_logging(log_level)

    # Validate coverage types
    valid_types = {"cov", "fw", "rv"}
    invalid_types = set(cov_types) - valid_types
    if invalid_types:
        console.print(f"[red]Invalid coverage types: {invalid_types}[/red]")
        console.print(f"[yellow]Valid types: {valid_types}[/yellow]")
        raise typer.Exit(1)

    try:
        # Load SRR data
        console.print(f"[blue]Loading SRR data from {srr_file}[/blue]")
        srr_data = pl.read_csv(srr_file)
        console.print(f"[green]Loaded {len(srr_data)} samples[/green]")

        # Initialize converter
        converter = ParquetConverter(table_path, parquet_path, srr_data)

        # Process each coverage type
        total_processed = 0
        total_failed = 0

        for cov_type in cov_types:
            console.print(
                f"\n[blue]Processing coverage type: {cov_type}[/blue]"
            )
            processed, failed = converter.process_samples_for_type_optimized(
                cov_fw_rv=cov_type,
                cluster_filter=cluster,
                num_workers=num_workers,
                skip_existing=skip_existing,
            )
            total_processed += processed
            total_failed += failed

            console.print(
                f"[green]Completed {cov_type}: {processed} processed, {failed} failed[/green]"
            )

        # Summary
        console.print("\n[bold green]Total Summary:[/bold green]")
        console.print(
            f"[green]Successfully processed: {total_processed}[/green]"
        )
        console.print(f"[red]Failed: {total_failed}[/red]")

        if total_failed > 0:
            raise typer.Exit(1)

    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def convert_fast(
    cov_types: List[str] = typer.Option(
        ["cov", "fw", "rv"],
        "--type",
        "-t",
        help="Coverage types to process (cov, fw, rv)",
    ),
    cluster: Optional[str] = typer.Option(
        "cell", "--cluster", "-c", help="Cluster type (cell, cluster, bulk)"
    ),
    num_workers: int = typer.Option(
        8, "--workers", "-w", help="Number of worker threads"
    ),
    batch_size: int = typer.Option(
        50, "--batch-size", "-b", help="Number of files to process per batch"
    ),
    srr_file: Path = typer.Option(
        DEFAULT_SRR_FILENAME,
        "--srr-file",
        "-s",
        help="Path to SRR CSV file",
        exists=True,
    ),
    table_path: Path = typer.Option(
        DEFAULT_TABLEPATH,
        "--table-path",
        help="Path to input CSV tables directory",
        exists=True,
    ),
    parquet_path: Path = typer.Option(
        DEFAULT_PARQUETPATH,
        "--parquet-path",
        "-o",
        help="Path to output Parquet directory",
    ),
    skip_existing: bool = typer.Option(
        True,
        "--skip-existing/--no-skip-existing",
        help="Skip files that already exist in parquet format",
    ),
    log_level: str = typer.Option(
        "WARNING",
        "--log-level",
        "-l",
        help="Logging level (DEBUG, INFO, WARNING, ERROR)",
    ),
) -> None:
    """Convert CSV files to Parquet format with fast processing (reduced logging)."""
    setup_logging(log_level)

    # Validate coverage types
    valid_types = {"cov", "fw", "rv"}
    invalid_types = set(cov_types) - valid_types
    if invalid_types:
        console.print(f"[red]Invalid coverage types: {invalid_types}[/red]")
        console.print(f"[yellow]Valid types: {valid_types}[/yellow]")
        raise typer.Exit(1)

    try:
        # Load SRR data
        console.print(f"[blue]Loading SRR data from {srr_file}[/blue]")
        srr_data = pl.read_csv(srr_file)
        console.print(f"[green]Loaded {len(srr_data)} samples[/green]")

        # Initialize converter
        converter = ParquetConverter(table_path, parquet_path, srr_data)

        # Process each coverage type
        total_processed = 0
        total_failed = 0

        for cov_type in cov_types:
            console.print(
                f"\n[blue]Processing coverage type: {cov_type}[/blue]"
            )
            processed, failed = converter.process_samples_for_type_optimized(
                cov_fw_rv=cov_type,
                cluster_filter=cluster,
                num_workers=num_workers,
                skip_existing=skip_existing,
                batch_size=batch_size,
            )
            total_processed += processed
            total_failed += failed

            console.print(
                f"[green]Completed {cov_type}: {processed} processed, {failed} failed[/green]"
            )

        # Summary
        console.print("\n[bold green]Total Summary:[/bold green]")
        console.print(
            f"[green]Successfully processed: {total_processed}[/green]"
        )
        console.print(f"[red]Failed: {total_failed}[/red]")

        if total_failed > 0:
            raise typer.Exit(1)

    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def list_existing(
    parquet_path: Path = typer.Option(
        DEFAULT_PARQUETPATH,
        "--parquet-path",
        "-o",
        help="Path to Parquet directory",
    ),
    cov_type: str = typer.Option(
        "cov", "--type", "-t", help="Coverage type to check (cov, fw, rv)"
    ),
) -> None:
    """List existing parquet files for a coverage type."""
    setup_logging("INFO")

    if cov_type not in TABLENAMES:
        console.print(f"[red]Invalid coverage type: {cov_type}[/red]")
        console.print(
            f"[yellow]Valid types: {list(TABLENAMES.keys())}[/yellow]"
        )
        raise typer.Exit(1)

    root_path = parquet_path / TABLENAMES[cov_type]

    if not root_path.exists():
        console.print(
            f"[yellow]No parquet directory found for {cov_type}[/yellow]"
        )
        return

    # Create a dummy converter to use the method
    converter = ParquetConverter(Path("."), parquet_path, pl.DataFrame())
    existing_files = converter.get_existing_parquet_files(root_path)

    console.print(
        f"[blue]Found {len(existing_files)} existing parquet files for {cov_type}:[/blue]"
    )
    for file_id in sorted(existing_files):
        console.print(f"  {file_id}")


@app.command()
def clean(
    parquet_path: Path = typer.Option(
        DEFAULT_PARQUETPATH,
        "--parquet-path",
        "-o",
        help="Path to Parquet directory",
    ),
    cov_type: Optional[str] = typer.Option(
        None,
        "--type",
        "-t",
        help="Coverage type to clean (cov, fw, rv). If not specified, cleans all.",
    ),
    confirm: bool = typer.Option(
        False,
        "--confirm/--no-confirm",
        help="Confirm deletion without prompting",
    ),
) -> None:
    """Clean (delete) existing parquet files."""
    setup_logging("INFO")

    types_to_clean = [cov_type] if cov_type else list(TABLENAMES.keys())

    for ct in types_to_clean:
        if ct not in TABLENAMES:
            console.print(f"[red]Invalid coverage type: {ct}[/red]")
            continue

        root_path = parquet_path / TABLENAMES[ct]

        if not root_path.exists():
            console.print(
                f"[yellow]No parquet directory found for {ct}[/yellow]"
            )
            continue

        if not confirm:
            response = typer.confirm(f"Delete all parquet files for {ct}?")
            if not response:
                console.print(f"[yellow]Skipping {ct}[/yellow]")
                continue

        import shutil

        shutil.rmtree(root_path)
        console.print(f"[green]Deleted parquet files for {ct}[/green]")


if __name__ == "__main__":
    app()
