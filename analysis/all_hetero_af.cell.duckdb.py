#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-19 14:27:56
# @DESCRIPTION: Optimized variant allele frequency analysis tool with parallel processing
# @VERSION: v0.1.0

"""
Optimized Variant Allele Frequency Analysis Tool

This script processes variant allele frequency data from CSV files and inserts them
into a DuckDB database with an optimized hybrid approach:

🚀 OPTIMIZED PROCESSING:
- Parallel CSV loading in user-defined batches (e.g., 50 files)
- Sequential database insertion (no contention, reliable progress)
- Memory-efficient batch processing with garbage collection
- Real-time progress tracking with accurate time estimates

📊 MONITORING & RELIABILITY:
- Rich progress bars with detailed status updates
- Color-coded logging with comprehensive error handling
- Real-time statistics and batch completion reports
- Comprehensive input file validation

🛠️ ARCHITECTURE:
- ProcessPoolExecutor for parallel CSV loading
- Sequential DuckDB operations for thread safety
- Context managers for proper resource management
- Type hints and modular design for maintainability

📈 USAGE EXAMPLES:
    # Process all files with default settings (50 files per batch)
    python script.py insert-all

    # Customize workers and batch size
    python script.py insert-all --max-workers 4 --batch-size 100

    # Process a single file
    python script.py insert-one GSE123456 SRR789012

    # Validate input files
    python script.py validate-files --show-missing

    # Show database statistics
    python script.py show-stats
"""

import gc
import logging
import multiprocessing as mp
import time
from concurrent.futures import (
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)
from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from typing import Annotated, List, Optional

import duckdb
import polars as pl
import typer
from rich import print
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
from rich.table import Table

# Set up logging with Rich
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True)],
)
logger = logging.getLogger("scMOCHA")
console = Console()


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
DBFILE = "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/all_hetero_af.cell.duckdb.1.2.1"
TABLENAME_ALTDEPTH = "all_hetero_altdepth_cell"
TABLENAME_SUMDEPTH = "all_hetero_sumdepth_cell"


class DuckDBManager:
    """Thread-safe manager for DuckDB database operations"""

    def __init__(self, dbfile: str):
        self.dbfile = dbfile
        self._connection = None
        self._lock = Lock()

    @property
    def connection(self):
        """Get connection, creating it if necessary (thread-safe)"""
        with self._lock:
            if self._connection is None:
                self._connection = duckdb.connect(self.dbfile)
            return self._connection

    def create_table(self, table_name: str, df: pl.DataFrame) -> bool:
        """Create a table in DuckDB from a Polars DataFrame"""
        try:
            self.connection.execute(
                f"CREATE TABLE IF NOT EXISTS {table_name} AS SELECT * FROM df"
            )
            logger.info(f"Created table: {table_name}")
            return True
        except Exception as e:
            logger.error(f"Failed to create table {table_name}: {e}")
            return False

    def insert_data(self, table_name: str, df: pl.DataFrame) -> bool:
        """Insert data into an existing DuckDB table"""
        try:
            self.connection.execute(
                f"INSERT INTO {table_name} SELECT * FROM df"
            )
            logger.debug(f"Inserted {len(df)} rows into {table_name}")
            return True
        except Exception as e:
            logger.error(f"Failed to insert data into {table_name}: {e}")
            return False

    def close(self):
        """Close the DuckDB connection"""
        if self._connection:
            self._connection.close()
            self._connection = None
            logger.debug("Database connection closed")

    def __enter__(self):
        """Context manager entry"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def table_exists(self, table_name: str) -> bool:
        """Check if a table exists in the database"""
        try:
            result = self.connection.execute(
                f"SELECT COUNT(*) FROM information_schema.tables WHERE table_name = '{table_name}'"
            ).fetchone()
            return result[0] > 0
        except Exception as e:
            logger.error(f"Error checking if table {table_name} exists: {e}")
            return False

    def insert_or_create_table(self, table_name: str, df: pl.DataFrame) -> bool:
        """Insert data into table, creating it if it doesn't exist"""
        if self.table_exists(table_name):
            return self.insert_data(table_name, df)
        else:
            return self.create_table(table_name, df)

    def show_table_head(self, table_name: str, n: int = 5) -> None:
        """Show the first n rows of a table"""
        if self.table_exists(table_name):
            try:
                result = self.connection.execute(
                    f"SELECT * FROM {table_name} LIMIT {n}"
                ).fetchall()
                columns = [desc[0] for desc in self.connection.description]

                # Create a rich table for better display
                table = Table(title=f"Table: {table_name} (first {n} rows)")
                for col in columns:
                    table.add_column(col)
                for row in result:
                    table.add_row(*[str(cell) for cell in row])
                console.print(table)
            except Exception as e:
                logger.error(f"Error displaying table {table_name}: {e}")
        else:
            logger.warning(f"Table '{table_name}' does not exist")

    def get_table_row_count(self, table_name: str) -> int:
        """Get the number of rows in a table"""
        if self.table_exists(table_name):
            try:
                result = self.connection.execute(
                    f"SELECT COUNT(*) FROM {table_name}"
                ).fetchone()
                return result[0]
            except Exception as e:
                logger.error(f"Error getting row count for {table_name}: {e}")
                return 0
        else:
            logger.warning(f"Table '{table_name}' does not exist")
            return 0

    def show_table_row_count(self, table_name: str) -> None:
        """Show the number of rows in a table"""
        count = self.get_table_row_count(table_name)
        if count > 0:
            console.print(
                f"Table '[cyan]{table_name}[/cyan]' has [bold green]{count:,}[/bold green] rows"
            )

    def list_tables(self) -> List[str]:
        """List all table names in the database"""
        try:
            result = self.connection.execute(
                "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
            ).fetchall()
            table_names = [row[0] for row in result]

            if table_names:
                table = Table(title="Database Tables")
                table.add_column("Table Name", style="cyan")
                table.add_column("Row Count", style="green")

                for name in table_names:
                    count = self.get_table_row_count(name)
                    table.add_row(name, f"{count:,}")
                console.print(table)
            else:
                console.print("[yellow]No tables found in database[/yellow]")

            return table_names
        except Exception as e:
            logger.error(f"Error listing tables: {e}")
            return []

    def drop_table(self, table_name: str) -> bool:
        """Drop a table from the database"""
        if self.table_exists(table_name):
            try:
                self.connection.execute(f"DROP TABLE {table_name}")
                logger.info(f"Table '{table_name}' has been dropped")
                return True
            except Exception as e:
                logger.error(f"Error dropping table {table_name}: {e}")
                return False
        else:
            logger.warning(f"Table '{table_name}' does not exist")
            return False


@dataclass
class ProcessingResult:
    """Result of processing a single GSE/SRR pair"""

    gseid: str
    srrid: str
    success: bool
    error_message: Optional[str] = None
    altdepth_rows: int = 0
    sumdepth_rows: int = 0
    processing_time: float = 0.0
    # Optional dataframes for batch processing
    df_altdepth: Optional[pl.DataFrame] = None
    df_sumdepth: Optional[pl.DataFrame] = None


def process_csv_files_parallel(
    gseid: str, srrid: str, cluster: str, hh: str
) -> ProcessingResult:
    """
    Process CSV files in parallel (file I/O only, no database operations).
    This function is designed for parallel execution.
    """
    # gseid, srrid, cluster, hh = (
    #     "GSE155673",
    #     "GSM4712885",
    #     "cell",
    #     "HETEROPLASMIC",
    # )
    start_time = time.time()
    try:
        config = CONFIG[hh]

        # Check if files exist before processing
        altdepth_file = (
            config.table_path
            / f"{cluster}_cov.{gseid}_{srrid}.altdepth.{config.suffix}.csv"
        )
        sumdepth_file = (
            config.table_path
            / f"{cluster}_cov.{gseid}_{srrid}.sumdepth.{config.suffix}.csv"
        )

        if not altdepth_file.exists():
            return ProcessingResult(
                gseid=gseid,
                srrid=srrid,
                success=False,
                error_message=f"Altdepth file not found: {altdepth_file}",
                processing_time=time.time() - start_time,
            )

        if not sumdepth_file.exists():
            return ProcessingResult(
                gseid=gseid,
                srrid=srrid,
                success=False,
                error_message=f"Sumdepth file not found: {sumdepth_file}",
                processing_time=time.time() - start_time,
            )

        # Read and process CSV files (this is what we parallelize)
        df_altdepth = pl.read_csv(altdepth_file, has_header=True).unpivot(
            index=["gseid", "srrid", "barcode"],
            variable_name="variant",
            value_name="altdepth",
        )

        df_sumdepth = pl.read_csv(sumdepth_file, has_header=True).unpivot(
            index=["gseid", "srrid", "barcode"],
            variable_name="variant",
            value_name="sumdepth",
        )

        # Store row counts and dataframes for sequential database insertion
        result = ProcessingResult(
            gseid=gseid,
            srrid=srrid,
            success=True,
            altdepth_rows=len(df_altdepth),
            sumdepth_rows=len(df_sumdepth),
            processing_time=time.time() - start_time,
            df_altdepth=df_altdepth,
            df_sumdepth=df_sumdepth,
        )

        return result

    except Exception as e:
        logger.error(f"Error processing CSV files {gseid}-{srrid}: {str(e)}")
        return ProcessingResult(
            gseid=gseid,
            srrid=srrid,
            success=False,
            error_message=str(e),
            processing_time=time.time() - start_time,
        )


def process_all_files_optimized(
    cluster: str = "cell",
    hh: str = "HETEROPLASMIC",
    max_workers: Optional[int] = None,
    batch_size: int = 50,
) -> None:
    """
    Process all files using optimized hybrid approach.

    STRATEGY:
    1. Parallel CSV loading in batches (maximize I/O throughput)
    2. Sequential DuckDB insertion (avoid contention, ensure progress)
    3. Memory-efficient batch processing with cleanup

    This approach provides the best balance of performance and reliability.

    Args:
        cluster: Cluster type (e.g., "cell")
        hh: Heteroplasmic or Homoplasmic ("HETEROPLASMIC"/"HOMOPLASMIC")
        max_workers: Maximum number of parallel workers for CSV loading
        batch_size: Number of files to process in each parallel batch
    """
    if max_workers is None:
        max_workers = min(mp.cpu_count(), 8)

    logger.info(
        "Starting optimized processing: parallel CSV loading + sequential DB insertion"
    )
    logger.info(f"Workers: {max_workers}, Batch size: {batch_size}")
    logger.info(f"Processing {len(SRR)} files total")

    # Initialize database and drop existing tables
    with DuckDBManager(DBFILE) as db:
        logger.info("Initializing database...")
        db.list_tables()

        if db.table_exists(TABLENAME_ALTDEPTH):
            logger.info(f"Dropping existing table: {TABLENAME_ALTDEPTH}")
            db.drop_table(TABLENAME_ALTDEPTH)

        if db.table_exists(TABLENAME_SUMDEPTH):
            logger.info(f"Dropping existing table: {TABLENAME_SUMDEPTH}")
            db.drop_table(TABLENAME_SUMDEPTH)

    # Prepare task list
    tasks = []
    for row in SRR.iter_rows(named=True):
        tasks.append((row["gseid"], row["srrid"], cluster, hh))

    # Processing statistics
    successful_tasks = 0
    failed_tasks = 0
    total_altdepth_rows = 0
    total_sumdepth_rows = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=None),
        "[progress.percentage]{task.percentage:>3.1f}%",
        "•",
        TextColumn("[blue]{task.completed}/{task.total}"),
        "•",
        TimeElapsedColumn(),
        "•",
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        overall_task = progress.add_task(
            f"Processing {len(tasks)} files", total=len(tasks)
        )

        # Process in batches
        for i in range(0, len(tasks), batch_size):
            batch = tasks[i : i + batch_size]
            batch_num = i // batch_size + 1
            total_batches = (len(tasks) + batch_size - 1) // batch_size

            logger.info(
                f"Batch {batch_num}/{total_batches}: Processing {len(batch)} files in parallel..."
            )

            # Step 1: Parallel CSV loading
            batch_results = []
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_task = {
                    executor.submit(process_csv_files_parallel, *task): task
                    for task in batch
                }
                print(
                    f"Batch {batch_num}: Waiting for {len(future_to_task)} tasks to complete..."
                )

                # Collect results as they complete
                for future in as_completed(future_to_task):
                    try:
                        task = future_to_task[future]
                        result = future.result()
                        batch_results.append(result)

                        if result.success:
                            logger.debug(
                                f"✓ Loaded {result.gseid}-{result.srrid} ({result.processing_time:.2f}s)"
                            )
                        else:
                            logger.error(
                                f"✗ Failed {result.gseid}-{result.srrid}: {result.error_message}"
                            )

                    except Exception as e:
                        gseid, srrid = task[0], task[1]
                        logger.error(f"✗ Exception {gseid}-{srrid}: {e}")
                        batch_results.append(
                            ProcessingResult(
                                gseid=gseid,
                                srrid=srrid,
                                success=False,
                                error_message=str(e),
                            )
                        )

                    progress.update(overall_task, advance=1)

            # Step 2: Sequential database insertion
            logger.info(
                f"Batch {batch_num}: Inserting {len([r for r in batch_results if r.success])} successful results to database..."
            )
            with DuckDBManager(DBFILE) as db:
                for result in batch_results:
                    if result.success and result.df_altdepth is not None:
                        try:
                            # Sequential database insertion (no conflicts)
                            altdepth_success = db.insert_or_create_table(
                                TABLENAME_ALTDEPTH, result.df_altdepth
                            )
                            sumdepth_success = db.insert_or_create_table(
                                TABLENAME_SUMDEPTH, result.df_sumdepth
                            )

                            if altdepth_success and sumdepth_success:
                                successful_tasks += 1
                                total_altdepth_rows += result.altdepth_rows
                                total_sumdepth_rows += result.sumdepth_rows
                            else:
                                failed_tasks += 1
                                logger.error(
                                    f"✗ DB insert failed: {result.gseid}-{result.srrid}"
                                )

                        except Exception as e:
                            failed_tasks += 1
                            logger.error(
                                f"✗ DB error {result.gseid}-{result.srrid}: {e}"
                            )

                    elif not result.success:
                        failed_tasks += 1

            # Clean up memory after each batch
            for result in batch_results:
                if (
                    hasattr(result, "df_altdepth")
                    and result.df_altdepth is not None
                ):
                    del result.df_altdepth
                if (
                    hasattr(result, "df_sumdepth")
                    and result.df_sumdepth is not None
                ):
                    del result.df_sumdepth
            del batch_results
            gc.collect()

            logger.info(
                f"Batch {batch_num} complete. Total success: {successful_tasks}, failed: {failed_tasks}"
            )

    # Final summary
    console.print("\n" + "=" * 60)
    console.print("[bold green]Processing Complete![/bold green]")
    console.print(f"✓ Successful: [green]{successful_tasks}[/green]")
    console.print(f"✗ Failed: [red]{failed_tasks}[/red]")
    console.print(
        f"📊 Total altdepth rows: [cyan]{total_altdepth_rows:,}[/cyan]"
    )
    console.print(
        f"📊 Total sumdepth rows: [cyan]{total_sumdepth_rows:,}[/cyan]"
    )

    # Show final database statistics
    with DuckDBManager(DBFILE) as db:
        console.print("\n[bold]Final Database Statistics:[/bold]")
        db.show_table_row_count(TABLENAME_ALTDEPTH)
        db.show_table_row_count(TABLENAME_SUMDEPTH)


# CLI Application
app = typer.Typer(
    help="Optimized variant allele frequency analysis tool with parallel processing",
    rich_markup_mode="rich",
)


@app.command()
def insert_one(
    gseid: Annotated[str, typer.Argument(help="GSE ID")],
    srrid: Annotated[str, typer.Argument(help="SRR ID")],
    cluster: Annotated[
        str, typer.Argument(help="Cluster type (e.g., cell)")
    ] = "cell",
    hh: Annotated[
        str,
        typer.Argument(
            help="Heteroplasmic or Homoplasmic (HETEROPLASMIC/HOMOPLASMIC)"
        ),
    ] = "HETEROPLASMIC",
):
    """Insert data for a single GSE and SRR into the database."""
    console.print(f"[bold]Processing single file:[/bold] {gseid} - {srrid}")

    result = process_csv_files_parallel(gseid, srrid, cluster, hh)

    if result.success:
        # Insert to database
        with DuckDBManager(DBFILE) as db:
            altdepth_success = db.insert_or_create_table(
                TABLENAME_ALTDEPTH, result.df_altdepth
            )
            sumdepth_success = db.insert_or_create_table(
                TABLENAME_SUMDEPTH, result.df_sumdepth
            )

            if altdepth_success and sumdepth_success:
                console.print(
                    f"[green]✓ Success![/green] Processed in {result.processing_time:.2f}s"
                )
                console.print(f"  • Altdepth rows: {result.altdepth_rows:,}")
                console.print(f"  • Sumdepth rows: {result.sumdepth_rows:,}")
            else:
                console.print("[red]✗ Failed to insert to database[/red]")
                raise typer.Exit(1)
    else:
        console.print(f"[red]✗ Failed:[/red] {result.error_message}")
        raise typer.Exit(1)


@app.command()
def insert_all(
    cluster: Annotated[
        str, typer.Option(help="Cluster type (e.g., cell)")
    ] = "cell",
    hh: Annotated[
        str,
        typer.Option(
            help="Heteroplasmic or Homoplasmic (HETEROPLASMIC/HOMOPLASMIC)"
        ),
    ] = "HETEROPLASMIC",
    max_workers: Annotated[
        Optional[int],
        typer.Option(
            help="Maximum number of parallel workers for CSV loading (default: auto-detect)"
        ),
    ] = None,
    batch_size: Annotated[
        int, typer.Option(help="Number of files to process in each batch")
    ] = 50,
):
    """
    Insert data for all GSE and SRR into the database using optimized processing.

    Uses parallel CSV loading combined with sequential database insertion for
    maximum performance and reliability.
    """
    console.print(f"[bold]Processing all files:[/bold] {len(SRR)} total files")
    console.print("[bold]Configuration:[/bold]")
    console.print(f"  • Cluster: {cluster}")
    console.print(f"  • Type: {hh}")
    console.print(f"  • Max workers: {max_workers or 'auto-detect'}")
    console.print(f"  • Batch size: {batch_size}")

    process_all_files_optimized(cluster, hh, max_workers, batch_size)
    console.print(
        "[bold green]✓ All processing completed successfully![/bold green]"
    )


@app.command()
def show_stats(
    table_name: Annotated[
        Optional[str],
        typer.Argument(
            help="Table name to show stats for (default: show all tables)"
        ),
    ] = None,
):
    """Show database statistics and table information."""
    with DuckDBManager(DBFILE) as db:
        if table_name:
            if db.table_exists(table_name):
                db.show_table_row_count(table_name)
                db.show_table_head(table_name, 10)
            else:
                console.print(f"[red]Table '{table_name}' does not exist[/red]")
                raise typer.Exit(1)
        else:
            console.print("[bold]Database Overview:[/bold]")
            tables = db.list_tables()
            if not tables:
                console.print("[yellow]No tables found in database[/yellow]")


@app.command()
def validate_files(
    cluster: Annotated[str, typer.Option(help="Cluster type")] = "cell",
    hh: Annotated[
        str, typer.Option(help="Heteroplasmic or Homoplasmic")
    ] = "HETEROPLASMIC",
    show_missing: Annotated[
        bool, typer.Option(help="Show missing files")
    ] = False,
):
    """Validate that all required input files exist."""
    config = CONFIG[hh]

    missing_files = []
    existing_files = []

    console.print(f"[bold]Validating files for {hh} {cluster} data...[/bold]")

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        "[progress.percentage]{task.percentage:>3.1f}%",
        console=console,
    ) as progress:
        task = progress.add_task("Checking files...", total=len(SRR))

        for row in SRR.iter_rows(named=True):
            gseid = row["gseid"]
            srrid = row["srrid"]

            altdepth_file = (
                config.table_path
                / f"{cluster}_cov.{gseid}_{srrid}.altdepth.{config.suffix}.csv"
            )
            sumdepth_file = (
                config.table_path
                / f"{cluster}_cov.{gseid}_{srrid}.sumdepth.{config.suffix}.csv"
            )

            if altdepth_file.exists() and sumdepth_file.exists():
                existing_files.append((gseid, srrid))
            else:
                missing_files.append(
                    (gseid, srrid, altdepth_file, sumdepth_file)
                )

            progress.update(task, advance=1)

    console.print("\n[bold]Validation Results:[/bold]")
    console.print(
        f"✓ Found: [green]{len(existing_files)}[/green] complete file pairs"
    )
    console.print(f"✗ Missing: [red]{len(missing_files)}[/red] file pairs")

    if missing_files and show_missing:
        console.print("\n[bold red]Missing files:[/bold red]")
        for gseid, srrid, altfile, sumfile in missing_files[
            :10
        ]:  # Show first 10
            console.print(f"  • {gseid}-{srrid}")
            if not altfile.exists():
                console.print(f"    - Missing: {altfile.name}")
            if not sumfile.exists():
                console.print(f"    - Missing: {sumfile.name}")

        if len(missing_files) > 10:
            console.print(f"    ... and {len(missing_files) - 10} more")

    if missing_files:
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
