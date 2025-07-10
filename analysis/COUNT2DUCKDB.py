#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-06-20 16:04:43
# @DESCRIPTION:
# @VERSION: v0.0.1


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
from typing import Annotated, List, Literal, Optional

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

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)


class DuckDBManager:
    """Thread-safe manager for DuckDB database operations"""

    def __init__(self, dbfile: str):
        self.dbfile = Path(dbfile)
        # Ensure the directory exists
        self.dbfile.parent.mkdir(parents=True, exist_ok=True)
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


TABLEPATH = Path("/home/liuc9/github/scMOCHA-data/analysis/zzz/db/TABLES")
DBPATH = Path("/home/liuc9/github/scMOCHA-data/analysis/zzz/db/DUCKDB")
DBFILE = DBPATH / "cov.duckdb"
TABLENAMES = {"cov": "covall", "fw": "covfw", "rv": "covrv"}


def load_data_pl(
    gseid: str, srrid: str, cluster: str, cov_fw_rv: Literal["cov", "fw", "rv"]
) -> pl.DataFrame:
    """Load data into DuckDB from Polars DataFrame"""
    # gseid, srrid, cluster, cov_fw_rv = "GSE147794", "GSM4446059", "cell", "cov"
    covall_filename = f"{cluster}_{cov_fw_rv}.{gseid}_{srrid}.csv"
    covall_df = (
        pl.read_csv(TABLEPATH / covall_filename)
        .lazy()
        .with_columns(
            [pl.lit(gseid).alias("gseid"), pl.lit(srrid).alias("srrid")]
        )
        .collect()
    )
    # covall_df.unpivot(
    #     index=["gseid", "srrid", "base", "barcode"],
    #     variable_name="pos",
    #     value_name="coverage",
    # )
    return covall_df


def insert_data_to_duckdb(
    gseid: str,
    srrid: str,
    cluster: str,
    cov_fw_rv: Literal["cov", "fw", "rv"],
) -> bool:
    """Insert data into DuckDB, creating table if it doesn't exist"""
    covall_df = load_data_pl(gseid, srrid, cluster, cov_fw_rv)
    table_name = TABLENAMES[cov_fw_rv]
    with DuckDBManager(DBFILE) as db_manager:
        db_manager.insert_or_create_table(table_name, covall_df)
    del covall_df  # Clean up memory
    gc.collect()  # Force garbage collection


# insert_data_to_duckdb(
#     gseid="GSE147794", srrid="GSM4446059", cluster="cell", cov_fw_rv="cov"
# )


def load_batch_data(batch_tasks: List[tuple], cov_fw_rv: str) -> pl.DataFrame:
    """Load a batch of data in parallel and concatenate into a single DataFrame"""
    batch_dfs = []

    with ThreadPoolExecutor(max_workers=min(20, len(batch_tasks))) as executor:
        # Submit all tasks in the batch
        future_to_task = {
            executor.submit(load_data_pl, gseid, srrid, cluster, cov_fw_rv): (
                gseid,
                srrid,
            )
            for gseid, srrid, cluster in batch_tasks
        }

        # Collect results with progress tracking
        for future in as_completed(future_to_task):
            gseid, srrid = future_to_task[future]
            try:
                df = future.result()
                batch_dfs.append(df)
                logger.debug(f"Loaded data for {gseid}_{srrid}")
            except Exception as e:
                logger.error(f"Failed to load data for {gseid}_{srrid}: {e}")

    # Concatenate all DataFrames in the batch
    if batch_dfs:
        combined_df = pl.concat(batch_dfs, how="vertical")
        logger.info(
            f"Combined {len(batch_dfs)} datasets into batch with {len(combined_df)} rows"
        )
        return combined_df
    else:
        return pl.DataFrame()


def process_data_in_batches(
    tasks: List[tuple],
    batch_size: int = 5,
    cov_types: List[str] = ["cov", "fw", "rv"],
) -> None:
    """Process all data in batches with parallel loading and sequential database insertion"""

    console.print(
        f"\n[bold blue]Processing {len(tasks)} tasks in batches of {batch_size}[/bold blue]"
    )

    # Create batches
    batches = [
        tasks[i : i + batch_size] for i in range(0, len(tasks), batch_size)
    ]
    console.print(f"Created {len(batches)} batches")

    # Process each coverage type
    for cov_type in cov_types:
        console.print(
            f"\n[bold green]Processing coverage type: {cov_type}[/bold green]"
        )
        table_name = TABLENAMES[cov_type]

        # Initialize progress bar
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
        ) as progress:
            task_id = progress.add_task(
                f"Processing {cov_type}", total=len(batches)
            )

            # Process each batch
            for batch_idx, batch in enumerate(batches, 1):
                progress.update(
                    task_id,
                    description=f"Processing {cov_type} - Batch {batch_idx}/{len(batches)}",
                )

                # Load batch data in parallel
                batch_df = load_batch_data(batch, cov_type)

                # Insert into database sequentially (thread-safe)
                if not batch_df.is_empty():
                    with DuckDBManager(DBFILE) as db_manager:
                        success = db_manager.insert_or_create_table(
                            table_name, batch_df
                        )
                        if success:
                            logger.info(
                                f"Inserted batch {batch_idx} into {table_name}"
                            )
                        else:
                            logger.error(
                                f"Failed to insert batch {batch_idx} into {table_name}"
                            )

                # Clean up memory
                del batch_df
                gc.collect()

                progress.update(task_id, advance=1)

        console.print(f"[green]Completed processing {cov_type}[/green]")


app = typer.Typer(help="scMOCHA Data Processing - Batch load data into DuckDB")


@app.command("process")
def process_command(
    batch_size: Annotated[
        int, typer.Option(help="Batch size for parallel processing")
    ] = 5,
    max_samples: Annotated[
        Optional[int],
        typer.Option(help="Maximum number of samples to process (for testing)"),
    ] = None,
    cov_types: Annotated[
        List[str], typer.Option(help="Coverage types to process")
    ] = ["cov", "fw", "rv"],
    show_tables: Annotated[
        bool,
        typer.Option(help="Show database tables before and after processing"),
    ] = True,
):
    """Process scMOCHA data in batches and load into DuckDB"""

    # Determine tasks to process
    sample_data = SRR.head(max_samples) if max_samples else SRR
    current_tasks = []
    cluster = "cell"

    for row in sample_data.iter_rows(named=True):
        current_tasks.append((row["gseid"], row["srrid"], cluster))

    console.print("[bold cyan]Starting scMOCHA Data Processing[/bold cyan]")
    console.print(
        f"Processing {len(current_tasks)} samples with batch size {batch_size}"
    )

    # Show initial database state
    if show_tables:
        with DuckDBManager(DBFILE) as db_manager:
            console.print("\n[bold]Current database state:[/bold]")
            db_manager.list_tables()

    # Process data in batches
    start_time = time.time()
    process_data_in_batches(
        current_tasks, batch_size=batch_size, cov_types=cov_types
    )
    end_time = time.time()

    # Show final database state
    if show_tables:
        with DuckDBManager(DBFILE) as db_manager:
            console.print("\n[bold]Final database state:[/bold]")
            db_manager.list_tables()

    console.print(
        f"\n[bold green]Processing completed in {end_time - start_time:.2f} seconds[/bold green]"
    )


@app.command("list-tables")
def list_tables_command():
    """List all tables in the DuckDB database"""
    with DuckDBManager(DBFILE) as db_manager:
        db_manager.list_tables()


@app.command("show-table")
def show_table_command(
    table_name: Annotated[
        str, typer.Argument(help="Name of the table to show")
    ],
    rows: Annotated[int, typer.Option(help="Number of rows to display")] = 5,
):
    """Show the first few rows of a specific table"""
    with DuckDBManager(DBFILE) as db_manager:
        db_manager.show_table_head(table_name, rows)


@app.command("drop-table")
def drop_table_command(
    table_name: Annotated[
        str, typer.Argument(help="Name of the table to drop")
    ],
    confirm: Annotated[
        bool, typer.Option("--confirm", help="Confirm table deletion")
    ] = False,
):
    """Drop a table from the database"""
    if not confirm:
        if not typer.confirm(
            f"Are you sure you want to drop table '{table_name}'?"
        ):
            console.print("[yellow]Operation cancelled[/yellow]")
            return

    with DuckDBManager(DBFILE) as db_manager:
        db_manager.drop_table(table_name)


@app.command("test")
def test_command(
    samples: Annotated[
        int, typer.Option(help="Number of samples to test with")
    ] = 2,
    batch_size: Annotated[int, typer.Option(help="Batch size for testing")] = 5,
):
    """Test the processing with a small subset of data"""
    console.print(f"[bold yellow]Testing with {samples} samples[/bold yellow]")
    process_command(
        batch_size=batch_size, max_samples=samples, show_tables=True
    )


def main_batch_processing():
    """Main function to run the batch processing"""
    console.print("[bold cyan]Starting scMOCHA Data Processing[/bold cyan]")

    # Create tasks from SRR data
    tasks = []
    cluster = "cell"

    for row in SRR.iter_rows(named=True):
        tasks.append((row["gseid"], row["srrid"], cluster))

    console.print(f"Processing {len(tasks)} samples with batch size 20")

    # Show initial database state
    with DuckDBManager(DBFILE) as db_manager:
        console.print("\n[bold]Current database state:[/bold]")
        db_manager.list_tables()

    # Process data in batches
    start_time = time.time()
    process_data_in_batches(tasks, batch_size=20)
    end_time = time.time()

    # Show final database state
    with DuckDBManager(DBFILE) as db_manager:
        console.print("\n[bold]Final database state:[/bold]")
        db_manager.list_tables()

    console.print(
        f"\n[bold green]Processing completed in {end_time - start_time:.2f} seconds[/bold green]"
    )


@app.command("insert-all")
def insert_all_command(
    batch_size: Annotated[
        int, typer.Option(help="Batch size for parallel processing")
    ] = 5,
    cov_types: Annotated[
        List[str], typer.Option(help="Coverage types to process")
    ] = ["cov", "fw", "rv"],
    drop_existing: Annotated[
        bool,
        typer.Option(
            "--drop-existing", help="Drop existing tables before inserting"
        ),
    ] = True,
    confirm: Annotated[
        bool,
        typer.Option("--confirm", help="Confirm dropping all existing tables"),
    ] = False,
):
    """Drop all existing tables and insert all data from SRR dataset"""

    # Prepare all tasks from SRR dataset
    all_tasks = []
    cluster = "cell"

    for row in SRR.iter_rows(named=True):
        all_tasks.append((row["gseid"], row["srrid"], cluster))

    console.print("[bold cyan]Insert All Data Command[/bold cyan]")
    console.print(f"Total samples to process: {len(all_tasks)}")
    console.print(f"Coverage types: {', '.join(cov_types)}")
    console.print(f"Batch size: {batch_size}")

    # Show current database state
    with DuckDBManager(DBFILE) as db_manager:
        console.print("\n[bold]Current database state:[/bold]")
        existing_tables = db_manager.list_tables()

    # Drop existing tables if requested
    if drop_existing and existing_tables:
        if not confirm:
            console.print(
                f"\n[bold red]This will drop {len(existing_tables)} existing tables:[/bold red]"
            )
            for table in existing_tables:
                console.print(f"  - {table}")

            if not typer.confirm("\nAre you sure you want to proceed?"):
                console.print("[yellow]Operation cancelled[/yellow]")
                return

        console.print("\n[bold red]Dropping existing tables...[/bold red]")
        with DuckDBManager(DBFILE) as db_manager:
            for table in existing_tables:
                if db_manager.drop_table(table):
                    console.print(f"[red]Dropped table: {table}[/red]")

        console.print("[green]All existing tables dropped[/green]")

    # Process all data in batches
    console.print("\n[bold blue]Starting full data processing...[/bold blue]")
    start_time = time.time()
    process_data_in_batches(
        all_tasks, batch_size=batch_size, cov_types=cov_types
    )
    end_time = time.time()

    # Show final database state
    with DuckDBManager(DBFILE) as db_manager:
        console.print("\n[bold]Final database state:[/bold]")
        db_manager.list_tables()

    console.print(
        f"\n[bold green]Full processing completed in {end_time - start_time:.2f} seconds[/bold green]"
    )
    console.print(
        f"[bold green]Successfully processed {len(all_tasks)} samples[/bold green]"
    )


@app.command("show-stats")
def show_stats_command(
    detailed: Annotated[
        bool,
        typer.Option(
            "--detailed", help="Show detailed statistics for each table"
        ),
    ] = False,
):
    """Show statistics for all tables in the database"""

    with DuckDBManager(DBFILE) as db_manager:
        console.print("[bold cyan]Database Statistics[/bold cyan]")

        # Get all tables
        table_names = db_manager.list_tables()

        if not table_names:
            console.print("[yellow]No tables found in database[/yellow]")
            return

        # Calculate total rows across all tables
        total_rows = 0
        table_stats = []

        for table_name in table_names:
            row_count = db_manager.get_table_row_count(table_name)
            total_rows += row_count
            table_stats.append((table_name, row_count))

        # Show summary statistics
        console.print("\n[bold]Summary Statistics:[/bold]")
        console.print(f"Total tables: [cyan]{len(table_names)}[/cyan]")
        console.print(
            f"Total rows across all tables: [green]{total_rows:,}[/green]"
        )

        if total_rows > 0:
            avg_rows = total_rows / len(table_names)
            console.print(
                f"Average rows per table: [yellow]{avg_rows:,.0f}[/yellow]"
            )

        # Show detailed statistics if requested
        if detailed:
            console.print("\n[bold]Detailed Table Statistics:[/bold]")

            # Sort tables by row count (descending)
            table_stats.sort(key=lambda x: x[1], reverse=True)

            detailed_table = Table(title="Table Details")
            detailed_table.add_column("Table Name", style="cyan")
            detailed_table.add_column(
                "Row Count", style="green", justify="right"
            )
            detailed_table.add_column(
                "Percentage", style="yellow", justify="right"
            )

            for table_name, row_count in table_stats:
                percentage = (
                    (row_count / total_rows * 100) if total_rows > 0 else 0
                )
                detailed_table.add_row(
                    table_name, f"{row_count:,}", f"{percentage:.1f}%"
                )

            console.print(detailed_table)

            # Show database file size if possible
            try:
                db_size = db_manager.dbfile.stat().st_size
                db_size_mb = db_size / (1024 * 1024)
                console.print(
                    f"\n[bold]Database file size:[/bold] [cyan]{db_size_mb:.2f} MB[/cyan]"
                )
            except Exception as e:
                logger.debug(f"Could not get database file size: {e}")


@app.command("help")
def help_command():
    """Show detailed help and usage examples"""

    console.print("[bold cyan]scMOCHA Data Processing Tool[/bold cyan]")
    console.print(
        "\nA comprehensive tool for processing scMOCHA data and loading it into DuckDB with batch processing capabilities.\n"
    )

    console.print("[bold]Available Commands:[/bold]")

    help_table = Table()
    help_table.add_column("Command", style="cyan", width=15)
    help_table.add_column("Description", style="white", width=60)
    help_table.add_column("Example", style="green", width=40)

    help_table.add_row(
        "process",
        "Process a subset of data with customizable options",
        "python COUNT2DUCKDB.py process --batch-size 10 --max-samples 50",
    )
    help_table.add_row(
        "insert-all",
        "Drop all tables and process the entire SRR dataset",
        "python COUNT2DUCKDB.py insert-all --batch-size 20 --confirm",
    )
    help_table.add_row(
        "test",
        "Test processing with a small subset of data",
        "python COUNT2DUCKDB.py test --samples 3 --batch-size 5",
    )
    help_table.add_row(
        "list-tables",
        "List all tables in the database",
        "python COUNT2DUCKDB.py list-tables",
    )
    help_table.add_row(
        "show-table",
        "Show the first few rows of a specific table",
        "python COUNT2DUCKDB.py show-table covall --rows 10",
    )
    help_table.add_row(
        "show-stats",
        "Show database statistics and table row counts",
        "python COUNT2DUCKDB.py show-stats --detailed",
    )
    help_table.add_row(
        "drop-table",
        "Drop a specific table from the database",
        "python COUNT2DUCKDB.py drop-table covall --confirm",
    )
    help_table.add_row(
        "help", "Show this help message", "python COUNT2DUCKDB.py help"
    )

    console.print(help_table)

    console.print("\n[bold]Key Features:[/bold]")
    console.print(
        "• [green]Batch Processing:[/green] Process data in configurable batches (default: 20)"
    )
    console.print(
        "• [green]Parallel Loading:[/green] Load multiple files simultaneously using ThreadPoolExecutor"
    )
    console.print(
        "• [green]Sequential Insertion:[/green] Thread-safe database operations"
    )
    console.print(
        "• [green]Progress Tracking:[/green] Real-time progress bars and detailed logging"
    )
    console.print(
        "• [green]Memory Management:[/green] Automatic garbage collection between batches"
    )
    console.print(
        "• [green]Rich Output:[/green] Beautiful console output with colors and tables"
    )

    console.print("\n[bold]Coverage Types:[/bold]")
    console.print(
        "• [cyan]cov:[/cyan] Overall coverage data → Table: [yellow]covall[/yellow]"
    )
    console.print(
        "• [cyan]fw:[/cyan] Forward strand coverage → Table: [yellow]covfw[/yellow]"
    )
    console.print(
        "• [cyan]rv:[/cyan] Reverse strand coverage → Table: [yellow]covrv[/yellow]"
    )

    console.print("\n[bold]Database Configuration:[/bold]")
    console.print(f"• Database file: [cyan]{DBFILE}[/cyan]")
    console.print(f"• Source data: [cyan]{SRR_FILENAME}[/cyan]")
    console.print(f"• Data path: [cyan]{TABLEPATH}[/cyan]")


if __name__ == "__main__":
    app()
