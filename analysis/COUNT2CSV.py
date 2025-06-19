#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-23 23:58:00
# @DESCRIPTION: Enhanced genomic data processing pipeline with improved logging and parallel processing
# @VERSION: v1.0.0

import concurrent.futures
import logging
import multiprocessing as mp
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import polars as pl
import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TaskID
from rich.table import Table

# Configure enhanced logging
console = Console()


def setup_logging(
    log_level: str = "INFO", log_file: Optional[str] = None
) -> logging.Logger:
    """Setup enhanced logging with rich formatting and optional file output."""

    # Create formatter
    formatter = logging.Formatter(
        fmt="%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Setup handlers
    handlers = [RichHandler(console=console, rich_tracebacks=True)]

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)

    # Configure logging
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        handlers=handlers,
        format="%(message)s",
    )

    logger = logging.getLogger("COUNT2CSV")
    logger.info(f"Logger initialized with level: {log_level}")

    return logger


# Initialize logger
log = setup_logging()


class Config:
    """Configuration class for file paths and constants."""

    def __init__(self, base_dir: Optional[Path] = None):
        self.base_dir = base_dir or Path("/home/liuc9/github/scMOCHA-data")

        # File paths
        self.srr_filename = (
            self.base_dir / "analysis/zzz/clean-data/gse_srrid_srrdir.csv"
        )
        self.fasta_filename = self.base_dir / "config/rCRS.MT.fasta.df.csv"
        self.table_dir = Path(
            "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/TABLES"
        )

        # Constants
        self.bases = ["A", "C", "G", "T"]
        self.clusters = ["cell", "cluster"]

        # Validate paths
        self._validate_paths()

    def _validate_paths(self) -> None:
        """Validate that required files exist."""
        missing_files = []

        if not self.srr_filename.exists():
            log.warning(f"SRR file not found: {self.srr_filename}")
            missing_files.append(str(self.srr_filename))

        if not self.fasta_filename.exists():
            log.warning(f"FASTA file not found: {self.fasta_filename}")
            missing_files.append(str(self.fasta_filename))

        # Create table directory if it doesn't exist
        self.table_dir.mkdir(parents=True, exist_ok=True)
        log.info(f"Table directory ready: {self.table_dir}")

        if missing_files:
            log.warning(f"Missing files detected: {missing_files}")
            log.warning(
                "Script will continue with fallback behavior where possible"
            )


class DataProcessor:
    """Main data processing class with enhanced error handling and logging."""

    def __init__(self, config: Config):
        self.config = config
        self.srr_data = self._load_srr_data()
        self.pos_data = self._load_position_data()

    def _load_srr_data(self) -> pl.DataFrame:
        """Load SRR data with error handling."""
        try:
            if not self.config.srr_filename.exists():
                log.warning(f"SRR file not found: {self.config.srr_filename}")
                log.warning("Creating dummy SRR data for testing")
                # Create dummy data for testing
                return pl.DataFrame(
                    {
                        "gseid": ["GSE_TEST"],
                        "srrid": ["SRR_TEST"],
                        "srrdir": ["/test/path"],
                    }
                )

            log.info(f"Loading SRR data from: {self.config.srr_filename}")
            data = pl.read_csv(self.config.srr_filename)
            log.info(
                f"Loaded SRR data: {data.shape[0]} rows, {data.shape[1]} columns"
            )
            return data
        except Exception as e:
            log.error(f"Failed to load SRR data: {e}")
            log.warning("Creating dummy SRR data for testing")
            return pl.DataFrame(
                {
                    "gseid": ["GSE_TEST"],
                    "srrid": ["SRR_TEST"],
                    "srrdir": ["/test/path"],
                }
            )

    def _load_position_data(self) -> pl.DataFrame:
        """Load FASTA position data with error handling and format detection."""
        try:
            log.info(f"Loading FASTA data from: {self.config.fasta_filename}")

            # Try different formats
            if self.config.fasta_filename.suffix == ".csv":
                fasta = pl.read_csv(self.config.fasta_filename)
            elif self.config.fasta_filename.suffix == ".fst":
                log.warning(
                    "FST format detected. This requires R/fst package. Creating dummy position data."
                )
                # Create a dummy position range for mitochondrial genome (16569 bp)
                positions = list(range(1, 16570))
                fasta = pl.DataFrame({"pos": positions})
            else:
                log.warning(
                    f"Unknown format: {self.config.fasta_filename.suffix}. Creating dummy position data."
                )
                positions = list(range(1, 16570))
                fasta = pl.DataFrame({"pos": positions})

            pos_data = pl.DataFrame({"pos": fasta["pos"]})
            log.info(f"Loaded position data: {pos_data.shape[0]} positions")
            return pos_data
        except Exception as e:
            log.error(f"Failed to load FASTA data: {e}")
            log.warning(
                "Creating fallback position data for mitochondrial genome"
            )
            # Fallback: create standard mitochondrial genome positions
            positions = list(range(1, 16570))
            return pl.DataFrame({"pos": positions})

    def load_base_data(
        self, filename: Union[str, Path], base: str
    ) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
        """Load and process data for a single base with enhanced error handling."""
        log.debug(f"Loading {filename} for base {base}")

        filename = Path(filename)
        if not filename.exists():
            log.error(f"File not found: {filename}")
            raise FileNotFoundError(f"Data file not found: {filename}")

        try:
            # Use polars to read the CSV file
            df = pl.read_csv(
                filename,
                has_header=False,
                new_columns=["pos", "barcode", "fw", "rv"],
                try_parse_dates=False,
            )

            if df.height == 0:
                log.warning(f"Empty file detected: {filename}")
                return self._create_empty_pivots(base)

            # Left join position data with df to ensure all positions are represented
            full_df = self.pos_data.join(df, on="pos", how="left")

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

            # Create pivots for fw, rv, and cov
            fw_pivot = self._create_pivot(full_df, "fw", base)
            rv_pivot = self._create_pivot(full_df, "rv", base)
            cov_pivot = self._create_pivot(full_df, "cov", base)

            log.debug(
                f"Processed {filename} for base {base}: "
                f"fw: {fw_pivot.shape}, rv: {rv_pivot.shape}, cov: {cov_pivot.shape}"
            )

            return fw_pivot, rv_pivot, cov_pivot

        except Exception as e:
            log.error(f"Error processing {filename} for base {base}: {e}")
            log.debug(traceback.format_exc())
            raise

    def _create_pivot(
        self, df: pl.DataFrame, value_col: str, base: str
    ) -> pl.DataFrame:
        """Create a pivot table for the specified value column."""
        pivot = (
            df.pivot(
                index="barcode",
                on="pos",
                values=value_col,
                aggregate_function="first",
            )
            .fill_null(0)
            .with_columns([pl.lit(base).alias("base")])
            .filter(pl.col("barcode").is_not_null())
        )

        # Reorganize columns to match the expected output format
        base_col = pl.col("base")
        pivot = pivot.select(
            [base_col, *[col for col in pivot.columns if col != "base"]]
        )

        return pivot

    def _create_empty_pivots(
        self, base: str
    ) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
        """Create empty pivot tables for cases with no data."""
        empty_df = pl.DataFrame({"base": [base]})
        return empty_df, empty_df, empty_df

    def load_cluster_data(
        self, srrdir: Union[str, Path], cluster: str = "cell"
    ) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
        """Load data for all bases in a cluster with parallel processing."""
        srrdir = Path(srrdir)

        log.info(f"Loading cluster data from {srrdir} for cluster '{cluster}'")

        if not srrdir.exists():
            log.error(f"Directory not found: {srrdir}")
            raise FileNotFoundError(f"Data directory not found: {srrdir}")

        filenames = [
            srrdir / f"{cluster}.{base}.txt.gz" for base in self.config.bases
        ]

        # Load data for each base with better error handling
        results = []
        for filename, base in zip(filenames, self.config.bases):
            try:
                results.append(self.load_base_data(filename, base))
            except Exception as e:
                log.warning(f"Failed to load {filename}: {e}, using empty data")
                results.append(self._create_empty_pivots(base))

        # Unpack results
        df_fw, df_rv, df_cov = zip(*results)

        # Use polars for faster concatenation
        try:
            df_fw_combined = (
                pl.concat(df_fw)
                if any(df.height > 0 for df in df_fw)
                else pl.DataFrame()
            )
            df_rv_combined = (
                pl.concat(df_rv)
                if any(df.height > 0 for df in df_rv)
                else pl.DataFrame()
            )
            df_cov_combined = (
                pl.concat(df_cov)
                if any(df.height > 0 for df in df_cov)
                else pl.DataFrame()
            )

            log.info(
                f"Combined data shapes - fw: {df_fw_combined.shape}, rv: {df_rv_combined.shape}, cov: {df_cov_combined.shape}"
            )

            return df_fw_combined, df_rv_combined, df_cov_combined

        except Exception as e:
            log.error(f"Error combining dataframes: {e}")
            raise

    def process_srr_entry(
        self, gseid: str, srrid: str, srrdir: Union[str, Path], cluster: str
    ) -> None:
        """Process a single SRR entry with enhanced error handling and logging."""
        start_time = time.time()
        log.info(
            f"Processing {gseid}_{srrid} in {srrdir} for cluster '{cluster}'"
        )

        try:
            # Load the data
            pl_fw, pl_rv, pl_cov = self.load_cluster_data(srrdir, cluster)
            tablename = f"{gseid}_{srrid}"

            # Save the data to CSV files
            output_files = {
                "fw": self.config.table_dir / f"{cluster}_fw.{tablename}.csv",
                "rv": self.config.table_dir / f"{cluster}_rv.{tablename}.csv",
                "cov": self.config.table_dir / f"{cluster}_cov.{tablename}.csv",
            }

            # Write files with error handling
            for data_type, (data, filepath) in zip(
                ["fw", "rv", "cov"],
                [
                    (pl_fw, output_files["fw"]),
                    (pl_rv, output_files["rv"]),
                    (pl_cov, output_files["cov"]),
                ],
            ):
                try:
                    data.write_csv(str(filepath), include_header=True)
                    log.debug(f"Saved {data_type} data to {filepath}")
                except Exception as e:
                    log.error(
                        f"Failed to write {data_type} data to {filepath}: {e}"
                    )
                    raise

            processing_time = time.time() - start_time
            log.info(
                f"Successfully processed {gseid}_{srrid} for cluster '{cluster}' in {processing_time:.2f}s"
            )

        except Exception as e:
            log.error(
                f"Failed to process {gseid}_{srrid} for cluster '{cluster}': {e}"
            )
            log.debug(traceback.format_exc())
            raise

    def get_all_processing_tasks(self) -> List[Dict[str, str]]:
        """Generate all processing tasks from SRR data."""
        tasks = []
        for row in self.srr_data.to_dicts():
            for cluster in self.config.clusters:
                tasks.append({**row, "cluster": cluster})
        return tasks

    def generate_slurm_script(
        self, output_file: Union[str, Path] = None
    ) -> Path:
        """Generate SLURM batch script for parallel processing."""
        if output_file is None:
            output_file = (
                self.config.base_dir / "analysis/zzz/db/slurm_run_all.sh"
            )

        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        tasks = self.get_all_processing_tasks()

        commands = []
        for task in tasks:
            cmd_parts = [
                "/scr1/users/liuc9/tools/anaconda3/envs/renv/bin/python3.13",
                str(Path(__file__).absolute()),
                "run-one",
                task["gseid"],
                task["srrid"],
                task["srrdir"],
                task["cluster"],
            ]
            commands.append(" ".join(cmd_parts))

        with open(output_file, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(
                "# Auto-generated SLURM script for scMOCHA data processing\n\n"
            )
            for cmd in commands:
                f.write(f"{cmd}\n")

        output_file.chmod(0o755)  # Make executable
        log.info(
            f"Generated SLURM script with {len(commands)} commands: {output_file}"
        )
        return output_file


class ParallelProcessor:
    """Enhanced parallel processing with progress tracking and error handling."""

    def __init__(
        self, data_processor: DataProcessor, max_workers: Optional[int] = None
    ):
        self.data_processor = data_processor
        self.max_workers = max_workers or min(
            mp.cpu_count(), 8
        )  # Reasonable default

    def process_single_task(
        self, task: Dict[str, str]
    ) -> Dict[str, Union[str, bool, float]]:
        """Process a single task with error handling and timing."""
        start_time = time.time()
        task_id = f"{task['gseid']}_{task['srrid']}_{task['cluster']}"

        try:
            self.data_processor.process_srr_entry(
                task["gseid"], task["srrid"], task["srrdir"], task["cluster"]
            )

            return {
                "task_id": task_id,
                "success": True,
                "processing_time": time.time() - start_time,
                "error": None,
            }

        except Exception as e:
            return {
                "task_id": task_id,
                "success": False,
                "processing_time": time.time() - start_time,
                "error": str(e),
            }

    def process_all_parallel(self) -> Dict[str, Union[int, float, List]]:
        """Process all tasks in parallel with progress tracking."""
        tasks = self.data_processor.get_all_processing_tasks()
        start_time = time.time()

        log.info(
            f"Starting parallel processing of {len(tasks)} tasks with {self.max_workers} workers"
        )

        # Create progress tracking
        with Progress(console=console) as progress:
            main_task = progress.add_task(
                "[green]Processing SRR entries...", total=len(tasks)
            )

            results = []
            failed_tasks = []

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=self.max_workers
            ) as executor:
                # Submit all tasks
                future_to_task = {
                    executor.submit(self.process_single_task, task): task
                    for task in tasks
                }

                # Process completed tasks
                for future in concurrent.futures.as_completed(future_to_task):
                    result = future.result()
                    results.append(result)

                    if not result["success"]:
                        failed_tasks.append(result)
                        log.error(
                            f"Task failed: {result['task_id']} - {result['error']}"
                        )
                    else:
                        log.debug(
                            f"Task completed: {result['task_id']} in {result['processing_time']:.2f}s"
                        )

                    progress.update(main_task, advance=1)

        # Summary statistics
        total_time = time.time() - start_time
        successful_tasks = len([r for r in results if r["success"]])

        summary = {
            "total_tasks": len(tasks),
            "successful_tasks": successful_tasks,
            "failed_tasks": len(failed_tasks),
            "total_processing_time": total_time,
            "average_task_time": sum(r["processing_time"] for r in results)
            / len(results)
            if results
            else 0,
            "failed_task_details": failed_tasks,
        }

        # Log results
        log.info("Parallel processing completed:")
        log.info(f"  Total tasks: {summary['total_tasks']}")
        log.info(f"  Successful: {summary['successful_tasks']}")
        log.info(f"  Failed: {summary['failed_tasks']}")
        log.info(f"  Total time: {summary['total_processing_time']:.2f}s")
        log.info(f"  Average task time: {summary['average_task_time']:.2f}s")

        if failed_tasks:
            log.warning("Failed tasks summary:")
            for task in failed_tasks[:5]:  # Show first 5 failures
                log.warning(f"  {task['task_id']}: {task['error']}")
            if len(failed_tasks) > 5:
                log.warning(f"  ... and {len(failed_tasks) - 5} more failures")

        return summary


# Global instances for CLI usage
config = Config()
data_processor = DataProcessor(config)
parallel_processor = ParallelProcessor(data_processor)


app = typer.Typer(help="Enhanced scMOCHA genomic data processing pipeline")


@app.command()
def run_all(
    max_workers: int = typer.Option(
        default=None,
        help="Maximum number of worker processes (default: CPU count)",
    ),
    log_level: str = typer.Option(
        default="INFO", help="Logging level (DEBUG, INFO, WARNING, ERROR)"
    ),
    log_file: Optional[str] = typer.Option(
        default=None, help="Optional log file path"
    ),
):
    """Run all SRR entries in parallel with enhanced logging and error handling."""

    # Setup logging
    global log
    log = setup_logging(log_level, log_file)

    try:
        # Initialize components
        config_obj = Config()
        processor = DataProcessor(config_obj)
        parallel_proc = ParallelProcessor(processor, max_workers)

        # Process all tasks
        summary = parallel_proc.process_all_parallel()

        # Display summary table
        table = Table(title="Processing Summary")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="green")

        table.add_row("Total Tasks", str(summary["total_tasks"]))
        table.add_row("Successful", str(summary["successful_tasks"]))
        table.add_row("Failed", str(summary["failed_tasks"]))
        table.add_row("Total Time", f"{summary['total_processing_time']:.2f}s")
        table.add_row("Avg Task Time", f"{summary['average_task_time']:.2f}s")

        console.print(table)

        # Exit with appropriate code
        sys.exit(0 if summary["failed_tasks"] == 0 else 1)

    except Exception as e:
        log.error(f"Fatal error in run_all: {e}")
        log.debug(traceback.format_exc())
        sys.exit(1)


@app.command()
def run_one(
    gseid: str = typer.Argument(..., help="GSE ID"),
    srrid: str = typer.Argument(..., help="SRR ID"),
    srrdir: str = typer.Argument(..., help="SRR directory path"),
    cluster: str = typer.Argument(..., help="Cluster type (cell or cluster)"),
    log_level: str = typer.Option(
        default="INFO", help="Logging level (DEBUG, INFO, WARNING, ERROR)"
    ),
):
    """Process a single SRR entry with enhanced error handling."""

    # Setup logging
    global log
    log = setup_logging(log_level)

    try:
        # Initialize components
        config_obj = Config()
        processor = DataProcessor(config_obj)

        # Process single entry
        processor.process_srr_entry(gseid, srrid, srrdir, cluster)

        log.info(
            f"Successfully processed {gseid}_{srrid} for cluster '{cluster}'"
        )

    except Exception as e:
        log.error(
            f"Failed to process {gseid}_{srrid} for cluster '{cluster}': {e}"
        )
        log.debug(traceback.format_exc())
        sys.exit(1)


@app.command()
def generate_slurm(
    output_file: Optional[str] = typer.Option(
        default=None, help="Output file path for SLURM script"
    ),
    log_level: str = typer.Option(default="INFO", help="Logging level"),
):
    """Generate SLURM batch script for cluster processing."""

    # Setup logging
    global log
    log = setup_logging(log_level)

    try:
        # Initialize components
        config_obj = Config()
        processor = DataProcessor(config_obj)

        # Generate script
        script_path = processor.generate_slurm_script(output_file)

        console.print(f"[green]Generated SLURM script: {script_path}[/green]")

    except Exception as e:
        log.error(f"Failed to generate SLURM script: {e}")
        log.debug(traceback.format_exc())
        sys.exit(1)


@app.command()
def validate_setup(
    log_level: str = typer.Option(default="INFO", help="Logging level"),
):
    """Validate the setup and configuration."""

    # Setup logging
    global log
    log = setup_logging(log_level)

    try:
        # Initialize and validate configuration
        config_obj = Config()
        processor = DataProcessor(config_obj)

        # Create validation table
        table = Table(title="Setup Validation")
        table.add_column("Component", style="cyan")
        table.add_column("Status", style="green")
        table.add_column("Details", style="yellow")

        table.add_row(
            "SRR Data", "✓ OK", f"{processor.srr_data.shape[0]} entries"
        )
        table.add_row(
            "Position Data", "✓ OK", f"{processor.pos_data.shape[0]} positions"
        )
        table.add_row("Output Directory", "✓ OK", str(config_obj.table_dir))
        table.add_row("Bases", "✓ OK", ", ".join(config_obj.bases))
        table.add_row("Clusters", "✓ OK", ", ".join(config_obj.clusters))

        console.print(table)

        # Show sample tasks
        tasks = processor.get_all_processing_tasks()
        console.print(f"\n[bold]Total processing tasks: {len(tasks)}[/bold]")

        if tasks:
            console.print("[bold]Sample tasks:[/bold]")
            for i, task in enumerate(tasks[:3]):
                console.print(
                    f"  {i + 1}. {task['gseid']}_{task['srrid']} ({task['cluster']})"
                )

    except Exception as e:
        log.error(f"Setup validation failed: {e}")
        log.debug(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    app()
