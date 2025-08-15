#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Author: Chunjie Liu
# Contact: chunjie.sam.liu.at.gmail.com
# Date: 2025-08-15
# Description: Efficiently merge 577 10x h5 files into one AnnData object with memory management
# Version: 0.2

import argparse
import gc
import logging
import os
import sys
from pathlib import Path
from typing import Optional

# Check if required packages are available
try:
    import anndata as ad
    import polars as pl
    import scanpy as sc
    from rich.console import Console
    from rich.logging import RichHandler
    from rich.progress import (
        BarColumn,
        Progress,
        SpinnerColumn,
        TaskProgressColumn,
        TextColumn,
    )
    from rich.traceback import install
except ImportError as e:
    print(f"Error importing required packages: {e}")
    print(
        "Please ensure all required packages are installed in your conda environment."
    )
    sys.exit(1)

# Install rich traceback handler
install()

# Setup rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console)],
)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_OUTDIR = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/EXPRMERGE"
)
DEFAULT_SRR_FILE = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)


def load_h5_safe(srrdir: Path, gseid: str, srrid: str) -> Optional[ad.AnnData]:
    """
    Safely load h5 file and return AnnData object with sample metadata.

    Args:
        srrdir: Path to the directory containing the h5 file
        gseid: GSE ID for the sample
        srrid: SRR ID for the sample

    Returns:
        AnnData object or None if loading failed
    """
    h5_filename = srrdir / "filtered_feature_bc_matrix.h5"
    sampleid = f"{gseid}_{srrid}"

    try:
        if not h5_filename.exists():
            logger.warning(f"H5 file not found: {h5_filename}")
            return None

        adata = sc.read_10x_h5(h5_filename)
        adata.var_names_make_unique()

        # Add sample metadata
        adata.obs["gseid"] = gseid
        adata.obs["srrid"] = srrid
        adata.obs["sampleid"] = sampleid

        # Make observation names unique by adding sample prefix
        adata.obs_names = [
            f"{sampleid}_{barcode}" for barcode in adata.obs_names
        ]

        logger.info(
            f"✓ Loaded {sampleid}: {adata.n_obs} cells, {adata.n_vars} genes"
        )
        return adata

    except Exception as e:
        logger.error(f"✗ Failed to load {sampleid}: {e}")
        return None


def process_batch(batch_data: list, batch_idx: int) -> Optional[ad.AnnData]:
    """
    Process a batch of samples and return concatenated AnnData object.

    Args:
        batch_data: List of (gseid, srrid, srrdir) tuples
        batch_idx: Batch index for logging

    Returns:
        Concatenated AnnData object or None if all samples failed
    """
    logger.info(
        f"Processing batch {batch_idx + 1} with {len(batch_data)} samples"
    )

    adatas = {}
    successful_loads = 0

    for gseid, srrid, srrdir in batch_data:
        srrdir = Path(srrdir)
        sampleid = f"{gseid}_{srrid}"

        adata = load_h5_safe(srrdir, gseid, srrid)
        if adata is not None:
            adatas[sampleid] = adata
            successful_loads += 1

    if not adatas:
        logger.warning(
            f"No samples successfully loaded in batch {batch_idx + 1}"
        )
        return None

    logger.info(
        f"Successfully loaded {successful_loads}/{len(batch_data)} samples in batch {batch_idx + 1}"
    )

    # Concatenate samples in this batch
    batch_adata = ad.concat(adatas, label="sample", keys=None)

    # Clear memory
    del adatas
    gc.collect()

    return batch_adata


def calculate_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """Calculate quality control metrics for the AnnData object."""
    logger.info("Calculating QC metrics...")

    # Gene annotations
    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # mitochondrial genes
    adata.var["ribo"] = adata.var_names.str.startswith(
        ("RPS", "RPL")
    )  # ribosomal genes
    adata.var["hb"] = adata.var_names.str.contains(
        "^HB[^(P)]"
    )  # hemoglobin genes

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    logger.info(
        f"QC metrics calculated: {adata.n_obs} cells, {adata.n_vars} genes"
    )
    return adata


def filter_and_normalize(adata: ad.AnnData) -> ad.AnnData:
    """Apply filtering and normalization to the AnnData object."""
    logger.info("Applying filtering and normalization...")

    # Store original counts
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)

    logger.info(
        f"Filtered: {n_cells_before} → {adata.n_obs} cells, {n_genes_before} → {adata.n_vars} genes"
    )

    # Save count data
    adata.layers["counts"] = adata.X.copy()

    # Normalize to median total counts
    sc.pp.normalize_total(adata)

    # Log-transform
    sc.pp.log1p(adata)

    logger.info("Normalization completed")
    return adata


def merge_h5_files(
    input_csv: Path,
    output_dir: Path,
    batch_size: int,
    max_samples: Optional[int],
) -> int:
    """
    Merge multiple 10x h5 files into one AnnData object efficiently.

    Returns:
        0 for success, 1 for failure
    """
    console.print(
        "[bold green]🧬 Starting scMOCHA data merging process[/bold green]"
    )

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    # Load sample information
    if not input_csv.exists():
        logger.error(f"Input CSV file not found: {input_csv}")
        return 1

    srr_data = pl.read_csv(input_csv)
    total_samples = len(srr_data)

    if max_samples:
        srr_data = srr_data.head(max_samples)
        total_samples = len(srr_data)
        logger.info(f"Limited to {max_samples} samples for testing")

    logger.info(f"Found {total_samples} samples to process")

    # Convert to list of tuples for batch processing
    sample_list = (
        srr_data.select(["gseid", "srrid", "srrdir"]).to_numpy().tolist()
    )

    # Process samples in batches
    batch_adatas = []
    num_batches = (total_samples + batch_size - 1) // batch_size

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(
            f"Processing {num_batches} batches", total=num_batches
        )

        for i in range(0, total_samples, batch_size):
            batch_data = sample_list[i : i + batch_size]
            batch_idx = i // batch_size

            batch_adata = process_batch(batch_data, batch_idx)
            if batch_adata is not None:
                batch_adatas.append(batch_adata)

            progress.update(task, advance=1)

            # Force garbage collection between batches
            gc.collect()

    if not batch_adatas:
        logger.error("No samples were successfully loaded!")
        return 1

    logger.info(f"Concatenating {len(batch_adatas)} batches...")

    # Concatenate all batches
    adata_raw = ad.concat(batch_adatas, label="batch", keys=None)
    adata_raw.obs_names_make_unique()

    # Clear batch data from memory
    del batch_adatas
    gc.collect()

    logger.info(
        f"Final raw data: {adata_raw.n_obs} cells, {adata_raw.n_vars} genes"
    )

    # Calculate QC metrics
    adata_raw = calculate_qc_metrics(adata_raw)

    # Save raw data
    raw_file = output_dir / "merged_raw_data.h5ad"
    logger.info(f"Saving raw data to: {raw_file}")
    adata_raw.write(raw_file)

    # Create normalized version
    adata_norm = adata_raw.copy()
    adata_norm = filter_and_normalize(adata_norm)

    # Save normalized data
    norm_file = output_dir / "merged_normalized_data.h5ad"
    logger.info(f"Saving normalized data to: {norm_file}")
    adata_norm.write(norm_file)

    # Print summary
    console.print(
        "\n[bold green]✅ Processing completed successfully![/bold green]"
    )
    console.print(
        f"[bold]Raw data:[/bold] {adata_raw.n_obs:,} cells × {adata_raw.n_vars:,} genes"
    )
    console.print(
        f"[bold]Normalized data:[/bold] {adata_norm.n_obs:,} cells × {adata_norm.n_vars:,} genes"
    )
    console.print(f"[bold]Files saved to:[/bold] {output_dir}")

    # Display memory usage info
    try:
        import psutil

        process = psutil.Process(os.getpid())
        memory_gb = process.memory_info().rss / 1024 / 1024 / 1024
        console.print(f"[dim]Peak memory usage: {memory_gb:.2f} GB[/dim]")
    except ImportError:
        logger.info("psutil not available for memory monitoring")

    return 0  # Success


def main():
    """Main CLI interface using argparse."""
    parser = argparse.ArgumentParser(
        description="""
Merge multiple 10x h5 files into one AnnData object efficiently.

This script processes 577 10x Genomics h5 files in batches to manage memory usage,
then saves both raw and normalized versions of the merged data.

Examples:
  # Test with 10 samples
  python %(prog)s --max-samples 10 --batch-size 5

  # Process all samples with default settings
  python %(prog)s

  # Custom batch size for high-memory machines
  python %(prog)s --batch-size 30
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--input-csv",
        "-i",
        type=Path,
        default=DEFAULT_SRR_FILE,
        help="CSV file containing GSE ID, SRR ID, and directory paths",
    )

    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory for merged files",
    )

    parser.add_argument(
        "--batch-size",
        "-b",
        type=int,
        default=25,
        help="Number of samples to process in each batch (default: 25)",
    )

    parser.add_argument(
        "--max-samples",
        "-m",
        type=int,
        default=None,
        help="Maximum number of samples to process (for testing)",
    )

    args = parser.parse_args()

    # Validate arguments
    if args.batch_size <= 0:
        parser.error("Batch size must be positive")

    if args.max_samples is not None and args.max_samples <= 0:
        parser.error("Max samples must be positive")

    exit_code = merge_h5_files(
        input_csv=args.input_csv,
        output_dir=args.output_dir,
        batch_size=args.batch_size,
        max_samples=args.max_samples,
    )

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
