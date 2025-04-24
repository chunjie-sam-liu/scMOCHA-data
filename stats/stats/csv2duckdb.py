#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-24 12:35:26
# @DESCRIPTION:
# @VERSION: v0.0.1

import concurrent.futures
import glob
import logging
import os

import duckdb
import polars as pl

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def import_csv_to_duckdb():
    """Import CSV files to DuckDB database."""
    # Path to CSV files
    csv_dir = "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/db/TABLES"
    # Output DuckDB file
    duckdb_path = (
        "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/db/cell_cov.duckdb"
    )

    # Find all cell_cov CSV files
    csv_files = glob.glob(os.path.join(csv_dir, "cell_cov.*.csv"))
    logging.info(f"Found {len(csv_files)} CSV files to import")

    def process_csv_file(csv_file):
        """Process a single CSV file and import it to DuckDB."""
        # Extract the table name from the filename
        filename = os.path.basename(csv_file)
        table_name = filename.replace("cell_cov.", "").replace(".csv", "")

        try:
            logging.info(f"Loading {filename} into memory")
            df = pl.read_csv(csv_file)

            # Open a connection for each file to avoid contention
            with duckdb.connect(duckdb_path) as conn:
                conn.register("temp_df", df)
                conn.execute(
                    f"CREATE TABLE IF NOT EXISTS {table_name} AS SELECT * FROM temp_df"
                )
                logging.info(f"Importing {filename} to table {table_name}")

            logging.info(f"Successfully imported {table_name}")
            return table_name
        except Exception as e:
            logging.error(f"Failed to import {filename}: {str(e)}")
            return None

    # Use ThreadPoolExecutor to process files in parallel
    max_workers = 10  # Limit to 8 or CPU count, whichever is smaller
    logging.info(f"Starting parallel import with {max_workers} workers")

    with concurrent.futures.ThreadPoolExecutor(
        max_workers=max_workers
    ) as executor:
        results = list(executor.map(process_csv_file, csv_files))

    successful_imports = [r for r in results if r]
    logging.info(
        f"Import process completed. Successfully imported {len(successful_imports)} tables."
    )


if __name__ == "__main__":
    import_csv_to_duckdb()
