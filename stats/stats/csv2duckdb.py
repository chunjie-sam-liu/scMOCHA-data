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

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def import_csv_to_duckdb():
    """Import CSV files to DuckDB database."""
    # Path to CSV files
    csv_dir = "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/db/TABLES"
    # Output DuckDB file
    duckdb_path = "/home/liuc9/github/scMOCHA-data/stats/stats/cell_cov.duckdb"

    # Connect to DuckDB with optimized settings
    logging.info(f"Connecting to DuckDB at {duckdb_path}")
    conn = duckdb.connect(duckdb_path)

    # Configure DuckDB for better performance
    conn.execute("PRAGMA threads=4")  # Adjust based on your CPU cores
    conn.execute("PRAGMA memory_limit='4GB'")  # Adjust based on available RAM

    # Find all cell_cov CSV files
    csv_files = glob.glob(os.path.join(csv_dir, "cell_cov.*.csv"))
    logging.info(f"Found {len(csv_files)} CSV files to import")

    # Create all tables in a single transaction
    conn.execute("BEGIN TRANSACTION")

    try:
        # Batch process the files
        for csv_file in csv_files:
            # Extract the table name from the filename
            filename = os.path.basename(csv_file)
            table_name = filename.replace("cell_cov.", "").replace(".csv", "")

            try:
                logging.info(f"Importing {filename} to table {table_name}")

                # Use COPY statement for faster imports
                conn.execute(
                    f"CREATE TABLE IF NOT EXISTS {table_name} AS SELECT * FROM read_csv_auto('{csv_file}', parallel=True)"
                )

                logging.info(f"Successfully imported {table_name}")
            except Exception as e:
                logging.error(f"Failed to import {filename}: {str(e)}")

        # Commit the transaction
        conn.execute("COMMIT")
        logging.info("All imports committed successfully")

    except Exception as e:
        # Roll back in case of error
        conn.execute("ROLLBACK")
        logging.error(f"Transaction rolled back due to error: {str(e)}")

    # Close connection
    conn.close()
    logging.info("Import process completed")


if __name__ == "__main__":
    import_csv_to_duckdb()
