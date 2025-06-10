#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-04-23 23:58:00
# @DESCRIPTION:
# @VERSION: v0.0.1

import concurrent.futures
import functools
import logging
import time
from pathlib import Path

import duckdb
import pandas as pd

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pd.read_csv(SRR_FILENAME)
SRR

FASTA_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/config/rCRS.MT.fasta.df.csv"
)
FASTA = pd.read_csv(FASTA_FILENAME)
FASTA
POS = FASTA["pos"].to_frame()
POS

BASES = ["A", "C", "G", "T"]
CLUSTERS = ["cell", "cluster"]


DUCKDB_DIR = Path("/home/liuc9/github/scMOCHA-data/analysis/zzz/db")
DUCKDB_VERSION = duckdb.__version__


# Configure logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("duckdb_operations")


def log_operation(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        start_time = time.time()
        logger.info(
            f"Starting {func.__name__} with args: {args}, kwargs: {kwargs}"
        )
        try:
            result = func(self, *args, **kwargs)
            elapsed_time = time.time() - start_time
            logger.info(
                f"Completed {func.__name__} in {elapsed_time:.2f} seconds"
            )
            return result
        except Exception as e:
            logger.error(f"Error in {func.__name__}: {str(e)}")
            raise

    return wrapper


class DuckDB:
    def __init__(self, path, dnname, version):
        self.path = path
        self.dbname = str(dnname)
        self.version = version
        logger.info(
            f"Connecting to database: {path}/{self.dbname}.duckdb.{self.version}"
        )
        self.connection = duckdb.connect(
            database=self.path / f"{self.dbname}.duckdb.{self.version}"
        )

    @log_operation
    def query(self, query):
        return self.connection.execute(query).fetchdf()

    @log_operation
    def create_table_from_df(self, tablename, df):
        self.connection.execute(f'DROP TABLE IF EXISTS "{tablename}"')
        self.connection.execute(
            f'CREATE TABLE "{tablename}" AS SELECT * FROM df'
        )

    @log_operation
    def create_table_from_csv(self, tablename, csvfile):
        print(f"Creating table {tablename} from {csvfile}")
        if not Path(csvfile).exists():
            print(f"File {csvfile} does not exist.")
            return
        self.connection.execute(f'DROP TABLE IF EXISTS "{tablename}"')

        # Check if the CSV file contains a 'chr' column before specifying types
        # Read the first line of the CSV file to get column names
        columns = []
        with open(csvfile, "r") as f:
            header = f.readline().strip()
            columns = header.split(",")

        print(f"Columns in {csvfile}: {columns}")
        if "chr" in columns:
            print("chr column found, specifying types")
        else:
            print("chr column not found, not specifying types")

        if "chr" in columns:
            # Specify column types in read_csv_auto
            self.connection.execute(f"""
                CREATE TABLE "{tablename}" AS
                SELECT * FROM read_csv_auto('{csvfile}', types={{'chr': 'VARCHAR'}})
            """)
        else:
            # Create table without specifying types for 'chr'
            self.connection.execute(f"""
                CREATE TABLE "{tablename}" AS
                SELECT * FROM read_csv_auto('{csvfile}')
            """)

    @log_operation
    def tables(self):
        return self.query("SHOW TABLES")

    @log_operation
    def table_info(self, tablename):
        return self.query(f"DESCRIBE {tablename}")

    @log_operation
    def table_count(self, tablename):
        return self.query(f"SELECT COUNT(*) FROM {tablename}").iloc[0, 0]

    @log_operation
    def table_head(self, tablename):
        return self.query(f"SELECT * FROM {tablename} LIMIT 5")

    @log_operation
    def table_tail(self, tablename):
        return self.query(
            f"SELECT * FROM {tablename} ORDER BY rowid DESC LIMIT 5"
        )

    def close(self):
        logger.info("Closing database connection")
        self.connection.close()


def show_tables(dbname):
    db = DuckDB(DUCKDB_DIR, dbname, DUCKDB_VERSION)
    print(f"DuckDB tables:{db.tables()}")
    db.close()


def show_table_info(dbname, tablename):
    db = DuckDB(DUCKDB_DIR, dbname, DUCKDB_VERSION)
    print(f"Table {tablename} info: {db.table_info(tablename)}")
    print(f"Table {tablename} count: {db.table_count(tablename)}")
    print(f"Table {tablename} head: {db.table_head(tablename)}")
    print(f"Table {tablename} tail: {db.table_tail(tablename)}")
    db.close()


def create_one_table(dbname, tablename, df):
    db = DuckDB(DUCKDB_DIR, dbname, DUCKDB_VERSION)
    db.create_table_from_df(tablename, df)
    db.close()


def load_df_from_base(filename, base):
    df = pd.read_csv(
        filename, sep=",", header=None, names=["pos", "barcode", "fw", "rv"]
    )
    df["cov"] = df["fw"] + df["rv"]

    # Left join POS with df to ensure all positions are represented
    full_df = POS.merge(df, on="pos", how="left")

    # Fill NaN values with 0 for counts
    full_df["fw"] = full_df["fw"].fillna(0).astype(int)
    full_df["rv"] = full_df["rv"].fillna(0).astype(int)
    full_df["cov"] = full_df["cov"].fillna(0).astype(int)

    # Create fw pivot df
    fw_pivot = (
        full_df[["barcode", "pos", "fw"]]
        .pivot(index="barcode", columns="pos", values="fw")
        .fillna(0)
        .astype(int)
    )
    fw_pivot = fw_pivot[fw_pivot.index.notna()]
    fw_pivot = fw_pivot.reset_index()
    # Create rv pivot df
    rv_pivot = (
        full_df[["barcode", "pos", "rv"]]
        .pivot(index="barcode", columns="pos", values="rv")
        .fillna(0)
        .astype(int)
    )
    rv_pivot = rv_pivot[rv_pivot.index.notna()]
    rv_pivot = rv_pivot.reset_index()
    # Create cov pivot df
    cov_pivot = (
        full_df[["barcode", "pos", "cov"]]
        .pivot(index="barcode", columns="pos", values="cov")
        .fillna(0)
        .astype(int)
    )
    cov_pivot = cov_pivot[cov_pivot.index.notna()]
    cov_pivot = cov_pivot.reset_index()
    # Create db
    # Add base information to pivot tables
    fw_pivot["base"] = base
    rv_pivot["base"] = base
    cov_pivot["base"] = base

    # For better organization, move the base column to the beginning
    fw_pivot_cols = ["base", "barcode"] + [
        col for col in fw_pivot.columns if col not in ["base", "barcode"]
    ]
    rv_pivot_cols = ["base", "barcode"] + [
        col for col in rv_pivot.columns if col not in ["base", "barcode"]
    ]
    cov_pivot_cols = ["base", "barcode"] + [
        col for col in cov_pivot.columns if col not in ["base", "barcode"]
    ]

    fw_pivot = fw_pivot[fw_pivot_cols]
    rv_pivot = rv_pivot[rv_pivot_cols]
    cov_pivot = cov_pivot[cov_pivot_cols]

    return fw_pivot, rv_pivot, cov_pivot


def load_df(srrdir, cluster="cell"):
    filenames = [f"{srrdir}/{cluster}.{base}.txt.gz" for base in BASES]
    df_bases = [
        load_df_from_base(filename, base)
        for filename, base in zip(filenames, BASES)
    ]
    df_fw, df_rv, df_cov = zip(*df_bases)
    df_fw = pd.concat(df_fw, ignore_index=True)
    df_rv = pd.concat(df_rv, ignore_index=True)
    df_cov = pd.concat(df_cov, ignore_index=True)
    return df_fw, df_rv, df_cov


def srrid2table(gseid, srrid, srrdir, cluster="cell"):
    # Load the data
    df_fw, df_rv, df_cov = load_df(srrdir, cluster)
    tablename = f"{gseid}_{srrid}"

    create_one_table(f"{cluster}_fw", tablename, df_fw)
    create_one_table(f"{cluster}_rw", tablename, df_rv)
    create_one_table(f"{cluster}_cov", tablename, df_cov)


def load_all_into_db(cluster="cell"):
    # Load the data
    for _, row in SRR.iterrows():
        gseid = row["gseid"]
        srrid = row["srrid"]
        srrdir = row["srrdir"]
        print(f"Loading {srrdir}")
        srrid2table(gseid, srrid, srrdir, cluster)


def load_all_clusters_in_parallel():
    """Load all data for both cell and cluster types in parallel."""
    logger.info("Starting parallel loading of all clusters")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Submit both cluster types for parallel processing
        futures = [
            executor.submit(load_all_into_db, cluster) for cluster in CLUSTERS
        ]

        # Process results as they complete
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logger.error(f"Error in parallel execution: {str(e)}")

    logger.info("Completed loading all clusters")


if __name__ == "__main__":
    # Execute the parallel loading function
    load_all_clusters_in_parallel()
    # load_all_into_db(cluster="cell")
    # load_all_into_db(cluster="cluster")

    # srrid2table(
    #     "GSE155673",
    #     "SRR11512399",
    #     "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE155673/final/GSM4712885",
    #     cluster="cell",
    # )
    # create_one_table(
    #     "/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE155673/final/GSM4712885/cell.C.txt.gz"
    # )
    # show_tables("cell_cov")
    # show_table_info("cell_fw", "GSE155673_SRR11512399")
