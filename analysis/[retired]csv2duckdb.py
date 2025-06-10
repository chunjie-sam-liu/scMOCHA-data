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
import polars as pl
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)

logger = logging.getLogger("rich")
logger.info("Hello, World!")


SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)
SRR
CLUSTERS = ["cell", "cluster"]
TABLEDIR = Path("/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/TABLES")

DUCKDB_DIR = Path("/home/liuc9/github/scMOCHA-data/analysis/zzz/db")
DUCKDB_VERSION = duckdb.__version__


def log_operation(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        start_time = time.time()
        logger.info(f"Starting {func.__name__}")
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

        logger.info(f"Columns in {csvfile}: {columns}")

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


def create_one_table_from_df(dbname, tablename, df):
    db = DuckDB(DUCKDB_DIR, dbname, DUCKDB_VERSION)
    db.create_table_from_df(tablename, df)
    db.close()


def create_one_table_from_csv(dbname, tablename, csvfile):
    db = DuckDB(DUCKDB_DIR, dbname, DUCKDB_VERSION)
    db.create_table_from_csv(tablename, csvfile)
    db.close()


if __name__ == "__main__":
    create_one_table_from_csv(
        dbname="cell_cov",
        tablename="GSE279945_GSM8583968",
        csvfile="/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/TABLES/cell_cov.GSE279945_GSM8583968.csv",
    )
    df = pl.read_csv(
        "/mnt/isilon/u01_project/large-scale/liuc9/raw/zzz/db/TABLES/cell_cov.GSE279945_GSM8583968.csv"
    )
    create_one_table_from_df("cell_cov", "GSE279945_GSM8583968", df)

    show_tables("cell_cov")
    show_table_info("cell_cov", "GSE279945_GSM8583968")
