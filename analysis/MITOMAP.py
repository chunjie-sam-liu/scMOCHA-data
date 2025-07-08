#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-05 18:53:48
# @DESCRIPTION:
# @VERSION: v0.0.1


import sqlite3

import polars as pl

FILENAME = "/mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3"


def connect_to_db():
    """Connect to the SQLite database and return the connection."""
    try:
        conn = sqlite3.connect(FILENAME)
        print(f"Successfully connected to {FILENAME}")
        return conn
    except sqlite3.Error as e:
        print(f"Error connecting to database: {e}")
        return None


def show_tables(conn):
    """Show all tables in the connected database."""
    if conn is None:
        return

    try:
        cursor = conn.cursor()
        # Query to get all table names from sqlite_master
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()

        if tables:
            print("Tables in the database:")
            for table in tables:
                print(f"- {table[0]}")
        else:
            print("No tables found in the database.")
    except sqlite3.Error as e:
        print(f"Error listing tables: {e}")


def show_table_info(conn, table_name):
    """Show the schema of a specific table."""
    if conn is None:
        return

    try:
        cursor = conn.cursor()
        # Query to get the schema of the specified table
        cursor.execute(f"PRAGMA table_info({table_name});")
        columns = cursor.fetchall()

        if columns:
            print(f"Schema for table '{table_name}':")
            for column in columns:
                print(f"- {column[1]}: {column[2]}")
        else:
            print(f"No information found for table '{table_name}'.")
    except sqlite3.Error as e:
        print(f"Error retrieving table info: {e}")


def extract_table_data(conn, table_name):
    """Extract data from a specific table and return it as a polars DataFrame."""
    if conn is None:
        print("No database connection available.")
        return None

    try:
        # Get column names from table info
        cursor = conn.cursor()
        cursor.execute(f"PRAGMA table_info({table_name});")
        columns = [column[1] for column in cursor.fetchall()]

        # Query all data from the table
        cursor.execute(f"SELECT * FROM {table_name};")
        rows = cursor.fetchall()

        # Create a polars DataFrame
        if rows:
            df = pl.DataFrame(rows, schema=columns)
            print(f"Successfully extracted {len(rows)} rows from {table_name}")
            print(f"DataFrame shape: {df.shape}")
            print(f"DataFrame schema: {df.schema}")
            return df
        else:
            print(f"No data found in table '{table_name}'.")
            return pl.DataFrame(schema=columns)
    except ImportError:
        print(
            "Error: polars library not installed. Install with 'pip install polars'"
        )
        return None
    except sqlite3.Error as e:
        print(f"Error extracting data from table: {e}")
        return None


if __name__ == "__main__":
    conn = connect_to_db()
    if conn:
        show_tables(conn)
        show_table_info(conn, "gnomad")
        show_table_info(conn, "conservation_rate")
        conservation_rate = extract_table_data(conn, "conservation_rate")
        conservation_rate.write_csv(
            "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/conservation_rate.csv"
        )
        extract_table_data(conn, "refseq")
        gnomad = extract_table_data(conn, "gnomad")
        gnomad.write_csv(
            "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/gnomad.csv"
        )

        conn.close()
