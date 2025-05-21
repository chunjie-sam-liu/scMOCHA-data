#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-16 14:40:50
# @DESCRIPTION:
# @VERSION: v0.0.1


import concurrent.futures
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
import typer
from rich import print
from tqdm import tqdm

OUTDIR = Path("/home/liuc9/github/scMOCHA-data/stats/stats/zzz/db/EXPR")
SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)

BARCODE_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/barcode_celltype.feather"
)

# Fix: Use BARCODE_FILENAME instead of SRR_FILENAME with read_ipc
BARCODE = pl.read_ipc(BARCODE_FILENAME)
celltypes = ["B", "CD4_T", "CD8_T", "DC", "Mono", "NK", "other", "other_T"]


"""
Prompt for Copilot to complete the function.
- filtered_feature_bc_matrix.h5 is cellranger output,
- BARCODE column are gseid, srrid, barcode, celltype

self.barcode = BARCODE.filter(pl.col("srrid") == self.srrid)
self.barcode matched with self.adata gene name

update the function celltype_expr to process adata, result in each gene mean expression for each cell type, the columns are: genename, *celltypes;

return polars dataframe
"""


class SCEXPR:
    """
    Single-cell expression (SCE) class for analyzing single-cell RNA-seq data.
    """

    def __init__(self, gseid: str, srrid: str, srrdir: Path):
        self.gseid = gseid
        self.srrid = srrid
        self.srrdir = srrdir
        self.barcode = BARCODE.filter(pl.col("srrid") == self.srrid)
        self.h5_file = Path(self.srrdir) / "filtered_feature_bc_matrix.h5"
        self.adata = sc.read_10x_h5(self.h5_file)
        self.adata.var_names_make_unique()

    def __repr__(self):
        return f"SCE(adata={self.adata})"

    def celltype_expr(
        self,
        normalize: bool = True,
        norm_method: str = "both",
        scale: bool = False,
    ):
        """
        For this sample, compute normalized mean expression per gene for each cell type.

        Parameters:
            normalize (bool): Whether to normalize data before computing means
            norm_method (str): Method for normalization ('log1p', 'cp10k', or 'both')
                - 'log1p': log(1+x) transformation
                - 'cp10k': counts per 10k normalization
                - 'both': both cp10k and log1p (default)
            scale (bool): Whether to z-score scale genes after normalization

        Returns a polars DataFrame: columns = genename, *celltypes.
        """
        # adata.obs_names: barcodes in AnnData
        # self.barcode: polars DataFrame with columns ['gseid', 'srrid', 'barcode', 'celltype']
        # Only keep barcodes present in adata
        adata_barcodes = set(self.adata.obs_names)
        barcode_df = self.barcode.filter(
            pl.col("barcode").is_in(list(adata_barcodes))
        )
        # Map barcode to celltype
        barcode2celltype = dict(
            zip(
                barcode_df["barcode"].to_list(),
                barcode_df["celltype"].to_list(),
            )
        )
        # Build celltype list for adata.obs
        celltype_list = [
            barcode2celltype.get(bc, "unknown") for bc in self.adata.obs_names
        ]
        # Add celltype to adata.obs
        self.adata.obs["celltype"] = celltype_list

        # Apply normalization if requested
        if normalize:
            # Work with a copy to avoid modifying the original data
            adata_norm = self.adata.copy()

            # Perform normalization based on the selected method
            if norm_method in ["cp10k", "both"]:
                # Normalize to counts per 10,000
                sc.pp.normalize_total(adata_norm, target_sum=1e4)

            if norm_method in ["log1p", "both"]:
                # Log transform
                sc.pp.log1p(adata_norm)

            if scale:
                # Scale each gene to unit variance and zero mean
                sc.pp.scale(adata_norm, max_value=10)

            # Use the normalized data
            data_matrix = adata_norm.X
        else:
            # Use raw data
            data_matrix = self.adata.X

        # Prepare result dict
        result = {"genename": self.adata.var_names}
        for ct in celltypes:
            # Get mask for cells of this celltype
            mask = np.array(self.adata.obs["celltype"] == ct)
            if mask.sum() == 0:
                # No cells of this type, fill with nan
                mean_expr = np.full(self.adata.shape[1], np.nan)
            else:
                # data_matrix: cells x genes
                mean_expr = np.asarray(data_matrix[mask].mean(axis=0)).ravel()
            result[ct] = mean_expr

        # # Add metadata
        # result["gseid"] = self.gseid
        # result["srrid"] = self.srrid

        # Build polars DataFrame
        df = pl.DataFrame(result)
        return df

    def celltype_expr_high_confident_genes(
        self,
        normalize: bool = True,
        norm_method: str = "both",
        scale: bool = False,
        min_cells_expr: int = 10,
    ):
        """
        For this sample, compute normalized mean expression per gene for each cell type.

        Parameters:
            normalize (bool): Whether to normalize data before computing means.
            norm_method (str): Method for normalization ('log1p', 'cp10k', or 'both').
                - 'log1p': log(1+x) transformation
                - 'cp10k': counts per 10k normalization
                - 'both': both cp10k and log1p (default)
            scale (bool): Whether to z-score scale genes after normalization.
            min_cells_expr (int): Minimum number of cells a gene must be expressed in
                                (count > 0) within any cell type to be retained.
                                If 0, no filtering is applied.

        Returns:
            A polars DataFrame: columns = genename, *celltypes
        """
        # Only keep barcodes present in adata
        adata_barcodes = set(self.adata.obs_names)
        barcode_df = self.barcode.filter(
            pl.col("barcode").is_in(list(adata_barcodes))
        )

        # Map barcode to celltype
        barcode2celltype = dict(
            zip(
                barcode_df["barcode"].to_list(),
                barcode_df["celltype"].to_list(),
            )
        )

        # Assign celltype labels to adata.obs
        celltype_list = [
            barcode2celltype.get(bc, "unknown") for bc in self.adata.obs_names
        ]
        self.adata.obs["celltype"] = celltype_list
        celltypes = sorted(set(celltype_list))

        # Keep raw matrix for gene filtering
        raw_X = self.adata.X

        # Normalize if requested
        if normalize:
            adata_norm = self.adata.copy()

            if norm_method in ["cp10k", "both"]:
                sc.pp.normalize_total(adata_norm, target_sum=1e4)

            if norm_method in ["log1p", "both"]:
                sc.pp.log1p(adata_norm)

            if scale:
                sc.pp.scale(adata_norm, max_value=10)

            data_matrix = adata_norm.X
        else:
            data_matrix = raw_X

        # Filter genes by number of expressing cells per cell type
        if min_cells_expr > 0:
            gene_mask = np.zeros(self.adata.shape[1], dtype=bool)

            for ct in celltypes:
                ct_mask = np.array(self.adata.obs["celltype"] == ct)
                if ct_mask.sum() == 0:
                    continue
                # Use raw counts to determine expression presence
                expr_sub = raw_X[ct_mask] > 0
                expr_counts = np.asarray(expr_sub.sum(axis=0)).ravel()
                gene_mask |= expr_counts >= min_cells_expr

            data_matrix = data_matrix[:, gene_mask]
            gene_names = self.adata.var_names[gene_mask]
        else:
            gene_names = self.adata.var_names

        # Compute mean expression per cell type
        result = {"genename": gene_names}
        for ct in celltypes:
            mask = np.array(self.adata.obs["celltype"] == ct)
            if mask.sum() == 0:
                mean_expr = np.full(data_matrix.shape[1], np.nan)
            else:
                mean_expr = np.asarray(data_matrix[mask].mean(axis=0)).ravel()
            result[ct] = mean_expr

        return pl.DataFrame(result)


def scexpr(gseid, srrid, srrdir):
    """
    Process a single SCE object and return a Polars DataFrame with gene-level expression data by cell type.

    Parameters:
        gseid (str): Gene set ID
        srrid (str): Sequence read archive ID
        srrdir (Path): Directory containing the SCE data

    Returns:
        pl.DataFrame: DataFrame with columns gseid, srrid, genename, and expression values for each cell type
    """
    # gseid, srrid, srrdir = SRR[0, "gseid"], SRR[0, "srrid"], SRR[0, "srrdir"]
    # srrdir = Path(srrdir)
    sce = SCEXPR(gseid, srrid, srrdir)
    df = sce.celltype_expr()
    outfile = OUTDIR / f"{gseid}_{srrid}_celltype_gene_expr.csv"
    df.write_csv(outfile)
    print(f"Saved expression data to {outfile}")

    return df


def sceexpr_all(max_works: int = 20):
    all_rows = list(SRR.iter_rows(named=True))

    def process_row(row):
        gseid = row["gseid"]
        srrid = row["srrid"]
        srrdir = Path(row["srrdir"])
        print(f"Processing {gseid} {srrid} {srrdir}")
        try:
            return scexpr(gseid, srrid, srrdir)
        except Exception as e:
            print(f"Error processing {gseid} {srrid} {srrdir}: {e}")
            return None

    exprs = []
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=max_works
    ) as executor:
        for expr in tqdm(
            executor.map(process_row, all_rows), total=len(all_rows)
        ):
            exprs.append(expr)

    # valid_exprs = [
    #     e
    #     if e is not None
    #     else pl.DataFrame(
    #         {
    #             "gseid": ["Unknown"],
    #             "srrid": ["Unknown"],
    #             "genename": ["Unknown"],
    #             **{ct: [np.nan] for ct in celltypes},
    #         }
    #     )
    #     for e in exprs
    # ]
    # if valid_exprs:
    #     exprs_out = pl.concat(valid_exprs)
    #     exprs_out.write_csv(
    #         OUTDIR / "gse_srrid_gene_expr.csv",
    #     )
    # else:
    #     print("No valid expression data found.")
    #     return None


def collect_expr(max_workers: int = 20):
    csvs = OUTDIR.glob("*.csv")
    dfs = []

    def process_csv(csv):
        df = pl.read_csv(csv)
        gseid, srrid = csv.stem.split("_")[:2]
        df = df.with_columns(
            pl.lit(gseid).alias("gseid"),
            pl.lit(srrid).alias("srrid"),
        )
        # print(f"Processing {csv} with shape {df.shape} {gseid} {srrid}")
        return df

    with concurrent.futures.ThreadPoolExecutor(
        max_workers=max_workers
    ) as executor:
        for df in tqdm(executor.map(process_csv, csvs), total=len(list(csvs))):
            if df is not None:
                dfs.append(df)

    # Combine all dataframes
    if dfs:
        dfs = pl.concat(dfs)
    else:
        dfs = pl.DataFrame()
    dfs.write_csv(
        OUTDIR / "gse_srrid_celltype_gene_expr.csv",
    )


app = typer.Typer()


@app.command()
def EXPRALL(max_works: int = 20):
    """
    Process all SCE data and save the expression data to a CSV file.

    Parameters:
        max_works (int): Maximum number of worker threads to use
    """
    sceexpr_all(max_works)


@app.command()
def EXPRONE(
    gseid_srrid_srrdir_csv: str = typer.Argument(
        ..., help="GSEID,SRRID,SRRDIR"
    ),
):
    gseid, srrid, srrdir = gseid_srrid_srrdir_csv.strip().split(",")
    srrdir = Path(srrdir)
    scexpr(gseid, srrid, srrdir)


@app.command()
def COLLECT(max_worksers: int = 20):
    collect_expr(max_workers=max_worksers)


if __name__ == "__main__":
    app()
