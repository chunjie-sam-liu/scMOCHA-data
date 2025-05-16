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

    def celltype_expr(self, normalize=True, norm_method="log1p", scale=False):
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


if __name__ == "__main__":
    app()
