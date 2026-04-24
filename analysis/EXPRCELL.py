#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-16 14:40:50
# @DESCRIPTION:
# @VERSION: v0.0.1


import concurrent.futures
from pathlib import Path

import numpy as np
import polars as pl
import scanpy as sc
import typer
from rich import print
from tqdm import tqdm

OUTDIR = Path("/home/liuc9/github/scMOCHA-data/analysis/zzz/db/EXPR")
SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)

BARCODE_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/clean-data/barcode_celltype.feather"
)

# Fix: Use BARCODE_FILENAME instead of SRR_FILENAME with read_ipc
BARCODE = pl.read_ipc(BARCODE_FILENAME)
celltypes = ["B", "CD4_T", "CD8_T", "DC", "Mono", "NK", "other", "other_T"]
CELLTYPES = celltypes


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
        self.adata = self.normalize()
        self.annotate_celltype()

    def __repr__(self):
        return f"SCE(adata={self.adata})"

    def normalize(self):
        adata = self.qc()
        # Saving count data
        adata.layers["counts"] = adata.X.copy()
        # Normalizing to median total counts
        sc.pp.normalize_total(adata)
        # Logarithmize the data
        sc.pp.log1p(adata)
        return adata

    def qc(self):
        adata = self.load_h5()
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=False
        )

        # filter cells and genes
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=10)

        return adata

    def load_h5(self):
        adata = sc.read_10x_h5(self.h5_file)
        adata.var_names_make_unique()

        # mitochondrial genes, "MT-" for human, "Mt-" for mouse
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        # ribosomal genes
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        # hemoglobin genes
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

        return adata

    def annotate_celltype(self):
        adata_barcodes = set(self.adata.obs_names)
        barcode_df = self.barcode.filter(
            pl.col("barcode").is_in(list(adata_barcodes))
        )
        barcode2celltype = dict(
            zip(
                barcode_df["barcode"].to_list(),
                barcode_df["celltype"].to_list(),
            )
        )
        celltype_list = [
            barcode2celltype.get(bc, "unknown") for bc in self.adata.obs_names
        ]
        self.adata.obs["celltype"] = celltype_list

    def celltype_mean_expr(self):
        result = {"genename": self.adata.var_names}
        data_matrix = self.adata.X
        for ct in CELLTYPES:
            # Get mask for cells of this celltype
            mask = np.array(self.adata.obs["celltype"] == ct)
            if mask.sum() == 0:
                # No cells of this type, fill with nan
                mean_expr = np.full(self.adata.shape[1], np.nan)
            else:
                # data_matrix: cells x genes
                mean_expr = np.asarray(data_matrix[mask].mean(axis=0)).ravel()
            result[ct] = mean_expr

        df = pl.DataFrame(result)
        # print(df.head())
        return df

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
    df = sce.celltype_mean_expr()
    outfile = OUTDIR / f"{gseid}_{srrid}_celltype_gene_expr.csv"
    df.write_csv(outfile)
    print(f"Saved expression data to {outfile}")

    return df


def sceexpr_all(max_workers: int = 20):
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
        max_workers=max_workers
    ) as executor:
        for expr in tqdm(
            executor.map(process_row, all_rows), total=len(all_rows)
        ):
            exprs.append(expr)


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
    dfs.write_ipc(
        OUTDIR / "gse_srrid_celltype_gene_expr.feather",
    )
    print(
        f"Saved combined expression data to {OUTDIR / 'gse_srrid_celltype_gene_expr.feather|csv'}"
    )


app = typer.Typer()


@app.command()
def EXPRALL(max_workers: int = 20):
    """
    Process all SCE data and save the expression data to a CSV file.

    Parameters:
        max_workers (int): Maximum number of worker threads to use
    """
    sceexpr_all(max_workers)


@app.command()
def EXPRONE(
    gseid_srrid_srrdir_csv: str = typer.Argument(
        ..., help="GSEID,SRRID,SRRDIR"
    ),
):
    """
    Generate expression data for a single-cell RNA-seq dataset.

    This function processes a single-cell RNA-seq dataset identified by a GEO Series ID,
    a Sample ID, and the directory containing the raw data files. It calls the `scexpr`
    function to process the data.
    Parameters
    ----------
    gseid_srrid_srrdir_csv : str
        A comma-separated string containing:
        - GSEID: GEO Series ID (e.g., "GSE155673")
        - SRRID: GEO Sample ID (e.g., "GSM4712885")
        - SRRDIR: Path to the directory containing the raw data files
        Example: "GSE155673,GSM4712885,/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE155673/final/GSM4712885"
    Returns
    -------
    None
        The function calls `scexpr` to process the data but doesn't return anything directly.
    """
    # gseid_srrid_srrdir_csv = "GSE155673,GSM4712885,/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE155673/final/GSM4712885"
    gseid, srrid, srrdir = gseid_srrid_srrdir_csv.strip().split(",")
    srrdir = Path(srrdir)
    scexpr(gseid, srrid, srrdir)


@app.command()
def COLLECT(max_workers: int = 20):
    """
    Collects expression data using the collect_expr function with parallelization.

    This function serves as a wrapper around the collect_expr function, providing
    a simplified interface for collecting expression data with configurable
    parallelization.

    Parameters
    ----------
    max_workers : int, default=20
        The maximum number of worker processes to use for parallel execution.
        Controls the level of parallelization when collecting expression data.

    Returns
    -------
    None
        This function does not return any value. The results of the data
        collection are typically saved to disk by the underlying collect_expr
        function.

    See Also
    --------
    collect_expr : The underlying function that performs the actual data collection.

    Examples
    --------
    >>> COLLECT()  # Use default 20 workers
    >>> COLLECT(max_workers=8)  # Use 8 workers for less resource usage
    """
    collect_expr(max_workers=max_workers)


@app.command()
def workflow(
    max_workers: int = 20,
):
    """
    Run the entire workflow: process all SCE data, collect expression data, and save to CSV.

    Parameters:
        max_workers (int): Maximum number of worker threads to use
    """
    EXPRALL(max_workers)
    COLLECT(max_workers)


if __name__ == "__main__":
    app()
