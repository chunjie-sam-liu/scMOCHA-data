#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-22 10:37:11
# @DESCRIPTION:
# @VERSION: v0.0.1


import concurrent.futures
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
CELLTYPES = ["B", "CD4_T", "CD8_T", "DC", "Mono", "NK", "other", "other_T"]


def load_h5(srrdir: Path) -> ad.AnnData:
    """
    Load h5 file and return AnnData object.
    """
    h5_filename = srrdir / "filtered_feature_bc_matrix.h5"
    print(h5_filename)
    try:
        adata = sc.read_10x_h5(h5_filename)
        adata.var_names_make_unique()
        print(f"Loaded {h5_filename}")
        return adata
    except Exception as e:
        print(f"Error loading {h5_filename}: {e}")
        return None


adatas = {}
for (
    gseid,
    srrid,
    srrdir,
) in SRR.head(100).iter_rows():
    # print(gseid, srrid, srrdir)
    srrdir = Path(srrdir)
    sampleid = f"{gseid}_{srrid}"
    adatas[sampleid] = load_h5(srrdir)

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)
