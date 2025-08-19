#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-08-18 17:23:08
# @DESCRIPTION:
# @VERSION: v0.0.1


# %%
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc

# %% core gene stress


core_genes_stress = pd.read_csv(
    "/home/liuc9/github/MTRNA-sc-cancer/auxiliary_data/coregene_df-FALSE-v3.csv"
)

list_core_genes = core_genes_stress[
    core_genes_stress["logFC"] > 0
].gene_symbol.to_numpy()

red_core_genes = core_genes_stress.head(40).gene_symbol.to_numpy()

# %% dissociation genes

dissociation_genes = (
    pl.read_csv(
        "/home/liuc9/github/scMOCHA-data/config/dissociation_genes-vanDenBrink2017.csv",
        has_header=False,
    )
    .select(pl.col("column_1"))
    .to_series()
    .str.to_uppercase()
    .to_numpy()
)

# %% dissociation genes machado
dissociation_genes_machado = (
    pl.read_csv(
        "/home/liuc9/github/scMOCHA-data/config/dissociation_Machado2021.csv",
        has_header=False,
    )
    .select(pl.col("column_1"))
    .to_series()
    .str.to_uppercase()
    .to_numpy()
)

# %% common dissociation genes
common_disso_genes = np.intersect1d(
    list_core_genes.astype(str),
    np.intersect1d(
        dissociation_genes.astype(str), dissociation_genes_machado.astype(str)
    ),
)
len(common_disso_genes)

# %% GOCC
gocc = {}
with open(
    "/home/liuc9/github/scMOCHA-data/config/GO_Cellular_Component_2013.txt",
    "r",
) as f:
    lines = f.readlines()
    for line in lines:
        vals = line.split("\t")
        # print(vals)
        gocc[vals[0]] = vals[2:-1]

pathways_gocc = ["mitochondrion (GO:0005739)", "cytoplasm (GO:0005737)"]

go_sigs = {}
for path in pathways_gocc:
    go_sigs[path] = gocc[path]

go_sigs["Dissociation stress"] = common_disso_genes
go_sigs.keys()


# %% load h5ad dataset
adata = sc.read_h5ad(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/EXPRMERGE/merged_normalized_data.h5ad"
)
adata.obs["Transcriptome variance"] = adata.to_df().var(axis=1)

# %%
for sig in go_sigs:
    sc.tl.score_genes(
        adata, gene_list=go_sigs[sig], score_name=sig.capitalize()
    )
# %%
adata

# %% save to h5ad
adata.write(
    "/home/liuc9/github/scMOCHA-data/analysis/zzz/db/EXPRMERGE/merged_normalized_data.disso_stress_score.h5ad"
)

# %%
