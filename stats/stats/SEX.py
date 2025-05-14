#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-05-14 12:14:10
# @DESCRIPTION:
# @VERSION: v0.0.1


from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import scanpy as sc
import typer

YGENES = ["RPS4Y1", "KDM5D", "DDX3Y"]
XGENES = ["XIST"]

SRR_FILENAME = Path(
    "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir.csv"
)
SRR = pl.read_csv(SRR_FILENAME)
SRR


class SCESEX:
    """
    Single-cell expression (SCE) class for analyzing single-cell RNA-seq data.
    """

    def __init__(self, h5_file: Path):
        """
        Initialize the SCE class with an AnnData object.

        Parameters
        ----------
        adata : sc.AnnData
            AnnData object containing single-cell RNA-seq data.
        """
        # h5_file = Path(
        #     "/home/liuc9/github/scMOCHA-data/data/GSE214865/final/GSM6616991/filtered_feature_bc_matrix.h5"
        # )
        self.h5_file = Path(h5_file)
        self.adata = sc.read_10x_h5(h5_file)
        self.adata.var_names_make_unique()
        self.mean_expr = None
        self.gene_expr = None

    def __repr__(self):
        return f"SCE(adata={self.adata})"

    def check_xy_genes(self):
        """
        Check if the X and Y genes are present in the AnnData object.
        """
        xgenes = self.adata.var_names.intersection(XGENES)
        ygenes = self.adata.var_names.intersection(YGENES)
        if len(xgenes) == 0:
            raise ValueError(f"X genes {XGENES} not found in AnnData object.")
        if len(ygenes) == 0:
            raise ValueError(f"Y genes {YGENES} not found in AnnData object.")
        return xgenes, ygenes

    def pseudo_bulk_ex_mean(self):
        """
        Create pseudo-bulk data from the AnnData object.
        """
        self.adata.layers["log1p"] = self.adata.X.copy()
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)

        gene_expr = pd.DataFrame(
            self.adata.X.toarray(),
            columns=self.adata.var_names,
        )

        mean_expr = gene_expr.mean(axis=0)
        self.gene_expr = gene_expr
        self.mean_expr = mean_expr
        return mean_expr

    def plot_genes_violin(self, keys):
        """
        Plot violin plots for the specified genes.

        Parameters
        ----------
        keys : list
            List of gene names to plot.
        """
        # outdir = self.h5_file.parent
        # outplotfile = outdir / "sex_estimate_violin.pdf"
        srrid = self.h5_file.parent.name
        gseid = self.h5_file.parent.parent.parent.name
        srrdir = self.h5_file.parent
        sc.settings.figdir = srrdir
        sc.pl.violin(
            self.adata,
            keys=keys,
            groupby=None,
            jitter=0.4,
            rotation=45,
            use_raw=False,
            show=False,
            save=f"_{srrid}_sex_markers.pdf",
        )

    def sex_score(self):
        self.pseudo_bulk_ex_mean()
        xgenes, ygenes = self.check_xy_genes()
        self.plot_genes_violin(keys=XGENES + YGENES)
        y_score = self.mean_expr[ygenes].mean() if len(ygenes) > 0 else 0
        x_score = self.mean_expr[xgenes].mean() if len(xgenes) > 0 else 0
        estimated_sex = "Unknown"
        if y_score > 0.5 and x_score < 0.5:
            estimated_sex = "Male"
        elif y_score < 0.2 and x_score > 0.5:
            estimated_sex = "Female"

        return pl.DataFrame(
            {
                "sex": estimated_sex,
                "x_score": x_score,
                "y_score": y_score,
            }
        )


# scesex = SCESEX(
#     # h5_file="/mnt/isilon/u01_project/large-scale/ting/raw/GSE161354/final/GSM4905214/filtered_feature_bc_matrix.h5"
#     h5_file="/mnt/isilon/u01_project/large-scale/liuc9/raw/GSE155223/final/GSM4697614/filtered_feature_bc_matrix.h5"
# )

# a = scesex.sex_score()


def estimate_sex(gseid: str, srrid: str, srrdir: Path):
    # gesid, srrid, srrdir = SRR[0, "gseid"], SRR[0, "srrid"], SRR[0, "srrdir"]
    # srrdir = Path(srrdir)
    h5_file = srrdir / "filtered_feature_bc_matrix.h5"
    scesex = SCESEX(h5_file)
    sexest = scesex.sex_score()
    return sexest


def estimate_sex_all_srr():
    sexests = []
    for row in SRR.iter_rows(named=True):
        gseid = row["gseid"]
        srrid = row["srrid"]
        srrdir = Path(row["srrdir"])
        print(f"Processing {gseid} {srrid} {srrdir}")
        sexest = estimate_sex(gseid, srrid, srrdir)
        sexests.append(sexest)
    sexests = pl.concat(sexests)
    SRR.with_columns(
        sexests,
    ).write_csv(
        "/home/liuc9/github/scMOCHA-data/stats/stats/zzz/clean-data/gse_srrid_srrdir_sex.csv",
        include_header=True,
    )


app = typer.Typer()


@app.command()
def SEX():
    estimate_sex_all_srr()


if __name__ == "__main__":
    SEX()
