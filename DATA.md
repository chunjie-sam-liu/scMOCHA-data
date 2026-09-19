# DATA.md

Input provenance for `scMOCHA-data`. One row per input that a track reads.
Pipeline outputs do not belong here; they belong to the stage that produces
them.

Every row carries: exact repository-relative or env-relative path, origin
(a person or the producing script, never "the cluster"), date received or
built, the env-file variable that exposes it, and the tracks or stages that
consume it. A coordinate-bearing input also carries its genome build.

This file starts with the inputs of `newplots/compare-with-mgatk`. It is not
yet a complete inventory of the repository; rows are added as each track's
inputs are verified. Do not treat an absent row as evidence that an input has
no provenance.

---

## Env variables that point at data

| Variable | Value | Covered below |
| --- | --- | --- |
| `ISILON_BASE` | `/home/cliu68/project/scmocha` | yes, for the `compare-with-mgatk` subtree only |
| `DATADIR`, `BASEDIR` | `${REPODIR}/data` (symlink) | not yet inventoried |
| `CLEANDATADIR` | `${REPODIR}/analysis/zzz/clean-data` | not yet inventoried |
| `DUCKDB_PATH` | `${ZZZDIR}/clean-data/all_hetero_af.cell.duckdb.1.2.1` | not yet inventoried |
| `DUCKDB_PATH_COV` | `${ZZZDIR}/db/DUCKDB/cov.duckdb` | not yet inventoried |

`data/`, `results/`, `logs/` and `tmp/` are symlinks to storage outside the
repository. Use `find -L`, `du -L`, `rsync -L` on them.

---

## newplots/compare-with-mgatk

Paired scMOCHA and original-mgatk variant-calling output, one archive per
sample, produced by Ting. Both arms were run on the same allele counts, which
is what makes the criterion-level comparison valid; that is confirmed for
GSE181279 and **not yet confirmed in writing for the other seven** (see
Unknowns).

Root: `${ISILON_BASE}/compare-with-mgatk`.

| Path | Origin | Date | Variable | Consumed by |
| --- | --- | --- | --- | --- |
| `${ISILON_BASE}/compare-with-mgatk/GSE181279-GSM5494116_5PPE.zip` | Ting | 2026-09-09 | `ISILON_BASE` | `newplots/compare-with-mgatk` steps 00 to 07 |
| `${ISILON_BASE}/compare-with-mgatk/GSE271107_GSM8369876_3PV3.zip` | Ting | 2026-09-11 | `ISILON_BASE` | same |
| `${ISILON_BASE}/compare-with-mgatk/GSE163314_GSM4976997_3PV2.zip` | Ting | 2026-09-11 | `ISILON_BASE` | same |
| `${ISILON_BASE}/compare-with-mgatk/GSE163668-GSM4995445_5PR2.zip` | Ting | 2026-09-11 | `ISILON_BASE` | same |
| `${ISILON_BASE}/compare-with-mgatk/GSE149689_GSM4509019_3PV3.zip` | Ting | 2026-09-11 | `ISILON_BASE` | same |
| `${ISILON_BASE}/compare-with-mgatk/GSE155673_GSM4712895_3PV3.zip` | Ting | 2026-09-18 | `ISILON_BASE` | same |
| `${ISILON_BASE}/compare-with-mgatk/GSE188632_GSM5687372_3PV3.zip` | Ting | 2026-09-18 | `ISILON_BASE` | same |
| `${ISILON_BASE}/compare-with-mgatk/GSE220189_GSM6793474_3PV3.zip` | Ting | 2026-09-18 | `ISILON_BASE` | same |

Origin stated by Chun-Jie Liu on 2026-09-18. Dates are the file modification
times read from the filesystem: the first five on 2026-09-17, the last three on
2026-09-18. All eight carry the identical 17-file layout and the identical
14-column `variant_stats` schema, verified per archive before use.

Each archive holds the same 17-file flat layout. `00-extract-archives.sh`
unpacks seven of them per sample into
`${ISILON_BASE}/compare-with-mgatk/samples/<sample_id>/`:

| File | Content | Genome build |
| --- | --- | --- |
| `cell.variant_stats.tsv.gz` | scMOCHA per-variant statistics, 14 columns | rCRS / chrM, positions 1-16569 |
| `cell.variant_stats_mgatk_original.tsv.gz` | original mgatk per-variant statistics, identical schema | rCRS / chrM |
| `cell.cell_heteroplasmic_df.tsv.gz` | scMOCHA cell-by-variant AF, post-C3 | rCRS / chrM |
| `cell.cell_heteroplasmic_df_mgatk_original.tsv.gz` | mgatk cell-by-variant AF, post-C3; its barcode set is mgatk's kept-cell set | rCRS / chrM |
| `cell.cell_heteroplasmic_df_raw.tsv.gz` | scMOCHA cell-by-variant AF, all candidates | rCRS / chrM |
| `cell.depthTable.txt` | per-cell mean MT coverage, rounded to two decimals | not applicable |
| `cell.coverage.txt.gz` | per-cell per-position coverage, `pos,barcode,cov` | rCRS / chrM |

Not extracted and not read by this stage: `cell.{A,T,C,G}.txt.gz`, `cell.rds`,
`cell.signac.rds`, `MT_refAllele.txt`, the two PNGs,
`cromwell_glob_control_file`.

Coordinates are mitochondrial positions on the rCRS reference shipped with the
caller (`config/rCRS.MT.fasta`), 1-16569. Nothing in this stage joins them to a
nuclear-genome coordinate, so no liftover is involved.

The flat `cell.*` files still present at the root of
`${ISILON_BASE}/compare-with-mgatk` are an earlier manual extraction of
`GSE181279-GSM5494116_5PPE.zip`. Nothing reads them since 2026-09-17; the stage
reads `samples/GSE181279_GSM5494116_5PPE/` instead.

---

## Unknowns

- **Pipeline versions behind the eight `compare-with-mgatk` archives.** Origin
  is recorded (Ting, 2026-09-18). Still needed: the scMOCHA and mgatk versions
  used, and written confirmation that the scMOCHA and original-mgatk arms ran
  on identical allele-count matrices for every sample. Decision `D1` of
  `newplots/compare-with-mgatk/DECISION.md` assumes they did, and that
  assumption is the basis of the whole comparison. Schema and candidate counts
  are consistent with it for all eight samples, but consistency is not
  confirmation.

## Superseded

None yet.
