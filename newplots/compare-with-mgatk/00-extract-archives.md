# 00-extract-archives.sh

## Purpose

Unpack each sample archive into its own directory under the stage input root,
extracting only the seven files the later steps read. Every step from 01
onwards reads `${ISILON_BASE}/compare-with-mgatk/samples/<sample_id>/` and
never the archives, so this step is what turns a delivered `.zip` into a
sample this stage can run.

It is idempotent: a sample whose seven files are already present and non-empty
is skipped, so re-running the driver costs nothing after the first extraction.

## Inputs

The five archives, all directly under `${ISILON_BASE}/compare-with-mgatk`:

| sample_id | Archive | Chemistry |
| --- | --- | --- |
| GSE149689_GSM4509019_3PV3 | `GSE149689_GSM4509019_3PV3.zip` | SC3Pv3 |
| GSE163314_GSM4976997_3PV2 | `GSE163314_GSM4976997_3PV2.zip` | SC3Pv2 |
| GSE163668_GSM4995445_5PR2 | `GSE163668-GSM4995445_5PR2.zip` | SC5P-R2 |
| GSE181279_GSM5494116_5PPE | `GSE181279-GSM5494116_5PPE.zip` | SC5P-PE |
| GSE271107_GSM8369876_3PV3 | `GSE271107_GSM8369876_3PV3.zip` | SC3Pv3 |

Two archive names separate the GSE from the GSM with `-`; `sample_id`
normalises that to `_` so one token is safe as a directory name, a factor
level and a file name. The same mapping is the `SAMPLES` registry in
`config.R`.

## Outputs

One directory per sample, `${ISILON_BASE}/compare-with-mgatk/samples/<sample_id>/`,
each holding exactly these seven files:

- `cell.variant_stats.tsv.gz`
- `cell.variant_stats_mgatk_original.tsv.gz`
- `cell.cell_heteroplasmic_df.tsv.gz`
- `cell.cell_heteroplasmic_df_mgatk_original.tsv.gz`
- `cell.cell_heteroplasmic_df_raw.tsv.gz`
- `cell.depthTable.txt`
- `cell.coverage.txt.gz`

About 1.3 GB in total across the five samples. The allele-count matrices,
`cell.rds`, `cell.signac.rds` and the PNG summaries stay inside the archives;
nothing in this stage opens them.

## Run

```bash
cd "$(git rev-parse --show-toplevel)"
bash newplots/compare-with-mgatk/00-extract-archives.sh
```

No arguments: it always walks all five samples and skips the ones already
extracted. `run-all.sh` calls it first.

## Verify

```bash
ls -1 ~/project/scmocha/compare-with-mgatk/samples/
for d in ~/project/scmocha/compare-with-mgatk/samples/*/; do
  echo "$(ls -1 "$d" | wc -l) $(du -sh "$d" | cut -f1) $d"
done
```

Expected: five sample directories, seven files in each, about 1.3 GB in total.
A second run of the script prints `[skip]` for every sample and writes
nothing.

## Notes

The script fails loudly rather than producing a partial sample. A missing
archive stops it before any extraction, and a member that `unzip` did not
produce, or produced empty, stops it right after the sample it belongs to. A
sample directory left behind by a failed run is therefore always detected as
incomplete on the next run and re-extracted.

The `members` array here and `ARCHIVE_MEMBERS` in `config.R` are the same list
in two languages. Adding a file to one without the other means the extraction
and the reader disagree about what a complete sample is, and the disagreement
only surfaces as a missing-input error inside step 01.

`unzip -j` is deliberate: the members are taken out of whatever directory the
archive nests them in and written flat into the sample directory, which is
what `stage_inputs()` expects.
