#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-11-26 17:09:06
# @DESCRIPTION:

# Number of input parameters
param=$#

gseids=(GSE149689 GSE155223 GSE155673 GSE157344 GSE163668 GSE166992 GSE171555 GSE181279 GSE226602)
basedir="/home/liuc9/github/scMOCHA-data/data"

update_gsmid() {
  local gseid=$1
  local gsmid=$2

  data_dir="$basedir/$gseid"
  targz_dir="$data_dir/targz"

  gsm_dir="$targz_dir/$gsmid"

  cd "$gsm_dir" || exit
  echo "Updating $gseid and $gsmid"

  source /home/liuc9/tools/anaconda3/etc/profile.d/conda.sh
  conda activate scmocha
  # cell cluster annotation
  Rscript /home/liuc9/github/scMOCHA/bin/azimuth.R \
    -h5file filtered_feature_bc_matrix.h5 \
    -npcs 10 \
    -reso 0.1 \
    -refname_celllevel refname=pbmcref celllevel=celltype.l1 \
    -nFeature_RNA_min 500 \
    -nFeature_RNA_max 8000 \
    -percent_mt_max 75.0 \
    -percent_ribo_max 50.0 \
    -percent_Lagest_Gene_max 50.0
  # cell level variant calling
  python /home/liuc9/github/scMOCHA/bin/variant_calling_cell_raw.py ./ cell 16569 10 MT
  # cluster level variant calling
  python /home/liuc9/github/scMOCHA/bin/variant_calling_cluster.py ./ cluster 16569 10 MT
  # plot scMOCHA results
  Rscript /home/liuc9/github/scMOCHA/bin/scMOCHA.R \
    -m cell_meta_data.tsv \
    -b barcode_cluster.tsv \
    -ceh cell.cell_heteroplasmic_df.tsv.gz \
    -cec cell.coverage.txt.gz \
    -clh cluster.cell_heteroplasmic_df.tsv.gz \
    -clc cluster.coverage.txt.gz \
    -chr cell.cell_heteroplasmic_df_raw.tsv.gz \
    -p /home/liuc9/github/scMOCHA/bin/get_variants_info.pl \
    -j /scr1/users/liuc9/tools/haplogrep3 \
    -s /mnt/isilon/xing_lab/liuc9/refdata/mitomaster/mitomap_sqlite_20230525.sqlite3 \
    -conda_root /home/liuc9/tools/anaconda3 \
    -conda_env scmocha

}

update_gse() {
  local gseid=$1
  echo "Updating $gseid"

  data_dir="$basedir/$gseid"
  targz_dir="$data_dir/targz"
  # find gsmid in the targz_dir

  gsmids=$(find "$targz_dir" -maxdepth 1 -type d -name 'GSM*' -exec basename {} \;)

  for gsmid in $gsmids; do
    echo "Updating $gseid and $gsmid"
    update_gsmid $gseid $gsmid &
  done

}

# update_gse GSE149689
update_all_gse() {
  for gseid in "${gseids[@]}"; do
    update_gse "$gseid"
  done
}

update_all_gse
