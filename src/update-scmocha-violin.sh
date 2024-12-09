#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-11-26 17:09:06
# @DESCRIPTION:

# Number of input parameters
param=$#

gseids=(GSE149689 GSE155223 GSE155673 GSE157344 GSE163668 GSE166992 GSE171555 GSE181279 GSE226602)
# gseids=(GSE149689 GSE155223 GSE155673 GSE157344 GSE163668 GSE166992 GSE171555 GSE181279)
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

  if [[ -f chemistry.csv ]]; then
    chem=$(sed -n '2p' chemistry.csv | cut -f1 -d,)
  else
    case "$gseid" in
    GSE163668 | GSE155223 | GSE171555)
      chem="SC5P-R2"
      ;;
    GSE149689 | GSE155673 | GSE157344)
      chem="SC3Pv3"
      ;;
    GSE166992 | GSE226602 | GSE181279)
      chem="SC5P-PE"
      ;;
    *)
      echo "Unknown GSE ID: $gseid"
      exit 1
      ;;
    esac
  fi

  # cell cluster annotation
  echo "$gseid $gsmid azimuth.R"
  Rscript /home/liuc9/github/scMOCHA/bin/azimuth.R \
    -h5file filtered_feature_bc_matrix.h5 \
    -refname_celllevel refname=pbmcref celllevel=celltype.l1 \
    -c ${chem}

  # cell level variant calling
  echo "$gseid $gsmid variant_calling_cell_raw.py"
  python /home/liuc9/github/scMOCHA/bin/variant_calling_cell_raw.py ./ cell 16569 10 MT
  # cluster level variant calling
  echo "$gseid $gsmid variant_calling_cluster.py"
  python /home/liuc9/github/scMOCHA/bin/variant_calling_cluster.py ./ cluster 16569 10 MT

  # plot scMOCHA results
  echo "$gseid $gsmid plot_scMOCHA.R"
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

update_all_gse() {
  for gseid in "${gseids[@]}"; do
    update_gse "$gseid"
  done
}

# update_gse GSE155673

# update_gsmid GSE226602 GSM7080058

update_all_gse
