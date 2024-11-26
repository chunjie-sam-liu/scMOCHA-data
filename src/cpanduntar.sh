#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-11-26 15:53:38
# @DESCRIPTION:

# Number of input parameters
param=$#
basedir="/home/liuc9/github/scMOCHA-data/data"
gseids=(GSE149689 GSE154567 GSE155223 GSE155673 GSE157344 GSE163668 GSE166992 GSE171555 GSE181279 GSE226602)

cp_targz_dir() {
  local gseid=$1
  echo "Updating $gseid"
  # Update the target directory
  data_dir="$basedir/$gseid"
  targz_dir="$data_dir/targz"

  # check if $targz_dir has any tar.gz file, if not run ${data_dir}/05.${gseid}.scmocha.cptargz.sh
  if [ -z "$(ls -A $targz_dir)" ]; then
    echo "No tar.gz files in $targz_dir"
    cmd="bash ${data_dir}/05.${gseid}.scmocha.cptargz.sh"
    echo $cmd
    eval $cmd
  else
    echo "There are tar.gz files in $targz_dir"
  fi
}

unzip_targz_dir() {
  local gseid=$1
  echo "Updating $gseid"
  # Update the target directory
  data_dir="$basedir/$gseid"
  targz_dir="$data_dir/targz"
  # unzip_dir="$data_dir/unzip"

  # check if $unzip_dir has any tar.gz file, if not run ${data_dir}/06.${gseid}.scmocha.unzip.sh
  if [ -z "$(ls -A $unzip_dir)" ]; then
    # echo "No tar.gz files in $unzip_dir"
    bash ${data_dir}/06.${gseid}.scmocha.clear.sh
    bash ${data_dir}/07.${gseid}.scmocha.untargz.sh

  else
    echo "There are tar.gz files in $unzip_dir"
  fi
}

# update_targz_dir "GSE171555"
cp_targz_dirs() {
  gseids_for_targz=(GSE157344 GSE171555)
  for gseid in ${gseids_for_targz[@]}; do
    update_targz_dir $gseid
  done
}

# update_targz_dirs

unzip_targz_dirs() {
  gseids_for_unzip=(GSE157344 GSE171555)
  for gseid in ${gseids_for_unzip[@]}; do
    unzip_targz_dir $gseid
  done
}

unzip_targz_dirs
