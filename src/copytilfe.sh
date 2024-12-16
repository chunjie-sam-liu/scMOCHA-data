#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-08-23 18:15:05
# @DESCRIPTION:

# Number of input parameters
param=$#

targzfile=$1
outdir=$2

filenames=$(cat ${targzfile})

for filename in ${filenames}; do
  # echo "Copying ${filename} to ${outdir}"
  cmd="cp ${filename} ${outdir} &"
  echo $cmd
  eval $cmd
done
