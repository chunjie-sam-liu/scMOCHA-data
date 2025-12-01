#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-11-30 23:12:35
# @DESCRIPTION:
# @VERSION: v0.0.1


if [ -f $HOME/.bashrc ]; then
  source ${HOME}/.bashrc
fi

source ~/tools/anaconda3/bin/activate


RSCRIPT="/scr1/users/liuc9/tools/anaconda3/envs/renv/bin/Rscript"
SCRIPT="/home/liuc9/github/scMOCHA-data/src/06.1-collect-variants-new.R"
BASEDIR="/mnt/isilon/u01_project/large-scale/liuc9/raw"

GSE_LIST="
GSE226602
GSE155673
GSE226598
GSE157344
GSE279945
GSE206283
GSE153421
GSE155223
GSE163314
GSE214865
GSE163633
GSE163668
GSE167825
GSE175499
GSE159117
GSE149689
GSE184703
GSE168453
GSE261140
GSE148215
GSE175524
GSE220189
GSE149313
GSE174125
GSE233844
GSE162117
GSE147794
GSE166992
GSE143353
GSE161354
GSE188632
GSE235050
GSE154386
GSE181279
GSE164690
GSE171555
"

echo "===== Start Running All GSE Jobs ====="
echo

for g in $GSE_LIST; do
    echo "▶️  Running $g ..."
    $RSCRIPT $SCRIPT -g "$g" -b "$BASEDIR"
    status=$?

    if [ $status -ne 0 ]; then
        echo "❌ $g failed with exit code $status — continue next"
    else
        echo "✔️ $g completed successfully"
    fi

    echo "--------------------------------------"
    echo
done

echo "===== All Tasks Finished ====="