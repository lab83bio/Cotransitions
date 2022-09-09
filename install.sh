#!/usr/bin/env bash
conda create  --name ct --file ct.txt
conda activate ct
Rscript -e 'install.packages("corpora")'
