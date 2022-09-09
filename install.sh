#!/usr/bin/env bash
conda create  --name cotr --file cotr.txt
conda activate cotrt
Rscript -e 'install.packages("corpora")'
