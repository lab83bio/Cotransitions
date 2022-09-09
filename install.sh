#!/usr/bin/env bash

conda create --name ct --file ct.txt -y
conda activate ct
Rscript -e 'install.packages("corpora")'
