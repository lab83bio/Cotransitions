#!/usr/bin/env bash
conda create  --name cotr --file cotr.txt
conda activate cotr
Rscript -e 'install.packages("corpora",repos="http://cran.us.r-project.org")'
