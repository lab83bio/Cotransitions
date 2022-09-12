#!/usr/bin/env bash

conda create  --name cotr --file cotr.txt
conda activate cotr
Rscript -e 'install.packages("corpora",repos="http://cran.us.r-project.org")'
pip3 install git+https://github.com/mojaie/pygosemsim.git
