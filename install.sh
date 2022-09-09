#!/usr/bin/env bash
conda create  --name cotr --file cotr.txt
conda activate cotr
<<<<<<< HEAD
Rscript -e 'install.packages("corpora",repos="http://cran.us.r-project.org")'
=======
Rscript -e 'install.packages("corpora")'
>>>>>>> 06c1b8e4f79196087a39a889e3d9bfe052b3b61e
