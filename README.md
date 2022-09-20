# Cotransitions
Statistical analysis of co-evolutionary transitions among genes.<br><br>
[![DOI](https://zenodo.org/badge/DOI/xxxx/zenodo.xxxx.svg)](https://doi.org/xxxx/zenodo.xxx) 
<br><br>


Pipeline in python3 and R to perform coevolutionary analysis
## Installation

### Requirements
Any version of [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

```{bash}
git clone https://github.com/lab83bio/Cotransitions.git
cd Cotransitions
source install.sh
```
The installation is successful if `(cotr)` in the command line prompt

```console
(cotr) user@pc:~/Cotransitions$ 
```

## Usage
To run the whole pipeline, please modify with your parameters in the `cotr_pipeline` VARAIBLES section. <br>
For example:
```bash
level="Eukaryota" #"Bacteria", "Archaea" (faster), "Mammalia", etc
tree="raxml" #raxml|ncbi|random
ladder=("RL" "LL" "NL") #tree orientation (RL=right-ladderized)
ncores=25
```
then you can run
```bash
./cotr_pipeline
```
## Notebook usage
[*Concordant.ipynb*](https://github.com/lab83bio/Cotransitions/blob/master/Notebook/Concordant.ipynb) and
[*Validation.ipynb*](https://github.com/lab83bio/Cotransitions/blob/master/Notebook/Validation.ipynb) files can be opened with `jupyter-lab` included in `cotr` conda environment <br>
The whole analysis can be performed by clicking ⏩ "Restart kernel and run all cells..."



## License

Cotr scripts are available under the MIT licence


