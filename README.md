# Cotransitions
Statistical analysis of co-evolutionary transitions among genes.<br><br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7578797.svg)](https://doi.org/10.5281/zenodo.7578797) 
<br><br>
![cotr_logo](./Figures/cotr_logo.png
Pipeline in python3 and R to perform coevolutionary analysis
## Citation:
Dembech, Elena, Marco Malatesta, Carlo De Rito, Giulia Mori, Davide Cavazzini, Andrea Secchi, Francesco Morandin, and Riccardo Percudani. “Identification of Hidden Associations among Eukaryotic Genes through Statistical Analysis of Coevolutionary Transitions.” Proceedings of the National Academy of Sciences 120, no. 16 (April 18, 2023): e2218329120. https://doi.org/10.1073/pnas.2218329120.

## Installation

### Requirements
Any version of [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

```{bash}
git clone https://github.com/lab83bio/Cotransitions.git
cd Cotransitions
source cotr_install
```
The installation is successful if `(cotr)` is added to the command line prompt

```console
(cotr) user@pc:~/Cotransitions$ 
```

## Usage
To run the whole pipeline with **OrthoDB**:

```bash
./cotr_pipeline
```

[you can modify with your parameters the `cotr_pipeline` VARIABLES section] <br>
```bash
level="Eukaryota" #"Bacteria", "Archaea" (faster), "Mammalia", etc
tree="raxml" #raxml|ncbi|random
ladder=("RL" "LL" "NL") #tree orientation (RL=right-ladderized)
ncores=25
```
### Usage with OMA HOGs
To run the whole pipeline with OMA HOGs, please install cotr_extra envinroment with:

```{bash}
source cotr_extra_install
```
```console
(cotr_extra) user@pc:~/Cotransitions$ 
```


then you can run
```bash
./cotr_pipeline_OMA.sh
```

## Notebook usage

[*Concordant.ipynb*](https://github.com/lab83bio/Cotransitions/blob/master/Notebook/Concordant.ipynb) and
[*Validation.ipynb*](https://github.com/lab83bio/Cotransitions/blob/master/Notebook/Validation.ipynb) files can be opened with `jupyter-lab` included in `cotr_extra` conda environment <br>
The whole analysis can be performed by clicking ⏩ "Restart kernel and run all cells..."

## How transitions are counted and scored

A "transition" is a change in the gene distribution vector. When two genes (orthogroups) share the same transition at the same position (e.g. 1/0, 1/0), a concordant transition is counted.

![fig1 concordant](./Figures/Figure1B/Concordant_main_figures_revision1B.png)

Conversely, when two genes (orthogroups) share the opposite transition at the same position (e.g. 1/0, 0/1), a discordant transition is counted (see [this figure](./Figures/Figure1B/Concordant_Supplementary_Figures_revision_FigS1.png)).

The co-transition (cotr) score and significance are comuputed form the number of concordant (c), discordant (d), and total transitions (t1, t2).

The values shown for "[Case1](./Figures/Figure1B/Figure1B_Case1.csv.num)" can be obtained with the following command:
```bash
./cotr_transitions.py Figures/Figure1B/Figure1B_Case1.csv.num | ./cotr_Fisher.r -
```
The values shown for "[Case2](./Figures/Figure1B/Figure1B_Case2.csv.num)" can be obtained with the following command:
```bash
./cotr_transitions.py Figures/Figure1B/Figure1B_Case2.csv.num | ./cotr_Fisher.r -
```

Note that the cotr_score and scignificance discriminate among presence/absence patterns with the same Jaccard or Pearson scores, but different coevolutionary information.
