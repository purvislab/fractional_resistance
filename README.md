# fractional_resistance

## Overview
<p>
  <img src="/doc/overview.png"/>
</p>

## Introduction
We performed multiplex, single-cell imaging to build a proteomic profile for individual tumor cells treated with palbociclib and identified fractionally resistant cells both in a cell culture model of ER+ breast cancer (T47D), as well as from live primary tumor cells resected from a patient. We show that ER+/HER2- tumor cells show subtle, cell-to-cell differences in core cell cycle regulators that allow a subset of tumor cells to escape CDK4/6 inhibitor therapy. For more details, please read the associated preprint. [Zikry TM†, Wolff SC†, Ranek JS, Davis H, Naugle A, Whitman AA, Kosorok MR, Spanheimer PM*, and Purvis JE*. Cell cycle plasticity underlies fractional resistance to palbociclib in ER+/HER2- breast tumor cells. _bioRxiv_. 2023.]()

## Installation
* You can clone the git repository by, 
```
git clone https://github.com/purvislab/fractional_resistance.git
```
* Once you've clone the repository, please change your working directory into this folder.
```
cd fractional_resistance
```
## Dependencies

Given that there are a number of required python and R packages, we recommend that you create a conda environment as follows. First, create the conda environment using the provided yml file. This contains most of the installation instructions.
```
conda env create -f venv_fr.yml
```
Once the environment is created, you can activate it by,
```
conda activate venv_fr
```
Then you can finish installing the required R packages by,
```
RScript install_R_dependencies.R
```
## Data access
You can download all of the preprocessed single-cell datasets (`.h5ad` files) from the [Zenodo](https://doi.org/10.5281/zenodo.7930054) repository. Each `.h5ad` object contains a sample profiled with [iterative indirect immunofluorescence imaging](https://pubmed.ncbi.nlm.nih.gov/30072512/) across three treatment conditions (untreated, 10 nM, and 100 nM palbociclib). Here, image preprocessing was performed using the [4i_analysis](https://github.com/purvislab/4i_analysis) pipeline to quantify the median intensity of core cell cycle regulators for individual cells. Datasets were subsequently preprocessed by performing feature selection and data standardization. Cells were annotated according to cell cycle state (proliferation, arrest) and phase (G0, G1, S, G2/M).  

## Description
To ascertain cell-cell differences in core cell cycle regulators that allow tumor cells to escape CDK4/6 inhibitor therapy, we performed a series of computational analyses as described below.

* `sketch_integrate.ipynb` - Downsamples data (2000 cells per treatment condition per sample) using [Kernel Herding sketching](https://dl.acm.org/doi/abs/10.1145/3535508.3545539). Then integrates T47D and primary tumor sample into a shared latent space using [TRANSACT](https://www.pnas.org/doi/10.1073/pnas.2106682118).
* `ci.R` - Computes 95% confidence intervals for the difference in mean expression for proliferative untreated cells vs. 10 nM of palbociclib and proliferative untreated cells vs. 100 nM palbociclib.
* `logistic_regression.ipynb` - Performs logistic regression to identify changes in phase-specific features for T47D and primary tumor samples upon treatment with palbociclib. 
* `ti.R` - Performs trajectory inference with [Slingshot](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0) for each sample (e.g. T47D, primary tumor) and treatment condition (e.g. untreated, 10 nM, 100 nM palbociclib), followed by trajectory alignment with [TrAGEDy](https://www.biorxiv.org/content/10.1101/2022.12.21.521424v1) to compare continuous proteomic expression profiles of the cell cycle across treatment conditions.

## License
This software is licensed under the MIT license (https://opensource.org/licenses/MIT).
