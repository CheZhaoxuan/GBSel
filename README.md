# GBSel

---

## Overview
`GBSel` is a new ensemble learning framework which combined gradient boosting strategy with stacking ensemble learning for genomic prediction, and this framework use ridge regression (RR), kernel ridge regression (KRR), and support vector regression (SVR) as base models.

## Installation
- The latest version:<br>
```
install.packages("devtools")
devtools::install_github("CheZhaoxuan/GBSel")
```
- If the installation fails, we recommend manually installing the dependency package:<br>
```
install.packages("reticulate");install.packages("parallel");install.packages("doParallel");
install.packages("foreach"); install.packages("elasticnet"); install.packages("Matrix"); 
devtools::install_github("xiaolei-lab/SIMER");
devtools::install_github("xiaolei-lab/rMVP");
```

## Data preparation
### Phenotype data
```
data(example_dataset)
head(phenotype,5)
|       |EarHT  |dpoll  |EarDia |
|:----: |:----: |:----: |:----: |
|33-16	|64.75  |64.5   | NaN   |
|38-11	|92.25  |68.5   | 37.897|
|4226  	|65.50  |59.5   |  32.21933 |
|4722 	|81.13  |71.5   |  32.42100 |
|A188   |27.50  |62.0   |   31.41900|
```

### Genotype data

### Genomic relationship matrix



## Usage
### load dependencies 


### start GBSel


### show result





