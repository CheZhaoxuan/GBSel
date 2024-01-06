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
```
|       |EarHT  |dpoll  |EarDia |
|:----: |:----: |:----: |:----: |
|33-16	|64.75  |64.5   |NaN    |
|38-11	|92.25  |68.5   |37.897 |
|4226  	|65.50  |59.5   |32.219 |
|4722 	|81.13  |71.5   |32.421 |
|A188   |27.50  |62.0   |31.419 |
- Note: The format is a matrix, and the example phenotype shows row and column names.

### Genotype data
```
data(example_dataset)
genotype[1:5,1:5]
```
|       |PZB00859.1 |PZA01271.1 |PZA03613.2 |PZA03613.1 | PZA03614.2 |
|:----: |:----:     |:----:     |:----:     |:----:     |:----:     |
|33-16	|2|0|0|2|2|
|38-11	|2|2|0|2|2|
|4226  	|2|0|0|2|2|
|4722 	|2|2|0|2|2|
|A188   |0|0|0|2|2|
- Note: The format is a matrix, and the example genotype shows row and column names.

## Usage
### load dependencies 
- Before running, make sure that the available **numpy** library and **scikit-learn** library in **python** are available
```


```
path <- reticulate::py_discover_config()
py.lib <- python.env(env="python", env.path=path$python, module.check=TRUE)	
```

### start GBSel


### show result





