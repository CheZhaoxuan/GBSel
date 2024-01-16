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
path <- reticulate::py_discover_config()
py.lib <- python.env(env="python", env.path=path$python, module.check=TRUE)	
```
- **Notice**: Note the environment and path of python, and `path` can be flexibly specified according to the  python you use. Python 3.6 or above is recommended.

### data input
- We used genomic relationship as the learning data, and the training set and test set were divided as examples
```
X <- genotype
X.FR <- Feature.reduce(X = X, method="GRM")
y.train <- phenotype[c(1:200),1]
X.train <- X.FR[c(1:200),]
X.test <- X.FR[-c(1:200),]
```
### start GBSel
```
Result <- GBSel(y_train=y.train,X_train=X.train,X_test=X.test,py.lib=py.lib,Parallel=FALSE)
```
### show result
```
str(Result)
List of 4
 $ MX_train   : num [1:200, 1:22] 56 58.8 55.2 57.8 53.1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:200] "33-16" "38-11" "4226" "4722" ...
  .. ..$ : chr [1:22] "Ridge_nIter0" "Ridge_nIter1" "Ridge_nIter2" "Ridge_nIter3" ...
 $ MX_test    : num [1:64, 1:22] 62.3 66.2 61.9 60.8 69.4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:22] "Ridge_nIter0" "Ridge_nIter1" "Ridge_nIter2" "Ridge_nIter3" ...
 $ GBSel_PV_tr: num [1:200(1d)] 58.1 57.1 56.6 56.2 55.1 ...
 $ GBSel_PV_te: num [1:64(1d)] 62.8 65.5 64.3 60.7 68.8 ...
```
- **MX_train**: Meta features for training set
- **MX_test**: Meta features for test set
- **GBSel_PV_tr**: Fit phenotype values of training set
- **GBSel_PV_te**: Prediction phenotype values of test set





