# MZINBVA: Variational approximation for multilevel zero-inflated negative-binomial models for association analysis in microbiome surveys

# Introduction
The microbiome data are over-dispersed with excess zeros. In addition, the observations from repeated measures in longitudinal studies are correlated. In order to address these challenges, we propose an effective variational approximation algorithm for fitting zero-inflated negative-binomial models for multilevel data and apply it to association analysis or differential abundance testing.

# Installation
You can install our MZINBVA package from Github
```r
install.packages("devtools")  
devtools::install_github("liudoubletian/MZINBVA")  
library(MZINBVA)  
```
# Basic Usage
## Data generation from three-level zero-inflated negative-binomial (ZINB) model
```r
sim_dat <- sim.data(sub.n, pos.n, vis.n, otu.n)
```
* `sub.n` : the number of subjects
* `pos.n` : the number of positions for each subject
* `vis.n` : the visit times for each subject
* `otu.n` : the number of OTUs
## Parameter estimation
```r
est_m <- step_alg(data, offset=NULL) 
```

## Association analysis or differential abundance testing for each taxon
```r
res.t <- data.pro(est_m,nodes=3) ### Wald test based on a variational sandwich covariance matrix
```
* `data` : a data frame that includs observed microbiome count data and covariates
* `est_m` : a list of estimated model parameters and variational parameters
* `nodes` : the required nodes for the parallel 

# Example
The following function shows that how to simulate a multilevel microbiome count data 
```r
sub.n <- 10
pos.n <- 3
vis.n <- 3
otu.n <- 100
sim_dat <- sim.data(sub.n, pos.n, vis.n, otu.n)
```
Next, we illustrate the association analysis or differential abundance testing procedure based on the example data
```r
sub.n <- 10
pos.n <- 3
vis.n <- 3
otu.n <- 100
sim_dat <- sim.data(sub.n, pos.n, vis.n, otu.n)
Y_mat <- sim_dat$Y_mat
data <- data.frame(Y=Y_mat[,1],ID=sim_dat$ID,cluster=sim_dat$cluster,x1=sim_dat$x1,x2=sim_dat$x2)
est_m <- step_alg(data,offset=NULL)
res.t <- data.pro(est_m,nodes=3)
```

[1] Tiantian Liu, Peirong Xu, Yueyao Du, Hui Lu, Hongyu Zhao, Tao Wang. (2021) MZINBVA: Variational approximation for multilevel zero-inflated negative-binomial models for association analysis in microbiome surveys.


