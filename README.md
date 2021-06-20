# MZINBVA: Variational approximation for multilevel zero-inflated negative-binomial models for association analysis in microbiome surveys

# Description
The microbiome count data are over-dispersed with excess zeros. In addition, the observations from repeated measures are correlated. In order to address these challenges, we propose an effective variational approximation method, MZINBVA, for fitting zero-inflated negative-binomial models for multilevel data and apply it to association analysis or differential abundance testing.

# Installation
You can install our MZINBVA package from Github
```r
install.packages("devtools")  
devtools::install_github("liudoubletian/MZINBVA")  
library(MZINBVA)  
```
# Basic Usage
## Data generation from three-level zero-inflated negative-binomial model
```r
sim_dat <- sim.data(sub.n, pos.n, vis.n, otu.n)
```
* `sub.n` : the number of subjects
* `pos.n` : the number of positions for each subject
* `vis.n` : the visit times for each subject
* `otu.n` : the number of OTUs

## Parameter estimation
```r
est_m <- step_alg(data).  ###parameter estimation
```
## Association analysis/differential abundance testing for each taxon
```r
res.t <- data.pro(est_m,nodes=3) ### Wald test based on a variational sandwich covariance matrix
```
* `data` : a list including observed microbiome data and covariates
* `est_m` : a list of the estimated model parameters and variational parameters
* `nodes` : the required nodes for the parallel 

# Example
The following step shows that how to generate multilevel microbiome count data.  
```r
sub.n <- 10
pos.n <- 3
vis.n <- 3
otu.n <- 100
sim_dat <- sim.data(sub.n, pos.n, vis.n, otu.n)

```

Run the *step_alg* and *data.pro* functions to conduct association analysis for a taxon and it returns a p-value.

```r
Y_mat <- sim_dat$Y_mat
sim_dat["Y_mat"] <- NULL
data <- c(sim_dat, list(Y=Y_mat[,1]))
est_m <- step_alg(data)  ###parameter estimation
res.t <- data.pro(est_m,nodes=3) ### Wald test based on a sandwich covariance structure
```

[1] Tiantian Liu, Peirong Xu, Yueyao Du, Hui Lu, Hongyu Zhao, Tao Wang. (2021) MZINBVA: Variational approximation for multilevel zero-inflated negative-binomial models for association analysis in microbiome surveys.



