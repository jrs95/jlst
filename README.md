# jlst <img src='man/figures/logo.png' align="right" height="139"/>

This package is used to fit joint mean (location) and variance (scale) tests, i.e. joint location-and-scale tests. The package also has functions to perform variability tests using the Breusch-Pagan and Brown-Forsythe methods.  

## Functions
* `jlssc`: joint location-and-scale score test. 
* `jlsp`: joint location-and-scale test using Fisher's method.  
* `vartest`: variability tests with Breusch-Pagan or Brown-Forsythe methods.  

## Installation
```
install.packages("remotes")
remotes::install_github("jrs95/jlst")
```

## Example
```
# Libraries
library(jlst)

# Data  
x <- rbinom(1000, 1, 0.5)
y <- 0.5 + 0.025 * x + rnorm(1000, 0, sqrt(0.005 * x)) + rnorm(1000, 0, 0.1)

# Variance test  
vartest(y, x = as.factor(x), type = 1) # Breusch-Pagan test
vartest(y, x = as.factor(x), type = 2) # Brown-Forsythe test

# Joint location-and-scale test using Fisher's method
jlsp(y, x = as.factor(x), var.type = 1) # Breusch-Pagan variance test
jlsp(y, x = as.factor(x), var.type = 2) # Brown-Forsythe variance test

# Joint location-and-scale score test   
jlssc(y, x = as.factor(x), type = 1) # Breusch-Pagan (BP) variance test
jlssc(y, x = as.factor(x), type = 2) # Brown-Forsythe (BF) variance test
jlssc(y, x = as.factor(x), type = 3) # Method of moments version of the test with the BP variance test
jlssc(y, x = as.factor(x), type = 4) # Method of moments version of the test with the BF variance test
```

## Citation
Staley JR, Windmeijer F, *et al.* A robust mean and variance test with application to high-dimensional phenotypes. [Eur J Epidemiol](https://pubmed.ncbi.nlm.nih.gov/34651232/) 2022;37(4):377-387.
