# jlst
This package is used to fit joint location-and-scale tests. This package is currently under development at the University of Bristol.

## Functions
* jlssc - joint location-and-scale score test. 
* jlsp - joint location-and-scale test using Fisher's method.  
* vartest - variability tests with Breusch-Pagan or Brown-Forsythe methods.  

## Installation
1. install.packages("devtools")
2. library(devtools)
3. install_github("jrs95/jlst")
4. library(jlst)  

## Example
\# Data  
x <- rbinom(1000, 1, 0.5)  
y <- 0.5 + rnorm(1000, 0.025, 0.025)*x + rnorm(1000, 0, 0.1)  

\# Variance test  
vartest(y, x=as.factor(x), type=1) # Breusch-Pagan test  
vartest(y, x=as.factor(x), type=2) # Brown-Forsythe test  

\# Joint location-and-scale test using Fisher's method   
jlsp(y, x=as.factor(x), var.type=1) # Breusch-Pagan variance test  
jlsp(y, x=as.factor(x), var.type=2) # Brown-Forsythe variance test  

\# Joint location-and-scale score test   
jlssc(y, x=as.factor(x), type=1) # Breusch-Pagan variance test  
jlssc(y, x=as.factor(x), type=2) # Brown-Forsythe variance test  
jlssc(y, x=as.factor(x), type=3) # Method of moments version of the test with the Breusch-Pagan variance test  
jlssc(y, x=as.factor(x), type=4) # Method of moments version of the test with the Brown-Forsythe variance test   

## Citation
Staley JR, Windmeijer F, et al. A robust mean and variance test with application to epigenome-wide association studies. TBC.
