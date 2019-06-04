# jlst
This package is used to fit joint location-and-scale tests. This package is currently under development at the University of Bristol.

## Functions
* jlssc - joint location-and-scale score test. 
* jlsp - joint location-and-scale test using Fisher's method.  
* vartest - variability tests with Breusch-Pagan or Brown-Forsythe methods.  

## Installation
1. install.packages("devtools")
2. library(devtools)
3. install_github("jrs95/jlst", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
4. library(hyprcoloc)
5. browseVignettes("jlst")  

