
#######################################################

###### Building the package and documentation #########

#######################################################

# This script sets up a run of the post-processing analysis 

# installing the R-package:

rm(list = ls())

setwd("~/pkg/paper/PostprocessingSST/")
options(max.print = 1e3)

#install.packages('devtools')
library(devtools)
devtools::load_all()
roxygen2::roxygenise()