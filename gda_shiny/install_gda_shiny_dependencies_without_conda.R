#!/usr/bin/env Rscript
# Script for installing the dependencies for the GDA Shiny app without using conda

list_of_packages <- c("devtools", "shiny", "ggplot2", "svglite", "gplots", "rjson", "reshape2", "gridExtra", "scales")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) {
  for(package in new_packages) {
    print(package)
    if(package=="devtools") {
      install.packages("devtools", dependencies=TRUE)
    } else if(package=="svglite") {
      devtools::install_github("r-lib/svglite")
    } else if(package=="gridExtra") {
      library("devtools")
      install_version("gridExtra", version="2.3", dependencies=TRUE)
    } else {
      install.packages(package, dependencies=TRUE)
    }
  }
} else {
  print("All required packages appear to be already installed")
}