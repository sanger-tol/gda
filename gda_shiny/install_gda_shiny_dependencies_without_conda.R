#!/usr/bin/env Rscript
# Script for installing the dependencies for the GDA Shiny app without using conda

# MIT License
# 
# Copyright (c) 2020-2021 Genome Research Ltd.
# 
# Author: Eerik Aunin (ea10@sanger.ac.uk)
# 
# This file is a part of the Genome Decomposition Analysis (GDA) pipeline.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
