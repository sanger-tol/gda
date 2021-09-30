# File for functions that are shared between modules of the GDA Shiny app

# MIT License
# 
# Copyright (c) 2020-2021 Genome Research Ltd.
# 
# Authors: Eerik Aunin (ea10@sanger.ac.uk) and Adam Reid (ar11@sanger.ac.uk)
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

get_clusters_palette <- function(df) {
  # Extracts colours for the UMAP scatter plot from the table that has been loaded from a CSV file
  colours_df <- df[,4:5]
  colours_df <- unique(colours_df)
  colours_df <- colours_df[order(colours_df$cluster),]
  clusters_palette <- colours_df$color
  return(clusters_palette)
}

txt_out <- function(input) {
  # _Function to print the name and the value of a variable to the console
  call_f <- match.call()
  call_f_list <- as.list(call_f)
  print(call_f_list[[-1]])
  print(input)
}
