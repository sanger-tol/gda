#!/usr/bin/env Rscript
# Main file of the GDA Shiny app

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

library("shiny")
library("ggplot2")
library("svglite")
library("gplots", warn.conflicts=FALSE)
library("rjson")
library("reshape2")
library("gridExtra")
library("scales")


check_input_files <- function(in_folder, selected_species) {
  # Function for checking if the required input files are present
  files_ok_flag <- TRUE
  if(in_folder == "" || selected_species == "") {
    files_ok_flag <- FALSE
    print("Please start this Shiny app by using command line options to select input folder and the species that will be analysed")
  } else {
    selected_species_folder <- paste0(in_folder, "/", selected_species)
    
    if(!dir.exists(in_folder)) {
      files_ok_flag <- FALSE
      print(paste0("Folder not found: ", in_folder))
    } else if(!dir.exists(selected_species_folder)) {
      files_ok_flag <- FALSE
      print(paste0("Folder not found: ", selected_species_folder))
    }
    
    clusterpos_json_file <- paste0(selected_species_folder, "/", "clusterpos.json")
    circos_json_file <- paste0(selected_species_folder, "/", "circos.json")
    feattable_path <- paste0(in_folder, "/feattable.csv")
    genomeprops_path <- paste0(in_folder, "/genomeprops.csv")
    umap_clusters_path <- paste0(in_folder, "/umap_clustering.csv")
    cluster_heatmap_data_path <- paste0(in_folder, "/cluster_heatmap.csv")
    paths_checklist <- c(clusterpos_json_file, circos_json_file, feattable_path, genomeprops_path, umap_clusters_path, cluster_heatmap_data_path)
    for(path in paths_checklist) {
      if(!file.exists(path)) {
        files_ok_flag <- FALSE
        print(paste0("File not found: ", path))
        break
      }
    }
  }
  
  stopifnot(files_ok_flag == TRUE)
}

gridextra_version <- packageVersion("gridExtra")
if((gridextra_version) < 2.3) {
  stop(paste0("gridExtra package version >= 2.3 is required for the GDA Shiny app but the version currently installed is ", gridextra_version, ". Please update gridExtra"))
}

in_folder <- getShinyOption("in_folder", "")
selected_species <- getShinyOption("selected_species", "")

check_input_files(in_folder, selected_species)


source("./modules/gda_umap_plots.R")
source("./modules/gda_heatmap.R")
source("./modules/gda_cluster_position_line_charts.R")
source("./modules/gda_cluster_locations.R")
source("./modules/gda_feature_tables.R")
source("./modules/gda_cluster_junctions.R")
source("./modules/gda_loading_message.R")
source("./modules/gda_shiny_shared_functions.R")


ui <- navbarPage("GDA",
  tabPanel("UMAP plots",
    fluidPage(
      umapPlotUI("umap1")
    )
          
 ),
 tabPanel("Cluster locations",
    fluidPage(
      clusterLocationsRasterPlotUI("rasterplot1"),
    )
 ),
 tabPanel("Cluster heatmaps",
    fluidPage(
      heatmapUI("features_heatmap1")
    )
  ),
 tabPanel("Feature tables",
    fluidPage(
      featureTablesUI("feature_tables1")
    )
  ),
 tabPanel("Cluster positions across chromosomes", 
    fluidPage(
      chromClusterPositionsUI("chrom_cluster_positions1")
    )),
 tabPanel("Chromosome cluster composition", 
    fluidPage(
      heatmapUI("chrom_cluster_composition_heatmap1")
    )),
 tabPanel("Cluster junction counts", 
  fluidPage(
    clusterJunctionsUI("cluster_junction_counts1")
  ))
)

check_input_files(in_folder, selected_species)

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  
  umapPlotServer("umap1", c(in_folder, selected_species))
  heatmapServer("features_heatmap1", in_folder, selected_species, "genomic_features")
  heatmapServer("chrom_cluster_composition_heatmap1", in_folder, selected_species, "cluster_locations")
  chromClusterPositionsServer("chrom_cluster_positions1", c(in_folder, selected_species))
  clusterLocationsRasterPlotServer("rasterplot1",  c(in_folder, selected_species)) 
  featureTablesServer("feature_tables1", c(in_folder, selected_species))
  clusterJunctionsServer("cluster_junction_counts1", c(in_folder, selected_species))
}

shinyApp(ui, server)
