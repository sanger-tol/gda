# Module for the cluster position line charts the GDA Shiny app

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

LINE_CHART_DEFAULT_WIDTH <- 1000

chromClusterPositionsUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        helpText("This page shows line charts of the relative positions of UMAP+HDBSCAN clusters in assembly scaffolds. The lengths of scaffolds have been normalised to be equal for this plot. The plot can be used for finding UMAP clusters that tend to be localised at particular regions of scaffolds (e.g. always at scaffold edges)."),
        sliderInput(ns("plot_width_slider"), "Plot width", min=1000, max=2000, value=LINE_CHART_DEFAULT_WIDTH, ticks=FALSE),
          sliderInput(ns("plot_height_slider"), "Plot height", min=1000, max=5000, value=2500, ticks=FALSE),
        uiOutput(ns("scaff_length_slider_ui")),
        selectInput(ns("plot_format_input"), "Plot file format", choices=c("png", "svg")),
        downloadButton(ns("save"), "Save plot")
      ),
      
      mainPanel(
        loadingMessageUI(ns("loading_message")),
        plotOutput(ns("cluster_positions_plot"))
      )
    )
    
  )
}



chromClusterPositionsServer <- function(id, server_input_args) {
  moduleServer(
    id, function(input, output, session) {
      
      get_plot_output_list <- function(cluster_df_list) {
        # Produces a list where each element is a line chart for one UMAP+HDBSCAN cluster
        
        input_n <- length(cluster_df_list)
        
        clusters_sequence <- NA
        
        if(is.element("cluster_-1", names(cluster_df_list))) {
          clusters_sequence <- seq(input_n-2,from=-1,by=1)
        } else {
          clusters_sequence <- seq(input_n-1,from=0,by=1)
        }

        plot_output_list <- lapply(clusters_sequence, function(i) {
          cluster_name <- paste0("cluster_", i)
          plotname <- paste("plot", i, sep="")
          p <- ggplot(cluster_df_list[[cluster_name]], aes(x=pos, y=abundance, group=scaff)) +
            geom_line(aes(color=scaff)) + ggtitle(paste0("Cluster ", i)) + xlab("Location relative to scaffold ends") + ylab("Cluster abundance") + 
            guides(color=guide_legend(title="Scaffold"))
        })
        
        return(plot_output_list)
      }
      
      
      load_cluster_df_list <- function(in_folder, selected_species) {
        # Reads cluster positions data from a json file and converts it to a list of data frames where each data frame contains the relative coordinates of one cluster
        clusterpos_json_path <- paste0(in_folder, "/", selected_species, "/clusterpos.json")
        clusterpos_json_data <- fromJSON(file=clusterpos_json_path)

        cluster_df_list <- list()
        for(cluster_name in names(clusterpos_json_data)) {
          cluster_data <- clusterpos_json_data[[cluster_name]]

          n_obs <- sapply(cluster_data, length)
          seq_max <- seq_len(max(n_obs))
          
          cluster_mat <- t(sapply(cluster_data, "[", i = seq_max)) #https://stackoverflow.com/questions/15201305/how-to-convert-a-list-consisting-of-vector-of-different-lengths-to-a-usable-data
          cluster_mat[is.na(cluster_mat)] <- 0
          
          cluster_df <- melt(cluster_mat)
          colnames(cluster_df) <- c("scaff", "pos", "abundance")

          cluster_df$pos <- as.numeric(cluster_df$pos)
          cluster_df$abundance <- as.numeric(cluster_df$abundance)
          
          cluster_df_list <- append(cluster_df_list, list(cluster_df))
        }
        names(cluster_df_list) <- paste0("cluster_", names(clusterpos_json_data))
        return(cluster_df_list)
      }
      
      
      get_scaff_lengths <- function(in_folder, selected_species) {
        # Extracts the information on scaffold lengths from input json file
        circos_json_path <- paste0(in_folder, "/", selected_species, "/circos.json")
        circos_graph_data <- fromJSON(file=circos_json_path)
        
        genome_df <- circos_graph_data$genome
        genome_df <- t(as.data.frame(matrix(unlist(genome_df), nrow=length(unlist(genome_df[1])))))
        colnames(genome_df) <- c("id", "label", "color", "len")
        genome_df <- as.data.frame(genome_df, stringsAsFactors=FALSE)
        genome_df$len <- as.numeric(genome_df$len)
        return(genome_df)
      }
      
    
      filter_data_by_scaff_len <- function(cluster_df_list, genome_df, len_cutoff) {
        # Filters cluster coordinates data and keeps only scaffolds with length above a cutoff
        ind <- which(genome_df$len < len_cutoff)
        too_short_scaffs <- genome_df$id[ind]
        filtered_cluster_df_list <- list()
        for(df in cluster_df_list) {
          ind <- which(is.element(df$scaff, too_short_scaffs))
          if(length(ind) > 0) {
            df <- df[-ind,]
          }
          filtered_cluster_df_list <- append(filtered_cluster_df_list, list(df))
        }
        names(filtered_cluster_df_list) <- names(cluster_df_list)
        return(filtered_cluster_df_list)
      }
      
      in_folder <- server_input_args[1]
      selected_species <- server_input_args[2]
      
      cluster_df_list <- load_cluster_df_list(in_folder, selected_species)
      genome_df <- get_scaff_lengths(in_folder, selected_species)
      max_scaff_len <- max(genome_df$len)
      
      
      observe({
        output$cluster_positions_plot <- renderPlot({
          do.call("grid.arrange", c(get_plot_output_list(filter_data_by_scaff_len(cluster_df_list, genome_df, input$scaff_length_slider)), ncol=1))
        }, height=input$plot_height_slider, width=input$plot_width_slider)
      })
      
      
      output$scaff_length_slider_ui <- renderUI({
        sliderInput(session$ns("scaff_length_slider"), "Minimum length filter for displaying scaffolds (bp)", min=0, max=max_scaff_len-1, round=2, ticks=FALSE, value=c(0))
      })
      
      
      output$save <- downloadHandler(
        filename = function() {
          plot_title <- paste0(format(Sys.Date(), "%Y%m%d"), "_", selected_species, "_cluster_position_line_charts")
          paste0(plot_title, ".", input$plot_format_input)
        },
        content = function(file) {
          width_value <- round((input$plot_width_slider/LINE_CHART_DEFAULT_WIDTH)*14, 2)
          height_value <- round((input$plot_height_slider/input$plot_width_slider)*width_value, 2)
          ggsave(file, plot=do.call("grid.arrange", c(get_plot_output_list(filter_data_by_scaff_len(cluster_df_list, genome_df, input$scaff_length_slider)), ncol=1)), device=input$plot_format_input, width=width_value, height=height_value, units="in", limitsize=FALSE)
        }
      )
      
    }
  )    
}
