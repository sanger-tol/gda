# Module for displaying cluster junction count results in the GDA Shiny app

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

clusterJunctionsUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarPanel(
      helpText("This page shows raster plots that are based on counts of junctions between windows belonging to each UMAP+HDBSCAN cluster. The junctions are the borders between clusters in the 'clusters.bed' file of GDA output"),
      selectInput(ns("selected_variable"), "Variable:",
                  c("Observed counts"="observed_counts",
                    "Expected counts"="expected_counts",
                    "Fold difference"="fold_difference",
                    "p value"="p_value")
      ),
      selectInput(ns("plot_palette"), "Plot palette:",
                  c("Blue-yellow"="blue_yellow",
                    "Red-green"="red_green")),
      sliderInput(ns("plot_width_slider"), "Plot width", min=400, max=3000, value=600, ticks=FALSE),
      sliderInput(ns("plot_height_slider"), "Plot height", min=400, max=2000, value=600, ticks=FALSE),
      selectInput(ns("pvalue_cutoff_selectinput"), "p value cutoff:",
                  c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-10, 1e-20, 1e-50, 1e-100), selected=1e-20),
      selectInput(ns("plot_format_input"), "Plot file format", choices=c("png", "svg")),
      downloadButton(ns("save_plot"), "Save plot"),
    ),
    mainPanel(
      fluidRow(
        loadingMessageUI(ns("loading_message")),
        textOutput(ns("error_message_output"), inline=TRUE),
        plotOutput(ns("tile_plot"), inline=TRUE),
        br(),
        htmlOutput(ns("variable_description"), inline=TRUE)
      )
    ),
  )
}


clusterJunctionsServer <- function(id, server_input_args) {
  moduleServer(
    id, function(input, output, session) {
      python_tuple_to_vector <- function(tuple) {
        # Input: a Python tuple that has been saved as a string e.g. '(0, -1)'
        # Output: tuple as an R vector of numbers
        tuple <- gsub("(", "", tuple, fixed=TRUE)
        tuple <- gsub(")", "", tuple, fixed=TRUE)
        tuple <- gsub(" ", "", tuple, fixed=TRUE)
        split_tuple <- strsplit(tuple, ",", fixed=TRUE)[[1]]
        split_tuple <- as.numeric(split_tuple)
        return(split_tuple)
      }
      
      get_pvalue_star <- function(pvalue, pvalue_cutoff) {
        pvalue_star <- ""
        if(pvalue < pvalue_cutoff) {
          pvalue_star <- "*"
        }
        return(pvalue_star)
      }
      
      
      counts_to_table <- function(df, selected_col, rounding_limits, pvalue_cutoff) {
        # Function for reformatting input data that has been read from a CSV so that it is ready to be plotted with ggplot2
        # Input: 1) data frame with cluster junctions data loaded from a CSV file, 2) name of a column in the data frame, 
        # 3) a vector with the rounding limits for numbers in different plots, 4) p value cutoff for determining statistical significance
        # Output: all vs all comparisons table of cluster junctions for the selected variable, ready to be plotted as a raster plot
        out_mat <- NULL
        for(i in 1:nrow(df)) {
          clusters_tuple <- df$X[i]
          clusters_vect <- python_tuple_to_vector(clusters_tuple)
          
          selected_font <- "plain"
          if(as.numeric(df$corrected_pvalue[i]) < as.numeric(pvalue_cutoff)) {
            selected_font <- "bold.italic"
          }
          
          out_vect <- c(clusters_vect[1], clusters_vect[2], df[[selected_col]][i], selected_font, df$corrected_pvalue[i])
          out_mat <- rbind(out_mat, out_vect)
          
          out_vect <- c(clusters_vect[2], clusters_vect[1], df[[selected_col]][i], selected_font, df$corrected_pvalue[i])
          out_mat <- rbind(out_mat, out_vect)
        }
        out_df <- as.data.frame(out_mat, stringsAsFactors=FALSE)
        colnames(out_df) <- c("x", "y", "value", "pvalue_font", "pv")
        out_df$x <- as.factor(out_df$x)
        out_df$y <- as.factor(out_df$y)
        out_df <- out_df[!duplicated(out_df), ]
        out_df$value <- as.numeric(out_df$value)
        out_df$value_text <- out_df$value
        
        
        if(selected_col == "corrected_pvalue") {
          raw_pvalues <- out_df$value
          out_df$value <- -log2(out_df$value)
          
          pvalue_stars <- sapply(raw_pvalues, get_pvalue_star, pvalue_cutoff=as.numeric(pvalue_cutoff))
          out_df$value_text <- paste0(round(out_df$value, rounding_limits[selected_col]), "\n(", round(raw_pvalues, rounding_limits[selected_col]), pvalue_stars, ")")
          
          ind <- which(is.infinite(out_df$value))
          if(!is.null(ind)) {
            ind2 <- which(!is.infinite(out_df$value))
            out_df$value[ind] <- max(out_df$value[ind2])
          }
        } else {
          out_df$value_text <- round(out_df$value_text, rounding_limits[selected_col])
        }
        return(out_df)
        
      }
      
      
      main <- function(infile_path, selected_species) {
        df <- read.csv(infile_path, sep=",", stringsAsFactors=FALSE)
        rounding_limits <- c(1, 1, 5, 5)
        names(rounding_limits) <- c("observed", "expected_int", "corrected_pvalue", "log2_observed_vs_expected_fold_diff")
        
        counts_table_observed <- reactive({counts_to_table(df, "observed", rounding_limits, input$pvalue_cutoff_selectinput)})
        counts_table_expected_int <- reactive({counts_to_table(df, "expected_int", rounding_limits, input$pvalue_cutoff_selectinput)})
        counts_table_pvalue <- reactive({counts_to_table(df, "corrected_pvalue", rounding_limits, input$pvalue_cutoff_selectinput)})
        counts_table_fold_change <- reactive({counts_to_table(df, "log2_observed_vs_expected_fold_diff", rounding_limits, input$pvalue_cutoff_selectinput)})
        
        plot_titles <- c("Observed counts of junctions between clusters", "Expected counts of junctions between clusters", "Fisher test p-value (observed vs expected junctions)", "Fold difference: log2(observed vs expected junctions)")
        names(plot_titles) <- c("observed_counts", "expected_counts", "p_value", "fold_difference")
        
        legend_titles <- c("Number of junctions", "Number of junctions", "-log2(p value)", "log2(fold difference)")
        names(legend_titles) <- c("observed_counts", "expected_counts", "p_value", "fold_difference")
        
        palette_colours <- list(c("deepskyblue3", "yellow"), c("red", "green"))
        names(palette_colours) <- c("blue_yellow", "red_green")
        
        variable_descriptions <- c("This plot shows the observed counts of junctions between windows belonging to each UMAP+HDBSCAN cluster. Junctions between windows belonging to the same type of cluster are included in the counts.",
                                   "This plot shows the expected counts of junctions between windows belonging to each cluster, if the junctions were distributed randomly (i.e. the number of windows per each cluster would be the same but the windows would appear in a randomly shuffled order). Junctions between windows belonging to the same type of cluster are included in the counts.",
                                   "This plot shows Fisher test p values for comparing whether observed counts of junctions differ from what is expected if they were distributed by chance. The p values have been corrected for multiple hypothesis testing using the BH method.<br>The upper row in each cell shows the -log2 transformed p value, and the lower rows show the untransformed p values, with star symbols that indicate statistical significance.\n",
                                   "This plot shows log2 of the observed counts divided by expected counts for each junction.")
        variable_descriptions <- paste0(variable_descriptions, "<br>Junctions with counts that are significantly different from what is expected by chance (based on Fisher test) are shown in <B><i>bold+italics</i></B> font.")
        names(variable_descriptions) <- c("observed_counts", "expected_counts", "p_value", "fold_difference")
        
        
        cluster_junctions_plot <- reactive({
          
          plotting_df_list <- list(counts_table_observed(), counts_table_expected_int(), counts_table_pvalue(), counts_table_fold_change())
          names(plotting_df_list) <- c("observed_counts", "expected_counts", "p_value", "fold_difference")
          
          ggplot(plotting_df_list[[input$selected_variable]], aes(x, y, fill=value)) + 
            geom_tile() + geom_text(aes(label = value_text, fontface=pvalue_font)) +
            scale_fill_gradient(low=palette_colours[[input$plot_palette]][1], high=palette_colours[[input$plot_palette]][2]) +
            xlab("Cluster") + ylab("Cluster") + ggtitle(plot_titles[input$selected_variable]) +
            theme(plot.title=element_text(hjust=0.5)) +
            labs(fill=legend_titles[input$selected_variable])
        })
        
        
        observe({
          output$tile_plot <- renderPlot({
            cluster_junctions_plot()
          }, height=input$plot_height_slider, width=input$plot_width_slider)
          
          output$variable_description <- renderText({
            variable_descriptions[input$selected_variable]
          })
        })
        
        
        output$save_plot <- downloadHandler(
          filename = function() {
            plot_title <- paste0(format(Sys.Date(), "%Y%m%d"), "_", selected_species, "_cluster_junctions_", input$selected_variable)  
            paste0(plot_title, ".", input$plot_format_input)
          },
          content = function(file){
            req(cluster_junctions_plot())
            width_value <- round((input$plot_width_slider/400)*7, 2)
            height_value <- round((input$plot_height_slider/input$plot_width_slider)*width_value, 2)
            ggsave(file, plot=cluster_junctions_plot(), device=input$plot_format_input, width=width_value, height=height_value, units="in")
          }
        )
      }
      
      in_folder <- server_input_args[1]
      selected_species <- server_input_args[2]
      
      infile_path <- paste0(in_folder, "/", selected_species, "/cluster_junctions_fisher_test.csv")
      
      if(file.exists(infile_path)) {
        main(infile_path, selected_species)
      } else {
        output$error_message_output <- renderText(paste0("Cannot display figures in this tab, as the input file for cluster junction counts was not found in the expected location (", infile_path, ")"))
      }
    }
  )
}
