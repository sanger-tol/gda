# Module for the features heat map of the GDA Shiny app

HEATMAP_DEFAULT_WIDTH <- 500

heatmapUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        uiOutput(ns("heatmap_helptext_output")),
        br(),
        sliderInput(ns("plot_width_slider"), "Plot width", min=HEATMAP_DEFAULT_WIDTH, max=2000, value=1200, ticks=FALSE),
        sliderInput(ns("plot_height_slider"), "Plot height", min=600, max=6000, value=800, ticks=FALSE),
        sliderInput(ns("row_labels_size_slider"), "Row labels size", min=0.25, max=3, value=1, ticks=FALSE),
        sliderInput(ns("column_labels_size_slider"), "Column labels size", min=0.25, max=3, value=1, ticks=FALSE),
        sliderInput(ns("legend_height_slider"), "Colour key height", min=0.25/2, max=0.5, value=0.2, ticks=FALSE),
        selectInput(ns("heatmap_palette_input"), "Heatmap palette", choices=c("red-yellow-green", "blue-yellow")),
        selectInput(ns("plot_format_input"), "Plot file format", choices=c("png", "svg")),
        downloadButton(ns("save_plot"), "Save plot"),
      ),
      mainPanel(
        loadingMessageUI(ns("loading_message")),
        plotOutput(ns("clusters_heatmap"), inline=TRUE),
      )
    )
    
  )
}



heatmapServer <- function(id, in_folder, selected_species, heatmap_data_type) {
  moduleServer(
    id, function(input, output, session) {
      
      get_heatmap_palette <- function(palette_name) {
        # Takes a palette name from the dropdown menu as an input and returns the corresponding palette
        out_palette <- NULL
        if(palette_name == "red-yellow-green") {
          out_palette <- colorRampPalette(c("red4", "red3", "red2", "orange", "yellow", "yellowgreen", "green2", "green3", "green4"))(n=100)
        } else if(palette_name == "blue-yellow") {
          out_palette <- colorRampPalette(c("blue", "yellow2"))(n=100)
        }
        return(out_palette)
      }
      
      draw_heatmap <- function(heatmap_mat, heatmap_colours, heatmap_title) {
        # Draws the heatmap
        heatmap.2(heatmap_mat,
                  main = heatmap_title,
                  col=heatmap_colours,
                  density.info="none",
                  trace="none",
                  cexRow=input$row_labels_size_slider,
                  cexCol=input$column_labels_size_slider,
                  margins=c(20, 30),
                  lhei=c(input$legend_height_slider, 1), lwid=c(1, 5))
      }
      
      crop_cluster_name <- function(x) {
        # Function for removing X characters from cluster names. The X characters were automatically added by R when loading the CSV
        x <- gsub("X.", "-", x, fixed=TRUE)
        x <- gsub("X", "", x, fixed=TRUE)
        x <- paste0("Cluster ", x)
        return(x)
      }
      
      heatmap_title <- NULL
      outfile_suffix <- NULL
      heatmap_mat <- NULL
      heatmap_helptext <- NULL
      heatmap_helptext2 <- NULL
      
      if(heatmap_data_type == "genomic_features") {
        heatmap_title <- "Features enriched in clusters"
        outfile_suffix <- "features_enriched_in_clusters"
        cluster_heatmap_data_path <- paste0(in_folder, "/cluster_heatmap.csv")
        cluster_heatmap_df <- read.csv(cluster_heatmap_data_path, stringsAsFactors=FALSE)
        rownames(cluster_heatmap_df) <- cluster_heatmap_df$X
        rownames(cluster_heatmap_df) <- paste0("Cluster ", rownames(cluster_heatmap_df))
        cluster_heatmap_df <- cluster_heatmap_df[,2:ncol(cluster_heatmap_df)]
        heatmap_mat <- t(as.matrix(cluster_heatmap_df))
        heatmap_helptext <- "This page shows which genomic features are enriched in each UMAP+HDBSCAN cluster (based on Kolmogorov-Smirnov test)."
        heatmap_helptext2 <- "If there is a high number of variables, you may need to increase plot height to display the names of all variables properly"
      } else if(heatmap_data_type == "cluster_locations") {
        heatmap_title <- "Chromosome cluster composition"
        outfile_suffix <- "chromosome_cluster_composition"
        chrom_cluster_heatmap_data_path <- paste0(in_folder, "/", selected_species, "/chrcompheat.csv")
        chrom_cluster_heatmap_df <- read.csv(chrom_cluster_heatmap_data_path, stringsAsFactors=FALSE)
        rownames(chrom_cluster_heatmap_df) <- chrom_cluster_heatmap_df$X
        colnames(chrom_cluster_heatmap_df) <- sapply(colnames(chrom_cluster_heatmap_df), crop_cluster_name)
        chrom_cluster_heatmap_df <- chrom_cluster_heatmap_df[,2:ncol(chrom_cluster_heatmap_df)]
        heatmap_mat <- as.matrix(chrom_cluster_heatmap_df)
        heatmap_helptext <- "This page shows a heat map of the distribution of UMAP+HDBSCAN clusters per scaffold."
        heatmap_helptext2 <- "If there is a high number of scaffolds, you may need to increase plot height to display the names of all variables properly"
      }
      
      heatmap_palette <- reactive({get_heatmap_palette(input$heatmap_palette_input)})
    
      observe({
        output$clusters_heatmap <- renderPlot({
          draw_heatmap(heatmap_mat, heatmap_palette(), heatmap_title)
        }, height=input$plot_height_slider, input$plot_width_slider)
      })
      
      output$heatmap_helptext_output <- renderUI({
        helpText(heatmap_helptext, br(), br(), heatmap_helptext2)
      })
      
      
      output$save_plot <- downloadHandler(
        filename = function() {
          plot_title <- paste0(format(Sys.Date(), "%Y%m%d"), "_", selected_species, "_", outfile_suffix)
          plot_title <- paste0(plot_title, ".", input$plot_format_input)
        },
        content = function(file){
          isolate({
            width_value <- round((input$plot_width_slider/HEATMAP_DEFAULT_WIDTH)*7, 2)
            height_value <- round((input$plot_height_slider/input$plot_width_slider)*width_value, 2)
            if(input$plot_format_input == "png") {
              png(file, width=width_value, height=height_value, units="in", res=350)
            } else if (input$plot_format_input == "svg") {
              svg(file, width=width_value, height=height_value)
            }
            draw_heatmap(heatmap_mat, heatmap_palette(), heatmap_title)
            dev.off()
          })
          
        }
      )
      
    }
  )    
}