# Module for the cluster locations plots of the GDA Shiny app

clusterLocationsRasterPlotUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      div(class="set1",
        sidebarPanel(
          helpText("This page shows a raster plot of the locations of UMAP+HDBSCAN clusters in assembly scaffolds.", br(), "x axis: scaffold coordinate.", br(), "y axis: scaffolds, sorted by size.", br(), "Colours: UMAP+HDBSCAN clusters. The colours for each cluster are the same as in the UMAP scatter plots"),
          sliderInput(ns("plot_width_slider"), "Plot width", min=400, max=2000, value=600, ticks=FALSE),
          sliderInput(ns("plot_height_slider"), "Plot height", min=600, max=3000, value=800, ticks=FALSE),
          uiOutput(ns("scaff_length_slider_ui")),
          uiOutput(ns("scaff_coord_slider_ui")),
          selectInput(ns("plot_format_input"), "Plot file format", choices=c("png", "svg")),
          downloadButton(ns("save"), "Save plot"),
          loadingMessageUI(ns("loading_message")),
        )
      ),
      
      mainPanel(
          plotOutput(ns("raster_plot"), inline = TRUE),
          br()
        )
    )
    
  )
}



clusterLocationsRasterPlotServer <- function(id, rasterplot_server_args) {
  moduleServer(
    id, function(input, output, session) {
    
      load_input_data <- function(in_folder, selected_species) {
        # Loads the data of cluster locations from a JSON file and converts it to a format that is suitable for plotting as a geom_tile() plot with ggplot2
        circos_json_path <- paste0(in_folder, "/", selected_species, "/circos.json")
        circos_graph_data <- fromJSON(file=circos_json_path)
        clusters_df <- circos_graph_data$clusters
        clusters_df <- t(as.data.frame(matrix(unlist(clusters_df), nrow=length(unlist(clusters_df[1])))))
        clusters_df <- as.data.frame(clusters_df, stringsAsFactors=FALSE)
        clusters_df <- clusters_df[2:6]
        colnames(clusters_df) <- c("chr", "start", "end", "cluster", "color")
        clusters_df$start <- as.numeric(clusters_df$start)
        clusters_df$end <- as.numeric(clusters_df$end)
        clusters_df$cluster <- as.numeric(clusters_df$cluster)
        
        genome_df <- circos_graph_data$genome
        genome_df <- t(as.data.frame(matrix(unlist(genome_df), nrow=length(unlist(genome_df[1])))))
        colnames(genome_df) <- c("id", "label", "color", "len")
        genome_df <- as.data.frame(genome_df, stringsAsFactors=FALSE)
        genome_df$len <- as.numeric(genome_df$len)
        
        
        clusters_df$scaff_size <- NA
        
        for(i in 1:nrow(genome_df)) {
          scaff_name <- genome_df$id[i]
          scaff_size <- genome_df$len[i]
          ind <- which(clusters_df$chr == scaff_name)
          if(length(ind) > 0){
            clusters_df$scaff_size[ind] <- scaff_size
          }
        }
        
        df2 <- clusters_df[order(clusters_df$scaff_size, clusters_df$chr, -clusters_df$start, decreasing = TRUE),]
        
        df2$chr <- factor(df2$chr, ordered=TRUE, levels=rev(unique(df2$chr)))
        df2$cluster <- as.factor(df2$cluster)
        return(df2)
      }
      
      in_folder <- rasterplot_server_args[1]
      selected_species <- rasterplot_server_args[2]
      
      df2 <- load_input_data(in_folder, selected_species)
      clusters_palette <- get_clusters_palette(df2)

      slider_scaler <- 1000000
      
      slider_max <- ceiling((max(df2$scaff_size) + 1)/slider_scaler)
      slider_min <- min(df2$scaff_size)/slider_scaler
      
      filtered_df2 <- reactive({df2[df2$scaff_size >= slider_scaler*input$scaff_length_slider[1] & df2$scaff_size <= slider_scaler*input$scaff_length_slider[2] & df2$start >= slider_scaler*input$scaff_coord_slider[1] & df2$end <= slider_scaler*input$scaff_coord_slider[2], ]})
      output$scaff_length_slider_ui <- renderUI({
        sliderInput(session$ns("scaff_length_slider"), "Length range filter for displaying scaffolds (Mb)", min=slider_min, max=slider_max, round=2, ticks=FALSE, value=c(slider_min, slider_max))
      })
      
      observe({
        if(!is.null(input$scaff_length_slider[2])) {
            output$scaff_coord_slider_ui <- renderUI({
              sliderInput(session$ns("scaff_coord_slider"), "Coordinate range filter for displaying scaffolds (Mb)", min=0, max=input$scaff_length_slider[2], round=2, ticks=FALSE, value=c(0, input$scaff_length_slider[2]))
          })
        }
      })
      
      cluster_locations_rasterplot <- reactive({
        p <- ggplot(data = filtered_df2(), aes(x = start, y = chr, fill = cluster)) +
          geom_tile() + scale_fill_manual(values=clusters_palette, drop=TRUE, name="UMAP cluster", limits=levels(df2$cluster)) + xlab("Scaffold coordinate (bp)") + ylab("Scaffold") +
        ggtitle("Raster plot of locations of UMAP clusters on scaffolds") + theme(plot.title=element_text(hjust=0.5)) + theme(axis.text.x = element_text(angle = 90)) +
          scale_x_continuous(labels=comma)
        scaff_count <- length(unique(filtered_df2()$chr))
        if(scaff_count >= 100) {
          p <- p + theme(axis.text.y=element_blank())
        }
        p
      })
      
      observe({
        output$raster_plot <- renderPlot({
          cluster_locations_rasterplot()
        }, height=input$plot_height_slider, width=input$plot_width_slider)
      })
      
      output$save <- downloadHandler(
        filename = function() {
          plot_title <- paste0(format(Sys.Date(), "%Y%m%d"), "_", selected_species, "_cluster_locations")  
          paste0(plot_title, ".", input$plot_format_input)
        },
        content = function(file){
          req(cluster_locations_rasterplot())
          width_value <- round((input$plot_width_slider/600)*7, 2)
          height_value <- round((input$plot_height_slider/input$plot_width_slider)*width_value, 2)
          ggsave(file, plot = cluster_locations_rasterplot(), device=input$plot_format_input, width=width_value, height=height_value, units="in")
        }
      )
      
    }
  )  
}