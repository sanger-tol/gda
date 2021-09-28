# Module for the UMAP plots of the GDA Shiny app

umapPlotUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      
      sidebarPanel(
        helpText("This page shows UMAP+HDBSCAN plot of clusters of genomic windows. Each dot in the plot is a genomic window. Cluster -1 is the noise cluster (unclassified windows)"),
        sliderInput(ns("plot_width_slider"), "Plot width", min=400, max=3000, value=600, ticks=FALSE),
        sliderInput(ns("plot_height_slider"), "Plot height", min=400, max=2000, value=600, ticks=FALSE),
        sliderInput(ns("umap_scatterplot_point_size"), "Point size:", min=0.0125/2, max=3, value=0.5, ticks=FALSE),
        selectInput(ns("plot_format_input"), "Plot file format", choices=c("png", "svg")),
        downloadButton(ns("save_hdbscan_plot"), "Save HDBSCAN clusters plot"),
        uiOutput(ns("per_species_plot_download_ui")),
        loadingMessageUI(ns("loading_message"))
      ),
      mainPanel(
        plotOutput(ns("umap_plot"), inline=TRUE),
        uiOutput(ns("umap_plot_by_species_ui"), inline=TRUE),
      )
    )
    
  )
}



umapPlotServer <- function(id, umap_server_args) {
  moduleServer(
    id, function(input, output, session) {
      in_folder <- umap_server_args[1]
      selected_species <- umap_server_args[2]
      umap_clusters_path <- paste0(in_folder, "/umap_clustering.csv")
      
      umap_df <- read.csv(umap_clusters_path, stringsAsFactors=FALSE)
      species_count <- reactive({length(unique(umap_df$species))})
      umap_palette <- get_clusters_palette(umap_df)
      
      umap_clusters_plot <- reactive({ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=as.factor(cluster))) + geom_point(size=input$umap_scatterplot_point_size) + scale_color_manual(values=umap_palette, name="Clusters") + guides(colour=guide_legend(override.aes=list(size=7))) +
          ggtitle("UMAP + HDBSCAN genomic window clusters") + theme(plot.title=element_text(hjust=0.5))})
      
      umap_species_plot <- reactive({ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=as.factor(species))) + geom_point(size=input$umap_scatterplot_point_size) + guides(colour=guide_legend(title="Species", override.aes=list(size=7))) +
          ggtitle("Genomic window clusters per species") + theme(plot.title=element_text(hjust=0.5))})
      
      observe({
        output$umap_plot <- renderPlot({
          umap_clusters_plot()
        }, height=input$plot_height_slider, width=input$plot_width_slider)
        
        output$umap_plot_by_species <- renderPlot({
          umap_species_plot()
        }, height=input$plot_height_slider, width=input$plot_width_slider)
      })
      
      output$save_hdbscan_plot <- downloadHandler(
        filename = function() {
          plot_title <- paste0(format(Sys.Date(), "%Y%m%d"), "_", selected_species, "_umap_hdbscan_clusters")  
          paste0(plot_title, ".", input$plot_format_input)
        },
        content = function(file){
          req(umap_clusters_plot())
          width_value <- round((input$plot_width_slider/400)*7, 2)
          height_value <- round((input$plot_height_slider/input$plot_width_slider)*width_value, 2)
          ggsave(file, plot=umap_clusters_plot(), device=input$plot_format_input, width=width_value, height=height_value, units="in")
        }
      )
      
      
      output$save_per_species_plot <- downloadHandler(
        filename = function() {
          plot_title <- paste0(format(Sys.Date(), "%Y%m%d"), "_", selected_species, "_umap_clusters_per_species")  
          paste0(plot_title, ".", input$plot_format_input)
        },
        content = function(file){
          req(umap_clusters_plot())
          width_value <- round((input$plot_width_slider/400)*7, 2)
          height_value <- round((input$plot_height_slider/input$plot_width_slider)*width_value, 2)
          ggsave(file, plot=umap_species_plot(), device=input$plot_format_input, width=width_value, height=height_value, units="in")
        }
      )
      
      observe({
        if(species_count()>1) {
          output$per_species_plot_download_ui <- renderUI({
            downloadButton(session$ns("save_per_species_plot"), "Save clusters per species plot")
          })
          
          output$umap_plot_by_species_ui <- renderUI({
            plotOutput(session$ns("umap_plot_by_species"))
          })
        }
      })
    }
  )    
}