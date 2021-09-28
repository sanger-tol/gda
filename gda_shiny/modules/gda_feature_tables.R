# Module for the features table of the GDA Shiny app

featureTablesUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    loadingMessageUI(ns("loading_message")),
    helpText("This table shows the features that were found enriched in each UMAP+HDBSCAN cluster of genomic window (based on Kolmogorov-Smirnov test)"),
    uiOutput(ns("datatables"))
  )
}


featureTablesServer <- function(id, feattable_server_args) {
  moduleServer(
    id, function(input, output, session) {

      in_folder <- feattable_server_args[1]
      selected_species <- feattable_server_args[2]
      
      generate_genomeprops_sentence <- function(cluster, genomeprops_df) {
        # Generates a sentence that says what percentage of windows is in a specified cluster
        species_name <- genomeprops_df$X
        cluster_name <- NA
        if(cluster == -1) {
          cluster_name <- "X.1"
        } else {
          cluster_name <- paste0("X", cluster)
        }
        cluster_percentage <- round(as.numeric(genomeprops_df[[cluster_name]]), 2)
        genomeprops_sentence <- paste0(cluster_percentage, "% of the ", species_name," genome is in cluster ", cluster)
        return(genomeprops_sentence)
      }
      
      feattable_path <- paste0(in_folder, "/feattable.csv")
      genomeprops_path <- paste0(in_folder, "/genomeprops.csv")
      
      df <- read.csv(feattable_path, stringsAsFactors=FALSE)
      genomeprops_df <- read.csv(genomeprops_path, stringsAsFactors=FALSE, header=TRUE)
      ind <- which(genomeprops_df$X == selected_species)
      if(length(ind) == 1) {
        genomeprops_df <- genomeprops_df[ind,]
      } else {
        print(paste0("Failed to parse the file of genomic proportions of clusters (",genomeprops_path, ")"))
      }
      
      clusters <- unique(df$cluster)
      
      output$datatables <- renderUI({
        lapply(1:length(clusters), function(i) {
          tagList(
            h3(paste0("Cluster ", clusters[i])),
            p(generate_genomeprops_sentence(clusters[i], genomeprops_df)),
            br(),
            renderDataTable(df[df$cluster == clusters[i],]),
            div(),
          )
        })
      })

    }
  )    
}