# File for functions that are shared between modules of the GDA Shiny app

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