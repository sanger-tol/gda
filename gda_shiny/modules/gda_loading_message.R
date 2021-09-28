# Module for displaying the 'Loading...' message
# based on: https://stackoverflow.com/questions/66825095/how-to-disconnect-the-loading-massage-to-the-content-when-the-app-is-fully-loade

loadingMessageUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$head(tags$style(type="text/css", "
           #loadmessage {
             position: fixed;
             top: 0px;
             left: 0px;
             width: 100%;
             padding: 5px 0px 5px 0px;
             text-align: center;
             font-weight: bold;
             font-size: 100%;
             color: #000000;
             background-color: #CCFF66;
             z-index: 105;
           }
        ")),
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     tags$div("Loading...",id=ns("loadmessage")))
  )

}