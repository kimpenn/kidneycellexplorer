





#' @export
strip_ui <- function(id) {
    ns <- NS(id)
    fluidRow(
        column(12, plotOutput(ns("hmap_strip"), width = "100%", height="30px"))
    )
}

#' @export
strip_server <- function(input, output, session, cur_g, pstats, pal = "RdYlBu", vlim) {
    output$hmap_strip <- renderPlot({
        heat_strip(pstats$gexpr, pal = pal, vlim = vlim)
    })
}



