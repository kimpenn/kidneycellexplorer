#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

# Define UI for application that draws a histogram
shinyUI(
    tagList(
        #tags$script(src = "selectize_override.js"),
        tags$head(includeHTML("www/google_analytics.html")),
        tags$head(tags$style(HTML("body{ background: #F7F7F7; }"))),
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
        ),
        tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 1100px; /* or 950px */}")),
        fluidPage(
            useShinyjs(),
            div(
                style = 'background-color: #990000;height:50px;',
                fluidRow(
                    column(4, tags$h3("Kidney Cell Explorer", class = "topbar-element", style = "color:white;"),
                           pivot_help_UI("choose_stats_info",title = NULL, label = NULL, icn="question-circle", type = "link", tooltip = T, style = "padding-left:10px;color:white;padding-top:5px;")),
                    column(5, 
                           uiOutput("input_gene_ui"),
                        actionLink("reset_gene", NULL, icon = icon("chevron-circle-left", "fa-1x"), class = "btn_leftAlign topbar-element", style = "padding-top:5px;"),
                        shinyBS::bsTooltip(
                            "reset_gene",
                            title = "Reset gene selection",
                            options = list(container = "body")
                        )
                    ),
                    column(3, 
                           tags$p("* Multiple genes accepted for heatmap. Feature plot use the first gene.", style = "font-size:11px; padding-right:10px; color: white;/*", class = "topbar-element")
                    )
                )
            ),
            tabPanel(
                title = NULL,
                #style = 'width:1000px;',
                div(
                    style= "padding: 5px 5px 5px 5px;",
                    fluidRow(
                        column(5, 
                               wellPanel(
                                   tags$b("tSNE of"),
                                   tags$select(id="tsne_id",
                                               class = "customDrop",
                                               tags$option(value = "global", "adult mouse kidney", selected=T),
                                               tags$option(value = "nephrogenic", "nephrogenic lineage"),
                                               tags$option(value = "ureteric", "ureteric lineage")
                                   ),
                                   tags$b("color by"),
                                   conditionalPanel("input.tsne_id == 'global'",
                                                    style="display: inline-block;",
                                       tags$select(id="colorBy1",
                                                   class = "customDrop",
                                                   tags$option(value = "Lineages", "grouping", selected=T),
                                                   tags$option(value = "Ontology_ID", "ontology ID"),
                                                   tags$option(value = "orig.ident", "zonation"),
                                                   tags$option(value = "sex", "gender")
                                       )
                                   ),
                                   conditionalPanel("input.tsne_id != 'global'",
                                                    style="display: inline-block;",
                                       tags$select(id="colorBy2",
                                                   class = "customDrop",
                                                   tags$option(value = "Ontology_ID", "ontology ID"),
                                                   tags$option(value = "orig.ident", "zonation"),
                                                   tags$option(value = "sex", "gender")
                                       )
                                   ),
                                   plotOutput("plot2d", height="270px")
                               )
                               #uiOutput("plot2d_tooltip")
                        ),
                        column(4,
                               wellPanel(
                                   tags$div(
                                       uiOutput("curg_text1", inline=T),
                                       tags$select(id="palette",
                                                   class="customDrop btn_rightAlign",
                                                   tags$option(value = "RedGrey", "RedGrey", selected=T),
                                                   tags$option(value = "RdBu", "RedBlue"),
                                                   tags$option(value = "RdYlBu", "RedYellowBlue"),
                                                   tags$option(value = "viridis", "viridis"),
                                                   tags$option(value = "inferno", "inferno")
                                       )
                                   ),
                                   plotOutput("exprPlot", height="270px")
                               )
                        ),
                        column(3,
                               wellPanel(
                                   uiOutput("featurePlot_ui")
                               )
                        )
                    )
                ),
                wellPanel(
                    fluidRow(
                        column(5, uiOutput("heat_strip_title")),
                        column(3, 
                               tags$br(),
                               uiOutput("bin_stats_ui")
                        ),
                        column(2, 
                               tags$br(),
                               tags$b("Palette"),
                               tags$select(id="bin_pal",
                                           class="customDrop",
                                           tags$option(value = "RdYlBu", "RdYlBu", selected=T),
                                           tags$option(value = "RdBu", "RedBlue"),
                                           tags$option(value = "viridis", "viridis"),
                                           tags$option(value = "inferno", "inferno")
                               )
                        ),
                        column(2, 
                               tags$br(),
                               uiOutput("bin_cut_ui")
                        )
                    ),
                    uiOutput("heat_strip_ui"),
                    tags$p("* numbers in boxes correspond to schematic. 'Rescaled' value is original value rescaled (gene-wise) to 0-1 range.", style = "font-size:12px; margin-bottom:0px;")
                ),
                
                wellPanel(
                fluidRow(
                    column(5,
                           img(src ="onto_plot.png", width = "100%")
                    ),
                    column(7,
                           DT::dataTableOutput("otg_show", width = "95%"),
                           tags$br(),
                           tags$p("For detailed information on the kidney anatomy please visit GUDMAP:", style="font-style: italic; line-height: 0.5;font-size:12px;"),
                           tags$a("https://www.gudmap.org/tutorials/kidney-dev/", href = "https://www.gudmap.org/tutorials/kidney-dev/", style="font-size:12px;"),
                           tags$br(),
                           tags$a("https://www.gudmap.org/tutorials/urogenital-dev/devmus.html#DMK", href = "https://www.gudmap.org/tutorials/urogenital-dev/devmus.html#DMK", style="font-size:12px;")
                    )
                )
                ),
                wellPanel(
                fluidRow(
                    column(12,  tags$h4("Number of genes detected in metacell", class = "panel-title"))
                ),
                plotlyOutput("gene_num_plot", height="200px")
                )
            ),
            br(),
            HTML("<li>Cite Kidney Cell Explorer: Ransick, A., Lindström, N.O., Liu, J., Zhu, Q., Guo, J.J., Alvarado, G.F., Kim, A.D., Black, H.G., Kim, J. and McMahon, A.P., 2019. Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney. Developmental cell, 51(3), pp.399-413.</li>"),
            HTML("<li>Source code for Kidney Cell Explorer has been deposited to <a href=https://github.com/qinzhu/kidneycellexplorer>https://github.com/qinzhu/kidneycellexplorer</li>"),
            # Sidebar with a slider input for number of bins 
            
            br(),
            br(),
            br(),
            br(),
            br(),
            br()
        )
    )
)
