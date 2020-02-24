#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    output$input_gene_ui <- renderUI({
        input$reset_gene
        div(
            selectizeInput("input_gene", NULL, choices = gene_tbl, multiple = T, width = "310px", options = list(placeholder = 'Gene Search: e.g. Flt1', maxItems = 50)),
            class = "topbar-element btn_leftAlign"
        )
    })
    
    tsne_proj <- reactive({
        req(input$tsne_id)
        tid <- input$tsne_id
        if(tid == "global") {
            proj <- cbind(kid_proj, pmeta)
        } else if(tid == "nephrogenic") {
            proj <- cbind(neph_proj, pmeta[match(rownames(neph_proj), rownames(pmeta)),])
        } else if(tid == "ureteric") {
            proj <- cbind(uret_proj, pmeta[match(rownames(uret_proj), rownames(pmeta)),])
        } else {
            return(NULL)
        }
        return(proj)
    })
    
    cBy <- reactive({
        tid <- input$tsne_id
        if(tid == "global") {
            req(input$colorBy1)
            cBy<-input$colorBy1
        } else {
            req(input$colorBy2)
            cBy<-input$colorBy2
        }
        return(cBy)
    })
    
    pp1 <- reactive({
        req(tsne_proj(), cBy())
        tid <- input$tsne_id
        onplotAnnot <- if(cBy() == "Ontology_ID") "repel" else NULL
        plotProj(tsne_proj(), dim_col =c(1,2), group.by=cBy(), pal=NULL, size = .8, plot_title=NULL, na.col = "lightgrey", legend=T, legend.title=cBy_title[cBy()], onplotAnnot = onplotAnnot, onplotAnnotSize = 4,legend.size = 5, legend.text.size = 11, 
                 keyheight = ifelse(cBy() %in% c("res.1", "Ontology_ID"), .5, 2.2),
                 ncol = ifelse(cBy() == "Ontology_ID", 2, 1)) + 
            theme(text=element_text(family = "Helvetica"),
                  legend.margin=margin(0,0,0,0))+
            theme_void()
    })
    
    output$plot2d <- renderPlot({
        pp1()
    })
    
    # 
    # output$plot2d_tooltip <- renderUI({
    #     hover <- input$plot2d_hover
    #     #assign("hover", hover, env=.GlobalEnv)
    #     x <- nearPoints(tsne_proj, hover, maxpoints = 1)
    #     req(nrow(x) > 0)
    #     y <- as.character(x[[input$colorBy]])
    #     tip <- paste0("<b>", cBy_title[input$colorBy], ": </b>", y, "<br/>")
    #     
    #     req(length(y) > 0)
    #     style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.65); ",
    #                     "left:", hover$coords_css$x + 2, "px; top:", hover$coords_css$y + 2, "px;",
    #                     "margin:5px; padding:5px 5px 0px 5px;")
    #     
    #     # actual tooltip created as wellPanel
    #     wellPanel(
    #         style = style,
    #         p(HTML(tip))
    #     )
    # })
    
    
    
    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x <- faithful[, 2] 
        bins <- seq(min(x), max(x), length.out = 10)
        
        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
    
    output$curg_text1 <- renderUI({
        tags$b(paste0(gvals$curg[1], " expression"))
    })
    
    # output$exprPlot_tooltip <- renderUI({
    #     hover <- input$exprPlot_hover
    #     req(gvals$gene_values)
    #     #assign("hover", hover, env=.GlobalEnv)
    #     x <- nearPoints(tsne_proj, hover, maxpoints = 1)
    #     req(nrow(x) > 0)
    #     y <- round(gvals$gene_values[rownames(x),, drop=F],3)
    #     tip <- paste0(sapply(1:ncol(y), function(i) paste0("<b>", colnames(y)[i], ": </b>", y[[i]], "<br/>")), collapse = "")
    #     
    #     req(length(y) > 0)
    #     style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.65); ",
    #                     "left:", hover$coords_css$x + 2, "px; top:", hover$coords_css$y + 2, "px;",
    #                     "margin:5px; padding:5px 5px 0px 5px;")
    #     
    #     # actual tooltip created as wellPanel
    #     wellPanel(
    #         style = style,
    #         p(HTML(tip))
    #     )
    # })
    
    gvals <- reactiveValues(curg=NULL)
    observe({
        req(tsne_proj())
        if(any(!input$input_gene %in% rownames(expr_data))) {
            showModal(modalDialog(
                title = "Gene not found",
                "Expression of this gene is below threshold.",
                easyClose = TRUE,
                size = "s",
                footer = NULL
            ))
            cur_selected <- input$input_gene[input$input_gene %in% rownames(expr_data)]
            updateSelectizeInput(session, "input_gene", selected = cur_selected)
            return()
        }
        gvals$curg <-if(is.null(input$input_gene))  "Flt1" else input$input_gene
    })
    
    output$exprPlot <- renderPlot({
        req(tsne_proj())
        req(input$palette)
        cur_g <- gvals$curg[1]
        gene_values<-t(as.matrix(expr_data[cur_g,rownames(tsne_proj()), drop=F]))
        assign("gene_values", gene_values, env=.GlobalEnv)
        ecut <- max(quantile(gene_values,.99),1)
        gene_values[gene_values == 0] <- NA
        visualize_gene_expression(gene_values,  tsne_proj()[,c(1,2)],
                                  limits=c(0,ecut),
                                  marker_size = .8, ncol=1,
                                  binary =F,
                                  pal=input$palette,
                                  na.col = "#CCCCCC", legend = F) + 
            theme_void()
    })
    
    output$featurePlot_ui <- renderUI({
        tagList(
            tags$div(
                tags$select(id="bpStype",
                            class="customDrop",
                            tags$option(value = "violin", "Violin plot", selected=T),
                            tags$option(value =  "box", "Box plot"),
                            tags$option(value =  "points", "Point plot")
                ),
                uiOutput("curg_text2", inline=T)
            ),
            plotOutput("featurePlot", height="270px")
        )
    })
    output$curg_text2 <- renderUI({
        req(cBy())
        tags$b(paste0("of ", gvals$curg[1], " by ", cBy_title[cBy()]))
    })
    
    output$featurePlot <- renderPlot({
        set.seed(2019)
        req(cBy())
        cur_g <- gvals$curg[1]
        bpGroup <-cBy()
        gene_values <- t(as.matrix(expr_data[cur_g,rownames(tsne_proj()), drop=F]))
        noise <- rnorm(n = length(x = gene_values[, cur_g]))/1e+05
        gene_values[,cur_g] <- gene_values[,cur_g] + noise
        ecut <- max(quantile(gene_values,1),1)
        g1<-feature_plot(gene_values, cur_g, 
                     group.by = bpGroup, 
                     meta = tsne_proj(), 
                     style = input$bpStype, log_scale = F,
                     text.size = 12, pointSize = .5, legend = F, 
                     order.by = "none", 
                     title=NULL,
                     pal=NULL,
                     legend.title = cBy_title[cBy()],
                     ylab.label =  "Log-Normalized\nExpression"
        ) +ylim(0, ecut)
        if(bpGroup %in% c("res.1", "Ontology_ID")) {
            g1<-g1 + coord_flip()
        } else {
            g1 <- g1
        }
        return(g1)
    })
    
    
    output$heat_strip_ui <- renderUI({
        if(length(gvals$curg) == 1) {
            tagList(
                fluidRow(
                    column(10,  plotOutput("hmap_top", width = "100%", height="140px")),
                    column(2)
                ),
                fluidRow(
                    column(12,  tags$h4("Juxtamedullary nephron and ureteric epithelium cell transcriptomes", class = "panel-title"))
                ),
                fluidRow(
                    column(12, plotOutput("hmap_bottom", width = "100%", height="140px"))
                )
            )
        } else {
            plot_output_list <- lapply(1:length(gvals$curg), function(i) {
                plotname <- paste("strip", i, sep="")
                strip_ui(plotname)
            })
            tagList(
                plotOutput("heat_strip_legend", width = "100%", height="30px"),
                plotOutput("heat_strip_text", width = "100%", height="20px"),
                do.call(tagList, plot_output_list),
                tags$br()
            )
        }
    })
    
    output$heat_strip_title <- renderUI({
        if(length(gvals$curg) == 1) {
            tags$h4("Cortical nephron and ureteric epithelium cell transcriptomes", class = "panel-title")
        } else {
            div(HTML("<h4>Juxtamedullary and cortical nephron,<br/> and collecting duct cell transcriptomes</h4>"), class = "panel-title")
        }
    })
    
    output$bin_stats_ui <- renderUI({
        req(gvals$curg)
        if(length(gvals$curg) > 1) {
            sel <- tags$select(id="bin_stats",
                        class="customDrop",
                        tags$option(value = "avg", "Average expression"),
                        tags$option(value = "avg_rs", "Average expression (rescaled)", selected =T),
                        tags$option(value = "frac", "Expressed proportion")
                        #tags$option(value = "frac_rs", "Expressed proportion (rescaled)")
            )
        } else {
            sel <- tags$select(id="bin_stats",
                        class="customDrop",
                        tags$option(value = "avg", "Average expression", selected=T),
                        tags$option(value = "avg_rs", "Average expression (rescaled)"),
                        tags$option(value = "frac", "Expressed proportion")
                        #tags$option(value = "frac_rs", "Expressed proportion (rescaled)")
            )
        }
        return(sel)
    })
    
    callModule(pivot_help, "choose_stats_info", title = "Kidney Cell Explorer Tutorial:", size = "l", content = list(
        tags$ol(
            tags$li("Input single gene to visualize its expression on the kidney model, e.g., Aqp2:"),
            img(src ="single_search_Aqp2_1.png", width = "100%"),
            img(src ="single_search_Aqp2_2.png", width = "100%"),
            tags$br(),
            tags$li("Alternatively, input multiple genes to visualize the data as heatmap:"),
            img(src ="multiple_g_1.png", width = "100%"),
            img(src ="multiple_g_2.png", width = "100%"),
            tags$br(),
            tags$li("Note that gene expression level can be visualized in various ways. You can choose to visualize the average expression, the average expression scaled across metacells, or the fraction of cells in each meta cell that express that gene. By default, red in the heatmap corresponds to the maximal level of expression (or 0.1 if expression is too low), but this can also be adjusted. For different genes, the expected expression levels are generally different. Therefore we recommend users to pay attention to these adjustable parameters. "),
            img(src ="stats_option.png", width = "60%"),
            tags$br(),
            tags$hr(),
            tags$p("Please kindly cite us if you used this tool for your study: "),
            tags$em("Ransick, A., Lindström, N.O., Liu, J., Qin, Z., Guo, J.J., Alvarado, G.F., Kim, A.D., Black, H.G., Kim, J. and McMahon, A.P., 2019. Single Cell Profiling Reveals Sex, Lineage and Regional Diversity in the Mouse Kidney. bioRxiv, p.673335.")
        )

    ))
    
    output$bin_cut_ui <- renderUI({
        req(gvals$curg)
        req(input$bin_stats)
        pstats <- get_pstats(input$bin_stats, gvals$curg)
        assign("pstats", pstats, env=.GlobalEnv)
        ccut<- if(input$bin_stats == "frac") 100 else if(input$bin_stats == "avg_rs") 1 else max(quantile(as.matrix(pstats$gexpr),1), .1)
        bin_c_ui <- tagList(
            tags$b("Max:"),
            tags$input(id = "bin_cut",
                               class="customDrop",
                               style="display:inline-block; ",
                               type = "number",
                               value = round(ccut,2),
                               step = .1,
                               min = 1e-2,
                               max = 100
        ))
        bin_c_ui <- if(grepl("rs", input$bin_stats)) {
            shinyjs::hidden(bin_c_ui)
        } else bin_c_ui
        tagList(
            bin_c_ui,
            shinyjs::hidden(textInput("bin_stats_fake", NULL, value =input$bin_stats))
        )
    })

    
    observe({
        gvals$curg
        isolate({
            cur_pal <- input$bin_pal
            updateSelectInput(session, "bin_pal", selected = cur_pal)
        })
    })
    
    get_pstats <- function(stats, cur_g) {
        if(grepl("avg", stats)) {
            df <- final_expr
            tt = "Average\nexpression"
            rd = 2
        } else if(grepl("frac", stats)) {
            df <- final_frac * 100
            tt = "Cell%\nexpressed"
            rd=1
        }
        gexpr <- df[cur_g,, drop=F]
        #assign("gexpr", gexpr, env = .GlobalEnv)
        if(grepl("_rs", stats)) {
            gexpr <- as.data.frame(t(apply(gexpr, 1, function(x) {
                if(all(x==0)) return(x)
                rescale(as.numeric(x), to=c(0,1))
            })))
            colnames(gexpr) <- colnames(final_expr)
            tt <- paste0(tt, "\n(rescaled)")
        }
        return(list(gexpr = gexpr, tt = tt, rd = rd))
    }
    
    output$hmap_top <- renderPlot({
        req(length(gvals$curg) == 1)
        req(input$bin_stats_fake == input$bin_stats)
        pstats <- get_pstats(input$bin_stats, gvals$curg)
        gexpr <- pstats$gexpr
        heat_top(gexpr, pal = input$bin_pal, vlim = c(0, input$bin_cut), legend.title = pstats$tt, round.digit = pstats$rd)
    })
    
    output$hmap_bottom <- renderPlot({
        req(length(gvals$curg) == 1)
        req(input$bin_stats_fake == input$bin_stats)
        pstats <- get_pstats(input$bin_stats, gvals$curg)
        gexpr <- pstats$gexpr
        heat_bottom(gexpr, pal = input$bin_pal, vlim = c(0, input$bin_cut), legend.title = pstats$tt, round.digit = pstats$rd)
    })
    
    output$heat_strip_legend <- renderPlot({
        #req(grepl("_rs", input$bin_stats))
        par(mar=c(0,0,0,0))
        plot(1:10,seq(0, 1,length.out=10), type='n', main="", xlab="x", ylab="y", xlim = c(0,10), ylim = c(0,1), bty="n",xaxt="n", yaxt="n", ann = F)
        rasterImage(as.raster(t(matrix(get_numeric_color("RdYlBu"), ncol=1))),10-1, 0.5, 10, 1)
        text(x = seq(9, 10,l=5), y= .55, pos=1, labels = round(seq(0,input$bin_cut,l=5),2), cex = .8)
        legend.title <- switch(input$bin_stats,
                      "avg"="Average expression",
                      "avg_rs"= "Average expression (rescaled)",
                      "frac"= "Expressed proportion",
                      "frac_rs"= "Expressed proportion (rescaled)")
        text(x = 9, y = .7,legend.title, cex = 1, pos = 2, font=2)
    })
    
    output$heat_strip_text <- renderPlot({
        #req(grepl("_rs", input$bin_stats))
        par(mar=c(0,0,0,0))
        max_x <- 10
        max_y <- 1
        min_x <- 0
        min_y <- 0
        
        # First assign correct position to the numbers
        strip_order <- c('1', '2', '3', '4', '5', '6', '7', '8', '9A', '9B', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19','22', '25','27','21','24', '20',  '23',  '26',  '28', '29', '30', '31', '32')
        bin_names <- strip_order
        bin_tbl <- data.frame(matrix(ncol = length(strip_order), nrow = 4), row.names = c("x_1","y_1","x_2","y_2"))
        colnames(bin_tbl) <- strip_order
        y_adj_min <- .01
        y_adj_max <- -.05
        ybinh <- max_y+y_adj_max-min_y-y_adj_min
        
        x_adj_min <- 0.25
        x_adj_max <- 0
        x_out <- length(strip_order)
        
        xbinw <- (max_x + x_adj_max- min_x-x_adj_min + 0.5)/x_out
        bin_tbl["x_1",] <- seq(min_x+x_adj_min, max_x+ x_adj_max, length.out = x_out)
        bin_tbl["x_2",] <- seq(min_x+x_adj_min, max_x+ x_adj_max, length.out = x_out)+ xbinw
        
        bin_tbl["y_1",] <- min_y
        bin_tbl["y_2",] <- min_y + ybinh
        plot(1:10,seq(0, max_y,length.out=10), type='n', main="", xlab="x", ylab="y", xlim = c(0,10), ylim = c(0,max_y), bty="n",xaxt="n", yaxt="n")
        for(i in 1:ncol(bin_tbl)){
            cur_nm <- colnames(bin_tbl)[i]
            if(cur_nm %in% c('3', '5', '7')) {
                col = "#df65b0"
            } else if (cur_nm %in% c('4', '6', '8')) {
                col = "#2b8cbe"
            } else {
                col = "black"
            }
            text(x = bin_tbl["x_1",cur_nm] + xbinw/2, y = bin_tbl["y_1",cur_nm] + ybinh/2, cur_nm, font=2, col = col)

        }
    })
    
    
    observe({
        req(input$bin_stats_fake == input$bin_stats)
        #req(grepl("_rs", input$bin_stats))
        input$bin_stats
        input$bin_pal
        input$bin_cut
        lapply(1:length(gvals$curg), function(i) {
            plotname <- paste("strip", i, sep="")
            callModule(strip_server, id = plotname, cur_g = gvals$curg[i], pstats = get_pstats(input$bin_stats, gvals$curg[i]), pal = input$bin_pal, vlim = c(0, input$bin_cut))
        })
    })

    
    
    output$otg_show <- DT::renderDataTable({
        DT::datatable(otg_tbl, rownames = F, selection = 'none',
                      style = 'bootstrap', class = 'table-condensed',
                      options = list(
                          searching=F, 
                          scrollX = T,
                          scrollY="450px",
                          paging = F,
                          lengthChange=F,
                          ordering=F,
                          info = F
                      )) %>%
            DT::formatStyle(columns = c(1, 2), fontSize = '80%')
    })
    
    output$gene_num_plot <- plotly::renderPlotly({
        gcount<-data.frame(id = colnames(final_expr), count=colSums(final_expr>0))
        gcount$id <- factor(gcount$id, levels = gcount$id)
        hues = seq(15, 375, length = length(gcount$id) + 1)
        gcolors<-hcl(h = hues, l = 65, c = 100)[1:length(gcount$id)]
        plot_ly(
            data=gcount,
            x = ~id,
            y = ~count,
            color = ~id, 
            colors = gcolors, 
            type = "bar"
        ) %>%
            layout(yaxis = list(title = '#Genes in metacell'),
                   xaxis = list(title = 'Ontology ID'),
                   showlegend = F)
    })

})
