



plotProj <- function (proj, dim_col = c(1,2), group.by=NULL, pal=NULL, size = 1, plot_title=NULL, na.col = "lightgrey", alpha=NULL, alpha_level=0.1, legend=T, trans = "identity", onplotAnnot=NULL, onplotAnnotSize = 2,  legend.size = 4, legend.text.size = 3, legend.position = "top", legend.title = waiver(), keywidth=0.1, keyheight=0.1, ncol = NULL, nudge_x = 0, nudge_y = 0, limits = NULL, breaks = waiver(), strwrap=T, ...) {
    plot_col <- colnames(proj)[dim_col]
    if(!is.null(alpha)) {
        proj$alpha <- alpha
    } else {
        proj$alpha <- rep("f", length(nrow(proj)))
    }
    alpha_manual <- c("f"=1,"t"=alpha_level)
    if(!is.null(limits)) {
        proj[[group.by]][proj[[group.by]] < limits[1]] <- limits[1]
        proj[[group.by]][proj[[group.by]] > limits[2]] <- limits[2]
    }
    pp<-ggplot(proj, aes_string(plot_col[1],plot_col[2])) +
        geom_point(aes_string(color=group.by, alpha="alpha"), size=size, stroke = 0) +
        scale_alpha_manual(values=alpha_manual) +
        theme_bw() +
        ggtitle(plot_title) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = legend.position)
    if(!is.null(onplotAnnot)) {
        label_data <- proj %>% group_by_at(group.by) %>% summarize_at(plot_col, median)
        if(length(breaks) > 0) {
            label_data <- label_data[label_data[[group.by]] %in% breaks,,drop=F]
        }
        if(onplotAnnot == "text") {
            geom_func <- geom_text
        } else if(onplotAnnot == "label") {
            geom_func <- geom_label
        } else {
            geom_func <- geom_label_repel
        }
        pp<- pp + geom_func(
            aes_string(
                label = group.by,
                color = group.by
            ),
            fontface = "bold",
            size = onplotAnnotSize,
            nudge_x = nudge_x,
            nudge_y = nudge_y,
            show.legend = F,
            data = label_data
        )
    }
    if(legend) {
        pp <- pp + guides(alpha=F, color = guide_legend(override.aes = list(size=legend.size),
                                                        title.theme = element_text(size = legend.text.size*1.2),
                                                        label.theme = element_text(size = legend.text.size, lineheight = 1),
                                                        title = legend.title, 
                                                        keywidth=keywidth, 
                                                        keyheight=keyheight, ncol = ncol))
    } else {
        pp <- pp + guides(alpha=F, color = F)
    }
    if(!is.null(names(pal))) {
        pp<-pp + scale_color_manual(values = pal, na.value=na.col, breaks = breaks)
    } else {
        if(!is.null(pal)) {
            if(is.factor(proj[[group.by]]) || is.character(proj[[group.by]])) {
                pp<-pp + scale_color_manual(values = get_factor_color(unique(na.omit(proj[[group.by]])), pal=pal, ...), na.value=na.col, breaks = breaks)
            } else {
                pp<-pp + scale_colour_gradientn(colours=get_numeric_color(pal), trans=trans)
                if(legend) {
                    pp<-pp + guides(color = guide_colorbar(barwidth = 10, barheight = 1))
                }
            }
        }
    }
    return(pp)
}



factor_color_opt <- function() {
    allowed_pals <- c("Set1", "Set2", "Paired", "Dark2", "Accent")
    return(allowed_pals)
}

get_factor_color <-function (labels, pal = "Set1", maxCol = 9, nogrey = T)
{
    unq <- unique(labels)
    maxCol <- min(length(unq), maxCol)
    hmcol <- RColorBrewer::brewer.pal(maxCol, pal)
    if(nogrey) {
        hmcol <- hmcol[!hmcol %in% c("#999999","#B3B3B3")]
    }
    colv <- rep(NA, length(labels))
    #if (length(unq) > maxCol) {
    cp <- colorRampPalette(hmcol)
    hmcol <- cp(length(unq))
    #}
    for (i in 1:length(unq)) {
        colv[labels == unq[i]] <- hmcol[i]
    }
    return(colv)
}


visualize_gene_expression <- function (gene_values, projection, limits = c(0, 10), marker_size = 0.1,
                                       title = NULL, ncol = NULL, title_size = 10, alpha=NULL, alpha_manual=NULL, binary = F, pal="rainbow2", na.col = "lightgrey", legend = T, legend_name = "Expression")
{
    projection_names <- colnames(projection)
    #colnames(projection) <- c("Component.1", "Component.2")
    if(!binary) {
        gene_values[gene_values < limits[1]] <- limits[1]
        gene_values[gene_values > limits[2]] <- limits[2]
        proj_gene <- data.frame(cbind(projection[c(1,2)], gene_values))
        proj_gene_melt <- melt(proj_gene, id.vars = colnames(projection))
        #idx_region <- which(proj_gene_melt$value > 0)
        use_color <- get_numeric_color(pal)
        #assign("proj_gene_melt", proj_gene_melt, env=.GlobalEnv)
        p <- ggplot(proj_gene_melt, aes_string(x=colnames(projection)[1], y=colnames(projection)[2])) +
            #geom_point(size=marker_size,color=na.col,show.legend=FALSE, stroke=0) +
            geom_point(aes_string(x=colnames(projection)[1], y=colnames(projection)[2], colour = "value"), size = marker_size, stroke=0)+
            scale_alpha_manual(values=alpha_manual) +
            guides(alpha=F)
        
        p <- p +
            #facet_wrap(~variable, ncol=ncol) +
            scale_color_gradientn(colours = use_color, na.value = na.col) +
            labs(x = projection_names[1], y = projection_names[2])
        if (!is.null(title)) {
            p <- p + ggtitle(title)
        }
        if(legend) {
            p <- p + 
                guides(colour = guide_colorbar(title = legend_name))
        } else {
            p <- p + guides(colour = F)
        }
        p <- p + theme_bw() + 
            theme(legend.position = c("top"), 
                  plot.title = element_text(hjust = 0.5), 
                  strip.text = element_text(size=title_size), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
    } 
    return(p)
}

numeric_color_opt <- function() {
    allowed_pals <- c('blue_green_gold', 'black_red_gold', 'black_red', 'red_yellow', 'black_yellow', 'black_yellow_gold' , 'rainbow2', 'rainbow', 'gg_color_hue', 'RdOgYl', 'RdYlBu', 'RdBu', 'viridis', 'magma', 'plasma', 'inferno', 'RedGrey')
    return(allowed_pals)
}

#' @export
get_numeric_color <- function(palette = NULL) {
    if(is.null(palette)) stop("please specify a palette")
    allowed_pals <- numeric_color_opt()
    if(!palette %in% allowed_pals) stop(paste0("Please specify one of '", paste(allowed_pals, collapse = "', '"), "'."))
    if(palette %in% get_brewer_set(c("sequential", "diverging"))) {
        colorRampPalette(rev(RColorBrewer::brewer.pal(9,palette)))(100)
    } else if(palette %in% list("viridis" = "viridis", "magma" = "magma", "plasma" = "plasma", "inferno" = "inferno")) {
        viridis(n = 100, option = palette)
    } else if(palette == "diverge_hcl") {
        colorRampPalette(colorspace::diverge_hcl(7))(100)
    } else if(palette == "redgreen") {
        rev(gplots::redgreen(75))
    } else if(palette == "rainbow") {
        colorRampPalette(rev(rainbow(10)))(100)
    } else if(palette == "rainbow2") {
        c("#CCCCCCCC",rainbow(500)[50:500])
    } else if(palette == "RedGrey") {
        c("#CCCCCC", "#CC0000")
    } else if(palette == "RdOgYl") {
        c("grey85", "red", "orange", "yellow")
    } else if(palette == "gg_color_hue") {
        gg_color_hue2(10)
    } else if(palette == "blue_green_gold"){
        c("grey85", "blue", "green", "#FFD200", "gold")
    } else if(palette == "black_red_gold"){
        c("grey85", "black", "red", "#FFD200")
    } else if(palette == "black_red") {
        c("grey85", "black", "red")
    } else if(palette == "red_yellow") {
        c("grey85",  "red", "yellow")
    } else if(palette == "black_yellow") {
        c("grey85",  "black", "yellow")
    } else if(palette == "black_yellow_gold") {
        c("grey85",  "black", "yellow", "gold")
    }
}


get_brewer_set <- function(palette = c("sequential", "diverging", "qualitative")) {
    match.arg(palette,
              several.ok = TRUE)
    
    sequential_palette <- c('Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys',
                            'Oranges', 'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
                            'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd')
    names(sequential_palette) <- sequential_palette
    diverging_palette <- c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
    names(diverging_palette) <- diverging_palette
    qualitative_palette <- c('Accent','Dark2','Paired', 'Pastel1', 'Pastel2','Set1', 'Set2', 'Set3')
    names(qualitative_palette) <- qualitative_palette
    return_palette = list()
    if("qualitative" %in% palette) {
        return_palette <- c(return_palette, as.list(qualitative_palette))
    }
    if("diverging" %in% palette) {
        return_palette <- c(return_palette, as.list(diverging_palette))
    }
    if("sequential" %in% palette) {
        return_palette <- c(return_palette, as.list(sequential_palette))
    }
    
    return(return_palette)
}



feature_plot <- function(df, selected_gene, group.by = "sample", meta = NULL, pal = "Set1", style = "box", log_scale = F, legend = T, legend.title = waiver(), legend.pos = "top", legend.point.size=4, text.size = 15, pointSize = 3, na.col = "lightgrey", breaks = waiver(),  axis.text.angle = 90, order.by = "none", ylab.label = "expression", title=NULL){
    if(is.null(df) || nrow(df) == 0 || is.null(meta) || is.null(palette)) {
        return()
    }
    colnames(df) <- "expression_level"
    df <- cbind(df, meta)
    
    if(order.by == "mean") {
        group_mean <- df %>% dplyr::group_by_at(group.by) %>% dplyr::summarize(mean = mean(expression_level))
        group_order <- group_mean[[group.by]][order(group_mean$mean, decreasing = T)]
        df[[group.by]] <- factor(as.character(df[[group.by]]), levels = group_order, ordered=T)
    }
    
    g1 <- ggplot(df, aes_string(x=group.by, y="expression_level"))
    
    if(style == "box") {
        g1 <- g1 + geom_boxplot(aes_string(fill = group.by, alpha = 0.2), outlier.size = 1, outlier.stroke = 0, outlier.alpha = .2)
    } else if(style == "violin") {
        g1 <- g1 + geom_violin(aes_string(fill = group.by, alpha = 0.2), trim = T, scale = "width")
    } else if(style == "points") {
        g1 <- g1 + geom_jitter(size = pointSize, aes_string(colour = group.by))
    }
    
    if(!is.null(names(pal))) {
        g1<-g1 + 
            scale_color_manual(values = pal, na.value=na.col, breaks = breaks) +
            scale_fill_manual(values = pal, na.value=na.col, breaks = breaks)
    } else {
        if(!is.null(pal)) {
            if(is.factor(df[[group.by]]) || is.character(df[[group.by]])) {
                pals <- get_factor_color(unique(na.omit(df[[group.by]])), pal=pal)
                g1<-g1 + 
                    scale_color_manual(values = pals, na.value=na.col, breaks = breaks) + 
                    scale_fill_manual(values = pals, na.value=na.col, breaks = breaks)
            } 
        }
    } 
    
    g1 <- g1 + 
        ggtitle(title) +
        xlab(legend.title) +
        ylab(ylab.label) + 
        theme(text = element_text(size=text.size), axis.text.x = element_text(angle=axis.text.angle, hjust=1), legend.position=legend.pos, plot.title = element_text(hjust = 0.5)) 
    if(legend) {
        g1 <- g1 + guides(alpha = F, fill=F, color = guide_legend(override.aes = list(size=legend.point.size), title = legend.title))
    } else {
        g1 <- g1 + guides(alpha = F, fill=F, color= F)
    }
    if(log_scale) {
        g1 <- g1 + scale_y_log10(breaks=c(25,100,400))
    }
    
    return(
        g1 + monocle_theme_opts()
    )
}


monocle_theme_opts <- function(){
    theme(strip.background = element_rect(colour = "white", fill = "white")) + 
        theme(panel.border = element_blank()) + theme(axis.line.x = element_line(size = 0.25, 
                                                                                 color = "black")) + theme(axis.line.y = element_line(size = 0.25, 
                                                                                                                                      color = "black")) + theme(panel.grid.minor.x = element_blank(), 
                                                                                                                                                                panel.grid.minor.y = element_blank()) + theme(panel.grid.major.x = element_blank(), 
                                                                                                                                                                                                              panel.grid.major.y = element_blank()) + theme(panel.background = element_rect(fill = "white")) + 
        theme(legend.key = element_blank())
}






# 
# SingleVlnPlot <- function(feature, data, cell.ident, do.sort, y.max, size.x.use, 
#                            size.y.use, size.title.use, adjust.use, point.size.use, cols.use, 
#                            gene.names, y.log, x.lab.rot, y.lab.rot, legend.position, 
#                            remove.legend) 
# {
#     feature.name <- colnames(data)
#     colnames(data) <- "feature"
#     feature <- "feature"
#     set.seed(seed = 42)
#     data$ident <- cell.ident
#     if (do.sort) {
#         data$ident <- factor(x = data$ident, levels = names(x = rev(x = sort(x = tapply(X = data[, 
#                                                                                                  feature], INDEX = data$ident, FUN = mean)))))
#     }
#     if (y.log) {
#         noise <- rnorm(n = length(x = data[, feature]))/200
#         data[, feature] <- data[, feature] + 1
#     }
#     else {
#         noise <- rnorm(n = length(x = data[, feature]))/1e+05
#     }
#     if (all(data[, feature] == data[, feature][1])) {
#         warning(paste0("All cells have the same value of ", feature, 
#                        "."))
#     }
#     else {
#         data[, feature] <- data[, feature] + noise
#     }
#     y.max <- SetIfNull(x = y.max, default = max(data[, feature]))
#     plot <- ggplot(data = data, mapping = aes(x = factor(x = ident), 
#                                               y = feature)) + geom_violin(scale = "width", adjust = adjust.use, 
#                                                                           trim = TRUE, mapping = aes(fill = factor(x = ident))) + 
#         guides(fill = guide_legend(title = NULL)) + xlab("Identity") + 
#         NoGrid() + ggtitle(feature) + theme(plot.title = element_text(size = size.title.use, 
#                                                                       face = "bold"), legend.position = legend.position, axis.title.x = element_text(face = "bold", 
#                                                                                                                                                      colour = "#990000", size = size.x.use), axis.title.y = element_text(face = "bold", 
#                                                                                                                                                                                                                          colour = "#990000", size = size.y.use))
#     if (point.size.use != 0) {
#         plot <- plot + geom_jitter(height = 0, size = point.size.use)
#     }
#     plot <- plot + ggtitle(feature.name)
#     if (y.log) {
#         plot <- plot + scale_y_log10()
#     }
#     else {
#         plot <- plot + ylim(min(data[, feature]), y.max)
#     }
#     if (feature %in% gene.names) {
#         if (y.log) {
#             plot <- plot + ylab(label = "Log Expression level")
#         }
#         else {
#             plot <- plot + ylab(label = "Expression level")
#         }
#     }
#     else {
#         plot <- plot + ylab(label = "")
#     }
#     if (!is.null(x = cols.use)) {
#         plot <- plot + scale_fill_manual(values = cols.use)
#     }
#     if (x.lab.rot) {
#         plot <- plot + theme(axis.text.x = element_text(angle = 45, 
#                                                         hjust = 1, size = size.x.use))
#     }
#     if (y.lab.rot) {
#         plot <- plot + theme(axis.text.x = element_text(angle = 90, 
#                                                         size = size.y.use))
#     }
#     if (remove.legend) {
#         plot <- plot + theme(legend.position = "none")
#     }
#     return(plot)
# }
# 
