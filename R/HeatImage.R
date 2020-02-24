

### Top heatmap for cortical ###

heat_top <- function(gexpr = NULL, vlim = c(0, max(quantile(gexpr,1), .1)), pal = "RdYlBu", legend.title=NULL, round.digit=2) {
    #hmap_img=readJPEG(empty_frame_path)
    max_x <- 10
    max_y <- 1
    min_x <- 0
    min_y <- 0
    
    # First assign correct position to the numbers
    bin_names <- colnames(gexpr)
    
    bin_list <- list(
        bin_set1 <- c("1","2"),
        bin_set2 <- c("9A", "15", "16", "17", "18", "19"),
        bin_set3 <- c("20", "23", "26", "28","29","30", "31", "32"),
        bin_set4 <- c("21", "24"),
        bin_set5 <- c("22", "25", "27"),
        bin_set6 <- c("3", "5", "7"),
        bin_set7 <- c("4","6", "8")
    )
    
    bin_tbl <- data.frame(matrix(ncol =  length(unlist(bin_list)), nrow = 4), row.names = c("x_1","y_1","x_2","y_2"))
    colnames(bin_tbl) <-unlist(bin_list)
    y_adj_min <- .05
    y_adj_max <- -.1
    y_int <- 0.01
    ybinh <- (max_y+y_adj_max-min_y-y_adj_min-y_int*2)/3
    
    x_adj_min <- .1
    x_adj_max <- -1.12
    x_out <- length(bin_set1) + length(bin_set2) + length(bin_set3) + length(bin_set6)
    
    xbinw <- (max_x + x_adj_max- min_x-x_adj_min + 0.5)/x_out
    x_1_series <- seq(min_x+x_adj_min, max_x+ x_adj_max, length.out = x_out)
    x_2_series <- seq(min_x+x_adj_min, max_x+ x_adj_max, length.out = x_out)+ xbinw
    
    
    for(i in unlist(bin_list[c(1,2)])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 1.5*ybinh
        bin_tbl["y_2",i] <- min_y + y_adj_min + 2.5*ybinh
    }
    for(i in unlist(bin_list[3])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 2*ybinh
        bin_tbl["y_2",i] <- min_y + y_adj_min + 3*ybinh
    }
    for(i in unlist(bin_list[4])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 1*ybinh
        bin_tbl["y_2",i] <- min_y + y_adj_min + 2*ybinh
    }
    for(i in unlist(bin_list[5])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min 
        bin_tbl["y_2",i] <- min_y + y_adj_min + 1*ybinh
    }
    
    for(i in unlist(bin_list[6])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 2.5*ybinh
    }
    for(i in unlist(bin_list[7])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 1.5*ybinh
    }
    
    bin_tbl["x_1",bin_list[[1]]] <- x_1_series[1:2]
    bin_tbl["x_2",bin_list[[1]]] <- x_2_series[1:2]
    bin_tbl["x_1",bin_list[[2]]] <- x_1_series[6:11]
    bin_tbl["x_2",bin_list[[2]]] <- x_2_series[6:11]
    bin_tbl["x_1",bin_list[[3]]] <- x_1_series[12:19]
    bin_tbl["x_2",bin_list[[3]]] <- x_2_series[12:19]
    bin_tbl["x_1",bin_list[[4]]] <- x_1_series[12:13]
    bin_tbl["x_2",bin_list[[4]]] <- x_2_series[12:13]
    bin_tbl["x_1",bin_list[[5]]] <- x_1_series[12:14] + xbinw/2
    bin_tbl["x_2",bin_list[[5]]] <- x_2_series[12:14] + xbinw/2
    bin_tbl["x_1",bin_list[[6]]] <- x_2_series[3:5] - xbinw/2
    bin_tbl["x_1",bin_list[[7]]] <- x_2_series[3:5] - xbinw/2
    
    bin_tbl[c("x_1", "x_2"),bin_list[[3]][-1]] <- bin_tbl[c("x_1", "x_2"),bin_list[[3]][-1]] + .05
    bin_tbl[c("x_1", "x_2"),bin_list[[4]][-1]] <- bin_tbl[c("x_1", "x_2"),bin_list[[4]][-1]] + .05
    bin_tbl[c("x_1", "x_2"),bin_list[[5]][-1]] <- bin_tbl[c("x_1", "x_2"),bin_list[[5]][-1]] + .05
    
    
    cor_len = 10
    glim <- c(vlim[1], vlim[2])
    if(pal %in% c("RdYlBu", "RdBu")) {
        orig_col<-rev(RColorBrewer::brewer.pal(9,pal))
        pals <- colorRampPalette(orig_col)(cor_len)
    } else if(pal == "RedGrey"){
        orig_col<- c("grey", "red")
        pals <- colorRampPalette(orig_col)(cor_len)
    } else if(pal %in% c("viridis", "inferno")) {
        pals <- viridis(n = cor_len, option =pal)
    }
    gexpr_col <- pals[sapply(gexpr, function(g) {
        thresh_pass<-which(g >= seq(glim[1],glim[2], length.out = cor_len))
        if(!length(thresh_pass))thresh_pass<-1
        max(thresh_pass)})]
    names(gexpr_col) <- colnames(gexpr)
    
    par(mar=c(0,0,0,0))

    plot(1:10,seq(0, max_y,length.out=10), type='n', main="", xlab="x", ylab="y", xlim = c(0,10), ylim = c(0,max_y), bty="n",xaxt="n", yaxt="n")
    #rasterImage(hmap_img, min_x, min_y, max_x, max_y, interpolate = F)
    #plot.new()
    for(i in 1:ncol(bin_tbl)){
        cur_nm <- colnames(bin_tbl)[i]
        #if(cur_nm %in% c(22)) lty=6 else lty=1
        lty = 1
        if(cur_nm %in% c(unlist(bin_list[c(1,2,6,7)]),20,21,22))  bcolor = "#7C7C7C" else bcolor <- "#3A3ACE"
        if(cur_nm %in% unlist(bin_list[c(6,7)])) {
            plotrix::draw.ellipse(x=bin_tbl["x_1",cur_nm], y=bin_tbl["y_1",cur_nm], a = xbinw/2, b = ybinh/2, angle = 0, border = "#7C7C7C", col = gexpr_col[cur_nm], lwd = 3)
            text(x = bin_tbl["x_1",cur_nm], y = bin_tbl["y_1",cur_nm], cur_nm, font=2)
        } else {
            rect(bin_tbl["x_1",cur_nm],bin_tbl["y_1",cur_nm],bin_tbl["x_2",cur_nm],bin_tbl["y_2",cur_nm], col=gexpr_col[cur_nm],border=bcolor, lwd=3, lty=lty)
            text(x = bin_tbl["x_1",cur_nm] + xbinw/2, y = bin_tbl["y_1",cur_nm] + ybinh/2, cur_nm, font=2)
        }
    }
    #text(min_x+.2, bin_tbl["y_1","1"]+1*ybinh, "Gene\nexpression", cex = 1.1, pos = 1)
    text(bin_tbl["x_1","5"], bin_tbl["y_1","5"]+1.25*ybinh, "female", cex = 1, pos = 1)
    text(bin_tbl["x_1","6"], bin_tbl["y_1","6"]-.57*ybinh, "male", cex = 1, pos = 1)
    
    bx0 = bin_tbl["x_1","3"]-.2
    bx1 =  bin_tbl["x_1","7"] +.2
    by0 = bin_tbl["y_1","3"]+.45*ybinh
    by1 = bin_tbl["y_1","3"]+.45*ybinh
    baj <- .05
    lines(x = c(bx0, bx0, bx1, bx1), y = c(by0, by0+baj, by0+baj, by0))
    lines(x = c(bx0, bx0, bx1, bx1), y = c(by0-1.9*ybinh, by0-1.9*ybinh-baj, by0-1.9*ybinh-baj, by0-1.9*ybinh))
    
    legend_image <- as.raster(matrix(rev(get_numeric_color(pal)), ncol=1))
    min_pos <- .13
    max_pos <- .73
    text(x= max_x-.1, y = seq(min_pos, max_pos,l=5), labels = round(seq(vlim[1],vlim[2],l=5),round.digit), cex = .8)
    rasterImage(legend_image,max_x-.4, min_pos, max_x-.3, max_pos)
    text(max_x-.13, bin_tbl["y_1","5"]+1.1*ybinh,legend.title, cex = .8, pos = 1, font=2)
}





### Bottom heatmap for jux ###

heat_bottom <- function(gexpr = NULL, vlim = c(0, max(quantile(gexpr,1), .1)), pal = "RdYlBu", legend.title=NULL, round.digit=2) {
    #hmap_img=readJPEG(empty_frame_path)
    max_x <- 10
    max_y <- 1
    min_x <- 0
    min_y <- 0
    
    # First assign correct position to the numbers
    bin_names <- colnames(gexpr)
    
    bin_list <- list(
        bin_set1 <- c("1","2"),
        bin_set2 <- c("14", "15", "16", "17", "18", "19"),
        bin_set3 <- c("20", "23", "26", "28", "29","30", "31", "32"),
        bin_set4 <- c("21", "24"),
        bin_set5 <- c("22", "25", "27"),
        bin_set6 <- c("3", "5", "7"),
        bin_set7 <- c("4","6", "8"),
        bin_set8 <- c("9B", "10", "11"),
        bin_set9 <- c("12"),
        bin_set10 <- c("13")
    )
    
    bin_tbl <- data.frame(matrix(ncol =  length(unlist(bin_list)), nrow = 4), row.names = c("x_1","y_1","x_2","y_2"))
    colnames(bin_tbl) <-unlist(bin_list)
    y_adj_min <- .05
    y_adj_max <- -.1
    y_int <- 0.01
    ybinh <- (max_y+y_adj_max-min_y-y_adj_min-y_int*2)/3
    
    x_adj_min <- .02
    x_adj_max <- -1
    x_out <- length(bin_set1) + length(bin_set2) + length(bin_set3) + length(bin_set6) + length(bin_set8) + length(bin_set9) 
    
    xbinw <- (max_x + x_adj_max- min_x-x_adj_min + 0.5)/x_out
    x_1_series <- seq(min_x+x_adj_min, max_x+ x_adj_max, length.out = x_out)
    x_2_series <- seq(min_x+x_adj_min, max_x+ x_adj_max, length.out = x_out)+ xbinw
    
    
    for(i in unlist(bin_list[c(1,2,8)])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 1.5*ybinh
        bin_tbl["y_2",i] <- min_y + y_adj_min + 2.5*ybinh
    }
    for(i in unlist(bin_list[c(3,9)])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 2*ybinh
        bin_tbl["y_2",i] <- min_y + y_adj_min + 3*ybinh
    }
    for(i in unlist(bin_list[c(4,10)])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 1*ybinh
        bin_tbl["y_2",i] <- min_y + y_adj_min + 2*ybinh
    }
    for(i in unlist(bin_list[5])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min 
        bin_tbl["y_2",i] <- min_y + y_adj_min + 1*ybinh
    }
    
    for(i in unlist(bin_list[6])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 2.5*ybinh
    }
    for(i in unlist(bin_list[7])) {
        bin_tbl["y_1",i] <- min_y + y_adj_min + 1.5*ybinh
    }
    
    bin_tbl["x_1",bin_list[[1]]] <- x_1_series[1:2]
    bin_tbl["x_2",bin_list[[1]]] <- x_2_series[1:2]
    bin_tbl["x_1",bin_list[[6]]] <- x_2_series[3:5] - xbinw/2
    bin_tbl["x_1",bin_list[[7]]] <- x_2_series[3:5] - xbinw/2
    bin_tbl["x_1",bin_list[[8]]] <- x_1_series[6:8]
    bin_tbl["x_2",bin_list[[8]]] <- x_2_series[6:8]
    bin_tbl["x_1",unlist(bin_list[c(9,10)])] <- x_1_series[9]
    bin_tbl["x_2",unlist(bin_list[c(9,10)])] <- x_2_series[9]
    bin_tbl["x_1",bin_list[[2]]] <- x_1_series[10:15]
    bin_tbl["x_2",bin_list[[2]]] <- x_2_series[10:15]
    bin_tbl["x_1",bin_list[[3]]] <- x_1_series[16:23]
    bin_tbl["x_2",bin_list[[3]]] <- x_2_series[16:23]
    bin_tbl["x_1",bin_list[[4]]] <- x_1_series[16:17]
    bin_tbl["x_2",bin_list[[4]]] <- x_2_series[16:17]
    bin_tbl["x_1",bin_list[[5]]] <- x_1_series[16:18] + xbinw/2
    bin_tbl["x_2",bin_list[[5]]] <- x_2_series[16:18] + xbinw/2
    
    bin_tbl[c("x_1", "x_2"),bin_list[[3]][-1]] <- bin_tbl[c("x_1", "x_2"),bin_list[[3]][-1]] + .05
    bin_tbl[c("x_1", "x_2"),bin_list[[4]][-1]] <- bin_tbl[c("x_1", "x_2"),bin_list[[4]][-1]] + .05
    bin_tbl[c("x_1", "x_2"),bin_list[[5]][-1]] <- bin_tbl[c("x_1", "x_2"),bin_list[[5]][-1]] + .05
    
    cor_len = 10
    glim <- c(vlim[1], vlim[2])
    if(pal %in% c("RdYlBu", "RdBu")) {
        orig_col<-rev(RColorBrewer::brewer.pal(9,pal))
        pals <- colorRampPalette(orig_col)(cor_len)
    } else if(pal == "RedGrey"){
        orig_col<- c("grey", "red")
        pals <- colorRampPalette(orig_col)(cor_len)
    } else if(pal %in% c("viridis", "inferno")) {
        pals <- viridis(n = cor_len, option =pal)
    }
    gexpr_col <- pals[sapply(gexpr, function(g) {
        thresh_pass<-which(g >= seq(glim[1],glim[2], length.out = cor_len))
        if(!length(thresh_pass))thresh_pass<-1
        max(thresh_pass)})]
    names(gexpr_col) <- colnames(gexpr)
    
    
    par(mar=c(0,0,0,0))
    
    plot(1:10,seq(0, max_y,length.out=10), type='n', main="", xlab="x", ylab="y", xlim = c(0,10), ylim = c(0,max_y), bty="n",xaxt="n", yaxt="n")
    #rasterImage(hmap_img, min_x, min_y, max_x, max_y, interpolate = F)
    #plot.new()
    #assign("bin_tbl",bin_tbl, env=.GlobalEnv)
    for(i in 1:ncol(bin_tbl)){
        cur_nm <- colnames(bin_tbl)[i]
        #if(cur_nm %in% c(20,21,24)) lty=2 else lty=1
        lty = 1
        if(cur_nm %in% c(unlist(bin_list[c(1,2,6,7,8,9,10)]), 20,21,22)) bcolor = "#7C7C7C" else bcolor <- "#3A3ACE"
        if(cur_nm %in% unlist(bin_list[c(6,7)])) {
            plotrix::draw.ellipse(x=bin_tbl["x_1",cur_nm], y=bin_tbl["y_1",cur_nm], a = xbinw/2, b = ybinh/2, angle = 0, border = "#7C7C7C", col = gexpr_col[cur_nm], lwd = 3)
            text(x = bin_tbl["x_1",cur_nm], y = bin_tbl["y_1",cur_nm], cur_nm, font=2)
        } else {
            rect(bin_tbl["x_1",cur_nm],bin_tbl["y_1",cur_nm],bin_tbl["x_2",cur_nm],bin_tbl["y_2",cur_nm], col=gexpr_col[cur_nm],border=bcolor, lwd=3, lty=lty)
            text(x = bin_tbl["x_1",cur_nm] + xbinw/2, y = bin_tbl["y_1",cur_nm] + ybinh/2, cur_nm, font=2)
        }
    }
    # text(min_x+.1, bin_tbl["y_1","1"]+1*ybinh, "Gene\nexpression", cex = 1.1, pos = 1)
    text(bin_tbl["x_1","5"], bin_tbl["y_1","5"]+1.25*ybinh, "female", cex = 1, pos = 1)
    text(bin_tbl["x_1","6"], bin_tbl["y_1","6"]-.57*ybinh, "male", cex = 1, pos = 1)
    
    bx0 = bin_tbl["x_1","3"]-.2
    bx1 =  bin_tbl["x_1","7"] +.2
    by0 = bin_tbl["y_1","3"]+.45*ybinh
    by1 = bin_tbl["y_1","3"]+.45*ybinh
    baj <- .05
    lines(x = c(bx0, bx0, bx1, bx1), y = c(by0, by0+baj, by0+baj, by0))
    lines(x = c(bx0, bx0, bx1, bx1), y = c(by0-1.9*ybinh, by0-1.9*ybinh-baj, by0-1.9*ybinh-baj, by0-1.9*ybinh))
    
    legend_image <- as.raster(matrix(rev(get_numeric_color(pal)), ncol=1))
    min_pos <- .13
    max_pos <- .73
    text(x= max_x-.1, y = seq(min_pos, max_pos,l=5), labels = round(seq(vlim[1],vlim[2],l=5),round.digit), cex = .8)
    rasterImage(legend_image,max_x-.4, min_pos, max_x-.3, max_pos)
    text(max_x-.2, bin_tbl["y_1","5"]+1.1*ybinh,legend.title, cex = .8, pos = 1, font=2)
}


# Creating a strip for selected gene
heat_strip <- function(gexpr = NULL, vlim = c(0, max(quantile(gexpr,1), .1)), pal = "RdYlBu", bin_tbl = NULL) {
    #hmap_img=readJPEG(empty_frame_path)
    max_x <- 10
    max_y <- 1
    min_x <- 0
    min_y <- 0
    
    # First assign correct position to the numbers
    strip_order <- c('1', '2', '3', '4', '5', '6', '7', '8', '9A', '9B', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19','22', '25','27','21','24', '20',  '23',  '26',  '28', '29', '30', '31', '32')
    gexpr <- gexpr[,strip_order]
    bin_names <- colnames(gexpr)
    
    bin_tbl <- data.frame(matrix(ncol = ncol(gexpr), nrow = 4), row.names = c("x_1","y_1","x_2","y_2"))
    colnames(bin_tbl) <- colnames(gexpr)
    y_adj_min <- .01
    y_adj_max <- -.05
    ybinh <- max_y+y_adj_max-min_y-y_adj_min
    
    x_adj_min <- .3
    x_adj_max <- -.1
    x_out <- ncol(gexpr)
    
    xbinw <- (max_x + x_adj_max- min_x-x_adj_min + 0.5)/x_out
    bin_tbl["x_1",] <- seq(min_x+x_adj_min, max_x+ x_adj_max, length.out = x_out)
    bin_tbl["x_2",] <- seq(min_x+x_adj_min, max_x+ x_adj_max, length.out = x_out)+ xbinw
    
    bin_tbl["y_1",] <- min_y
    bin_tbl["y_2",] <- min_y + ybinh
    
    cor_len = 10
    glim <- c(vlim[1], vlim[2])
    if(pal %in% c("RdYlBu", "RdBu")) {
        orig_col<-rev(RColorBrewer::brewer.pal(9,pal))
        pals <- colorRampPalette(orig_col)(cor_len)
    } else if(pal == "RedGrey"){
        orig_col<- c("grey", "red")
        pals <- colorRampPalette(orig_col)(cor_len)
    } else if(pal %in% c("viridis", "inferno")) {
        pals <- viridis(n = cor_len, option =pal)
    }
    gexpr_col <- pals[sapply(gexpr, function(g) {
        thresh_pass<-which(g >= seq(glim[1],glim[2], length.out = cor_len))
        if(!length(thresh_pass))thresh_pass<-1
        max(thresh_pass)})]
    names(gexpr_col) <- colnames(gexpr)
    
    cur_g <- rownames(gexpr)
    par(mar=c(0,0,0,0))
    plot(1:10,seq(0, max_y,length.out=10), type='n', main="", xlab="x", ylab="y", xlim = c(0,10), ylim = c(0,max_y), bty="n",xaxt="n", yaxt="n")
    #assign("bin_tbl",bin_tbl, env=.GlobalEnv)
    for(i in 1:ncol(bin_tbl)){
        cur_nm <- colnames(bin_tbl)[i]
        #if(cur_nm %in% c(20,21,24)) lty=2 else lty=1
        lty = 1
        if(cur_nm %in% c(1:19,'9A','9B',20,21,22)) bcolor = "#7C7C7C" else bcolor <- "#3A3ACE"
        rect(bin_tbl["x_1",cur_nm],bin_tbl["y_1",cur_nm],bin_tbl["x_2",cur_nm],bin_tbl["y_2",cur_nm], col=gexpr_col[cur_nm],border=bcolor, lwd=2, lty=lty)
        #text(x = bin_tbl["x_1",cur_nm] + xbinw/2, y = bin_tbl["y_1",cur_nm] + ybinh/2, cur_nm, font=2)
    }
    text(x = 0.3, y = bin_tbl["y_1",1] + ybinh/2, cur_g, font=2, pos = 2)
}




