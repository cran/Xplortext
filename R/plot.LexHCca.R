#' @rawNamespace import(ape, except = c(rotate, ladderize))
#' @rawNamespace import(ggdendro, except = c(theme_dendro))
#' @import ggplot2
#' @import ggrepel
#' @import ggforce
#' @importFrom ggforce geom_mark_ellipse geom_mark_hull
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom grDevices rainbow
#' @rawNamespace import(dendextend, except=c(cutree))
#' @export
plot.LexHCca <- function(x, type=c("map", "tree", "phylo", "clado", "radial", "fan"), 
                         plot=c("points","labels","centers"), selClust="ALL", selInd="ALL", axes=c(1,2),
                         theme=theme_bw(), palette=NULL, title=NULL, axis.title=NULL, axis.text=NULL, 
                         xlim=NULL,ylim=NULL, hvline=NULL, points=NULL, labels=NULL, centers=NULL, 
                         traject=NULL, hull=NULL, rotate=FALSE, branches=NULL,...)
{
  
  object <- x
  
# Reading information  -------------------------
  dots <- list(...)
   #  isnothing <- function(x) {is.null(x) || is.na(x) || is.nan(x) }
   # isnothing <- function(x) {is.null(x) || is.na(x) || is.nan(x) || length(x)=0}
   
  # Selection type
  str.type <- c("map","tree", "phylo", "clado", "radial", "fan")
  str.type.2 <- c("phylo", "clado", "radial", "fan")
  # Selection for type error
  
  if(!all(type %in% str.type))
    stop(paste0("type must be one of: ", str.type))
  type <- match.arg(type)
  if(type %in% str.type.2) bphylo <- TRUE else bphylo <- FALSE 
  
  nclust <- nrow(object$clust.count)
  ncases <- nrow(object$data.clust)
  xlim.seg <- xlim
  ylim.seg <- ylim
  
  if("horiz" %in% names(dots)) 
  {
    warning("horiz argument is deprecated for tree type. Please, use rotate argument")
    horiz <- dots$horiz
    rotate <- horiz
  } 
  
  if(is.null(hvline)){
    if(type=="map") hvline <- TRUE else hvline <- FALSE
  }



# Selection --------------------------
  selection <- function(sel1,xobj,bType,axx,axy) {
    xx<-""
    if(length(sel1)==1) if(sel1=="ALL") sel1 <- c(1:dim(xobj$coord)[1])
    if(length(sel1)==1) {
      if(sel1=="meta") sel1 <- "meta 3"
      if(sel1=="char") sel1 <- "char 0.05"
      xx <- gregexpr(pattern =' ',sel1)[[1]][1]-1
      xx <- substr(sel1, 1, xx)
    }
    
    
    if(xx=="coord" | xx=="cos2" | xx=="contrib" | xx=="meta" | xx=="char") 
    { 
      nc <- nchar(sel1)
      if(xx=="coord") {
        # Selection by coordinates"
        sel1 <- as.numeric(substr(sel1, 7, nc))
        # dft Coordinates for elements x two selected dimensions
        dft <- data.frame(xobj$coord[,c(axx,axy),drop=FALSE])
        # fval maximum value of element from two dimensions
        fval <- apply(dft, 1, function(z) max(abs(z)))
        ordmax <- rank(fval)
        posic <- which(ordmax > (length(fval)-sel1))
        sel1 <- rownames(xobj$coord)[posic]
      }
      
      if(xx=="contrib") {
        # Selection by contrib
        sel1 <- as.numeric(substr(sel1, 9, nc))
        dft <- data.frame(xobj$contrib[,c(axx,axy),drop=FALSE])
        fval <- apply(dft, 1, function(z) max(abs(z)))
        posic <- which(fval>= sel1)
        sel1 <- rownames(xobj$contrib)[posic]
      }
      if(xx=="cos2") {
        # Selection by cos2
        sel1 <- as.numeric(substr(sel1, 5, nc))
        dft <- data.frame(xobj$cos2[,c(axx,axy),drop=FALSE])
        fval <- apply(dft, 1, function(z) sum(z))
        posic <- which(fval>= sel1)
        sel1 <- rownames(xobj$cos2)[posic]
      }
      
      if(xx=="meta") {
        sel1 <- as.numeric(substr(sel1, 5, nc))
        sMeta <-  rownames(xobj$coord)[which(xobj$contrib[,axx]>mean(xobj$contrib[,axx])*sel1)]
        sMeta <-  c(sMeta,rownames(xobj$coord)[which(xobj$contrib[,axy]>mean(xobj$contrib[,axy])*sel1)])
        sel1 <- unique(sMeta)
      }
    } # End if(xx=="coord" ....
    if(is.character(sel1)) sel1 <- which(rownames(xobj$coord) %in% sel1)
    sel1 <- rownames(xobj$coord)[sel1]
    sel1  <- sel1[!is.na(sel1)]  
    return(sel1)
  }
  
  
  
  
# palette definition ----------
  if(is.null(palette))  palette <- rainbow(nclust)
  if(length(palette)==1) palette<- rep(palette, times=nclust)
  if(length(palette)<nclust) palette <- c(palette, rep("black", nclust-length(palette)))
  if(length(palette)>nclust) palette <- palette[1:nclust]
  
  
  
  

# Type no map -------------  
  if(type!="map") 
  {
    # if(type=="tree" || type=="phylo"|| type=="clado"|| type=="radial"|| type=="fan") 
    if(is.null(object$call$t)) stop("No dendrogram is built for kmeans or consolidation")
  palette.branches <- palette
  
  # Title for the dendrogram -------------
  if(is.null(title)) {
    if(type=="tree") title <- paste0("Cluster dendrogram. ", toupper(x$type))
  }

  
  
  
 
  # Building Dendrogram -----------
  hc <- as.hclust(x$call$t$tree)  # diana or agnes dendrogram
  groups <- stats::cutree(hc, k = nclust)  # Cut the dendrogram into k groups
  # Convert hclust object to dendrogram
  dend <- as.dendrogram(hc)  # "dendrogram"
  dend.D <- dend
  
  grL <- FALSE  # And not labels in group labels
  # Add numbers or roman letters to each cluster
  # If groupLabels=TRUE then numeric group labels will be added to each cluster. 
  # If a vector is supplied then these entries will be used as the group labels. 
  # If a function is supplied then it will be passed a numeric vector of groups (e.g. 1:5) 
  # and must return the formatted group labels. 
  

  # Height. To remove
  # if(is.null(axis.text)) axis.text <- FALSE
  

   #### Branches   -----
   # palette.seg <- "black"

  if(!is.null(branches)) if(!is.list(branches)) stop("Please, use branches argument with a list structure") 
  
  if(length(branches[["color"]])>0) {
    palette.branches <- branches[["color"]]
   
    
    if(length(palette.branches)==1) palette.branches <- rep(palette.branches , times=nclust)
    

    if(length(palette.branches)<nclust) palette.branches  <- c(palette.branches, rep("black", nclust-length(palette.branches )))
    if(length(palette.branches)>nclust) {
      palette.branches <- palette.branches[1:nclust] 
    }
  }

  
  if(length(branches[["linesize"]])>0)
    dend.D<-  dendextend::assign_values_to_branches_edgePar(dend.D,k=nclust, value = branches[["linesize"]], edgePar = "lwd")
  if(length(branches[["linetype"]])>0)
    dend.D<-  dendextend::assign_values_to_branches_edgePar(dend.D,k=nclust, value = branches[["linetype"]], edgePar = "lty")
  
 
 
   #### Labels ---------------------------------
  if(!is.null(labels)) if(!is.list(labels)) stop("Please, use labels argument with a list structure") 
  
  grL <- NULL
  if(!is.null(labels)) {
    if(!is.null(labels[["groupLabels"]])) {
      # If number of elements equal to 1
      if(length(labels[["groupLabels"]])==1){
        if(is.logical(labels[["groupLabels"]])) grL <- labels[["groupLabels"]]
        if(labels[["groupLabels"]]=="as.roman")  grL <- as.roman
        if(labels[["groupLabels"]]=="letters")  grL <- letters[1:nclust]
        if(labels[["groupLabels"]]=="LETTERS")  grL <- LETTERS[1:nclust]
      } else {
        grL <- unname(labels[["groupLabels"]])  
        # For two labels in the same box:
        # labels=list(groupLabels=c(paste0("A1","\n", "Word 2"), "B2", "C3", "D4")))

        if(length(grL)>nclust) grL <- grL[1:nclust]
        if(length(grL)<nclust) grL <- c(grL, rep("", (nclust-length(grL))))
  
      }
    } 
  } # End null(labels)
  

  

  if(!is.null(labels[["color"]])) {
    if(length(labels[["color"]])>0) {
      palette.labels <- labels[["color"]]
      if(length(palette.labels)==1) palette.labels <- rep(palette.labels , times=ncases)
      
      if(length(palette.labels)<nclust) palette.labels  <- c(palette.labels, rep("black", ncases-length(palette.labels)))
      if(length(palette.labels)>nclust) palette.labels  <- palette.labels[1:ncases]
    }
  } 
  else  palette.labels <- palette
  

  # Size of the labels

  nlcex <- 1
  if(!is.null(labels[["cex"]])) {
    if(length(labels[["cex"]])==1)  {
      if(labels[["cex"]]==0) palette.labels  <- rep("transparent",length(palette.labels)) 
        else nlcex <- labels[["cex"]]
      } else {
        nlcex <- labels[["cex"]]
        }

  }

  dend.D <- dendextend::set(dend.D, "labels_cex", nlcex)  # Adjust the size of labels  

  
  ### 
  if(length(branches[["color"]])==1)
    dend.D <- dendextend::set(dend.D, "branches_col", branches[["color"]])
    dend.D <- dendextend::color_labels(dend.D, k=nclust, col = palette.labels)
    dend.D <- dendextend::color_branches(dend.D, k=nclust, col = palette.branches, groupLabels = grL)
  

    
    
    
    
    
    
    
    
    
    
    
    
    
  
  
  #-----------------  Cluster selection ---------------------  
  if(selClust[1]!="ALL") {
    dend.D.sel <- dend.D
    dend_list <- dendextend::get_subdendrograms(dend.D, nclust) 
    str.dend <- lapply(dend_list, labels)
    remov.labels <- labels(dend.D)[!labels(dend.D) %in% unlist(str.dend[selClust])]
    dend.D <- dendextend::prune(dend= dend.D, leaves=remov.labels)
  } else dend.D.sel <- NULL
  
  if(is.null(axis.title)) axis.title <- TRUE 
  
  
  
  
# type tree -------  
  if(type=="tree") {     
    if(rotate==FALSE) rotate<-0 ; if(rotate==TRUE) rotate<- 90
    if(rotate==0) horiz <- FALSE else horiz <- TRUE
    if(rotate==0 | rotate==90)  plot(dend.D , axes=axis.title,  horiz=horiz, main=title)
    if(rotate==-90) dendextend::plot_horiz.dendrogram(dend.D, side = TRUE, main=title)
    
    ### Rectangles hull ------------------
    if(!is.null(hull)) {
      
      if(is.null(hull[["which"]]))  str.which <- 1:nclust
      else str.which <- unlist(hull[["which"]])
      
      if(is.null(hull[["color"]])) hull[["color"]] <- palette.branches[str.which]
      
      # Color selection
      str.col <- hull[["color"]]
      
      # prop_k_height  a (scalar) value (should be between 0 to 1), 
      # indicating what proportion of the height our rect will be between the height needed
      # for k and k+1 clustering. By default 0.5
      if(is.null(hull[["prop_k_height"]])) prop_k_height <- 0.5 else 
        prop_k_height <- hull[["prop_k_height"]]
      
      if(is.null(hull[["lower_rect"]])) lower_rect<- 0 else
        lower_rect <- hull[["lower_rect"]]
      
      if(is.null(hull[["upper_rect"]])) upper_rect<- 0 else
        upper_rect <- hull[["upper_rect"]]
      
      dendextend::rect.dendrogram(dend.D, k=nclust,which =str.which, border=str.col, horiz=horiz,
                                  prop_k_height= prop_k_height, lower_rect=lower_rect,
                                  upper_rect=upper_rect)
    } # End hull
    
    
    
    
    
    #### Add line ----------------------
    if(selClust[1]!="ALL") hvline <- NULL
    
    if(!is.null(hvline)) if(hvline[1]!=FALSE) {
      hv.lwd <- 2 ; hv.lty <- 2 ; hv.col <- "black"
      hv.pos <- dendextend::heights_per_k.dendrogram(dend.D)[nclust]
      
      if(!is.logical(hvline[1])) {
        if(!is.null(hvline[["linesize"]])) hv.lwd <- hvline[["linesize"]]
        if(!is.null(hvline[["linetype"]])) hv.lty <- hvline[["linetype"]]
        if(!is.null(hvline[["color"]])) hv.col <- hvline[["color"]]
        
        if(!is.null(hvline[["pos"]])) hv.pos <- hvline[["pos"]]
      }
      if(horiz==TRUE) abline(v = hv.pos, lwd = hv.lwd, lty = hv.lty, col = hv.col)
      else abline(h = hv.pos, lwd = hv.lwd, lty = hv.lty, col = hv.col)
    }
    
    
    
    
    
    
    
  } # End type tree
  
  


  #### ape required ------------
  if(bphylo) {

    if(system.file(package='ape')=="") install.packages("ape")
    
    if(is.null(dend.D.sel))
      dend.DP <- ape::as.phylo(dend.D)
    else
      dend.DP <- ape::as.phylo(dend.D.sel)
    
    # Computing color for branches in phylo graph
    df <- data.frame(dend.DP$edge)
    colnames(df) <- c("node", "case")
    
    # Adding names and groups
    df$names <- dend.DP$tip.label[df$case]
    df$gr <- x$clust.content[df$names, "Cluster"] # node case names gr
    
    # Selecting possible nodes to join
    df$sel <- TRUE
    df$sel[is.na(df$gr)] <- FALSE
    
    
    ##############################
    df$gr2 <- df$gr
    df$order <- 1:nrow(df)
    
    # Number of groups
    ngr <- sum(!is.na(unique(df$gr)))
    
    for(i in 1:ngr) {
      # Only duplicated nodes can be joined
      if(sum(df$sel & !is.na(df$gr))>ngr) {
        
        df.sel <- df[df$sel==TRUE,]
        dupl <- rownames(df.sel[duplicated(df.sel[c("node")],fromLast = F), ])
        dupl <- c(dupl,rownames(df.sel[duplicated(df.sel[c("node")],fromLast = T), ]))
        df.sel <- df[dupl,]
        u.sel <- unique(df.sel[c("node","gr")])
        
        # Put in gr.y the groups of the new nodes
        # rownames case node names gr.x   sel gr.y
        T1 <- merge(df, u.sel, by.x="case", by.y="node", all.x=TRUE)
        T1 <- T1[order(T1$order),]
        # Put in gr.x the new nodes
        T1[is.na(T1$gr.x),"gr.x"] <- T1[is.na(T1$gr.x),"gr.y"]
        # sel variable for joined cases must be FALSE now
        # Put TRUE in sel variable for not NA values in gr.y
        # u.sel$node
        
        # If case is in u.sel node sel TRUE.
        T1[T1$case %in% u.sel$node, "sel"] <- TRUE
        # If node está en u.sel node sel FALSE
        T1[T1$node %in% u.sel$node, "sel"] <- FALSE
        T1$gr.y <- NULL
        colnames(T1)[4] <- "gr"
        df <- T1
      }
    } # End i
    
    df$gr[is.na(df$gr)] <- ngr+1
    if(length(palette.branches) < ngr+1)
      palette.branches <- c(palette.branches, rep("black", (ngr+1-length(palette.branches))))
    ecn <-df$gr
    
    if(!is.null(dend.D.sel)) { # selecting cases in phylo
      # dend_list <- dendextend::get_subdendrograms(dend.D.sel, nclust)
      # str.dend <- lapply(dend_list, labels)
      # select.labels <- labels(dend.D.sel)[!labels(dend.D.sel) %in% unlist(str.dend[selClust])]
      palette.branches[!c(1:length(palette.branches)) %in% selClust] <- "transparent"
      palette.labels[!c(1:length(palette.labels)) %in% selClust] <- "transparent"
    }
    
    if(type=="phylo") ptype <- "u"
    if(type=="radial") ptype <- "r"
    if(type=="clado") ptype <- "c"
    if(type=="fan") ptype <- "fan"
    
    # rotate.tree
    # for "fan", "unrooted", or "radial" trees: the rotation of the whole tree in degrees (negative values are accepted).
    if(rotate==FALSE) rotate <- 0
    if(rotate==TRUE) rotate <- 90
    
    ape::plot.phylo(dend.DP, type=ptype, main=title,
                    use.edge.length = FALSE,
                    tip.col=palette.labels[groups],
                    edge.color=palette.branches[ecn], cex=nlcex,
                    rotate.tree=rotate,
                    edge.width=branches[["linesize"]],
                    edge.lty= branches[["linetype"]],
                    font=1)
    # nodelabels() # blue
    # tiplabels() # yellow
    # edgelabels() # green
  } # End bphylo
} # End no  map
  
  
  

#### Type map
  if(type=="map") {
    #---- Check what is plotted and limits ----
    if(plot[1]=="ALL") bplot <- TRUE else bplot <-FALSE
    bPoints <- bLabels <- bCenters <- bplot 
    if("points" %in% plot) bPoints <- TRUE
    if("labels" %in% plot) bLabels <- TRUE
    if("centers" %in% plot) bCenters <- TRUE
    if("traject" %in% plot) bTrajectory <- TRUE else bTrajectory <- FALSE
    if("hull" %in% plot) bHull <- TRUE else bHull <- FALSE
    
    if(!is.null(xlim))
      xlim <- c(min(object$coord.clust[axes[1]]), max(object$coord.clust[axes[1]]))
    if(!is.null(ylim)) 
      ylim <- c(min(object$coord.clust[axes[2]]), max(object$coord.clust[axes[2]]))
    

    
    
    #---- Main title ----
    if(!is.null(title)) if(!is.list(title)) stop("Please, use title argument with a list structure") 
    
    if(!is.list(title)) {
      s.title <- title
      title <- list(text=NULL, color=NULL, size= NULL, family=NULL,face=NULL,just=NULL)
      if(!is.null(s.title)) title[["text"]] <- s.title
    } 
    
    if(is.null(title[["text"]])) title.text <- "Clusters on the CA map" else title.text <- title[["text"]]
    
    if(is.null(title[["color"]])) title.color <- "black" else title.color <- title[["color"]]
    if(is.null(title[["size"]])) title.size <- 18 else title.size <- title[["size"]]
    if(is.null(title[["family"]])) title.family <- 'serif' else title.family <- title[["family"]]
    if(is.null(title[["face"]])) title.face <- 'plain' else title.face <- title[["face"]]  
    if(is.null(title[["just"]])) title.hjust <- "c" else title.hjust <- title[["just"]] 
    if(title.hjust=="c"| title.hjust=="centered") title.hjust <- 0.5
    if(title.hjust=="l"| title.hjust=="left") title.hjust <- 0
    if(title.hjust=="r"| title.hjust=="right") title.hjust <- 1
    
    
    
    
    #---- Axestitle ----
    if(!is.list(axis.title)) title <- list(text.x= NULL, text,y=NULL, color=NULL, size= NULL, family=NULL,
                                           face=NULL,just=NULL)
    
    if(is.null(axis.title[["text.x"]]))
      labx <- paste0("Dim ", axes[1], " (", round(object$call$CA$eig[axes[1],2],2), "%)")
    else labx <- axis.title[["text.x"]]
    
    if(is.null(axis.title[["text.y"]]))
      laby <- paste0("Dim ", axes[2], " (", round(object$call$CA$eig[axes[2],2],2), "%)")
    else laby <- axis.title["text.y"]
    
    
    if(is.null(axis.title[["color"]])) axis.title.color <- "black" else axis.title.color <- axis.title[["color"]]
    if(is.null(axis.title[["size"]]))  axis.title.size <- 12 else axis.title.size <- axis.title[["size"]]
    if(is.null(axis.title[["family"]])) axis.title.family <- 'serif' else axis.title.family <- title[["family"]]
    if(is.null(axis.title[["face"]])) axis.title.face <- 'plain' else axis.title.face <- axis.title[["face"]]
    if(is.null(axis.title[["just"]])) axis.title.hjust <- "c" else axis.title.hjust <- axis.title[["just"]]  
    if(axis.title.hjust=="c"| axis.title.hjust=="centered") axis.title.hjust <- 0.5
    if(axis.title.hjust=="l"| axis.title.hjust=="left") axis.title.hjust <- 0
    if(axis.title.hjust=="r"| axis.title.hjust=="right") axis.title.hjust <- 1
    axis.title <- element_text(color = axis.title.color, family=axis.title.family, face=axis.title.face,
                               hjust =axis.title.hjust, size=axis.title.size)
    
    
    # Axes text/numbers
    if(is.null(axis.text[["color"]])) axis.text.color <- "black" else axis.text.color <- axis.text[["color"]]
    if(is.null(axis.text[["size"]])) axis.text.size <- 8 else axis.text.size <- as.numeric(axis.text[["size"]])
    if(is.null(axis.text[["family"]])) axis.text.family <- 'serif' else axis.text.family <- axis.text[["family"]]
    if(is.null(axis.text[["face"]])) axis.text.face <- 'plain' else axis.text.face <- axis.text[["face"]]
    
    
 
    #---- Starting plot --------------
    
    nclust <- nrow(object$clust.count)
    ncases <-  nrow(object$data.clust)
    if(max(axes[1],axes[2]) > ncol(object$centers)-1) stop("The number of axes is bigger than saved dimensions")
    
    dataT  <- object$coord.clust                # Total data
    dataT$ID_ <- rownames(dataT)
    dataC  <-  object$centers               # Centers
    dataC$ID_ <- paste0("cl.",dataC$Clust_)
    dataC$Clust_C <- dataC$Clust_
    dataC <- dataC[,-1] 
    
    # if(is.null(palette)) palette<- rainbow(nclust)
    if(selClust[1]!="ALL") palette[!(1:nclust) %in% selClust] <- "transparent"
    
    if(deparse(substitute(theme))=="theme_bw()") 
      theme <- theme+ theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank())
    
    theme <- theme+
      theme(axis.text= element_text(color=axis.text.color, size=axis.text.size,family=axis.text.family,
                                    face=axis.text.face),
            axis.title.x = axis.title,
            axis.title.y = axis.title,
            plot.title = element_text(hjust = title.hjust, colour=title.color, size=title.size, 
                                      family=title.family, face=title.face)
      )
    p <- ggplot(data=dataT, aes(x=dataT[,axes[1]], y=dataT[,axes[2]])) +  theme +
      coord_equal() + labs(x=labx)+labs(y=laby) + ggtitle(title.text)
    
    
    
    
    
    
    #---- hvline -----------------
    if(!is.null(hvline)) if(hvline[[1]]!=FALSE) {
      if(!is.list(hvline)) {
        hvline <- list(intercept.x=NULL, intercept.y=NULL, linetype=NULL, color=NULL, linesize= NULL,
                       family=NULL,face=NULL,just=NULL, alpha.t=NULL)
        intercept.y <- intercept.x <- 0
        intercept.linetype <- 2 # Dashed
        intercept.color <- "grey"
        intercept.size <- .5
        intercept.alpha <- 1
      } else {  # Is list
        if(is.null(hvline[["intercept.y"]])) intercept.y<-0 else  intercept.y <- as.numeric(hvline[["intercept.y"]])
        if(is.null(hvline[["intercept.x"]])) intercept.x<-0 else  intercept.x <- as.numeric(hvline[["intercept.x"]])
        if(is.null(hvline[["linetype"]])) intercept.linetype<-"dashed" else  intercept.linetype <- hvline[["linetype"]]
        if(is.null(hvline[["color"]])) intercept.color <-"gray" else  intercept.color <- hvline[["color"]]
        if(is.null(hvline[["linesize"]])) intercept.size <- .5 else  intercept.size <- as.numeric(hvline[["linesize"]])
        if(is.null(hvline[["alpha.t"]])) intercept.alpha <- 1 else  intercept.alpha <- as.numeric(hvline[["alpha.t"]])
      }
      
      p <- p + geom_hline(yintercept=intercept.y, linetype=intercept.linetype, color = intercept.color,
                          linewidth=intercept.size, alpha=intercept.alpha)
      p <- p + geom_vline(xintercept=intercept.x, linetype=intercept.linetype, color = intercept.color,
                          linewidth=intercept.size, alpha=intercept.alpha)
    }
    
    
    
    
    #---- Selection of cases. by cos2, contrib,... ----    
    if(is.null(selInd[1])) {
      bPoints <- bLabels <- FALSE
    } else {
      if(object$call$cluster.CA=="docs")  selInd <-  selection(selInd, object$call$CA$row,"Doc", axes[1],axes[2])
      else selInd <-  selection(selInd, object$call$CA$col,"Word", axes[1],axes[2])
    } # End if(is.null(selInd[1]))
    
    ncases.sel <- length(selInd)    
    if(ncases.sel==0) stop("No selected cases")
    # df.sel <- dataK$data.clust[selInd,] 
    df.sel <- dataT[selInd,] 
    
    
    
    
    #---- Only print the selected cases, do not fix them transparent  ---------------
    if(is.null(points[["shape"]])) points.shape<- 21 else points.shape<- as.numeric(points[["shape"]])
    # Shape between 0 to 20 do not fill the symbol    
    if(is.null(points[["size"]])) points.size <- 2 else points.size<- as.numeric(points[["size"]])
    if(is.null(points[["fill"]])) points.fill <- palette[df.sel$Clust_] else points.fill<- points[["fill"]]
    if(is.null(points[["stroke"]])) points.stroke<-0 else points.stroke<- as.numeric(points[["stroke"]])
    
    

    if(is.null(points[["border"]])) points.border <-"transparent" # palette[df.sel$Clust_] 
    else 
    {
      points.border<- points[["border"]]
      if(points.stroke==0) points.stroke<-1
    }
    
    if(is.null(points[["alpha.t"]])) points.alpha<-1 else points.alpha <- as.numeric(points[["alpha.t"]])
    if(is.null(points.size) | points.size==0 )  bPoints <- FALSE
    if(bPoints==FALSE) points.size <- NA
    
    
    if(bPoints)
      p <- p + geom_point(data=df.sel, aes(x=df.sel[,axes[1]], y=df.sel[,axes[2]]),
                          shape=points.shape, size=points.size, fill=points.fill , color=points.border,
                          stroke=points.stroke, alpha=points.alpha, show.legend=FALSE)
    
    
    
    
    #---- Labels  ----  
    if(!is.null(labels)) if(!is.list(labels)) stop("Please, use labels argument with a list structure") 
    
    if(bLabels) {
      if(is.null(labels[["numbers"]])) bNumbers <- FALSE else bNumbers=labels[["numbers"]]
      if(is.null(labels[["size"]])) labels.size <- 4 else labels.size<- as.numeric(labels[["size"]])
      if(labels.size==0) bLabels<-FALSE
      if(is.null(labels[["family"]])) labels.family <- 'serif' else labels.family <- labels[["family"]]
      if(is.null(labels[["face"]])) labels.face <- 'plain' else labels.face <- labels[["face"]]
      if(is.null(labels[["hjust"]])) labels.hjust <-1  else labels.hjust <- labels[["hjust"]]
      if(is.null(labels[["vjust"]])) labels.vjust <-1  else labels.vjust <- labels[["vjust"]]
      if(is.null(labels[["max.overlaps"]])) max.overlaps <-10  else max.overlaps <- labels[["max.overlaps"]]
      
      if(is.null(labels[["color"]])) labels.color.text <- points.fill else labels.color.text<-labels[["color"]]
      if(is.null(labels[["alpha.t.text"]])) labels.alpha.text<-1 else labels.alpha.text <- as.numeric(labels[["alpha.t.text"]])
      if(labels.alpha.text!=1) {
        labels.color.text <- apply(sapply(labels.color.text, col2rgb)/255, 2,  function(x) 
          rgb(x[1], x[2], x[3], alpha=labels.alpha.text))} 
      
      if(is.null(labels[["color.fill"]])) labels.color.fill <- "transparent" else labels.color.fill=labels[["color.fill"]]
      if(labels.color.fill[1]==FALSE)  labels.color.fill <- "transparent"
      if(labels.color.fill[1]==TRUE)  labels.color.fill <- points.fill
      
      if(is.null(labels[["alpha.t.fill"]])) labels.alpha.fill<-1 else labels.alpha.fill <- as.numeric(labels[["alpha.t.fill"]])
      if(labels.alpha.fill!=1) {
        labels.color.fill <- apply(sapply(labels.color.fill, col2rgb)/255, 2,  function(x) 
          rgb(x[1], x[2], x[3], alpha=labels.alpha.fill))} 
      
      if(is.null(labels[["rect"]])) labels.rect <- FALSE else labels.rect<- labels[["rect"]]
      if(is.null(labels[["force"]])) labels.force <-1  else labels.force <- as.numeric(labels[["force"]])
      if(labels.force==FALSE) labels.force <-0
      if(labels.force>0) { labels.hjust <- NULL ; labels.vjust <- NULL}
      if(is.null(labels[["set.seed"]])) labels.seed <- set.seed(Sys.time())  else labels.seed <- labels[["set.seed"]]
    } # End bLabels
    
    if(bLabels) {
      if(bNumbers) IDL <- df.sel$Clust_ else IDL <- df.sel$ID_
      if(labels.force!=0)
        if(labels.rect) {
          p <- p + ggrepel::geom_label_repel(data=df.sel, aes(x=df.sel[,axes[1]], y=df.sel[,axes[2]],label = IDL,
                                                              hjust=labels.hjust, vjust=labels.vjust, ),
                                             size=labels.size, family = labels.family, fontface = labels.face,
                                             colour=labels.color.text, fill= labels.color.fill,force=labels.force,
                                             max.overlaps = max.overlaps,
                                             seed=labels.seed) 
        } else {
          # Label without border
          p <- p + ggrepel::geom_text_repel(data=df.sel, aes(x=df.sel[,axes[1]], y=df.sel[,axes[2]],label = IDL,
                                                             hjust=labels.hjust, vjust=labels.vjust),
                                            size=labels.size, family = labels.family, fontface = labels.face,
                                            force=labels.force, colour=labels.color.text,
                                            max.overlaps = max.overlaps,
                                            seed=labels.seed)
        }
      else { # No repel
        if(labels.rect) {
          p <- p + geom_label(data=df.sel, aes(x=df.sel[,axes[1]], y=df.sel[,axes[2]],label = IDL,
                                               hjust=labels.hjust, vjust=labels.vjust, ),
                              size=labels.size, family = labels.family, fontface = labels.face,
                              colour=labels.color.text, fill= labels.color.fill) 
          
        } else {
          p <- p + geom_text(data=df.sel, aes(x=df.sel[,axes[1]], y=df.sel[,axes[2]],label = IDL,
                                              hjust=labels.hjust, vjust=labels.vjust),
                             size=labels.size, family = labels.family, fontface = labels.face,
                             colour=labels.color.text)
        }
        
      }
    } # End bLabels  
    
    
    
    
    
    
    #---- trajectory ------------
    if(bTrajectory) {
      if(is.null(traject[["color"]])) traject.color <- "blue" else traject.color=traject[["color"]]
      if(is.null(traject[["space"]])) traject.space <- 0 else traject.space <- as.numeric(traject[["space"]])
      if(is.null(traject[["size"]])) traject.size <- 1 else traject.size<-as.numeric(traject[["size"]])
      if(is.null(traject[["linetype"]])) traject.linetype <- 'solid' else traject.linetype<- traject[["linetype"]]
      if(is.null(traject[["arrow.length"]])) arrow.length <- .3 else   arrow.length<- as.numeric(traject[["arrow.length"]])
      if(is.null(traject[["arrow.type"]])) arrow.type <- "closed" else   arrow.type<-traject[["arrow.type"]]
      if(is.null(traject[["arrow.angle"]])) arrow.angle <- 30 else arrow.angle<-as.numeric(traject[["arrow.angle"]])
      if(is.null(traject[["alpha.t"]])) traject.alpha<-1 else traject.alpha <-as.numeric(traject[["alpha.t"]])

       ntr <- nrow(df.sel)
      dftr <- data.frame("x1"=df.sel[-ntr,axes[1]],"y1"=df.sel[-ntr,axes[2]],
                         "xend"=df.sel[-1,axes[1]],"yend"=df.sel[-1,axes[2]] )
      sc.x <-traject.space * max(df.sel[,axes[1]])- min(df.sel[,axes[1]])/10000
      sc.y <-traject.space * max(df.sel[,axes[2]])- min(df.sel[,axes[2]])/10000
      modul <- sqrt((dftr$xend-dftr$x1)^2+(dftr$yend-dftr$y1)^2)
      
      dftrj <- data.frame("x1"=dftr$x1 +sc.x* (dftr$xend-dftr$x1)/modul)
      dftrj$y1 <- dftr$y1 +sc.y* (dftr$yend-dftr$y1)/modul
      v.xend <- dftr$xend -sc.x*(dftr$xend-dftr$x1)/modul
      v.yend <- dftr$yend -sc.y* (dftr$yend-dftr$y1)/modul
      
      if(traject.alpha!=1) {
        traject.color <- apply(sapply(traject.color, col2rgb)/255, 2,  function(x) 
          rgb(x[1], x[2], x[3], alpha=traject.alpha ))  
      } 
      # dftrj <- fortify(dftrj)
      xinit <- dftr$x1
      yinit <- dftr$y1
      
      p<- p+ ggplot2::geom_segment(data=dftrj, aes(x=xinit, y=yinit, xend=v.xend, yend=v.yend), 
                                   arrow=arrow(length=unit(arrow.length,"cm"), type=arrow.type, 
                                               angle=arrow.angle), arrow.fill =traject.color,
                                   linewidth=traject.size, colour=traject.color,
                                   linetype=traject.linetype, linejoin = "mitre") 
      }
    
    
    
    #---- Centers ----
    if(bCenters) {
      if(!is.null(centers)) if(!is.list(centers)) stop("Please, use centers argument with a list structure") 
      
      centers.labels <- dataC$ID_
      if(!is.null(centers[["labels"]])) {
        nc  <- length(centers[["labels"]])
        nl  <- length(centers.labels)
        if(nl > nc) centers.labels[1:nc] <-  centers[["labels"]] else
          centers.labels <- centers[["labels"]][1:nl] }
      
      centers.fill <- rep("transparent", nclust)
      if(!is.null(centers[["fill"]])) {
        nc  <- length(centers[["fill"]])
        if(nc==1) { centers.fill[1:nclust] <- centers[["fill"]][1] ; nc <- nclust } else
          centers.fill[1:nclust] <- centers[["fill"]][1:nclust]} 
      

      centers.color <- NULL
      if(!is.null(centers[["color"]])) {
          nc  <- length(centers[["color"]])
         if(nc==1) centers.color[1:nclust] <- centers[["color"]]
        else
        { centers.color[1:nclust] <- centers[["color"]][1:nclust]
        centers.color[is.na(centers.color)] <- "black"}
      } else centers.color <- palette
 

      if(is.null(centers[["size"]])) centers.size <- 5 else 
      {
        centers.size <- centers[["size"]]  
        
      }

      
      if(!is.null(centers[["alpha.t"]])) {
        centers.alpha <- centers[["alpha.t"]]
        centers.fill <- apply(sapply(centers.fill, col2rgb)/255, 2,  function(x) 
          rgb(x[1], x[2], x[3], alpha=centers.alpha)) 
      } 
      
      if(is.null(centers[["family"]])) centers.family <- 'serif' else centers.family <- centers[["family"]]
      
      if(is.null(centers[["face"]])) centers.face <- 3 else centers.face <- centers[["face"]]
      if(centers.face=="italic") centers.face <- 3
      
      p <- p +  ggplot2::geom_label(data=dataC,
                                    aes(x=dataC[,axes[1]], y=dataC[,axes[2]]), 
                                    fill=centers.fill,
                                    size = centers.size,
                                    family = centers.family,
                                    fontface = centers.face,
                                    colour= centers.color,
                                    label = centers.labels)
      
    }  # End centers
    
     
    
    
    #---- Hulls ----
    if(bHull) {
       if(is.null(hull[["type"]])) hull.type <- "ellipse" else hull.type <-hull["type"]
      if(is.null(hull[["alpha.t"]])) hull.alpha <- .1 else hull.alpha <- as.numeric(hull[["alpha.t"]])
      if(is.null(hull[["color"]])) hull.color <- "black"  else hull.color<-hull[["color"]]  
      if(is.null(hull[["linetype"]])) hull.linetype <- "dotted"  else hull.linetype<- hull[["linetype"]]
      
      
      v.Clust_ <- as.factor(dataT$Clust_)
      if (hull.type =="ellipse") 
        p <- p + ggforce::geom_mark_ellipse(data=dataT, aes(x= dataT[,axes[1]], y= dataT[,axes[2]],
                                                            fill=v.Clust_), alpha= hull.alpha, color=hull.color, linetype=hull.linetype,
                                            show.legend =FALSE)
      if (hull.type =="hull")   
        
        p <- p + ggforce::geom_mark_hull(data=dataT, aes(x= dataT[,axes[1]], y= dataT[,axes[2]],
                                                         fill=v.Clust_), alpha= hull.alpha, color=hull.color, linetype=hull.linetype,
                                         show.legend =FALSE)   
      if (is.null(xlim))  xlim <-ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
      if (!is.null(ylim)) ylim(ylim) else ylim <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
    } # End hull
    
    
    
    
    
  }#  End map ---------------
    
  #---- return
  if(type=="map") {
  suppressMessages(p <- p + ggplot2::coord_cartesian(xlim = xlim.seg, ylim = ylim.seg))
  p <- p + theme(legend.position = "none") 
  return(p)
}
} # End function
  
  
 
  
  
  
  
