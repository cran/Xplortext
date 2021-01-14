#' import ggplot2
#' import ggrepel
#' import ggforce
#' @importFrom ggforce geom_mark_ellipse geom_mark_hull
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' #' @export
plot.LexHCca <- function(x, type="map", plot=c("points","labels","centers"),
                             selClust="ALL", selInd="ALL", axes=c(1,2), theme=theme_bw(), palette=NULL,
                             title=NULL, axis.title=NULL, axis.text=NULL,
                             points=NULL, labels=NULL, centers=NULL, traject=NULL, hull=NULL,
                               xlim=NULL,ylim=NULL,hvline=TRUE,...) 
{
  object <- x
  if(type=="tree")
  {
    if(is.null(object$call$t)) stop("No dendrogram is built for kmeans or consolidation")
    nclust<- object$call$t$nb.clust
    if(is.null(palette)) palette<- rainbow(nclust)
    plot(object$call$t$tree, hang=-1)
    rect.hclust(object$call$t$tree, k=nclust, border=palette)
  } else {
  # draw map   
 
  
  
  
  
  #---- Check what is plotted ----
  if(plot[1]=="ALL") bplot <- TRUE else bplot <-FALSE
  bPoints <- bLabels <- bCenters <- bplot 
  if("points" %in% plot) bPoints <- TRUE
  if("labels" %in% plot) bLabels <- TRUE
  if("centers" %in% plot) bCenters <- TRUE
  if("traject" %in% plot) bTrajectory <- TRUE else bTrajectory <- FALSE
  if("hull" %in% plot) bHull <- TRUE else bHull <- FALSE

  
 # stop(paste0(bPoints, bLabels, bCenters, bTrajectory))
  
  isnothing = function(x) {is.null(x) || is.na(x) || is.nan(x) }
  
  #---- Main title ----
  if(isnothing(title["text"])) title.text <- "Clusters on the CA map" else title.text <- title["text"]
  if(isnothing(title["color"])) title.color <- "black" else title.color <- title["color"]
  if(isnothing(title["size"])) title.size <- 18 else title.size <- as.numeric(title["size"])
  if(isnothing(title["family"])) title.family <- 'serif' else title.family <- title["family"]
  if(isnothing(title["face"])) title.face <- 'plain' else title.face <- title["face"]  
  if(isnothing(title["just"])) title.hjust <- "c" else title.hjust <- title["just"] 
  if(title.hjust=="c"| title.hjust=="centered") title.hjust <- 0.5
  if(title.hjust=="l"| title.hjust=="left") title.hjust <- 0
  if(title.hjust=="r"| title.hjust=="right") title.hjust <- 1


  # Axes. axis.title
  if(isnothing(axis.title["text.x"])) labx <- paste0("Dim ", axes[1], " (", round(object$call$CA$eig[axes[1],2],2), "%)")
  else labx <- axis.title["text.x"]
  if(isnothing(axis.title["text.y"])) laby <- paste0("Dim ", axes[2], " (", round(object$call$CA$eig[axes[2],2],2), "%)")
  else laby <- axis.title["text.y"]
  if(isnothing(axis.title["color"])) axis.title.color <- "black" else axis.title.color <- axis.title[["color"]]
  if(isnothing(axis.title["size"]))  axis.title.size <- 12 else axis.title.size <- as.numeric(axis.title["size"])
  if(isnothing(axis.title["family"])) axis.title.family <- 'serif' else axis.title.family <- title["family"]
  if(isnothing(axis.title["face"])) axis.title.face <- 'plain' else axis.title.face <- axis.title["face"]
    if(isnothing(axis.title["just"])) axis.title.hjust <- "c" else axis.title.hjust <- axis.title["just"]  
    if(axis.title.hjust=="c"| axis.title.hjust=="centered") axis.title.hjust <- 0.5
    if(axis.title.hjust=="l"| axis.title.hjust=="left") axis.title.hjust <- 0
    if(axis.title.hjust=="r"| axis.title.hjust=="right") axis.title.hjust <- 1
    axis.title <- element_text(color = axis.title.color, family=axis.title.family, face=axis.title.face,
                            hjust =axis.title.hjust)
  
  
  # Axes text/numbers
   if(isnothing(axis.text["color"])) axis.text.color <- "black" else axis.text.color <- axis.text["color"]
   if(isnothing(axis.text["size"])) axis.text.size <- 8 else axis.text.size <- as.numeric(axis.text["size"])
   if(isnothing(axis.text["family"])) axis.text.family <- 'serif' else axis.text.family <- axis.text["family"]
   if(isnothing(axis.text["face"])) axis.text.face <- 'plain' else axis.text.face <- axis.text["face"]

    
    # Starting plot 
    
    nclust <- nrow(object$clust.count)
    ncases <-  nrow(object$data.clust)
        if(max(axes[1],axes[2]) > nclust) stop("The number of axes is bigger than saved dimensions")


     dataT  <- object$coord.clust                # Total data
     dataT$ID_ <- rownames(dataT)
     dataC  <-  object$centers               # Centers
     dataC$ID_ <- paste0("cl.",dataC$Clust_)
     dataC$Clust_C <- dataC$Clust_
     dataC <- dataC[,-1] 

    
    
    if(is.null(palette)) palette<- rainbow(nclust)
    if(selClust[1]!="ALL") palette[!(1:nclust) %in% selClust] <- "transparent"
    
    
    if(!is.null(hvline)) if(hvline[[1]]==FALSE) hvline<-NULL
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
      coord_equal() + labs(x=labx)+labs(y=laby) + ggtitle (title.text) 
    
    
    if(!isnothing(hvline))  {
      if(isnothing(hvline["intercept.y"])) intercept.y<-0 else  intercept.y <- as.numeric(hvline["intercept.y"])
      if(isnothing(hvline["intercept.x"])) intercept.x<-0 else  intercept.x <- as.numeric(hvline["intercept.x"])
      if(isnothing(hvline["linetype"])) intercept.linetype<-"dashed" else  intercept.linetype <- hvline["linetype"]
      if(isnothing(hvline["color"])) intercept.color <-"gray" else  intercept.color <- hvline["color"]
      if(isnothing(hvline["size"])) intercept.size <- .5 else  intercept.size <- as.numeric(hvline["size"])
      if(isnothing(hvline["alpha.t"])) intercept.alpha <- 1 else  intercept.alpha <- as.numeric(hvline["alpha.t"])
      
      p <- p + geom_hline(yintercept=intercept.y, linetype=intercept.linetype, color = intercept.color,
                          size=intercept.size, alpha=intercept.alpha)
      p <- p + geom_vline(xintercept=intercept.x, linetype=intercept.linetype, color = intercept.color,
                          size=intercept.size, alpha=intercept.alpha)
    }
  
    
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


#### Only print the selected cases, do not fix them transparent
if(isnothing(points["shape"])) points.shape<- 21 else points.shape<- as.numeric(points["shape"])
# Shape between 0 to 20 do not fill the symbol    
if(isnothing(points["size"])) points.size <- 2 else points.size<- as.numeric(points["size"])
if(isnothing(points["fill"])) points.fill <- palette[df.sel$Clust_] else points.fill<- points["fill"]
if(isnothing(points["border"])) points.border <- palette[df.sel$Clust_] else points.border<- points["border"]
if(isnothing(points["stroke"])) points.stroke<-0 else points.stroke<- as.numeric(points["stroke"])
if(isnothing(points["alpha.t"])) points.alpha<-1 else points.alpha <- as.numeric(points["alpha.t"])
if(isnothing(points.size) | points.size==0 )  bPoints <- FALSE
if(bPoints==FALSE) points.size <- NA

if(bPoints)
p <- p + geom_point(data=df.sel, aes(x=df.sel[,axes[1]], y=df.sel[,axes[2]]),
                    shape=points.shape, size=points.size, fill=points.fill , color=points.border,
                    stroke=points.stroke, alpha=points.alpha, show.legend=FALSE)


#---- Labels  ----  
if(bLabels) {
if(isnothing(labels["numbers"])) bNumbers <- FALSE else bNumbers=labels["numbers"]
if(isnothing(labels["size"])) labels.size <- 4 else labels.size<- as.numeric(labels["size"])
if(labels.size==0) bLabels<-FALSE
if(isnothing(labels["family"])) labels.family <- 'serif' else labels.family <- labels["family"]
if(isnothing(labels["face"])) labels.face <- 'plain' else labels.face <- labels["face"]
if(isnothing(labels["hjust"])) labels.hjust <-1  else labels.hjust <- labels["hjust"]
if(isnothing(labels["vjust"])) labels.vjust <-1  else labels.vjust <- labels["vjust"]
if(isnothing(labels["max.overlaps"])) max.overlaps <-10  else max.overlaps <- labels["max.overlaps"]



if(isnothing(labels["color.text"])) labels.color.text <- points.fill else labels.color.text<-labels["color.text"]
if(isnothing(labels["alpha.t.text"])) labels.alpha.text<-1 else labels.alpha.text <- as.numeric(labels["alpha.t.text"])
if(labels.alpha.text!=1) {
  labels.color.text <- apply(sapply(labels.color.text, col2rgb)/255, 2,  function(x) 
    rgb(x[1], x[2], x[3], alpha=labels.alpha.text))} 

if(isnothing(labels["color.fill"])) labels.color.fill <- "transparent" else labels.color.fill=labels["color.fill"]
if(labels.color.fill[1]==FALSE)  labels.color.fill <- "transparent"
if(labels.color.fill[1]==TRUE)  labels.color.fill <- points.fill

if(isnothing(labels["alpha.t.fill"])) labels.alpha.fill<-1 else labels.alpha.fill <- as.numeric(labels["alpha.t.fill"])
if(labels.alpha.fill!=1) {
  labels.color.fill <- apply(sapply(labels.color.fill, col2rgb)/255, 2,  function(x) 
    rgb(x[1], x[2], x[3], alpha=labels.alpha.fill))} 

if(isnothing(labels["rect"])) labels.rect <- FALSE else labels.rect<- labels["rect"]
if(isnothing(labels["force"])) labels.force <-1  else labels.force <- as.numeric(labels[["force"]])
if(labels.force==FALSE) labels.force <-0
if(labels.force>0) { labels.hjust <- NULL ; labels.vjust <- NULL}
if(isnothing(labels["set.seed"])) labels.seed <- set.seed(Sys.time())  else labels.seed <- labels["set.seed"]

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


# trajectory
if(bTrajectory) {
  if(isnothing(traject["color"])) traject.color <- "blue" else traject.color=traject["color"]
  if(isnothing(traject["space"])) traject.space <- 0 else traject.space <- as.numeric(traject["space"])
  if(isnothing(traject["size"])) traject.size <- 1 else traject.size<-as.numeric(traject["size"])
  if(isnothing(traject["linetype"])) traject.linetype <- 'solid' else traject.linetype<- traject["linetype"]
  if(isnothing(traject["arrow.length"])) arrow.length <- .3 else   arrow.length<- as.numeric(traject["arrow.length"])
  if(isnothing(traject["arrow.type"])) arrow.type <- "closed" else   arrow.type<-traject["arrow.type"]
  if(isnothing(traject["arrow.angle"])) arrow.angle <- 30 else arrow.angle<-as.numeric(traject["arrow.angle"])
  if(isnothing(traject["alpha.t"])) traject.alpha<-1 else traject.alpha <-as.numeric(traject["alpha.t"])
}


# Starting trajectory    
if(bTrajectory) {
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
                   size=traject.size, colour=traject.color,
                  linetype=traject.linetype, linejoin = "mitre")  
}  # End trajectory




if(bCenters) {
  if(isnothing(centers["size"])) centers.size <- 5 else centers.size <- as.numeric(centers["size"])    
  if(isnothing(centers["family"])) centers.family <- 'serif' else centers.family <- centers["family"]
  if(isnothing(centers["face"])) centers.face <- 'italic' else centers.face <- centers["face"] 

  if(isnothing(centers["color"])) centers.color <- palette[dataC$Clust_C] else centers.color <- centers["color"] 
  
  if(isnothing(centers["alpha.t"])) centers.alpha<-1 else centers.alpha <-as.numeric(centers["alpha.t"])

  if(isnothing(centers["labels1"])) centers.labels <- dataC$ID_ else
  {centers.labels<- vector()
    for(i in 1:nclust) {
      namecl <- paste0("labels",i)
      centers.labels[i] <- centers[namecl]
    }
  }
  
  

  if(isnothing(centers["fill"])) { centers.fill <- rep("transparent",times=nrow(dataC));
                                                       centers.alpha <-1} else 
  {
    if(centers["fill"][1]==TRUE) centers.fill <-  palette else  centers.fill <- centers["fill"] 
  }

  if(centers.alpha!=1) {
     centers.fill <- apply(sapply(centers.fill, col2rgb)/255, 2,  function(x) 
      rgb(x[1], x[2], x[3], alpha=centers.alpha))  
  }
  p <- p +  ggplot2::geom_label(data=dataC,
                                aes(x=dataC[,axes[1]], y=dataC[,axes[2]]), fill=centers.fill,
                                size = centers.size,
                                family = centers.family,
                                fontface = centers.face,
                                colour=centers.color,
                                label = centers.labels)
                                
}  # End centers


#---- Hulls ----
if(bHull) {
  if(isnothing(hull["type"])) hull.type <- "ellipse" else hull.type<-hull["type"]
  if(isnothing(hull["alpha.t"])) hull.alpha <- .1 else hull.alpha <- as.numeric(hull["alpha.t"])
  if(isnothing(hull["color"])) hull.color <- "black"  else hull.color<-hull["color"]  
  if(isnothing(hull["linetype"])) hull.linetype <- "dotted"  else hull.linetype<- hull["linetype"]
  

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





suppressMessages(p <- p + coord_cartesian(xlim = xlim, ylim = ylim))
p <- p + theme(legend.position = "none") 
return(p)
}
}

 