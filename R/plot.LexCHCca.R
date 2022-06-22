#' @importFrom grDevices palette
#' @importFrom graphics layout legend par plot text
#' @export
plot.LexCHCca <-function(x, axes=c(1, 2), type=c("tree","map","bar"), rect=TRUE, 
      title=NULL, ind.names=TRUE, new.plot=FALSE, max.plot=15, 
      tree.barplot=TRUE, ...) 
{
  if(is.null(x)) stop("Missing argument for x object")
  if (!inherits(x,"LexCHCca")) stop("x object should be LexCHCca class")   
  dots <- list(...)
 if(!is.null(dots$choice)) stop("Change choice argument by type argument")
 if(length(type)>1) type <- type[1]
  res <- x
  X <- res$call$X
  max<-res$call$max
  min<-res$call$min
  max.plot <- max(res$call$max,max.plot)
  nb.clust <- length(levels(X$clust))
  levs <- levels(X$clust)

  par(mfrow=c(1,1))
  
  if (type == "tree") {
    if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
      dev.new()
    
    if (is.null(title)) 
      title = "Constrained hierarchical clustering"
    if (tree.barplot) {
      def.par <- par(no.readonly = TRUE)
      par(mar = c(0.5, 2, 0.75, 0))
      lay = matrix(ncol = 5, nrow = 5, c(2, 4, 4, 4, 4, 
                                         2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 
                                         1, 3, 3, 3, 3))
      layout(lay, respect = TRUE)
      inv.height = rev(res$call$t$height)
      vec<-inv.height[1:max.plot]
      barplot(height = vec, col = c(rep("black", nb.clust - 1), rep("grey", max(max, max.plot) - nb.clust + 1)), 
              space = 0.7)
      plot(x = 1, xlab = "", ylab = "", main = "", col = "white", 
           axes = FALSE)
      text(1, 1, title, cex = 2)
      plot(x = 1, xlab = "", ylab = "", main = "", col = "white", 
           axes = FALSE)
      legend("top", "Agg. criterion", box.lty = NULL, cex = 1)
    }
    

    plot(res$call$t$tree, hang = -1, xlab = "", sub = "", main=title,...)
    
    if (rect) {
      y = (res$call$t$tree$height[length(res$call$t$tree$height) - nb.clust + 2] + 
             res$call$t$tree$height[length(res$call$t$tree$height) - nb.clust + 1])/2
      ordColo <- unique(res$call$X$clust[res$call$t$tree$order])
      rect = rect.hclust(res$call$t$tree, h = y, border = ordColo)
    }
  }
  
  if (type == "map") {
    if (!res$call$vec) {
      if (is.null(title)) 
        title = "CA map"
      Y = X[, -ncol(X)]
      leg.map = NULL
      for (p in 1:nrow(X)) leg.map[p] = paste("cluster", 
                                              X$clust[p], " ", sep = " ")
      Y = cbind.data.frame(Y, as.factor(leg.map))
      res2 = PCA(Y, quali.sup = ncol(Y), scale.unit = FALSE, 
                 row.w = res$call$t$res$call$row.w, ncp = Inf, 
                 graph = FALSE)
      res2$eig <- res$call$t$res$eig
      if (ind.names) 
        plot.PCA(res2, title = title, habillage = ncol(Y), 
                 cex = 0.8, axes = axes, new.plot = new.plot, 
                 palette = palette(), graph.type = "classic", ...)
      else plot.PCA(res2, title = title, habillage = ncol(Y), 
                    cex = 0.8, axes = axes, label = "none", new.plot = new.plot, 
                    palette = palette(), graph.type = "classic", ...)
    }
  }
  if (type == "bar") {
    if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
      dev.new()
    inv.height = rev(res$call$t$height)
    vec<-inv.height[1:max.plot]
    names.arg = NULL
    if (is.null(title)) 
      title = "Aggregation criterion"
    for (i in 1:(max.plot)) names.arg = c(names.arg, paste(i, 
                                                           "-", i + 1, sep = ""))
    barplot(vec, names.arg = names.arg, col = c(rep(1, nb.clust - 1),
                                                rep(18, length(vec) - nb.clust + 1)), main = title, 
            xlab = "level of cutting", ylab = "Criterion")
  }
  par(mfrow=c(1,1))
#  if (type == "tree" & tree.barplot) 
#    par(def.par)
#  invisible()
  
  
}


