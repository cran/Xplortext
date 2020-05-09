#' @export
LexHCca <- function(object, cluster.CA="docs", nb.clust="click",  min=3, max=NULL, kk=Inf, 
                    consol=FALSE, iter.max=10, order=TRUE, graph=TRUE, description=TRUE, proba=0.05, 
                    nb.desc=5, size.desc=80, seed=12345,...)
  
{

if(is.character(nb.clust)) 
  nb.clust <- ifelse(nb.clust == "auto", -1,0) 
  

    
if(nb.clust>0) if(kk!=Inf) stop("For k-means nb.clust must be 0 (user selection) or -1 (automatic selection)")
# graph.scale = add top roght the barplot
  graph.scale="inertia"
  options(stringsAsFactors = FALSE)

##### Function select ############
  select <- function(Y, default.size, method, coord.centers) {
    clust <- Y[1, ncol(Y)]	
    Y <- Y[, -ncol(Y), drop = FALSE]	
    colnames(coord.centers) <- colnames(Y)
    Z <- rbind(Y, coord.centers)	     
    if (nrow(Y) == 1) {				
      distance <- data.frame(0, row.names = "")				
      colnames(distance) <- rownames(Z[1, ])				
    } else {				
      distance <- as.matrix(dist(Z, method = method))				
      distance <- distance[(nrow(Y) + 1):nrow(distance), 				
                           -((nrow(Y) + 1):ncol(distance))]		
      distance <- sort(distance[clust, ], decreasing = FALSE)			
    }
    if (length(distance) > default.size) 				
      distance <- distance[1:default.size]				
    else distance <- distance				
  }	 # Final select function
  
  
###### Function distinctivness  #################
  distinctivness <- function(Y, default.size, method, coord.centers) {	
    clust <- as.numeric(Y[1, ncol(Y)])		
    Y <- Y[, -ncol(Y), drop = FALSE]		
    colnames(coord.centers) <- colnames(Y)
    Z <- rbind(Y, coord.centers)		
      if (nrow(Y) == 1) {		
      distance <- as.matrix(dist(Z, method = method))		
      ind.car <- vector(length = 1, mode = "numeric")		
      ind.car <- min(distance[-c(1, (clust + 1)), 1])		
      names(ind.car) <- rownames(Z[1, ])		
    }		
    else {		
      distance <- as.matrix(dist(Z, method = method))		
      distance <- distance[(nrow(Y) + 1):nrow(distance), 		
                           -((nrow(Y) + 1):ncol(distance))]		
      if (nrow(distance) == 2) 		
        center.min <- distance[-clust, ]		
      else center.min <- apply(distance[-clust, ], 2, min)		
      ind.car <- sort(center.min, decreasing = TRUE)		
    }		
    if (length(ind.car) > default.size) 		
      ind.car = ind.car[1:default.size]		
    else ind.car = ind.car		
  } # Final distinctivness



################  Function auto.cut.tree  ###########################
  auto.cut.tree <- function(res, min, max, metric, method, weight = NULL, cla = NULL, ...)				
  {				
    if (order) {				
      sss = cbind.data.frame(res$ind$coord, res$call$X, res$call$row.w, res$call$row.w.init)				
      weight <- weight[order(sss[, 1], decreasing = FALSE)]				
      sss = sss[order(sss[, 1], decreasing = FALSE), ]				
      res$ind$coord = sss[, 1:ncol(res$ind$coord), drop = FALSE]	
      res$ind$contrib = sss[, 1:ncol(res$ind$contrib), drop = FALSE]	
      res$ind$cos2 = sss[, 1:ncol(res$ind$cos2), drop = FALSE]	      
      # Las contribuciones y los cos2 no estaban ordenadas				
      res$call$X = sss[, (ncol(res$ind$coord) + 1):(ncol(sss) -  2)]				
      res$call$row.w = sss[, ncol(sss) - 1]				
      res$call$row.w.init = sss[, ncol(sss)]				
      names(res$call$row.w.init) <- names(res$call$row.w) <-  rownames(res$call$X)				
    }				
    set.seed(seed)				
    X <- as.data.frame(res$ind$coord)				
    do <- dist(X, method = metric)^2	# Matriz de distancias			
    eff <- outer(weight, weight, FUN = function(x, y, n) {x * y/n/(x + y)}, n = sum(weight))				
    dissi <- do * eff[lower.tri(eff)]				
    hc <- flashClust::hclust(dissi, method = method, members = weight)				
    inert.gain <- rev(hc$height)				
    intra <- rev(cumsum(rev(inert.gain)))				
    quot <- intra[min:(max)]/intra[(min - 1):(max - 1)]				
    nb.clust <- which.min(quot) + min - 1				
    rm(.Random.seed,envir=globalenv()) # Restaura set.seed	
    return(list(res = res, tree = hc, nb.clust = nb.clust, 				
                within = intra, inert.gain = inert.gain, quot = quot))				
  }				# Final auto.cut.tree
#################################################################################
################ Consolidation functions ======				
    consolidationW <- function(X, clust,  weights =NULL, iter.max = iter.max, ...) {			# iter.max = 10,	
    iter.max <- iter.max+1				
    cla <- cbind(X*weights,weights)	
    centers1 <- aggregate(cla, by=list(cla.cluster=clust), FUN=sum)		
    w_centers <- sweep(centers1[, c(2:(ncol(centers1)-1))], 1, centers1$weight, FUN = "/")				
    w_centers.before <- w_centers		
    centers.weight.before <- centers1$weight
      
    iter <- cluster_i <- 0				
    cluster_iminus1 <- 1				
    while(identical(cluster_i, cluster_iminus1) == F && iter < iter.max){				
      # update iteration  				
      iter <- iter + 1				
      suppressWarnings(				
        cluster_kmeans <- kmeans(x = X, centers = w_centers, iter.max = 1))#$cluster)				
      centers1 <- aggregate(cla, by=list(cla.cluster=cluster_kmeans$cluster), FUN=sum)				
      w_centers <- sweep(centers1[, c(2:(ncol(centers1)-1))], 1, centers1$weight, FUN = "/")				
      # update cluster_i and cluster_iminus1				
      if(iter == 1) {cluster_iminus1 <- 0} else{cluster_iminus1 <- cluster_i}				
      cluster_i <- cluster_kmeans$cluster				
    } # Final while				
    # X$cluster <- cluster_i				
    return(list(cluster=cluster_i,centers.before=w_centers.before, centers.after=w_centers,
                centers.weight=centers1$weights, centers.weight.before=centers.weight.before))				
  } # end  consolidationW
#######################################################
############### descfreq_New function #############################################
  descfreq_New <- function (donnee, by.quali = NULL, proba = 0.05, row.INIT, col.INIT) 				
  {
    lab.sauv <- lab <- colnames(donnee)				
    for (i in 1:length(lab)) lab[i] = gsub(" ", ".", lab[i])				
    if (!is.null(by.quali)) {				
      donnee <- as.data.frame(matrix(unlist(by(donnee, by.quali, 				
                                               apply, 2, sum)), ncol = ncol(donnee), byrow = TRUE))				
      rownames(donnee) <- levels(by.quali)
      marge.li <- aggregate(row.INIT, by=list(by.quali), FUN="sum")  # dataframe 4 x 2 columnas				
      marge.li <- as.vector(marge.li$x)  # vector 4 observaciones grupos con freq iniciales				
    }
    colnames(donnee) = lab				
    
    old.warn = options("warn")				
    options(warn = -1)				
    # if (is.null(by.quali))  marge.li = apply(donnee, 1, sum)				
    nom = tri = structure(vector(mode = "list", length = nrow(donnee)), 				
                          names = rownames(donnee))				
    marge.col=col.INIT # marge.col = apply(donnee, 2, sum)				
    
    SumTot=sum(marge.li)				
    nr <- nrow(donnee)  				
    nc <- ncol(donnee) 				

    for (j in 1:nr) {				
      aux3 = marge.li[j]/SumTot				
      for (k in 1:nc) {				
        aux2 = donnee[j, k]/marge.col[k]				
        
        if (aux2 > aux3) # sobrerepresentada				
          aux4 = phyper(donnee[j, k] - 1, marge.col[k], 				
                        SumTot - marge.col[k], marge.li[j],lower.tail = FALSE) * 2				
        else aux4 = phyper(donnee[j, k], marge.col[k], SumTot - marge.col[k], marge.li[j]) * 2				
        
        if (aux4 > 1) 				
          aux4 <- 2 - aux4				
        if (aux4 < proba) {				
          aux5 = (1 - 2 * as.integer(aux2 > aux3)) * qnorm(aux4/2)				
          aux1 = donnee[j, k]/marge.li[j]				
          tri[[j]] = rbind(tri[[j]], c(aux1 * 100, sum(marge.col[k])/SumTot * 				
                                         100, donnee[j, k], marge.col[k], aux4, aux5))				
          nom[[j]] = rbind(nom[[j]], c(colnames(donnee)[k], colnames(donnee)))				
        }				
      }				
    }				
    
    for (j in 1:nrow(donnee)) {				
      if (!is.null(tri[[j]])) {				
        oo = rev(order(tri[[j]][, 6]))				
        tri[[j]] = tri[[j]][oo, ]				
        nom[[j]] = nom[[j]][oo, ]				
        if (nrow(matrix(tri[[j]], ncol = 6)) > 1) 				
          rownames(tri[[j]]) = nom[[j]][, 1]				
        else {				
          tri[[j]] = matrix(tri[[j]], ncol = 6)				
          rownames(tri[[j]]) = nom[[j]][1]				
        }				
        colnames(tri[[j]]) = c("Intern %", "glob %", "Intern freq", 				
                               "Glob freq ", "p.value", "v.test")				
      }				
      
      
      
    }				
   names(tri) <-paste("cluster",1:length(tri),sep="_")
    res = tri				
    options(old.warn)				
    class(res) <- c("descfreq", "list ")				
    return(res)				
  }		  # End   descfreq_New
######################################################
############  create_dendrog function #######################   
  create_dendrog <- function(X, weights, clust, centers, centers.weight) {
    X.consol <- cbind(X,weights, clust) 
    
    only.dissi <- function(X, n=NULL) {				
      weight <- X[,(ncol(X))]				
      if(is.null(n)) n <- sum(weight)				
      X <- X[,c(1:(ncol(X)-1))] 				
      do <- dist(X, method = metric)^2	# Matriz de distancias			
      eff <- outer(weight, weight, FUN = function(x, y, n) {x * y/n/(x + y)}, n=n) 				
      dissi <- do * eff[lower.tri(eff)]
      return(dissi)
    }

      
    flash.clust<- function(X, n=NULL) {				
      weight <- X[,(ncol(X))]				
      if(is.null(n)) n <- sum(weight)				
      X <- X[,c(1:(ncol(X)-1))] 				
      do <- dist(X, method = metric)^2	# Matriz de distancias			
      eff <- outer(weight, weight, FUN = function(x, y, n) {x * y/n/(x + y)}, n=n) 				
      dissi <- do * eff[lower.tri(eff)]				
      hc <- flashClust::hclust(dissi, method = method, members = weight)				
      #  hc <- hclust(dissi, method = method, members = weight)				
      #   rdo <- list(dissi=dissi,eff=eff, do=do, hc=hc)			
      #   return(rdo)
      return(hc)
    }  # End flash.clust
    

    n <- sum(weights) 
    X.tot <-  cbind(X,weights)
    cl.tot <- flash.clust( X.tot, n)
    # plot(cl.tot$hc)
   # plot(cl.tot)
    nb.clust <- length(levels(X.consol$clust))
    
    X.centers <- cbind(centers,centers.weight)
    cl.centers <- flash.clust(X.centers, n)
    # plot(cl.centers$hc)
    # plot(cl.centers)
    
    
    # Intentamos determinar a partir de los centros y las ponderaciones las agrupaciones
    clustList <- hcAc <- vector('list', length(nb.clust))
    for(i in 1:nb.clust) {
      clustList[[i]] <- X.consol[X.consol$clust==i, -ncol(X.consol)] 
      if(nrow(clustList[[i]])>1) 
      clustList[[i]] <- flash.clust(clustList[[i]],n)   # el height esta bien
            hcAc[[i]] <- clustList[[i]] 
    }  

    
    ncolc <- ncol(centers)
    nrowc <- nrow(centers)


    for(i in 1:(nrow(centers)-1)){
     #  return(X.centers) centros de los clusteres (1.2..num.cluster x Dims) + centers.weight
      newdissi <- as.matrix(only.dissi(X.centers, n)) # Recalcula el height
      # stop(n) n es la frecuencia total
      # newdissi matriz 3x3 con las distancias entre los centros
      # Añade el maximo a la parte superior de la matriz
       newdissi[upper.tri(newdissi,diag=TRUE)] <- max(newdissi) 
      # return(newdissi) Matriz en la que la parte superior tiene el maximo de similaridad
       irow <- min(which(newdissi == min(newdissi), arr.ind = TRUE))
      # stop(irow) Fila que obtiene el minimo de la la tabla  ¿¿¿¿???? El resultado del a fila es 1
       jcol <- max(which(newdissi == min(newdissi), arr.ind = TRUE))
      # subcluster.add.height <- as.matrix(newdissi)[irow,jcol]
      newdissi <- only.dissi(X.centers[c(irow,jcol),], n) # Recalcula el height
      # newdissi es un valor con el height
   
            
    
      if(class(clustList[[irow]])=="hclust" & class(clustList[[irow]])=="hclust"){
        newhc <- as.hclust(merge(as.dendrogram(clustList[[irow]]),as.dendrogram(clustList[[jcol]]),
                       height=newdissi))
      } else {
          # Al menos uno de los dos contiene un cluster con solo 1 elemento
        if(class(clustList[[irow]])=="data.frame" & class(clustList[[jcol]])=="data.frame") {
          # lOS DOS CONTIENEN SOLO UN elemento
        
          newhc <- list()  # initialize empty object
          newhc$merge <-  matrix(c(-1, -2), ncol=2, byrow=TRUE)
          newhc$order <- c(1,2)
          newhc$labels <- c(rownames(clustList[[irow]]),rownames(clustList[[jcol]]) )
          newhc$height <- newdissi[1]
          newhc$method <- "ward"
          newhc$dist.method <- "euclidean"
               } else {
          i_row <- irow; j_col <- jcol
          if(class(clustList[[irow]])=="data.frame") {
            i_row <- jcol; j_col <- irow 
          }
          # UNo de ellos es cluster y el otro no, ahora el primero i_row es el cluster
          newhc <- list()  # initialize empty object
          newhc$merge <- rbind(clustList[[i_row]]$merge,
                               c(nrow(clustList[[i_row]]$merge),
                                 -(length(clustList[[i_row]]$order)+1)))
          newhc$order <- c(clustList[[i_row]]$order, 
                           (length(clustList[[i_row]]$order))+1)
          newhc$labels <- c(clustList[[i_row]]$labels, 
                            rownames(clustList[[j_col]]) )
          newhc$height <- c(clustList[[i_row]]$height,newdissi[1])
          newhc$method <- "ward"
          newhc$dist.method <- "euclidean"
        } # Final de uno de ellos es cluster y el otro no
          
        }
    class(newhc) <- "hclust" 

      # Calculo el flashClust del nuevo con los dos objetos unidos. 2x Dims + centers.weight
     clanew <- cbind(centers*centers.weight, centers.weight)[c(irow,jcol),]
     centers.weight <- c(centers.weight[-c(irow,jcol)] , sum(clanew[ncolc+1]))
     # Nuevas ponderaciones, las primeras las iniciales sin las eliminadas + ultim fila la conjunta
     
     # clanew contiene el nuevo centro de los unidos
          clanew <- apply(clanew/centers.weight[length(centers.weight)],2,FUN=sum)[1:ncolc]
          # Recalcula centros del cluster de agregacion
          centers <- rbind(centers[-c(irow,jcol),], clanew ) # Antiguos centros + nuevo
          X.centers <- cbind(centers,centers.weight)

          # Quitamos el cluster con mayor numero de los unidos de la lista
          # Quitamos el menor cluster de los unidos de la lista
          #  return(X.centers)
          clustList <- clustList[-c(irow,jcol)]
          clustList[[nrow(centers)]]  <- newhc
          
          
          
    } # End for
        clustList[[1]]$call <- NULL
    return(as.hclust(clustList[[1]]))
  }

  ###################### final functions #####################

  
# description ------------------------------------------------
# If marg.doc ==before las marginales fila (documentos activosxpalabras) son los originales, 
# con todas las palabras, hayan sido seleccionadas o no. Solamente los documentos activos son eliminados
# Si True las marginales fila (documentos activos) son el numero de palabras de los documentos 
# despues del umbral, es decir, de la tabla DocTerm
############# Falta edit.par  
  #if TRUE, the literal text of the parangon and specific documents are listed in the results 
  # (by default FALSE)
  ######## Falta description
  ### data.clust	
  ## the original active lexical table used in LexCA plus a new column called clust containing the partition
  
  # iter.max consolidacion 
  # nb.clust es el number of cluster
# If 0, the tree is cut at the level the user clicks on.
# If -1, the tree is automatically cut at the suggested level (see details).
#  If a (positive) integer, the tree is cut with nb.cluters clusters
# kk Inf y dendrograma/jerarquico (aunque consol=TRUE)
#  A boolean. If TRUE, clusters are ordered following their center coordinate on the first axis
  # graph.scale	
  # A character string. By default "inertia" and the height of the tree corresponds to the inertia gain, else "sqrt-inertia" the square root of the inertia gain.
  # nb.par	
  # An integer. The number of edited paragons.
#  proba	
#  The probability used to select axes and variables in catdes (see catdes for details
  

  
  
# Starting checks -----------------------------------------------------------------------------  
  if (!inherits(object, "LexCA")) stop("object should be LexCA class")
  if (min <=1) min<-2
  metric = "euclidean"
  method = "ward"
  # kk = Inf
  # graph.scale = "inertia"
  if (cluster.CA=="rows" ) cluster.CA <- "docs"
  if (cluster.CA=="columns" ) cluster.CA <- "words"

    if ((kk != Inf) & (consol == TRUE)) {				
    warning("No consolidation has been done after the hierarchical clustering since kk is different from Inf")				
    consol <- FALSE				
  }				
  
###############################################
# object$call$X
#  RESUMEN	DOCS	PALABRAS	TOTAL PAL
#  object$call$X	297	37	1282
# write.csv2(object$call$X, file="C:/Xplortext/X.csv")  # suma 1282
# Los número de los docs están bien
###############################################


  doc.w.after <- apply(object$call$X, 1, sum)    # Initial weighting for documents after selection
  nm <- names(doc.w.after)
  doc.w.before <- object$rowINIT[which(rownames(object$rowINIT) %in% nm) ,]
  
  
  names(doc.w.before) <- nm
  word.w <- apply(object$call$X, 2, sum) 
  
  
#   if(marg.doc=="before")      doc.w <- doc.w.before # before word's selection
  
  if (cluster.CA == "docs") {
    weight.dw.cl <- doc.w.after
    coord.dw <- object$row$coord
    
  } else {
    weight.dw.cl <- word.w
    coord.dw <- object$col$coord
  }
    res.save <- object
    cla <- NULL	
    vec <- FALSE
    
    # Operations with kk ----------------
    #  If kk = Inf => hierarchical cluster
    #  If kk is not Inf => k-means cluster. If k-means coordinates * sqrt(eigenvalue)   
       
    
    
    ################# Selecting cluster of rows (documents) or columns (words) ####################
    ## If hierarchical kk = Inf
      if (kk < nrow(coord.dw)) {
        res <- as.data.frame(sweep(coord.dw, 2, sqrt(object$eig[1:ncol(coord.dw), 1]), FUN = "*"))
        kk <- min(kk, nrow(unique(coord.dw)))
      } else kk=Inf				
    ######################################
# If kk >= nrow(word) now a hierarchical cluster is performed, res don't exist
    


# KM ############ Weighted kmeans    
    if (exists("res")) {	# It is k-means. res has a list actif documents (or words) x dimensions	
      # res <- res[, unlist(lapply(res, is.numeric)), drop = FALSE] 	
      if (kk < nrow(res)) kk <- min(kk, nrow(unique(res)))	
      if (kk <= nrow(res)) {	
        # Weighted k-means (no with Factoclass) Do it with flexclust   
        set.seed(seed)

         cla <- flexclust::cclust(res,k=kk, method="hardcl",weights= weight.dw.cl)
         cla <- data.frame(cla@cluster) # documents/words win names x 1 column whith number of cluster	
         cla <- cbind(res, weight.dw.cl ,cla[rownames(res),,drop=FALSE])
         cla1 <- cbind(cla[,c(1:(ncol(cla)-2))]*cla[,(ncol(cla)-1)],cla[,c((ncol(cla)-1):ncol(cla))])	
         cla.agg <- aggregate(cla1[,c(1:(ncol(cla1)-1))], by=list(cla.cluster=cla$cla.cluster), FUN=sum)	
         cla.centers <-  cla.agg[,c(2:(ncol(cla.agg)-1))] /cla.agg[,ncol(cla.agg)]    
         res <- PCA(cla.centers, row.w = cla.agg[,ncol(cla.agg)], scale.unit = FALSE, ncp = Inf, graph = FALSE)
         rm(.Random.seed,envir=globalenv()) # Retrieve set.seed
      }          # end      if (kk <= nrow(res)) 
    } # Final of k-means
    # res has one PCA object if k means has been computed. res is NULL if is hierarchical cluster


    if(kk==Inf){ ## It is hierarchical cluster, no k-means
        aux <- object$eig			
        res <- FactoMineR::PCA(coord.dw, scale.unit = FALSE, ncp = Inf, graph = FALSE, row.w = weight.dw.cl)	
        res$eig <- aux			
    }
############## Now, all res object are PCA objects type
#======================================================================
    
    
    
### For build the tree    
    #  Checking the maximum number of clusters. "ind" is documents if cluster.CA="docs", words if cluster.CA="words"
    if (is.null(max))
    max <- min(10, round(nrow(res$ind$coord)/2))			
    max <- min(max, nrow(res$ind$coord) - 1)  


    # Retrieve the tree with cases if hierarchical or groups if k-means	----------
    # row.w.init has the weighting from rows from res object of PCA
    t <- auto.cut.tree(res, min = min, max = max, metric = metric,method = method,
                       weight = res$call$row.w.init, cla = cla, order = order, ...)	
    # t$nb.clust is the number of auto selected clusters
    # If nb.clust is -1 (automatic) we take the number of clusters computed in t (t$nb.clust)
    
    if(nb.clust==-1) nb.clust <- t$nb.clust
  
    # inertia.height not exists 
    if (graph.scale == "inertia") {							
      nb.ind <- nrow(t$res$ind$coord)			## Son individuos (docs o words) 				
      inertia.height <- rep(0, nb.ind - 1)							
      for (i in 1:(nb.ind - 1)) inertia.height[i] <- t$inert.gain[(nb.ind -i)]							
      inertia.height <- sort(inertia.height, decreasing = FALSE)	
      t$tree$height <- inertia.height							
    }			# end of graph.scale =="inertia"
   
   
    # Checking if the selection of the number of clusters is doing cliking on the dendrogram     
    if (nb.clust == 0)  graph <- TRUE     # If is 0 the user click on the dendrogram
    if (nb.clust == 1)  #
    {
      warning("The number of clusters can not be 1. Automatic change: nb.clust=",-1)				
      nb.clust <- -1
    }  

    if(nb.clust>1)
    if(kk<=nb.clust) {
      warning("The number of clusters used in a Kmeans preprocessing before the hierarchical clustering
              is less than the number of clusters nb.clust. Automatic change: nb.clust=",kk-1)				
      nb.clust <- kk-1
      #  kk <- kk+1
    }    

    if(nb.clust>max) {
      warning("The number of selected clusters is bigger than the number of maximum clusters.
              Automatic change: max=", min(nb.clust, nrow(res$ind$coord) - 1))				
      max<- min(nb.clust,nrow(res$ind$coord) - 1)
      nb.clust <- max  # New instruction
      #  kk <- kk+1
    }    

    auto.haut <- ((t$tree$height[length(t$tree$height) - t$nb.clust + 2]) + 							
                    (t$tree$height[length(t$tree$height) - t$nb.clust + 1]))/2							
    # Break point propositon


    # Print title and Inertia gain graph							
    if (graph) {
          if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 							
      dev.new()							
      old.mar <- par()$mar							
      par(mar = c(0.5, 2, 0.75, 0))							
      lay <- matrix(ncol = 5, nrow = 5, c(2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4,
                                         2, 4, 4, 4, 4, 1, 3, 3, 3, 3))							
      layout(lay, respect = TRUE)							
      barplot(t$inert.gain[1:max(15, max)], col = c(rep("black",t$nb.clust - 1),							
                                                    rep("grey", max(max, 15) - t$nb.clust + 1)), 							
              rep(0.1, max(max, 15)), space = 0.9)							
      plot(x = 1, xlab = "", ylab = "", main = "", col = "white", axes = FALSE)							
      text(1, 1, "Hierarchical Clustering", cex = 2)							
      plot(x = 1, xlab = "", ylab = "", main = "", col = "white", axes = FALSE)							
      legend("top", "inertia gain  ", box.lty = NULL, cex = 1)			
    } # End graph
    else {			
      if (nb.clust == 0 | nb.clust == 1) 			
        nb.clust <- -1			
    }			

    if (nb.clust == 0| (nb.clust == 1))  {							
      print("Click on the graph to cut the tree")							
      flush.console()							
      plot(t$tree, hang = -1, main = "Click to cut the tree", xlab = "", sub = "")							
      abline(h = auto.haut, col = "black", lwd = 3)							
      coupe <- locator(n = 1)							
      while (coupe$y < min(t$tree$height)) {							
        cat("No class \n")							
        coupe <- locator(n = 1)							
      }						
      # Fin de while	
      y <- coupe$y							
      # y is a value of height, value of xy axis							
      # If it is not manual selections inertia gain graph is not plotted
    } # Final of nb.clust==0  
    else {
          if (nb.clust < 0) 	
        y = auto.haut	
      else y = (t$tree$height[length(t$tree$height) - nb.clust + 2] + 
                  t$tree$height[length(t$tree$height) - nb.clust + 1])/2	
    }
  
    if(graph){
    layout(lay, respect = TRUE)							
    barplot(t$inert.gain[1:max(15, max)], col = c(rep("black",t$nb.clust - 1),							
                                                  rep("grey", max(max, 15) - t$nb.clust + 1)), 							
            rep(0.1, max(max, 15)), space = 0.9)							
    plot(x = 1, xlab = "", ylab = "", main = "", col = "white", axes = FALSE)							
    text(1, 1, "Hierarchical Clustering", cex = 2)							
    plot(x = 1, xlab = "", ylab = "", main = "", col = "white", axes = FALSE)							
    legend("top", "inertia gain  ", box.lty = NULL, cex = 1)							
    plot(t$tree, hang = -1, main="", xlab = "", sub = "")							
    } # Enf if(graph)
    
    clust <- cutree(as.hclust(t$tree), h = y)	
    
    # integer with the number of the cluster
    # The name is the label of the doc/word 
    nb.clust <- max(clust)	
    # Order of clusters in dendrogram (duplicates are removed)
    ordColo <- unique(clust[t$tree$order])	

    
    X <- as.data.frame(t$res$ind$coord)	   # Docs or words x dimensions

    
# Add rectangles to graph					
    if (graph) {							
      rect <- rect.hclust(t$tree, h = y, border = ordColo)							
      clust <- NULL							
      for (j in 1:nb.clust) clust <- c(clust, rep(j, length(rect[[j]])))							
      clust <- as.factor(clust)							
      belong <- cbind.data.frame(t$tree$order, clust)							
      belong <- belong[do.call("order", belong), ]							
      clust <- as.factor(belong$clust)							
      names(clust) <- rownames(X)							
      if (nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 							
        layout(matrix(nrow = 1, ncol = 1, 1), respect = TRUE)							
    }							

    
   if (kk != Inf)  consol <- FALSE		# Consolidation only for hierarchical, not kmeans		
    
    # If it is k means all ss are computed later, here only for hierarchical					
    if(kk==Inf) {	
      ss.between <- sum(rev(t$tree$height)[1:(nb.clust -1)])
      ss.tot <- t$within[1] 
      ss.within <- t$within[nb.clust]		
    }	
    
    

# if(consol==TRUE) t_tree_before <- t$tree    # Only for hierarchical cluster
    if (consol==TRUE) {		
    t_tree_before <- t$tree    # Only for hierarchical cluster
    # Only if not kmeans  # X cases by dimensions			
    # Order of X is the same as dendrogram ; res$call$row.w.init is ordered as initial docs or alphab. if words
      weights_init <- res$call$row.w.init[rownames(X)]   # same order than dendrogram and clust variable
      clust.afterconsol <-  consolidationW(X, clust=clust,  weights=weights_init, iter.max = iter.max, ...) 
      clust.after <- as.factor(clust.afterconsol$cluster)		
      clust.before <- clust	
      centers.after <- clust.afterconsol$centers.after					# Coordinates		
      centers.before <- clust.afterconsol$centers.before	
      centers.weight <- clust.afterconsol$centers.weight	      # Weights after
      
      centr.dim  <- colSums(centers.after * centers.weight) / sum(centers.weight)
      ss.between.afterconsol <- sum((sweep(centers.after,2,centr.dim ,FUN="-"))^2 * 
                                      centers.weight)/sum(centers.weight)
      
      if(ss.between.afterconsol < ss.between) { # Check if the new solutions is worst					
        clust <- clust.before							
        centers <- centers.before							
        centers.weight <- clust.afterconsol$centers.weight.before							
        ss.between.afterconsol <- ss.between							
      } else {							
        centers<-  centers.after							
        clust <- clust.after							
      } # end if(ss.between.afterconsol < ss.between)							
  
      if(order){							
        aux <- names(clust) # Gz86...Zp08							
        ord <- order(centers[, 1, drop = FALSE])							
        centers <- centers[ord, , drop = FALSE]							
        centers.before <- centers.before[ord, , drop = FALSE]	
        centers.weight <- centers.weight[ord]							
        clust <- (order(ord))[clust]							
        names(clust) <- aux							
        X <- X[aux,]							
      } else {							
        clust <- clust[rownames(X)]							
      }  # Final de if(order)	
    } # Final de if(consol)  



    
    if (!consol | kk==Inf) {	
        weights_init <- res$call$row.w.init[rownames(X)]			
        cla2 <- cbind(X*weights_init,weights_init)		
        centers1 <- aggregate(cla2, by=list(cla2.cluster=clust), FUN=sum)							
        centers.weight <- centers1$weight							
        centers <- sweep(centers1[, c(2:(ncol(centers1)-1))], 1, centers1$weight, FUN = "/")							
        if(order) {							
          #aux <- names(clust) <- rownames(X)							
          aux <- names(clust)	# Gz86Gz89CS81Gz82Gz93Su79Rj11Zp04Az96Zp08				
          ord <- order(centers[, 1, drop = FALSE])	# 43251
          centers <- centers[ord, , drop = FALSE]	
          clust <- (order(ord))[clust]		
          centers.weight <- centers.weight[ord,drop = FALSE]	
          names(clust) <- aux	
          clust <- as.factor(clust)
        }		# End If order	
        
        
        
        if(kk<Inf) {
          clust <- clust[order(as.integer(names(clust)))]	
          rr <- clust[cla$cla.cluster]
          names(rr) <- rownames(cla)	
          clust<- cla$cla.cluster <- rr	
        } # End kk<Inf
    } # End if !consol
    clust <- as.factor(clust)

    
    
    # If consolidation = TRUE dendrogram must be re calculated again							
    if(consol) {							
      dendroConsol <- create_dendrog(X,weights_init, clust, centers, centers.weight)
      #t$tree <- as.hclust(dendroConsol)	
      t$tree <- dendroConsol
      t$tree$method <- "ward"							
      t$tree$dist.method <- "euclidean"							
      t$inert.gain <- rev(t$tree$height)							
      t$within <- rev(cumsum(t$tree$height))							
      t$quot <- t$within[min:(max)]/t$within[(min - 1):(max - 1)]	
      if(graph==TRUE) {							
        if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 							
          dev.new()							
        plot(t$tree, hang=-1, main="Reconstructed dendrogram after consolidation",							
             xlab = "", sub = "", cex=0.7, cex.main = 0.8, cex.axis=0.7)							
        rect.hclust(t$tree, k=nb.clust,  border = ordColo)							
      }							
    } # end if(consol) 							

 if(exists("t_tree_before")) t$tree$call <- t_tree_before$call
    
  
    


# Computing statistical results for clusters    
# res.save$call$X contiene la tabla lexica/agregada  docs x words
# dfdw has the frequencies of active docs(or words) + number of cluster)  (antes dfword)  	
    if(cluster.CA=="docs") {	
      dfdw <- data.frame(cbind(res.save$call$X,clust_=clust[rownames(res.save$call$X)]))	
      colnames(dfdw) <- c(colnames(res.save$call$X),"clust_")
    #  rownames(dfdw) <- rownames(res.save$call$X)
      coord.clust <- res.save$row$coord # data.clust <- res.save$row$coord	
    } else {					
      dfdw <- data.frame(cbind(t(res.save$call$X),clust_=clust[colnames(res.save$call$X)]))	
      colnames(dfdw) <- c(rownames(res.save$call$X),"clust_")
#      rownames(dfdw) <- colnames(res.save$call$X)
#      colnames(dfdw)[1:(ncol(dfdw)-1)] <- rownames(res.save$call$X)
      coord.clust <- res.save$col$coord  #  data.clust <- res.save$col$coord							
    }

  
    # coord.clust Coordenadas docs x (dim + cluster + ponderacion) si docs							
    coord.clust <- data.frame(cbind(coord.clust, weight= weight.dw.cl[rownames(coord.clust)]),							
                              clust=clust[rownames(coord.clust)])
    
        # Computing sum of squares (sctot) for each group		
    sum.weights <- sum(weight.dw.cl)							
    dimen <- ncol(coord.clust)-2				
    cla2 <- cbind(coord.clust[,c(1:dimen)]*coord.clust$weight, 							
                  coord.clust[,c((dimen+1):ncol(coord.clust))])	

  
    mediadim <- apply(cla2[,c(1:(ncol(cla2)-2))],2,FUN=sum)/sum.weights # Near 0							
    sctot <- apply((coord.clust[,c(1:dimen)] - mediadim)^2*(coord.clust$weight/sum.weights),1,FUN="sum")	
    sc.tot <- sum(sctot)  # <---------- pendiente de utilizacion posterior
  
    # sum of squares for each cluster
    sctot <- aggregate(sctot, by=list(cla.cluster=coord.clust$clust), FUN=sum)$x
    # sctot <- sctot$x   # suma de cuadrados totales de cada cluster	
    # sctot is the sum of squares for each document (word)							
    # Computing centers		
  
    
        
    centers1 <- aggregate(cla2[,c(1:(dimen+1))], by=list(cla.cluster=coord.clust[,(dimen+2)]), FUN=sum)	
    w_centers <- sweep(centers1[, c(2:(dimen+1))], 1, centers1$weight, FUN = "/") 		
    

    scinter <- apply((w_centers - mediadim)^2*(centers1$weight/sum.weights),1,FUN="sum")							
    sc.inter <- sum(scinter)


  
    # Computing sum squares for each cluster    							
    # dfclust <- as.data.frame(table(dfword$clust_))
    # Number of elements of each cluster
    dfclust <- as.data.frame(table(dfdw$clust_))
    
    dfclust$weights <- centers1$weight
    dfclust$sstot <- sctot							
    dfclust$ssinter <- scinter							
    dfclust$ssintra <- dfclust$sstot-dfclust$ssinter 							
    dfclust$Explained <- dfclust$ssinter*100/ dfclust$sstot							
    
    if (cluster.CA == "docs") {							
      colnames(dfclust) <- c("Cluster", "Docs", "Occurrences", "ss.tot", "ss.inter", "ss.intra", "%Explained")							
    } else {							
      colnames(dfclust) <- c("Cluster", "Words", "Occurrences", "ss.tot", "ss.inter", "ss.intra", "%Explained")							
    }							

    compo.cluster<-list()							
    for (iclus in 1:nb.clust)  compo.cluster[[iclus]]<- rownames(dfdw)[dfdw$clust_ == iclus]							
    
    
         desc.ind <- descaxes <- descword <- descwordsup <- descdoc <- descquali <- descquanti <- NULL

    if(description) {
      if(cluster.CA=="docs") {
        # Cluster Description by active words
        descword <- descfreq_New(dfdw[,c(1:length(word.w))],dfdw[,ncol(dfdw)],prob=proba, 
                                 row.INIT=doc.w.before[rownames(dfdw)],col.INIT= word.w) 
      } else { 
        descdoc <- descfreq_New(dfdw[,c(1:length(doc.w.before))],dfdw[,ncol(dfdw)],prob=proba,
                                 row.INIT=word.w[rownames(dfdw)],
                                  col.INIT=doc.w.after[colnames(dfdw)])
      } 
    ######################## FINAL descword, Cluster description by words and descdoc

      
      
    
  ######################## Palabras suplementarias    
      if(!is.null(res.save$col.sup$coord)) {
        # Hay palabras suplementarias
        #descwordsup <-  data.clust[,c(rownames(res.save$col.sup$coord), "clust_")]
        # stop(rownames(res.save$col.sup$coord))
        if(cluster.CA=="docs") {
          descwordsup <- data.frame(cbind(res.save$call$Xtot[rownames(res.save$call$X), 
                                                             rownames(res.save$col.sup$coord)],
                                          clust_=clust[rownames(res.save$call$X)])  )
          colnames(descwordsup) <- c(rownames(res.save$col.sup$coord), "clust_")
          WSUP <- apply(descwordsup[,c(1:(ncol(descwordsup)-1))],2,"sum")
          descwordsup <- descfreq_New(descwordsup[,c(1:(ncol(descwordsup)-1))],descwordsup[,ncol(descwordsup)],
                                      prob=proba, row.INIT=doc.w.before[rownames(descwordsup)], 
                                      col.INIT=WSUP)
        }}
      ######################## FINAL descwordsup, descripcion de los clusteres de documentos por palabras suplement.
  

      ############ Descripcion de documentos por cualitativas suplementarias         
      if(!is.null(res.save$quali.sup$coord)) {
        if(cluster.CA=="docs") {
          a1 <- res.save$call$Xtot[rownames(res.save$call$X), rownames(res.save$quali.sup$eta2), drop=FALSE]
          a2 <- cbind(a1, clust_=clust[rownames(a1)])
          descquali <- FactoMineR::catdes(a2, num.var=ncol(a2), proba = proba, row.w=doc.w.before[rownames(a2)])
          names(descquali$category) <-paste("cluster",1:length(descquali$category),sep="_")
        }}
      ############ Final cualitativas suplementarias         
      
      
      

    ############# Inicio cuantitativas suplementarias
    if(!is.null(res.save$quanti.sup$coord)) {
      if(cluster.CA=="docs") {
        a1 <- res.save$call$Xtot[rownames(res.save$call$X), rownames(res.save$quanti.sup$coord),drop=FALSE]
        a2 <- cbind(a1, clust_=clust[rownames(a1)])
 #       descquanti <- catdes.New(a2, num.var=ncol(a2), proba = proba, row.w=doc.w.before[rownames(a2)])
#     return(descquanti)
        
          descquanti <- FactoMineR::catdes(a2, num.var=ncol(a2), proba = proba, row.w=doc.w.before[rownames(a2)])

          
      if(!is.null(descquanti$quanti)>0)  
        names(descquanti$quanti) <-paste("cluster",1:length(descquanti$quanti),sep="_")
      }}
    ############# Final cuantitativas suplementarias

      
 
      
      
      
    ############## descripcion de los ejes    
    if(cluster.CA=="docs") {
      descaxes <- data.frame(cbind(res.save$row$coord,clust_=as.factor(clust[rownames(res.save$row$coord)])))
      descaxes$clust_ <- as.factor(descaxes$clust_)
      descaxes <- FactoMineR::catdes(descaxes, num.var=ncol(descaxes), proba = proba, row.w=doc.w.before[rownames(descaxes)])  
    } else { # columnas
      descaxes <- data.frame(cbind(res.save$col$coord,clust_=as.factor(clust[rownames(res.save$col$coord)])))
      descaxes$clust_ <- as.factor(descaxes$clust_)
      descaxes <- FactoMineR::catdes(descaxes, num.var=ncol(descaxes), proba = proba, row.w=word.w[rownames(descaxes)]) 
    }
    } # Final descriptions

         
    # tabInd coordinates + cluster if active docs/words
    if (cluster.CA == "docs"){ 
      tabInd <- cbind.data.frame(res.save$row$coord, clust=dfdw[rownames(res.save$row$coord),ncol(dfdw)])
    }else{ 
      tabInd <- cbind.data.frame(res.save$col$coord, clust=dfdw[rownames(res.save$col$coord),ncol(dfdw)])
    }
 

    if(description) {  
      if(nb.desc>0){
        # Calculo de los centros
        # list.centers <- by(tabInd[, -ncol(tabInd), drop = FALSE], tabInd[, ncol(tabInd)], colMeans)
        # centers <- matrix(unlist(list.centers), ncol = ncol(tabInd) - 1, byrow = TRUE)		
        m_centers <- data.matrix(w_centers)

        Cluster <- tabInd[,ncol(tabInd), drop=FALSE]
        para <- by(tabInd, Cluster, simplify = FALSE, select, default.size = nb.desc, 
                   method = metric, coord.centers = m_centers)
        
        
        dist <- by(tabInd, Cluster, simplify = FALSE, distinctivness, 
                   default.size = nb.desc, method = metric, coord.centers = m_centers)
        desc.ind <- list(para = para, dist = dist)
        } # End nb.par
    } # End description
   

         
         
    ss <- list(info=dfclust, ss.tot=sum(dfclust$ss.tot), ss.between=sum(dfclust$ss.inter),
              ss.intra=sum(dfclust$ss.intra), inertia.ratio=sum(dfclust$ss.inter)*100/sum(dfclust$ss.tot))  
    Xnew <- coord.clust[,-(ncol(coord.clust)-1)]
    
    newcall <- list(t=t, min=min, max=max, X=Xnew, vec=FALSE, call = match.call())							

    newres <- list(data.clust=dfdw, centers=w_centers, clust.count=ss, clust.content=compo.cluster, 
                      call=newcall)	
   
     
        if(!is.null(descaxes)) {
          class(descaxes) <- NULL
          newres$description$desc.axes <- descaxes
        }
  
  
      if(cluster.CA =="docs") {
        if(!is.null(descword)) 
          newres$description$desc.cluster.doc <- list(words=descword, wordsup=descwordsup, qualisup=descquali, quantisup=descquanti)
      #   if(!is.null(desc.ind))  
      #   newres$docslabels <- desc.ind  # antes $desc.doc
      } else{
        if(!is.null(descdoc)) 
          newres$description$desc.cluster.word <- list(docs=descdoc)
      #  if(!is.null(desc.ind))  
      #   newres$wordslabels <- desc.ind   # antes $desc.word
      }
  
    

    if(nb.desc>0){ 
      #### description of the clusters by the characteristic documents
     
      if(is.null(desc.ind$para)) stop("Please, use description=TRUE \n  to obtain a description of the clusters by the characteristic documents ")
      
      #### creation of the structure of corpus
      var.text <- object$info$var.text[[1]]
      str.base <- object$info$base[[1]]
      str.envir <- object$info$menvir[[1]]
      base <- get(str.base, envir=str.envir)
      corpus <- base[, var.text[1]]
      if(length(var.text) > 1) {					
        for (i in 2:length(var.text)){					
          corpus <- paste(corpus, base[, var.text[i]], sep = ".")
        }}
      corpus <- data.frame(corpus, stringsAsFactors = FALSE)
      rownames(corpus) <- rownames(base)
      #### corpus[rownames(object$call$X),] #########?
      lispara <- vector(mode="list")

      for(iclus in 1:nb.clust)
      {
        doctot <- length(desc.ind$para[[iclus]])
        ntdoc <- min(nb.desc, doctot)
        lispara[[iclus]]<-data.frame()

      for (i in 1:ntdoc)
        {
          lispara[[iclus]][i,1] <- names(desc.ind$para[[iclus]][i])
          lispara[[iclus]][i,2] <- desc.ind$para[[iclus]][i]
          if(cluster.CA=="docs") lispara[[iclus]][i,3] <- strtrim(corpus[lispara[[iclus]][i,1],1], size.desc)
        } # End for (i in 1:ntdoc)
        if(cluster.CA=="docs") colnames(lispara[[iclus]]) <- c("DOCUMENT", "CRITERION", "TEXT")
         else colnames(lispara[[iclus]])  <- c("WORD", "CRITERION")
       } # Final for iclus
        names(lispara) <-paste("cluster",1:nb.clust,sep="_")
 #     newres$docspara <- lispara

      if(cluster.CA=="docs") newres$description$desc.cluster.doc$para <- lispara else newres$description$desc.cluster.word$para <- lispara
      lisdist <- vector(mode="list")
      for(iclus in 1:nb.clust)
      {
        doctot <- length(desc.ind$dist[[iclus]])
        ntdoc <- min(nb.desc, doctot)
        lisdist[[iclus]]<-data.frame()
         for (i in 1:ntdoc)
        {
          lisdist[[iclus]][i,1] <- names(desc.ind$dist[[iclus]][i])
          lisdist[[iclus]][i,2] <- desc.ind$dist[[iclus]][i]
          if(cluster.CA=="docs") lisdist[[iclus]][i,3] <- strtrim(corpus[lisdist[[iclus]][i,1],1], size.desc)
        } # End for (i in 1:ntdoc)
        if(cluster.CA=="docs") colnames(lisdist[[iclus]]) <- c("DOCUMENT", "CRITERION", "TEXT")
        else colnames(lisdist[[iclus]])  <- c("WORD", "CRITERION")
      } # Final for iclus
      names(lisdist) <-paste("cluster",1:nb.clust,sep="_")
 #    newres$docsdist <- lisdist
      if(cluster.CA=="docs") newres$description$desc.cluster.doc$dist <- lisdist else newres$description$desc.cluster.word$dist <- lisdist
       } # End nb.par >0 
      
   class(newres) <- c("LexHCca")  # , "HCPC")			
    return(newres)
   } # Final correspondiente a la funcion