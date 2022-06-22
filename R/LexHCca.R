#' @import flexclust
#' @export
LexHCca <- function(x, cluster.CA="docs", nb.clust="click",  min=2, max=NULL, kk=Inf, 
                    consol=FALSE, iter.max=500, graph=TRUE, description=TRUE, proba=0.05, 
                    nb.desc=5, size.desc=80, seed=12345,...)
  
{
  dots <- list(...)
  if(is.null(dots$method)) method.sel <- "ward"
    else method.sel <-  dots$method
    

    object <- x
if (!inherits(object, "LexCA")) stop("object should be LexCA class")
  options(stringsAsFactors = FALSE)
  
  marg.doc<-"before"
 if(description==FALSE)  nb.desc <-0
  
#---- Initial checks -----------------------------------------
  ## Type selection number of clusters
if(is.character(nb.clust))
  nb.clust <- ifelse(nb.clust == "auto", -1,0)
 # if(nb.clust=="click") nb.clust<- 0
 #  if(nb.clust=="auto") nb.clust<- -1
  
  
# if(nb.clust>0) if(kk!=Inf) stop("For k-means nb.clust must be 0 (user selection) or -1 (automatic selection)")
# graph.scale = add top right the barplot
  graph.scale="inertia"

  if ((kk != Inf) & (consol == TRUE)) {
    warning("No consolidation has been done after the hierarchical clustering since kk is different from Inf")				
    consol <- FALSE
  }	
  ### nb.clust can be a number; 0 is click; -1 is automatic
  
  ##
  if (min <2) min<-2
  metric = "euclidean"
  if (cluster.CA=="rows" ) cluster.CA <- "docs"
  if (cluster.CA=="columns" ) cluster.CA <- "words"
 #---- End Initial checks -----------------------------------------


  
#---- F1. auto.cut.tree function -----
  auto.cut.tree <- function(res, min, max, metric, method=method.sel, weight = NULL, cla = NULL, ...)
  {
    set.seed(seed)
    X <- as.data.frame(res$ind$coord)
    do <- dist(X, method = metric)^2			
    eff <- outer(weight, weight, FUN = function(x, y, n) {x * y/n/(x + y)}, n = sum(weight))
    dissi <- do * eff[lower.tri(eff)]
    hc <- flashClust::hclust(dissi, method = method, members = weight)
    inert.gain <- rev(hc$height)
    intra <- rev(cumsum(rev(inert.gain)))
    quot <- intra[min:(max)]/intra[(min-1):(max-1)]
    nb.clust <- which.min(quot) + min-1
    rm(.Random.seed,envir=globalenv()) # Retrieve set.seed
    return(list(res = res, tree = hc, nb.clust = nb.clust,
                within = intra, inert.gain = inert.gain, quot = quot))
  }	  # End Function auto.cut.tree  
  
  

  
#---- F2. print.inertia function ----
  print.inertia <- function(title) {
    if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
      dev.new()
    old.mar <- par()$mar
    par(mar = c(0.5, 2, 0.75, 0))
    lay = matrix(ncol = 5, nrow = 5, c(2, 4, 4, 4, 4, 
                                       2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4,1, 3, 3, 3, 3))
    layout(lay, respect = TRUE)
    barplot(rev(t$tree$height)[1:max(15, max)], col = c(rep("black", 
                                                            t$nb.clust - 1), rep("grey", max(max, 15) - t$nb.clust + 1)), rep(0.1, max(max, 15)), space = 0.9)
    plot(x = 1, xlab = "", ylab = "", main = "", col = "white",axes = FALSE)
    text(1,1,title,cex=2)
    plot(x = 1, xlab = "", ylab = "", main = "", col = "white",axes = FALSE)
    legend("top", "Aggreg. Crit", box.lty = NULL, cex = 1)
  }    # End function print.inertia
  
  


#---- F3. sum of squares ----
  ssquares <- function(coord, weight_,clust_){ # }, weight_, clust_) {
    dfw <- data.frame(coord,"weight_"=weight_, "clust_"=clust_)
        # write.csv2(dfw,file="C:/Xplortext/dfw.csv")
        sum.weights <- sum(dfw$weight_)
        pct.weight <- dfw$weight_/ sum.weights  
        dfclust <- as.data.frame(table(clust_))
        dimen <- ncol(coord)
    mediadim <- apply(dfw[,c(1:dimen)]*pct.weight,2,FUN=sum)
    dfw$sstot_ <- apply(sweep(dfw[,c(1:dimen)],2,mediadim,"-")^2  *pct.weight,1,FUN="sum")
    cl.weight  <- aggregate(pct.weight, by=list(dfw$clust_), FUN=sum)$x
    coord_centers <- aggregate(dfw[,c(1:dimen)]*pct.weight/cl.weight[dfw$clust_],by=list(dfw$clust_), FUN=sum) 
    dfw$ssintra_ <- apply(((dfw[,c(1:dimen)] - coord_centers[clust_, c(2:ncol(coord_centers))])^2),1,sum) *pct.weight
    dfw$ssinter_ <- dfw$sstot_-dfw$ssintra_
    dfw$ssexplain_ <- 100*dfw$ssinter_ / dfw$sstot_ 
    dfclust$weight <- cl.weight*sum.weights 
    dfclust$sstot <- aggregate( dfw$sstot_,by=list(dfw$clust_), FUN=sum)$x
    dfclust$ssintra <- aggregate( dfw$ssintra_,by=list(dfw$clust_), FUN=sum)$x
    dfclust$ssinter <- dfclust$sstot-dfclust$ssintra
    dfclust$ssexplain <- 100*dfclust$ssinter / dfclust$sstot 
    ss.tot <- sum(dfclust$sstot)  
    ss.intra <- sum(dfclust$ssintra)
    ss.inter <- sum(dfclust$ssinter)
    pct.explained <- 100* ss.inter / ss.tot

       if (cluster.CA == "docs") {							
      colnames(dfclust) <- c("Cluster", "Docs", "Occurrences", "ss.tot", "ss.intra", "ss.between", "%Explained")							
    } else {							
      colnames(dfclust) <- c("Cluster", "Words", "Occurrences", "ss.tot", "ss.intra", "ss.between", "%Explained")							
    }	
     dfw <- dfw[,c((ncol(dfw)-5):ncol(dfw))]
     colnames(dfw) <- c("Weight","Cluster", "ss.tot", "ss.intra", "ss.between","%Explained")
    
    res.ss <- list("ss.tot"=ss.tot, "ss.intra"=ss.intra,"ss.inter"= ss.inter, "coord.centers"=coord_centers, "dfw"=dfw,
                   "dfclust"=dfclust)
    return(res.ss)
  }  # End function F3. sum of squares
  
  
  
  #---- F4. Consolidation function ----
  consolidationW <- function(coordinates, cluster_i,  weights =NULL, method="hardcl", iter.max = iter.max, ...) {		
    # X son las coordenadas documents x dims
    # iter.max <- iter.max+1
     # init_centers <- as.matrix(centers[, c(2:ncol(centers))])
    iter <- 0
    cluster_iminus1 <- 1

    while(identical(cluster_i, cluster_iminus1) == F && iter < iter.max){
      iter <- iter + 1			 # update iteration	
      suppressWarnings(
        ck_flex <- flexclust::cclust(coordinates,k=cluster_i, method="hardcl",weights= weights))
       # w_centers <- as.matrix(ck_flex@centers)
      cluster_iminus1 <- cluster_i
      cluster_i <- ck_flex@cluster
    } # End while
    return(ck_flex@cluster)
  } 
  
  #---- F4 End Consolidation function ----
  
  #---- F5 descfreq_New ----
  descfreq_New <- function (donnee, by.quali = NULL, proba = 0.05, row.INIT, col.INIT)
  {
    lab <- colnames(donnee)
    for (i in 1:length(lab)) lab[i] = gsub(" ", ".", lab[i])
    if (!is.null(by.quali)) {
      donnee <- as.data.frame(matrix(unlist(by(donnee, by.quali,apply, 2, sum)), ncol = ncol(donnee), byrow = TRUE))				
      rownames(donnee) <- levels(by.quali)
      marge.li <- aggregate(row.INIT, by=list(by.quali), FUN="sum")
      marge.li <- as.vector(marge.li$x)
    }
    colnames(donnee) = lab
    old.warn = options("warn")
    options(warn = -1)				
    nom = tri = structure(vector(mode = "list", length = nrow(donnee)), names = rownames(donnee))				
    marge.col=col.INIT
    SumTot=sum(marge.li)
    nr <- nrow(donnee)
    nc <- ncol(donnee)
    for (j in 1:nr) {
      aux3 = marge.li[j]/SumTot
      for (k in 1:nc) {
        aux2 = donnee[j, k]/marge.col[k]
        
        if (aux2 > aux3) # Over representation
          aux4 <- phyper(donnee[j, k] - 1, marge.col[k],SumTot - marge.col[k], marge.li[j],lower.tail=FALSE)*2
        else aux4 = phyper(donnee[j, k], marge.col[k], SumTot - marge.col[k], marge.li[j])*2
        
        if (aux4 > 1)
          aux4 <- 2 - aux4
        if (aux4 < proba) {
          aux5 <- (1-2 * as.integer(aux2 > aux3)) * qnorm(aux4/2)
          aux1 = donnee[j, k]/marge.li[j]
          tri[[j]] = rbind(tri[[j]], c(aux1 * 100, sum(marge.col[k])/SumTot *100, donnee[j, k], marge.col[k], aux4, aux5))
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
        colnames(tri[[j]]) = c("Intern %", "glob %", "Intern freq","Glob freq ", "p.value", "v.test")
      }
    }
    names(tri) <-paste("cluster",1:length(tri),sep="_")
    res <- tri
    options(old.warn)
    class(res) <- c("descfreq", "list ")
    return(res)				
  }	
  #---- End descfreq_New function -------------------------------
  
  
  
  ########################################################################################  
  # Starting analysis ------------------------------------------------------------  
  ######################################################################################## 


   #---- 1. Weighting documents (after/before) and words ----
  doc.w.after <- apply(object$call$X,1,sum)    # Initial weighting for documents after selection
  #nm <- names(doc.w.after)
  # doc.w.before before word's selection
  #doc.w.before <- object$rowINIT[which(rownames(object$rowINIT) %in% nm) ,]
  #names(doc.w.before) <- nm
  word.w <- apply(object$call$X, 2, sum) 
    #if(marg.doc=="before") doc.w <- doc.w.before else doc.w <- doc.w.after # before word's selection


  #---- 2. Selecting words or documents to cluster  ----
  if (cluster.CA == "docs") {
    weight.dw <- doc.w.after # 
    coord.dw <- object$row$coord
  } else {
    weight.dw <- word.w
    coord.dw <- object$col$coord
  }   
 
  
  #---- 3. Building PCA object ----
  # Check if kmeans will be performed
  # Operations with kk
  #  If kk = Inf => hierarchical cluster
  #  If kk is not Inf => k-means cluster. If k-means coordinates * sqrt(eigenvalue)   
  ################# Selecting cluster of rows (documents) or columns (words) 
  cla <- res <- NULL	
  
  
  
  if (kk < nrow(coord.dw)) {  # Kmeans
    res <- as.data.frame(sweep(coord.dw, 2, sqrt(object$eig[1:ncol(coord.dw), 1]), FUN = "*"))
    kk <- min(kk, nrow(unique(coord.dw)))
  } else kk=Inf				


   if (!is.null(res)) {	# It is k-means. res is a dataframe coordinates active documents (or words) x dimensions	
    if (kk < nrow(res)) kk <- min(kk, nrow(unique(res)))	
    if (kk <= nrow(res)) {	
      # Weighted k-means (no with Factoclass) Do it with flexclust   
      set.seed(seed)
      cla <- flexclust::cclust(res,k=kk, method="hardcl",weights= weight.dw)
      cla <- data.frame(cla@cluster) # documents/words win names x 1 column whith number of cluster	
      cla <- cbind(res, weight.dw ,cla[rownames(res),,drop=FALSE])
      cla1 <- cbind(cla[,c(1:(ncol(cla)-2))]*cla[,(ncol(cla)-1)],cla[,c((ncol(cla)-1):ncol(cla))])
      cla.agg <- aggregate(cla1[,c(1:(ncol(cla1)-1))], by=list(cla.cluster=cla$cla.cluster), FUN=sum)	
      cla.centers <-  cla.agg[,c(2:(ncol(cla.agg)-1))] /cla.agg[,ncol(cla.agg)]  
      res <- FactoMineR::PCA(cla.centers, row.w = cla.agg[,ncol(cla.agg)], scale.unit = FALSE, ncp = Inf, graph = FALSE)
      rm(.Random.seed,envir=globalenv()) # Retrieve set.seed
    }          # end      if (kk <= nrow(res)) 
  } # Final of k-means
  

  # If kk >= nrow(word) now a hierarchical cluster is performed, res don't exist
  if(kk==Inf){ ## It is hierarchical cluster, no k-means
    aux <- object$eig			
    res <- FactoMineR::PCA(coord.dw, scale.unit = FALSE, ncp = Inf, graph = FALSE, row.w = weight.dw)	
    res$eig <- aux			
  }  # Now, all res object are PCA objects type;  End Building PCA object



    
  
  
  #---- 4. Checking the minimum and maximum number of clusters ----
  # max min  "ind" is documents if cluster.CA="docs", words if cluster.CA="words"
  nb.ind <- nrow(coord.dw)			## Cases (docs o words) 

  if (is.null(max))    max <- min(10, round(nb.ind/2))			
  max <- min(max, nb.ind  - 1)  
  if(max<2) {
    stop("The maximum number of clusters (",max,") can not be lower than 2")
  }
  if(max < min) {
    warning("The minimum number of clusters (",min,") can not be bigger than the maximum (",max,"). 
              Automatic change: min=max")				
    min<- max
  }   # End max min
  

  

  
  

  
  #---- 5. Building the tree to select number clusters -----   
  # Retrieve the tree with cases if hierarchical or groups if k-means	
  # row.w.init has the weighting from rows from res object of PCA
  # method=ward, metric=euclidean , res$call$row.w.init is absolute frequency after
  
    t <- auto.cut.tree(res, min = min, max = max, metric = metric,method.sel = method.sel,
                     weight = res$call$row.w.init, cla = cla, order = order, ...)	
 
   # tree can be plotted plot(t$tree)
  # We need to select the number of clusters
  # Only is necessary the plot if nb.clust = click (0 value)

  if (graph.scale == "inertia") {
    nb.ind.clI <- nrow(t$res$ind$coord)
    inertia.height <- rep(0, nb.ind.clI - 1)							
    for (i in 1:(nb.ind.clI - 1)) inertia.height[i] <- t$inert.gain[(nb.ind.clI -i)]							
    inertia.height <- sort(inertia.height, decreasing = FALSE)	
    t$tree$height <- inertia.height							
  }			# end of graph.scale =="inertia"
  
  
  
  
  #---- 6. Check number clusters (max,min) before selection ----
  options(warn=0)
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
    } 
  
  if(nb.clust>max) {
    warning("The number of selected clusters is bigger than the number of maximum clusters.
              Automatic change: max=", min(nb.clust, nrow(res$ind$coord) - 1))				
    max<- min(nb.clust,nrow(res$ind$coord) - 1)
    nb.clust <- max  # New instruction
  }    
  
 
  
  #---- 7. Height for proposed number of clusters, necessary to plot for selection ----
  auto.haut <- ((t$tree$height[length(t$tree$height) - t$nb.clust + 2]) +		
                  (t$tree$height[length(t$tree$height) - t$nb.clust + 1]))/2 
  
  #---- 8. Print Inertia gain graph for selection case ----
  if (nb.clust==0 ) {
    par(mfrow=c(1,1))			
    print.inertia("Hierarchical Clustering")
    print("Click on the graph to cut the tree")	
    flush.console()							
    plot(t$tree, hang = -1, main = "Click to cut the tree", xlab = "", sub = "")
    text(1, 1, "Hierarchical Clustering", "cex" = 2)	
    abline(h = auto.haut, col = "black", lwd = 3)	
    coupe <- locator(n = 1)							
    while (coupe$y < min(t$tree$height)) {							
      cat("No class \n")							
      coupe <- locator(n = 1)							
    }						    # End while	
    y <- coupe$y							
    # y is a value of height, value of xy axis							
    # If it is not manual selections inertia gain graph is not plotted
  } else {
    if (nb.clust < 0) 	
      y <- auto.haut	
    else     y <- (t$tree$height[length(t$tree$height) - nb.clust + 2] + 
                t$tree$height[length(t$tree$height) - nb.clust + 1])/2	
  }

  # y has the point of cut
  # composition of clusters (if kmeans not the cases)
  clust <- cutree(as.hclust(t$tree), h = y)	

  # integer with the number of the cluster
  # The name is the label of the doc/word 
  nb.clust <- max(clust)	
  
  
### REVISAR SI ordColo es necesario después
  # Order of clusters in dendrogram (duplicates are removed)
  ordColo <- unique(clust[t$tree$order])
  
  
  ### REVISAR SI X es necesario después
    X <- as.data.frame(t$res$ind$coord)	   # Coordinates Docs or words x dimensions, If kmeans groups x dimensions
  #---- End Print Inertia gain graph for selection case ----
  
  
  #---- 9. df.coord.cla docs x  (Dim.1;...;Dim k; weight cluster ----
  if(kk==Inf) {  # Jerárquico
    df.coord.cla <- data.frame(X, weight= res$call$row.w.init[rownames(X)], cluster=clust[rownames(X)])
  }
  if(kk!=Inf) {   # Kmedias
    df.coord.cla <-   data.frame(coord.dw, weight= cla[rownames(coord.dw),"weight.dw"], cluster=clust[cla$cla.clust])                               #cluster=cla[rownames(coord.dw),"cluster"])
  }


    
  #---- 10 Print dendrogram - Selection of the number of the clusters has been done, cla ----
  if(graph) if(kk!=Inf) {
    warning("Warning: Dendrogram is not plotted for kmeans method")
   # graph <- FALSE
  }
  if(graph) if(consol) {
    warning("Warning: Dendrogram is not plotted for consolidation method")
  #  graph <- FALSE
  }

  if(graph & kk==Inf) {
    par(mfrow=c(1,1))
    print.inertia("Hierarchical Clustering")
    plot(t$tree, hang = -1, main = "", xlab = "", sub = "")
    text(1, 1, "Hierarchical Clustering", "cex" = 2)
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
  
  
    
    
  #---- 11. Consolidation if hierarchical ----			
 if(consol==TRUE){
   clust.afterconsol <- consolidationW(df.coord.cla[,c(1:(ncol(df.coord.cla)-2))], cluster_i=df.coord.cla[,"cluster"],
                                               weights=df.coord.cla[,"weight"] , iter.max = iter.max)
   df.after <- as.data.frame( clust.afterconsol ) 
   dftemp<- data.frame(df.coord.cla[,c(1:(ncol(coord.dw)))],  weight=df.coord.cla[,"weight"], cluster=df.coord.cla[,"cluster"])
   dftemp$cluster <- df.after[rownames(dftemp),"clust.afterconsol"]
   df.coord.cla <- dftemp
}

    
  #---- 12. Computing centers. Building RDO object  ----			
  RDO <-   ssquares(df.coord.cla[,c(1:(ncol(coord.dw)))],  df.coord.cla[,"weight"], df.coord.cla[,"cluster"]) 
  
  RDO <- list("data.clust"=data.frame(df.coord.cla[,c(1:(ncol(coord.dw)))], "clust"=df.coord.cla[,"cluster"]),   
              "centers"=RDO$coord.centers, "clust.count"=RDO$dfclust, "clust.content"=RDO$dfw,
              "global"= list("ss.tot"=RDO$ss.tot, "ss.intra"=RDO$ss.intra,"ss.inter"= RDO$ss.inter,			
                             "%Explained"= 100*RDO$ss.inter/RDO$ss.tot), cluster.CA=cluster.CA)
  RDO$CA <- object
  class(RDO) <- "LexHCca"
  # End 12. Computing centers. Building RDO object 
  

  
  
  
  
  

  #---- 14. Weighting documents (after/before) and words ----
  #  doc.w.after <- apply(object$call$X,1,sum)    # Initial weighting for documents after selection
  nm <- names(doc.w.after)
  # doc.w.before before word's selection
  doc.w.before <- object$rowINIT[which(rownames(object$rowINIT) %in% nm) ,]
  names(doc.w.before) <- nm
  # word.w <- apply(object$call$X, 2, sum) 
  if(marg.doc=="before") doc.w <- doc.w.before else doc.w <- doc.w.after # before word's selection
  # num.docs <- nrow(object$call$X)
  # num.words <- ncol(object$call$X)
  # Final 14 weighting docs
 
#  doc.w.after <- apply(object$call$X,1,sum)    # Initial weighting for documents after selection
  
  
  
  
  # tabInd coordinates + cluster if active docs/words
  if (cluster.CA == "docs"){ 
    descr.data <- object$call$X   # words in columns
    descr.cl <- data.frame("cluster_"=df.coord.cla[rownames(descr.data),"cluster"])  # drop=FALSE
    rownames(descr.cl) <- rownames(descr.data)
    tabInd <- cbind.data.frame(object$row$coord, clust= as.factor(df.coord.cla[rownames(descr.data),"cluster"]))
  }else{ 
    descr.data <- t(object$call$X)
    descr.cl <- data.frame("cluster_"=df.coord.cla[rownames(descr.data),"cluster"] ) # drop=FALSE
    rownames(descr.cl) <- rownames(descr.data)
    tabInd <- cbind.data.frame(object$col$coord, clust=as.factor(df.coord.cla[rownames(descr.data),"cluster"]))
  }

  
  
  
  
  
  
  desc.ind <- descaxes <- descword <- descwordsup <- descdoc <- descquali <- descquanti <- NULL
  #---- 15. descword, Cluster description by words and descdoc ----
  if(description) { 
    
    if(cluster.CA=="docs") {
      descr.w.rows <- doc.w[rownames(descr.data)]
      descr.w.cols <- word.w[colnames(descr.data)]
      # Cluster Description by active words
      descword <- descfreq_New(descr.data, descr.cl$cluster_,prob=proba, row.INIT= descr.w.rows, col.INIT=descr.w.cols) 
    } else { 
      descr.w.rows <- word.w[rownames(descr.data)]
     # descr.w.cols <- doc.w[colnames(descr.data)]
      # Cluster Description by active documents # Using doc.w.after
      descdoc <- descfreq_New(descr.data, descr.cl$cluster_ ,prob=proba, row.INIT= descr.w.rows, col.INIT=doc.w.after) 
    } 
    #----  End descword, Cluster description by words and descdoc ----    


    
    
    #---- Document description by qualitative supplementary variables ----
    ### Por defecto es before...
    if(!is.null(object$quali.sup$coord)) {
      if(cluster.CA=="docs") {
    #    descr.cl <- data.frame("cluster_"=df.coord.cla[rownames(descr.data),"cluster"])  # drop=FALSE
        a1 <- object$call$Xtot[rownames(descr.data), rownames(object$quali.sup$eta2), drop=FALSE]
        a2 <- cbind(a1, "clust_"=as.factor(descr.cl[rownames(a1),"cluster_"]))
         descquali <- FactoMineR::catdes(a2, num.var=ncol(a2), proba = proba, row.w=doc.w[rownames(a2)])
         names(descquali$category) <-paste("cluster",1:length(descquali$category),sep="_")
      }}
    #---- End Document description by qualitative supplementary variables ----    
   
    #---- Document description by quantitative supplementary variables ----
    if(!is.null(object$quanti.sup$coord)) {
    
      if(cluster.CA=="docs") {
        a1 <- object$call$Xtot[rownames(descr.data), rownames(object$quanti.sup$coord), drop=FALSE]
        a2<- cbind(a1, "clust_"=as.factor(descr.cl[rownames(a1),"cluster_"]))
        descquanti <- FactoMineR::catdes(a2, num.var=ncol(a2), proba = proba, row.w=doc.w[rownames(a2)])
        if(!is.null(descquanti$quanti)>0)  
          names(descquanti$quanti) <-paste("cluster",1:length(descquanti$quanti),sep="_")
      }}
    #---- End Document description by quantitative supplementary variables ----
    
    
    #---- Description of axes ----
    if(cluster.CA=="docs") {
      descaxes <- data.frame(cbind(object$row$coord,clust_=as.factor(descr.cl[rownames(object$row$coord),"cluster_"])))
      descaxes$clust_ <- as.factor(descaxes$clust_)
      descaxes <- FactoMineR::catdes(descaxes, num.var=ncol(descaxes), proba = proba, row.w=doc.w[rownames(descaxes)])  
    } else { # columnas
      descaxes <- data.frame(cbind(object$col$coord,clust_=as.factor(clust[rownames(object$col$coord)])))
      descaxes$clust_ <- as.factor(descaxes$clust_)
      descaxes <- FactoMineR::catdes(descaxes, num.var=ncol(descaxes), proba = proba, row.w=word.w[rownames(descaxes)]) 
    }
    
 
  
    if(nb.desc>0){
      # Calculo de los centros
      # list.centers <- by(tabInd[, -ncol(tabInd), drop = FALSE], tabInd[, ncol(tabInd)], colMeans)
      # centers <- matrix(unlist(list.centers), ncol = ncol(tabInd) - 1, byrow = TRUE)
      #---- Function select --------------------------------
      select <- function(Y, default.size, method, coord.centers) {
        clust <- as.numeric(Y[1, ncol(Y)])
        Y <- Y[, -ncol(Y), drop = FALSE]
        colnames(coord.centers) <- colnames(Y)
        Z <- rbind(Y, coord.centers)

        if (nrow(Y) == 1) {  #No elements
          distance <- data.frame(0, row.names = "")
          colnames(distance) <- rownames(Z[1, ])
        }
        else {
          distance <- as.matrix(dist(Z, method = method))
          distance <- distance[(nrow(Y) + 1):nrow(distance), 
                               -((nrow(Y) + 1):ncol(distance))]
          distance <- sort(distance[clust, ], decreasing = FALSE)
        }
        if (length(distance) > default.size) 
          distance <- distance[1:default.size]
        else distance <- distance
      }

            #---- Function distinctivness  ---------------------
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
          distance <- distance[(nrow(Y) + 1):nrow(distance),-((nrow(Y) + 1):ncol(distance))]
          if (nrow(distance) == 2)
            center.min <- distance[-clust, ]
          else center.min <- apply(distance[-clust, ], 2, min)
          ind.car <- sort(center.min, decreasing = TRUE)		
        }	
        if (length(ind.car) > default.size)
          ind.car = ind.car[1:default.size]
        else ind.car = ind.car
      } # Final distinctivness
      #---- End Function distinctivness  ---------------------
      
      m_centers <- data.matrix(RDO$centers)[,-1]
      Cluster <- tabInd[,ncol(tabInd), drop=FALSE]
      para <- by(tabInd, Cluster, simplify = FALSE, select, 
                 default.size = nb.desc, method = metric, coord.centers = m_centers)
      dist <- by(tabInd, Cluster, simplify = FALSE, distinctivness, 
                 default.size = nb.desc, method = metric, coord.centers = m_centers)
      desc.ind <- list(para = para, dist = dist)

    } # End nb.par
  }  
  #---- Final Descriptions ----
  

  newcall <- list(t = t, min = min, max = max, X = tabInd[,-ncol(tabInd)], vec = FALSE, call = match.call(), 
                  cluster.CA=cluster.CA)  
  
  df.tmp <- cbind(descr.data,descr.cl)
  colnames(df.tmp)[ncol(df.tmp)] <- "Clust_"

  colnames(RDO$centers)[1] <- "Clust_" 
  newres <- list(data.clust = df.tmp, centers =as.data.frame(RDO$centers), clust.count = RDO$clust.count, 
                 clust.content = RDO$clust.content, ss=RDO$global, call = newcall)
  

  if (!is.null(descaxes)) {
    class(descaxes) <- NULL
    newres$description$desc.axes <- descaxes
  }
  if (cluster.CA == "docs") {
    if (!is.null(descword)) 
      newres$description$desc.cluster.doc <- list(words = descword, 
                                                  wordsup = descwordsup, qualisup = descquali, 
                                                  quantisup = descquanti)
  }
  else {
    if (!is.null(descdoc)) 
      newres$description$desc.cluster.word <- list(docs = descdoc)
  }
  
  

  
  
  if (!is.null(object$info)) {
    if (nb.desc > 0) {
      if (is.null(desc.ind$para)) 
        stop("Please, use description=TRUE \n  to obtain a description of the clusters by the characteristic documents ")
     

      if(!is.null(object$var.agg)) vbaggr <- TRUE else vbaggr <- FALSE
      
      if(vbaggr==FALSE) {
      var.text <- object$info$var.text[[1]]
      str.base <- object$info$base[[1]]
      str.envir <- object$info$menvir[[1]]
      base <- get(str.base, envir = str.envir)
      corpus <- base[, var.text[1]]
      if (length(var.text) > 1) {
        for (i in 2:length(var.text)) {
          corpus <- paste(corpus, base[, var.text[i]], sep = ".")
        }
      }
      }
      if(vbaggr==FALSE) {
      corpus <- data.frame(corpus, stringsAsFactors = FALSE)
      rownames(corpus) <- rownames(base)
      }
      
     
      lispara <- vector(mode = "list")
   
      for (iclus in 1:nb.clust) {
        doctot <- length(desc.ind$para[[iclus]])
        ntdoc <- min(nb.desc, doctot)
        lispara[[iclus]] <- data.frame()
        for (i in 1:ntdoc) {
          lispara[[iclus]][i, 1] <- names(desc.ind$para[[iclus]][i])
          lispara[[iclus]][i, 2] <- desc.ind$para[[iclus]][i]

          if (cluster.CA == "docs") { if(vbaggr==FALSE)
            lispara[[iclus]][i, 3] <- strtrim(corpus[lispara[[iclus]][i,1], 1], size.desc)
            else lispara[[iclus]][i, 3] <- ""
          }
        }

        if (cluster.CA == "docs") {
          colnames(lispara[[iclus]]) <- c("DOCUMENT", "CRITERION", "TEXT") 
        }
        else colnames(lispara[[iclus]]) <- c("WORD", "CRITERION")
      }
      names(lispara) <- paste("cluster", 1:nb.clust, sep = "_")
      if (cluster.CA == "docs") 
        newres$description$desc.cluster.doc$para <- lispara
      else newres$description$desc.cluster.word$para <- lispara
      lisdist <- vector(mode = "list")
      for (iclus in 1:nb.clust) {
        doctot <- length(desc.ind$dist[[iclus]])
        ntdoc <- min(nb.desc, doctot)
        lisdist[[iclus]] <- data.frame()
        for (i in 1:ntdoc) {
          lisdist[[iclus]][i, 1] <- names(desc.ind$dist[[iclus]][i])
          lisdist[[iclus]][i, 2] <- desc.ind$dist[[iclus]][i]
          if (cluster.CA == "docs") 
          { if(vbaggr==FALSE)
            lisdist[[iclus]][i, 3] <- strtrim(corpus[lisdist[[iclus]][i, 1], 1], size.desc)
          else lisdist[[iclus]][i, 3] <-""
        }}
        if (cluster.CA == "docs") 
          colnames(lisdist[[iclus]]) <- c("DOCUMENT","CRITERION", "TEXT")
        else colnames(lisdist[[iclus]]) <- c("WORD", "CRITERION")
      }
      names(lisdist) <- paste("cluster", 1:nb.clust, sep = "_")
      if (cluster.CA == "docs") 
        newres$description$desc.cluster.doc$dist <- lisdist
      else newres$description$desc.cluster.word$dist <- lisdist
    }
  }
  newres$call$CA <- object
  newres$coord.clust <- cbind(newres$call$X, "Clust_"=newres$data.clust[rownames(newres$call$X),"Clust_"])
  
  class(newres) <- c("LexHCca")
  #---- 13. Plot Factorial Plane with centers ----
  if(graph){
    p<-plot.LexHCca(newres, type="map", draw=c("points","centers"))
    plot(p)
  }
  # 13. Plot Factorial Plane with centers
  
return(newres)

}
