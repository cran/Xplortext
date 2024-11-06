#' @importFrom graphics abline layout legend locator par plot text
#' @importFrom grDevices palette
#' @importFrom vegan mantel
#' @export


LexCHCca <- function (object, ncp=5, nb.clust=0, min=2, 
                      max=NULL, nb.par=5, graph=TRUE, proba=0.05, cut.test = FALSE, alpha.test =0.05,
                      description=FALSE, nb.desc=5, size.desc=80) 
{
  x <- object
  ######################  New. Modification of x removing no factors used ######################
  ncp.max <- ncol(x$row$coord)
  if(ncp>ncp.max)   {
    warning(paste0("Number of components ncp= ", ncp ," is bigger than number of components in ",
                   deparse(substitute(object)), "= ", ncp.max),"\nncp is changed to ",ncp.max)			
    ncp <- ncp.max
  } 

  x$row$coord <- x$row$coord[, 1:ncp]
  x$row$contrib <- x$row$contrib[, 1:ncp]
  x$row$cos2 <- x$row$cos2[, 1:ncp]
  x$row$inertia <- x$row$inertia[1:ncp]
  
  x$col$coord <- x$col$coord[, 1:ncp]
  x$col$contrib <- x$col$contrib[, 1:ncp]
  x$col$cos2 <- x$col$cos2[, 1:ncp]
  x$col$inertia <- x$col$inertia[1:ncp]
  
  if(!is.null(x$quanti.sup$coord)) {
    x$quanti.sup$coord <- x$quanti.sup$coord[, 1:ncp]
    x$quanti.sup$cos2 <- x$quanti.sup$cos2[, 1:ncp]   }
  
  
  if(!is.null(x$quali.sup$coord)) {
    x$quali.sup$coord <- x$quali.sup$coord[, 1:ncp]
    x$quali.sup$cos2 <- x$quali.sup$cos2[, 1:ncp] 
    x$quali.sup$v.test <- x$quali.sup$v.test[, 1:ncp] 
    x$quali.sup$eta2 <- x$quali.sup$eta2[, 1:ncp] }
  
  if(!is.null(x$meta$Word)) x$meta$Word <- x$meta$Word[as.numeric(x$meta$Word$Dim)<(ncp+1),]
  if(!is.null(x$meta$Doc)) x$meta$Doc <- x$meta$Doc[as.numeric(x$meta$Doc$Dim)<(ncp+1),]
  ##############################################################################################
  
  # cut.stat. NULL, percent, median
  if(cut.test == FALSE) cut.stat <- NULL else cut.stat <- "percent"
  
  options(warn=0)
  marg.doc<- "after"
  cor.stat <- "pearson"                      # Or spearman, not included in the function
  if (!inherits(x,"LexCA")) stop("Object x should be LexCA class")
  options(stringsAsFactors = FALSE)
  
  if(proba<0 | proba >1) stop("No valid values for proba argument")
  
 # Check if Legendre Test is applied
  if(is.null(cut.stat)) bTest <- FALSE
  else {
    if(!is.null(alpha.test)) {
      if(alpha.test<0 | alpha.test >1) stop("No valid values for alpha.test argument")
    }
    bTest<- TRUE
    if(cut.stat!="percent" & cut.stat!="median")
      warning(paste0("Not allowed value ", cut.stat, " for cut.stat argument, changed to 'percent'"))
      cut.stat <- "percent"
  }
  bAddTest <- FALSE


#  stop("*****  FASE A")
  ## Type selection number of clusters
  # If Legendre test is applied, the number of clusters is the number obtained by the test.
  if(is.character(nb.clust)) {
    if(nb.clust == "auto") nb.clust <- -1
    if(nb.clust == "click") nb.clust <- 0
  }
      #  nb.clust <- ifelse(nb.clust == "auto", -1,0)
  # if(nb.clust=="click") nb.clust<- 0
  #  if(nb.clust=="auto") nb.clust<- -1
  if(nb.clust==0) graph <-TRUE  
  
  
  if (nb.clust == 1)  #
  {
    stop("The number of clusters can not be 1. Automatic change: nb.clust= 'auto'")				
    }  
  nb.cluster.test <- 0  
  
  
  
  if(bTest==TRUE & nb.clust != -1) {
    warning("Only automatic results are allowed using Test in cut.stat.\n nb.clust is changed to 'auto' value")
    nb.clust < -1
  }
  

  
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
  # Starting functions
   #---- F1. sum of squares ----
  ssquares <- function(coord, weight_,clust_){ # }, weight_, clust_) {
    dfw <- data.frame(coord,"weight_"=weight_, "clust_"=clust_)
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
    ###############################  REVISAR
      colnames(dfclust) <- c("Cluster", "Docs", "Occurrences", "ss.tot", "ss.intra", "ss.between", "%Explained")							
    dfw <- dfw[,c((ncol(dfw)-5):ncol(dfw))]
    colnames(dfw) <- c("Weight","Cluster", "ss.tot", "ss.intra", "ss.between","%Explained")
    
    res.ss <- list("ss.tot"=ss.tot, "ss.intra"=ss.intra,"ss.inter"= ss.inter, "coord.centers"=coord_centers, "dfw"=dfw,
                   "dfclust"=dfclust)
    return(res.ss)
  }  # End function F1. sum of squares
  
  
  ### function F2 select
  select <- function(Y, default.size, method, coord.centers) {
    clust <- Y[1, ncol(Y)]
    Y <- Y[, -ncol(Y), drop = FALSE]
    Z <- rbind(Y, coord.centers)
    if (nrow(Y) == 1) {
      distance <- data.frame(0, row.names = "")
      colnames(distance) <- rownames(Z[1, ])
    }
    else {
      distance <- as.matrix(dist(Z, method = method))
      distance <- distance[(nrow(Y) + 1):nrow(distance),-((nrow(Y) + 1):ncol(distance))]
      distance <- sort(distance[clust, ], decreasing = FALSE)
    }
    if (length(distance) > default.size) 
      distance <- distance[1:default.size]
    else distance <- distance
  }
  
  ### function F3 distinctivness
  distinctivness = function(Y, default.size, method, coord.centers) {
    clust <- as.numeric(Y[1, ncol(Y)])
    Y <- Y[, -ncol(Y), drop = FALSE]
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
  }
  
  
  ### function F4 auto.cut.constree
  auto.cut.constree <- function(res, min, max, alpha.test) 
  {
    X <- as.data.frame(res$ind$coord)  # Documents x Dimensions
    
    hcc <- hcclust(X, alpha.test)
    height <- hcc$height
    inv.height <-rev(height)
    quot<-inv.height[(min-1):(max-1)]/inv.height[(min):(max)]
    
    if(is.null(cut.stat))
    { nb.clust <- which.max(quot) + min - 1 } else nb.clust <-  nb.cluster.test
    #  stop(nb.clust)
    return(list(res = res, tree = hcc, nb.clust = nb.clust,height=height,quot = quot))
  }  # End auto.cut.constree function
  
  
  ### function F5 hcclust
  hcclust <- function(X, alpha.test)  
  {
    options(warn=-1)
    ## PHASE 1 ##
    ndoc<-nrow(X)                  # Number of documents to cluster
    d<-dist(X)                     # Euclidean distance of LexCA coordinates of documents as matrix format. Mdist
    mDist.0 <- mDist<-as.matrix(d)

    # mSim is the original similarity matrix 
    maxd<-max(mDist)+1e-10
    mSim<-(maxd-mDist)
    mCont<-rep(0,ndoc*ndoc)
    attr(mCont,"dim") <- c(ndoc, ndoc)   # matrix, array
    # If only restricted contiguity clusters
    for(i in 1:(nrow(mCont)-1)) mCont[i,i+1] <- 1
    for(i in 2:(nrow(mCont))) mCont[i,i-1] <- 1
    rownames(mCont) <- colnames(mCont) <- colnames(mSim)
    ###   Introduction of the constriction into the similarity and distance matrices
    mSimCont<-mSim*mCont
    mDistCont<-mDist*mCont
    maxMsim <- max(mSim)
    # detail. For explaining who is joined in each step
    detail <- data.frame(n.step=character(),joined=logical(),cluster.Left=character(), cluster.Right=character())
    ## END PHASE 1 ##
    

    
    ## PHASE 2 ##
    ###    Groups existing at Phase 0. Empty list with as many items as individuals
    groups <- list()                           
    for (i in 1:nrow(mSim)) groups[[i]]<-i
    ###    Critagreg will file the successive values of the agregation criterium
    Critagreg<-numeric()                     # value of the distance in the (n-1) fusions of nodes
    ###  clust will archive the successive pairs of nodes that are merged
    clust<-list()
    ###    Building the hierarchy: (n-1) intermediate nodes have to be created
    i<-1
    # indice<-nrow(Mdist)-1
    
  
    #    if(bTest) {
    ### ONLY WHEN alpha.test !=0. mSingletons save if one individual is a singleton
    #    mSingletons <- matrix(0, nrow = nrow(mCont), ncol = 2)
    #    colnames(mSingletons) <- c("Left", "Right")
    #    rownames(mSingletons) <- rownames(mCont)
    #    mSingletons[1,1] <- mSingletons[nrow(mDist),2] <-1
    #  write.csv2(mSingletons,file= "C:/LexCHCCA/mSingletons.csv")
    ### End ONLY WHEN alpha.test !=0
    #    } # End mSingletons
    
    ## END PHASE 2 ##
    
    ############ Pendiente esto, solamente si bTest 
    indice <- 1

    #####################################################################################################################################################
    # creating the n-1 nodes 
    #  while(indice>0) {	 
    while(sum(mCont)>0) {	                                 ### PHASE 3. DO UNTIL  AGGREGATION IS POSSIBLE
      #        Indentifying the maximum similarity
      maxsim<- max(mSimCont)  
      mat_pos <- which(mSimCont == maxsim,  arr.ind = TRUE)
      minpos <- min(mat_pos[1,])
      maxpos <- max(mat_pos[1,])
      # minpos is the row and maxps is the column to group with the similarities
      

      ##############  Checking Legendre test ib bTest=TRUE ##############
      if(bTest) {                                                 ### PHASE 4 ########
        # n.minpos is the number of cases in the cluster int the position minpos
        # n.maxpos is the number of cases in the cluster int the position maxpos
        n.minpos <- length(groups[[minpos]])
        n.maxpos <- length(groups[[maxpos]])
        

        
        if(n.minpos+n.maxpos==2) {   # Each one of the two checked clusters has only one case.
          bAddTest <- TRUE   # Aggregation of two single clusters
          detail <- rbind(detail,data.frame(n.step=indice,bAddTest, cluster.Left=colnames(mDist)[minpos], cluster.Right=colnames(mDist)[maxpos] ) )
                 }
        else{
          mDistExtrem <- cbind(matrix(0,nrow=n.minpos, ncol=n.minpos),
                               matrix(1,nrow=n.minpos, ncol=n.maxpos))
          mDistExtrem <- rbind(mDistExtrem, cbind(t(matrix(1,nrow=n.minpos, ncol=n.maxpos)),
                                                  matrix(0,nrow=n.maxpos, ncol=n.maxpos)))
          
          mDist.CL <- mDist.0[c(groups[[minpos]],groups[[maxpos]]),
                              c(groups[[minpos]],groups[[maxpos]])]
          
          if(cut.stat =="percent") {
            n.lower <-  ((( n.minpos+n.maxpos)^2)-( n.minpos+ n.maxpos))/2
            n.ones <- n.minpos*n.maxpos
            cut.stat.value <- quantile(mDist.CL[row(mDist.CL)>col(mDist.CL)],
                                       probs = 1-n.ones/n.lower)
          } else { # End "percent
            cut.stat.value <- median(mDist.CL[row(mDist.CL)!=col(mDist.CL)])
            
          } # End if(cut.stat =="percent")
          ind01 <- mDist.CL >= cut.stat.value 
          ind01 <- +ind01
          
          suppressMessages(Mantel.res <- vegan::mantel(mDistExtrem, ind01, method=cor.stat, permutations=9999))
          
          bAddTest <- ifelse(Mantel.res$signif < alpha.test,FALSE,TRUE)  ### Only one tail greater
          
          detail <- rbind(detail,data.frame(n.step=indice,bAddTest, cluster.Left=colnames(mDist)[minpos], 
                                            cluster.Right=colnames(mDist)[maxpos] ) ) 
          if(bAddTest==FALSE) 
            mCont[minpos, maxpos] <- mCont[maxpos, minpos] <- 0 
          
        } # End  if(n.minpos+n.maxpos==2) End  #### Hay test  
      } # End bTest is TRUE
      
      

      if(bTest==FALSE | bAddTest==TRUE) {   # Adding clusters to dendrogram
        Critagreg[i]<- mDistCont[minpos,maxpos]
        clust[[i]] <- vector(mode="list", length=2)
        clust[[i]][[1]] <- groups[[minpos]] 
        clust[[i]][[2]] <- groups[[maxpos]]

    
        
        colnames(mSim)[minpos]<-rownames(mSim)[minpos] <- colnames(mDist)[minpos]<-rownames(mDist)[minpos] <- 
        colnames(mCont)[minpos]<-rownames(mCont)[minpos] <- paste0(rownames(mSim)[minpos],"-",colnames(mSim)[maxpos])
        ###### Lance & Williams. Values for complete4 link (CLCA) alpha=0.5 ; gamma=.5
        # Actualization of Msim, Mdist and Mcont
        # When for distances alfa=0.5 and gamma  = 1/2 (Complete linkage). Joining of maximum distances
        # When for similarities alfa =0.5 and gamma = -1/2 (single linkage). Joining minimum similarities
        #   fila <- minpos; columna <- maxpos
        #  if (fila!=1) {
        #      Msim[fila,fila-1]<-Msim[fila-1,fila]<- 0.5*Msim[fila,fila-1]+0.5*Msim[columna,fila-1]-0.5*abs(Msim[fila,fila-1]-Msim[columna,fila-1])
        #      Mdist[fila,fila-1]<-Mdist[fila-1,fila]<-0.5*Mdist[fila,fila-1]+0.5*Mdist[columna,fila-1]+0.5*abs(Mdist[fila,fila-1]-Mdist[columna,fila-1])
        #      Mcont[fila-1,fila]<-Mcont[fila,fila-1]<-1
        #   }
        #    if (columna!=nrow(MSimCont)) {
        #      Msim[columna+1,fila]<-Msim[fila,columna+1]<-0.5*Msim[fila,columna+1]+0.5*Msim[columna,columna+1]-0.5*abs(Msim[fila,columna+1]-Msim[columna,columna+1])
        #     Mdist[columna+1,fila]<-Mdist[fila,columna+1]<-0.5*Mdist[fila,columna+1]+0.5*Mdist[columna,columna+1]+0.5*abs(Mdist[fila,columna+1]-Mdist[columna,columna+1])
        #      Mcont[columna+1,fila]<-Mcont[fila,columna+1]<-1
        #    }
        
        # New algorithm, only when method=complete is equal to max distance and minimum similarity
        # Compute the similarities and distances for all old clusters with the new cluster 
        mSim[minpos,]<- mSim[,minpos]<- apply(X=mSim[,c(minpos,maxpos)], MARGIN=1, FUN=min)
        mSim[minpos,minpos]<- maxMsim
        mDist[minpos,]<- mDist[,minpos]<- apply(X=mDist[,c(minpos,maxpos)], MARGIN=1, FUN=max)
        mDist[minpos,minpos]<- 0
   
        
 # Check if there is a cluster on the left and on the right       
        if(maxpos!=nrow(mCont)) mCont[maxpos+1,minpos] <- mCont[minpos,maxpos+1] <-1
        if(minpos!=1) mCont[minpos-1,minpos] <- mCont[minpos,minpos-1] <-1
        

        groups[[minpos]] <- c(groups[[minpos]], groups[[maxpos]])
        # Removing rows and columns not necessary
        groups <- groups[-maxpos]
        mSim <- mSim[-maxpos,-maxpos]
        mDist <- mDist[-maxpos,-maxpos]
        mCont <- mCont[-maxpos,-maxpos]
        i <- i + 1         
      }  # End  if(bTest==FALSE | bAddTest==TRUE)
      
      mSimCont <- mSim * mCont
      mDistCont <- mDist * mCont
      indice <- indice + 1
      

      if(bTest) if(sum(mCont)==0) {    # No more possible aggregations
        nb.cluster.test <<- ncol(mCont)
        max.Dist <- max(mDist)
        Critagreg.old <- Critagreg
        if(length(Critagreg.old)== ndoc-1) stop(paste0("No clusters for alpha.test level= ", alpha.test,
                                               "\nIncrease this value"))
        Critagreg[i:(ndoc-1)]<- max.Dist

        for(k in i:(ndoc-1)) {
          clust[[k]] <- vector(mode="list", length=2)
          clust[[k]][[1]] <- groups[[1]] 
          clust[[k]][[2]] <- groups[[2]]
          groups[[1]] <- c(groups[[1]], groups[[2]])
          groups <- groups[-2]
        }
        
      }
      
    } # End     while(sum(mCont)>0)  
    
    height <- Critagreg
    
   

    ######### Eliminar los singletons después #################
    order <- as.integer(1:nrow(X))
    ##################
    grups_blocs <- list()
    grups_blocs[[1]] <- rep(0, nrow(X))
    ############
    
    
    for (i in 1:(length(clust) - 1)) {
      grups_blocs[[i + 1]] <- grups_blocs[[i]]
      grups_blocs[[i + 1]][c(clust[[i]][[1]], clust[[i]][[2]])] <- i
    }
    merge <- matrix(nrow = (ndoc - 1), ncol = 2, 0)
    
    for (i in 1:(ndoc - 1)) {
      if (length(clust[[i]][[1]]) == 1 & length(clust[[i]][[2]]) == 1) {
        merge[i, 1] <- (-clust[[i]][[1]])
        merge[i, 2] <- (-clust[[i]][[2]])
      }
      else {
        if (length(clust[[i]][[1]]) == 1) {
          merge[i, 1] <- (-clust[[i]][[1]])
        }
        else {
          merge[i, 1] <- grups_blocs[[i]][clust[[i]][[1]][1]]
        }
        if (length(clust[[i]][[2]]) == 1) {
          merge[i, 2] <- (-clust[[i]][[2]])
        }
        else {
          merge[i, 2] <- grups_blocs[[i]][clust[[i]][[2]][1]]
        }
      }
    }
    
    #  if(bTest==FALSE) cut.names <- stats::cutree()
    #  else cut.names <- rownames(mDist)
    #  
    
    if(bTest==FALSE) detail <- NULL
    
    hcccall<-call("hcclust",match.call()$res)		
    hcc<-structure(list(merge = merge,height=height,
                        order=order, labels = attr(d, "Labels"), method = "complete", 
                        call = hcccall, dist.method = "euclidean",clust=clust), class="hclust")
    options(warn=0)
      return(hcc)

    # stop("FINAL DE OBTENCION DEL DENDROGRAMA")
    
    #    dhcc <- as.dendrogram(hcc)
    #    plot(cut(dhcc, h = Critagreg[[ndoc-1]]-.000000001)$lower[[1]])
    #    plot(dhcc)
  } # End hclust
  
  
  
  
  

  
  ### main part of the function
  res.ca<- x           ##################  <- Pendiente ¿Necesario?
  res.sauv <- res.ca            ##################  <- Pendiente ¿Necesario?
  metric="euclidean"
  method="complete"             ##################  <- Pendiente
  if (min <= 1) min<-2
  

  ### PCA equivalent to CA. It is necessary to know the max number of clusters
  aux <- res.ca$eig
  weight<- row.w <- res.ca$call$marge.row * sum(res.ca$call$X)
  res <- PCA(res.ca$row$coord, scale.unit = FALSE, ncp = Inf, graph = FALSE, row.w = row.w)
  res$eig <- aux
  

    if (is.null(max)) max <- min(10, round(nrow(res$ind$coord)/2))
     max <- min(max, nrow(res$ind$coord) - 1)
  ###


  ######## auto.cut.constree call to auto.cut.constree function
  ######## auto.cut.constree call to hcclust function
  t <- auto.cut.constree(res=res,min = min, max = max, alpha.test=alpha.test) 
  # It is possible to have inversions in t$tree
  # t$tree contains the constrained dendrogram. índependentely of nb.clust
  nb.ind <- nrow(t$res$ind$coord)
  
  auto.haut <- ((t$tree$height[length(t$tree$height) - 
                                 t$nb.clust + 2]) + (t$tree$height[length(t$tree$height) - t$nb.clust + 1]))/2

  graph.1 <- function() {
    if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
      dev.new()
    old.mar <<- par()$mar
    par(mar = c(0.5, 2, 0.75, 0))
    lay <- matrix(ncol = 5, nrow = 5, c(2, 4, 4, 4, 4, 
                                        2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 
                                        1, 3, 3, 3, 3))
    layout(lay, respect = TRUE)

    barplot(rev(t$tree$height)[1:max(15, max)], col = c(rep("black", 
          t$nb.clust - 1), rep("grey", max(max, 15) - t$nb.clust +1)),
            rep(0.1, max(max, 15)), space = 0.9)
    plot(x = 1, xlab = "", ylab = "", main = "", col = "white", axes = FALSE)
    text(1, 1, "Constrained Hierarchical Clustering", cex = 2)
    plot(x = 1, xlab = "", ylab = "", main = "", col = "white", axes = FALSE)
    legend("top", "Aggreg. Crit", box.lty = NULL, cex = 1)
  }

  
  if(bTest==FALSE & nb.clust == 0) {  # click
        old.mar <- par()$mar
        graph.1()
        plot(t$tree, hang = -1, main = "Click to cut the tree", xlab = "", sub = "")
        abline(h = auto.haut, col = "black", lwd = 3)
        coupe <- locator(n = 1)
        while(coupe$y < min(t$tree$height) | coupe$y > max(t$tree$height)) {
          cat("No class selected \n")
          coupe <- locator(n = 1)           
        }
        y <- coupe$y
        par(mar = old.mar)
  } 
  
  if(bTest==FALSE & nb.clust==-1) y <- auto.haut # No test, automatic clusters
  if(bTest==FALSE & nb.clust>1 )  # No test, fixed number clusters
    y <- (t$height[length(t$height)-nb.clust+1]+t$height[length(t$height)-nb.clust+2])/2

  if(bTest==TRUE) {
    # nb.clust <- nb.cluster.test
    #  y <- (t$height[length(t$height)]+t$height[length(t$height)-nb.clust+1])/2
      clust <- stats::cutree(as.hclust(t$tree), k = nb.cluster.test)
  } else {
    clust <- stats::cutree(as.hclust(t$tree), h = y)
  } 

  

  
  nb.clust <- max(clust)
  X <- as.data.frame(t$res$ind$coord)
  ordColo <- unique(clust[t$tree$order])


  
  if (graph) {
    graph.1()
    plot(t$tree, hang = -1, main = "",xlab = "", sub = "")
    if(bTest==TRUE) 
      rect <- rect.hclust(t$tree, h = t$tree$height[length(t$tree$height)], border = ordColo)
    
    else rect <- rect.hclust(t$tree, h = y, border = ordColo)
    clust <- NULL
    

    for (j in 1:nb.clust) clust <- c(clust, rep(j, length(rect[[j]])))
    
    clust <- as.factor(clust)
    belong <- cbind.data.frame(t$tree$order, clust)
    belong <- belong[do.call("order", belong), ]
    clust <- as.factor(belong$clust)
    if (nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
      layout(matrix(nrow = 1, ncol = 1, 1), respect = TRUE)
  }



  # Ver si aux no es necesario
  aux <- names(clust) <- rownames(X)
  # names(clust) <- aux
  clust <- as.factor(clust)
  X <- cbind.data.frame(X, clust)
  
  #  data.clust <- cbind.data.frame(res.sauv$call$Xtot[rownames(t$res$call$X),], clust)
  #  data.clust <- data.clust[rownames(res.sauv$row$coord), ]
  # return(data.clust)  # <-   PROBLEMA   Contiene vocabulario + variable cuantitativa + cluster
  # <----------- Xtot contiene las contextuales cuantitativas, lo cambio por S
  ###########   REVISAR EN LexHCca
  
  data.clust <- cbind.data.frame(res.sauv$call$X[rownames(t$res$call$X),], "_clust"=clust)
  data.clust <- data.clust[rownames(res.sauv$row$coord), ]
  #  return(data.clust)  # <-   YA NO PROBLEMA   Contiene vocabulario + cluster
  # 
  
  ################################################################
  # Checking description marg.doc before or after
  ################################################################
  doc.w.after <- apply(x$call$X,1,sum)    # Initial weighting for documents after selection
  nm <- names(doc.w.after)
  
  # doc.w.before before word's selection
  doc.w.before <- object$rowINIT[which(rownames(object$rowINIT) %in% nm) ,]
  names(doc.w.before) <- nm
  
  
  #@ Si marg.doc = before o no (after)
  
  # doc.w.before before word's selection
  doc.w.before <- x$rowINIT[which(rownames(x$rowINIT) %in% nm) ,]
  names(doc.w.before) <- nm
  
  if(marg.doc=="before") doc.w <- doc.w.before else doc.w <- doc.w.after # before word's selection
  doc.w <- as.data.frame(doc.w)
  
  
  #  stop("insertar weight en la siguiente")
  df.coord.cla<- data.frame(X[,c(1:(ncol(X)-1))], weight=doc.w[rownames(X),], cluster=X[,"clust"])
  

  
  #---- 12. Computing centers. Building RDO object  ----			
  
  RDO <-   ssquares(df.coord.cla[,c(1:(ncol(df.coord.cla)-2))],  df.coord.cla[,"weight"], df.coord.cla[,"cluster"])
  

  #  desc.var <- descfreq(data.clust[, -which(sapply(data.clust,is.factor))], 
  #                       data.clust[, ncol(data.clust)], proba = proba)

if(description) {
  desc.var <- descfreq(data.clust[,1:(ncol(data.clust)-1)], 
                       data.clust[, ncol(data.clust)], proba = proba)
  ### Modificar en el siguientes las salidas de los objetos:
  # Link between the cluster and the quantitative variables
  # Description of each cluster by quantitative variables
#   desc.axe <- FactoMineR::catdes(X, num.var=ncol(X), proba = proba)
} # End description

  
  tabInd <- cbind.data.frame(res.sauv$row$coord, data.clust[rownames(res.sauv$row$coord), 
                                                            ncol(data.clust)])
  colnames(tabInd)[ncol(tabInd)] <- "Cluster"
  
  list.centers <- by(tabInd[, -ncol(tabInd), drop = FALSE], 
                     tabInd[, ncol(tabInd)], colMeans)
  
  centers <- matrix(unlist(list.centers), ncol = ncol(tabInd) - 1, byrow = TRUE)
  colnames(centers) <- colnames(tabInd)[-ncol(tabInd)]
  cluster <- tabInd[, ncol(tabInd), drop = FALSE]
  
if(description) {
  para <- by(tabInd, cluster, simplify = FALSE, select, default.size = nb.par, 
             method = metric, coord.centers = centers)
  
  dist <- by(tabInd, cluster, simplify = FALSE, distinctivness, 
             default.size = nb.par, method = metric, coord.centers = centers)
  desc.doc <- list(para = para, dist = dist)
} # End description  

  colnames(data.clust)[ncol(data.clust)] <- "Clust_"  
  colnames(df.coord.cla)[ncol(df.coord.cla)] <- "Clust_"  
  
if(description ==FALSE)   desc.var <- desc.axe <- desc.doc <- NULL
  call <- list(t = t, min = min, max = max, X = X, vec=FALSE, call = match.call())
  description <- list(des.word = desc.var, des.doc = desc.doc)
  res.LexCHCca <- list(data.clust = data.clust, coord.clust= df.coord.cla, centers=centers, description=description,  
                        call = call, dendro=t$tree, phases=t$tree$clust, sum.squares=RDO)
  class(res.LexCHCca) <- "LexCHCca"

  if (graph) par(mar = old.mar)
  
  return(res.LexCHCca)
 
}
