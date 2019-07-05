#' @importFrom graphics abline layout legend locator par plot text
#' @importFrom grDevices palette
#' @export


######### Mínimo en 3
######### Ver para y dist con marginales nuevas.. y descripciones
LexCHCca = function (object, nb.clust=0, min=3, 
    max=NULL, nb.par=5, graph=TRUE, proba=0.05) 
{
#### Esta se repite después
if (!inherits(object,"LexCA")) stop("Object should be LexCA class")

###  
hcclust= function(X)  
  {
   ndoc<-nrow(X)
   d<-dist(X)
# Mdist
   Mdist<-as.matrix(d)
# Msim
   maxd<-max(Mdist)+1e-10
   Msim<-(maxd-Mdist)
#Mcont
   Mcont<-rep(0,nrow(Msim)*ncol(Msim))
   attr(Mcont,"dim") <- c(nrow(Msim), ncol(Msim))
   Mcont[1,2]<-1
   for (i in 2:(nrow(Mcont)-1))
     {
	  Mcont[i,i+1]<-1
	  Mcont[i,i-1]<-1
	}
    Mcont[nrow(Mcont),nrow(Mcont)-1]<-1
    rownames(Mcont)<-rownames(Msim)
    colnames(Mcont)<-colnames(Msim)



###   Introduction of the constriction into the similarity and distance matrices
    MSimCont<-Msim*Mcont
    MDistCont<-Mdist*Mcont

###    Groups existing at Phasis 0
    groups<-list()                           
    for (i in 1:nrow(Msim)) groups[[i]]<-i
###    Critagreg will file the successive values of the agregation criterium
    Critagreg<-numeric()                     # valor of the distance in the (n-1) fusions of nodes
###  clust will archiv the successive pairs of nodes that are merged
    clust<-list()
###    Building the hierarchy: (n-1) intermediate nodes have to be created
    i<-1
    indice<-nrow(Mdist)-1

#####################################################################################################################################################
# creating the n-1 nodes 
     while(indice>0) {	                 
#        Indentifying the maximum similarity
       maxsim<-max(MSimCont)
       posmaxsim<-which(MSimCont==maxsim)			
       if (posmaxsim[1]%%nrow(MSimCont)==0) {
	   fila<-posmaxsim[1]%/%nrow(MSimCont)
	   col<-nrow(MSimCont)
       }else{
	   fila<-posmaxsim[1]%/%nrow(MSimCont)+1
	   col<-posmaxsim[1]%%nrow(MSimCont)   }

	 maxfc<-max(fila,col)
       minfc<-min(fila,col)     # maxfc and minfc correspond to the row and col that are grouped
       Critagreg[i]<-MDistCont[fila,col]
	 clust[[i]]<-vector(mode="list",length=2)
       clust[[i]][[1]]<-groups[[minfc]]
       clust[[i]][[2]]<-groups[[maxfc]]

       rownames(Msim)[minfc]<-colnames(Msim)[minfc]<-rownames(Mdist)[minfc]<-colnames(Mdist)[minfc]<-rownames(Mcont)[minfc]<-
       colnames(Mcont)[minfc]<-paste(rownames(Msim)[minfc],"-",rownames(Msim)[maxfc])

# Actualization of Msim, Mdist and Mcont
       if (minfc!=1) {
         Msim[minfc,minfc-1]<-Msim[minfc-1,minfc]<-0.5*Msim[minfc,minfc-1]+0.5*Msim[maxfc,minfc-1]-0.5*abs(Msim[minfc,minfc-1]-      Msim[maxfc,minfc-1])
         Mdist[minfc,minfc-1]<-Mdist[minfc-1,minfc]<-0.5*Mdist[minfc,minfc-1]+0.5*Mdist[maxfc,minfc-1]+0.5*abs(Mdist[minfc,minfc-1]-      Mdist[maxfc,minfc-1])
         Mcont[minfc-1,minfc]<-Mcont[minfc,minfc-1]<-1
                     }
       if (maxfc!=nrow(MSimCont)) {
         Msim[maxfc+1,minfc]<-Msim[minfc,maxfc+1]<-0.5*Msim[minfc,maxfc+1]+0.5*Msim[maxfc,maxfc+1]-0.5*abs(Msim[minfc,maxfc       +1]-Msim[maxfc,maxfc+1])
         Mdist[maxfc+1,minfc]<-Mdist[minfc,maxfc+1]<-0.5*Mdist[minfc,maxfc+1]+0.5*Mdist[maxfc,maxfc+1]+0.5*abs(Mdist       [minfc,maxfc+1]-Mdist[maxfc,maxfc+1])
         Mcont[maxfc+1,minfc]<-Mcont[minfc,maxfc+1]<-1
                                  }

#		
	 groups[[minfc]]<-c(groups[[minfc]],groups[[maxfc]])
	 groups<-groups[-maxfc]
	 Msim<-Msim[-maxfc,-maxfc]
	 Mdist<-Mdist[-maxfc,-maxfc]
	 Mcont<-Mcont[-maxfc,-maxfc]	
	 MSimCont<-Msim*Mcont
	 MDistCont<-Mdist*Mcont
	 indice<-indice-1
	 i<-i+1
    }   # endwhile # endwhile 



 
# ############################################################################################################################################
# The hierarchy is completed, the "type hclust" has to be built
#
#    hccc<-list()
    height<-Critagreg
    order<-as.integer(1:nrow(X))


#### begin of grups_blocs
    grups_blocs<-list()
    grups_blocs[[1]]<-rep(0,nrow(X))
    for (i in 1:(length(clust)-1)) {
      grups_blocs[[i+1]]<-grups_blocs[[i]]
	grups_blocs[[i+1]][c(clust[[i]][[1]],clust[[i]][[2]])]<-i
                                   }
     merge<-matrix(nrow=(ndoc-1),ncol=2,0)


### begin of the loop of type for that builds merge
    for(i in 1:(ndoc-1)) {
      if (length(clust[[i]][[1]])==1&length(clust[[i]][[2]])==1)  {
	  merge[i,1]<-(-clust[[i]][[1]])
	  merge[i,2]<-(-clust[[i]][[2]])
      }else{
	  if (length(clust[[i]][[1]])==1) {
	     merge[i,1]<-(-clust[[i]][[1]])
	  }else{
	     merge[i,1]<-grups_blocs[[i]][clust[[i]][[1]][1]]
	     }		
	   if (length(clust[[i]][[2]])==1) {
	     merge[i,2]<-(-clust[[i]][[2]])
	   }else{
	   merge[i,2]<-grups_blocs[[i]][clust[[i]][[2]][1]]
	                                    }
                                                                 } 	  
                                      }
### end for

    hcccall<-call("hcclust",match.call()$res)		
    hcc<-structure(list(merge = merge,height=height,
		       order=order, 
                   labels = attr(d, "Labels"), method = "complete", 
                   call = hcccall, dist.method = "euclidean",clust=clust),class="hclust")

    return(hcc)
}


######
    auto.cut.constree = function(res, min, max) 
    {
        X = as.data.frame(res$ind$coord)
        hcc <- hcclust(X)
        height=hcc$height
        inv.height<-rev(height)
        quot<-inv.height[(min-1):(max-1)]/inv.height[(min):(max)]
        nb.clust = which.max(quot) + min - 1
        return(list(res = res, tree = hcc, nb.clust = nb.clust,height=height,quot = quot))
    }


### function select
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
            distance <- distance[(nrow(Y) + 1):nrow(distance), 
                -((nrow(Y) + 1):ncol(distance))]
            distance <- sort(distance[clust, ], decreasing = FALSE)
        }
        if (length(distance) > default.size) 
            distance <- distance[1:default.size]
        else distance <- distance
    }
### function distinctivness
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

### main part of the function
   if (!inherits(object, "LexCA"))  stop ("object should be LexCA class")
   res.ca<- object
   res.sauv <- res.ca
   metric="euclidean"
   method="complete"
   if (min <= 1) min<-2



### PCA equivalent to CA
      aux <- res.ca$eig
      weight=row.w = res.ca$call$marge.row * sum(res.ca$call$X)
      res <- PCA(res.ca$row$coord, scale.unit = FALSE, ncp = Inf, graph = FALSE, 
      row.w = row.w)
      res$eig <- aux

###
   if (is.null(max)) 
        max <- min(10, round(nrow(res$ind$coord)/2))
    max <- min(max, nrow(res$ind$coord) - 1)
###
   t <- auto.cut.constree(res=res,min = min, max = max) 
   nb.ind <- nrow(t$res$ind$coord)


   auto.haut <- ((t$tree$height[length(t$tree$height) - 
            t$nb.clust + 2]) + (t$tree$height[length(t$tree$height) - 
            t$nb.clust + 1]))/2


   if (graph) {
            if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
                dev.new()
            old.mar <- par()$mar
            par(mar = c(0.5, 2, 0.75, 0))
            lay = matrix(ncol = 5, nrow = 5, c(2, 4, 4, 4, 4, 
                2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 
                1, 3, 3, 3, 3))
            layout(lay, respect = TRUE)
            barplot(rev(t$tree$height)[1:max(15, max)], col = c(rep("black", 
                t$nb.clust - 1), rep("grey", max(max, 15) - t$nb.clust + 
                1)), rep(0.1, max(max, 15)), space = 0.9)
            plot(x = 1, xlab = "", ylab = "", main = "", col = "white", 
                axes = FALSE)
            text(1, 1, "Constrained Hierarchical Clustering", cex = 2)
            plot(x = 1, xlab = "", ylab = "", main = "", col = "white", 
                axes = FALSE)
            legend("top", "Aggreg. Crit", box.lty = NULL, cex = 1)
        } else {
            if (nb.clust == 0 | nb.clust == 1) 
                nb.clust <- -1
        }

        if ((nb.clust == 0) | (nb.clust == 1)) {
            plot(t$tree, hang = -1, main = "Click to cut the tree", 
                xlab = "", sub = "")
            abline(h = auto.haut, col = "black", lwd = 3)
            coupe <- locator(n = 1)
            while (coupe$y < min(t$tree$height)) {
                cat("No class \n")
                coupe <- locator(n = 1)
            }
            y <- coupe$y
        } else {
            if (graph) 
                plot(t$tree, hang = -1, main = "Constrained Hierarchical Classification", 
                  xlab = "", sub = "")
            if (nb.clust < 0) 
                y = auto.haut
            else y = (t$tree$height[length(t$tree$height) - nb.clust + 
                2] + t$tree$height[length(t$tree$height) - nb.clust + 
                1])/2
        }

t$tree$height[t$tree$height < 0.0000000000001]<- 0
    clust <- cutree(as.hclust(t$tree), h = y)
    nb.clust <- max(clust)

    X = as.data.frame(t$res$ind$coord)
    ordColo = unique(clust[t$tree$order])
    if (graph) {
        rect <- rect.hclust(t$tree, h = y, border = ordColo)
        clust <- NULL
        for (j in 1:nb.clust) clust <- c(clust, rep(j, length(rect[[j]])))
        clust <- as.factor(clust)
        belong <- cbind.data.frame(t$tree$order, clust)
        belong <- belong[do.call("order", belong), ]
        clust <- as.factor(belong$clust)
        if (nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
            layout(matrix(nrow = 1, ncol = 1, 1), respect = TRUE)
    }
    aux <- names(clust) <- rownames(X)
    names(clust) <- aux
    clust <- as.factor(clust)

    X <- cbind.data.frame(X, clust)
    data.clust <- cbind.data.frame(res.sauv$call$Xtot[rownames(t$res$call$X), 
                  ], clust)
    data.clust <- data.clust[rownames(res.sauv$row$coord), ]
    desc.var <- descfreq(data.clust[, -which(sapply(data.clust, 
        is.factor))], data.clust[, ncol(data.clust)], proba = proba)
    desc.axe <- catdes(X, ncol(X), proba = proba)
    tabInd <- cbind.data.frame(res.sauv$row$coord, data.clust[rownames(res.sauv$row$coord), 
            ncol(data.clust)])
    colnames(tabInd)[ncol(tabInd)] = "Cluster"


    list.centers <- by(tabInd[, -ncol(tabInd), drop = FALSE], 
        tabInd[, ncol(tabInd)], colMeans)
    centers <- matrix(unlist(list.centers), ncol = ncol(tabInd) - 
        1, byrow = TRUE)
    colnames(centers) = colnames(tabInd)[-ncol(tabInd)]
    cluster <- tabInd[, ncol(tabInd), drop = FALSE]
    para <- by(tabInd, cluster, simplify = FALSE, select, default.size = nb.par, 
        method = metric, coord.centers = centers)
    dist <- by(tabInd, cluster, simplify = FALSE, distinctivness, 
        default.size = nb.par, method = metric, coord.centers = centers)
    desc.doc <- list(para = para, dist = dist)

   call <- list(t = t, min = min, max = max, X = X, vec=FALSE, call = match.call())
   res.LexCHCca <- list(data.clust = data.clust, desc.word = desc.var, 
        desc.axes = desc.axe, call = call, desc.doc = desc.doc, dendro=t$tree$clust)

    res = x = res.LexCHCca
    X = res$call$X
    max = res$call$max
    min = res$call$min
    max.plot = max(res$call$max, 15)
    nb.clust = length(levels(X$clust))
    levs = levels(X$clust)
    label = "ind"
    new.plot = TRUE
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
    axes=c(1,2)
    if(graph) plot.PCA(res2, title = title, habillage = ncol(Y), 
                  cex = 0.8, axes = axes, new.plot = new.plot, 
                  palette = palette())
 
    class(res.LexCHCca) = "LexCHCca"
    if (graph) par(mar = old.mar)
    return(res.LexCHCca)
}

