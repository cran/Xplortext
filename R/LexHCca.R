#' @export
LexHCca<-function (object, nb.clust=0, consol=TRUE, iter.max=10, min=3, 
    max=NULL, order=TRUE, nb.par=5, edit.par=FALSE, graph=TRUE, proba=0.05, ...) 
{


if (!inherits(object, "LexCA")) stop("object should be LexCA class")
if (min <=1) min<-2
cluster.CA = "rows"
metric = "euclidean"
method = "ward"
kk = Inf

graph.scale = "inertia"
res.hcpc<-HCPC(object, nb.clust=nb.clust, consol=consol, iter.max=iter.max, min=min,max=max , metric=metric, 
               method=method , order=order, graph.scale=graph.scale, nb.par=nb.par, graph=graph, proba=proba, 
               cluster.CA=cluster.CA, kk=kk ) 
names(res.hcpc)[which(names(res.hcpc)=="desc.ind")]<-"desc.doc"
names(res.hcpc)[which(names(res.hcpc)=="desc.var")]<-"desc.wordvar"
res.hcpc$call$call<-match.call()

#      Computing the number of documents in each cluster
nb.clust<-nlevels(as.factor(res.hcpc$data.clust$clust))
count.cluster<-matrix(nrow=nb.clust,ncol=1)
partition<-as.factor(res.hcpc$data.clust$clust)
for (nc in 1:nb.clust) count.cluster[nc,1]<-length(which(partition==nc))
row.names(count.cluster)<-paste("cluster",1:nb.clust,sep="-")
colnames(count.cluster)<-c("counts")

#      Content of the clusters
data.clust<-res.hcpc$data.clust
compo.cluster<-list()
for (iclus in 1:nb.clust) compo.cluster[[iclus]]<- rownames(data.clust)[data.clust$clust == iclus]
names(compo.cluster) <- row.names(count.cluster)



#      Description of the clusters by the axes
   if ("quanti" %in% names(res.hcpc$desc.axes))
    {
    for(ilev in 1:nb.clust)
       {
       if(!is.null(res.hcpc$desc.axes$quanti[[ilev]]))  {
          colnames(res.hcpc$desc.axes$quanti[[ilev]])[colnames(res.hcpc$desc.axes$quanti[[ilev]])=="Mean in category"]<-"Mean in cluster"
          colnames(res.hcpc$desc.axes$quanti[[ilev]])[colnames(res.hcpc$desc.axes$quanti[[ilev]])=="sd in category"]<-"sd in cluster"
                                                       }
       }
     }

#### description of the clusters by the characteristic documents
if (  (is.null(object$var.agg)) & (edit.par == TRUE) )
{
#### creation of the structure corpus
 var.text <- object$info$var.text[[1]]
 str.base <- object$info$base[[1]]
 str.envir <- object$info$menvir[[1]]
 base <- get(str.base, envir=str.envir)
 corpus <- base[, var.text[1]]

 if(length(var.text) > 1) {					
   for (i in 2:length(var.text)){					
    corpus <- paste(corpus, base[, var.text[i]], sep = ".")
                          }}

    corpus <- data.frame(corpus)
    rownames(corpus) <- rownames(base)
    corpus[rownames(object$call$X),]

    nb.clus<-nlevels(as.factor(res.hcpc$data.clust[,"clust"]))
    vclus<-c(1:nb.clus)

    lispara <- vector(mode="list")
    for(iclus in 1:nb.clus)
    {
        doctot <- which(res.hcpc$data.clust[rownames(object$call$X),"clust",drop=FALSE] == iclus)
        ntdoc <- min(nb.par, length(doctot))
        lispara[[iclus]]<-data.frame()
        for (i in 1:ntdoc)
        {
      lispara[[iclus]][i,1] <- rownames(as.data.frame(res.hcpc$desc.doc$para[[iclus]]) )[[i]]
      lispara[[iclus]][i,2] <- res.hcpc$desc.doc$para[[iclus]][[i]]
      lispara[[iclus]][i,3] <- as.character(corpus[rownames
                            (as.data.frame(res.hcpc$desc.doc$para[[iclus]]) )[[i]],1])
        }
       colnames(lispara[[iclus]]) <- c("DOCUMENT", "CRITERION", "---------------------TEXT---------------------")
   }
   names(lispara) <-paste("cluster",1:nb.clus,sep="_")

    lisdist <- vector(mode="list",length=nb.clus)
    for(iclus in 1:nb.clus)
    {
        doctot <- which(res.hcpc$data.clust[rownames(object$call$X),"clust",drop=FALSE] == iclus)
        ntdoc <- min(nb.par, length(doctot))
        lisdist[[iclus]]<-data.frame()
        for (i in 1:ntdoc)
       {
   lisdist[[iclus]][i,1] <- rownames(as.data.frame(res.hcpc$desc.doc$dist[[iclus]]) )[[i]]
   lisdist[[iclus]][i,2] <- res.hcpc$desc.doc$dist[[iclus]][[i]]
   lisdist[[iclus]][i,3] <- as.character(corpus[rownames
                            (as.data.frame(res.hcpc$desc.doc$dist[[iclus]]) )[[i]],1])
       }
colnames(lisdist[[iclus]]) <- c("DOCUMENT", "CRITERION", "---------------------TEXT---------------------")
    }
 names(lisdist) <- paste("cluster",1:nb.clus,sep="_")
}

if (is.null(object$var.agg) )  {
          listcom <- list(count.cluster,compo.cluster)
          names(listcom)<-c("clust.count","clust.content")
          if (edit.par ==  TRUE)   {
                listcom$docspara<-lispara
                listcom$docsdist<-lisdist 
                                   }
          res<-c(res.hcpc,listcom)
                               }
if (!is.null(object$var.agg)) { 
           listcom <- list(count.cluster,compo.cluster)
           names(listcom)<-c("clust.count","clust.content")
           res<-c(res.hcpc,listcom)
                              }
    class(res) = c("LexHCPC", "HCPC")
    return(res)
}



