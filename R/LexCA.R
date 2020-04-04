#' @import FactoMineR
#' @importFrom methods hasArg
#' @export
LexCA <- function(object,ncp=5, context.sup="ALL", doc.sup=NULL, word.sup=NULL, 
 segment=FALSE, graph=TRUE, axes=c(1, 2), lmd=3, lmw=3)
{
if(is.null(object)) stop("Missing argument for object")
if (!inherits(object,"TextData")) stop("Object should be TextData class")
 options(stringsAsFactors = FALSE)
 
# Starting CA_New  
 CA_New <- function (X, ncp = 5, row.sup = NULL, col.sup = NULL, quanti.sup = NULL, 
                     quali.sup = NULL, graph = TRUE, axes = c(1, 2), row.w = NULL, 
                     excl = NULL) 
 {
   
   fct.eta2 <- function(vec, x, weights) {
     VB <- function(xx) {
       return(sum((colSums((tt * xx) * weights)^2)/ni))
     }
     tt <- tab.disjonctif(vec)
     ni <- colSums(tt * weights)
     unlist(lapply(as.data.frame(x), VB))/colSums(x * x *  weights)
   }
   
   
   if (is.table(X)) X <- matrix(as.vector(X), nrow(X), dimnames = dimnames(X))
   
   if (is.null(rownames(X))) rownames(X) <- 1:nrow(X)
   
   if (is.null(colnames(X))) colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "V")
   
   X <- as.data.frame(X)
   is.quali <- which(!unlist(lapply(X, is.numeric)))
   
   X[, is.quali] <- lapply(X[, is.quali, drop = FALSE], as.factor)
   for (i in is.quali) X[, i] = as.factor(X[, i])
   X <- droplevels(X)
   Xtot <- X
   
   
   if (any(!sapply(X, is.numeric))) {
     auxi = NULL
     for (j in (1:ncol(X))[!((1:ncol(X)) %in% quali.sup)]) if (!is.numeric(X[, j])) 
       auxi = c(auxi, colnames(X)[j])
     if (!is.null(auxi)) 
       stop(paste("\nThe following variables are not quantitative: ", auxi))
   }
   if (!inherits(X, "data.frame")) stop("X is not a data.frame")
   if (!is.null(row.sup)) X <- as.data.frame(X[-row.sup, ])
   if ((!is.null(col.sup)) || (!is.null(quanti.sup)) || (!is.null(quali.sup))) 
     X <- as.data.frame(X[, -c(col.sup, quanti.sup, quali.sup)])
   
   if (any(apply(X, 1, sum) == 0)) {
     warning(paste0("The rows ", paste(rownames(X)[which(apply(X, 1, sum) == 0)], collapse = ", "),
                    " sum at 0. They were suppressed from the analysis"))
     X <- X[-which(apply(X, 1, sum) == 0), , drop = FALSE]
   }
   
   
   if (any(apply(X, 2, sum) == 0)) {
     warning(paste0("The columns ", paste(colnames(X)[which(apply(X, 
                                                                  2, sum) == 0)], collapse = ", "), " sum at 0. They were suppressed from the analysis"))
     X <- X[, -which(apply(X, 2, sum) == 0), drop = FALSE]
   }
   
   
   if (is.null(row.w)) row.w = rep(1, nrow(X))
   row.w.init <- row.w
   
   if (length(row.w) != nrow(X)) stop("length of vector row.w should be the number of active rows")
   total <- sum(X * row.w)
   
   F <- as.matrix(X) * (row.w/total)
   marge.col <- colSums(F)
   marge.row <- rowSums(F)
   ncp <- min(ncp, (nrow(X) - 1), (ncol(X) - 1))
   Tc <- t(t(F/marge.row)/marge.col) - 1
   if (!is.null(excl)) marge.col[excl] <- 1e-15
   tmp <- svd.triplet(Tc, row.w = marge.row, col.w = marge.col, ncp = ncp)
   
   if (!is.null(excl))  marge.col[excl] <- 0
   eig <- tmp$vs^2
   vp <- matrix(NA, length(eig), 3)
   rownames(vp) <- paste("dim", 1:length(eig))
   colnames(vp) <- c("eigenvalue", "percentage of variance", 
                     "cumulative percentage of variance")
   vp[, "eigenvalue"] <- eig
   vp[, "percentage of variance"] <- (eig/sum(eig)) * 100
   vp[, "cumulative percentage of variance"] <- cumsum(vp[,"percentage of variance"])
   V <- tmp$V
   U <- tmp$U
   
   eig <- eig[1:ncol(U)]
   coord.col <- t(t(V) * sqrt(eig))
   coord.row <- t(t(U) * sqrt(eig))
   dist2.col <- colSums(Tc^2 * marge.row)
   contrib.col <- t(t(coord.col^2 * marge.col)/eig)
   cos2.col <- coord.col^2/dist2.col
   colnames(coord.col) <- colnames(contrib.col) <- colnames(cos2.col) <- paste("Dim",1:length(eig))
   rownames(coord.col) <- rownames(contrib.col) <- rownames(cos2.col) <- attributes(X)$names
   dist2.row <- rowSums(t(t(Tc^2) * marge.col))
   contrib.row <- t(t(coord.row^2 * marge.row)/eig)
   cos2.row <- coord.row^2/dist2.row
   colnames(coord.row) <- colnames(contrib.row) <- colnames(cos2.row) <- paste("Dim",1:length(eig))
   rownames(coord.row) <- rownames(contrib.row) <- rownames(cos2.row) <- attributes(X)$row.names
   
   
   inertia.row = marge.row * dist2.row
   inertia.col = marge.col * dist2.col
   names(inertia.col) <- attributes(coord.col)$row.names
   names(inertia.row) <- attributes(coord.row)$row.names
   res.call <- list(X = X, marge.col = marge.col, marge.row = marge.row, 
                    ncp = ncp, row.w = row.w, excl = excl, call = match.call(), 
                    Xtot = Xtot, N = sum(row.w * rowSums(X)))
   
   res.col <- list(coord = as.matrix(coord.col[, 1:ncp]), contrib = as.matrix(contrib.col[, 1:ncp] * 100), 
                   cos2 = as.matrix(cos2.col[, 1:ncp]), inertia = inertia.col)
   res.row <- list(coord = coord.row[, 1:ncp], contrib = contrib.row[, 1:ncp] * 100,
                   cos2 = cos2.row[, 1:ncp], inertia = inertia.row)
   res <- list(eig = vp[1:min(nrow(X) - 1, ncol(X) - 1), , drop = FALSE], 
               call = res.call, row = res.row, col = res.col, svd = tmp)
   
   
   
   
   if (!is.null(row.sup)) {
     X.row.sup <- as.data.frame(Xtot[row.sup, ])
     
     if ((!is.null(col.sup)) || (!is.null(quanti.sup)) || 
         (!is.null(quali.sup))) 
       X.row.sup <- as.data.frame(X.row.sup[, -c(col.sup, quanti.sup, quali.sup)])
     somme.row <- rowSums(X.row.sup)
     X.row.sup <- X.row.sup/somme.row
     coord.row.sup <- crossprod(t(as.matrix(X.row.sup)), V)
     dist2.row <- rowSums(t((t(X.row.sup) - marge.col)^2/marge.col))
     cos2.row.sup <- coord.row.sup^2/dist2.row
     coord.row.sup <- coord.row.sup[, 1:ncp, drop = FALSE]
     cos2.row.sup <- cos2.row.sup[, 1:ncp, drop = FALSE]
     colnames(coord.row.sup) <- colnames(cos2.row.sup) <- paste("Dim", 1:ncp)
     rownames(coord.row.sup) <- rownames(cos2.row.sup) <- rownames(X.row.sup)
     res.row.sup <- list(coord = coord.row.sup, cos2 = cos2.row.sup)
     res$row.sup <- res.row.sup
     res$call$row.sup <- row.sup
   }
   
   if (!is.null(col.sup)) {
     X.col.sup <- as.data.frame(Xtot[, col.sup])
     if (!is.null(row.sup)) 
       X.col.sup <- as.data.frame(X.col.sup[-row.sup, ])
     X.col.sup <- X.col.sup * row.w
     colnames(X.col.sup) <- colnames(Xtot)[col.sup]
     somme.col <- colSums(X.col.sup)
     X.col.sup <- t(t(X.col.sup)/somme.col)
     coord.col.sup <- crossprod(as.matrix(X.col.sup), U)
     dist2.col <- colSums((X.col.sup - marge.row)^2/marge.row)   
     coord.col.sup <- as.matrix(coord.col.sup[, 1:ncp, drop = FALSE])
     cos2.col.sup <- coord.col.sup^2/dist2.col
     cos2.col.sup <- cos2.col.sup[, 1:ncp, drop = FALSE]
     colnames(coord.col.sup) <- colnames(cos2.col.sup) <- paste("Dim",1:ncp)   
     rownames(coord.col.sup) <- rownames(cos2.col.sup) <- colnames(X.col.sup)
     res.col.sup <- list(coord = coord.col.sup, cos2 = cos2.col.sup)
     res$col.sup <- res.col.sup
     res$call$col.sup <- col.sup
   }
   
   
   if (!is.null(quanti.sup)) {
     coord.quanti.sup <- matrix(NA, length(quanti.sup), ncp)
     
     if (is.null(row.sup)) 
       coord.quanti.sup <- cov.wt(cbind.data.frame(res$row$coord, 
                                                   Xtot[, quanti.sup, drop = FALSE]), cor = TRUE, 
                                  wt = marge.row, method = "ML")$cor[-(1:ncp), 1:ncp, drop = FALSE]
     else coord.quanti.sup <- cov.wt(cbind.data.frame(res$row$coord, 
                                                      Xtot[-row.sup, quanti.sup, drop = FALSE]), wt = marge.row, 
                                     cor = TRUE, method = "ML")$cor[-(1:ncp), 1:ncp, drop = FALSE]
     dimnames(coord.quanti.sup) <- list(colnames(Xtot)[quanti.sup], paste("Dim", 1:ncp, sep = "."))   
     res$quanti.sup$coord <- coord.quanti.sup
     res$quanti.sup$cos2 <- coord.quanti.sup^2
     res$call$quanti.sup <- quanti.sup
   }
   
   
   if (!is.null(quali.sup)) {
     if (!is.null(row.sup)) 
       X.del <- as.data.frame(Xtot[-row.sup, -c(col.sup,quanti.sup, quali.sup)])
     else X.del <- Xtot[, -c(col.sup, quanti.sup, quali.sup)]
     
     X.quali.sup <- NULL
     Xtot2 <- Xtot
     if (!is.null(row.sup))  Xtot2 <- Xtot[-row.sup,] 
     
     for (j in 1:length(quali.sup)) {
       Xtot2[,quali.sup[j]] <- droplevels(Xtot2[,quali.sup[j]] , reorder=FALSE)  
       X.quali.sup <- rbind(X.quali.sup, matrix(unlist(by(X.del, 
                                                          Xtot2[, quali.sup[j]], colSums)), ncol = ncol(X.del), byrow = T))
     }
     
     somme.quali <- rowSums(X.quali.sup)
     X.quali.sup <- X.quali.sup/somme.quali
     coord.quali.sup <- crossprod(t(as.matrix(X.quali.sup)), V)
     dist2.quali <- rowSums(t((t(X.quali.sup) - marge.col)^2/marge.col))
     cos2.quali.sup <- coord.quali.sup^2/dist2.quali
     coord.quali.sup <- coord.quali.sup[, 1:ncp, drop = FALSE]
     cos2.quali.sup <- cos2.quali.sup[, 1:ncp, drop = FALSE]
     
     rownames(coord.quali.sup) <- rownames(cos2.quali.sup) <- paste(rep(colnames(Xtot2)[quali.sup], 
                                                                        lapply(Xtot2[, quali.sup, drop = FALSE], nlevels)), 
                                                                    unlist(lapply(Xtot2[, quali.sup, drop = FALSE], levels)),sep = ".") 
     
     colnames(coord.quali.sup) <- colnames(cos2.quali.sup) <- paste("Dim", 1:ncp)
     res$quali.sup <- list(coord = coord.quali.sup, cos2 = cos2.quali.sup)
     
     Zqs <- tab.disjonctif(Xtot2[, quali.sup])
     
     Nj <- colSums(Zqs * row.w)
     Nj <- colSums(Zqs * marge.row) * total
     if (total > 1) 
       coef <- sqrt(Nj * ((total - 1)/(total - Nj)))
     else coef <- sqrt(Nj)
     res$quali.sup$v.test <- res$quali.sup$coord * coef
     
     eta2 = matrix(NA, length(quali.sup), ncp)
     eta2 <- sapply(as.data.frame(Xtot2[, quali.sup, drop = FALSE]), 
                    fct.eta2, res$row$coord, weights = marge.row)
     
     eta2 <- t(as.matrix(eta2, ncol = ncp))
     colnames(eta2) = paste("Dim", 1:ncp)
     rownames(eta2) = colnames(Xtot)[quali.sup]
     res$quali.sup$eta2 <- eta2
     res$call$quali.sup <- quali.sup
   }
   class(res) <- c("CA", "list")
   if (graph & (ncp > 1)) {
     print(plot(res, axes = axes))
     if (!is.null(quanti.sup)) 
       print(plot(res, choix = "quanti.sup", axes = axes, 
                  new.plot = TRUE))
   }
   return(res)
 }
 
 
##### Final CA_New 
 
# Functions -------------------
plotLexCA <- function()
{

# if(dev.interactive()) dev.new()
plot.LexCA(res, selDoc=NULL, selWord=NULL, eigen=TRUE, axes=axes, graph.type = "classic", new.plot=TRUE)


if(!is.null(res$quanti.sup$coord)) {
#  if(dev.interactive()) dev.new()
 plot.LexCA(res, selDoc=NULL, selWord=NULL, quanti.sup=rownames(res$quanti.sup$coord), eigen=FALSE, axes=axes, graph.type = "classic",
            new.plot=TRUE)
}

#if(dev.interactive()) dev.new()
  plot.LexCA(res, selDoc= paste("meta ", lmd), selWord= paste("meta ", lmw), eigen=FALSE, axes=axes, graph.type = "classic",
             new.plot=TRUE)

if(!is.null(rownames(res$segment$coord))& segment==TRUE){
#  if(dev.interactive()) dev.new()
  plot.LexCA(res, selDoc= NULL, selWord= NULL, selSeg="ALL", axes=axes, graph.type = "classic",
             new.plot=TRUE)}

if(!is.null(res$quali.sup)) {
# if(dev.interactive()) dev.new()
  plot.LexCA(res, selDoc= NULL, selWord= NULL, quali.sup="ALL", axes=axes, graph.type = "classic",
             new.plot=TRUE)
 }

if(!is.null(res$row.sup$coord)){
   ident <- identical(dimnames(res$quali.sup$coord)[1], dimnames(res$row.sup$coord)[1])
   if(ident==FALSE) {
#  if(dev.interactive()) dev.new()
   plot.LexCA(res, selDoc= NULL, selWord= NULL, selDocSup="ALL", axes=axes, graph.type = "classic",
              new.plot=TRUE)
 }
}

if(!is.null(res$col.sup$coord)){
   ident <- identical(dimnames(res$col.sup$coord)[1], dimnames(res$segment$coord)[1])
   if(ident==FALSE) {
#  if(dev.interactive()) dev.new()
   plot.LexCA(res, selDoc= NULL, selWord= NULL, selWordSup="ALL", selSeg=NULL, axes=axes, graph.type = "classic",
              new.plot=TRUE) }
}
}


# Put documents at the end
remdocs <- function(z,doc_sup){
  df <- data.frame(z[,, drop=FALSE])
  tmprow <- rownames(df)[doc_sup]
  rowqq <- data.frame(df[doc_sup,, drop=FALSE])	
  df <- rbind(df,rowqq)	
  df <- data.frame(df[-doc_sup,,drop=FALSE])	
  rownames(df)[(ndoc-ndocsup+1):ndoc] <- tmprow	
  z <- data.frame(df[,, drop=FALSE])
 return(z)}

# --------------- Check context quanti-quali
checkcontext <- function(context_, TDContext) {
 if (is.character(context_)) 	
        context_ <- which(colnames(TDContext) %in% context_)	
        context_ <- colnames(TDContext)[context_]	
        context_ <- context_[!is.na(context_)]	
        context_ <- which(colnames(TDContext) %in% context_)	
  if(length(context_ )==0) context_ <- NULL
 return(context_) }


# --------------- Keys function
 Keys <- function(res,lmd,lmw,naxes, axes) {			
 title <- "Metakeys & Dockeys"
 col.word <- "red"; col.doc<-"blue"; cex=1
 axx <- axes[1]
 axy <- axes[2]			
 Ccontr <- res$col$contrib			
 Ccoor <- res$col$coord			
 Rcontr <-res$row$contrib			
 Rcoor <- res$row$coord			
 eigen <- res$eig			
if(naxes>=min(nrow(Rcontr),nrow(Ccontr)))			
 naxes <- min((nrow(Rcontr)-1),(nrow(Ccontr)-1))			
# Computing the metakeys: words with contributions over lmw*average_contribution of words  		
Metakeys <-vector(mode="list",length=naxes)		
for (i in 1:naxes){		
	Metakeys[[i]]=vector(mode="list",length=2)}	

for(i in 1:naxes){		
	Metakeys[[i]][[1]]<-sort(Ccontr[which(Ccontr[,i]>lmw*mean(Ccontr[,i])&Ccoor[,i]>0),i],decreasing=TRUE)	
	if (length(Metakeys[[i]][[1]])==1) 	
          names(Metakeys[[i]][[1]])<-rownames(Ccontr)[which(Ccontr[,i]%in%sort(Ccontr[which(Ccontr[,i]>lmw*mean(Ccontr[,i])&Ccoor[,i]>0),i],decreasing=TRUE))]		
	Metakeys[[i]][[2]]<-sort(Ccontr[which(Ccontr[,i]>lmw*mean(Ccontr[,i])&Ccoor[,i]<0),i],decreasing=TRUE)	
	if (length(Metakeys[[i]][[2]])==1) 	
          names(Metakeys[[i]][[2]])<-rownames(Ccontr)[which(Ccontr[,i]%in%sort(Ccontr[which(Ccontr[,i]>lmw*mean(Ccontr[,i])&Ccoor[,i]<0),i],decreasing=TRUE))]		
}		
# Computing the Keydocs : documents/answers with contributions over lmd*average_contribution of documents 			
#lmd <- 3; lmw<-3


Keydocs<-vector(mode="list",length=naxes)			
for (i in 1:naxes){			
	Keydocs[[i]]=vector(mode="list",length=2)}		
for(i in 1:naxes){			
	Keydocs[[i]][[1]]<-sort(Rcontr[which(Rcontr[,i]>lmd*mean(Rcontr[,i])&Rcoor[,i]>0),i],decreasing=TRUE)		
	if (length(Keydocs[[i]][[1]])==1) 		
       names(Keydocs[[i]][[1]])<-rownames(Rcontr)[which(Rcontr[,i]%in%sort(Rcontr[which(Rcontr[,i]>lmd*mean(Rcontr[,i])&Rcoor[,i]>0),i],decreasing=TRUE))]			
	Keydocs[[i]][[2]]<-sort(Rcontr[which(Rcontr[,i]>lmd*mean(Rcontr[,i])&Rcoor[,i]<0),i],decreasing=TRUE)		
	if (length(Keydocs[[i]][[2]])==1) 		
        names(Keydocs[[i]][[2]])<-rownames(Rcontr)[which(Rcontr[,i]%in%sort(Rcontr[which(Rcontr[,i]>lmd*mean(Rcontr[,i])&Rcoor[,i]<0),i],decreasing=TRUE))]			
}			
# metakeys and Keydocs				
       Metakeys.Keydocs <-vector(mode="list",naxes)				
       names(Metakeys.Keydocs)<-paste("DIM",1:naxes,sep="")				
	for (i in 1:naxes){			
                X <-list(names(Metakeys[[i]][[1]]),names(Keydocs[[i]][[1]]),				
                      names(Metakeys[[i]][[2]]),names(Keydocs[[i]][[2]]))				
            names(X)<-c("Metakeys+","Keydocs+","Metakeys-","Keydocs-")				
           Metakeys.Keydocs[[i]]<-X       }				
 #Dimension of the words				
    mkeys<-numeric()				
	for (i in 1:naxes){			
		for (j in 1:2){		
			mkeys<-c(mkeys,names(Metakeys[[i]][[j]]))	}}
mkeys<-as.factor(mkeys)	
levels(mkeys)	
matmkeys<-t(rbind(summary(mkeys,maxsum=10000),rep(naxes,length(levels(mkeys))),round(summary(mkeys,maxsum=10000)*100/naxes,2)))	
colnames(matmkeys)<-c("Dim","Total.Dim","%Dim")	
DimensionWord<-matmkeys[order(matmkeys[,1],decreasing = T),]	

 # Graph of metakeys and Keydocs on the "axes" of the CA
 axes<-c(axx,axy)
 # Identify metakeys of the dimensions "axes"
 ClusPal1<-which(Ccontr[,axes[1]]%in%Metakeys[[axes[1]]][[1]])
 ClusPal2<-which(Ccontr[,axes[1]]%in%Metakeys[[axes[1]]][[2]])
 ClusPal3<-which(Ccontr[,axes[2]]%in%Metakeys[[axes[2]]][[1]])
 ClusPal4<-which(Ccontr[,axes[2]]%in%Metakeys[[axes[2]]][[2]])
 ClusPal<-c(ClusPal1,ClusPal2,ClusPal3,ClusPal4)
 ClusPal<-ClusPal[!duplicated(ClusPal)]

 
 
 # Identify Keydocs of the dimensions "axes"
 ClusDoc1<-which(Rcontr[,axes[1]]%in%Keydocs[[axes[1]]][[1]])
 ClusDoc2<-which(Rcontr[,axes[1]]%in%Keydocs[[axes[1]]][[2]])
 ClusDoc3<-which(Rcontr[,axes[2]]%in%Keydocs[[axes[2]]][[1]])
 ClusDoc4<-which(Rcontr[,axes[2]]%in%Keydocs[[axes[2]]][[2]])
 ClusDoc<-c(ClusDoc1,ClusDoc2,ClusDoc3,ClusDoc4)
 ClusDoc<-ClusDoc[!duplicated(ClusDoc)]

 if(length(ClusPal)>0&length(ClusDoc)>0){
 xl=c(min(Ccoor[ClusPal,axes[1]],Rcoor[ClusDoc,axes[1]]),max(Ccoor[ClusPal,axes[1]],Rcoor[ClusDoc,axes[1]]))
 yl=c(min(Ccoor[ClusPal,axes[2]],Rcoor[ClusDoc,axes[2]]),max(Ccoor[ClusPal,axes[2]],Rcoor[ClusDoc,axes[2]]))}

 if(length(ClusPal)>0&length(ClusDoc)==0){
 xl=c(min(Ccoor[ClusPal,axes[1]]),max(Ccoor[ClusPal,axes[1]]))
 yl=c(min(Ccoor[ClusPal,axes[2]]),max(Ccoor[ClusPal,axes[2]]))}

if(length(ClusPal)==0&length(ClusDoc)>0){
 xl=c(min(Rcoor[ClusDoc,axes[1]]),max(Rcoor[ClusDoc,axes[1]]))
 yl=c(min(Rcoor[ClusDoc,axes[2]]),max(Rcoor[ClusDoc,axes[2]]))}

 


if(length(ClusPal)>0|length(ClusDoc)>0){
dfW<-dfD<-NULL
if(length(ClusPal)>0) {
 dfW <- data.frame(Words=character(),Dim=integer(), Coordinates=integer(), stringsAsFactors = FALSE)
for(i in 1:naxes) {
for(k in 1:2){
Q <- ifelse(k==1,1,-1)
if(length(Metakeys[[i]][[k]])>0)  {
   for(j in 1: length(Metakeys[[i]][[k]]))  {
      dfW[nrow(dfW)+1,] <- c(names(Metakeys[[i]][[k]])[j], i, Q*Metakeys[[i]][[k]][j])}
}}}


if(length(ClusDoc)>0) {
 dfD <- data.frame(Docs=character(),Dim=integer(), Coordinates=integer(), stringsAsFactors = FALSE)
for(i in 1:naxes) {
 for(k in 1:2){
Q <- ifelse(k==1,1,-1)
if(length(Keydocs[[i]][[k]])>0) {
  for(j in 1: length(Keydocs[[i]][[k]]))  {
 dfD[nrow(dfD)+1,] <- c(names(Keydocs[[i]][[k]])[j], i, Q* Keydocs[[i]][[k]][j])}
}}}
}
keysR <- list(Word=dfW, Doc = dfD, lmd=lmd, lmw=lmw)
 return(keysR)
}

}else{
 cat("\nThere are no elements to plot, lower lmd and/or lmw values\n")
}
}

 Xtab.disjonctif <- function(tab)
 {
   tab <- as.data.frame(tab,  stringsAsFactors = TRUE)
   modalite.disjonctif <- function(i) {
     moda <- as.factor(tab[, i])
     n <- length(moda)
     x <- matrix(0L, n, nlevels(moda))
     x[(1:n) + n * (unclass(moda) - 1L)] <- 1L
     return(x)
   }
   if (ncol(tab) == 1) {
     res <- modalite.disjonctif(1)
     dimnames(res) <- list(attributes(tab)$row.names, levels(tab[,1]))
   }
   else {
     variable <- rep(attributes(tab)$names, lapply(tab, nlevels))
     listModa <- unlist(lapply(tab, levels))
     wlistModa <- which((listModa) %in% c("y", "n", 
                                          "Y", "N"))
     if (!is.null(wlistModa)) 
       listModa[wlistModa] <- paste(variable[wlistModa], 
                                    listModa[wlistModa], sep = ".")
     numlistModa <- which(unlist(lapply(listModa, is.numeric)))
     if (!is.null(numlistModa)) 
       listModa[numlistModa] <- paste(variable[numlistModa], 
                                      listModa[numlistModa], sep = ".")
     res <- lapply(1:ncol(tab), modalite.disjonctif)
     res <- as.matrix(data.frame(res, check.names = FALSE))
     dimnames(res) <- list(attributes(tab)$row.names, listModa)
   }
   return(res)
 }
 
 
 # End Functions -------------------------

 var.text <- object$info$var.text[[1]]
 str.base <- object$info$base[[1]]
 str.envir <- object$info$menvir[[1]]
 
 

info <- list(var.text= object$info$var.text, base=object$info$base,
  menvir= object$info$menvir, catquali = object$context$quali)

naxes <- ncp
context <- context.sup
var.agg <- object$info$name.var.agg[[1]]
quali.sup <- NULL; quanti.sup<-NULL	
ndocsup <-0; nsegm <-0	
nwordsup <- 0; nquali<-0; nquanti<-0	



if(segment==TRUE)	
 if(is.null(object$DocSeg))stop("No segment selected in TextData object")


if(var.agg=="") {
context.quanti<-context.quali<-NULL

if(!is.null(context)) {
 if(length(context)==1)
  if(context=="ALL") {
  context<- NULL
 if(!is.null(object$context$quali)) context <- as.character(colnames(object$context$quali))
 if(!is.null(object$context$quanti)) context <- c(context,colnames(object$context$quanti))
   }
  if(!is.null(object$context$quanti)){
   context.quanti <- checkcontext(context,object$context$quanti)
   context.quanti <- colnames(object$context$quanti)[context.quanti]

  fval <- apply(object$context$quanti[context.quanti], 2, function(x) max(x)-min(x))
  fval <- context.quanti[which(fval==0)]
  if(length(fval)>0){
   warning("Quantitative variables with the same value are removed: ",fval)
   context.quanti <- context.quanti[- which(fval==0)]
  }
  if(length(context.quanti)==0) context.quanti<-NULL
    } else {context.quanti<-NULL}

  if(!is.null(as.character(colnames(object$context$quali)))){
   context.quali <- checkcontext(context,object$context$quali)
   context.quali <- colnames(object$context$quali)[context.quali]

   if(length(context.quali)==0) context.quali<-NULL

} else { context.quali<- NULL}
} # Final if(!is.null(context))



# ========== Check quali variables $context$quali no aggregate
if(!is.null(context.quali))	
object$context$quali <- data.frame(object$context$quali[,context.quali,drop=FALSE] )

# =======  Check quanti variables $context$quanti no aggregate
if(!is.null(context.quanti)) 	
object$context$quanti <- data.frame(object$context$quanti[,context.quanti,drop=FALSE])	

LT <- as.matrix(object$DocTerm)	
ndoc <- nrow(LT)	
if(segment==TRUE){  LTS <- as.matrix(object$DocSeg)
 rownames(LTS) <- rownames(LT)}	

# =================== Selecting supplementary documents docs by rownumber or rowname	
if(!is.null(doc.sup)) {	
      if (is.character(doc.sup)) 	
         doc.sup <- which(rownames(LT) %in% doc.sup)	
         doc.sup <- rownames(LT)[doc.sup]	 
         doc.sup <- doc.sup[!is.na(doc.sup)]	
         doc.sup <- which(rownames(LT)%in% doc.sup)	
         ndocsup <- length(doc.sup)	
 if(ndocsup==0) doc.sup <- NULL	
}	
ndocact <- ndoc - ndocsup


# Put supplementary docs at the end	
if(ndocsup>0) {	
  tmpdf <- LT[doc.sup,,drop=FALSE]
  LT <- LT[-doc.sup,,drop=FALSE] 
  LT <- rbind(LT,tmpdf)
 if(!is.null(context.quali)) if(length(context.quali)>0)  	
     object$context$quali<- remdocs(object$context$quali,doc.sup)	
 if(!is.null(context.quanti)) if(length(context.quanti)>0)	
  object$context$quanti<- remdocs(object$context$quanti,doc.sup)

  #_____________________________________________
if(segment==TRUE){	
 tmpdf <- LTS[doc.sup,,drop=FALSE]
  LTS <- LTS[-doc.sup,,drop=FALSE] 
  LTS <- rbind(LTS,tmpdf)
}} # final de if(!is.null(doc.sup))	


#_____________________________________________	
nword <- ncol(LT)
nwordact <- nword
# Put supplementary words at the right of LT	
if(!is.null(word.sup)) {	
      if (!is.character(word.sup)) 					
      word.sup <- colnames(LT)[word.sup]
      word.sup <- which(colnames(LT) %in% word.sup)	
      nwordsup <- length(word.sup)

   if(nwordsup==0){ word.sup <- NULL} else {
    LT <- cbind(LT,LT[,word.sup])	
    wcolnames <- colnames(LT)[word.sup]
     LT <- LT[,-word.sup]	
   colnames(LT)[(ncol(LT)-nwordsup+1):ncol(LT)] <- wcolnames
   nwordact <- ncol(LT) - nwordsup}	}


#_____________________________________________	
# Remove supplementary words with zero sum	
if(!is.null(word.sup)) if(nwordsup>0) {
  pos.elim<-which(colSums(LT[c(1:ndocact),c((nwordact+1):nword),drop=FALSE])==0)	
if(length(pos.elim)>0)	
    LT <- LT[,-pos.elim] 
  nwordsup <- nwordsup - length(pos.elim)
   nword <- ncol(LT)}

#_____________________________________________	
# Remove words with zero sum	
  pos.elim<-which(colSums(LT[c(1:ndocact),,drop=FALSE])==0)	
  if(length(pos.elim)>0)	
    LT <- LT[,-pos.elim] 	
  nword <- ncol(LT)	
  nwordact <- nword - nwordsup

#_____________________________________________	
if(segment==TRUE){	
# Remove segments with zero sum	
pos.elim<-which(colSums(LTS[c(1:ndocact),,drop=FALSE])==0)	
  if(length(pos.elim)>0)	
    LTS <- LTS[,-pos.elim] 	
  nsegm <- ncol(LTS)	
   LT <- cbind(LT,LTS)	}	

  
# ==== Yuxtapostition of qualitative and quantitative variables 	
 if(length(context.quali)>0) {
      yQL <- data.frame(object$context$quali[rownames(LT),context.quali])
      colnames(yQL) <- colnames(object$context$quali)
	LT <- cbind(LT,yQL)
	nquali <- ncol(yQL)}
 if(length(context.quanti)>0) {	
      yQ <- data.frame(object$context$quanti[rownames(LT),context.quanti])
      colnames(yQ) <- colnames(object$context$quanti)
	LT <- cbind(LT,yQ)
      nquanti <- ncol(yQ)}	

  
#_____________________________________________
# ==== Remove supplementary empty docs 
if(ndocsup>0) {
pos.elim<-which(rowSums(LT[c((ndocact+1):ndoc),c(1:nwordact),drop=FALSE])==0)
  if(length(pos.elim)>0){
  LT <- LT[-pos.elim,] 
  ndocsup <- ndocsup - length(pos.elim)
}}


#_____________________________________________
# ==== Remove empty docs 
  pos.elim<-which(rowSums(LT[c(1:nrow(LT)),c(1:nwordact),drop=FALSE])==0)
if(length(pos.elim)>0)
  LT <- LT[-pos.elim,] 
 ndoc <- nrow(LT)
 ndocact <- ndoc-ndocsup

 

 
 
#_____________________________________________
# Fill NA quantitative values with the average
if(nquanti >0){
posquanti <- ncol(LT)-nquanti+1
for(i in (posquanti:ncol(LT))){
if(any(is.na(LT[,i]))) warning("\n", names(LT)[i], " variable has missing values. They will be replaced by the mean\n") 
  LT[is.na(LT[,i]), i] <- mean(LT[,i], na.rm = TRUE)
}}


if(ndocact <3)
    stop("At least three active documents have to be selected")
if(nwordact <3)
    stop("At least three active words have to be selected")
minrowcol <- min((ndocact-1), (nwordact-1))
ncp <- min(ncp, minrowcol)


if(ndocsup >0){ row.sup <- c(c(ndocact+1):ndoc)} else {row.sup <-NULL}

pos.word.sup <- pos.seg.sup <- NULL
ncoltemp <- nwordact
if(nwordsup >0){ 
  pos.word.sup <- c((ncoltemp +1):(ncoltemp+nwordsup))   
  ncoltemp <- ncoltemp + nwordsup
}

if(nsegm >0){
  pos.seg.sup <- c(c(ncoltemp +1):c(ncoltemp+nsegm))       # 1081 al 2713  ; Hay 1487 segmentos
  ncoltemp <- ncoltemp + nsegm
}
pos.col.sup <- c(pos.word.sup,pos.seg.sup)


if(nquali >0){ quali.sup <- c((ncoltemp+1):(ncoltemp+nquali)); ncoltemp <- ncoltemp+nquali}
if(nquanti >0){ quanti.sup <- c((ncoltemp+1):(ncoltemp+nquanti)); ncoltemp <- ncoltemp+nquanti}

# res <- FactoMineR::CA(LT, ncp, row.sup=row.sup, col.sup=word.seg.sup, quali.sup=quali.sup, quanti.sup=quanti.sup, graph=FALSE)
# res <- FactoMineR::CA(LT, ncp, row.sup=row.sup, col.sup=pos.col.sup, quali.sup=quali.sup, quanti.sup=quanti.sup, graph=FALSE)
res <- CA_New(LT, ncp, row.sup=row.sup, col.sup=pos.col.sup, quali.sup=quali.sup, quanti.sup=quanti.sup, graph=FALSE)
res$meta <- Keys(res,lmd,lmw,naxes, axes)



#################################
############ Modificada la siguiente:
##### res$VCr <- round(sqrt(sum(res$eig[, 1])/ (minrowcol-1) ), 4)
 res$VCr <- round(sqrt(sum(res$eig[, 1])/ minrowcol ), 4)
 res$Inertia <- round(sum(res$eig[, 1]), 4)
 res$info <- info

 
 if(segment==TRUE) if(ncol(LTS)==0) segment<-NULL
 
 if(segment==TRUE) {
  nsegm <- ncol(LTS)	

 if(nwordsup>0){
 res$segment <- res$col.sup
 res$segment$coord <- res$segment$coord[-c(1:nwordsup),,drop=FALSE]
 res$segment$cos2 <- res$segment$cos2[-c(1:nwordsup),,drop=FALSE]
 res$col.sup$coord <- res$col.sup$coord[c(1:nwordsup),,drop=FALSE]
 res$col.sup$cos2 <- res$col.sup$cos2[c(1:nwordsup),,drop=FALSE] 
} else {
 res$segment <- ""
 suppressWarnings(res$segment$coord <-res$col.sup$coord[(1:nsegm),,drop=FALSE]) 
 suppressWarnings(res$segment$cos2 <- res$col.sup$cos2[(1:nsegm),,drop=FALSE])
}
  }


if(!is.null(res$quali.sup$coord)) {
  res$quali.sup$cos2 <- res$quali.sup$coord[,,drop=FALSE]
  if(is.null(row.sup))
    LTnosup <- LT
  else LTnosup <- LT[-row.sup,, drop=FALSE]

  
   Xact <- as.matrix(LTnosup[,1:nwordact, drop=FALSE]) 
   PJ <- colSums(Xact)/sum(Xact)
   tt <-as.matrix( t(Xtab.disjonctif(as.matrix(LTnosup[,quali.sup,drop=FALSE]))))
   PIJ <- tt %*% (Xact/sum(Xact))
  disto2 <- ((t(PIJ/rowSums(PIJ)) -PJ)^2)/PJ
  disto2 <- colSums(disto2)
  cos2 <- res$quali.sup$coord^2 
  res$quali.sup$cos2 <- cos2 / disto2
}

res$rowINIT <- object$rowINIT 
 class(res) <- c("LexCA","CA", "list")
 
if(graph==TRUE) plotLexCA()
return(res)

 
 
 
 
 } else {


## Aggregate analysis
context.quanti<-context.quali<-NULL
if(!is.null(context)) {
ctmp <- NULL
if(length(context)==1)
if(context=="ALL") {
 context<- NULL
 if(!is.null(object$context$quali$qualivar[,1])) context <- as.character(object$context$quali$qualivar[,1])
 if(!is.null(object$context$quanti)) context <- c(context,colnames(object$context$quanti))
}

  if(!is.null(object$context$quanti)){
   context.quanti <- checkcontext(context,object$context$quanti)
   context.quanti <- colnames(object$context$quanti)[context.quanti]
  if(length(context.quanti)==0) context.quanti<-NULL
    } else {context.quanti<-NULL}

  if(!is.null(as.character(object$context$quali$qualivar[,1]))){
   context.quali <- as.character(object$context$quali$qualivar[,1])


if(is.null(context)) context<-""
 if(is.character(context)) 
        context_ <- which(context %in% context.quali)
        context_ <- as.character(object$context$quali$qualivar[context_,1])
        context_ <-  context_[!is.na(context_)]
        context.quali <- which(context_ %in% context.quali)
 if(length(context.quali )==0) context.quali <- NULL
} else { context.quali<- NULL}


} # Final if(!is.null(context))


# =======  Check quanti variables $context$quanti aggregate	
if(!is.null(context.quanti))  
        data.frame(object$context$quanti[,context.quanti,drop=FALSE])

# Aggregate lexical table	
LT <- as.matrix(object$DocTerm)	
ndoc <- nrow(LT)	
if(segment==TRUE) {
 LTS <- data.frame(as.matrix(object$DocSeg)[,,drop=FALSE])
 rownames(LTS) <- rownames(LT)}

# =================== Selecting supplementary documents docs by rownumber or rowname					
if(!is.null(doc.sup)) {					
      if (is.character(doc.sup)) 					
         doc.sup <- which(rownames(LT) %in% doc.sup)					
         doc.sup <- rownames(LT)[doc.sup]		
         doc.sup <- doc.sup[!is.na(doc.sup)]					
         doc.sup <- which(rownames(LT)%in% doc.sup)					
 ndocsup <- length(doc.sup)					
 if(ndocsup==0) doc.sup <- NULL					
}					
ndocact <- ndoc - ndocsup

LTQ <- NULL; LTQL<- NULL
if(!is.null(context.quanti)){
 LTQ <- data.frame(object$context$quanti[,context.quanti,drop=FALSE])
 nquanti <- length(LTQ)}

# Put supplementary docs at the end
if(ndocsup>0) {
  tmpdf <- LT[doc.sup,,drop=FALSE]
  LT <- LT[-doc.sup,,drop=FALSE] 
  LT <- rbind(LT,tmpdf)		

  if(segment==TRUE){	
  tmpdf <- LTS[doc.sup,,drop=FALSE]
  LTS <- LTS[-doc.sup,,drop=FALSE] 
  LTS <- rbind(LTS,tmpdf)	
}

 if(!is.null(context.quanti)){ 
  tmpdf <- LTQ[doc.sup,,drop=FALSE]
  LTQ <- LTQ[-doc.sup,,drop=FALSE] 
  LTQ <- rbind(LTQ,tmpdf)
}

ndocact <- ndoc-ndocsup
if(ndocact <3)
    stop("At least three active documents have to be selected")
} # Final if(ndocsup>0) 

 
if(!is.null(context.quali))
   LTQL <- data.frame(object$context$quali$qualitable[,,drop=FALSE])

nwordact <- ncol(LT); nword<- nwordact	



#==================================================================
# Put supplementary words at the right word.sup	
if(!is.null(word.sup)) {	
   if (!is.character(word.sup)) 					
    word.sup <- colnames(LT)[word.sup]
    word.sup <- which(colnames(LT) %in% word.sup)
    nwordsup <- length(word.sup)	

 if(nwordsup==0){ word.sup <- NULL} else {
    LTSUP <- data.frame(LT[,word.sup,drop=FALSE])	
    LT <- cbind(LT,LTSUP)	
    LT <- data.frame(LT[,-word.sup,drop=FALSE])	
 if(!is.null(context.quali)){ 	
    LTQLSUP <- data.frame(LTQL[,word.sup,drop=FALSE])	
    LTQL <- cbind(LTQL,LTQLSUP)	
    LTQL <- data.frame(LTQL[,-word.sup,drop=FALSE])  	
}}	
nwordact <- nword - nwordsup	
}	


########### Remove supplementary words with zero sum			
if(nwordsup>0) {			
  pos.elim<-which(colSums(LT[c(1:ndocact),c((nwordact+1):nword),drop=FALSE])==0)			
if(length(pos.elim)>0) {			
 if(!is.null(LTQL))  LTQL <- LTQL[,-pos.elim] 			
  LT <- LT[,-pos.elim] 			
nword <- ncol(LT)			
nwordsup <- nword-nwordact			
}}			

if(nwordact <3)
    stop("At least three active words have to be selected")
		
# Removing active words with zero margins			
if(ndocact>0) {			
  pos.elim<-which(colSums(LT[c(1:ndocact),c(1:nwordact),drop=FALSE])==0)			
  if(length(pos.elim)>0) {			
  LT <- LT[,-pos.elim, drop=FALSE] 			
 if(!is.null(LTQL))			
  LTQL <- LTQL[,-pos.elim, drop=FALSE] 			
nword <- ncol(LT)			
nwordact <- nword-nwordsup 			
}}			


# Fill NA quantitative values with the average	
if(!is.null(LTQ)){	
for(i in (1:ncol(LTQ))){
if(any(is.na(LTQ[,i]))) warning("\n", names(LTQ)[i], " variable has missing values. They will be replaced by the mean\n") 
  LTQ[is.na(LTQ[,i]), i] <- mean(LTQ[,i], na.rm = TRUE)	
}}	

if(ndocact <3)		
    stop("At least three active documents have to be selected")	
if(nwordact <3)		
    stop("At least three active words have to be selected")
	
minrowcol <- min((ndocact-1), (nwordact-1))		
ncp <- min(ncp, minrowcol)		

if(segment==FALSE) nsegm<-0 else nsegm <- ncol(LTS)
if(ndocsup >0){ row.sup <- c((ndocact+1):(ndocact+ndocsup))} else {row.sup <-NULL}
word.sup <-NULL
ncoltemp <- nwordact
if(nwordsup >0){ word.sup <- c((ncoltemp+1):(nwordact+nwordsup)); ncoltemp <- ncoltemp+nwordsup}
if(nsegm >0){ word.sup <- c(word.sup,c((ncoltemp+1):(ncoltemp+nsegm))); ncoltemp <- ncoltemp+nsegm
LT <- cbind(LT,LTS)}


if(!is.null(LTQL)) nquali <- nrow(LTQL) else nquali<-0

if(nquali >0){
  doc.sup.q <- c((ndocact+ndocsup+1):(ndocact+ndocsup+nquali))
  ndifcol <- ncol(LT)-ncol(LTQL) 
if(ndifcol>0) {
  tempcol <- matrix(c(rep.int(NA,(nquali*ndifcol))),nrow=nquali,ncol=ndifcol)
  newcol <- data.frame(tempcol)
  colnames(newcol) <- colnames(LT)[(ncol(LTQL)+1):ncol(LT)]
  LTQL <- cbind(LTQL,newcol)
 }
  colnames(LTQL) <- colnames(LT)
  LT <- rbind(LT,LTQL)
}

quanti.sup.q <- NULL
if(nquanti >0){ 
 quanti.sup.q <- c((ncol(LT)+1):(ncoltemp+nquanti))
 ndif <- nrow(LT)-nrow(LTQ)
 if(ndif >0) {
 temprow <- matrix(c(rep.int(NA,(ndif*ncol(LTQ)))),nrow=ndif,ncol=ncol(LTQ))
 newrow <- data.frame(temprow)[,,drop=FALSE]
 rownames(newrow) <- rownames(LT)[(nrow(LT)-ndif+1):nrow(LT)]
 colnames(newrow) <- colnames(LTQ)
 LTQ <- rbind(LTQ,newrow)
 } 
 LT <- cbind(LT,LTQ)
}



if((ndocsup+nquali)>0) doc.sup <- c((ndocact+1):(ndocact+ndocsup+nquali))
res <- CA_New(LT, ncp, row.sup=doc.sup, col.sup=word.sup, quanti.sup=quanti.sup.q, graph=FALSE)
res$var.agg <- var.agg
res$meta <- Keys(res,lmd,lmw,naxes, axes)
####### Modificar lo siguiente
############################
######res$VCr <- round(sqrt(sum(res$eig[, 1])/ (minrowcol-1) ), 4)
res$VCr <- round(sqrt(sum(res$eig[, 1])/ minrowcol ), 4)
res$Inertia <- round(sum(res$eig[, 1]), 4)
res$info <- info



if(segment==TRUE){
if(nwordsup>0){
res$segment <- res$col.sup
res$segment$coord <- res$segment$coord[-(1:nwordsup),]
res$segment$cos2 <- res$segment$cos2[-(1:nwordsup),]
} else {
res$segment <- ""
 suppressWarnings(res$segment$coord <-res$col.sup$coord[(1:nsegm),]) 
suppressWarnings(res$segment$cos2 <- res$col.sup$cos2[(1:nsegm),])
}
 if(!is.null(var.agg))
{
   rownames(res$segment$coord)<- substring(rownames(res$segment$coord),2)
   aa <- sub('\\.', ':', rownames(res$segment$coord)) 
   rownames(res$segment$coord) <- gsub("\\."," ",aa)
   rownames(res$segment$cos2) <-  rownames(res$segment$coord)
   rownames(res$col.sup$coord)[(nwordsup+1):nrow(res$col.sup$coord)] <- rownames(res$segment$coord)
   rownames(res$col.sup$cos2)<-  rownames(res$col.sup$coord)
}
}


if(nquali >0) {
 res$quali.sup <- res$row.sup
 if(ndocsup>0){
 res$quali.sup$coord <- res$quali.sup$coord[-(1:ndocsup),,drop=FALSE]
 res$quali.sup$cos2 <- res$quali.sup$cos2[-(1:ndocsup),,drop=FALSE]
 res$row.sup$coord  <- res$row.sup$coord[(1:ndocsup),,drop=FALSE]
 res$row.sup$cos2  <- res$row.sup$cos2[(1:ndocsup),,drop=FALSE]}
}
res$rowINIT <- object$rowINIT 
 class(res) <- c("LexCA","CA", "list")

 if(graph==TRUE) plotLexCA()
 return(res)
}
}