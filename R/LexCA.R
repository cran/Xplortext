#' @import FactoMineR
#' @export
LexCA <- function(object,ncp=5, context.sup="ALL", doc.sup=NULL, word.sup=NULL, 
 segment=FALSE, graph=TRUE, axes=c(1, 2), lmd=3, lmw=3)
{
if(is.null(object)) stop("Missing argument for object")
if (!inherits(object,"TextData")) stop("Object should be TextData class")

# Functions -------------------
plotLexCA <- function()
{
# if(dev.interactive()) dev.new()
plot.LexCA(res, selDoc=NULL, selWord=NULL, eigen=TRUE, axes=axes)

if(!is.null(res$quanti.sup$coord)) {
#  if(dev.interactive()) dev.new()
 plot.LexCA(res, selDoc=NULL, selWord=NULL, quanti.sup=rownames(res$quanti.sup$coord), eigen=FALSE, axes=axes)
}

#if(dev.interactive()) dev.new()
  plot.LexCA(res, selDoc= paste("meta ", lmd), selWord= paste("meta ", lmw), eigen=FALSE, axes=axes)

if(!is.null(rownames(res$segment$coord))& segment==TRUE){
#  if(dev.interactive()) dev.new()
  plot.LexCA(res, selDoc= NULL, selWord= NULL, selSeg="ALL", axes=axes)}

if(!is.null(res$quali.sup)) {
# if(dev.interactive()) dev.new()
 plot.LexCA(res, selDoc= NULL, selWord= NULL, quali.sup="ALL", axes=axes)
 }

if(!is.null(res$row.sup$coord)){
   ident <- identical(dimnames(res$quali.sup$coord)[1], dimnames(res$row.sup$coord)[1])
   if(ident==FALSE) {
#  if(dev.interactive()) dev.new()
   plot.LexCA(res, selDoc= NULL, selWord= NULL, selDocSup="ALL", axes=axes)
 }
}

if(!is.null(res$col.sup$coord)){
   ident <- identical(dimnames(res$col.sup$coord)[1], dimnames(res$segment$coord)[1])
   if(ident==FALSE) {
#  if(dev.interactive()) dev.new()
   plot.LexCA(res, selDoc= NULL, selWord= NULL, selWordSup="ALL", selSeg=NULL, axes=axes)
 }
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
object$context$quali <- data.frame(object$context$quali[,context.quali,drop=FALSE])	

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








#_____________________________________________	
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
word.sup <-NULL
ncoltemp <- nwordact
if(nwordsup >0){ word.sup <- c((nwordact+1):(nwordact+nwordsup)); ncoltemp <- ncoltemp+nwordsup}
if(nsegm >0){ word.supseg <- c(c(c(ncoltemp+1):c(ncoltemp+nsegm))); ncoltemp <- ncoltemp+nsegm

word.sup <- c(word.sup,word.supseg) }
if(nquali >0){ quali.sup <- c((ncoltemp+1):(ncoltemp+nquali)); ncoltemp <- ncoltemp+nquali}
if(nquanti >0){ quanti.sup <- c((ncoltemp+1):(ncoltemp+nquanti)); ncoltemp <- ncoltemp+nquanti}


res <- CA(LT, ncp, row.sup=row.sup, col.sup=word.sup, quali.sup=quali.sup, quanti.sup=quanti.sup, graph=FALSE)
 res$meta <- Keys(res,lmd,lmw,naxes, axes)
res$VCr <- round(sqrt(sum(res$eig[, 1])/ (minrowcol-1) ), 4)
res$Inertia <- round(sum(res$eig[, 1]), 4)
res$info <- info

if(segment==TRUE) if(ncol(LTS)==0) segment<-NULL
if(segment==TRUE) {
 nsegm <- ncol(LTS)	
 if(nwordsup>0){
 res$segment <- res$col.sup
 res$segment$coord <- res$segment$coord[-c(1:nwordsup),,drop=FALSE]
 res$segment$cos2 <- res$segment$cos2[-c(1:nwordsup),,drop=FALSE]
} else {
 res$segment <- ""
 suppressWarnings(res$segment$coord <-res$col.sup$coord[(1:nsegm),,drop=FALSE]) 
 suppressWarnings(res$segment$cos2 <- res$col.sup$cos2[(1:nsegm),,drop=FALSE])
}}

if(!is.null(res$quali.sup$coord)) {
  res$quali.sup$cos2 <- res$quali.sup$coord[,,drop=FALSE]
  Xact <- as.matrix( LT[,1:nwordact, drop=FALSE])
  PJ <- colSums(Xact)/sum(Xact)
  tt <-as.matrix( t(tab.disjonctif(as.matrix(LT[,quali.sup,drop=FALSE]))))
  PIJ <- tt %*% (Xact/sum(Xact))
  disto2 <- ((t(PIJ/rowSums(PIJ)) -PJ)^2)/PJ
  disto2 <- colSums(disto2)
  cos2 <- res$quali.sup$coord^2 
  res$quali.sup$cos2 <- cos2 / disto2
}

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
res <- CA(LT, ncp, row.sup=doc.sup, col.sup=word.sup, quanti.sup=quanti.sup.q, graph=FALSE)
res$var.agg <- var.agg
res$meta <- Keys(res,lmd,lmw,naxes, axes)
res$VCr <- round(sqrt(sum(res$eig[, 1])/ (minrowcol-1) ), 4)
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

 class(res) <- c("LexCA","CA", "list")

 if(graph==TRUE) plotLexCA()
 return(res)
}
}