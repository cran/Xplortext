#' @importFrom methods hasArg
#' @importFrom graphics barplot
#' @export
plot.LexCA <- function(x, selDoc="ALL", selWord="ALL", selSeg=NULL,
  selDocSup=NULL, selWordSup=NULL, quanti.sup=NULL, quali.sup=NULL, maxDocs=20,  
  eigen=FALSE, title=NULL, axes=c(1,2), col.doc="blue", col.word="red",
  col.doc.sup="darkblue", col.word.sup="darkred", col.quanti.sup = "blue",
  col.quali.sup="darkgreen", col.seg="cyan4",col="grey", cex=1, 
  xlim=NULL, ylim=NULL, shadowtext=FALSE, habillage="none", unselect=1,
  label="ALL", autoLab=c("auto", "yes", "no"), new.plot=TRUE, 
  graph.type = c("classic", "ggplot"),...) 
{
  
 if (!inherits(x, "LexCA"))  stop("x object should be LexCA class")
  options(stringsAsFactors = FALSE)
  if(is.null(label)) label <- "none"
  if(label=="ALL") label<-"all"
  if(length(graph.type)>1) graph.type <- "classic"

if(eigen==TRUE) {selDoc<-selWord<-NULL}
if(!is.null(quali.sup)) {
 if(!is.null(x$var.agg)) { 
 if(length(quali.sup)==1) if(quali.sup=="ALL") 
 quali.sup <- as.character(x$info$catquali$qualivar[,1])
 df1 <- data.frame(x$info$catquali$qualivar)
 df1$inic <- cumsum(df1$qualincat) 
 df1$inic <- df1$inic - df1$qualincat+1 
 posic <- which(df1[,1] %in% quali.sup)
 df2 <- df1[posic,]  
 if(nrow(df2)>0)
  for(i in 1:nrow(df2)) {
  quali.sup <- c(quali.sup,dimnames(x$info$catquali$qualitable)[[1]][df2$inic[i]: (df2[i,"inic"]+df2[i,"qualincat"]-1)])
  if(length(quali.sup)==0) quali.sup <- NULL}
} else {
# Not aggregate
 if(length(quali.sup)==1) if(quali.sup=="ALL") 
  quali.sup <- as.character(rownames(x$quali.sup$coord))
  namevar <- colnames(x$info$catquali)
  posic <- which(namevar %in% quali.sup)
  if(length(posic)>0) {  # There are variables
    temp1 <- sapply(x$info$catquali, levels)[posic]
    temp2 <- unlist(temp1, use.names=FALSE)
    temp3 <- rep(names(temp1), lengths(temp1))
    temp3 <- paste0(temp3,".",temp2)
    quali.sup <- c(quali.sup, unlist(temp1),unlist(temp3))
    posic <- which(rownames(x$quali.sup$coord) %in% quali.sup)
    quali.sup <- rownames(x$quali.sup$coord)[posic]}
} # Final not aggregate
}

selection <- function(sel1, xobj, bType, axx, axy)
{
if(is.null(sel1)) return(sel1)
 xx <- ""
  if(length(sel1)==1){

     if(sel1=="ALL") sel1 <- c(1:dim(xobj$coord)[1])
  }
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
   if(xx=="cos2") {
 # Selection by cos2
 sel1 <- as.numeric(substr(sel1, 5, nc))
 dft <- data.frame(xobj$cos2[,c(axx,axy),drop=FALSE])
 fval <- apply(dft, 1, function(z) sum(z))
 posic <- which(fval>= sel1)
 sel1 <- rownames(xobj$cos2)[posic]
   }
   if(xx=="contrib") {
 # Selection by contrib
 if(bType=="Seg") stop("Segments can not be selected by contribution")
 if(bType=="Quali") stop("Contextual categorical variables can not be selected by contribution")
 if(bType=="Quanti") stop("Contextual quantitative variables can not be selected by contribution")
 if(bType=="Dsup") stop("Supplementary documents can not be selected by contribution")
 if(bType=="Wsup") stop("Supplementary words can not be selected by contribution")
 sel1 <- as.numeric(substr(sel1, 9, nc))
 dft <- data.frame(xobj$contrib[,c(axx,axy),drop=FALSE])
 fval <- apply(dft, 1, function(z) max(abs(z)))
 posic <- which(fval>= sel1)
 sel1 <- rownames(xobj$contrib)[posic]
  }
 if(xx=="meta") {
 if(bType=="Seg") stop("Segments can not be selected by meta")
 if(bType=="Quali") stop("Contextual categorical variables can not be selected by meta")
 if(bType=="Quanti") stop("Contextual quantitative variables can not be selected by meta")
 if(bType=="Dsup") stop("Supplementary documents can not be selected by meta")
 if(bType=="Wsup") stop("Supplementary words can not be selected by meta")
 sel1 <- as.numeric(substr(sel1, 5, nc))
 sMeta <-  rownames(xobj$coord)[which(xobj$contrib[,axx]>mean(xobj$contrib[,axx])*sel1)]
 sMeta <-  c(sMeta,rownames(xobj$coord)[which(xobj$contrib[,axy]>mean(xobj$contrib[,axy])*sel1)])
 sel1 <- unique(sMeta)
}
 if(xx=="char") {
 if(bType!="Word") stop("char only can be used to select active words (selWord)")
  proba <- as.numeric(substr(sel1, 5, nc))
 if(proba<0|proba>1) stop("proba should be between 0 and 1")
 if(is.null(selDoc)) stop("selDoc can not be NULL when words are selected by char")
 ntdoc <- nrow(x$call$X) 
 if(ntdoc>maxDocs) stop("Computing the characteristic words when more than ", maxDocs, " documents is not allowed")
 resCharWord <- descfreq(x$call$X,proba=proba)
 sW <- NULL
 for (idoc in 1:ntdoc)
    {
 if(rownames(x$row$coord)[idoc] %in% selDoc){
     df <- data.frame(resCharWord[idoc])
     df <-rownames(df[df[6]>0,])
     sW <- c(sW,df)}
     }
 sel1 <- unique(sW)
 }
 } else {
 if(is.character(sel1)) sel1 <- which(rownames(xobj$coord) %in% sel1)
 sel1 <- rownames(xobj$coord)[sel1]
 sel1  <- sel1[!is.na(sel1)]
}
  return(sel1)
}
# Final functions




if(eigen) {

 if(is.null(title)) titleE <- "Eigenvalues" else titleE <- title
  args <- list(...)
  exist <- "names.arg" %in% names(args)

 if(hasArg(names.arg)) {
     valarg <- args$names.arg
     if(length(valarg)==1) {
      if(valarg=="") names.arg <- rep(valarg,nrow(x$eig))
      else names.arg <- paste0(valarg,1:nrow(x$eig))
       }
          else names.arg <- args$names.arg
 } else  names.arg <- paste("dim",1:nrow(x$eig))

  

 if(new.plot & graph.type!="interact") dev.new()
 if(graph.type == "classic") { 
   barplot(x$eig[, 1], main = titleE, col=col, names.arg = names.arg) 
      }  else  {   # No es classic
      x<- names.arg ; y <- x$eig[,1]
      df <- data.frame(x,y)
     # df <- data.frame(x=names.arg, y=x$eig[, 1])
      df$x <- factor(df$x,levels=unique(df$x))
      bplot<- ggplot(data=df, aes(x=x, y=y)) + ggtitle(titleE)+ 
        labs(x = NULL)+ labs(y = NULL)+
        geom_bar(stat="identity", width=0.5, fill=col)+  
        theme_minimal()+
        theme(plot.title = element_text(hjust = 0.5))
      if(graph.type=="ggplot") {print(bplot); return(bplot)}

   #   if(graph.type == "interact") {
       my.text <-  paste(rownames(x$eig), "<br>", "Eigenvalue",round(x$eig[,1],4),
                    "<br>", "Pct.Variance", round(x$eig[,2],4),
                    "<br>", "AcumPct.Variance", round(x$eig[,3],4))
       bplot <- plotly::plotly_build(bplot)
       bplot$x$data[[1]]$text <- my.text
       
       fplot <- function(x) {
        #  default_opts <- callr::r(function(){options()}); options(default_opts)
         oo.viewer <- getOption("viewer")
         on.exit({ # print("Restoring viewer...")
           Sys.sleep(0)
           options(viewer = oo.viewer)
         }, add = TRUE)
         options(viewer = NULL)
         print(x)
       }
       if(new.plot) fplot(bplot) else print(bplot)
return(bplot)
  }

  
 stemp <- c(selDoc, selWord,selSeg,selDocSup,selWordSup,quanti.sup,quali.sup)
 stemp <- unique(stemp)
 if(!is.null(stemp)) {
   cat("When eigen=TRUE other elements cannot be plotted at the same time")
 return("")
                     }
} else {
invisib <- c("quali.sup")


if(!is.null(quanti.sup))
 {
  if (is.null(x$quanti.sup)) 
    stop("No quantitative supplementary variables in LexCA")
  qsn <- rownames(x$quanti.sup$coord)
  
  if(is.numeric(quanti.sup))   quanti.sup <-  na.omit(qsn[quanti.sup])
  sel1 <- selection(quanti.sup,x$quanti.sup,"Quanti",axes[1],axes[2])
 #  sel1  <- which(sel1 %in% qsn)  Changed version  1.5.4
  sel1  <-  which(qsn %in% sel1)

  
if(length(sel1)==0)
    stop("No quantitative supplementary variables in plot.LexCA")
 
objqs <- x$quanti.sup
drow <- unlist(dimnames(objqs[[1]])[1])
dcol <- unlist(dimnames(objqs[[1]])[2])
coord.qs <- matrix(x$quanti.sup[[1]],length(drow),length(dcol))
rownames(coord.qs) <- drow
colnames(coord.qs) <- dcol


coord.qs <- coord.qs[sel1,,drop=FALSE]
x$quanti.sup$coord <- coord.qs
 if(is.null(title)) titleS <- "Supplementary quantitative variables on the CA map"
 else titleS <- title


 res.quanti <-  FactoMineR::plot.CA(x, axes=axes, choix = c("quanti.sup"), col.quanti.sup =col.quanti.sup, title= titleS, 
         cex=cex,graph.type =graph.type)
 stemp <- c(selDoc, selWord,selSeg,selDocSup,selWordSup,quali.sup)

 stemp <- unique(stemp)
 if(!is.null(stemp)) {
   cat("When quanti.sup is not NULL other elements cannot be plotted at the same time")
 return("")
 }
 return(res.quanti)
} else {

if(!is.null(selDoc)) { 
 selDoc <- selection(selDoc,x$row,"Doc",axes[1],axes[2])
 if(length(selDoc)==0) selDoc <- NULL 
}

if(!is.null(selWord)) {
 selWord <- selection(selWord, x$col,"Word", axes[1],axes[2])
 if(length(selWord)==0) selWord <- NULL 
}


#====================================
if(!is.null(quali.sup)) {
  if(!is.null(x$quali.sup)) 
    quali.sup <- selection(quali.sup,x$quali.sup,"Quali", axes[1],axes[2])
    else stop("No categorical supplementary variables in LexCA object")}
if(!is.null(quali.sup)) {
  if(!is.null(x$var.agg)) 
   x$quali.sup$coord <- data.frame(x$quali.sup$coord[,,drop=FALSE])
  posic <- which(rownames(x$quali.sup$coord) %in% quali.sup)
   x$quali.sup$coord <- x$quali.sup$coord[posic,,drop=FALSE]
   x$quali.sup$coord <- x$quali.sup$coord[apply(!is.na(x$quali.sup$coord), 1, any),,drop=FALSE ] 

}

# Supplementary documents
if(!is.null(selDocSup)) 
  if(!is.null(x$row.sup)) 
    {selDocSup <- selection(selDocSup,x$row.sup,"Dsup", axes[1],axes[2])
     if(length(selDocSup)==0) selDocSup <- NULL
    } else {stop("No supplementary documents in LexCA object")}

# Qualitative variables
if(!is.null(quali.sup)) 
if(!is.null(x$quali.sup$coord)) {
 rdo <- which(rownames(x$quali.sup$coord) %in% rownames(x$row$coord))
 if(length(rdo)>0) rownames(x$quali.sup$coord)[rdo] <- paste0("_",rownames(x$quali.sup$coord)[rdo])

 if(!is.null(x$row.sup$coord)) {
 rdo <- which(rownames(x$row.sup$coord) %in% rownames(x$row$coord))
 if(length(rdo)>0) rownames(x$row.sup$coord)[rdo] <- paste0("_",rownames(x$row.sup$coord)[rdo])}
 ifelse(!is.null(x$row.sup), nrowsup <- nrow(x$row.sup$coord),nrowsup <-0)

if(is.null(selDocSup)){ 
    x$row.sup$coord <- x$quali.sup$coord
     selDocSup <- rownames(x$quali.sup$coord)
    col.doc.sup <- col.quali.sup
  } else {
 colnames(x$quali.sup$coord) <- colnames(x$row$coord)
 x$row.sup$coord <- rbind(x$row.sup$coord, x$quali.sup$coord)
 selDocSup <- c(selDocSup, rownames(x$quali.sup$coord))
 col.doc.sup <- c(rep(col.doc.sup,nrowsup), rep(col.quali.sup,nrow(x$quali.sup$coord)))
}

 if(length(selDocSup)==0) selDocSup <- NULL
}
if(is.null(selDocSup)) invisib <- c(invisib,"row.sup")
else selDoc <- c(selDoc,selDocSup)


#============================================================
# Supplementary Segments & words
#====================================
nsegm <-0; nwordsup <-0
 if(!is.null(x$segment$coord)) nsegm <- nrow(x$segment$coord)
 if(!is.null(x$col.sup$coord)) nwordsup <- (nrow(x$col.sup$coord))-nsegm


if(!is.null(selSeg)) {
  if(nsegm==0) stop("No segments in LexCA object")
   selSeg <- selection(selSeg,x$segment,"Seg", axes[1],axes[2])
   nsegm <- length(selSeg)
   if(nsegm==0) selSeg <- NULL } else {nsegm<-0}


# Supplementary words
if(!is.null(selWordSup)) {
  if(nwordsup==0) stop("No supplementary words in LexCA object")
    tmp <- x$col.sup
    tmp$coord <- tmp$coord[1:nwordsup,,drop=FALSE]
    tmp$cos2 <- tmp$cos2[1:nwordsup,,drop=FALSE]
    selWordSup <- selection(selWordSup,tmp,"Wsup", axes[1],axes[2])
    nwordsup <- length(selWordSup)
    if(nwordsup==0) selWordSup <- NULL} else {nwordsup<-0} 

if(!is.null(x$col.sup$coord)) x$col.sup$coord <- x$col.sup$coord[selWordSup,,drop=FALSE]
if(!is.null(x$segment$coord)) x$segment$coord <- x$segment$coord[selSeg,,drop=FALSE]

if(nsegm>0) {
 if(nwordsup==0) x$col.sup$coord <- x$segment$coord 
  else x$col.sup$coord <- rbind(x$col.sup$coord,x$segment$coord)
}

 col.word.sup <- c(rep(col.word.sup,nwordsup), rep(col.seg,nsegm))

 if(length(col.word.sup)==0) { col.sup <- NULL; invisib <- c(invisib,"col.sup")}
 invisib <- c(invisib,"quanti.sup")
 if(is.null(selDoc))  invisib <- c(invisib,"row")
 if(is.null(selWord))  invisib <- c(invisib,"col")

if(!is.null(x$row.sup$coord)) colnames(x$row.sup$coord) <- colnames(x$row$coord)
if(!is.null(selWord)){
 if(!is.null(selWordSup)) selWord <- c(selWord, selWordSup)
 if(!is.null(selSeg)) selWord <- c(selWord, selSeg)
}


invisib <- unique(invisib)
if(length(invisib)==6) stop("No selected elements to plot")


if(eigen!=TRUE)
if(new.plot==TRUE){if(dev.interactive()) dev.new()}
 plot.CA(x, axes=axes, invisible=invisib,
    choix=c("CA"), title= title, cex=cex, selectCol=selWord, selectRow=selDoc, 
    xlim=xlim, ylim=ylim, shadowtext=shadowtext,habillage=habillage, unselect=unselect,
    autoLab=autoLab, col.row=col.doc, col.row.sup=col.doc.sup, label=label, 
    col.col=col.word, col.col.sup=col.word.sup, graph.type =graph.type)
}}}
