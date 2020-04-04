#' @importFrom graphics lines
#' @export
 ellipseLexCA <- function (object, selWord="ALL",selDoc="ALL", 
    nbsample = 100, level.conf = 0.95, axes = c(1, 2), ncp=NULL, xlim = NULL, ylim = NULL, title=NULL,
    col.doc = "blue", col.word = "red", col.doc.ell = col.doc, 
    col.word.ell = col.word, cex=1) 
{
    if (!inherits(object, "LexCA"))  stop("Object should be LexCA class")
   graph.type <- "classic"
   options(stringsAsFactors = FALSE)
   
if(is.null(ncp)) ncp <- ncol(object$col$coord) else
ncp <- min(ncol(object$col$coord),ncp)
if(max(axes)> ncp) stop("Axis=" , max(axes), " is bigger than ncp=", ncp, ". Change ncp value in LexCA" )
   
   
selectionDW <- function(sel1, xobj,bType, axx, axy)
{
if(is.null(sel1)) {  sel1<-c()
                     return(sel1)}
xx <- ""
if(length(sel1)==1){
     if(sel1=="ALL") sel1 <- c(1:dim(xobj$coord)[1])
                   }
if(length(sel1)==1) {
   if(sel1=="meta") sel1 <- "meta 3"
   xx <- gregexpr(pattern =' ',sel1)[[1]][1]-1
   xx <- substr(sel1, 1, xx)
                     }
if(xx=="coord" | xx=="cos2" | xx=="contrib" | xx=="meta" ) 
 { 
  nc <- nchar(sel1)
  if(xx=="coord") {
 # Selection by coordinates"
 sel1 <- as.numeric(substr(sel1, 7, nc))
 dft <- data.frame(xobj$coord[,c(axx,axy),drop=FALSE])
 fval <- apply(dft, 1, function(x) max(abs(x)))
 ordmax <- rank(fval)
 posic <- which(ordmax > (length(fval)-sel1))
 sel1 <- rownames(xobj$coord)[posic]
                   }
  if(xx=="cos2") {
 # Selection by cos2
 sel1 <- as.numeric(substr(sel1, 5, nc))
 dft <- data.frame(xobj$cos2[,c(axx,axy),drop=FALSE])
 fval <- apply(dft, 1, function(x) sum(x))
 posic <- which(fval>= sel1)
 sel1 <- rownames(xobj$cos2)[posic]
                  }
  if(xx=="contrib") {
 # Selection by contrib
 sel1 <- as.numeric(substr(sel1, 9, nc))
 dft <- data.frame(xobj$contrib[,c(axx,axy),drop=FALSE])
 fval <- apply(dft, 1, function(x) max(abs(x)))
 posic <- which(fval>= sel1)
 sel1 <- rownames(xobj$contrib)[posic]
  }
 if(xx=="meta") {
 sel1 <- as.numeric(substr(sel1, 5, nc))
 sKeys <-  rownames(xobj$coord)[which(xobj$contrib[,axx]>mean(xobj$contrib[,axx])*sel1)]
 sKeys <-  c(sKeys,rownames(xobj$coord)[which(xobj$contrib[,axy]>mean(xobj$contrib[,axy])*sel1)])
 sel1 <- unique(sKeys)
}
 } else {
 if(is.character(sel1)) sel1 <- which(rownames(xobj$coord) %in% sel1)
 sel1 <- rownames(xobj$coord)[sel1]
 sel1  <- sel1[!is.na(sel1)]
}
 return(sel1)
}
# Final functions

X = object$call$X
if(!is.null(selDoc)) selDoc <- selectionDW(selDoc,object$row,"Doc",axes[1],axes[2])
if(!is.null(selWord)) selWord <- selectionDW(selWord, object$col,"Word",axes[1],axes[2])
if ( (length(selDoc) + length(selWord))  == 0)
    stop("No elements are selected to be plotted;", "\n",
         "  verify the expression of both arguments selWord and/or selDoc")

    concCol = X
    concRow = X
    sampcol=0
    samprow=0
    method <- "multinomial" 
    proba = unlist(c(X))
    N = sum(proba)
    proba = proba/N
    aa = rmultinom(nbsample, size = N, prob = proba)
    for (i in 1:nbsample) {
            aux = matrix(aa[, i, drop = FALSE], nrow = nrow(X))
            dimnames(aux) = dimnames(X)
            if(!is.null(selWord))
               {
               aaa<-apply(aux,2,sum)  
               if (!(0 %in% aaa)) 
                      {
                      concCol = cbind.data.frame(concCol, aux)
                      } else {
                      sampcol<-sampcol+1 
                      } 
               }
            if(!is.null(selDoc)) 
               {
               bbb<-apply(aux,1,sum)  
               if (!(0 %in% bbb)) 
                      {
                      concRow = rbind.data.frame(concRow, aux)
                      } else {
                      samprow<-samprow+1 
                      } 
                }
        }

    
  
  if ((!is.null(selWord)) &  (sampcol/nbsample)>0.1 )  {
    stop("over 10% of the replicated samples present words with null frequency","\n",
          "  it is advisable to increase the threshold on the word frequency in TextData") }

  if ( (!is.null(selDoc)) & (samprow/nbsample)>0.1 )   {
    stop("over 10% of the replicated samples present documents with null length","\n",
         "  it is advisable to increase the threshold on the word frequency in TextData") }

    vdvword<-rep("transparent",ncol(X))
    vdvword[which(colnames(X) %in% selWord)]<-col.word.ell
    vdvdoc<-rep("transparent",nrow(X))
    vdvdoc[which(rownames(X) %in% selDoc)]<-col.doc.ell

    
    
    
    
    Xcoord.ellipse <- function (coord.simul, centre = NULL, axes = c(1, 2), level.conf = 0.95, 
              npoint = 100, bary = FALSE) 
    {
  
      
      Xellipse <- function(x, scale = c(1, 1), centre = c(0, 0), 
                          level = 0.95, t = sqrt(qchisq(level, 2)), which = c(1, 2), npoints = 100) {
        names <- c("x", "y")
        if (is.matrix(x)) {
          xind <- which[1]
          yind <- which[2]
          r <- x[xind, yind]
          if (missing(scale)) {
            scale <- sqrt(c(x[xind, xind], x[yind, yind]))
            if (scale[1] > 0) 
              r <- r/scale[1]
            if (scale[2] > 0) 
              r <- r/scale[2]
          }
          if (!is.null(dimnames(x)[[1]])) 
            names <- dimnames(x)[[1]][c(xind, yind)]
        }
        else r <- x
        r <- min(max(r, -1), 1)
        d <- acos(r)
        a <- seq(0, 2 * pi, len = npoints)
        matrix(c(t * scale[1] * cos(a + d/2) + centre[1], t * 
                   scale[2] * cos(a - d/2) + centre[2]), npoints, 2, 
               dimnames = list(NULL, names))
      }
      nbre.fact <- nlevels(coord.simul[, 1])
      res <- label <- NULL
      lev <- levels(coord.simul[, 1])
    
      
      for (f in 1:nbre.fact) {
        x <- coord.simul[which(coord.simul[, 1] == lev[f]), axes[1] + 1]
        y <- coord.simul[which(coord.simul[, 1] == lev[f]), axes[2] + 1]
        if (is.null(centre)) 
          center <- c(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE))
        else {
          if (ncol(coord.simul) != ncol(centre)) 
            stop("ncol de centre incorrect")
          if (!all.equal(lev, levels(centre[, 1]))) 
            stop("Levels of centre are not corrects")
          center <- as.numeric(centre[which(centre[, 1] ==levels(centre[, 1])[f]), 
                                      c(axes[1] + 1, axes[2] + 1)])
        }
        
       # tab <- data.frame(x = x, y = y)
        tab <- data.frame(x = x, y = y, stringsAsFactors = TRUE)
   
        if (nrow(tab) > 1) 
          mat.cov <- cov(tab)
        else mat.cov <- matrix(0, 4, ncol = 2)
        if (bary) 
          mat.cov = mat.cov/nrow(tab)
      elli.tmp <- Xellipse(mat.cov, centre = center, 
                                     level = level.conf, npoints = npoint)
        res <- rbind(res, elli.tmp)
      }
      label <- factor(rep(lev, each = npoint), levels = lev)
      result <- data.frame(facteur = label, res)
      colnames(result)[1] = "facteur"
      colnames(result) <- colnames(coord.simul)[c(1, axes + 1)]
      return(list(res = result, call = npoint))
    }
    

    if (!is.null(selWord)) {
        colCA <- FactoMineR::CA(concCol, col.sup = (ncol(X) + 1):ncol(concCol), graph = FALSE, ncp=ncp)
        aux3 <- colCA$col.sup$coord[, axes]
        rownames(aux3) <- paste("r", 1:nrow(aux3), sep = "")
        aux1 <- cbind.data.frame(label=paste("word", 1:ncol(X), sep = ""), aux3)
        aux1$label <- as.factor(aux1$label)
        ellCol <- Xcoord.ellipse(aux1, level.conf = level.conf)$res
    }
    if (!is.null(selDoc)) {
        rowCA = FactoMineR::CA(concRow, row.sup = (nrow(X) + 1):nrow(concRow), graph = FALSE, ncp=ncp)
        aux2 = cbind.data.frame(label=paste("doc", 1:nrow(X), sep = ""), 
            rowCA$row.sup$coord[, axes])
        aux2$label <- as.factor(aux2$label)
        ellRow = Xcoord.ellipse(aux2, level.conf = level.conf)$res
    }
    if (is.null(xlim)) {
        if (  (!is.null(selWord)) & (!is.null(selDoc))  ) 
            xlim <- c(min(ellCol[, 2], ellRow[, 2]), max(ellCol[, 2], ellRow[, 2]))
        else {
            if (!is.null(selWord))
                xlim <- c(min(ellCol[, 2]), max(ellCol[, 2]))
            if (!is.null(selDoc)) 
                xlim <- c(min(ellRow[, 2]), max(ellRow[, 2]))
        }
    }
    if (is.null(ylim)) {
        if (  (!is.null(selWord)) & (!is.null(selDoc))  ) 
            ylim <- c(min(ellCol[, 3], ellRow[, 3]), max(ellCol[, 3], ellRow[, 3]))
        else {
            if (!is.null(selWord)) 
                ylim <- c(min(ellCol[, 3]), max(ellCol[, 3]))
            if (!is.null(selDoc))
                ylim <- c(min(ellRow[, 3]), max(ellRow[, 3]))
        }
    }
    if (  (!is.null(selWord)) & (!is.null(selDoc))  )
 plot.LexCA(object, axes = axes, selWord=selWord,selDoc=selDoc,xlim = xlim, ylim = ylim, col.doc = col.doc, col.word = col.word,
            cex=cex,title=title, graph.type =graph.type, plot.new=TRUE)
    else {
          if (!is.null(selWord))  
         plot.LexCA(object, axes = axes, selWord=selWord,selDoc=NULL,xlim = xlim, ylim = ylim, col.word = col.word,cex=cex,title=title,
                    graph.type =graph.type, plot.new=TRUE) 
          if (!is.null(selDoc)) 
         plot.LexCA(object, axes = axes, selWord=NULL,selDoc=selDoc,xlim = xlim, ylim = ylim,  col.doc = col.doc,,cex=cex,title=title,
                    graph.type =graph.type, plot.new=TRUE)
         }

    if (!is.null(selDoc)) {
        lev <- paste("doc", 1:nlevels(ellRow[, 1]), sep = "")
        for (e in 1:nlevels(ellRow[, 1])) {
            data.elli <- ellRow[ellRow[, 1] == lev[e], -1]
            lines(x = data.elli[, 1], y = data.elli[, 2], col = vdvdoc[e])
        }
    }
    if (!is.null(selWord)) {
        lev <- paste("word", 1:nlevels(ellCol[, 1]), sep = "")
        for (e in 1:nlevels(ellCol[, 1])) {
            data.elli <- ellCol[ellCol[, 1] == lev[e], -1]
            lines(x = data.elli[, 1], y = data.elli[, 2], col = vdvword[e])
        }
    }
}

