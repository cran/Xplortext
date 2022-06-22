#' @import stringr
#' @export
LexChar <- function(object,  proba = 0.05, maxDocs=20, maxCharDoc = 10, maxPrnDoc=100, marg.doc="before",
                    context=NULL, correct=TRUE,  nbsample=500, seed=12345, ...)
 {
  set.seed(seed)
  options(stringsAsFactors = FALSE)
  # maxDocs: Maximum number of input documents or aggregate categories
  # maxCharDoc: Maximum number of characteristic documents to show
  # correct = TRUE and Correction for pvalue of Characteristic words (FactoMineR, Agresti...)
  # marg.doc; after/before/before.RW

  # Version 1.4.1 used context.sup, version 1.4.2. use context argument
  varnames<- lapply(substitute(list(...))[-1], deparse)
  if(!is.null(varnames$context.sup)) {
    context <- gsub("[[:punct:]]", "", varnames$context.sup)   # Only for Compatibility version 1.4.1
    warning("Xplortext Versions > 1.4.1 use context, no context.sup argument") 
  }
  ###################################################################
  
  # Checking object type, proba limits and maxCharDoc values
  if (!inherits(object,"TextData") & !inherits(object,"DocumentTermMatrix") & !inherits(object,"matrix")
      & !inherits(object,"data.frame"))
    stop("Object should be TextData, DocumentTermMatrix, matrix or data frame")
  if(proba<0|proba>1) stop("proba should be between 0 and 1")
  
  if(is.null(maxCharDoc) |  maxCharDoc <0) maxCharDoc <- 0 

  
  ###################################################################
  ## Functions
  ############### descfreq_New function #############################################
  descfreq_NewChar <- function (data, proba = 0.05, marge.li, marge.col) 				
  {
    lab.sauv <- lab <- colnames(data)				
    for (i in 1:length(lab)) lab[i] = gsub(" ", ".", lab[i])				
    colnames(data) <- lab
    
    old.warn = options("warn")			
    options(warn = -1)			# suppress warnings globally

    nom <- tri <- structure(vector(mode = "list", length = nrow(data)), names = rownames(data))		
    if(is.data.frame(marge.li)) marge.li <- as.vector(marge.li[,1]) else marge.li <- as.vector(marge.li)
    if(is.data.frame(marge.col)) marge.col <- as.vector(marge.col[,1]) else marge.col <- as.vector(marge.col)
    

    SumTot<-sum(marge.li)				
    nr <- nrow(data)  				
    nc <- ncol(data) 				

    for (j in 1:nr) {
      aux3 <- marge.li[j]/SumTot				# % Occurrences before or after
      
      for (k in 1:nc) {
        aux2 <- data[j, k]/marge.col[k]
        daux <- dhyper(data[j, k], marge.col[k], SumTot - marge.col[k], marge.li[j])
        if(correct==FALSE) daux <- daux*2  
        #  P[X <= OCURR.W-1]
        A <- phyper(data[j, k]-1, marge.col[k], 				
                    SumTot - marge.col[k], marge.li[j],lower.tail = TRUE)
        #  P[X > OCURR.W]
        B <- phyper(data[j, k], marge.col[k], SumTot - marge.col[k], marge.li[j],lower.tail = FALSE)
        aux4 <- min(A,B)*2+daux
        if (aux4 <= proba) {
          aux5 = (1 - 2 * as.integer(aux2 > aux3)) * qnorm(aux4/2)  # vtest			
          aux1 = data[j, k]/marge.li[j]				
          tri[[j]] = rbind(tri[[j]], c(aux1 * 100, sum(marge.col[k])/SumTot * 				
                                         100, data[j, k], marge.col[k], aux4, aux5))	
          nom[[j]] = rbind(nom[[j]], c(colnames(data)[k], colnames(data)))				
        }
      }
    }
    
    
    for (j in 1:nr) {				
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
    res = tri				
    options(old.warn)				
    class(res) <- c("descfreq", "list ")				
    return(res)				
  }		  # End   descfreq_New
 
   ##############  Function to compute Vocab$quali$stats  
  chi.quali<- function(X, QL)  {
    list.Agg <- lapply(seq_along(QL),FUN=function(i) t(t(as.matrix(X))%*% as.matrix(QL[[i]])))
    
    Nq <- sum(X)
    
    dfq<-NULL
    old.warn = options("warn")				
    options(warn = -1)			# suppress warnings globally
    for(i in 1:length(list.Agg)) {
      ch <- chisq.test(list.Agg[[i]], correct=FALSE)[1:3]
      ph2 <- ch$statistic/Nq
      ph <- sqrt(ph2)
      minrowcol <- min(ncol(list.Agg[[i]]),nrow(list.Agg[[i]]))-1
      VCr <-sqrt(ph2/ minrowcol )
      d <- data.frame(ch,ph,VCr)
      dfq <- rbind(dfq,d)
    }
    colnames(dfq) <- c("Chi.Squared", "df", "p.value", "phi", "Cramer's V")
    rownames(dfq) <-  names(QL)
    options(old.warn)	
    return(dfq)
  }  
  chi.quali.single <- function(X)  {
    old.warn = options("warn")				
    options(warn = -1)			# suppress warnings globally
    ch <- chisq.test(X, correct=FALSE)[1:3]
    ph2 <- ch$statistic/sum(X)
    ph <- sqrt(ph2)
    minrowcol <- min(ncol(X),nrow(X))-1
    VCr <-sqrt(ph2/ minrowcol )
    d <- data.frame(ch,ph,VCr)    
    colnames(d) <- c("Chi.Squared", "df", "p.value", "phi", "Cramer's V")
    rownames(d) <-  "LexicalTable"
    options(old.warn)	
    return(d)
  }

  ##############  Function to compute Vocab$quanti$stats    
  mean.p <- function(V,poids) res<-sum(V*poids,na.rm=TRUE)/sum(poids[!is.na(V)])
  var.p <- function(V,poids) res<-sum(V^2*poids,na.rm=TRUE)/sum(poids[!is.na(V)])
  
  
  vocabQuanti <- function(vdt,vX, vrow.INIT ) {
    vrow.INIT <- vrow.INIT[,1]
    vX <- as.data.frame(vX)

    for(i in 1:ncol(vX)){
      if(any(is.na(vX[,i]))) warning("\n", names(vX)[i], " variable has missing values. They will be replaced by the mean\n") 
      vX[is.na(vX[,i]), i] <- mean(vX[,i], na.rm = TRUE)
    }
    
    # Weight documents
    Wi <- vrow.INIT/sum(vrow.INIT)
    col.INIT <- colSums(vdt)
    y <- sum(col.INIT)
    #  Weighted average of quantitative variables  (1 x 2 variables)
    
    Aver.X <- apply(vX,2,mean.p,Wi)
    Var.X<-apply(sweep(vX,2,Aver.X,"-"),2,var.p,Wi)

    #Average of quantitative variables for each word
    Wij.cj <- as.matrix(sweep(vdt,2,col.INIT ,"/"))
    
    MeanWord <- t(vX) %*% Wij.cj 

    # Variance of words
    coef <- as.matrix((y/col.INIT-1)/(y-1))
    Var.Y <- coef %*% Var.X
    sc.inter <-   apply(sweep(t(sweep(t(MeanWord), 2, Aver.X,"-")^2),2,col.INIT,"*"),1,sum)/y
    pct.expl <- sc.inter/Var.X
    Desv.Y <- sqrt(Var.Y)
    dj <- t(sweep(MeanWord,1,Aver.X,"-")) / Desv.Y
    
    
    ### Permutations
    # Doing list with vocabulary elements
    nWords <- ncol(vdt)                                                       
    nDocs <- nrow(vdt)                                                          
    
    relcq.perm <-vector(mode='list',length=nWords)
    for (i in 1:nWords) relcq.perm[[i]]<-matrix(nrow=nbsample,ncol=ncol(vX),0) 
    relcq.SS.perm <-vector(mode='list',length=ncol(vX))
    for (i in 1:ncol(vX)) relcq.SS.perm[[i]]<- matrix(nrow=nbsample,ncol=1,0) 	 
    

    for (i in 1:nbsample){                                                      
      X.perm <- vX[sample(1:nDocs),,drop=FALSE]
      Aver.X.perm<-apply(X.perm,2,mean.p,Wi)
      MeanWord.perm <- t(X.perm) %*% Wij.cj
      Var.X.perm<-apply(sweep(X.perm,2,Aver.X.perm,"-"),2,var.p,Wi)
      Var.Y.perm <- coef %*% Var.X.perm
      # dj words x variables
      dj.perm <- t(sweep(MeanWord.perm,1,Aver.X.perm,"-")) / sqrt(Var.Y.perm)   
      

      for (iw in 1:nWords) 
        relcq.perm[[iw]][i,]<-dj.perm[iw,]    
      sc.inter.perm <-   apply(sweep(t(sweep(t(MeanWord.perm), 2, Aver.X.perm,"-")^2),2,col.INIT,"*"),1,sum)/y
      pct.expl.perm <- sc.inter.perm/Var.X.perm
      
      for (iw in 1:ncol(vX)) relcq.SS.perm[[iw]][i,] <- pct.expl.perm[iw] # /  SC.X.perm[iw]
    }
  
    
    #  names(relcq.perm) <- colnames(vdt)
    a <- matrix(nrow=ncol(vdt),ncol=ncol(vX),0)
    b <- matrix(nrow=1,ncol=ncol(vX),0)
    

    for (k in 1:ncol(vX)){
      for (i in 1:ncol(vdt)){
        if (dj[i,k]>0) a[i,k] <- sum(relcq.perm[[i]][,k]>dj[i,k])/nbsample
        if (dj[i,k]<0) a[i,k] <- sum(relcq.perm[[i]][,k]<dj[i,k])/nbsample
      }
      b[1,k] <- sum(relcq.SS.perm[[k]]> pct.expl[k]) /nbsample
    }
    

    
    
    res.mat <- matrix(nrow=ncol(vX),ncol=4,0) 	
    rownames(res.mat)<-colnames(vX)
    colnames(res.mat )<- c("GlobalAverage","AverageWord","Differ.","pvalue")
    #    names(relcq.perm) <- colnames(vdt)
    for (i in 1:nWords) {
      relcq.perm[[i]]<-res.mat
      relcq.perm[[i]][,1] <- t(Aver.X)
      relcq.perm[[i]][,2] <- MeanWord[,i]
      relcq.perm[[i]][,3] <-  relcq.perm[[i]][,2]- relcq.perm[[i]][,1]
      relcq.perm[[i]][,4] <-  t(a[i,])
      names(relcq.perm)[[i]] <- colnames(vdt)[i]
      relcq.perm[[i]] <- subset(relcq.perm[[i]], proba > relcq.perm[[i]][,4])
    }
    
    relcq.perm <- relcq.perm[sapply(relcq.perm, function(x) ifelse(nrow(x)==0,F,T))]
    res.mat.SS <- matrix(nrow=ncol(vX),ncol=4,0) 	
    rownames(res.mat.SS)<-colnames(vX)
    colnames(res.mat.SS)<- c("TotSumSquares","BetweenSquares","%Explained","pvalue")
    for (i in 1:ncol(vX)) {
      res.mat.SS[i,1] <- Var.X[i]
      res.mat.SS[i,2] <- sc.inter[i]
      res.mat.SS[i,3] <- pct.expl[i]
      res.mat.SS[i,4] <- b[1,i]
    }
    
    quanti <- list(CharWord=relcq.perm, stats=res.mat.SS)
    return(quanti)
  }
  ##### End of functions #############################################
  
  

  #### step 1. Selecting the type of object ##############
  ### Verifying if it is a TextDataObject and aggregated TextData
  bTextData <- ifelse(inherits(object,"TextData"),TRUE,FALSE)
  ###  To know if it is an aggregate analysis for an TextData object
  # If it is an aggregated table
  bvaragg <- FALSE
  if(inherits(object,"TextData"))
    if(!is.null(object$info)) 
      if(object$info$name.var.agg[[1]]!="") bvaragg <- TRUE 
  
  
  #=================================================================	
  # Computing DocTerm and row.INIT
  #### step 2. Detecting names contextual variables ##############
  if(bTextData) {
    # Is bTextData
   # df.qual<- strquali <- strquanti <- NULL
    strquali <- strquanti <- NULL
    if(!is.null(context))
    if(length(context)==1) if(context=="ALL") {
      if(bvaragg) {
        context  <- colnames(object$SourceTerm.qual)
        context  <- c(context,colnames(object$SourceTerm.quant))  
      } else {
        context  <- colnames(object$context$quali)
        context  <- c(context,colnames(object$context$quanti))
      } }
    
    
    # context, quanti and quali names of variables
    strquali<- colnames(object$SourceTerm.qual)[which(colnames(object$SourceTerm.qual) %in% context)]
    strquanti<- colnames(object$SourceTerm.quant)[which(colnames(object$SourceTerm.quant) %in% context)]
    if(!is.null(strquali)) if(length(strquali)==0) strquali <- NULL
    if(!is.null(strquanti)) if(length(strquanti)==0) strquanti <- NULL
    
  } else {
    # Is not bTextData
   # df.context.quali <- df.context.quanti <-NULL
    if("DocumentTermMatrix" %in% class(object)) object <- as.matrix(object)
    if(is.matrix(object))      object <- as.data.frame(object)
    if(!is.data.frame(object)) stop("Error: object must be a dataframe")

    if(marg.doc!="after") {
      options(warn=0)
      warning("Only marg.doc==after is allowed ; is changed to after")
      marg.doc <- "after"
    }

    # DocTerm <- object
    sel.context <- colnames(object)[which(colnames(object) %in% context)]
    
    strquali <- strquanti <- NULL
    if(length(sel.context)!=0) {
      df.context <- object[,sel.context,drop=FALSE]
      context.type <- sapply(df.context, class)
      strquali <- names(context.type[ which(context.type %in% c("factor","character", "logical", "Date"))])
      strquanti <- names(context.type[!names(context.type) %in% strquali])
    }
  }
  

  
  #=================================================================	
  
  ##### step 3. Building object to descfreq_NewChar function ##############
  if(bTextData) { # There is TextData object
    # 1. Computing DocTerm
    DocTerm <- as.matrix(object$DocTerm)
    if(marg.doc=="after") {
      DocTerm <- DocTerm[rowSums(DocTerm)!=0,]
      row.INIT <- data.frame(Occurrences.after=rowSums(DocTerm))
    }
    if(marg.doc=="before") {
      DocTerm <- DocTerm[rowSums(DocTerm)!=0,]
      row.INIT <- as.data.frame(object$rowINIT)        # Frequency Occurrences.before
      # before0 for remove rows with 0 frequency. Occurrences.before
      # row.INIT <- row.INIT[row.INIT$Occurrences.before>0,,drop=FALSE] 
      row.INIT <- row.INIT[rownames(DocTerm),,drop=FALSE]
    }
    
    if(marg.doc=="before.RW") {
      row.INIT<- as.data.frame(object$summDoc[, "Occurrences.before",drop=FALSE])
      rownames(row.INIT) <- object$summDoc$DocName
      row.INIT <- row.INIT[row.INIT$Occurrences.before!=0,,drop=FALSE] 
      
      # NoNullBefore rows with no null documents before but after threshold
      DT2 <- DocTerm
      NoNullBefore <- row.INIT[!rownames(row.INIT) %in% rownames(DocTerm),,drop=FALSE]
      
      if(nrow(NoNullBefore)>0){
        DT2 <- cbind(DT2,as.data.frame(row.INIT[rownames(DocTerm),]))
        nrDT2 <- nrow(DT2)
        strnames <- rownames(NoNullBefore)
        # Add NoNullBefore rows to dataframe and fill them with 0
        DT2[( nrDT2+1):( nrDT2+  nrow(NoNullBefore)),] <- 0
        # Add rownames to new rows from strnames
        rownames(DT2)[(nrDT2+1):(nrDT2+nrow(NoNullBefore))]  <-  strnames  
        DT2[(nrDT2+1):( nrDT2+  nrow(NoNullBefore)),ncol(DT2)] <- NoNullBefore[,1]
        DT2 <- DT2[rownames(row.INIT),]
      } 
      else {
        DT2 <- cbind(DT2,as.data.frame(row.INIT[rownames(DocTerm),]))
      }# End if(nrow(NoNullBefore)>0)
      
      freq.after <- apply(DT2[,c(1:(ncol(DT2)-1))],1,sum)
      freq.before <- DT2[,ncol(DT2)]      
      # Adding a column with the name of RemovedWords
      DT2 <- cbind(DT2,as.data.frame(freq.before-freq.after))
      colnames(DT2)[ncol(DT2)] <- "RemovedWords"
      DT2 <- DT2[,-(ncol(DT2)-1)]
      #  colnames(DT2)[ncol(DT2)] <- "RemovedWords"
      DocTerm <- DT2
    }  # End before.RW
    # DocTerm can have some document/aggregate document with margin zero
    pos.0 <- which(rowSums(DocTerm)==0)
    if(length(pos.0)!=0) {
      # Remove empty documents or aggregate documents 
      DocTerm <- DocTerm[!rownames(DocTerm) %in% names(pos.0),,drop=FALSE]
      row.INIT <- row.INIT[!rownames(row.INIT) %in% names(pos.0),,drop=FALSE]
    }
    
  } else {
    # Is not bTextData

    DocTerm <- object
    if(!is.null(strquali))  DocTerm <- DocTerm[,!colnames(DocTerm) %in% strquali,drop=FALSE] 
    if(!is.null(strquanti))  DocTerm <- DocTerm[,!colnames(DocTerm) %in% strquanti,drop=FALSE]
    
    context.type <- sapply(DocTerm, class)

    strquali.rem <- names(context.type[ which(context.type %in% c("factor","character", "logical", "Date"))])
    DocTerm <- DocTerm[,!colnames(DocTerm) %in% strquali.rem,drop=FALSE]
    if(!is.null(strquali)) df.context.quali <- object[,colnames(object) %in% strquali,drop=FALSE]
    if(!is.null(strquanti)) df.context.quanti <- object[,colnames(object) %in% strquanti,drop=FALSE]
    pos.0 <- which(rowSums(DocTerm)==0)

    if(length(pos.0)!=0) {
      # Remove empty documents of non TextData object 
      DocTerm <- DocTerm[!rownames(DocTerm) %in% names(pos.0),,drop=FALSE]
      if(!is.null(strquali))  df.context.quali<- df.context.quali[!rownames(df.context.quali) %in% names(pos.0),,drop=FALSE] 
      if(!is.null(strquanti))  df.context.quanti<- df.context.quanti[!rownames(df.context.quanti) %in% names(pos.0),,drop=FALSE]  
    }

    row.INIT <- data.frame(rowSums(DocTerm))
#    stop("Is not bTextData")
  }
  col.INIT <- as.vector(colSums(DocTerm))
  row.INIT.B <- row.INIT

  #--------------------------------------
  #  Computing Lexical Table: Chi.Squared  df   p.value       phi Cramer's V. No weigthing is applied
  n.words<- ncol(DocTerm)
  d.single <- chi.quali.single(DocTerm[rowSums(DocTerm)>0,c(1:n.words)]) 
  #--------------------------------------------------------------------------------------------
  resCharWord <- descfreq_NewChar(DocTerm, proba = proba, row.INIT, col.INIT) 	
  # Return words for active categories
  res <- list(CharWord=resCharWord, stats=d.single)
  # Chi square not depend on row.INIT
  ###########################################################################################################



  
  #=================================================================	  
  ##### step 4. Quali contextual variables ##############
  if(length(strquali)==0) strquali <- NULL

  
  if(!is.null(strquali)){
   if(bTextData) {  # Is bTextData
    df.context.quali <- object$SourceTerm.qual[,strquali,drop=FALSE]
    
       if(bvaragg) {   ## si bvaragg add aggr variable
         df.aggr <- object$SourceTerm.var.agg
         df.context.quali <- cbind(df.aggr,df.context.quali)
         }  
    
    if(ncol(df.context.quali)>0) {
      df.context.quali$new.Cat <- ""
      df.context.quali$new.Cat  <- paste0(df.context.quali$new.Cat, paste0(colnames(df.context.quali)[1],"."), df.context.quali[,1])
    
      if(ncol(df.context.quali)>2)
      for(i in 2:(ncol(df.context.quali)-1)) {
        df.context.quali$new.Cat  <- paste0(df.context.quali$new.Cat, 
                                            paste0("@",colnames(df.context.quali)[i],"."), df.context.quali[,i])
      } }

    if(!is.null(strquanti)) {
       X.quanti.agg <- data.frame(object$SourceTerm.quant[,strquanti,drop=FALSE],"new.Cat"=df.context.quali$new.Cat)
    } # End if(!is.null(strquanti))
    
    
        if(bvaragg) {
      DT3 <- as.matrix(object$SourceTerm.dtm)   # original documents x total vocabulary
      row.INIT <- data.frame("row.INIT"=rowSums(DT3))
      DT4 <-  DT3[,colnames(DT3) %in% colnames(DocTerm),drop=FALSE]
      
      if(marg.doc=="before.RW"){
             DT4 <- cbind(DT4, row.INIT - data.frame(rowSums(DT4))); colnames(DT4)[ncol(DT4)] <- "RemovedWords"
             df.temp <- data.frame(DT4, row.INIT,df.context.quali[rownames(DT4),ncol(df.context.quali),drop=FALSE])      
      } else {
        if(marg.doc=="after") row.INIT <- data.frame("row.INIT"=rowSums(DT4))
         df.temp <- data.frame(DT4,row.INIT, df.context.quali[rownames(DT4),ncol(df.context.quali),drop=FALSE])
      }  # End marg.doc


      df.temp.agg <- aggregate(.~new.Cat, df.temp, sum)
      rownames(df.temp.agg) <- df.temp.agg[,1]
      df.temp.agg <- df.temp.agg[,-1]
      df.temp.agg <- df.temp.agg[rowSums(df.temp.agg[1:(ncol(df.temp.agg)-1)])!=0,]
      row.init.QL <- df.temp.agg[, "row.INIT", drop=FALSE]
      DT3 <- df.temp.agg[,-ncol(df.temp.agg)]
      } else {   # no aggr
        df.temp <- data.frame(DocTerm,row.INIT, df.context.quali[rownames(DocTerm),ncol(df.context.quali),drop=FALSE])
        df.temp.agg <- aggregate(.~new.Cat, df.temp, sum)
        rownames(df.temp.agg) <- df.temp.agg[,1]
        df.temp.agg <- df.temp.agg[,-1]
        df.temp.agg <- df.temp.agg[rowSums(df.temp.agg[1:(ncol(df.temp.agg)-1)])!=0,]
        row.init.QL <- df.temp.agg[, ncol(df.temp.agg), drop=FALSE]
        DT3 <- df.temp.agg[,-ncol(df.temp.agg)]
    }# End no aggregated

    
  } else {    # No bTextData

    
    df.context.quali <- object[,strquali,drop=FALSE]
    
    if(ncol(df.context.quali)>0) {
      df.context.quali$new.Cat <- ""
      df.context.quali$new.Cat  <- paste0(df.context.quali$new.Cat, paste0(colnames(df.context.quali)[1],"."), df.context.quali[,1])
      
     if(ncol(df.context.quali)>2)      
      for(i in 2:(ncol(df.context.quali)-1)) {
        df.context.quali$new.Cat  <- paste0(df.context.quali$new.Cat, 
                                            paste0("@",colnames(df.context.quali)[i],"."), df.context.quali[,i])
      } }
  
    df.temp <- data.frame(DocTerm,row.INIT, df.context.quali[rownames(DocTerm),ncol(df.context.quali),drop=FALSE])
    df.temp.agg <- aggregate(.~new.Cat, df.temp, sum)
    rownames(df.temp.agg) <- df.temp.agg[,1]
    df.temp.agg <- df.temp.agg[,-1]
    df.temp.agg <- df.temp.agg[rowSums(df.temp.agg[1:(ncol(df.temp.agg)-1)])!=0,]
    row.init.QL <- df.temp.agg[, ncol(df.temp.agg), drop=FALSE]
    row.init.QL <-  row.init.QL[,1,drop=TRUE]
    DT3 <- df.temp.agg[,-ncol(df.temp.agg)]
    if(length(strquanti)==0) strquanti <- NULL
    if(!is.null(strquanti)){
      df.context.quanti <- data.frame(object[rownames(df.context.quali),strquanti,drop=FALSE],df.context.quali$new.Cat)
    }
  }
    col.INIT <- colSums(DT3)
    resCharWord.QL <- descfreq_NewChar(DT3, proba = proba, row.init.QL, col.INIT) 
    
    d.single <- chi.quali.single(DT3)
    res$Vocab$quali$stats <- d.single 
    res$Vocab$quali$CharWord <- resCharWord.QL
  } # End (!is.null(strquali))
    


  
  
  #=================================================================	  
  ##### step 5. Quanti contextual variables ##############
  if(!is.null(strquanti)){
    
    if(bTextData) {
        if(is.null(strquali)) {

          if(bvaragg) {
          X.quanti.TEMP <- cbind(object$SourceTerm.quant[rownames(object$SourceTerm.var.agg),strquanti,drop=FALSE],
                                 object$SourceTerm.var.agg)
          colnames(X.quanti.TEMP)[ncol(X.quanti.TEMP)] <- "new.Cat"
          X.quanti.TEMP <- aggregate(.~new.Cat, X.quanti.TEMP, mean)
          rownames(X.quanti.TEMP) <- X.quanti.TEMP$new.Cat
          X.quanti.TEMP <- X.quanti.TEMP[rownames(row.INIT.B),-1,drop=FALSE]
          res$Vocab$quanti <- vocabQuanti(DocTerm[rownames(row.INIT.B),,drop=FALSE],
                                          X.quanti.TEMP,row.INIT.B)
          
         } else {  # No aggregated
           X.quanti.TEMP <- object$SourceTerm.quant[rownames(DocTerm),strquanti,drop=FALSE]
           res$Vocab$quanti <- vocabQuanti(DocTerm[rownames(row.INIT.B),,drop=FALSE],
                                           X.quanti.TEMP,row.INIT.B)
       
         }
      } else {    # There are quali variables
        X.quanti.agg <- aggregate(.~new.Cat, X.quanti.agg, mean)
        rownames(X.quanti.agg) <- X.quanti.agg[,1]
        X.quanti.agg <- X.quanti.agg[rownames(row.init.QL),-1,drop=FALSE]
        res$Vocab$quanti.aggr <- vocabQuanti(DT3[rownames(row.init.QL),,drop=FALSE],
                                             X.quanti.agg[rownames(row.init.QL),,drop=FALSE],row.init.QL) 
      }
    }  # End bTextData==TRUE
    
    if(bTextData==FALSE) {

      row.INIT.B <- data.frame(rowSums(DocTerm)) 

      # Relationships with original documents, no qualitative contextual variables are considered
      res$Vocab$quanti <- vocabQuanti(DocTerm[rownames(row.INIT.B),,drop=FALSE],
                                      df.context.quanti[rownames(row.INIT.B),strquanti,drop=FALSE],row.INIT.B)
  
      stop("0000000000000000000000000")
      
          
      
      if(!is.null(strquali)) {  # quali and quanti contextual
        df.temp.quali <- data.frame(DocTerm,row.INIT, df.context.quali[rownames(DocTerm),ncol(df.context.quali),drop=FALSE])
        df.temp.quali.agg <- aggregate(.~new.Cat, df.temp.quali, sum)
        rownames(df.temp.quali.agg) <- df.temp.quali.agg[,1]
        df.temp.quali.agg <- df.temp.quali.agg[-1]
        row.INIT.B <- df.temp.quali.agg[,ncol(df.temp.quali.agg),drop=FALSE]
        DocTerm.quant <- df.temp.quali.agg[,-ncol(df.temp.quali.agg),drop=FALSE]
        colnames(df.context.quanti)[ncol(df.context.quanti)] <- "new.Cat"
        df.temp.quanti.agg <- aggregate(.~new.Cat, df.context.quanti, mean)
        rownames(df.temp.quanti.agg) <- df.temp.quanti.agg$new.Cat
        df.temp.quanti.agg <- df.temp.quanti.agg[,-1, drop=FALSE]
        res$Vocab$quanti.aggr <- vocabQuanti(DocTerm.quant[rownames(row.INIT.B),,drop=FALSE],
                                             df.temp.quanti.agg[rownames(row.INIT.B),,drop=FALSE],row.INIT.B)
      }
  }}

  

#============================= 

  if(bTextData & !bvaragg & is.null(strquali)) maxDocs <- 0
  if(!bTextData & is.null(strquali)) maxDocs <- 0

  
  if(maxDocs>0)  {
    if(!bTextData) { # Not TextData objext
      df.QUAL<- df.context.quali[,"new.Cat", drop=FALSE]
      motsval <- sapply(resCharWord.QL,function(x) if(!is.null(x))  x[order(rownames(x)),6,drop=FALSE],simplify = FALSE)
      DT3 <- DocTerm
    }
    
    if(bTextData)
    if(bvaragg | !is.null(strquali)) {   # Yes Extraction of the characteristic documents
      #  stop("Tiene sentido extraer ducumentos")     
      var.text <- object$info$var.text[[1]]
      str.base <- object$info$base[[1]]
      str.envir <- object$info$menvir[[1]]
      base <- get(str.base, envir=str.envir)
      corpus <- base[, var.text[1]]
      #--------- Save corpus var.text  -------------------						
      if(length(var.text) > 1) {					
        for (i in 2:length(var.text)){					
          corpus <- paste(corpus, base[, var.text[i]], sep = ".")}}
      
      corpus <- data.frame(corpus, stringsAsFactors = FALSE)
      rownames(corpus) <- rownames(base)
      # corpus <- corpus[rownames(dtm4),]

       if(bvaragg) DT3<- as.matrix(object$SourceTerm) else DT3 <- DocTerm
      corpus <- corpus[rownames(DT3),,drop=FALSE]

       if(is.null(strquali)) { # only aggregated
         df.QUAL<- object$var.agg[rownames(DT3),, drop=FALSE]
         colnames(df.QUAL) <- "new.Cat"
         resCharWord.QL <- resCharWord
         motsval <- sapply(resCharWord,function(x) if(!is.null(x))  x[order(rownames(x)),6,drop=FALSE],simplify = FALSE)
      } else {
        df.QUAL <- df.context.quali[rownames(DT3),"new.Cat",drop=FALSE]
        motsval <- sapply(resCharWord.QL,function(x) if(!is.null(x))  x[order(rownames(x)),6,drop=FALSE],simplify = FALSE)
      }
}

      ########################################
      # motsval vtest of words for each group
      # Intern %     glob % Intern freq Glob freq       p.value    v.test
      # For each group, alphabetical order of words and |vtest|>1.96 
      # motsval <- sapply(resCharWord.QL,function(x) if(!is.null(x))  x[order(row.names(x)),6,drop=FALSE],simplify = TRUE)
      # nlev Nombre of categories
      vlev <- names(resCharWord.QL)
      nlev <- length(vlev)	
      lisresult <- vector(mode="list",length=nlev)

      ######### -----------------------------      
      for (ilev in 1:nlev)     {
        # docpos doc position for each group, starting the firs ilev=1)
        docpos <- which(df.QUAL$new.Cat == vlev[ilev])  # vlev vector with the names of aggregated categories

        # ntrep maximum nombre of documents to print, minumum between maxPrnDoc and the size of the group
        ntrep <- min(maxCharDoc, length(docpos))
        lisresult[[ilev]] <- data.frame(nlev,3)  # Name of the group
        SourceTermcurs <- DT3[docpos,,drop=FALSE ] 
        # SourceTermcurs the selected rows for non aggregated table

        #lisresult[[ilev]] <- SourceTermcurs
        ly<-as.matrix(rep(0,ncol(DT3)))
        # Rows are the vocabulary of after words
        rownames(ly)<-colnames(DT3)

        
        if(is.null(motsval[[ilev]])) {  
          lisresult[ilev] <- list(NULL)
        } else {
          # There are significant words for this category
          motsvalcurs<- as.matrix(motsval[[ilev]])
          motsvalcurs<-	motsvalcurs[which(rownames(motsvalcurs)!="RemovedWords"),,drop=FALSE ]
          
          # ly has the words in rows + vtest =0 if not significant, else the vtest (+ or -)
          ly[rownames(motsvalcurs),1] <- motsvalcurs[,1]
          
      
          
          # SourceTermcurs has the docs of the category in rows and the vocabulary in columns
          # a use vtest + and -
          a <- crossprod(ly,t(SourceTermcurs))
          b <- rowSums(SourceTermcurs)          # sum of vocabulary after of category, categories in columns
          repvaltest <- a
          repvaltest[b > 0] <- a[b > 0]/b[b > 0]

          # Order of the docs into category
          ordrep <- order(repvaltest, decreasing = "TRUE")

          # ntrep is the number of docs to print
          for (i in 1:ntrep)
          {
            lisresult[[ilev]][i,1] <- rownames(DT3)[docpos[ordrep[i]]]
             lisresult[[ilev]][i,2] <- repvaltest[ordrep[i]]
             if(bTextData)
            lisresult[[ilev]][i,3] <- substr(corpus[docpos[ordrep[i]],], start = 1, stop = maxPrnDoc)
             else lisresult[[ilev]][i,3] <- ""
          }
          if(bTextData)
          colnames(lisresult[[ilev]]) <- c("DOCUMENT", "CRITERION", "---------------------TEXT---------------------")	
          else colnames(lisresult[[ilev]]) <- c("DOCUMENT", "CRITERION", "")	
        } # End is.null(motsval[[ilev]]
        
      } # End ilev
      
      names(lisresult) <- vlev
      res$CharDoc <- lisresult
  } # End if(maxDocs>0) 

  res$Proba <- proba
  class(res) <- c("LexChar", "list")  
  return(res)
  
}
