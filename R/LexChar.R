#' @import stringr
#' @export
LexChar <- function(object,  proba = 0.05, maxCharDoc = 10, maxPrnDoc=100, marg.doc="before",
                    context=NULL, correct=TRUE,  nbsample=500, seed=12345, ...)
 { 
  set.seed(seed)
  options(stringsAsFactors = FALSE)
  # maxCharDoc: Maximum number of characteristic documents to show
  # correct = TRUE and Correction for pvalue of Characteristic words (FactoMineR, Agresti...)
  # marg.doc; after/before/before.RW

  # To compute vtest from pvalue when pvalue=0.000
  dots <- list(...)
  if("eps" %in% names(dots)) eps <- dots$eps else eps <- 0.0001

  # Version 1.4.1 used context.sup, version 1.4.2. use context argument
  varnames<- lapply(substitute(list(...))[-1], deparse)
  if(!is.null(varnames$context.sup)) {
    context <- gsub("[[:punct:]]", "", varnames$context.sup)   # Only for Compatibility version 1.4.1
    warning("Xplortext Versions > 1.4.1 use context, no context.sup argument") 
  }
  ###################################################################

  if(marg.doc!="before"& marg.doc!="before.RW" & marg.doc!="after") stop("Error in marg.doc argument")
  
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
    vX <- as.data.frame(vX)  # Quantitative values of variables

    for(i in 1:ncol(vX)){
      if(any(is.na(vX[,i]))) warning("\n", names(vX)[i], " variable has missing values. They will be replaced by the mean\n") 
      vX[is.na(vX[,i]), i] <- mean(vX[,i], na.rm = TRUE)
    }
 
    # Weight documents. Different from marg.doc
    Wi <- vrow.INIT/sum(vrow.INIT)
    col.INIT <- colSums(vdt)
  
    y <- sum(col.INIT)
    #  Weighted average of quantitative variables  (1 x 2 variables)
    
    Aver.X <- apply(vX,2,mean.p,Wi)
    Var.X<-apply(sweep(vX,2,Aver.X,"-"),2,var.p,Wi)

    #Average of quantitative variables for each word
    Wij.cj <- as.matrix(sweep(vdt,2,col.INIT ,"/"))
    
    nk <- colSums(vdt)

    ntot <- sum(vrow.INIT)
    MeanWord <- t(vX) %*% Wij.cj 

    # Variance of words
    coef <- as.matrix((y/col.INIT-1)/(y-1))
    Var.Y <- coef %*% Var.X
    sc.inter <-   apply(sweep(t(sweep(t(MeanWord), 2, Aver.X,"-")^2),2,col.INIT,"*"),1,sum)/y
    pct.expl <- sc.inter/Var.X
    Desv.Y <- sqrt(Var.Y)
    dj <- t(sweep(MeanWord,1,Aver.X,"-")) / Desv.Y
    
    # colSums(vdt)[colnames(vdt)[i]] 
    

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
    
    # res.mat <- matrix(nrow=ncol(vX),ncol=4,0) 
    res.mat <- matrix(nrow=ncol(vX),ncol=5,0) 	
    rownames(res.mat)<-colnames(vX)
   # colnames(res.mat )<- c("GlobalAverage","AverageWord","Differ.","pvalue")
    colnames(res.mat )<- c("GlobalAverage","AverageWord","Differ.","vtest", "pvalue")
    #    names(relcq.perm) <- colnames(vdt)


    
    for (i in 1:nWords) {
      relcq.perm[[i]] <- res.mat
      relcq.perm[[i]][,1] <- t(Aver.X)
      relcq.perm[[i]][,2] <- MeanWord[,i]
      relcq.perm[[i]][,3] <-  relcq.perm[[i]][,2]- relcq.perm[[i]][,1]
      # vtest non weighted
      # relcq.perm[[i]][,4] <- (relcq.perm[[i]][,2]- relcq.perm[[i]][,1] )/ 
      #    sqrt( (ntot- nk[colnames(MeanWord)[i]]) * t(Var.X)/ ((ntot-1)* nk[colnames(MeanWord)[i]]) )
      # 1.- Numerator relcq.perm[[i]][,2]- relcq.perm[[i]][,1]
      # 2.- Variance quantitative variable t(Var.X)
      # 3. n is ntot (total occurrences)
      # 4. nk occurrences word k
      relcq.perm[[i]][,4] <- relcq.perm[[i]][,5] <- t(a[i,])
       relcq.perm[[i]][,4][relcq.perm[[i]][,4]<eps ] <- eps 
     # relcq.perm[[i]][,4] <- sign(relcq.perm[[i]][,3]) * abs(qnorm(relcq.perm[[i]][,4]))
       relcq.perm[[i]][,4] <- sign(relcq.perm[[i]][,3]) * abs(qnorm(relcq.perm[[i]][,4]/2))

      names(relcq.perm)[[i]] <- colnames(vdt)[i]
      
      relcq.perm[[i]] <- subset(relcq.perm[[i]], relcq.perm[[i]][,5] < proba)
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
  # Computing strquanti and strquali
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
    strquali<- context[which(context %in% colnames(object$SourceTerm.qual))]
    # strquanti<- colnames(object$SourceTerm.quant)[which(colnames(object$SourceTerm.quant) %in% context)]
    strquanti<- context[which(context %in% colnames(object$SourceTerm.quant))]
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
    if(length(context)==1)
      if(context=="ALL") stop("You must define context variables by index or name of object")

    if(is.numeric(context)) context <- colnames(object)[context]
   #  sel.context <- colnames(object)[which(context %in% colnames(object))]
    sel.context <- context[which(context %in% colnames(object))]

    strquali <- strquanti <- NULL
    

    if(length(sel.context)!=0) {
      df.context <- object[,sel.context,drop=FALSE]
      context.ordered <-  names(df.context)[sapply(df.context, is.ordered)]
      context.type <- sapply(df.context, class)
      
      strquali <-  names(context.type[which(context.type %in% c("factor","character", "logical", "Date"))])
      strquali <-  c(strquali,context.ordered)
      strquanti <- sel.context[!sel.context %in% strquali]

      if(length(strquanti)==0) strquanti <- NULL
      if(length(strquali)==0) strquali <- NULL
    }
    
  }  # End step 2
  
 

  # STEP 3. 
  #=================================================================	
  # if after/before, frequencies after/before TextData selection are used as document weighting (by default "before");
  #  if before.RW all words under threshold in TextData function are included as a new word named RemovedWords
  ##### step 3. Building object to descfreq_NewChar function ##############
  if(bTextData) { # There is TextData object
    # 1. Computing DocTerm

    DT2 <- as.matrix(object$DocTerm)
    if(marg.doc=="after") {
      DT2 <- DT2[rowSums(DT2)!=0,]
      row.INIT.2 <- data.frame(Occurrences.after=rowSums(DT2))
    }

    if(marg.doc=="before") {
      DT2 <- DT2[rowSums(DT2)!=0,]
      row.INIT.2 <- as.data.frame(object$rowINIT)        # Frequency Occurrences.before
      row.INIT.2 <- row.INIT.2[rownames(DT2),,drop=FALSE]
    }

     if(marg.doc=="before.RW") {
      row.INIT.2 <- as.data.frame(object$summDoc[, "Occurrences.before",drop=FALSE])  # All documents
      rownames(row.INIT.2) <- object$summDoc$DocName
      row.INIT.2 <- row.INIT.2[row.INIT.2$Occurrences.before!=0,,drop=FALSE] 
      # NoNullBefore rows with no null documents before but after threshold
      NoNullBefore <- row.INIT.2[!rownames(row.INIT.2) %in% rownames(DT2),,drop=FALSE]
  
      if(nrow(NoNullBefore)>0){
        DT2 <- cbind(DT2,as.data.frame(row.INIT.2[rownames(DT2),]))
        nrDT2 <- nrow(DT2)
        strnames <- rownames(NoNullBefore)
        # Add NoNullBefore rows to dataframe and fill them with 0
        DT2[( nrDT2+1):( nrDT2+  nrow(NoNullBefore)),] <- 0
        # Add rownames to new rows from strnames
        rownames(DT2)[(nrDT2+1):(nrDT2+nrow(NoNullBefore))]  <-  strnames 
        DT2[(nrDT2+1):( nrDT2+  nrow(NoNullBefore)),ncol(DT2)] <- NoNullBefore[,1]
        DT2 <- DT2[rownames(row.INIT.2),]
      } 
      else {
        DT2 <- cbind(DT2,as.data.frame(row.INIT.2[rownames(DT2),]))
      }# End if(nrow(NoNullBefore)>0)
      
      # Is before.RW"
      freq.after <- apply(DT2[,c(1:(ncol(DT2)-1))],1,sum)
      freq.before <- DT2[,ncol(DT2)]  
      # Adding a column with the name of RemovedWords
      DT2 <- cbind(DT2,as.data.frame(freq.before-freq.after))
      colnames(DT2)[ncol(DT2)] <- "RemovedWords"
      DT2 <- DT2[,-(ncol(DT2)-1)]
      #  colnames(DT2)[ncol(DT2)] <- "RemovedWords"
    }  # End before.RW
    # DocTerm can have some document/aggregate document with margin zero
    
    row.INIT.2 <- row.INIT.2[rownames(DT2),,drop=FALSE]
    pos.0 <- which(rowSums(DT2)==0)
    if(length(pos.0)!=0) {
      # Remove empty documents or aggregate documents 
      DT2 <- DT2[!rownames(DT2) %in% names(pos.0),,drop=FALSE]
      row.INIT.2 <- row.INIT.2[!rownames(row.INIT.2) %in% names(pos.0),,drop=FALSE]
    }  
    # dbTextData 
  } # End bTextData for step3
  else{   # Is not TextData
    # Is not bTextData
    # We must have DocTerm without quali and quanti variables
    
     strQLQN <- c(strquali, strquanti)
    DT2.init <- data.frame(object)
    DT2 <- object[,!colnames(DT2.init) %in% strQLQN,drop=FALSE]  # without context variables
    
    df.context.quali <- df.context.quanti <- NULL
    #  pos.0 <- which(rowSums(DT2)==0)
    # Remove empty documents of non TextData object 
    #   if(length(pos.0)!=0)   DT2 <- DT2[!rownames(DT2) %in% rownames(DT2)[pos.0],,drop=FALSE]

    row.INIT.2 <- data.frame(rowSums(DT2))

    if(!is.null(strquali)) df.context.quali<- DT2.init[rownames(DT2),strquali,drop=FALSE]
    if(!is.null(strquanti)) df.context.quanti<- DT2.init[rownames(DT2),strquanti,drop=FALSE]
    #    stop("Is not bTextData")
  }
  #--------------------------------------
  # This is valid for any kind of objects
  col.INIT <- as.vector(colSums(DT2))
  # row.INIT.B <- row.INIT
  #  Computing Lexical Table: Chi.Squared  df   p.value       phi Cramer's V. No weigthing is applied
  n.words<- ncol(DT2)
  d.single <- chi.quali.single(DT2[rowSums(DT2)>0,c(1:n.words)]) 
  #--------------------------------------------------------------------------------------------
  resCharWord <- descfreq_NewChar(DT2, proba = proba, row.INIT.2, col.INIT) 	
  # Return words for active categories
  res <- list(CharWord=resCharWord, stats=d.single)
  res$Proba <- proba
  # Chi square not depend on row.INIT
  


  ###########################################################################################################
  if(length(strquali)==0) strquali <- NULL
  if(length(strquanti)==0) strquanti <- NULL
  strQNQL <- c(strquali, strquanti)
  DT.Z <- NULL
  

  
  ##########################  Contextual variables. DT3
  new.QUALI <- strquali 
  if(bTextData & bvaragg)  new.QUALI <- c(colnames(object$var.agg), new.QUALI)
  # new.QUALI <- unique(c(strQNQL, new.QUALI))  # Quanti, Quali and var.agg
  


  ###################  NO QUALI BUT QUANTI CONTEXTUAL VARIABLES ###################
  
  ##### step XX. Quali / quanti contextual variables ##############
  if(is.null(strquali) & !is.null(strquanti)) {
    # Only quanti
    if(bTextData) {
      if(bvaragg) { # TextData, bvaragg, no quali, yes quanti
        df.context.quanti <- as.matrix(object$context$quanti)
        ROW.Q <- row.INIT.2
      } else { # "TextData, no bvaragg, no quali, yes quanti"
        ROW.Q <- row.INIT.2[rownames(DT2),,drop=FALSE]
        df.context.quanti <- object$SourceTerm.quant
        df.context.quanti <-  df.context.quanti[rownames(DT2),,drop=FALSE]
      }
      #  ROW.W <- ROW.W[rownames(DT2),,drop=FALSE ]
      res$quanti <- vocabQuanti(DT2, df.context.quanti, ROW.Q)
    } else { # 
      # No TextData
      res$quanti <- vocabQuanti(DT2, df.context.quanti, row.INIT.2)
    }
  }
  
  if(length(new.QUALI)==0)  new.QUALI <- NULL
  
  

  
  ##################  If there are context qualitative variables
  if(!is.null(new.QUALI)){  # new.QUALI has the names of var.agg and strquali
    DT.Z2 <- NULL
    if(bvaragg){  # Es agregado
      DT.Z <- data.frame(as.matrix(object$SourceTerm))       # 292 x 24
      df.AGGR <- data.frame(object$var.agg)                 # 292
      df.AGGR <- df.AGGR[rownames(DT.Z),,drop=FALSE]
      
      if(!is.null(strquali)) {
        df.context.quali <- object$SourceTerm.qual[,strquali,drop=FALSE]
        df.context.quali <- df.context.quali[rownames(DT.Z),,drop=FALSE]
        df.AGGR <- cbind(df.AGGR,df.context.quali)               # 292
      }
      
      df.context.quali <- df.AGGR
      rm(df.AGGR)
      strquali <- colnames(df.context.quali)
      
      if(marg.doc=="after") {
        ROW.W <- data.frame(rowSums(DT.Z))
        rownames(ROW.W) <- rownames(DT.Z)  # 292
      }
      
      if(marg.doc=="before") {
        ROW.W <- object$SourceTerm.freq # 292
        rownames(ROW.W) <- ROW.W$DocName
        ROW.W <- ROW.W[,-1,drop=FALSE]
      }
      
      if(marg.doc=="before.RW") {
        DT.Temp <- as.matrix(object$SourceTerm.dtm)
        DT.Temp <- DT.Temp[rowSums(DT.Temp)!=0,]   # 298
        
        ROW.W <- data.frame(rowSums(DT.Temp))
        rownames(ROW.W) <- rownames(DT.Temp)       # 298
        rm(DT.Temp)
        NoNullBefore <- ROW.W[!rownames(ROW.W) %in% rownames(DT.Z),,drop=FALSE]
        
        if(nrow(NoNullBefore)>0){
          DT.Z2 <- cbind(DT.Z,as.data.frame(ROW.W[rownames(DT.Z),]))  # 292
          nrDT2 <- nrow(DT.Z2)
          strnames <- rownames(NoNullBefore)
          
          # Add NoNullBefore rows to dataframe and fill them with 0
          DT.Z2[( nrDT2+1):( nrDT2+  nrow(NoNullBefore)),] <- 0     # 298
          rownames(DT.Z2)[(nrDT2+1):(nrDT2+nrow(NoNullBefore))]  <-  strnames 
          DT.Z2[(nrDT2+1):( nrDT2+  nrow(NoNullBefore)),ncol(DT.Z2)] <- NoNullBefore[,1]
          DT.Z <- DT.Z2[rownames(ROW.W),]
          # Add rownames to new rows from strnames
        } 
        else {
          DT.Z2 <- cbind(DT.Z,as.data.frame(ROW.W[rownames(DT.Z),]))
        }# End if(nrow(NoNullBefore)>0)
        
        # Is before.RW"
        freq.after <- apply(DT.Z2[,c(1:(ncol(DT.Z2)-1))],1,sum)
        freq.before <- DT.Z2[,ncol(DT.Z2)]  
        # Adding a column with the name of RemovedWords
        DT.Z2 <- cbind(DT.Z2,as.data.frame(freq.before-freq.after))
        colnames(DT.Z2)[ncol(DT.Z2)] <- "RemovedWords"
        DT.Z <- DT.Z2[,-(ncol(DT.Z2)-1)]
        #  colnames(DT2)[ncol(DT2)] <- "RemovedWords"
        
        df.context.quali <- object$SourceTerm.var.agg[rownames(DT.Z),,drop=FALSE]
        if(!is.null(object$SourceTerm.qual))
          df.context.quali <- cbind(df.context.quali,object$SourceTerm.qual[rownames(DT.Z),,drop=FALSE] )
      }
      if(!is.null(strquanti)) {
        df.context.quanti <- data.frame(object$SourceTerm.quant)  # <-------------------------------300.
        df.context.quanti <- df.context.quanti[rownames(DT.Z),,drop=FALSE]  # 292 o 298
      }
    } # End bvaragg aggregate
    
    
  
  
    
    
    if(!bvaragg & bTextData){  # Is not aggregate
      #  df.context.quali <- object$context$quali   # 292
      df.context.quali <- object$SourceTerm.qual[,strquali,drop=FALSE]  # 300
      DT.Z <- data.frame(as.matrix(object$DocTerm))      # 292 (always) x 24
      
      
      if(marg.doc=="after") {
        ROW.W <- data.frame(rowSums(DT.Z))
        rownames(ROW.W) <- rownames(DT.Z)  # 292  
        DT.Z2 <- DT.Z[rownames(ROW.W),,drop=FALSE]
        if(!is.null(strquanti)) {
          df.context.quanti <- as.matrix(object$context$quanti)  # 292
          df.context.quanti <- df.context.quanti[rownames(DT.Z),,drop=FALSE]
        }
      }
      
      if(marg.doc=="before"){
        ROW.W <- object$rowINIT                 # 300  
        ROW.W <- ROW.W[rownames(DT.Z),,drop=FALSE]  # 292
        DT.Z2 <- DT.Z[rownames(ROW.W),,drop=FALSE]
        if(!is.null(strquanti)) {
          df.context.quanti <- as.matrix(object$context$quanti)  # 292
          df.context.quanti <- df.context.quanti[rownames(DT.Z),,drop=FALSE]
        }
      } 
      
      if(marg.doc=="before.RW") {  
        ROW.W <- object$rowINIT                 # 300  
        ROW.W <- ROW.W[ROW.W$Occurrences.before!=0,,drop=FALSE]  # 298
        
        # NoNullBefore rows with no null documents before but after threshold
        NoNullBefore <- ROW.W[!rownames(ROW.W) %in% rownames(DT.Z),,drop=FALSE]
        
        if(nrow(NoNullBefore)>0){
          DT.Z2 <- cbind(DT.Z,as.data.frame(ROW.W[rownames(DT.Z),]))  # 292
          nrDT2 <- nrow(DT2)
          strnames <- rownames(NoNullBefore)
          # Add NoNullBefore rows to dataframe and fill them with 0
          DT.Z2[( nrDT2+1):( nrDT2+  nrow(NoNullBefore)),] <- 0     # 304
          # Add rownames to new rows from strnames
          rownames(DT.Z2)[(nrDT2+1):(nrDT2+nrow(NoNullBefore))]  <-  strnames 
          DT.Z2[(nrDT2+1):( nrDT2+  nrow(NoNullBefore)),ncol(DT.Z2)] <- NoNullBefore[,1]
          DT.Z2 <- DT.Z2[rownames(ROW.W),]
          DT.Z <- DT.Z2
        } 
        else {
          DT.Z2 <- cbind(DT.Z2,as.data.frame(ROW.W[rownames(DT.Z2),]))
        }# End if(nrow(NoNullBefore)>0)
        # Is before.RW"
        freq.after <- apply(DT.Z2[,c(1:(ncol(DT.Z2)-1))],1,sum)
        freq.before <- DT.Z2[,ncol(DT.Z2)]  
        # Adding a column with the name of RemovedWords
        DT.Z2 <- cbind(DT.Z2,as.data.frame(freq.before-freq.after))
        
        colnames(DT.Z2)[ncol(DT.Z2)] <- "RemovedWords"
        DT.Z2 <- DT.Z2[,-(ncol(DT.Z2)-1)]
        DT.Z <- DT.Z2
      }
      
      df.context.quali <- df.context.quali[rownames(DT.Z),,drop=FALSE]
      
      if(!is.null(strquanti)) {
        df.context.quanti <- data.frame(object$SourceTerm.quant)
        df.context.quanti <- df.context.quanti[rownames(DT.Z),,drop=FALSE]
      }
    } # End !bvaragg
    
    if(!is.null(DT.Z2)) rm(DT.Z2)
    
   } ####### End quali
  
  

  
  ####### 
  if(!is.null(new.QUALI)) {
 
    if(!bTextData)  {
      DT.Z <- DT2
      ROW.W <- data.frame(rowSums(DT.Z))
      rownames(ROW.W) <- rownames(DT.Z)
    }

    
    # df.context.quali si es mayor que uno unir con @
    # if df.context.quanti
    # ROW.W
    # return(DT.Z)
    df.context.quali$new.Cat <- ""
    df.context.quali$new.Cat  <- paste0(df.context.quali$new.Cat, paste0(colnames(df.context.quali)[1],"."), df.context.quali[,1])
    if(ncol(df.context.quali)>2)
      for(i in 2:(ncol(df.context.quali)-1)) {
        df.context.quali$new.Cat  <- paste0(df.context.quali$new.Cat, 
                                            paste0("@",colnames(df.context.quali)[i],"."), df.context.quali[,i])
        df.context.quali <- df.context.quali[rownames(DT.Z),,drop=FALSE]          
      }
   

    if(length(strquanti)==0) strquanti <- NULL
    ROW.W <- ROW.W[rownames(DT.Z),,drop=FALSE]
    df.context.quali <- df.context.quali[rownames(DT.Z),,drop=FALSE]
    ROW.W$new.Cat <- df.context.quali$new.Cat
    DT.Z$new.Cat <-  df.context.quali$new.Cat

    # -------  Hacer agregada
    DT.Z.aggr <- aggregate(.~new.Cat, DT.Z, sum) 
    rownames(DT.Z.aggr) <- DT.Z.aggr$new.Cat
    DT.Z.aggr <- DT.Z.aggr[,-1,drop=FALSE]
  

    ROW.W.aggr <- aggregate(.~new.Cat, ROW.W, sum) 
    rownames(ROW.W.aggr) <- ROW.W.aggr$new.Cat
    ROW.W.aggr <- ROW.W.aggr[,-1,drop=FALSE]
    

    if(!is.null(strquanti)) {
      df.context.quanti <- df.context.quanti[rownames(DT.Z),,drop=FALSE]
      df.context.quanti$new.Cat <- df.context.quali$new.Cat
      df.context.quanti.aggr <- aggregate(.~new.Cat, df.context.quanti, mean) 
      rownames(df.context.quanti.aggr) <- df.context.quanti.aggr$new.Cat
      df.context.quanti.aggr <- df.context.quanti.aggr[,-1,drop=FALSE]
      res$quanti <- vocabQuanti(DT.Z.aggr,
                                df.context.quanti.aggr[rownames(DT.Z.aggr),strquanti,drop=FALSE],
                                ROW.W.aggr[rownames(DT.Z.aggr),,drop=FALSE])
    }
    res$quali$CharWord  <- descfreq_NewChar(DT.Z.aggr, proba = proba, ROW.W.aggr[rownames(DT.Z.aggr),,drop=FALSE],
                                            col.INIT) 
    res$quali$stats <- chi.quali.single(DT.Z.aggr)
  }
  
  
 
  
  ############################################3
  fCharDoc  <- function(CharWord.z, df.QUAL.z, DT.z) {
    
    motsval <- sapply(CharWord.z ,function(x) if(!is.null(x))  x[order(rownames(x)),6,drop=FALSE],simplify = FALSE)
    ########################################
    # motsval vtest of words for each group
    # Intern %     glob % Intern freq Glob freq       p.value    v.test
    # For each group, alphabetical order of words and |vtest|>1.96 
    # motsval <- sapply(resCharWord.QL,function(x) if(!is.null(x))  x[order(row.names(x)),6,drop=FALSE],simplify = TRUE)
    # vlev Name of supplementary categories 
    vlev <- names(CharWord.z)
    nlev <- length(vlev)	
    lisresult <- vector(mode="list",length=nlev)
    DT.z <- DT.z[,-ncol(DT.z),drop=FALSE]
    
    
    ######### -----------------------------  ´
    for (ilev in 1:nlev)     {
      # docpos doc position for each group, starting the firs ilev=1)
      docpos <- which(df.QUAL.z$new.Cat == vlev[ilev])  # vlev vector with the docs names of aggregated categories

      
      # ntrep maximum number of documents to print, minimum between maxPrnDoc and the size of the group
      ntrep <- min(maxCharDoc, length(docpos))   # length(docpos) is the number of documents in each group
      
      lisresult[[ilev]] <- data.frame(nlev,3)  # Name of the group

      SourceTermcurs <- DT.z[docpos,,drop=FALSE ] 
      # SourceTermcurs the selected rows for non aggregated table
      #lisresult[[ilev]] <- SourceTermcurs
      ly<-as.matrix(rep(0,ncol(DT.z)),drop=FALSE)
      # Rows are the vocabulary of after words
      rownames(ly)<-colnames(DT.z)

      
      if(is.null(motsval[[ilev]])) {  
        lisresult[ilev] <- list(NULL)
      } else {
        # There are significant words for this category
        motsvalcurs<- as.matrix(motsval[[ilev]])
        motsvalcurs<-	motsvalcurs[which(rownames(motsvalcurs)!="RemovedWords"),,drop=FALSE ]
        # ly has the words in rows + vtest =0 if not significant, else the vtest (+ or -)
        #ly[rownames(motsvalcurs),1] <- motsvalcurs[,1]
        # str.motsvalcurs <- motsvalcurs %in% 
        ly[rownames(motsvalcurs),1] <- motsvalcurs[,1]
        
        if(sum(motsvalcurs[,1])==0) {
          lisresult[ilev] <- list(NULL)
        } else
        {          
          
          # SourceTermcurs has the docs of the category in rows and the vocabulary in columns
          # a use vtest + and 
          a <- crossprod(ly,t(SourceTermcurs))
          b <- rowSums(SourceTermcurs)          # sum of vocabulary after of category, categories in columns
          repvaltest <- a
          repvaltest[b > 0] <- a[b > 0]/b[b > 0]
          # New the following line
          # The next is new. To remove CRITERIO zero
          repvaltest <- repvaltest[,colSums(repvaltest)>0,drop=FALSE]
          # Order of the docs into category
          ordrep <- order(repvaltest, decreasing = "TRUE")

         ntrep2 <- min(ntrep, length(ordrep))

          ## return(length(ordrep))  30
          
          # ntrep is the number of docs to print
          for (i in 1:ntrep2)
          {
            lisresult[[ilev]][i,1] <- rownames(DT.z)[docpos[ordrep[i]]]
            lisresult[[ilev]][i,2] <- repvaltest[ordrep[i]]
            
            if(bTextData) {
              lisresult[[ilev]][i,3] <- substr(corpus[docpos[ordrep[i]],], start = 1, stop = maxPrnDoc)
              #  if(length(ordrep)==0) lisresult[[ilev]][i,3] <- NULL
            } else {lisresult[[ilev]][i,3] <- ""}
          }


          if(bTextData) {
            colnames(lisresult[[ilev]]) <- c("DOCUMENT", "CRITERION", "---------------------TEXT---------------------")	
          } else { colnames(lisresult[[ilev]]) <- c("DOCUMENT", "CRITERION", "")	}
        } # End is.null(motsval[[ilev]]
      } # New is. sum(ly[,])==0)
      if(length(ordrep)==0) lisresult[ilev] <- list(NULL)

    } # End for
    names(lisresult) <- vlev
    return(lisresult)
  }  

  
  if(bTextData & !bvaragg & is.null(strquali)) maxCharDoc <- 0 # maxDocs <- 0
  if(!bTextData & is.null(strquali)) maxCharDoc <- 0# maxDocs <- 0
  



 #  if(maxDocs>0)  {
  if(maxCharDoc>0)  {

    if(bTextData) {
      base.new <- object$info$base[[1]]
      var.text <- object$info$var.text[[1]]   #  Important  Relaunch
     # str.base <- object$info$base[[1]]    # open.question, character
    #  str.envir <- object$info$menvir[[1]]

        base.new[, var.text[1]] <- as.character(base.new[, var.text[1]])
        base.new[is.na(base.new[var.text[1]]), var.text[1]] <- ""
        corpus <-  base.new[, var.text[1]]

      if(length(var.text) > 1) {					
        for (i in 2:length(var.text)){	
          corpus2 <- as.character(base.new[, var.text[i]])
          dos <-which(corpus!="" & corpus2!="")
          corpus[dos] <- paste(corpus[dos], corpus2[dos], sep=". ")
          uno <-which(corpus=="" & corpus2!="")
          corpus[uno] <- corpus2[uno]
        }
        rm(corpus2)
      }					
      corpus <- data.frame(corpus, stringsAsFactors = FALSE)	
      rownames(corpus) <- rownames(base.new)					# 300
      
      corpus <- corpus[rownames(DT.Z),,drop=FALSE]   # 292 o 298 for Fmin
      DT6 <- DT.Z[rowSums(DT.Z[,1:(ncol(DT.Z)-1)])!=0,,drop=FALSE]
      corpus <- corpus[rownames(DT6),,drop=FALSE]
      df.QUAL<- df.context.quali[rownames(DT6),"new.Cat",drop=FALSE]
      res$quali$CharDoc <- fCharDoc(res$quali$CharWord, df.QUAL, DT6)
    }
    
     if(!bTextData & !is.null(strquali)){
       df.QUAL<- df.context.quali[rownames(DT.Z),"new.Cat",drop=FALSE]
       res$quali$CharDoc <- fCharDoc(res$quali$CharWord, df.QUAL, DT.Z)
     } 
    
    
  } # End maxDocs maxCharDoc
  

  
  res$Proba <- proba
  class(res) <- c("LexChar", "list")    
  return(res)
}


