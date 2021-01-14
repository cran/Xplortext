#' @import stringr
#' @export
LexChar <- function(object, proba = 0.05, maxDocs=20, maxCharDoc = 10, maxPrnDoc=100, marg.doc="before", correct=TRUE,
                  context.sup="ALL", nbsample=500, seed=12345)
 {
  set.seed(seed)
  options(stringsAsFactors = FALSE)
# maxDocs: Maximum number of input documents or aggregate categories
# maxCharDoc: Maximum number of characteristic documents to show
# correct = TRUE and Correction for pvalue of Characteristic words (FactoMineR, Agresti...)

  if (!inherits(object,"TextData") & !inherits(object,"DocumentTermMatrix") & !inherits(object,"matrix")
      & !inherits(object,"data.frame"))
    stop("Object should be TextData, DocumentTermMatrix, matrix or data frame")
  if(proba<0|proba>1) stop("proba should be between 0 and 1")
  
  if(is.null(maxCharDoc) |  maxCharDoc <0) maxCharDoc <-0 
 
  ############### descfreq_New function #############################################
  descfreq_NewChar <- function (data, proba = 0.05, marge.li, marge.col) 				
  {
    lab.sauv <- lab <- colnames(data)				
    for (i in 1:length(lab)) lab[i] = gsub(" ", ".", lab[i])				
    colnames(data) = lab				
  
    old.warn = options("warn")				
    options(warn = -1)			# suppress warnings globally
    
    nom = tri = structure(vector(mode = "list", length = nrow(data)), names = rownames(data))				
    marge.li <- as.vector(marge.li)
    marge.col <- as.vector(marge.col)
    SumTot<-sum(marge.li)				
    nr <- nrow(data)  				
    nc <- ncol(data) 				
    
    
    for (j in 1:nr) {				
      aux3 <- marge.li[j]/SumTot				# % Ocurrences before or after
      
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
  ######################################################
  ##################  main function  ##################
  ###      Step 1. Verifying the correctness of the arguments
  bTextData <- ifelse(inherits(object,"TextData"),TRUE,FALSE)
  ###### Contextual data
    
  
  
##############  Function to compute Vocab$quali$stats  
  chi.quali<- function(X, QL)  {

    list.Agg <- lapply(seq_along(QL),FUN=function(i) t(t(as.matrix(X))%*% as.matrix(QL[[i]])))

    
    Nq <- sum(X)
    dfq<-NULL
    old.warn = options("warn")				
    options(warn = -1)			# suppress warnings globally
    
    for(i in 1:length(list.Agg)) {
      ch <- chisq.test(list.Agg[[i]], correct=FALSE)[1:3]
      ph <- sqrt(ch$statistic/Nq)
      d <- data.frame(ch,ph)
      dfq <- rbind(dfq,d)
    }
    colnames(dfq) <- c("Chi.Squared", "df", "p.value", "phi")
    rownames(dfq) <-  names(QL)
    options(old.warn)	
    return(dfq)
  }
##########################################
  chi.quali.single <- function(X)  {
    old.warn = options("warn")				
    options(warn = -1)			# suppress warnings globally
    ch <- chisq.test(X, correct=FALSE)[1:3]
    ph <- sqrt(ch$statistic/sum(X))
    d <- data.frame(ch,ph)
    colnames(d) <- c("Chi.Squared", "df", "p.value", "phi")
    rownames(d) <-  "LexicalTable"
    options(old.warn)	
    return(d)
  }
 
  
  vocabQuanti <- function(vdt,vX, vrow.INIT ) {
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

    Aver.X <- apply(X,2,mean.p,Wi)
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
    nWords <- ncol(vdt)                                                           # 1090 words
    nDocs <- nrow(vdt)                                                            # 11 discourses
    
    relcq.perm <-vector(mode='list',length=nWords)
    for (i in 1:nWords) relcq.perm[[i]]<-matrix(nrow=nbsample,ncol=ncol(vX),0) 	  # 500 x (1090*2varcuanti)
        relcq.SS.perm <-vector(mode='list',length=ncol(vX))
    for (i in 1:ncol(vX)) relcq.SS.perm[[i]]<- matrix(nrow=nbsample,ncol=1,0) 	 # vX matriz de variables cuantitat, 500x2
        
    for (i in 1:nbsample){                                                       # nbsample, 500 by default
      X.perm <- X[sample(1:nDocs),,drop=FALSE]
      Aver.X.perm<-apply(X.perm,2,mean.p,Wi)
      MeanWord.perm <- t(X.perm) %*% Wij.cj
      Var.X.perm<-apply(sweep(X.perm,2,Aver.X.perm,"-"),2,var.p,Wi)
      Var.Y.perm <- coef %*% Var.X.perm
      # dj words x variables
      dj.perm <- t(sweep(MeanWord.perm,1,Aver.X.perm,"-")) / sqrt(Var.Y.perm)    # 1090 words x 2 average variables
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
    colnames(res.mat )<- c("GlobalAverage","AverageWord","Difer.","pvalue")
  
  
    
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
 #   res$Vocab$quanti$CharWord <- relcq.perm
 #   res$Vocab$quanti$stats <- res.mat.SS
  }
  
  
  
  
  
  
  ###  To know if it is an aggregate analysis
  # If it is an aggregated table
  bvaragg <- FALSE
  if(inherits(object,"TextData"))
    if(!is.null(object$info)) 
      if(object$info$name.var.agg[[1]]!="") bvaragg <- TRUE 
  
  
  #=================================================================		
  ###### Step 2. Extraction of the characteristic words #### 				     
  
  if(bTextData) { # There is TextData
    DocTerm <- as.matrix(object$DocTerm)
    # Pendiente de las frecuencias row.INIT dependiendo del marg
    if(marg.doc=="before") {
      row.INIT <- object$rowINIT 
      DocTerm <- cbind(DocTerm,as.data.frame(row.INIT -apply(DocTerm,1,sum)))
      colnames(DocTerm)[ncol(DocTerm)] <- "RestOfWords"
    }
    else row.INIT <- data.frame(Ocurrences.After=rowSums(DocTerm))
    row.INIT <- row.INIT[,1]
  } else{  # "CharWord extern docterm"
    DocTerm <- as.matrix(object)
    row.INIT <- rowSums(DocTerm)
    col.INIT <- as.vector(colSums(DocTerm))
  }
  col.INIT <- as.vector(colSums(DocTerm))
  resCharWord <- descfreq_NewChar(DocTerm, proba = 0.05, row.INIT, col.INIT) 		
  d.single <- chi.quali.single(DocTerm) 
  res <- list(CharWord=resCharWord, stats=d.single)
  
  
  
  ##################3###      Step 3. Verifying context variables and extracting information
  # --------------- Check context quanti-quali
  context <- context.sup
  context.quanti <-  context.quali <- NULL 
  #_____________________________________________
  # Fill NA quantitative values with the average
  mean.p <- function(V,poids) res<-sum(V*poids,na.rm=TRUE)/sum(poids[!is.na(V)])
  var.p <- function(V,poids) res<-sum(V^2*poids,na.rm=TRUE)/sum(poids[!is.na(V)])
  
  
  if(!is.null(context)) # Check that names of context.quanti exist		
  if(bTextData) { # There is TextData
    if(bvaragg) { # Is agreggated table
      if(length(context)==1)  if(context=="ALL") {
        context <- as.character(object$context$quali$qualivar[,1])
        context <- c(context,colnames(object$context$quanti))
      }# eND  if(context=="ALL")

      if(!is.null(object$context$quali$qualivar))
        strquali<- object$context$quali$qualivar[,1][which(object$context$quali$qualivar[,1] %in% context)]
        else strquali<- NULL
      if(length(strquali)>0)  {
          numb.qual <- length(strquali)
          dfq <- NULL
          if(marg.doc=="after") {
            df <- object$context$quali$qualivar
            df <- cbind(df, ac = cumsum(df$qualincat))
            df$qualincat <- df$ac - df$qualincat+1
            for(i in 1:numb.qual) {
              df2 <- df[df$qualivar==strquali[i],]
              rowdf2 <- c(df2$qualincat:df2$ac)
              Z2 <- object$context$quali$qualitable[rowdf2,]
              d.single <- chi.quali.single(Z2) 
              rownames(d.single) <- strquali[i]
              dfq <- rbind(dfq,d.single)     # stats for qualivariables after
              resCharWordQL <- descfreq_NewChar(Z2, proba = 0.05, rowSums(Z2), col.INIT) 
              res$Vocab$quali$CharWord[[i]] <-resCharWordQL
            } # end for
          } else {  # before aggregation
            SourceTerm <-as.matrix(object$SourceTerm)
            occ.after <- rowSums(SourceTerm)
            occ.before <- object$SourceTerm.freq$Occurrences.before
            df <- cbind(SourceTerm, "RestOfWords"= (occ.before-occ.after))
            for(i in 1:numb.qual) {
              Zqs <- tab.disjonctif(object$SourceTerm.qual[ rownames(SourceTerm), i])  # 300 x 2
              Z2 <- t(Zqs) %*% as.matrix(df )
              d.single <- chi.quali.single(Z2) 
              rownames(d.single) <- strquali[i]
              dfq <- rbind(dfq,d.single)     # stats for qualivariables before
              resCharWordQL <- descfreq_NewChar(Z2, proba = 0.05, rowSums(Z2), col.INIT) 
              res$Vocab$quali$CharWord[[i]] <-resCharWordQL
            } # End for
                    } # End if(marg.doc=="after")
          res$Vocab$quali$stats <- dfq
          names(res$Vocab$quali$CharWord) <- strquali
                  } # End context.quali
 

        if(!is.null(colnames(object$context$quanti))) 
          strquanti<- colnames(object$context$quanti)[which(colnames(object$context$quanti) %in% context)]
        else strquanti <- NULL
        
        
        if(length(strquanti)>0)  {
            X <- object$context$quanti[,strquanti, drop=FALSE]
            res$Vocab$quanti <- vocabQuanti(DocTerm,X, row.INIT)
          } # End strquanti

 
        
    } else { # Is not agreggated table
      
      if(length(context)==1)  if(context=="ALL") {
        context <- colnames(object$context$quali)
        context <- c(context,colnames(object$context$quanti)) }

      if(!is.null(object$context$quali))
       strquali <-  colnames(object$context$quali)[which(colnames(object$context$quali) %in% context)]
   
         if(!is.null(colnames(object$context$quanti))) 
          strquanti<- colnames(object$context$quanti)[which(colnames(object$context$quanti) %in% context)]
        else strquanti <- NULL
      
      
      
         if(length(strquali)==0) str.quali<-NULL  else {
          #df.context.quali <- context[,context.quali,drop=FALSE]
          numb.qual <- length(strquali)
          dfq <- NULL
          df.context.quali <- object$context$quali[,strquali,drop=FALSE]
          # Hago las tablas de contingencia agregada en una lista
          dis.X <- lapply(colnames(df.context.quali), function(x) FactoMineR::tab.disjonctif(df.context.quali[,x,drop=FALSE])) 
          names(dis.X) <- context.quali
          res.chi.quali <- chi.quali(DocTerm, dis.X)
          rownames(res.chi.quali ) <- strquali
          res$Vocab$quali$stats <- res.chi.quali   
          
          for(i in 1:length(strquali)) {
            AggTable <- t(t(as.matrix(DocTerm))%*% dis.X[[i]])
            resCharWordQL <- descfreq_NewChar(AggTable, proba = 0.05, row.INIT, col.INIT) 
            res$Vocab$quali$CharWord[[i]] <-resCharWordQL
          }
          names(res$Vocab$quali$CharWord) <- strquali
         } # End if(length(strquali)==0)
  
         if(length(strquanti)>0)  {
           X <- object$context$quanti[,strquanti, drop=FALSE]
           res$Vocab$quanti <- vocabQuanti(DocTerm,X, row.INIT)
         } # End strquanti
         
         
         
    }
  } else { # There is not TextData
    
    if(!is.data.frame((context.sup)))
    if(length(context.sup)==1) if(context.sup=="ALL") { context <- context.sup <- NULL }

    if(!is.null(context.sup)) 
      if(is.data.frame(context.sup))  context <- colnames(context.sup) else
          context.sup <- as.data.frame(context.sup, drop=FALSE)

    dt <- as.matrix(object) # DocTerm Matrix
    if(!is.null(context.sup)) {
      
      context.quanti <- colnames(context.sup)[unlist(lapply(context.sup, is.numeric))]
      context.quali  <- colnames(context.sup)[!unlist(lapply(context.sup, is.numeric))]


      
      if(length(context.quali)==0) context.quali<-NULL  else {
        
        df.context.quali <- context.sup[,context.quali,drop=FALSE]
        # To make contingency aggregated tables into a list
       dis.X <- lapply(colnames(df.context.quali), function(x) FactoMineR::tab.disjonctif(df.context.quali[,x,drop=FALSE])) 
        names(dis.X) <- context.quali
       res.chi.quali <- chi.quali(dt, dis.X)
       res$Vocab$quali$stats <- res.chi.quali 

       
       #  list.Agg <- lapply(seq_along(QL),FUN=function(i) t(t(as.matrix(X))%*% QL[[i]]))
       for(i in 1:length(context.quali)) {
         AggTable <- t(t(as.matrix(dt))%*% as.matrix(dis.X[[i]]))
         resCharWordQL <- descfreq_NewChar(AggTable, proba = 0.05, row.INIT, col.INIT) 
         res$Vocab$quali$CharWord[[i]] <-resCharWordQL
       }
#     if(length(res$Vocab$quali$CharWord)> 1)
         names( res$Vocab$quali$CharWord) <- context.quali
        } # End context.quali

      
      
      if(length(context.quanti)==0) context.quanti<-NULL  else {
        
        
        X <- context.sup[,context.quanti,drop=FALSE]
      #  DocTerm <- as.matrix(object)
      #  row.INIT <- rowSums(DocTerm)
        row.INIT <- rowSums(dt)

        
        
        ### Inicio de pasar a función   
        for(i in 1:ncol(X)){
          if(any(is.na(X[,i]))) warning("\n", names(X)[i], " variable has missing values. They will be replaced by the mean\n") 
               X[is.na(X[,i]), i] <- mean(X[,i], na.rm = TRUE)
        }
       res$Vocab$quanti <- vocabQuanti(dt,X, row.INIT)
      } # End strquanti
       }  # End is.null(context)
  } # End !bTextData

  
  if(!bvaragg) maxDocs <- 0 # Not allowed if there is not source documents
    
  ###### Step  3. Extraction of the characteristic documents in the case of an aggregate table	
  if(maxDocs>0)   	
  {	
    vlev <- rownames(DocTerm)			
    nlev <- nrow(DocTerm)
    
    SourceTerm<-as.matrix(object$SourceTerm)
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
    corpus <- corpus[rownames(object$var.agg),]
    # motsval vtest of words for each group
    # Intern %     glob % Intern freq Glob freq       p.value    v.test
    # For each group, alphabetical order of words and |vtest|>1.96 
    motsval <- sapply(resCharWord,function(x) if(!is.null(x))  x[order(row.names(x)),6,drop=FALSE],simplify = TRUE)
    # nlev Nombre of categories
    lisresult <- vector(mode="list",length=nlev)
    if(marg.doc=="before")  DT <- DocTerm[,-ncol(DocTerm)] else DT <- DocTerm
    
    for (ilev in 1:nlev)     {
      # repsel doc position for each group, starting the firs ilev=1)
      respsel <- which(object$var.agg == vlev[ilev])  # vlev vector with the names of aggregated categories
      
      # ntrep maximum nombre of documents to print, minumum between maxPrnDoc and the size of the group
      ntrep <- min(maxCharDoc, length(respsel))
      lisresult[[ilev]] <- data.frame(nlev,3)  # Nombre of the group
      # SourceTerm the rows for non aggregated table
      SourceTermcurs <- SourceTerm[respsel,,drop=FALSE ]   
      ly<-as.matrix(rep(0,ncol(DT)))
      # Rows are the vocabulary of after words
      rownames(ly)<-colnames(DT)
      if(is.null(motsval[[ilev]])) {  lisresult[ilev] <- NULL } else {
        # There are significant words for this category
        motsvalcurs<-as.matrix(motsval[[ilev]])
        motsvalcurs<-	  motsvalcurs[which(rownames(motsvalcurs)!="RestOfWords"),,drop=FALSE ]
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
          lisresult[[ilev]][i,1] <- rownames(object$SourceTerm)[respsel[ordrep[i]]]
          lisresult[[ilev]][i,2] <- repvaltest[ordrep[i]]
          lisresult[[ilev]][i,3] <- substr(corpus[respsel[ordrep[i]]], start = 1, stop = maxPrnDoc)
        }
        
        colnames(lisresult[[ilev]]) <- c("DOCUMENT", "CRITERION", "---------------------TEXT---------------------")	
        
      } # End  if(is.null(motsval[[ilev]]))
    } # End for ilev in 1:nlev
    names(lisresult) <- vlev
    res$CharDoc <- lisresult
  } # End Extraction of documents
  
  res$Proba <- proba
  class(res) <- c("LexChar", "list")  
return(res)
} # End function
 