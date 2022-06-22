#' @importFrom MASS ginv
#' @export
LexGalt <- function (object, context="ALL", conf.ellip=FALSE, 
                     nb.ellip = 100, graph= TRUE, axes = c(1, 2), label.group=NULL) 
{

  # ncp= number of dimensions kept to compute (by default NULL to indicate all dimensions)
  ncp <- NULL
  
  # LexGalt CAN'T work with Aggregated Lexical Tables (ALT)
  # level.ventil, a proportion corresponding to the level under which the category is ventilated; by default, 0 and no ventilation is done
  # If level.ventil, seed is a random value by defect and results may differ when running the same script twice.
  # Use the same values for the set.seed() function.
  
  # scale. Variables are are scaled to unit variance (standardized) (by default TRUE)
  scale <- TRUE
  level.ventil = 0; set.seed(1234)
  
  if(!is.null(object$var.agg)) stop("LexGalt needs non aggregate TextData objects")
  options(stringsAsFactors = FALSE)
  
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


    #######################################################################
  # 1.- Checking if it is a list (multiple) or only one object (simple)
  num.group <- ifelse(is.vector(object), length(object),1)
  

  if(num.group==1) {
    # Simple, only one TextData object
    # If it is not a TextData object: stop
    if (!inherits(object,"TextData")) stop("Object should be TextData class")
    if(is.null(context)) stop("No contextual variables are selected")
  } else {
    # Multiple. Check all objects in list are TextData objects
    for(i in 1:num.group) 
      if (!inherits(object[[i]],"TextData")) stop("Object should be TextData class")
  }
  #######################################################################
  
  #######################################################################
  ##2. Contextual variables
  context.quanti <- context.quali <- NULL
  #############  Start selecting contextual variables	
  if(length(context)==1)	# Only one contextual variable
    if(context=="ALL") {	# Contextual variable can be "ALL" (defect) (quanti+quali)
      context<- NULL	# New object context to separe quanti and quali variables
      if(num.group==1) {      # It is Simple  LexGalt
        if(!is.null(object$context$quali))     # Check if there are quali variables in TextData object
          context <- as.character(colnames(object$context$quali))
        if(!is.null(object$context$quanti)) # Check if there are quanti variables in TextData object
          context <- c(context,colnames(object$context$quanti))	# context has quanti and quali variables
      } else {    # It is Multiple LexGalt; Only read contextual variables from first TextData object
        if(!is.null(object[[1]]$context$quali)) 	
          context <- as.character(colnames(object[[1]]$context$quali))	
        if(!is.null(object[[1]]$context$quanti)) 
          context <- c(context,colnames(object[[1]]$context$quanti))	
      } # Final num.group
    }  # Final context=1	
  ##########################################################################3
  
  if(num.group ==1) {  # If it is simple
    if(!is.null(colnames(object$context$quanti))) # Selection quanti variables 
      context.quanti <- colnames(object$context$quanti)[which(colnames(object$context$quanti) %in% context)]
    if(!is.null(colnames(object$context$quali)))  # Selection quali variables
      context.quali <- colnames(object$context$quali)[which(colnames(object$context$quali) %in% context)]
  } else { # If is multiple
    context.quanti <- context.quali <- context
    # Only selection contextual variables if they are in all the groups
    for(i in 1:num.group) {
      if(!is.null(colnames(object[[i]]$context$quanti)))
        context.quanti <- colnames(object[[i]]$context$quanti)[which(colnames(object[[i]]$context$quanti) %in% context.quanti)]
      else   context.quanti <- NULL
      if(!is.null(colnames(object[[i]]$context$quali)))
        context.quali <- colnames(object[[i]]$context$quali)[which(colnames(object[[i]]$context$quali) %in% context.quali)]
      else   context.quali <- NULL
    } # End For
  } # Final  if(num.group ==1) 
  ##### Final of contextual variables
  # nQN is the number of selected quanti contextual variables
  # nQL is the number of selected quali contextual variables
  nQN <- ifelse(length(context.quanti)==0, 0, length(context.quanti))
  nQL <- ifelse(length(context.quali)==0, 0, length(context.quali))
  ################################################################################
  
  
  
  #######################################################################
  ##3. Functions
  mean.p<-function(V, poids) res <- sum(V * poids, na.rm = TRUE)/sum(poids[!is.na(V)])
  sd.p<-function(V, poids) res <- sqrt(sum(V^2 * poids, na.rm = TRUE)/sum(poids[!is.na(V)]))
  

  #######################################################################
  ## 4. Simple LexGalt 
  ####### If LexGalt is simple, not multiple
  if(num.group==1) {  # It is simple, not multiple
    
    #### There are Quantitative variables
    if(nQN!=0) {
      # df.quanti is a dataframe with quantitative variables in columns
      df.quanti <- data.frame(object$context$quanti[,context.quanti,drop=FALSE])
      # Check if one variable has all cases the same value (variance=0)
      # The following line returns a value for each variable with the range of each variable
      # If the range is 0 then all the vaues of the variable are the same
      fval <- apply(df.quanti, 2, function(x) max(x)-min(x))
      fval <- colnames(df.quanti)[which(fval==0)]
      if(length(fval)>0){
        warning("Quantitative variable with the same value are removed: ",fval)
        context.quanti <- context.quanti[!context.quanti %in% fval]
        nQN <- length(context.quanti) # Compute the new number of quanti variables
      }
      if(nQN<2) stop("The number of quantitative variables must be > 1")
      XQN <- data.frame(object$context$quanti[,context.quanti,drop=FALSE])	
      # minN is the maximum of factors that can be extracted
      minN <- min(nrow(XQN)- 1, ncol(XQN)) 
      # Si se especificó en la llamada el número de factores ncp
      ncpQN <- ifelse(is.null(ncp), minN, min(ncp, minN))
      # Save in two dataframes (XQN and XQN.initial) quanti variables
      XQN.initial<- XQN
    }
    #################### Final quantitative variables
  
        
    ######### Categorical variables
    if(nQL!=0){
      XQL <- data.frame(object$context$quali[,context.quali,drop=FALSE]) 
      minL <- min(nrow(XQL) - 1, ncol(Xtab.disjonctif(XQL)) - ncol(XQL))
      ncpQL <- ifelse(is.null(ncp), minL, min(ncp, minL))
      if(nQL<2) stop("The number of factors for qualitative variables ncp must be > 1")
      XQL.initial <- XQL
    }
    
    if(nQN+nQL==0) stop("There are not contextual variables")
    
    ##############################################
    # 5.  Starting 
    Y <- as.matrix(object$DocTerm) # (Documents x words) (Respondents x words)
    P <- as.matrix(Y/sum(Y))       # Proportion matrix 
    PI. <- apply(P, 1, sum)        # Margin vector % words for each document (respondent)
    P.J <- apply(P, 2, sum)        # Margin vector % words for each word
    # I respondents (docs), J words, K quantitative variables (nQN)
    ##############################################  
    
    
    ########### Start funct1 ############					
    funct1 <- function(phi.stand,ncpf,XP,su) {					
      # su: scale. TRUE if quantitative scaled 
      # su: is FALSE if qualitative or quantitative not scaled (centered)
      # P is not passed because is always the same 
      # P is % words in matrix documents x words
      # phi.stand have the eigenvectors (docs x dimensiones)
      # ncpf is the number of factors
      # XP is the quantitative data matrix centered (and scaled if su=TRUE)
      # L Matrix (words x dimensions)
      # L <- t(P) *phi.stand/P.J  
      L <- sweep(crossprod(P, phi.stand), 1, P.J, "/")					
      # mT (Words x categories):  (J words x I docs) * (I docs x variables o categorías) 
      mT <- crossprod(P, XP)	
      C <- crossprod(sweep(XP, 1, PI.^(1/2), "*"), sweep(XP, 1, PI.^(1/2),"*"))	
      # J words x L variables or categories
      W <- sweep(crossprod(t(mT), MASS::ginv(C)), 1, P.J, "/")	
      colnames(W) <- colnames(mT)
      # PCA using L as actif and W as supplementary
      diag.L <- PCA(cbind(L, W), quanti.sup = (ncpf + 1):(ncpf + ncol(mT)),
                    scale.unit = FALSE, ncp = ncpf, row.w = P.J, graph = FALSE)
      # Document coordinates					
      coord.ind <- sweep(crossprod(t(P), diag.L$svd$U), 1, PI., "/")
      # Cos2 of documents
      cos2.ind <- sweep(coord.ind^2, 1, apply(coord.ind^2, 1, sum),"/")
      # Building an object res.doc with coordinates and cos2 for the documents
      res.doc <- list(coord = coord.ind, cos2 = cos2.ind)
      # Building an object res.word with coordinates and cos2 for the docs
      res.word <- list(coord = diag.L$ind$coord, cos2 = diag.L$ind$cos2, 
                       contrib = diag.L$ind$contrib)
      res <- list(eig = diag.L$eig, doc = res.doc, word = res.word, diag.L=diag.L, L=L, W=W)
    }
    ########### Final of function funct1 ############					

    
    
    # =======  Call when LexGalt with quantitative variables			
    if(nQN!=0) {			
      diag.XQN <- PCA(XQN, scale.unit = scale, ncp = ncpQN, row.w = PI., graph = FALSE)
      # phi.stand.QN has the eigenvectors (Docs x dimensions)
      phi.stand.QN <- diag.XQN$svd$U
      # XQN matriz docs x quanti variables (centered or scale)  
      # XQN always is centered:
      XQN <- as.matrix(sweep(XQN, 2, apply(XQN, 2, mean.p, PI.), "-"))		
      # If scale=TRUE divided for standard deviation
      if (scale == "TRUE") XQN <- sweep(XQN, 2, apply(XQN, 2, sd.p, PI.), "/")	
      resQN <- funct1(phi.stand.QN, ncpQN, XQN, scale)
      resQN$quanti.var <- resQN$diag.L$quanti.sup	
      #  }
      class(resQN) <- c("LexGalt", "list")
    }			
    # =======  Final Call when LexGalt with quantitative variables	
    
    
    # =======  Call when LexGalt with qualitative variables			
    if(nQL!=0) {			
      diag.XQL <- MCA(XQL, row.w = PI., ncp = ncpQL, level.ventil = level.ventil, graph = FALSE)
      if (ncol(XQL) > 1) 
        XQL <- sweep(Xtab.disjonctif(XQL), 2, apply(Xtab.disjonctif(XQL), 2, mean.p, PI.), "-")
      else XQL <- Xtab.disjonctif(XQL)
      phi.stand.QL <- diag.XQL$svd$U
      # El resultado de funct1 var.contexto cualitativo se almacena en resQL
      resQL <- funct1(phi.stand.QL, ncpQL, XQL, FALSE)			
      resQL$quali.var <- list(coord = resQL$diag.L$quanti.sup$coord, cos2 = resQL$diag.L$quanti.sup$cos2)
      class(resQL) <- c("LexGalt", "list")
    }
    # ======= FInal  Call when LexGalt with qualitative variables			
    
    

    #______________________________ Ellipses ---------------------
    if(conf.ellip==TRUE)
      if(nb.ellip>0) {
        PJ <- vector(mode = "list", length = nb.ellip)
        #     stop("Hace elipses")
        
        for(n in 1:nb.ellip) {
          samp <- sample(1:nrow(Y), replace = TRUE)
          while (sum(apply(Y[samp, ], 2, sum) > 0) != ncol(Y)) samp <- sample(1:nrow(Y), replace = TRUE)
          Y.samp <- Y[samp, ]
          P.samp <- as.matrix(Y.samp/sum(Y.samp))
          PI.samp <- apply(P.samp, 1, sum)
          P.J.samp <- apply(P.samp, 2, sum)
          PJ[[n]] <- P.J.samp
          
          if(nQL>0) {  X.samp <- sweep(Xtab.disjonctif(XQL.initial[samp, ]),2, 
                                       apply(Xtab.disjonctif(XQL.initial[samp, ]), 2, mean.p, PI.samp), "-")
          phi.stand.samp <- phi.stand.QL[samp, ]
          T.samp <- crossprod(P.samp, X.samp)
          C.samp <- crossprod(sweep(X.samp, 1, PI.samp^(1/2), "*"), sweep(X.samp, 1, PI.samp^(1/2), "*"))
          
          if (n == 1) {  
            LQ.samp <- sweep(crossprod(P.samp, phi.stand.samp), 1, P.J.samp, "/")
            WQ.samp <- sweep(crossprod(t(T.samp), MASS::ginv(C.samp)), 1, P.J.samp, "/")
            WQ.stand.samp <- sweep(sweep(crossprod(t(T.samp), 
                                                   MASS::ginv(C.samp)), 1, P.J.samp, "/"), 2, apply(sweep(crossprod(t(T.samp), 
                                                   MASS::ginv(C.samp)), 1, P.J.samp, "/"), 2, sd.p, P.J.samp), "/")
          }
          else {       
            LQ.samp <- rbind(LQ.samp, sweep(crossprod(P.samp, phi.stand.samp), 1, P.J.samp, "/"))
            WQ.samp <- cbind(WQ.samp, sweep(crossprod(t(T.samp),MASS::ginv(C.samp)), 1, P.J.samp, "/"))
            WQ.stand.samp <- rbind(WQ.stand.samp, sweep(sweep(crossprod(t(T.samp), 
                          MASS::ginv(C.samp)), 1, P.J.samp, "/"), 2, apply(sweep(crossprod(t(T.samp),
                          MASS::ginv(C.samp)), 1, P.J.samp, "/"), 2, sd.p, P.J.samp), "/"))
          }
          }  # End nQL>0
          
          if(nQN>0) {
            X.samp <- as.matrix(sweep(XQN.initial[samp, ], 2, apply(XQN.initial[samp,], 2, 
                                                                    mean.p, PI.samp), "-"))
            if(scale==TRUE) 
              X.samp <- sweep(X.samp, 2, apply(X.samp, 2, sd.p, PI.samp), "/")
            phi.stand.samp <- phi.stand.QN[samp, ]
            T.samp <- crossprod(P.samp, X.samp)
            C.samp <- crossprod(sweep(X.samp, 1, PI.samp^(1/2), "*"), sweep(X.samp, 1, PI.samp^(1/2), "*"))
            if (n == 1) {
              LN.samp <- sweep(crossprod(P.samp, phi.stand.samp), 1, P.J.samp, "/")
              WN.samp <- sweep(crossprod(t(T.samp), MASS::ginv(C.samp)), 1, P.J.samp, "/")
              WN.stand.samp <- sweep(sweep(crossprod(t(T.samp), 
                              MASS::ginv(C.samp)), 1, P.J.samp, "/"), 2, apply(sweep(crossprod(t(T.samp), 
                              MASS::ginv(C.samp)), 1, P.J.samp, "/"), 2, sd.p, P.J.samp), "/")
            }
            else {
              LN.samp <- rbind(LN.samp, sweep(crossprod(P.samp, phi.stand.samp), 1, P.J.samp, "/"))
              WN.samp <- cbind(WN.samp, sweep(crossprod(t(T.samp),MASS::ginv(C.samp)), 1, P.J.samp, "/"))
              WN.stand.samp <- rbind(WN.stand.samp, sweep(sweep(crossprod(t(T.samp), 
                              MASS::ginv(C.samp)), 1, P.J.samp, "/"), 2, apply(sweep(crossprod(t(T.samp),
                                                                                                                                           MASS::ginv(C.samp)), 1, P.J.samp, "/"), 2, sd.p, P.J.samp), "/"))
            }
          } # End if(nQN>0) 
        } ### End for nb.ellip  
          
          
        
        if(nQL>0)  {
          rownames(LQ.samp) <- paste(rep(rownames(resQL$L), nb.ellip), rep(1:nb.ellip, each = nrow(resQL$L)), sep = "")
          freq.ellip.coord <- as.data.frame(PCA(rbind(resQL$L, LQ.samp), ncp = ncpQL, 
                                                ind.sup = (ncol(Y) + 1):((nb.ellip + 1) * ncol(Y)), row.w = P.J, 
                                                scale.unit = FALSE, graph = FALSE)$ind.sup$coord)
          colnames(WQ.samp) <- paste(rep(colnames(resQL$W), nb.ellip), rep(1:nb.ellip, each = ncol(resQL$W)), sep = "")
          
          var.ellip.coord <- as.data.frame(PCA(cbind(resQL$L, WQ.samp), 
                                               quanti.sup = (ncpQL + 1):(ncpQL + ncol(WQ.samp)), scale.unit = FALSE, 
                                               ncp = ncpQL, row.w = P.J, graph = F)$quanti.sup$coord)
          freq.ellip.coord$FREQ <- rep(rownames(resQL$L), nb.ellip)
          var.ellip.coord$VAR <- rep(colnames(XQL), nb.ellip)
          resQL$ellip <- list(word = freq.ellip.coord, var = var.ellip.coord)
          class(resQL) <- c("LexGalt", "list")
        } # End if(nQL>0) 
        
        if(nQN>0)  {
          rownames(LN.samp) <- paste(rep(rownames(resQN$L), nb.ellip), rep(1:nb.ellip, each = nrow(resQN$L)), sep = "")
          freq.ellip.coord <- as.data.frame(PCA(rbind(resQN$L, LN.samp), ncp = ncpQN, 
                                                ind.sup = (ncol(Y) + 1):((nb.ellip + 1) * ncol(Y)), row.w = P.J, 
                                               scale.unit = FALSE, graph = FALSE)$ind.sup$coord)
            for (n in 1:nb.ellip) {
            aux <- as.matrix(freq.ellip.coord[(((n - 1) * ncol(Y)) + 1):(n * ncol(Y)), 1:ncpQN])
            aux.cent <- sweep(aux, 2, apply(aux, 2, mean.p, PJ[[n]]), "-")
            aux.stand <- as.matrix(sweep(aux.cent, 2, apply(aux.cent, 2, sd.p, PJ[[n]]), "/"))
            if (n == 1) 
              var.ellip.coord <- as.data.frame(crossprod(WN.stand.samp[(((n -1) * ncol(Y)) + 1):(n * ncol(Y)), ],
                                                         sweep(aux.stand, 1, PJ[[n]], "*")))
            else var.ellip.coord <- rbind(var.ellip.coord, as.data.frame(crossprod(WN.stand.samp[(((n -1)
                                 * ncol(Y)) + 1):(n * ncol(Y)), ], sweep(aux.stand, 1, PJ[[n]], "*"))))
                         }  # End For
          
          freq.ellip.coord$FREQ <- rep(rownames(resQN$L), nb.ellip)
          var.ellip.coord$VAR <- rep(colnames(XQN), nb.ellip)
          resQN$ellip <- list(word = freq.ellip.coord, var = var.ellip.coord)
        } # End if(nQN>0)
        
        #______________________________ Final Ellipses ---------------------   
      }

    if(nQL==0) resQL <- ""
    if(nQN==0) resQN <- ""
    
    resNew <- list(SQL = resQL, SQN=resQN)
    class(resNew) <- c("LexGalt","list") 
    #______________________________ Plots simple ---------------------
    if (graph) {
      if(nQN!=0) { 
        plot.LexGalt(resNew, type="QN", eigen=TRUE, new.plot = TRUE, title="Eigenvalues from Quantitative variables")
        plot.LexGalt(resNew, type="QN", selDoc="ALL",axes = axes, new.plot = TRUE, title="Documents from Quantitative variables")
        plot.LexGalt(resNew, type="QN", selQuantiVar="ALL", axes = axes, conf.ellip = conf.ellip, new.plot = TRUE)
        
        if (nrow(resNew$SQN$word$coord) < 50) {
          plot.LexGalt(resNew, type="QN", selWord="ALL", autoLab="yes", axes = axes, conf.ellip = conf.ellip, new.plot = TRUE,
                       title="Words from Quantitative variables")
        }
        else {
          plot.LexGalt(resNew, type="QN", selWord="coord 49", autoLab="yes", axes = axes, conf.ellip = conf.ellip, new.plot = TRUE,
                       title="Words from Quantitative variables")
          warning("The first 50 words that have the highest coordinates on the 2 dimensions of your plot are drawn.")
        }
      } # End nQN
      if(nQL>0) { 
        plot.LexGalt(resNew, type="QL", selDoc="ALL", axes = axes, new.plot = TRUE, title="Documents from Qualitative variables")
        plot.LexGalt(resNew, type="QL", selQualiVar="ALL", axes = axes, new.plot = TRUE, title="Qualitative variables")
        
        if (nrow(resNew$SQL$word$coord) < 50) {
          plot.LexGalt(resNew, type="QL", selWord="ALL", conf.ellip = conf.ellip, new.plot = TRUE,
                       title="Words from Qualitative variables")
        }
        else {
          plot.LexGalt(resNew, type="QL", selWord="coord 49",autoLab="yes", conf.ellip = conf.ellip, new.plot = TRUE,
                       title="Words from Qualitative variables")
          warning("The first 50 words that have the highest coordinates on the 2 dimensions of your plot are drawn.")     
        }
      } # End nQL
    } # End graph
    
    #___________________________ Final plots
    
    return(resNew)
    
  } # Final of if(num.group==1), simple
  
   else {
          ############################  MULTIPLE

     Y<-vector(mode='list',length=num.group)
     for(i in 1:num.group) Y[[i]] <- as.matrix(object[[i]]$DocTerm)
     
     if(is.null(label.group)) label.group <- paste("GROUP",1:num.group,sep = ".")	
     if(length(label.group) != num.group) stop("The name of groups is != number of groups")
     name.group <- label.group
     
     N<-sum(unlist(Y))
     Nl<-sapply(Y,sum)
     num.ind<-sapply(Y,nrow)
     num.freq<-sapply(Y,ncol)
     names(num.freq)<-names(num.ind)<-names(Nl)<-label.group
     Pl<-mapply(function(y,nl) y/nl,Y,Nl,SIMPLIFY = FALSE)
     pi.l<-lapply(Pl,function(p) apply(p,1,sum))
     p.jl<-lapply(Pl,function(p) apply(p,2,sum))
     wl<-Nl/N
     
     DIl<-mapply(function(pi,wl) pi*wl,pi.l,wl,SIMPLIFY = FALSE)
     DJl<-mapply(function(pj,wl) pj*wl,p.jl,wl,SIMPLIFY = FALSE)
     DI<-unlist(DIl)
     DJ<-unlist(DJl)
     
     nXQL <- nXQN <- 0 

          
     if(nQL!=0) {
       XQL<-vector(mode='list',length=num.group)
       nXQL <- length(context.quali) 
       if(nXQL<2) stop("The number of qualitative variables must be > 1")
       for(i in 1:num.group) XQL[[i]] <- as.matrix(object[[i]]$context$quali[context.quali])
     }
     Xl<-vector(mode='list',length=num.group)
     XG<-data.frame()
     resQL <- resQN <- NULL 
     
     if(nXQL!=0) {
       for(i in 1:num.group) Xl[[i]]<-as.matrix(Xtab.disjonctif(XQL[[i]]))
       if(unique(sapply(XQL,ncol))!=1) Xl<-mapply(function(x,pi) sweep(x,2,apply(x,2,mean.p,pi),"-"),Xl,pi.l,SIMPLIFY = FALSE)
       
       for (i in 1:num.group) XG<-as.matrix(rbind(XG,Xl[[i]]))
       Yal<-mapply(function(y,x) crossprod(y,x),Y,Xl,SIMPLIFY = FALSE)
       Pal<-mapply(function(ya) ya/N,Yal,SIMPLIFY = FALSE)
       C<-crossprod(sweep(XG,1,DI^(1/2),"*"),sweep(XG,1,DI^(1/2),"*"))
       Ml<-mapply(function(di,dj,w,x) diag(dj)%*%matrix(nrow=length(dj),ncol=length(di),1/w)%*%diag(di)%*%x,DIl,DJl,wl,Xl,SIMPLIFY = FALSE)				
       Zl<-mapply(function(pa,m,dj) sweep(pa-m,1,dj,"/")%*%MASS::ginv(C),Pal,Ml,DJl,SIMPLIFY=FALSE)
       ponde<-mapply(function(z,dj) Re(eigen(t(z)%*%diag(dj)%*%z%*%C)$values[1]),Zl,DJl,SIMPLIFY = FALSE)				
       row.w<-unlist(mapply(function(dj,pond) dj/pond,DJl,ponde,SIMPLIFY = FALSE))				
       Z<-data.frame()				
       
       for (i in 1:num.group) Z<-as.matrix(rbind(Z,Zl[[i]]))				
       diag.MFAGALT<-eigen(t(Z)%*%diag(row.w)%*%Z%*%C)				
       eig<-Re(diag.MFAGALT$values)				
       eig<-eig[1:(length(eig)-unique(sapply(XQL,ncol)))]				
       vp<-as.data.frame(matrix(NA,length(eig),3))				
       rownames(vp)<-paste("comp",1:length(eig))				
       colnames(vp)<- c("eigenvalue","percentage of variance","cumulative percentage of variance")				
       vp[,"eigenvalue"]<-eig				
       vp[,"percentage of variance"]<-(eig/sum(eig))*100				
       vp[,"cumulative percentage of variance"]<-cumsum(vp[,"percentage of variance"])	
       U<-sweep(Re(diag.MFAGALT$vectors[,1:length(eig)]),2,sqrt(diag(t(Re(diag.MFAGALT$vectors[,1:length(eig)]))%*%C%*%Re(diag.MFAGALT$vectors[,1:length(eig)]))),"/")				
       coord.freq<-Z%*%C%*%U				
       contrib.freq<-sweep(sweep(coord.freq^2,1,row.w,FUN="*")*100,2,eig,FUN="/")				
       dist2.freq<-diag(Z%*%C%*%t(Z))				
       cos2.freq<-sweep(as.matrix(coord.freq^2),1,dist2.freq,FUN="/")				
       coord.var<-sweep(U,2,sqrt(eig),"*")				
       dist2.var<-apply(sweep(Z,1,sqrt(row.w),FUN="*")^2,2,sum)				
       cos2.var<-sweep(as.matrix(coord.var^2),1,dist2.var,FUN="/")	
       coord.freq.partiel<-mapply(function(z) z%*%C%*%U,Zl,SIMPLIFY=FALSE)				
       coord.var.partiel<-mapply(function(z,f,dj,pond) sweep(t(z)%*%sweep(f,1,dj/pond,"*")*num.group,2,sqrt(eig),"/"),Zl,coord.freq.partiel,DJl,ponde,SIMPLIFY=FALSE)				
       tab<-(Z%*%C%*%t(Z))^2				
       tab<-sweep(sweep(tab,2,row.w,"*"),1,row.w,"*")				
       Lg<-matrix(0,num.group+1,num.group+1)				
       ind.gl<-0				

     
       for (gl in 1:num.group) {				
         ind.gc<-0			
         for (gc in 1:num.group) {			
           Lg[gl,gc]<-Lg[gc,gl]<-sum(tab[(ind.gl+1):(ind.gl+num.freq[gl]),(ind.gc+1):(ind.gc+num.freq[gc])])		
           ind.gc<-ind.gc+num.freq[gc]			
         }			
         Lg[num.group+1,gl]<-Lg[gl,num.group+1]<-sum(tab[(ind.gl+1):(ind.gl+num.freq[gl]),1:ncol(tab)])/eig[1]			
         ind.gl<-ind.gl+num.freq[gl]			
       }	
       Lg[num.group+1,num.group+1]<-sum(tab[1:ncol(tab),1:ncol(tab)])/(eig[1]^2)				
       RV<-sweep(Lg,2,sqrt(diag(Lg)),"/")				
       RV<-sweep(RV,1,sqrt(diag(Lg)),"/")				
       contrib.group<-matrix(NA,num.group,ncol(U))				
       dist2.group<-vector(length=num.group)				
       freq.gr<-0
       
       for (g in 1:num.group) {				
         if (g==1) contrib.group[g,]<-apply(contrib.freq[1:num.freq[g],],2,sum)			
         else contrib.group[g,]<-apply(contrib.freq[(freq.gr+1):(freq.gr+num.freq[g]),],2,sum)			
         ponde.tot<-eigen(t(Zl[[g]])%*%diag(DJl[[g]])%*%Zl[[g]]%*%C)$values	
         dist2.group[g]<-sum((ponde.tot/ponde[[g]])^2)			
         freq.gr<-freq.gr+num.freq[g]			
       }
       coord.group<-sweep(contrib.group/100,2,eig,"*")				
       cos2.group<-sweep(coord.group^2,1,dist2.group,"/")				
       rownames(coord.var)<-rownames(cos2.var)<-colnames(XG)					
       colnames(coord.freq)<-colnames(contrib.freq)<-colnames(cos2.freq)<-colnames(coord.var)<-colnames(cos2.var)<-paste("Dim",c(1:ncol(U)),sep = ".")
       
       names(coord.var.partiel)<-label.group	
       
       for(i in 1:num.group){				
         rownames(coord.var.partiel[[i]])<-colnames(XG)			
         colnames(coord.var.partiel[[i]])<-paste("Dim",c(1:ncol(U)),sep = ".")			
       }
       rownames(coord.group)<-rownames(contrib.group)<-rownames(cos2.group)<-label.group				
       colnames(coord.group)<-colnames(contrib.group)<-colnames(cos2.group)<-paste("Dim",c(1:ncol(U)),sep = ".")				
       rownames(Lg)<-colnames(Lg)<-rownames(RV)<-colnames(RV)<-c(label.group,"MFA-GALT")				
       res<-list()				
       res$eig<-vp				
       res$word<-list(coord=coord.freq,cos2=cos2.freq,contrib=contrib.freq)

       type <- "n"
       res$quali.var<-list(coord=coord.var,cos2=cos2.var,coord.partial=coord.var.partiel)	
       res$group<-list(Lg=Lg,RV=RV,coord=coord.group,contrib=contrib.group,cos2=cos2.group)				
       res$call<-list(num.groups=num.group,name.groups=label.group,num.freq=num.freq,type=type)				
       
       class(res) <- c("LexGalt", "list")	# Remove MFAGALT to the end	
       resQL <- res	
       class(resQL) <- c("LexGalt", "list")		
     } # FInal     if(nXQL!=0) 
     
     if(nQN!=0) {
       XQN<-vector(mode='list',length=num.group)
       nXQN <- length(context.quanti)
       if(nXQN<2) stop("The number of quantitative variables must be > 1")
       for(i in 1:num.group) {
         XQN[[i]] <- as.matrix(object[[i]]$context$quanti[context.quanti])
         Xl[[i]]<-as.matrix(XQN[[i]])}
       Xl<-mapply(function(x,pi) sweep(x,2,apply(x,2,mean.p,pi),"-"),Xl,pi.l,SIMPLIFY = FALSE)
       XG<-data.frame()
       
       for (i in 1:num.group) XG<-as.matrix(rbind(XG,Xl[[i]]))
       Yal<-mapply(function(y,x) crossprod(y,x),Y,Xl,SIMPLIFY = FALSE)
       Pal<-mapply(function(ya) ya/N,Yal,SIMPLIFY = FALSE)
       C<-crossprod(sweep(XG,1,DI^(1/2),"*"),sweep(XG,1,DI^(1/2),"*"))
       Ml<-mapply(function(di,dj,w,x) diag(dj)%*%matrix(nrow=length(dj),ncol=length(di),1/w)%*%diag(di)%*%x,DIl,DJl,wl,Xl,SIMPLIFY = FALSE)				
       Zl<-mapply(function(pa,m,dj) sweep(pa-m,1,dj,"/")%*%MASS::ginv(C),Pal,Ml,DJl,SIMPLIFY=FALSE)
       ponde<-mapply(function(z,dj) Re(eigen(t(z)%*%diag(dj)%*%z%*%C)$values[1]),Zl,DJl,SIMPLIFY = FALSE)				
       row.w<-unlist(mapply(function(dj,pond) dj/pond,DJl,ponde,SIMPLIFY = FALSE))				
       Z<-data.frame()				
       for (i in 1:num.group) Z<-as.matrix(rbind(Z,Zl[[i]]))				
       diag.MFAGALT<-eigen(t(Z)%*%diag(row.w)%*%Z%*%C)				
       eig<-Re(diag.MFAGALT$values)				
       vp<-as.data.frame(matrix(NA,length(eig),3))				
       rownames(vp)<-paste("comp",1:length(eig))				
       colnames(vp)<- c("eigenvalue","percentage of variance","cumulative percentage of variance")				
       vp[,"eigenvalue"]<-eig				
       vp[,"percentage of variance"]<-(eig/sum(eig))*100				
       vp[,"cumulative percentage of variance"]<-cumsum(vp[,"percentage of variance"])
       
       U<-sweep(Re(diag.MFAGALT$vectors[,1:length(eig)]),2,sqrt(diag(t(Re(diag.MFAGALT$vectors[,1:length(eig)]))%*%C%*%Re(diag.MFAGALT$vectors[,1:length(eig)]))),"/")				
       coord.freq<-Z%*%C%*%U				
       contrib.freq<-sweep(sweep(coord.freq^2,1,row.w,FUN="*")*100,2,eig,FUN="/")				
       dist2.freq<-diag(Z%*%C%*%t(Z))				
       cos2.freq<-sweep(as.matrix(coord.freq^2),1,dist2.freq,FUN="/")				
       coord.var<-sweep(U,2,sqrt(eig),"*")				
       dist2.var<-apply(sweep(Z,1,sqrt(row.w),FUN="*")^2,2,sum)				
       cos2.var<-sweep(as.matrix(coord.var^2),1,dist2.var,FUN="/")	
       coord.freq.partiel<-mapply(function(z) z%*%C%*%U,Zl,SIMPLIFY=FALSE)				
       coord.var.partiel<-mapply(function(z,f,dj,pond) sweep(t(z)%*%sweep(f,1,dj/pond,"*")*num.group,2,sqrt(eig),"/"),Zl,coord.freq.partiel,DJl,ponde,SIMPLIFY=FALSE)				
       # Si es numerica
       cor.var<-sweep(sweep(t(Z)%*%sweep(coord.freq,1,row.w,"*"),1,apply(Z,2,function(V, poids) res <- sqrt(sum(V^2 * poids, na.rm = TRUE)),row.w),"/"),2,sqrt(eig),"/")
       tab<-(Z%*%C%*%t(Z))^2				
       tab<-sweep(sweep(tab,2,row.w,"*"),1,row.w,"*")				
       Lg<-matrix(0,num.group+1,num.group+1)				
       ind.gl<-0				
       for (gl in 1:num.group) {				
         ind.gc<-0			
         for (gc in 1:num.group) {			
           Lg[gl,gc]<-Lg[gc,gl]<-sum(tab[(ind.gl+1):(ind.gl+num.freq[gl]),(ind.gc+1):(ind.gc+num.freq[gc])])		
           ind.gc<-ind.gc+num.freq[gc]			
         }			
         Lg[num.group+1,gl]<-Lg[gl,num.group+1]<-sum(tab[(ind.gl+1):(ind.gl+num.freq[gl]),1:ncol(tab)])/eig[1]			
         ind.gl<-ind.gl+num.freq[gl]			
       }				
       
       Lg[num.group+1,num.group+1]<-sum(tab[1:ncol(tab),1:ncol(tab)])/(eig[1]^2)				
       RV<-sweep(Lg,2,sqrt(diag(Lg)),"/")				
       RV<-sweep(RV,1,sqrt(diag(Lg)),"/")				
       contrib.group<-matrix(NA,num.group,ncol(U))				
       dist2.group<-vector(length=num.group)				
       freq.gr<-0				

       
       for (g in 1:num.group) {				
         if (g==1) contrib.group[g,]<-apply(contrib.freq[1:num.freq[g],],2,sum)			
         else contrib.group[g,]<-apply(contrib.freq[(freq.gr+1):(freq.gr+num.freq[g]),],2,sum)			
         ponde.tot<-eigen(t(Zl[[g]])%*%diag(DJl[[g]])%*%Zl[[g]]%*%C)$values			
         dist2.group[g]<-sum((ponde.tot/ponde[[g]])^2)			
         freq.gr<-freq.gr+num.freq[g]			
       }
       
       
       coord.group<-sweep(contrib.group/100,2,eig,"*")				
       cos2.group<-sweep(coord.group^2,1,dist2.group,"/")				
       rownames(coord.var)<-rownames(cos2.var)<-colnames(XG)					
       colnames(coord.freq)<-colnames(contrib.freq)<-colnames(cos2.freq)<-colnames(coord.var)<-colnames(cos2.var)<-paste("Dim",c(1:ncol(U)),sep = ".")
       rownames(cor.var)<-colnames(XG)
       colnames(cor.var)<-paste("Dim",c(1:ncol(U)),sep = ".")
       
       
       names(coord.var.partiel)<-label.group				
       for(i in 1:num.group){				
         rownames(coord.var.partiel[[i]])<-colnames(XG)			
         colnames(coord.var.partiel[[i]])<-paste("Dim",c(1:ncol(U)),sep = ".")			
       }				

       rownames(coord.group)<-rownames(contrib.group)<-rownames(cos2.group)<-label.group				
       colnames(coord.group)<-colnames(contrib.group)<-colnames(cos2.group)<-paste("Dim",c(1:ncol(U)),sep = ".")				
       rownames(Lg)<-colnames(Lg)<-rownames(RV)<-colnames(RV)<-c(label.group,"MFA-GALT")				
       res<-list()				
       res$eig<-vp				
       res$word<-list(coord=coord.freq,cos2=cos2.freq,contrib=contrib.freq)	    		
       res$quanti.var<-list(coord=coord.var,cor=cor.var,cos2=cos2.var,coord.partial=coord.var.partiel)
       res$group<-list(Lg=Lg,RV=RV,coord=coord.group,contrib=contrib.group,cos2=cos2.group)
       #   \item{type}{the type of variables: "c" or "s" for quantitative variables and "n" for categorical variables. The difference is that for "s" variables are scaled to unit variance (by default, variables are scaled to unit variance)}
       
       type <- ifelse(scale==TRUE,"c","s")
       res$call<-list(num.groups=num.group,name.groups=label.group,num.freq=num.freq,type=type)
       
       class(res) <- c("LexGalt", "list")	# Remove MFAGALT	
       resQN <- res	
       class(resQN) <- c("LexGalt", "list")
     }
     
     resNew <- list(MQL=resQL,MQN=resQN)
     class(resNew) <- c("LexGalt", "list") 
 
     
     if(graph) {
       if(nQL!=0){
         plot.LexGalt(resNew, type="QL", eigen=TRUE, new.plot = TRUE, title="Eigenvalues from Quantitative variables")
         plot.LexGalt(resNew, type="QL",plot.group="TRUE",new.plot = TRUE, title="Groups representation (LexGalt) from qualitative variables")
  #       plot.LexGalt(resNew, type="QL", selDoc="ALL",axes = axes, new.plot = TRUE, title="Documents from Qualitative variables")
         plot.LexGalt(resNew, type="QL", selQualiVar="ALL", axes = axes, new.plot = TRUE, title="Qualitative variables")
         plot.LexGalt(resNew, type="QL", selWord="ALL", new.plot = TRUE,
                      title="Words from Qualitative variables")
       }
       
       if(nQN!=0){
         plot.LexGalt(resNew, type="QN", eigen=TRUE, new.plot = TRUE, title="Eigenvalues from Quantitative variables")
                  plot.LexGalt(resNew, type="QN",plot.group="TRUE", new.plot = TRUE,
                      title="Groups representation (LexGalt) from quantitative variables")
         plot.LexGalt(resNew, type="QN", selDoc="ALL",axes = axes, new.plot = TRUE, title="Documents from Quantitative variables")
         plot.LexGalt(resNew , type="QN", selQuantiVar="ALL", axes = axes, new.plot = TRUE, partial=TRUE, title="Quantitative variables")
         plot(resNew, type="QN", selWord="ALL", new.plot=TRUE, title="Words from Quantitative variables")
       } # FINAL nQN
     } # Final graph
     
     return(resNew)
   } # Final multiple or simple, if(num.group==1)
  
} # Final LexGalt