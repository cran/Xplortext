#' @importFrom grDevices adjustcolor col2rgb rgb
#' @importFrom graphics arrows points
#' @export
plot.LexGalt <- function(x, type="QL", selDoc=NULL, selWord=NULL, selQualiVar=NULL, 
  selQuantiVar=NULL, conf.ellip=FALSE, selWordEllip=NULL, selQualiVarEllip=NULL,
  selQuantiVarEllip=NULL, level.conf=0.95, eigen=FALSE, title = NULL, axes = c(1, 2), 
  xlim = NULL, ylim = NULL, col.eig="grey", col.doc = "black", col.word = NULL, 
  col.quali = "blue", 
  col.quanti = "blue",col="grey", pch = 20, label = TRUE, 
  autoLab = c("auto", "yes", "no"),  palette = NULL, unselect = 1,  
  selCov=FALSE, selGroup="ALL", partial=FALSE, plot.group=FALSE,
  col.group=NULL, label.group=NULL, legend=TRUE, pos.legend="topleft", new.plot = TRUE,cex=1, ...)
{
  

    if(!inherits(x,"LexGalt")) stop("non convenient data")
    if(is.null(selDoc) & is.null(selWord) & is.null(selQualiVar)&
     is.null(selQuantiVar) & conf.ellip==FALSE & is.null(selWordEllip) & is.null(selQualiVarEllip) &
     is.null(selQuantiVarEllip) & eigen==FALSE) stop("There is nothing to plot\nSelect some argument")
  
     if(unselect==TRUE) unselect <- 1
#######################  Previous control
    if (is.numeric(unselect)) 
    if ((unselect > 1) | (unselect < 0)) 
      stop("unselect should be between 0 and 1")
  
  somethingsel <- all(sapply(list(selDoc, selWord, selQualiVar, selQuantiVar, selWordEllip, selQualiVarEllip, selQuantiVarEllip,
                                    selGroup, eigen), is.null))
  if(somethingsel) stop("There is nothing selected to plot")

  
####################### Functions
  selX <- function(sel1)
  {
    xx <- gregexpr(pattern =' ',sel1)[[1]][1]-1
    xx <- substr(sel1, 1, xx)
  }
  autoLab <- match.arg(autoLab, c("auto", "yes", "no"))
  if (autoLab == "yes") 
    autoLab = TRUE
  if (autoLab == "no") 
    autoLab = FALSE
  auto.Lab <- autoLab
  

    PALETTE <- palette(c("black", "red", "green3", "blue", "cyan", 
                         "magenta", "darkgray", "darkgoldenrod", "darkgreen", 
                         "violet", "turquoise", "orange", "lightpink", 
                         "lavender", "yellow", "lightgreen", "lightgrey", 
                         "lightblue", "darkkhaki", "darkmagenta", "darkolivegreen", 
                         "lightcyan", "darkorange", "darkorchid", "darkred", 
                         "darksalmon", "darkseagreen", "darkslateblue", 
                         "darkslategray", "darkslategrey", "darkturquoise", 
                         "darkviolet", "lightgray", "lightsalmon", "lightyellow", 
                         "maroon")) 
    if (!is.null(palette)) PALETTE[1:length(palette)] <- palette
                   
  
  selectionX <- function(sel1, xobj, bType, axx, axy)
  {
    if(is.null(sel1)) return(sel1)
    xx <- ""
    if(length(sel1)==1){
      if(sel1=="ALL") sel1 <- c(1:dim(xobj$coord)[1])
    }
    if(length(sel1)==1) {
      if(sel1=="meta") sel1 <- "meta 3"
      # if(sel1=="char") sel1 <- "char 0.05"
      xx <- gregexpr(pattern =' ',sel1)[[1]][1]-1
      xx <- substr(sel1, 1, xx)
    }  }
########## End Functions

  xlimT <- xlim ; ylimT <- ylim


  # Check if it is a Multiple analysis
  typeSM <- ifelse((length(x$SQL)+ length(x$SQN)) > 0, "typeS", "typeM")
  
  if(typeSM=="typeS") {
    if(type=="QL") {        # CATEGORICAL VARIABLES
      if(length(x$SQL)==1) stop("There are not categorical variables selected in LexGalt object")
      res.cagalt <- x$SQL }
    if(type=="QN") {
      if(length(x$SQN)==1)stop("There are not quantitative variables selected in LexGalt object")
      res.cagalt <- x$SQN }
    lab.x <- paste("Dim ", axes[1], " (", format(res.cagalt$eig[axes[1],2],
                                                 nsmall = 2, digits = 2), "%)", sep = "")
    lab.y <- paste("Dim ", axes[2], " (", format(res.cagalt$eig[axes[2],2],
                                                 nsmall = 2, digits = 2), "%)", sep = "")
    

  
###################################### Plot eigenvalues   
  if(eigen) {
    #  nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")) # TRUE  when inside RStudio ,  # FALSE when outside RStudio
    # Cambiado
    if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {
      dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
    } else dev.new()
    
    if(is.null(title)) titleE <- "Eigenvalues" else titleE <- title
    barplot(res.cagalt$eig[, 1], main = titleE, col=col.eig, cex.axis=cex, cex.names=cex,
            names.arg = paste("dim",1:nrow(res.cagalt$eig)))}
###################################### End Plot eigenvalues 


    
###################################### Document selection selDoc. Only for separate analysis. Answers if open question
  if (!is.null(selDoc)) {
    if (is.null(title)) 
      titre <- "Document factor map (LexGalt)"
    else titre <- title
    coord.doc <- res.cagalt$doc$coord[, axes, drop = FALSE]
    num.doc <- nrow(coord.doc)
    xlim <- xlimT; ylim <- ylimT
    if (is.null(xlim)) xlim <- c(min(coord.doc[, 1]), max(coord.doc[, 1])) * 1.2
    if (is.null(ylim)) ylim <- c(min(coord.doc[, 2]), max(coord.doc[, 2])) * 1.2
    selection <- NULL
    rdo<-""
    
    if (mode(selDoc) == "numeric") selection <- selDoc
    else if(length(selDoc)==1) {
      if (selDoc=="ALL")  selection <- c(1:num.doc)
      else rdo <- selX(selDoc)
    } else { selection <-  which(rownames(coord.doc) %in% selDoc)}
  
    
    
    if (rdo=="coord") 
      selection <- (rev(order(apply(coord.doc^2, 1, max))))[1:min(num.doc,
                 sum(as.integer(unlist(strsplit(selDoc, "coord"))), na.rm = T))]
    if (rdo=="cos2") {
      selcos2 <- as.numeric(substr(selDoc, 5, nchar(selDoc)))
      selection <- which(apply(res.cagalt$doc$cos2[, axes], 1, sum)> selcos2)
    }
    if (rdo=="contrib") stop("There is not possible to use contrib selection for selDoc")
    
    if(is.null(selection)) if(rdo=="") if(length(selDoc)==1) selection <- which(rownames(coord.doc) %in% selDoc)
    if(length(selection)==0) stop("There are not selected elements to plot")
    
    if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {  
      dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
   } else dev.new() 

    plot(0, 0, main = titre, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, col = "white", asp = 1, ...)
    abline(v = 0, lty = 2, ...)
    abline(h = 0, lty = 2, ...)
 ipch <- 20
    if(!is.null(pch)) {                                     # symbols in graphs
      ipch <- rep(pch[1], num.doc)   # By defect symbol 20
    } 
    coo <- labe <- coll <- fonte <- NULL
    
  #  coo <- labe <- coll <- ipch <- fonte <- NULL
    coo <- coord.doc
    if (label) 
      labe <- rownames(coord.doc)
    else labe <- rep("", num.doc)

        
    if(length(col.doc)==1)  coll <- rep(col.doc, num.doc) else{
     if(length(col.doc)< num.doc) col.doc <- c(col.doc, rep("grey", num.doc-length(col.doc)))
      if(length(col.doc)> num.doc) col.doc <- col.doc[1:num.doc] 
      coll <- col.doc
    }

    fonte <- rep(1, num.doc)

    if (!is.null(selection)) 
        coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% selection)])), 
                                                       alpha = 255 * (1 - unselect), maxColorValue = 255)

    pos.text <- NULL
    if (any(labe != "")) {
      autoLab <- FALSE
      if (auto.Lab == "auto") autoLab = (length(selection) < 50)
      if (auto.Lab == TRUE) autoLab <- TRUE
      if (autoLab == TRUE) {
        if(is.null(pch))    warning("autoLab must draw points")
        autoLab(coo[selection, 1], y = coo[selection,2],
                labels = labe[selection], col = coll[selection], 
                font = fonte[selection], cex=cex,...)
        points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
      }  # End of autolab ==TRUE
      else {  
        # It is not autolab
        if(!is.null(pch)) { 
          points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
          pos.text <- 3
        }
        text(coo[labe != "", 1], y = coo[labe != "",2],
             labels = labe[labe != ""], col = coll[labe != ""], 
             font = fonte[labe != ""], pos=pos.text , cex=cex,...)
      }
    } else {
      # Without labels
      if(is.null(pch)) stop("You must plot points and/or text")
      points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
    }
   }  # End plot documents
###################################### Document selection selDoc

    
    
##############################  Selection of words selWord
  if (!is.null(selWord)) {
    if (is.null(title)) 
      titre <- "Word factor map (LexGalt)"
    else titre <- title
    coord.word <- res.cagalt$word$coord[, axes, drop = FALSE]	
    num.word <- nrow(coord.word)
    
    xlim <- xlimT; ylim <- ylimT
    if (is.null(xlim)) 
      xlim <- c(min(coord.word[, 1]), max(coord.word[, 1])) * 1.2
    if (is.null(ylim)) 
      ylim <- c(min(coord.word[, 2]), max(coord.word[, 2])) * 1.2
    selection <- NULL
    rdo<-""
    
    # Selection of words
    if (mode(selWord) == "numeric") 
      selection <- selWord
    else if(length(selWord)==1) {
      if (selWord=="ALL")  selection <- c(1:num.word)
      else rdo <- selX(selWord)
    } else { selection <-  which(rownames(coord.word) %in% selWord)}
  
    
    if (rdo=="coord") 	
      selection <- (rev(order(apply(coord.word, 1, max))))[1:min(nrow(coord.word), 	
            sum(as.integer(unlist(strsplit(selWord, "coord"))), na.rm = TRUE))]	
    if (rdo=="cos2") {	
      selcos2 <- as.numeric(substr(selWord, 5, nchar(selWord)))	
      selection <- which(apply(res.cagalt$word$cos2[, axes], 1, sum)> selcos2) }	
    if (rdo=="contrib") {	
      selcontrib <- as.numeric(substr(selWord, 8, nchar(selWord)))	
      dft <- data.frame(res.cagalt$word$contrib[,c(axes[1],axes[2]),drop=FALSE])
      fval <- apply(dft, 1, function(z) max(abs(z)))
      selection <- which(fval>= selcontrib) }	
    if (rdo=="meta") {	
      selmeta <- as.numeric(substr(selWord, 5, nchar(selWord)))	
      sMeta <- which(res.cagalt$word$contrib[,1] > selmeta*mean(res.cagalt$word$contrib[,1]))
      sMeta <- c(sMeta, which(res.cagalt$word$contrib[,2] > selmeta*mean(res.cagalt$word$contrib[,2])))
      selection <- unique(sMeta) }	
    
    if(is.null(selection)) if(rdo=="") if(length(selWord)==1) selection <- which(rownames(coord.word) %in% selWord)
    if(length(selection)==0) stop("There are not selected words to plot")	
    
    if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")))   
      dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
    else dev.new()
    plot(0, 0, main = titre, xlab = lab.x, ylab = lab.y, 	
         xlim = xlim, ylim = ylim, col = "white", asp = 1, ...)	
    abline(v = 0, lty = 2, ...)	
    abline(h = 0, lty = 2, ...)	
    
    ipch <- 20
    if(!is.null(pch)) {                                     # symbols in graphs
      ipch <- rep(pch[1], nrow(coord.word))   # By defect symbol 20
    } 
    
    coo <- labe <- coll <- fonte <- NULL	
    coo <- coord.word	

    if (label) 	
      labe <- rownames(coord.word)	
    else labe <- rep("", nrow(coord.word))	


    if(length(col.word)==1)  coll <- rep(col.word, nrow(coord.word)) else{
      if(length(col.word)< nrow(coord.word)) col.word <- c(col.word, rep("darkred", nrow(coord.word)-length(col.word)))
      if(length(col.word)> nrow(coord.word)) col.word <- col.word[1:nrow(coord.word)] 
      coll <- col.word
    }

  #  coll <- rep(col.word, nrow(coord.word))	
  #  ipch <- rep(15, nrow(coord.word))	
    fonte <- rep(1, nrow(coord.word))	
    
    
    if (!is.null(selection)) 
      coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% selection)])), 
                                                     alpha = 255 * (1 - unselect), maxColorValue = 255)
    pos.text <- NULL
    
    if (any(labe != "")) {
      autoLab <- FALSE
      if (auto.Lab == "auto") autoLab = (length(selection) < 50)
      if (auto.Lab == TRUE) autoLab <- TRUE
      if (autoLab == TRUE) {
        if(is.null(pch))    warning("autoLab must draw points")
        autoLab(coo[selection, 1], y = coo[selection,2],
                labels = labe[selection], col = coll[selection], 
                font = fonte[selection], cex=cex,...)
        points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
      }  # Final de autolab ==TRUE
      else {  
        # NO es autolab
        if(!is.null(pch)) { 
          points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
          pos.text <- 3
        }
        text(coo[labe != "", 1], y = coo[labe != "",2],
             labels = labe[labe != ""], col = coll[labe != ""], 
             font = fonte[labe != ""], pos=pos.text , cex=cex,...)
      }
    } else {
      # Sin etiquetas
      if(is.null(pch)) stop("You must plot points and/or text")
      points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
    }
  }  # Final plot words

       
# Final selection of words    
#################################################3

    
    
    
  # Word ellipses selected
  if(!is.null(selWordEllip)){
    selectionEll <- NULL
    rdoEllip<-""
    if (mode(selWordEllip) == "numeric") 
      selWordEllip <- rownames(coord.word)[selWordEllip]
    else if(length(selWordEllip)==1) {
      if (selWordEllip=="ALL") selWordEllip <- rownames(coord.word)
      else rdoEllip <- selX(selWordEllip)
    }
    if (rdoEllip=="coord") 
      selectionEll  <- (rev(order(apply(coord.word[,axes]^2, 1, max))))[1:min(num.word, 
                  sum(as.integer(unlist(strsplit(selWordEllip, "coord"))), na.rm = TRUE))]
    if (rdoEllip=="cos2") {
      selcos2 <- as.numeric(substr(selWordEllip, 5, nchar(selWordEllip)))
      selectionEll  <- which(apply(res.cagalt$word$cos2[, axes], 1, sum)> selcos2) }
    if (rdoEllip=="contrib") {
      selcontrib <- as.numeric(substr(selWordEllip, 8, nchar(selWordEllip)))
      dft <- data.frame(res.cagalt$word$contrib[,c(axes[1],axes[2]),drop=FALSE])
      fval <- apply(dft, 1, function(z) max(abs(z)))
      selectionEll  <- which(fval>= selcontrib) }
    if (rdoEllip=="meta") {
      selmeta <- as.numeric(substr(selWordEllip, 5, nchar(selWordEllip)))
      sMeta <- which(apply(res.cagalt$word$contrib[, axes], 1, function(a) a[axes]> selmeta*mean(a[axes])))
      selectionEll  <- unique(sMeta) }
    if(!is.null(selectionEll)) selWordEllip <- rownames(coord.word)[selectionEll]
    selWordEllip <-selWordEllip[which(selWordEllip %in% rownames(coord.word)[selection])]
    if(length(selWordEllip)==0) stop("There are not selected words to plot ellipses")

            
    dfEll <- res.cagalt$ellip$word[,c(ncol(res.cagalt$ellip$word), axes),drop=FALSE]
    dfEll$FREQ <- as.factor(dfEll$FREQ)
    coord.ellip <- FactoMineR::coord.ellipse(dfEll,level.conf=level.conf,bary = FALSE)
    

    col.word <- coll[which(rownames(coord.word) %in% selWordEllip)]      
    for(i in 1:length(selWordEllip)){
      lines(coord.ellip$res[coord.ellip$res$FREQ==selWordEllip[i],2], 
            coord.ellip$res[coord.ellip$res$FREQ==selWordEllip[i],3], lty = 2, lwd = 2, col = col.word[i])
    }
}# Final Words and ellipses of words

    
    
  
# Selection of quantitative variables  
  if (!is.null(selQuantiVar)) {
    if (is.null(res.cagalt$quanti.var)) 
      stop("Variables are not quantitative")
    
}
# Final quantitative variables  
  

##############################################  Qualitative variables
  if (!is.null(selQualiVar)) { 
    if (is.null(res.cagalt$quali.var)) 
      stop("Variables are not categorical")
    if (is.null(title)) 
      titre <- "Categories factor map (CaGalt)"
    else titre <- title
    coord.var <- res.cagalt$quali.var$coord[, axes, drop = FALSE]
    num.var <-  nrow(coord.var)
    
    xlim <- xlimT; ylim <- ylimT
    
    if (is.null(xlim)) 
      xlim <- c(min(coord.var[, 1]), max(coord.var[, 1])) * 1.2
    if (is.null(ylim)) 
      ylim <- c(min(coord.var[, 2]), max(coord.var[, 2])) * 1.2
    
    selection <- NULL
    rdo<-""
    
    if (mode(selQualiVar) == "numeric") selection <- selQualiVar
    else if(length(selQualiVar)==1) {
      if (selQualiVar=="ALL")  selection <- c(1:num.var)
      else rdo <- selX(selQualiVar)
    } else { selection <-  which(rownames(coord.var) %in% selQualiVar)}
    
    if (rdo=="coord") 	
      selection <- (rev(order(apply(coord.var^2, 1, max))))[1:min(num.var, 	
           sum(as.integer(unlist(strsplit(selQualiVar, "coord"))), na.rm = TRUE))]	
    if (rdo=="cos2") {	
      selcos2 <- as.numeric(substr(selQualiVar, 5, nchar(selQualiVar)))	
      selection <- which(apply(res.cagalt$quali.var$cos2[, axes], 1, sum)> selcos2)	
    }	
    
    if (rdo=="contrib") stop("There is not possible to use contrib selection for selQualiVar")
    if (rdo=="meta") stop("There is not possible to use meta selection for selQualiVar")
    
    
    
    if(is.null(selection)) if(rdo=="") if(length(selQualiVar)==1) selection <- which(rownames(coord.var) %in% selQualiVar)
    if(length(selection)==0) stop("There are not selected elements to plot")
    
    if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")))   
      dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
    else dev.new()
    plot(0, 0, main = titre, xlab = lab.x, ylab = lab.y, 
         xlim = xlim, ylim = ylim, col = "white", asp = 1, ...)
    abline(v = 0, lty = 2, ...)
    abline(h = 0, lty = 2, ...)

    ipch <- 20      
    if(!is.null(pch)) {                                     # symbols in graphs
      ipch <- rep(pch[1], nrow(coord.var))   # By defect symbol 20
    } 
    
    
    coo <- labe <- coll <-  fonte <- NULL
    coo <- coord.var
    if (label) 
      labe <- rownames(coord.var)
    else labe <- rep("", nrow(coord.var))
    coll <- rep(col.quali, nrow(coord.var))

    if(length(col.quali)==1)  coll <- rep(col.quali, nrow(coord.var)) else{
      if(length(col.quali)< nrow(coord.var)) col.quali <- c(col.quali, rep(col, nrow(coord.var)-length(col.quali)))
      if(length(col.quali)> nrow(coord.var)) col.quali <- col.quali[1:nrow(coord.var)] 
      coll <- col.quali
    }
    
    fonte <- rep(1, nrow(coord.var))
    
    if (!is.null(selQualiVar)) {
        coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% 
           selection)])), alpha = 255 * (1 - unselect), maxColorValue = 255)
    }
   
    pos.text <- NULL
      if (any(labe != "")) {
      
      autoLab <- FALSE
      if (auto.Lab == "auto") autoLab = (length(selection) < 50)
      if (auto.Lab == TRUE) autoLab <- TRUE
      
      if (autoLab == TRUE) {
        if(is.null(pch))    warning("autoLab must draw points")
        autoLab(coo[selection, 1], y = coo[selection,2],
                labels = labe[selection], col = coll[selection], 
                font = fonte[selection], cex=cex,...)
        points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
      }  # Final de autolab ==TRUE
      else {  
        # NO es autolab
        if(!is.null(pch)) { 
          points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
          pos.text <- 3
        }
        text(coo[labe != "", 1], y = coo[labe != "",2],
             labels = labe[labe != ""], col = coll[labe != ""], 
             font = fonte[labe != ""], pos=pos.text , cex=cex,...)
      }
     
      } else {
        # Sin etiquetas
        if(is.null(pch)) stop("You must plot points and/or text")
        points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
      }

   
    
    
#########################################  Final qualitative variables


  # Qualitative categories ellipses
    if(!is.null(selQualiVarEllip)){
      selectionEll <- NULL
      rdoEllip<-""
      if(selQualiVarEllip[1]=="ALL")  selQualiVarEllip <- "ALL"

      if (mode(selQualiVarEllip) == "numeric") 
        selQualiVarEllip <- rownames(coord.var)[selQualiVarEllip]
      else if(length(selQualiVarEllip)==1) {
        if (selQualiVarEllip=="ALL") selQualiVarEllip <- rownames(coord.var)
        else rdoEllip <- selX(selQualiVarEllip)
      }# 
      
      if (rdoEllip=="coord") 	
        selectionEll  <- (rev(order(apply(coord.var[,axes]^2, 1, max))))[1:min(num.var, 	
                        sum(as.integer(unlist(strsplit(selQualiVarEllip, "coord"))), na.rm = TRUE))]	
      if (rdoEllip=="cos2") {	
        selcos2 <- as.numeric(substr(selQualiVarEllip, 5, nchar(selQualiVarEllip)))	
        selectionEll  <- which(apply(res.cagalt$quali.var$cos2[, axes], 1, sum)> selcos2) }	
      if (rdoEllip=="contrib") {	
        selcontrib <- as.numeric(substr(selQualiVarEllip, 8, nchar(selQualiVarEllip)))	
        dft <- data.frame(res.cagalt$quali.var$contrib[,c(axes[1],axes[2]),drop=FALSE])	
        fval <- apply(dft, 1, function(z) max(abs(z)))	
        selectionEll  <- which(fval>= selcontrib) }	
      if (rdoEllip=="meta") {
        selmeta <- as.numeric(substr(selQualiVarEllip, 5, nchar(selQualiVarEllip)))
        sMeta <- which(apply(res.cagalt$quali.var$contrib[, axes], 1, function(a) a[axes]> selmeta*mean(a[axes])))
        selectionEll  <- unique(sMeta) }
      
      if(!is.null(selectionEll)) selQualiVarEllip <- rownames(res.cagalt$quali.var$coord)[selectionEll]
      selQualiVarEllip <- selQualiVarEllip[which(selQualiVarEllip %in% rownames(res.cagalt$quali.var$coord)[selection])]
      if(length(selQualiVarEllip)==0) stop("There are not selected elements to plot ellipses")

      # Coordinates of all simulations. Last column has the word
      dfEll <- res.cagalt$ellip$var[,c(ncol(res.cagalt$ellip$var), axes),drop=FALSE]
      dfEll$VAR <- as.factor(dfEll$VAR)
      coord.ellip <- FactoMineR::coord.ellipse(dfEll,level.conf=level.conf,bary = FALSE)


col.quali <- coll[which(rownames(res.cagalt$quali.var$coord) %in% selQualiVarEllip)]      
      for(i in 1:length(selQualiVarEllip )){
        lines(coord.ellip$res[coord.ellip$res$VAR==selQualiVarEllip[i],2], 
              coord.ellip$res[coord.ellip$res$VAR==selQualiVarEllip[i],3], lty = 2, lwd = 2, col = col.quali[i])
      }
    } # Final ellipses for qualitative variables
  } # Final qualitative variables    


    
    
    
          
########################  Quantitative variables
  if (!is.null(selQuantiVar)) {
    if (is.null(res.cagalt$quanti.var)) 
      stop("Variables are not quantitative")
    if (is.null(title)) 
      titre <- "Variables factor map (CaGalt)"
    else titre <- title
 
    # Review if selCov is necessary   
    if(selCov==TRUE) {
      coord.Q <- res.cagalt$quanti.var$coord[, axes, drop = FALSE]
      xlim <- xlimT; ylim <- ylimT
      if (is.null(xlim)) 
        xlim <- c(min(coord.Q[, 1]), max(coord.Q[, 1])) * 1.2
      if (is.null(ylim)) 
        ylim <- c(min(coord.Q[, 2]), max(coord.Q[, 2])) * 1.2
    } else {
      coord.Q <- res.cagalt$quanti.var$cor[, axes, drop = FALSE]
      xlim <- ylim <- c(-1, 1)
    }
    selection <- NULL
    rdo<-""
    
    if (mode(selQuantiVar) == "numeric") 
      selection <- selQuantiVar
    else if(length(selQuantiVar)==1) {
      if (selQuantiVar=="ALL")  selection <- c(1:nrow(coord.Q))
      else rdo <- selX(selQuantiVar)
    } else { selection <-  which(rownames(coord.Q) %in% selQuantiVar)}
    
    if (rdo=="coord") 
      selection <- (rev(order(apply(coord.Q^2, 1, max))))[1:min(nrow(coord.Q), 
                                                                       sum(as.integer(unlist(strsplit(selQuantiVar, "coord"))), na.rm = T))]
    if (rdo=="cos2") {
      selcos2 <- as.numeric(substr(selQuantiVar, 5, nchar(selQuantiVar)))
      selection <- which(apply(res.cagalt$quanti.var$cos2[, axes], 1, sum)> selcos2)
    }
    if (rdo=="contrib") stop("There is not possible to use contrib selection for selQuantiVar")
    
    
    if(is.null(selection)) if(rdo=="") if(length(selQuantiVar)==1) selection <- which(rownames(coord.Q) %in% selQuantiVar)
    if(length(selection)==0) stop("There are not selected elements to plot")
 
#    if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")))   
#      dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
#    else if(new.plot==TRUE) dev.new()

    if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
      dev.new()

   
        plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, 
         col = "white", asp = 1, main = titre, ...)
    if(selCov==FALSE){
      x.cercle <- seq(-1, 1, by = 0.01)
      y.cercle <- sqrt(1 - x.cercle^2)
      lines(x.cercle, y = y.cercle, ...)
      lines(x.cercle, y = -y.cercle, ...)}
    
    abline(v = 0, lty = 2, ...)
    abline(h = 0, lty = 2, ...)
    coll <- coo <- labe <- posi <- NULL
    
    coll <- rep(col.quanti, nrow(coord.Q))
    coo <- coord.Q                                 # Coordinates quantitative variables without simulation

    
    if (label) 	
      labe <- rownames(coord.Q)	
    else labe <- rep("", nrow(coord.Q))	
    
    if (!is.null(selection)) {	
      if (is.numeric(unselect)) 	
        coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% 	
                                                                          selection)])), alpha = 255 * (1 - unselect),maxColorValue = 255)	
      else coll[!((1:length(coll)) %in% selection)] = unselect	
      labe[!((1:length(coll)) %in% selection)] <- ""	
    }	
    
    for (v in 1:nrow(coord.Q)) {	
      arrows(0, 0, coord.Q[v, 1], coord.Q[v, 2], 	
             length = 0.1, angle = 15, code = 2, col = coll[v])	
      if (label) {	
        if (abs(coord.Q[v, 1]) > abs(coord.Q[v,2])) {	
          if (coord.Q[v, 1] >= 0) 	
            posi <- c(posi, 4)	
          else posi <- c(posi, 2)	
        }	
        else {	
          if (coord.Q[v, 2] >= 0) 	
            posi <- c(posi, 3)	
          else posi <- c(posi, 1)	
        }	
      }	
    }
    

    if (any(labe != "")) {	
      autoLab <- FALSE	
      if (auto.Lab == "auto")  autoLab = (length(which(labe != "")) < 50)	
      if (auto.Lab == TRUE) autoLab <- TRUE
      if (autoLab == TRUE){ 	
        autoLab(coo[labe != "", 1], y = coo[labe != "", 2], labels = labe[labe != ""],	
      #          col = coll[labe != ""], shadotext = shadowtext,cex=cex, ...)	}
      col = coll[labe != ""], cex=cex, ...)	}
      else	
        text(coo[labe != "", 1], y = coo[labe != "", 2], labels = labe[labe != ""], 	
             pos = posi[labe != ""], col = coll[labe != ""], cex=cex,...)	
    }	
    
    # Ellipses selection for Quantitative analysis
    if(!is.null(selQuantiVarEllip)){	
      selectionEll <- NULL	
      rdoEllip<-""	
      if (mode(selQuantiVarEllip) == "numeric") 	
        selQuantiVarEllip <- rownames(res.cagalt$quanti.var$coord)[selQuantiVarEllip]	
      else if(length(selQuantiVarEllip)==1) {	
        if (selQuantiVarEllip=="ALL") selQuantiVarEllip <- rownames(res.cagalt$quanti.var$coord)	
      }# Final no 1	
      

      selQuantiVarEllip <- selQuantiVarEllip[which(selQuantiVarEllip %in% rownames(res.cagalt$quanti.var$coord)[selection])]
      if(length(selQuantiVarEllip)==0) stop("There are not selected elements to plot ellipses")
      

   
      # Contiene las coordenadas de todas las simulaciones para los ejes seleccionados,
      # La primera columna contiene la palabra
      # Contiene las correlaciones en todos los casos
      dfEll <- res.cagalt$ellip$var[,c(ncol(res.cagalt$ellip$var), axes),drop=FALSE]
      dfEll$VAR <- as.factor(dfEll$VAR)

      if(selCov==TRUE)
      {
        factorBoot  <- data.frame(res.cagalt$quanti.var$coord[,1]/res.cagalt$quanti.var$cor[,1],
                                  rownames(res.cagalt$quanti.var$cor))
        colnames(factorBoot) <- c("Value","VAR")
        tmp <- merge(x =dfEll, y = factorBoot, by = "VAR", all.x = TRUE)
        tmp[,2:3] <- tmp[,2:3]*tmp[,"Value"]
        dfEll <- tmp[,1:3] 
      }
      
   
      coord.ellip <- coord.ellipse(dfEll,level.conf=level.conf,bary = FALSE)

      if(selCov==FALSE) {
        tmp <- sqrt(coord.ellip$res[,2]^2+coord.ellip$res[,3]^2)
        coord.ellip$res[which(tmp>1),2] <- coord.ellip$res[which(tmp>1),2]/tmp[which(tmp>1)]
        coord.ellip$res[which(tmp>1),3] <- coord.ellip$res[which(tmp>1),3]/tmp[which(tmp>1)]
      }
      
      


      for(i in 1:length(selQuantiVarEllip )){
        lines(coord.ellip$res[coord.ellip$res$VAR==selQuantiVarEllip[i],2], 
              coord.ellip$res[coord.ellip$res$VAR==selQuantiVarEllip[i],3], 
              lty = 2, lwd = 2, col = col.quanti)
      }
    } # Final selQuantiVarEllip
    
  } # Final Quantitative variables    

  } else { 
# This is the multiple case
    
    
    
    
    if(type=="QL") {		
      if(is.null(x$MQL) | length(x$MQL)==1) stop("There are not categorical variables selected or in LexGalt object")		
      res.mfagalt <- x$MQL }		
    if(type=="QN") {			
      if(is.null(x$MQN) | length(x$MQN)==1)stop("There are not quantitative variables selected or in LexGalt object")		
      res.mfagalt <- x$MQN }

    # Selection of groups
    ## Problema. Tenemos en res.mfagalt$call el nombre de los grupos name.groups y el número num.groups
    ## Podemos tener seleccionados grupos en selGroup
    ## Podemos tener nuevas etiquetas en label.group
    name.groups <- res.mfagalt$call$name.groups
    num.groups <-res.mfagalt$call$num.groups
    if(is.null(selGroup)) selGroup <- c(1:num.groups)	
    if(length(selGroup)==1) if(selGroup=="ALL") selGroup <- c(1:num.groups)
    if (is.numeric(selGroup)) selGroup <- name.groups[selGroup]
    selGroup <- selGroup[which(selGroup %in%  name.groups)]
    if(length(selGroup)==0) stop("There is not selected groups")
    noselGroup <- name.groups[-which(name.groups %in% selGroup)]		
  
    # Labels of groups
    # Tienen que ser el mismo número que name.groups, estén seleccionadas o no

    
   old.name.groups <- name.groups
        if(is.null(label.group)) label.group <- name.groups
    if(length(label.group)==num.groups) name.groups <- label.group 
    else if(length(label.group)==length(selGroup)) name.groups[which(name.groups %in% selGroup)] <- label.group
         else      warning("Number of label.group is not the same as number of groups or selected selGroup) \n     label.group will be ignorated")
     label.group <- name.groups 
  
    
      
    # Color selection		
    if(is.null(col.group)) col.group <- PALETTE[1:res.mfagalt$call$num.groups]	
    if (is.numeric(unselect)) if ((unselect > 1) | (unselect < 0)) stop("unselect should be between 0 and 1")		
    # Labels 
    lab.x <- paste("Dim ", axes[1], " (", format(res.mfagalt$eig[axes[1], 2], nsmall = 2, digits = 2), "%)", sep = "")		
    lab.y <- paste("Dim ", axes[2], " (", format(res.mfagalt$eig[axes[2], 2], nsmall = 2, digits = 2), "%)", sep = "")		

    ###################################### Plot eigenvalues   ########   Same as simple
    if(eigen) {
      #  nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")) # TRUE  when inside RStudio ,  # FALSE when outside RStudio
      
      if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")))   
        dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
      else if(new.plot==TRUE) dev.new()
      
      if(is.null(title)) titleE <- "Eigenvalues" else titleE <- title
      barplot(res.mfagalt$eig[, 1], main = titleE, col=col.eig, cex.axis=cex, cex.names=cex,
              names.arg = paste("dim",1:nrow(res.mfagalt$eig)))}
    ###################################### End Plot eigenvalues 


    
    # selection <- NULL  # Number of non selected groups
    if(length(noselGroup)>0)
      selection <- which(!rownames(res.mfagalt$group$coord[, axes, drop = FALSE]) %in% noselGroup)
        else selection <- c(1:length(selGroup))

    
    
#######################################################################
#### Gráfico en donde aparecen como modalidades los grupos
# If plot groups					

    
if(plot.group==TRUE){	
  coord.actif <- res.mfagalt$group$coord[, axes, drop = FALSE]	
  xlim <- xlimT; ylim <- ylimT	
  if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")))   
    dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
  else dev.new()
  
  coo <- labe <- coll <- ipch <- fonte <- NULL			
  if (is.null(xlim)) xlim <- c(0,1)			
  if (is.null(ylim)) ylim <- c(0,1)			
  if (is.null(title)) title <- "Groups representation (MfaGalt)"			
  plot(0, 0, main = title, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, asp = 1, col = "white", ...)			
  abline(v = 0, lty = 2, ...)			
  abline(h = 0, lty = 2, ...)			
  # coo <- rbind(coo, coord.actif)	
  coo <- coord.actif
  if(length(col.group)==1) col.group[1:num.groups] <- col.group
  coll <- col.group

 # label.group
    if(length(noselGroup)>0) 
       coll[!((1:length(coll)) %in% selection)] <- rgb(t(col2rgb(coll[!((1:length(coll)) %in% selection)])), 
                                               alpha = 255 * (1 - unselect), maxColorValue = 255)


  # Posibilidad de que el número de grupos no sea igual que el de símbolos
  
  if(is.null(pch)) ipch <- c(15:(15+num.groups-1))
  if(length(pch)==1)  ipch <- c(pch:(pch+num.groups-1))	else ipch <- pch	
  if(length(pch)>length(ipch)) maxgr <- pch else maxgr <- ipch
     pch <- ipch <- maxgr

  fonte <- rep(1, nrow(coord.actif))

  

  
    if(label) labe <- label.group
  pos.text <- NULL
    if (any(labe != "")) {
    autoLab <- FALSE
    if (auto.Lab == "auto") autoLab = (length(selection) < 50)
    if (auto.Lab == TRUE) autoLab <- TRUE
    if (autoLab == TRUE) {
      if(is.null(pch))    warning("autoLab must draw points")
      autoLab(coo[, 1], y = coo[,2],
              labels = labe,  col=coll, 
              font = fonte, cex=cex,...)
      points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
    }  # Final de autolab ==TRUE
    else {  
      # NO es autolab
      if(!is.null(pch)) { 
        points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
        pos.text <- 3
      }
      # text(coo[labe != "", 1], y = coo[labe != "",2],
      #     labels = labe[labe != ""], col = coll[labe != ""], 
      #     font = fonte[labe != ""], pos=pos.text , cex=cex,...)
      text(coo[, 1], y = coo[,2], labels = labe, col = coll, 
               font = fonte, pos=pos.text , cex=cex,...)
           
    }
  } else {
    # Sin etiquetas
    if(is.null(pch)) stop("You must plot points and/or text")
    points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
  }
} # Final if plot.group
#### Finalizado el plot.group


    
    
    
    
    
    
if(!is.null(selWord)){		
  #############  Plot selected words		
  titre <- ifelse(is.null(title), "Word factor map (M LexGalt)", title)		
  coord.word <- res.mfagalt$word$coord[, axes, drop = FALSE]	  # Cambiar después para simplificar
  num.word <- nrow(coord.word)
  selection.word <- NULL                                                    #  <-------------- PDte 

    rdo <-""                                                              # Para la selección cos2...
  sel <- nosel <- NULL
  xlim <- xlimT; ylim <- ylimT		
  if (is.null(xlim)) xlim <- c(min(coord.word[,1]),max(coord.word[,1])) * 1.2		
  if (is.null(ylim)) ylim <- c(min(coord.word[,2]),max(coord.word[,2])) * 1.2		

#  if(length(col.word)==1)  coll <- rep(col.word, nrow(coord.word)) else{
#    if(length(col.word)< nrow(coord.word)) coll <- c(col.word, rep(col.word, nrow(coord.word)-length(col.word)))
#    if(length(col.word)> nrow(coord.word)) coll <- col.word[1:nrow(coord.word)] 
#    #  coll <- col.word
#  }
  
    # First of all we must to remove not selected words by groups
  # Selecting rows (words)	from groups	
  pos <- matrix(, nrow = num.groups, ncol = 3)		
  colnames(pos) <- c("group","start","end")		
  pos[,1] <- name.groups	
  ncont <- 0			
  for(i in 1:num.groups) {			
    pos[i,2] <- ncont+1			
    ncont <- ncont+res.mfagalt$call$num.freq[i]			
    pos[i,3] <- ncont    }			
  # pos has the position of each group
  # group     start end  
  # [1,] "GROUP.1" "1"   "169"
  # [2,] "GROUP.2" "170" "316"
  numnoselgroup <- length(noselGroup)	
    if(numnoselgroup>0) {
     for(i in 1:num.groups) {						
      tmp <- c(pos[i,2]:pos[i,3])						
      tmp2 <- which(name.groups[i] %in% noselGroup)
      if(length(tmp2)==0) sel <- c(sel,tmp) # else nosel <- c(nosel,tmp)	
    }						
  } else { sel <- c(1:num.word)}	 # Final numnoselgroup>0
  # sel contiene las posiciones de los individuos seleccionados por grupos, nosel los no seleccionados
nosel <- which(!c(1:num.word) %in% sel)
noselGROUP <- nosel

if(is.null(pch)) ipch <- c(15:(15+num.groups-1))
if(length(pch)==1)  ipch <- c(pch:(pch+num.groups-1))	else ipch <- pch	


   ########################### Selection of words
 if(!is.list(selWord)) {
     rdo <- selX(selWord)
    if(length(rdo)==1){
    if (rdo=="cos2") {						
      selcos2 <- as.numeric(substr(selWord, 5, nchar(selWord)))	
      # Calcula las palabras que no tienen una suma de cos2 entre los dos factores mayor que la seleccionada
      noselcos2 <- which(apply(res.mfagalt$word$cos2[, axes], 1, sum) < selcos2)
      nosel <- unique(c(nosel,noselcos2)) }		
    
     if (rdo=="coord") {						
       selcoord <- as.numeric(substr(selWord, 6, nchar(selWord)))	
       C <- res.mfagalt$word$coord						
       rownames(C) <- c(1:nrow(C))		
       C <- data.frame(C[,axes], c(1:nrow(C)), rownames(coord.word))	# Añade a la última columna la palabra
       C2 <- data.frame(apply(C^2, 1, max), c(1:nrow(C)), rownames(coord.word))[sel,]						
       selcoord <- rownames(C2[ order(-C2[,1]), ]  [c(1:selcoord),])		
       temp <- 1:nrow(C)						
       noselcoord <- temp[!temp %in% selcoord]						
       nosel <- unique(c(nosel,noselcoord)) }
     
     if (rdo=="contrib") {						
       selcontrib <- as.numeric(substr(selWord, 8, nchar(selWord)))						
       C <- res.mfagalt$word$contrib						
       rownames(C) <- c(1:nrow(C))		
       C <- data.frame(C[,axes], c(1:nrow(C)), rownames(res.mfagalt$word$contrib))	
       
       C2 <- data.frame(apply(C[,c(1,2)], 1, max), c(1:nrow(C)), rownames(res.mfagalt$word$coord))[sel,]
       noselcontrib <- which(C2[,1]< selcontrib)						
       nosel <- unique(c(nosel,noselcontrib))
       }   
     
     if (rdo=="meta") {						
       selmeta <- as.numeric(substr(selWord, 5, nchar(selWord)))						
       C <- res.mfagalt$word$contrib						
       rownames(C) <- c(1:nrow(C))						
       C <- data.frame(C[,axes], c(1:nrow(C)), rownames(res.mfagalt$word$contrib))						
       sMeta <- which(C[,1] > selmeta*mean(C[,1]))						
       sMeta <- c(sMeta, which(C[,2] > selmeta*mean(C[,2])))						
       #       sMeta <- which(apply(C[, axes], 1, function(a) a[axes]> selmeta*mean(a[axes])))						
       sMeta <- unique(sMeta)   						
       temp <- 1:nrow(C)						
       nosel <- unique(c(nosel,temp[!temp %in% sMeta]))}   

      sel <- which(!c(1:nrow(res.mfagalt$word$coord)) %in% nosel)
     } # Final de si es solo un valor

     # nosel contiene los no seleccionados por palabras


     if(selWord[1]!="ALL"){
      if(rdo[1]=="") {

       if(is.numeric(selWord)) selWord <- rownames(coord.word)[selWord] 
       nosel2 <-  which(!rownames(coord.word) %in% selWord)
       nosel <- unique(c(nosel2,nosel))
       sel <- which(!c(1:num.word) %in% nosel)
     }} # Final rdo==""
    #  stop("NO Es una lista")
  #stop("Final de no es una lista")
  } else {
    # is.list(selWord)
    
    numselgroup <- num.groups-numnoselgroup
    
    if(numselgroup != length(selWord)) 
      stop(paste0("Error, number of selected groups is " , numselgroup, 
                  "\n  but the number of groups in selWord is ", length(selWord) ))
    
      sel <- NULL
    cont <- 1
    
    
    for(i in 1:num.groups) {						
      if(name.groups[i] %in% selGroup) {
        strtmp <- rownames(coord.word)[pos[i,2]:pos[i,3]]  # Palabras totales del grupo i
        sel <-  c(sel, (which(strtmp %in%  selWord[[cont]])-1+as.numeric(pos[i,2])))  # MOdificar quitar primer sel
        cont <- cont+1
      }
      #            
    } 
      #         stop("Es una lista")
    nosel <- which(!c(1:num.word) %in% sel)
  }
  
  ############## Graficos
# 1.- Poner el color de la palabra dependiendo del grupo al que pertenece
  if(length(sel)==0) stop("There are no words to plot")
  coo <- labe <- coll <- fonte <- NULL			
  # coord.word is a matrix, with words in rows and dimensions in columns
  coo <- coord.word
  

  # label.group
  if(length(col.group)!= num.groups) {
    warning("col.group is bad defined, it will be ignorated")
    col.group <- PALETTE[1:res.mfagalt$call$num.groups]
  }

   ipch <- NULL
  # Posibilidad de que el número de grupos no sea igual que el de símbolos
  if(is.null(pch)) 
    ipch <- rep(26,num.word) else {
    pch <- rep(pch,num.groups)
  if(length(pch)==1)  {  ipch <- rep(as.numeric(pch),num.word) 
  }else {
  if(length(pch)>num.groups) pch <- pch[1:num.groups] 
  if(length(pch)<num.groups) stop("length pch is < than num.groups")
  if(is.null(ipch))
  for(i in 1:num.groups) {						
      itmp <- as.numeric(pos[i,3]) - as.numeric(pos[i,2]) +1 
      ipch <- c(ipch, rep(pch[i], itmp))
  } 
  }}

   
      coll <- rep(col.group,res.mfagalt$call$num.freq)	 
  if(!is.null(col.word)) {
    if(length(col.word)!=length(coll)) stop(paste0("length.col must have ",length(coll), " values, not ", length(col.word) ))
   coll <- col.word 
  }
 coll[c(1:num.word) %in% nosel] <- rgb(t(col2rgb(coll[c(1:num.word) %in% nosel])), 
              alpha = 255 * (1 - unselect), maxColorValue = 255)	
 coll[((1:num.word) %in% noselGROUP)] <- rgb(t(col2rgb(coll[((1:num.word) %in% noselGROUP)])), 
                                         alpha = 0, maxColorValue = 255)

 if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")))   
   dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
 else dev.new()
 
  plot(0, 0, main = titre, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, col = "white", asp = 1, ...)		
  abline(v = 0, lty = 2, ...)		
  abline(h = 0, lty = 2, ...)		


  if(label)
    labe <-  rownames(coord.word)
   else labe <- rep("",nrow(coord.word))

  fonte <- rep(1, nrow(coord.word))	
#  ipch <- 20
#  if(!is.null(pch)) {                                     # symbols in graphs
#    ipch <- rep(pch[1], nrow(coord.word))   # By defect symbol 20
#  } 
  #if(length(nosel)!=0)
  #  text(coo[nosel, 1], y = coo[nosel, 2], labels = rownames(coord.word)[nosel], cex=cex, col=coll[nosel], ...)	 
  

  if (any(labe != "")) {
    autoLab <- FALSE
    if (auto.Lab == "auto") autoLab = (length(sel) < 50)
    if (auto.Lab == TRUE) autoLab <- TRUE
    if (autoLab == TRUE) {
      if(is.null(pch))    warning("autoLab must draw points")
      autoLab(coo[sel, 1], y = coo[sel,2],
              labels = labe[sel], col = coll[sel], 
              font = fonte[sel], cex=cex,...)
      points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
    }  # Final de autolab ==TRUE
    else {  
      # NO es autolab
      if(!is.null(pch)) { 
        points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
        pos.text <- 3
      } else pos.text <- NULL
            text(coo[labe != "", 1], y = coo[labe != "",2],
           labels = labe[labe != ""], col = coll[labe != ""], 
           font = fonte[labe != ""], pos=pos.text , cex=cex,...)
    }
  } else {
    # Sin etiquetas
    if(is.null(pch)) stop("You must plot points and/or text")
    points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
  }
  
#  text(coo[sel, 1], y = coo[sel, 2], labels = labe, cex=cex, col=coll[sel], ...)	
  
# Controlar antes los pch si solo es uno o hay varios

  

pch2 <- pch
if(is.null(pch2)) pch2 <- rep(NA, num.groups)

   if(legend) if(length(label.group)==num.groups){						
#    legend(pos.legend,label.group,pch=15:(15+res.mfagalt$call$num.group-1),cex=cex*1.5,text.col=col.word,col=col.word)				
  #  legend(pos.legend,label.group,pch=pch2,cex=cex*1.1,text.col=col.group,col=col.group)	
   
     legend(pos.legend,label.group[which(old.name.groups %in% selGroup)],pch=pch2,cex=cex*1.1,text.col=col.group,col=col.group)			

    } else { warning("The number of groups of name.groups is ",  res.mfagalt$call$num.groups, " not " , length(label.group))}						
  } # Final if(!is.null(selWord)){		

    
###########################################################
    if(!is.null(selQualiVar)){				
      if (is.null(title)) titre <- "Categories factor map (MfaGalt)"					
        else titre <- title					
      coord.var <- res.mfagalt$quali.var$coord[, axes]					
      xlim <- xlimT; ylim <- ylimT						
      # Si tiene parcial o no, calculo de los limites		
      # col.var.partial <- col.group	
      ncateg <- nrow(coord.var)
      col.var.partial <- rep(col.quali[1], ncateg)
      
 
     # ipch <- 20
     if(!is.null(pch)) {                                     # symbols in graphs
         ipch <- rep(pch[1],ncateg)   # By defect symbol 20
      } else ipch <- rep(26,ncateg) 


      ############ Seleccion de variables y/o modalidades					
      selection <- NULL						
      if (is.numeric(selQualiVar)) selQualiVar <- rownames(coord.var)[selQualiVar]						
      if(selQualiVar[1]=="ALL") selQualiVar <- rownames(coord.var)						
      selection  <- which(rownames(coord.var) %in% selQualiVar)					
      
      if(length(selQualiVar)==1){						
        rdo <- selX(selQualiVar)						
        if (rdo=="cos2") {						
          selcos2 <- as.numeric(substr(selQualiVar, 5, nchar(selQualiVar)))						
          selection <- which(apply(res.mfagalt$quali.var$cos2[, axes], 1, sum) >= selcos2) }		
        if (rdo=="coord") {						
          selcoord <- as.numeric(substr(selQualiVar, 6, nchar(selQualiVar)))						
          selection <- rev(order(apply(res.mfagalt$quali.var$coord[,axes]^2, 1, max)))[1:selcoord] }	
      }  # Final selQualiVar
      
  #    return(selection)
      
      if(partial==TRUE | partial=="ALL"){				
        strcateg <- rownames(coord.var)
        for(i in 1:res.mfagalt$call$num.groups){	
          new.coord <- res.mfagalt$quali.var$coord.partial[[i]][, axes, drop = FALSE]
          rownames(new.coord) <- paste0(name.groups[i],".", strcateg)
          coord.var <- rbind(coord.var,new.coord)
          if(is.null(col.group[i])) tmp <- rep(col.quali[1], ncateg) else tmp<- rep(col.group[i], ncateg)
          col.var.partial <- c( col.var.partial, tmp)
          if(length(pch) ==i) 
              {pch <- c(pch,pch[i])
            tmp <- rep(pch[1], ncateg)} else tmp<- rep(pch[i+1], ncateg)
                    ipch <- c(ipch, tmp)
                            }
       } # Final partial TRUE
      if (is.null(xlim)) xlim <- c(min(coord.var[,1]),max(coord.var[,1])) * 1.2						
      if (is.null(ylim)) ylim <- c(min(coord.var[,2]),max(coord.var[,2])) * 1.2

      coll <- col.var.partial
      

      #       if (is.null(xlim)) xlim <- c(min(min(coord.var[,1]),xPmin),max(max(coord.var[,1]),xPmax)) * 1.2						
      #      if (is.null(ylim)) ylim <- c(min(min(coord.var[,2]),yPmin),max(max(coord.var[,2]),yPmax)) * 1.2						
      #      } else {						
      #        if (is.null(xlim)) xlim <- c(min(coord.var[,1]),max(coord.var[,1])) * 1.2						
      #        if (is.null(ylim)) ylim <- c(min(coord.var[,2]),max(coord.var[,2])) * 1.2			
      
      if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
        dev.new(width = min(14, max(8, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1]))), height = 8)	
      plot(0, 0, main = title, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, col = "white", asp = 1, ...)						
      abline(v = 0, lty = 2, ...)					
      abline(h = 0, lty = 2, ...)					
    
      coo <- labe <- fonte <- NULL	
      labe <- trimws(rownames(coord.var))
      coo <-  coord.var
   
      # ipch <- rep(20,nrow(coord.var))
 
  #    if(!is.null(col.word)) {
  #      if(length(col.word)!=length(coll)) stop(paste0("length.col must have ",length(coll), " values, not ", length(col.word) ))
  #      coll <- col.word 
  #    }
      # selection
      
      nosel <- which(!c(1:ncateg) %in% selection)
      coll[c(1:ncateg) %in% nosel] <- rgb(t(col2rgb(coll[c(1:ncateg) %in% nosel])), 
                                            alpha = 255 * (1 - unselect), maxColorValue = 255)	
      
      if(partial==TRUE | partial=="ALL"){		
        for(i in 1:res.mfagalt$call$num.groups){	
          if(name.groups[i] %in% noselGroup) 
            coll[((i*ncateg) + c(1:ncateg))] <- rgb(t(col2rgb(coll[((i*ncateg) + c(1:ncateg))])), 
                                                alpha = 255 * (1 - unselect), maxColorValue = 255)	
        }
      }
    
      sel <- selection
      if (any(labe != "")) {
        autoLab <- FALSE
        if (auto.Lab == "auto") autoLab = (length(sel) < 100)
        if (auto.Lab == TRUE) autoLab <- TRUE
        if (autoLab == TRUE) {
          if(is.null(pch))    warning("autoLab must draw points")
          autoLab(coo[, 1], y = coo[,2],
                  labels = labe, col = coll, 
                  font = fonte, cex=cex,...)
          points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
        }  # Final de autolab ==TRUE
        else {  
          # NO es autolab
          if(!is.null(pch)) { 
            points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
            pos.text <- 3
          } else pos.text <- NULL
          text(coo[labe != "", 1], y = coo[labe != "",2],
               labels = labe[labe != ""], col = coll[labe != ""], 
               font = fonte[labe != ""], pos=pos.text , cex=cex,...)
        }
      } else {
        # Sin etiquetas
        if(is.null(pch)) stop("You must plot points and/or text")
        points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, cex=cex,...)
      }

      
      pch2 <- pch
      if(is.null(pch2)) pch2 <- rep(NA, num.groups)
      if(length(label.group)==num.groups){						
        #    legend(pos.legend,label.group,pch=15:(15+res.mfagalt$call$num.group-1),cex=cex*1.5,text.col=col.word,col=col.word)				
        #  legend(pos.legend,label.group,pch=pch2,cex=cex*1.1,text.col=col.group,col=col.group)				
        legend(pos.legend,label.group,pch=pch2,cex=cex*1.1,text.col=col.group,col=col.group)				
      } else { warning("The number of groups of name.groups is ",  res.mfagalt$call$num.groups, " not " , length(label.group))}						
    } # Final is.null(selQualiVar)						
  
 
    
    
 
    
    

    #####################################################################3			
    if(!is.null(selQuantiVar)){			
      if (is.null(res.mfagalt$quanti.var)) 			
        stop("Variables are not quantitative")			
      if (is.null(title)) 			
        titre <- "Variables factor map (LexGalt)"			
      else titre <- title			
      selection <- NULL			
      rdo<-""			
      
      

      xlim <- xlimT; ylim <- ylimT		
      if(selCov==TRUE) {		
        coord.Q <- res.mfagalt$quanti.var$coord[, axes, drop = FALSE]		
        if (is.null(xlim)) 		
          xlim <- c(min(coord.Q[, 1]), max(coord.Q[, 1])) * 1.2		
        if (is.null(ylim)) 		
          ylim <- c(min(coord.Q[, 2]), max(coord.Q[, 2])) * 1.2		
      } else {		
        coord.Q <- res.mfagalt$quanti.var$cor[, axes, drop = FALSE]		
        xlim <- ylim <- c(-1, 1)		
      } # Final selCov		
      

      
      if (mode(selQuantiVar) == "numeric") 			
        selection <- selQuantiVar			
      else if(length(selQuantiVar)==1) {			
        if (selQuantiVar=="ALL")  selection <- c(1:nrow(coord.Q))			
        else rdo <- selX(selQuantiVar)			
      } else { selection <-  which(rownames(coord.Q) %in% selQuantiVar)}			


      if (rdo=="coord") 			
        selection <- (rev(order(apply(coord.Q^2, 1, max))))[1:min(nrow(coord.Q), 			
                                                                         sum(as.integer(unlist(strsplit(selQuantiVar, "coord"))), na.rm = T))]			
      if (rdo=="cos2") {			
        selcos2 <- as.numeric(substr(selQuantiVar, 5, nchar(selQuantiVar)))			
        selection <- which(apply(res.mfagalt$quanti.var$cos2[, axes], 1, sum)> selcos2)			
      }			
      if (rdo=="contrib") stop("There is not possible to use contrib selection for selQuantiVar")			
      
      if(is.null(selection)) if(rdo=="") if(length(selQuantiVar)==1) selection <- which(rownames(coord.Q) %in% selQuantiVar)		
      if(length(selection)==0) stop("There are not selected elements to plot")		
  

      
      
      if(new.plot) if(!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY")))   
        dev.new(width = min(14, 8 * (xlim[2] - xlim[1])/(ylim[2] - ylim[1])), height = 8) # Outside RStudio 
      else if(new.plot==TRUE) dev.new()
      
      
      coo <- labe <- coll <- ipch <- fonte <- NULL	
    
      
      
      
      
      if (label) labe <- rownames(coord.Q)		
      else labe <- rep("", nrow(coord.Q))	
      #		if(is.null(col.quanti)) col.quanti="blue"	
      coll <- rep(col.quanti, nrow(coord.Q))			
      #      	coo <- coord.Q		
      plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, col = "white", asp = 1, main = title, ...)
      
      if(selCov==FALSE) {				
        x.cercle <- seq(-1, 1, by = 0.01)	
        y.cercle <- sqrt(1 - x.cercle^2)	
        lines(x.cercle, y = y.cercle, ...)		
        lines(x.cercle, y = -y.cercle, ...)		
      } # Final selCov==FALSE
      abline(v = 0, lty = 2, ...)		
      abline(h = 0, lty = 2, ...)		
      
      

      col.var.partial <- col.group					

    
   
      if(partial==FALSE | partial=="ALL"){				
        posi <- NULL
        
        if(length(col.quanti)== 1)
        coll <- rep(col.quanti, nrow(coord.Q))		
        else {
          coll <- rep("blue", nrow(coord.Q))
          coll[1:length(col.quanti)] <- col.quanti
        }
        
        coo <- coord.Q					
        
        if (label) 					
          labe <- rownames(coord.Q)					
        else labe <- rep("", nrow(coord.Q))					
        

        if (!is.null(selection)) {					
          if (is.numeric(unselect)) 					
            coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% 					
                                                                              selection)])), alpha = 255 * (1 - unselect),maxColorValue = 255)					
          else coll[!((1:length(coll)) %in% selection)] = unselect					
          labe[!((1:length(coll)) %in% selection)] <- ""					
        }				
        

        
        for (v in 1:nrow(coord.Q)) {					
          arrows(0, 0, coord.Q[v, 1], coord.Q[v, 2], 					
                 length = 0.1, angle = 15, code = 2, col = coll[v])					
          if (label) {					
            if (abs(coord.Q[v, 1]) > abs(coord.Q[v,2])) {					
              if (coord.Q[v, 1] >= 0) 					
                posi <- c(posi, 4)					
              else posi <- c(posi, 2)					
            }					
            else {					
              if (coord.Q[v, 2] >= 0) 					
                posi <- c(posi, 3)					
              else posi <- c(posi, 1)					
            }					
          }		
        }


          if (any(labe != "")) {					
            autoLab <- FALSE					
            if (auto.Lab == "auto")  autoLab = (length(which(labe != "")) < 50)					
            if (auto.Lab == TRUE) autoLab <- TRUE					
            if (autoLab == TRUE){ 					
              autoLab(coo[labe != "", 1], y = coo[labe != "", 2], labels = labe[labe != ""],					
                      col = coll[labe != ""], cex=cex, ...)	}				
            else					
              text(coo[labe != "", 1], y = coo[labe != "", 2], labels = labe[labe != ""], 					
                   pos = posi[labe != ""], col = coll[labe != ""], cex=cex,...)					
         #  }		
      }	# Final partial				
  
} # End Temporal      

      
      



      #==============================================================	
      if(partial==TRUE | partial=="ALL"){	
        n.group <- res.mfagalt$call$num.groups	
        n.var <- nrow(res.mfagalt$quanti.var$coord)	
        nam.group <- res.mfagalt$call$name.groups	
        coord.Q <- NULL			
        
        for(i in 1:n.group) {
          coord.Q <- rbind(coord.Q, res.mfagalt$quanti.var$coord.partial[[i]])
        }
        if(selCov==FALSE) {			
          coord.Q  <- coord.Q/sqrt(apply(coord.Q^2,1,sum))
          xlim <- ylim <- c(-1, 1)	
        } # End selCov			
        
        coord.Q <- cbind(coord.Q,rep(nam.group, each=n.var))		
        coord.Q <-  coord.Q[,c(axes,ncol(coord.Q))]	

        posi <- NULL			
        coll <- rep(col.var.partial, each=n.var)			# Basic colors of variables by groups
        coo <- coord.Q	
        if (label) 			
          labe <- rownames(coord.Q)			
        else labe <- rep("", nrow(coord.Q))			
       
       selection <-  which(rownames(coord.Q) %in% rownames(coord.Q)[selection] & coord.Q[,3] %in% selGroup)

        
        
        if (!is.null(selection)) {			
          if (is.numeric(unselect)) 			
            coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% selection)])),
                                                           alpha = 255 * (1 - unselect),maxColorValue = 255)			
          else coll[!((1:length(coll)) %in% selection)] = unselect			
          labe[!((1:length(coll)) %in% selection)] <- ""			
        }			
        ipch <- rep(25, nrow(coord.Q))				
        fonte <- rep(1, nrow(coord.Q))				
        
        for (v in 1:nrow(coord.Q)) {		
          arrows(0, 0, as.numeric(coord.Q[v, 1]), as.numeric(coord.Q[v, 2]), 		
                 length = 0.1, angle = 15, code = 2, col = coll[v])		
          #   if (shadowtext) points(coo[, 1], y = coo[, 2], pch = ipch, col = coll,cex=cex, ...)				
          
        } # End for

        posi <- NULL			
        coll <- rep(col.var.partial, each=n.var)			
        coo <- coord.Q		
        if (label) 			
          labe <- rownames(coord.Q)			
        else labe <- rep("", nrow(coord.Q))			
        
        if (!is.null(selection)) {			
          if (is.numeric(unselect)) 			
            coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% 			
                                                                              selection)])), alpha = 255 * (1 - unselect),maxColorValue = 255)			
          else coll[!((1:length(coll)) %in% selection)] = unselect		
          
          # labe[!((1:length(coll)) %in% selection)] <- ""
          # labe[((1:length(coll)) %in% selection)] <- rownames(coord.Q)[selection]
        }			
        ipch <- rep(25, nrow(coord.Q))				
        fonte <- rep(1, nrow(coord.Q))				
        
        for (v in 1:nrow(coord.Q)) {		
          arrows(0, 0, as.numeric(coord.Q[v, 1]), as.numeric(coord.Q[v, 2]), 		
                 length = 0.1, angle = 15, code = 2, col = coll[v])		
       #   if (shadowtext) points(coo[, 1], y = coo[, 2], pch = ipch, col = coll,cex=cex, ...)				
          
        } # End for
        if (any(labe != "")) {				
          labe <- rownames(coord.Q)
          autoLab(as.numeric(coord.Q[labe != "", 1]), y = as.numeric(coord.Q[labe != "", 2]), 
                  labels = labe[labe != ""], col = coll[labe != ""], font = fonte[labe != ""], cex=cex, ...)					
        }


             

        nselGr <- which(old.name.groups %in% selGroup)
       
        
        if(legend) 
         legend(pos.legend,label.group[nselGr],pch=NULL,cex=cex*1.5,text.col=col.var.partial[ nselGr],
                          col=col.var.partial[ nselGr])	  
                 
      } # End partial
      #==============================================================	

      
      
    } ########################## Final selQuantiVar
    
    

  } # Final multiple case
  
}


