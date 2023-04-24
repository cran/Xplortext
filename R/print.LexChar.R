#' @export
print.LexChar <- function (x,file = NULL, sep = ";", dec=".",...) 
{
  
  if (!inherits(x, "LexChar")) 
    stop("x object must be LexChar class")

  df.res <- data.frame(name = character(),
                       description = character(),
                       stringsAsFactors = FALSE)

  if("CharWord" %in% names(x)) df.res[nrow(df.res)+1,] <- c("$CharWord", "Words characterizing the documents")
  if("stats" %in% names(x)) df.res[nrow(df.res)+1,] <- c("$stats", "Association statistics of the lexical table")
  if("CharDoc" %in% names(x)) df.res[nrow(df.res)+1,] <- c("$CharDoc", "Documents characterizing categories of aggregate table")
  if("Proba" %in% names(x)) df.res[nrow(df.res)+1,] <- c("$Proba", "Threshold on the p-value used when selecting the characteristic words
  (by default 0.05")
                                                                         
  if(!is.null(x$quali$CharWord)) df.res[nrow(df.res)+1,] <-c("$quali$CharWord", "Words characterizing the qualitative variables and their categories")
  if(!is.null(x$quali$stats)) df.res[nrow(df.res)+1,] <- c("$quali$stats", "Association statistics for vocabulary and qualitative variables")
  if(!is.null(x$quali$CharDoc)) df.res[nrow(df.res)+1,] <-c("$quali$CharDoc", "Documents characterizing categories of aggregate table from qualitative variables and their categories")
 
  if(!is.null(x$quanti$CharWord)) df.res[nrow(df.res)+1,] <- c("$quanti$CharWord", "Provides characteristic quantitative variables for each word") 
  if(!is.null(x$quanti$stats)) df.res[nrow(df.res)+1,] <- c("$quanti$stats", "Association statistics for vocabulary and quantitative variables")
  
  df.res <- as.matrix(df.res) 
  index <- nrow(df.res)
  res <- array("", c(index, 2), list(1:index, c("name", "description")))
  for(i in 1:index) res[i,] <-df.res[i,]
  cat("\n*The results are available in the following objects:\n\n")
  print(res[c(1:index),])

  sink.reset <- function(){
    for(i in seq_len(sink.number())){sink()}}
  sink.reset
  
  if (is.null(file))
    cat("(*) To obtain more detailed information use print to file (print.LexChar)")
  else 
    {
      sink(file)
  
  fCharWord <- function(z)  {
    nobjects <- length(z) 

    cat("\n\nCHARACTERISTIC WORDS FOR EACH DOCUMENT\n(DETAILED INFORMATION)\n\n")
    cat("\nOver_used_words are ordered from the most characteristic word to the less one\n")
    cat("Infra_used_words are ordered from the most anti-characteristic word to the less one\n\n")
    stit <- t(c("Word","Intern %","glob %","Intern freq","Glob freq","p.value", "v.test"))
  
    
    for(i in 1:nobjects){	
      out1 <- c(names(z[i]), "Over_used_words")	
      out2 <- c(names(z[i]), "Infra_used_words")	
      
      #    write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
      temp <- z[[i]]	
      posic <- temp[,"v.test"]>0	
      
      out.over <- data.frame(temp[posic,,drop=FALSE])	
      out.infra <- data.frame(temp[!posic,,drop=FALSE])	

      if(nrow(out.over)>0){	
        write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
        write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
        write.table(out.over, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec) 
        cat("\n")	
      }
      
      if(nrow(out.infra)>0){	
        write.table(out2, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
        write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
        write.table(out.infra, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec) 
        cat("\n")	
      }
    }
  }
  
  fCharQuanti<- function(x)  {
    tac <- NULL
    strcolnames<- c("GlobalAverage", "AverageWord","Differ.", "pvalue", "Word", "Variable")
    for(i in 1:length(x)) {
      t1<- as.data.frame(x[i,drop=FALSE])
      t2 <- data.frame(t1,rep(names(x[i]),length(x[i])), rownames(t1))
      colnames(t2) <- strcolnames

      if(i==1) tac <- t2 else tac <- rbind(tac,t2)
      SP <- split(tac,f=tac$Variable, drop=FALSE)
      str.colnames<- c("Word", "GlobalAverage", "AverageWord","Differ.", "pvalue")
      empty_list = structure(vector(mode = "list", length = length(SP)), names = names(SP))		
      

      for(i in 1:length(SP)) {
        t1<- as.data.frame(SP[i,drop=FALSE])
        t2.pos <- t1[t1[,3]>0, ,drop=FALSE]
        t2.neg <- t1[t1[,3]<0, ,drop=FALSE]
        if(nrow(t2.pos)>0) {
          t2.pos <- t2.pos[order(t2.pos[,4]),,drop=FALSE]
          rownames(t2.pos) <- paste0("P", c(1:nrow(t2.pos)))
          t3.pos <- t2.pos[,c(5,1:4)]
          colnames(t3.pos) <- str.colnames
          empty_list[[i]]$posit <- t3.pos
        }
        if(nrow(t2.neg)>0) {
          t2.neg <- t2.neg[order(t2.neg[,4]),,drop=FALSE]
          rownames(t2.neg) <- paste0("N", c(1:nrow(t2.neg)))
          t3.neg <- t2.neg[,c(5,1:4)]
          colnames(t3.neg) <- str.colnames
          #    print(t3.neg)
          empty_list[[i]]$negat <- t3.neg
        }
      } # End for
    }
    return(empty_list)
  } # End fCharQuanti
  
  
  #################  Check CharWord
  if(!is.null(x$CharWord)) {
    
  fCharWord(x$CharWord)
  nobjects <- length(x$CharWor)
  
  cat("\n\nDOCUMENTS (OR AGGREGATE DOCS) RELATED WITH A WORD\n")	
  tw <- NULL	
  for(i in 1:nobjects){	
    tw <- c(tw,rownames(x$CharWord[[i]]))	
  }	
  tw <- paste(unique(tw))	
  tw <- sort(tw, decreasing = FALSE)	
  ntw <- length(tw)	
  
  
  for(i in 1:ntw){	
    vsearch <- tw[i]	
    vt_ <- NULL; posit <- NULL; negat <- NULL  	
    for(j in 1:nobjects){	
      cat_ <- which(rownames(x$CharWord[[j]]) %in% vsearch)	
      if(length(cat_)>0){	
        vt  <- x$CharWord[[j]][cat_,6]	
        if(vt>0) posit <- c(posit,sep,names(x$CharWord[j])) else	
          negat <- c(negat,sep,names(x$CharWord[j]))	
      }	
    }	
    cat("\n",vsearch)	
    if(length(posit)>0)  cat("\n",sep, "Over", posit) 	
    if(length(negat)>0)  cat("\n",sep, "Infra", negat) 	
    cat("\n")	
  }	
    }
  
  
  fCharWord2 <- function(z,i)  {
    nobjects <- length(z) 
    
  #  cat("\n\nCHARACTERISTIC WORDS FOR EACH DOCUMENT\n(DETAILED INFORMATION)\n\n")
  #  cat("\nOver_used_words are ordered from the most characteristic word to the less one\n")
  #  cat("Infra_used_words are ordered from the most anti-characteristic word to the less one\n\n")
    stit <- t(c("Word","Intern %","glob %","Intern freq","Glob freq","p.value", "v.test"))
    
    
  #  for(i in 1:nobjects){	
      out1 <- c(names(z[i]), "Over_used_words")	
      out2 <- c(names(z[i]), "Infra_used_words")	
      
      #    write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
      temp <- z[[i]]	
      posic <- temp[,"v.test"]>0	
      
      out.over <- data.frame(temp[posic,,drop=FALSE])	
      out.infra <- data.frame(temp[!posic,,drop=FALSE])	
      
      if(nrow(out.over)>0){	
        write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
        write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
        write.table(out.over, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec) 
        cat("\n")	
      }
      
      if(nrow(out.infra)>0){	
        write.table(out2, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
        write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
        write.table(out.infra, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec) 
        cat("\n")	
      }
   # }
  }
  

  
  
    
    #################  Check of stats exists  
  if(!is.null(x$stats)) {
    cat("\n\n\nSTATISTIC ASSOCIATION LEXICAL TABLE DOCUMENTS (OR AGGREGATE DOCS) AND VOCABULARY\n\n")
      stit <- t(c("",colnames(x$stat)))
      write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      write.table(x$stats, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec)
    }

  
  
  #################  Check CharDoc
  if(!is.null(x$CharDoc)) {
    cat("\n\n\nCHARACTERISTIC SOURCE-DOCUMENTS\n")
    #     print(object$CharDoc)
    stit <- t(names(x$CharDoc[[1]]))
    for(i in 1:length(x$CharDoc)){	
      cat("\n")
      write.table(names(x$CharDoc[i]), quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      write.table(x$CharDoc[[i]][1:3], quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
    }	
  }  
  
  ###########  Qualitative contextual variables
  if(!is.null(x$quali$stats)){	
    cat("\n\nASSOCIATION STATISTICS OF AGGREGATE CONTEXTUAL QUALITATIVE TABLE\n")
    stit <- t(c("",colnames(x$quali$stats)))
    write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
    write.table(x$quali$stats, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec)
  } # End Vocab$quali$CharWord$stats
  
  if(!is.null(x$quali$CharWord)){	
    for(i in 1:length(x$quali$CharWord)) {
      cat("\n\nCharacteristic of aggregate-documents for qualitative variable: ", names(x$quali$CharWord[i]), "\n")
      cat("=====================================================================================")
      fCharWord2(x$quali$CharWord,i)
    }
    
    #################  Check CharDoc
    if(!is.null(x$quali$CharDoc)) {
      cat("\n\n\nCHARACTERISTIC SOURCE-DOCUMENTS\n")
      #     print(object$CharDoc)
      stit <- t(names(x$quali$CharDoc[[1]]))
      for(i in 1:length(x$quali$CharDoc)){	
        cat("\n")
        write.table(names(x$quali$CharDoc[i]), quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
        write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
        write.table(x$quali$CharDoc[[i]][1:3], quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      }	
    }   
  } # End CharWord

  

  ###########  Quantitative contextual variables
  if(!is.null(x$quanti$stats)){	
    cat("\nRELATIONSHIPS BETWEEN VOCABULARY AND CONTEXTUAL QUANTITATIVE VARIABLES\n")
    cat("\nStatistics for quantitative variables")
    cat("\n-------------------------------------\n")
    stit <- t(c("",colnames(x$quanti$stats)))
    write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
    write.table(x$quanti$stats, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec)
    

  #  print(object$quanti$CharWord)
    res.quanti <- fCharQuanti(x$quanti$CharWord)
    num.var.quanti <- length(res.quanti)
    stit<- t(c("Word", "GlobalAverage", "AverageWord","Differ.", "pvalue"))
    
       for(i in 1:num.var.quanti) {
      str.name.var <- rownames(x$quanti$stats)[i]
      cat("\n\nCharacteristic quantitative variables for each word. Variable:", str.name.var)
      cat("\n----------------------------------------------------------------------------------------\n")
      
      
      if(length(res.quanti[[str.name.var]]$posit)>0) {
        cat(paste0("\n$",str.name.var, "$posit\n"))
        write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
        write.table(res.quanti[[str.name.var]]$posit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      }
   
      if(length(res.quanti[[str.name.var]]$negat)>0) {
        cat(paste0("\n$",str.name.var, "$negat\n"))
        write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
        write.table(res.quanti[[str.name.var]]$negat, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      }    
      
    }
    
  } # End quanti
  
  cat("\n------- Proba threshold -----------\n")
  cat(paste("Proba=", x$Proba))
  
  
  sink()
    }  # End sink file    
  if (!is.null(file)) cat(paste("(*) All the results are in file", file))
}
