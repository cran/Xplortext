###' @export
print.LexChar <- function (x,file = NULL, sep = ";", dec=".",...) 
{

 if (!inherits(x, "LexChar")) 
        stop("x object should be LexChar class")


sink.reset <- function(){
    for(i in seq_len(sink.number())){sink()}}
sink.reset

indice <- 8
cat("*The results are available in the following objects:\n\n")
 res.LexChar <- x
 res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
 res[1,] <- c("$CharWord", "Words characterizing the documents")
 res[2,] <- c("$CharDoc", "Documents with a characterizing word")
 res[3,] <- c("$stats", "Association statistics of the lexical table")
 res[4,] <- c("$Vocab", "Information about the vocabulary")
 res[5,] <- c("$Vocab$quali$CharWord", "Provides the qualitative variables and their categories")
 res[6,] <- c("$Vocab$quali$stats", "Provides association statistics for vocabulary and qualitative variables")
 res[7,] <- c("$Vocab$quanti$CharWord", "Provides characteristic quantitative variables for each word") 
 res[8,] <- c("$Vocab$quali$stats", "Provides association statistics for vocabulary and quantitative variables")
  print(res[c(1:indice),])
 cat("\n")
 

if (is.null(file))
 print("(*) To obtain more detailed information use print to file")
 
if (!is.null(file)) {
 sink(file)

  
fCharWord <- function(z)  {
  nobjects <- length(z) 
  cat("\nOver_used_words are ordered from the most characteristic word to the less one\n")
  cat("Infra_used_words are ordered from the most anti-characteristic word to the less one\n\n")
  stit <- t(c("Word","Intern %","glob %","Intern freq","Glob freq","p.value", "v.test"))

  for(i in 1:nobjects){	
    
    out1 <- c(names(z[i]), "Over_used_words")	
    out2 <- c(names(z[i]), "Infra_used_words")	

    
#    write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
    temp <- z[[i]]	
    posic <- temp[,"v.test"]>0	
    out3 <- data.frame(temp[posic,,drop=FALSE])	

    
    if(nrow(out3)>0){	
      write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
      
      write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
      write.table(out3, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec)

    
    cat("\n")	
 #   if(nrow(out3)>0){	
    write.table(out2, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
    out3 <- data.frame(temp[!posic,,drop=FALSE])	
    out3 <- out3[order(out3["v.test"]), ]	
 
    
  #  if(nrow(out3)>0){	
      write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      write.table(out3, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec) 
    cat("\n")	
} # End if(nrow(out3)>0)
  }
}
  
 if(!is.null(x$CharWord)) fCharWord(x$CharWord)

if(!is.null(x$CharDoc)){	
  nobjects <- length(x$CharDoc) 
  
  # Print documents	
  stit <- t(c("Classification criterion","Document","Characteristic documents"))	
  cat("\nCharacteristic documents (Words frequency criterion)\n")	
  
  for(i in 1:nobjects){	
    write.table(names(x$CharWord[i]), quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)	
    if(nrow(x$CharDoc[[i]][2])>0){	
      out3 <- c(format(round(x$CharDoc[[i]][2],3),nsmall=3),x$CharDoc[[i]][1],
                x$CharDoc[[i]][3])
      write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      write.table(out3, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
      cat("\n")}
  }


cat("\n", "Documents related to a word")	
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

if(!is.null(x$stats)){	
  cat("\n\nAssociation statistics of the lexical table\n")
  stit <- t(c("",colnames(x$stat)))
  write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
  write.table(x$stats, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec)
}

if(!is.null(x$Vocab$quali$CharWord)){	
  for(i in 1:length(x$Vocab$quali$CharWord)) {
      fCharWord(x$Vocab$quali$CharWord[[i]])
}
  } # End CharWord

if(!is.null(x$Vocab$quali$stats)){	
  cat("\n\nAssociation statistics of the aggregated contextual qualitative  table\n")
  stit <- t(c("",colnames(x$Vocab$quali$stats)))
  write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
  write.table(x$Vocab$quali$stats, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec)
} # End Vocab$quali$CharWord$stats

if(!is.null(x$Vocab$quanti$CharWord)){	
  cat("\n\nContextual quantitative variables related with words\n")
  stit <- t(c("",colnames(x$Vocab$quanti$CharWord[[1]])))
  for(i in 1:length(x$Vocab$quanti$CharWord)) {
    cat("\nWord: ", names(x$Vocab$quanti$CharWord)[i],"\n" )
    write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
    write.table(x$Vocab$quanti$CharWord[[1]], quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec)
  }
  cat("\n")
  
  if(!is.null(x$Vocab$quanti$stats)){	
    cat("\n\nRelationship between vocabulary and contextual quantitative variables\n")
    stit <- t(c("",colnames(x$Vocab$quanti$stats)))
    write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE,dec=dec)
    write.table(x$Vocab$quanti$stats, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE,dec=dec)
  } # End Vocab$quali$CharWord$stats
  
  
}






sink()
if (!is.null(file)) print(paste("All the results are in file", file))
}
}
