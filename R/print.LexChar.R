###' @export
print.LexChar <- function (x,file = NULL, sep = ";", ...) 
{

 if (!inherits(x, "LexChar")) 
        stop("x object should be LexChar class")


sink.reset <- function(){
    for(i in seq_len(sink.number())){sink()}}
sink.reset

indice <- 2
cat("*The results are available in the following objects:\n\n")
 res.LexChar <- x
 res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
 res[1,] <- c("$CharWord", "Words characterizing the documents")
 res[2,] <- c("$CharDoc", "Documents with a characterizing word")
  print(res[c(1:indice),])
 cat("\n")
if (is.null(file))
 print("(*) To obtain more detailed information use print to file")

if (!is.null(file)) {
 sink(file)

# Print words
  nobjects <- length(x$CharWord)
  cat("\nOver_used_words are ordered from the most characteristic word to the less one\n")
  cat("Infra_used_words are ordered from the most anti-characteristic word to the less one\n\n")
  stit <- t(c("Word","Intern %","glob %","Intern freq","Glob freq","p.value", "v.test"))
for(i in 1:nobjects){	
 out1 <- c(names(x$CharWord[i]), "Over_used_words")	
 out2 <- c(names(x$CharWord[i]), "Infra_used_words")	
 write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)	
 temp <- x$CharWord[[i]] 	
 posic <- temp[,"v.test"]>0	
 out3 <- data.frame(temp[posic,,drop=FALSE])	
	
 if(nrow(out3)>0){	
 write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)	
 write.table(out3, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE)}	
	
 cat("\n")	
 write.table(out2, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)	
 out3 <- data.frame(temp[!posic,,drop=FALSE])	
 out3 <- out3[order(out3["v.test"]), ]	
  if(nrow(out3)>0){	
	write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
	write.table(out3, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE)}
cat("\n")	
}	
	
if(!is.null(x$CharDoc)){	
# Print documents	
 stit <- t(c("Classification criterion","Document","Characteristic documents"))	
 cat("\nCharacteristic documents (Words frequency criterion)\n")	
	
for(i in 1:nobjects){	
 write.table(names(x$CharWord[i]), quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)	
 if(nrow(x$CharDoc[[i]][2])>0){	
	out3 <- c(format(round(x$CharDoc[[i]][2],3),nsmall=3),x$CharDoc[[i]][1],
       x$CharDoc[[i]][3])
	write.table(stit, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
	write.table(out3, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
	cat("\n")}
}}	
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
sink()
if (!is.null(file)) print(paste("All the results are in file", file))
}
}
