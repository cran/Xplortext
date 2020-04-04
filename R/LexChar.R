#' @import stringr
#' @export
LexChar <- function(object, proba = 0.05, maxDocs=20, maxCharDoc = 10, maxPrnDoc=100, marg.doc="before")
 {
  options(stringsAsFactors = FALSE)
# maxDocs: Maximum number of input documents or aggregate categories
# maxCharDoc: Maximum number of characteristic documents to show

##################  main function  ##################
###      Step 1. Verifying the correctness of the arguments
if (!inherits(object,"TextData")) stop("Object should be TextData class")
if(proba<0|proba>1) stop("proba should be between 0 and 1")
  

###  To know if it is an aggregate analysis
bvaragg <- ifelse(object$info$name.var.agg[[1]]=="",FALSE,TRUE)
bCharDoc <- ifelse(is.null(maxCharDoc),FALSE,TRUE)
bextract <- ifelse((bvaragg & bCharDoc), TRUE,FALSE)

if(bCharDoc & !bvaragg) 
	warning("The documents are not aggregate, \nextracting characteristic documents has no sense")

###   To control if there is not too many documents
DocTerm <- as.matrix(object$DocTerm)
if(nrow(DocTerm)>maxDocs) stop("Computing the characteristic words when more than ", maxDocs, " documents is not allowed")
vlev <- rownames(DocTerm)			
nlev <- nrow(DocTerm)



############### descfreq_New function #############################################
descfreq_NewChar <- function (data, proba = 0.05, marge.li, marge.col) 				
{
  lab.sauv <- lab <- colnames(data)				
  for (i in 1:length(lab)) lab[i] = gsub(" ", ".", lab[i])				
  colnames(data) = lab				
  
  old.warn = options("warn")				
  options(warn = -1)				
  # if (is.null(by.quali))  marge.li = apply(donnee, 1, sum)
 # marge.li <- row.INIT
  nom = tri = structure(vector(mode = "list", length = nrow(data)), names = rownames(data))				
#  marge.col=col.INIT # marge.col = apply(donnee, 2, sum)			
  marge.li <- as.vector(marge.li)
  marge.col <- as.vector(marge.col)
  SumTot=sum(marge.li)				
  nr <- nrow(data)  				
  nc <- ncol(data) 				

  for (j in 1:nr) {				
    aux3 <- marge.li[j]/SumTot				# % Ocurrences before
    
    for (k in 1:nc) {				
      aux2 <- data[j, k]/marge.col[k]	
      if (aux2 > aux3) # sobrerepresentada				
        aux4 = phyper(data[j, k] - 1, marge.col[k], 				
                      SumTot - marge.col[k], marge.li[j],lower.tail = FALSE) * 2				
      else aux4 = phyper(data[j, k], marge.col[k], SumTot - marge.col[k], marge.li[j]) * 2				

      if (aux4 > 1) 				
        aux4 <- 2 - aux4				
      if (aux4 < proba) {				
        aux5 = (1 - 2 * as.integer(aux2 > aux3)) * qnorm(aux4/2)				
        aux1 = data[j, k]/marge.li[j]				
        tri[[j]] = rbind(tri[[j]], c(aux1 * 100, sum(marge.col[k])/SumTot * 				
                                       100, data[j, k], marge.col[k], aux4, aux5))				
        nom[[j]] = rbind(nom[[j]], c(colnames(data)[k], colnames(data)))				
      }				
    }				
  }				
  
  for (j in 1:nrow(data)) {				
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








#=================================================================		
###### Step 2. Extraction of the characteristic words #### 				     
# resCharWord <- descfreq(DocTerm,proba=proba)
if(marg.doc=="before") row.INIT <- object$rowINIT else row.INIT <- data.frame(Ocurrences.After=rowSums(DocTerm))
row.INIT <- row.INIT[,1]
col.INIT <- as.vector(colSums(DocTerm))
resCharWord <- descfreq_NewChar(DocTerm, proba = 0.05, row.INIT, col.INIT) 				



###### Step  3. Extraction of the characteristic documents in the case of an aggregate table		
	
if(bextract & maxCharDoc>0)   	
 {	
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

 
	for (ilev in 1:nlev)
	{
	  # repsel doc position for each group, starting the firs ilev=1)
	respsel <- which(object$var.agg == vlev[ilev])  # vlev vector con los nombres de las modalidades agregadas
	# ntrep maximum nombre of documents to print, minumum between maxPrnDoc and the size of the group
	ntrep <- min(maxCharDoc, length(respsel))
	lisresult[[ilev]] <- data.frame(nlev,3)  # Nombre of the group
	# SourceTerm the rows for non aggregated table
	SourceTermcurs <- SourceTerm[respsel,,drop=FALSE ]   
	ly<-as.matrix(rep(0,ncol(DocTerm)))
	# Rows are the vocabulary of after words
	rownames(ly)<-colnames(DocTerm)
  
	
  if(is.null(motsval[[ilev]])) {  lisresult[ilev] <- NULL } else {
    # There are significant words for this category
	  motsvalcurs<-as.matrix(motsval[[ilev]])
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
}
	}
 names(lisresult) <- vlev
 }


###### Content of the result of the function according to non-aggregate or aggregate documents
bres <- FALSE


if(bextract==TRUE) if(maxCharDoc>0) bres<- TRUE
	ifelse(bres, 
	res <- list(CharWord=resCharWord,CharDoc=lisresult,Proba=proba),
	res <- list(CharWord=resCharWord,Proba=proba))
  class(res) <- c("LexChar", "list")


###### Printing characteristic words on the screen      			
cat("\n\nCHARACTERISTIC WORDS\n(DETAILED INFORMATION)\n")
 for (ilev in 1:nlev)
 {			
 # Intern % glob % Intern freq Glob freq   p.value v.test
 cat(paste("\n","Group",ilev,": ",rownames(DocTerm)[ilev],sep="", "\n"))
 cat(paste(rep("-",40)))

if(is.null(resCharWord[[ilev]])) { cat(paste("\n"))} else { 
 cat("\n      Word          Intern %  glob % Intern freq Glob freq  p.value    v.test \n")
 cat(paste(rep("-",40)))
	for(j in 1:nrow(resCharWord[[ilev]]))
	{
	cat(paste("\n",stringr::str_pad(rownames(resCharWord[[ilev]])[j],width=20, side="left",pad=" "),"  ",sep=""))
	cat(paste(stringr::str_pad(format(round(resCharWord[[ilev]][j,1],3),nsmall=3),width=6,side="left",pad=" "),"  ",sep=""))
	cat(paste(stringr::str_pad(format(round(resCharWord[[ilev]][j,2],3),nsmall=3),width=6,side="left",pad=" "),"  ",sep=""))
	cat(paste(stringr::str_pad(resCharWord[[ilev]][j,3],width=9,side="left",pad=" "),"  ",sep=""))
	cat(paste(stringr::str_pad(resCharWord[[ilev]][j,4],width=9,side="left",pad=" "),"  ",sep=""))
	cat(paste(stringr::str_pad(format(round(resCharWord[[ilev]][j,5],5),nsmall=5,digits=5, format="fg"),width=7,side="left",pad=" "),"  ",sep=""))
	cat(paste(stringr::str_pad(format(round(resCharWord[[ilev]][j,6],6),nsmall=4, format="%9.5f"),width=10)))
	}
cat("\n")			
 }	}		


###### Printing characteristic documents on the screen	
if(bres) 	
 {	
  cat("\n\nCHARACTERISTIC DOCUMENTS\n(WORDS FREQUENCY CRITERION)\n")	
  for(i in 1:length(lisresult))	
  {	
  cat(paste("\n","GROUP ",i,": ",names(lisresult)[i],sep="", "\n"))	
  cat(paste(rep("-",35)))	
if(is.null(resCharWord[[i]])) { cat(paste("\n"))} else { 
  cat("\nCLASSIFICATION       DOCUMENT     CHARACTERISTIC\n")	
  cat("CRITERION                         DOCUMENT\n")	
  cat(paste(rep("-",35)))	
  for(j in 1:nrow(lisresult[[i]]))	
  {	

if(lisresult[[i]][j,2]>0){
  cat("\n",paste(stringr::str_pad(format(round(lisresult[[i]][j,2],3),	
  nsmall=3),width=6,side="left",pad=" "),"  ----  ",	
  stringr::str_pad(as.character(lisresult[[i]][j,1]), width=13, side="left", pad=" "), "  ",	
  stringr::str_pad(as.character(substr(lisresult[[i]][j,3], 1,maxPrnDoc) ),width=maxPrnDoc,side="right", pad=" "),sep=""))	
}
  }	
  cat("\n\n")	
  }	
} }	
return(res)
}

