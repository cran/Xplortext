#' @import stringr
#' @export
LexChar <- function(object, proba = 0.05, maxDocs=20, maxCharDoc = 10, maxPrnDoc=100)
 {
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

#=================================================================		
###### Step 2. Extraction of the characteristic words #### 				     
resCharWord <- descfreq(DocTerm,proba=proba)


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
 corpus <- data.frame(corpus)
 rownames(corpus) <- rownames(base)
 corpus <- corpus[rownames(object$var.agg),]
 motsval <- sapply(resCharWord,function(x) if(!is.null(x))  x[order(row.names(x)),6,drop=FALSE],simplify = TRUE)
lisresult <- vector(mode="list",length=nlev)

	for (ilev in 1:nlev)
	{
	respsel <- which(object$var.agg == vlev[ilev])
	ntrep <- min(maxCharDoc, length(respsel))
	lisresult[[ilev]] <- data.frame(nlev,3)
	SourceTermcurs <- SourceTerm[respsel,,drop=FALSE ]
	ly<-as.matrix(rep(0,ncol(DocTerm)))
	rownames(ly)<-colnames(DocTerm)

if(is.null(motsval[[ilev]])) {  lisresult[ilev] <- NULL } else {
	motsvalcurs<-as.matrix(motsval[[ilev]])
	ly[rownames(motsvalcurs),1] <- motsvalcurs[,1]
      a <- crossprod(ly,t(SourceTermcurs))
	b <- rowSums(SourceTermcurs)
	repvaltest <- a
	repvaltest[b > 0] <- a[b > 0]/b[b > 0]
	ordrep <- order(repvaltest, decreasing = "TRUE")
	   for (i in 1:ntrep)
	   {
	   lisresult[[ilev]][i,1] <- rownames(object$SourceTerm)[respsel[ordrep[i]]]
	   lisresult[[ilev]][i,2] <- repvaltest[ordrep[i]]
	   lisresult[[ilev]][i,3] <- corpus[respsel[ordrep[i]]]
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
	cat(paste("\n",str_pad(rownames(resCharWord[[ilev]])[j],width=20, side="left",pad=" "),"  ",sep=""))
	cat(paste(str_pad(format(round(resCharWord[[ilev]][j,1],3),nsmall=3),width=6,side="left",pad=" "),"  ",sep=""))
	cat(paste(str_pad(format(round(resCharWord[[ilev]][j,2],3),nsmall=3),width=6,side="left",pad=" "),"  ",sep=""))
	cat(paste(str_pad(resCharWord[[ilev]][j,3],width=9,side="left",pad=" "),"  ",sep=""))
	cat(paste(str_pad(resCharWord[[ilev]][j,4],width=9,side="left",pad=" "),"  ",sep=""))
	cat(paste(str_pad(format(round(resCharWord[[ilev]][j,5],5),nsmall=5,digits=5, format="fg"),width=7,side="left",pad=" "),"  ",sep=""))
	cat(paste(str_pad(format(round(resCharWord[[ilev]][j,6],6),nsmall=4, format="%9.5f"),width=10)))
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
  cat("\n",paste(str_pad(format(round(lisresult[[i]][j,2],3),	
  nsmall=3),width=6,side="left",pad=" "),"  ----  ",	
  str_pad(as.character(lisresult[[i]][j,1]), width=13, side="left", pad=" "), "  ",	
  str_pad(as.character(substr(lisresult[[i]][j,3], 1,maxPrnDoc) ),width=maxPrnDoc,side="right", pad=" "),sep=""))	
}
  }	
  cat("\n\n")	
  }	
} }	
return(res)
}

