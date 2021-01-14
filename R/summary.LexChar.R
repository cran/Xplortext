###' @export
summary.LexChar <- function (object, CharWord=TRUE, stats=TRUE, CharDoc=TRUE, Vocab=TRUE,  file = NULL, ...) 
{
if (!inherits(object, "LexChar")) stop("non convenient object")
  options(stringsAsFactors = FALSE)
 
  
  
  fCharWord<- function(x)  {
    nlev <- length(object$CharWord)
    for (ilev in 1:nlev)
    {			
      # Intern % glob % Intern freq Glob freq   p.value v.test
      cat(paste("\n",names(object$CharWord)[ilev],sep="", "\n"))
      cat(paste(rep("-",40)))
      
      if(is.null(object$CharWord[[ilev]])) { cat(paste("\n"))} else { 
        bposit <- TRUE
        cat("\n      Word          Intern %  glob % Intern freq Glob freq  p.value    v.test \n")
        cat(paste(rep("-",40)))
        for(j in 1:nrow(object$CharWord[[ilev]]))
        {
          if(bposit & object$CharWord[[ilev]][j,6] <0) {cat(paste("\n")) ; bposit<-FALSE }
          cat(paste("\n",stringr::str_pad(rownames(object$CharWord[[ilev]])[j],width=20, side="left",pad=" "),"  ",sep=""))
          cat(paste(stringr::str_pad(format(round(object$CharWord[[ilev]][j,1],3),nsmall=3),width=6,side="left",pad=" "),"  ",sep=""))
          cat(paste(stringr::str_pad(format(round(object$CharWord[[ilev]][j,2],3),nsmall=3),width=6,side="left",pad=" "),"  ",sep=""))
          cat(paste(stringr::str_pad(object$CharWord[[ilev]][j,3],width=9,side="left",pad=" "),"  ",sep=""))
          cat(paste(stringr::str_pad(object$CharWord[[ilev]][j,4],width=9,side="left",pad=" "),"  ",sep=""))
          cat(paste(stringr::str_pad(format(round(object$CharWord[[ilev]][j,5],5),nsmall=5,digits=5, format="fg"),width=7,side="left",pad=" "),"  ",sep=""))
          cat(paste(stringr::str_pad(format(round(object$CharWord[[ilev]][j,6],6),nsmall=6, format="fg"),width=10)))
          
        }
        cat("\n")			
      }	}		
    
  } # End fCharWord
    
  
  
  if(!is.null(file)) sink(file) else file=""
  

  
  #################  Comprobar que hay CharWord
  if("CharWord" %in% names(object)) if(CharWord==TRUE)   {
  ###### Printing characteristic words on the screen      	
  cat("\n\nCHARACTERISTIC WORDS\n(DETAILED INFORMATION)\n")
  fCharWord(object$CharWord) 
  } # end charword
  

  
  #################  Check of stats exists  
  if("stats" %in% names(object)) if(stats==TRUE) {
    cat("\n\n\nSTATISTIC ASSOCIATION LEXICAL TABLE DOCUMENTS (OR AGGREGATED DOCS,) AND VOCABULARY\n\n")
    print(object$stats)
  }
  
  #################  Check CharDoc
  if("CharDoc" %in% names(object)) if(CharDoc==TRUE)   {
    cat("\n\n\nCHARACERISTIC SOURCE-DOCUMENTS\n\n")
    print(object$CharDoc)
  }  
  
  if("Vocab" %in% names(object)) if(Vocab==TRUE)   {

      if(!is.null(object$Vocab$quali))  {
      cat("\n\nCHARACTERISTIC QUALITATIVE VARIABLES\n")
      cat("\nStatistics for qualitative variables\n")
      cat("\n-------------------------------------\n")
      print(object$Vocab$quali$stats)
      for(i in 1:length(object$Vocab$quali$CharWord)) {
        cat("\n\nCharacteristic words for qualitative variable: ", names(object$Vocab$quali$CharWord[i]), "\n")
        fCharWord(object$Vocab$quali$CharWord[i]) 
      }
      } # End quali

    
    if(!is.null(object$Vocab$quanti))  {
      cat("\n\nCHARACTERISTIC QUANTITATIVE VARIABLES\n")
      cat("\nStatistics for quantitative variables\n")
      cat("\n-------------------------------------\n")
      print(object$Vocab$quanti$stats)

      cat("\n\nCharacteristic quantitative variables for each word\n")
      cat("\n-----------------------------------------------------\n")
      print(object$Vocab$quanti$CharWord)
      
      cat("\n\nCharacteristic words for each quantitative variable\n ")
      cat("\n---------------------------------------------------\n")
      
      fCharQuanti<- function(x)  {
        tac <- NULL
        strcolnames<- c("GlobalAverage", "AverageWord","Difer.", "pvalue", "Word", "Variable")
        for(i in 1:length(x)) {
          t1<- as.data.frame(x[i,drop=FALSE])
          t2 <- data.frame(t1,rep(names(x[i]),length(x[i])), rownames(t1))
          colnames(t2) <- strcolnames
          if(i==1) tac <- t2 else tac <- rbind(tac,t2)
          SP <- split(tac,f=tac$Variable, drop=FALSE)
          str.colnames<- c("Word", "GlobalAverage", "AverageWord","Difer.", "pvalue")
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
      
      res <- fCharQuanti(object$Vocab$quanti$CharWord)
      print(res)

   #   for(i in 1:length(object$Vocab$quanti$CharWord)) {
    #    cat("\n\nCharacteristic words for each quantitative variable: ", names(object$Vocab$quanti$CharWord[i]), "\n")
    #    print(object$Vocab$quanti$CharWord$stats)
    #  }

    } # End quanti
    
    
  } 
  
 
  
 



if(file!=""){
  sink()
cat("\nAll the results are in file: ",file,"\n")
}
}
