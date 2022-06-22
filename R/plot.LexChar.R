#' @import gridExtra
#' @import stringr 
#' @importFrom utils modifyList
#' @export
plot.LexChar <- function (x, char.negat=TRUE, col.char.posit="blue", col.char.negat="red",
  col.lines="black", theme=theme_bw(), text.size=12,numr=1,numc=2, top=NULL, max.posit=15, max.negat=15, 
 # type=c("CharWord","quanti","quali"),context=NULL, ...) 
 type=c("CharWord","quanti","quali"), ...) 
   # 1. Eliminar de la ayuda context, controlar de las versiones anteriores
   
{
   
   
  options(stringsAsFactors = FALSE)
   varnames<- lapply(substitute(list(...))[-1], deparse)
 #  if(!is.null(varnames$context.sup)) {
 #     context <- gsub("[[:punct:]]", "", varnames$context.sup)   # Only for Compatibility version 1.4.1
 #     warning("Versions > 1.4.1 use context, no context.sup argument") 
 #  }
   
  type <- match.arg(type[1], c("CharWord", "quanti", "quali"))
 if(max.negat<1) char.negat<-FALSE 
  if(char.negat==FALSE) max.negat<-0
 # context.sup <- context

  
  marrangeGrob2<- function(grobs, ncol, nrow, ..., top ) 
{
    n <- length(grobs)
    nlay <- nrow * ncol
    pages <- n%/%nlay + as.logical(n%%nlay)
    groups <- split(seq_along(grobs), gl(pages, nlay, n))
    pl <- vector(mode = "list", length = pages)

    for (page in seq_along(groups)) {
        g <- page
        params <- modifyList(list(...), list(top = eval(top), 
            nrow = nrow, ncol = ncol))
        pl[[g]] <- do.call(gridExtra::arrangeGrob, c(grobs[groups[[g]]],params))
    }
    class(pl) <- c("arrangelist", class(pl))
    pl
}

 if (!inherits(x, "LexChar"))  stop("x object should be LexChar class")


 words <-vtest <- NULL
 

 
 fCharWord <- function(x, top)  {
   icont<-0
   ldoc<-names(x)
   ntdoc<-length(ldoc)
   pword <-list()
 #  proba<-Proba

   theme$text$size <- text.size
   for (idoc in 1:ntdoc)
   {
     if(!is.null(x[[idoc]])) {
       df <- data.frame(x[idoc])
       df[,1] <- rownames(df)
       df[,2] <- df[,6] 
       df[,3:6] <- NULL
       rownames(df) <- NULL
       colnames(df) <- c("words", "vtest")
       df$words <- reorder(df$words,df$vtest)
       
       if(!char.negat)      df <- df[-df$vtest<0,,drop=FALSE]
       numposit<-nrow(df[df$vtest>0,])
       numnegat <- nrow(df[df$vtest<0,])
       
       if(numposit > max.posit) {
         df <- df[-((max.posit+1):numposit),,drop=FALSE]
         numposit <- max.posit
       }

       if(numnegat > max.negat) {
         df <- df[-((numposit+1):(nrow(df)-max.negat)),,drop=FALSE]
         numnegat <- max.negat
       }
       
       colorXX <- c(rep(col.char.posit,numposit),rep(col.char.negat,nrow(df)- numposit)) 
       
       subtitle <- names(x[idoc]) 
       icont <- icont+1
       
       pword[[icont]] <- ggplot(df) + geom_bar(aes(x=words,y=vtest),stat = "identity", color = col.lines,
                                               fill = colorXX)+ coord_flip() +
         labs(title = subtitle)+ ylab("vtest") + xlab("") + 
         # ggtitle(title) + 
         theme(axis.text = element_text(size = text.size))+ theme

     }
   }
   return(pword)
 }

 
 
 
 
if(type=="CharWord")  {
  if(is.null(top)) top <- paste0("Characteristic words. Proba= ", x$Proba)
 pword <- fCharWord(x$CharWord,  top)
} # End CharWord
 

 
 fChar<- function(x)  {
   tac <- NULL
   strcolnames<- c("GlobalAverage", "AverageWord","Difer.", "pvalue", "Word", "Variable")
   for(i in 1:length(x)) {
     t1<- as.data.frame(x[i,drop=FALSE])
     t2 <- data.frame(t1,rep(names(x)[i],length(x[i])), rownames(t1))
     colnames(t2) <- strcolnames
     if(is.null(tac)) tac <- t2 else  tac <- rbind(tac,t2)
    }
     
     SP <- split(tac,f=tac$Variable, drop=FALSE)
     str.colnames<- c("Word", "GlobalAverage", "AverageWord","Differ.", "pvalue")
     empty_list = structure(vector(mode = "list", length = length(SP)), names = names(SP))			
     for(i in 1:length(SP)) {
       t1<- as.data.frame(SP[i,drop=FALSE])
       t2.pos <- t1[t1[,3]>0, ,drop=FALSE]
       t2.neg <- t1[t1[,3]<0, ,drop=FALSE]

         if(nrow(t2.pos)>0) {
         t2.pos <- t2.pos[order(-t2.pos[,3]),,drop=FALSE]
         rownames(t2.pos) <- paste0("P", c(1:nrow(t2.pos)))
         t3.pos <- t2.pos[,c(5,1:4)]
         colnames(t3.pos) <- str.colnames
         empty_list[[i]]$posit <- t3.pos
                             }
       if(nrow(t2.neg)>0) {
         t2.neg <- t2.neg[order(t2.neg[,3]),,drop=FALSE]
         rownames(t2.neg) <- paste0("N", c(1:nrow(t2.neg)))
         t3.neg <- t2.neg[,c(5,1:4)]
         colnames(t3.neg) <- str.colnames
         empty_list[[i]]$negat <- t3.neg
       }
     } # End for
  
   return(empty_list)
 }
 
 if(type=="quanti")  {
   if(is.null(x$Vocab$quanti$CharWord)) stop("There is not quantitative variables to plot in LexChart object")
   res <- fChar(x$Vocab$quanti$CharWord)
   ldoc<-names(res)
   ntdoc<-length(ldoc)
   pword <-list()
   proba<-x$Proba

   if(is.null(top)) top <- paste0("Characteristic words. Proba= ",proba)
   theme$text$size <- text.size
   icont<-0


   for (idoc in 1:ntdoc)
   {
     if(!is.null(res[[idoc]])) {
       t1<- res[idoc,drop=FALSE][[1]]

       if(!is.null(t1$posit)) numposit<-nrow(t1$posit) else numposit <- 0
       if(!is.null(t1$negat)) numnegat<-nrow(t1$negat) else numnegat <- 0
       
        df<-NULL
       if(numposit > max.posit) {
         t1$posit <- t1$posit[-((max.posit+1):numposit),,drop=FALSE]
         numposit <- max.posit
        }

       if(numnegat > max.negat) {
         t1$negat <- t1$negat[-((max.negat+1):numnegat),,drop=FALSE]
         
         numnegat <- max.negat
       }
  if(!is.null(t1$posit)) df <- rbind(df,t1$posit)
  if(!is.null(t1$negat)) df <- rbind(df,t1$negat)
     #  colorXX <- c(rep(col.char.negat,numnegat),rep(col.char.posit,numposit)) 
       colorXX <- c(rep(col.char.posit,numposit),rep(col.char.negat,nrow(df)- numposit))     
      # title <- names(res[idoc])

       df <- df[order(df[,"Differ."]),c("Word", "Differ."),drop=FALSE]     
       colnames(df) <- c("Word", "Difference")
       df$Word <- reorder(df$Word,df$Difference)
       
       subtitle <- ldoc[idoc]
       icont  <- icont+1

    pword[[icont]] <- ggplot(df) + geom_bar(aes(x=df$Word,y=df$Difference),stat = "identity", 
                                               color = col.lines, fill = colorXX)+ coord_flip() +
        labs(title = subtitle)+ 
           ylab("Difference over the average") + xlab("") +
        theme(axis.text = element_text(size = text.size))+ theme
     } # End if
   
     } # End for
# return(pword)
 }

 
 if(type=="quali") {
   if(is.null(x$Vocab$quali$CharWord)) stop("There are not qualitative variables in LexChart object")
         pword <- fCharWord(x$Vocab$quali$CharWord,  top)
 } # End Type quali

  

 suppressWarnings(marrangeGrob2(grobs=pword, nrow = numr, ncol = numc, top=top))
}
