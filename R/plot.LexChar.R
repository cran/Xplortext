#' @import gridExtra
#' @import stringr 
#' @importFrom utils modifyList
#' @export
plot.LexChar <- function (x, char.negat=TRUE, col.char.posit="blue", col.char.negat="red",
  col.lines="black", theme=theme_bw(), text.size=12,numr=1,numc=2, top=NULL, 
  max.posit=15, max.negat=15, type=c("CharWord","quanti","quali"), sel.var.cat="ALL", 
  txt.var.cat=NULL, sel.words="ALL",  ...)   
{
 ## Pendiente ayudas:
  # quanti, variables en orden alfabético
  # xlab
  
  dots <- list(...)
## Pendiente hjust
 # if("hjust" %in% names(dots)) hjust<- dots$hjust else hjust <- 0  ### <------------------------ Pendiente
 if("ylab" %in% names(dots)) ylab <- dots$ylab else ylab <- "vtest"
  # text.size

  
  if(is.null(type)) type <- "CharWord"
  
  options(stringsAsFactors = FALSE)
  if (!inherits(x, "LexChar"))  stop("x object should be LexChar class")
  words <-vtest <- NULL
  
  
  
  type <- match.arg(type[1], c("CharWord", "quanti", "quali"))
  if(max.negat<1) char.negat<-FALSE 
  if(char.negat==FALSE) max.negat<-0
  
    if(is.null(top)) top <- paste0("Characteristic words. Proba= ", x$Proba)

  theme$text$size <- text.size
  icont<-0
  
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
  
  
  # -------------------------------------------------------------------------------------
  
  fCharWord <- function(x, top)  {
    svc <- NULL
    icont<-0
    ldoc<-names(x)
    ntdoc<-length(ldoc)
    pword <-list()

    theme$text$size <- text.size
    
    if(length(sel.var.cat)==1) if(sel.var.cat=="ALL") sel.var.cat <- ldoc
    
    if(is.numeric(sel.var.cat)) 
      svc <- ldoc[sel.var.cat] else svc <- sel.var.cat[sel.var.cat %in% ldoc]
    # if(is.null(svc)) svc <- ldoc
    n.svc <- length(svc)
    
    for(i.svc in 1:n.svc)
      #   for (idoc in 1:ntdoc)
    {
      name.svc <- svc[i.svc]
      

      if(!is.null(x[[name.svc]])) {
        df <- data.frame(x[name.svc])
        df[,1] <- rownames(df)
        df[,2] <- df[,6] 
        df[,3:6] <- NULL
        rownames(df) <- NULL
        colnames(df) <- c("words", "vtest")
        df$words <- reorder(df$words,df$vtest)
       
      
        if(length(sel.words)==1) if(sel.words=="ALL") sel.words <- df$words
        if(is.null(sel.words)) sel.words <- df$words
        df <- df[df$words %in% sel.words,,drop=FALSE]
        
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
        
      
        # icont <- icont+1
        
    
        if(!is.null(txt.var.cat))  subtitle <- txt.var.cat[i.svc] else
        subtitle <- names(x[name.svc])
        
        pword[[i.svc]] <- ggplot(df) + geom_bar(aes(x=words,y=vtest),stat = "identity", color = col.lines,
                                                   fill = colorXX)+ coord_flip() +
          labs(title = subtitle)+ ylab(ylab) + xlab("") + 
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
  
  if(type=="quali")  {
    if(is.null(x$quali$CharWord)) stop("There are not qualitative variable to plot in LexChar object")
    if(is.null(top)) top <- paste0("Characteristic words. Proba= ", x$Proba)
    pword <- fCharWord(x$quali$CharWord,  top)
  } # End CharWord
  
  
  
  
  fCharWord.Quanti <- function(x, top, p.Proba)  {
    tac <- NULL
    pword <-list()
    str.colnames<- c("Word", "GlobalAverage", "AverageWord","Differ.", "vtest", "pvalue")
    if(is.null(top)) top <- paste0("Characteristic words. Proba= ", p.Proba)
    theme$text$size <- text.size
    
    # Variable selection
    if(length(sel.var.cat)==1) if(sel.var.cat=="ALL") sel.var.cat <- names.quanti
        if(is.numeric(sel.var.cat)) 
      svc <- names.quanti[sel.var.cat] else svc <- sel.var.cat[sel.var.cat %in% names.quanti]
    if(length(svc)==0) stop("There are no quantitative variables to select")
    n.svc <- length(svc)

    strcolnames<- c("GlobalAverage", "AverageWord","Difer.", "vtest", "pvalue", "Word", "Variable")

  

    for(i in 1:length(x)) {
      t1<- as.data.frame(x[i,drop=FALSE])
      t2 <- data.frame(t1,rep(names(x)[i],length(x[i])), rownames(t1))
      colnames(t2) <- strcolnames
      if(is.null(tac)) tac <- t2 else  tac <- rbind(tac,t2)
    }

    total.words <- unique(tac$Word)
    
    if(length(sel.words)==1) if(sel.words=="ALL") sel.words <- total.words
    sel.words <- sel.words[sel.words %in% total.words]
    tac <- tac[tac$Word %in% sel.words, , drop=FALSE]

    
    for(i.svc in 1:n.svc) {
      SP <- tac[tac$Variable==svc[i.svc],,drop=FALSE]
      t1 <- as.data.frame(tac[tac$Variable==svc[i.svc],,drop=FALSE])
     # t1 <- t1[t1$Word %in% sel.words,,drop=FALSE]
      t2.pos <- t1[t1[,3]>0, ,drop=FALSE]
      t2.neg <- t1[t1[,3]<0, ,drop=FALSE] 
      
      if(nrow(t2.pos)>0) {
        t2.pos <- t2.pos[order(-t2.pos[,4]),,drop=FALSE]
        rownames(t2.pos) <- paste0("P", c(1:nrow(t2.pos)))  
        t2.pos <- t2.pos[,c(6,4,5)]
        t2.pos <- t2.pos[1:min(max.posit,nrow(t2.pos)) ,,drop=FALSE]
      }
      if(nrow(t2.neg)>0) {
        t2.neg <- t2.neg[order(t2.neg[,4]),,drop=FALSE]
        rownames(t2.neg) <- paste0("N", c(1:nrow(t2.neg)))  
        t2.neg <- t2.neg[,c(6,4,5)]
        t2.neg <- t2.neg[1:min(max.negat,nrow(t2.neg)) ,,drop=FALSE]
        t2.neg <- t2.neg[order(-t2.neg[,2]),,drop=FALSE]
      } 
      
      numposit <- nrow(t2.pos)
      numnegat <- nrow(t2.neg)
      colorXX <- c(rep(col.char.posit,numposit),rep(col.char.negat,numnegat)) 

      df <- rbind(t2.pos, t2.neg)

      if(!is.null(txt.var.cat))  subtitle <- txt.var.cat[i.svc] else
        subtitle <- svc[i.svc]
      df$word <- factor(df$Word,                                    # Factor levels in decreasing order
                       levels = df$Word[order(df$vtest, decreasing = FALSE)])

      pword[[i.svc]] <- ggplot(df) + geom_bar(aes(x=word,y=vtest),stat = "identity", color = col.lines,
                                              fill = colorXX)+ coord_flip() + ylab(ylab) + xlab("") +
        labs(title = subtitle)+ 
        # ggtitle(title) + 
        theme(axis.text = element_text(size = text.size))+ theme 
    } # End for
    
    return(pword)

  } # End fCharWord.Quanti
  
  
  
  
  if(type=="quanti")  {
    if(is.null(x$quanti$CharWord)) stop("There are not quantitative variables to plot in LexChar object")
    if(is.null(top)) top <- paste0("Characteristic words. Proba= ", x$Proba)
    names.quanti <- rownames(x$quanti$stats)
    pword <- fCharWord.Quanti(x$quanti$CharWord, top, x$Proba)
    
  } # End CharWord
  
  
  

if(length(pword)==1) numr <- numc <- 1

 suppressWarnings(marrangeGrob2(grobs=pword, nrow = numr, ncol = numc, top=top))
}
