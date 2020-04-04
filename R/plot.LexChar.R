#' @import gridExtra
#' @importFrom utils modifyList
#' @export
plot.LexChar <- function (x, char.negat=TRUE, col.char.posit="blue", col.char.negat="red",
  col.lines="black", theme=theme_bw(), text.size=12,numr=1,numc=2, top=NULL, max.posit=15, max.negat=15,...) 
{
  options(stringsAsFactors = FALSE)
  
  marrangeGrob2<- function (grobs, ncol, nrow, ..., top ) 
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
 ldoc<-names(x$CharWord)
 ntdoc<-length(ldoc)
 pword <-list()
 proba<-x$Proba
if(is.null(top)) top <- paste0("Characteristic words. Proba= ",proba)
theme$text$size <- text.size

icont<-0
 for (idoc in 1:ntdoc)
    {
    if(!is.null(x$CharWord[[idoc]])) {
     df <- data.frame(x$CharWord[idoc])
     df[,1] <- rownames(df)
     df[,2] <- df[,6] 
     df[,3:6] <- NULL
     rownames(df) <- NULL
     colnames(df) <- c("words", "vtest")
     df$words <- reorder(df$words,df$vtest)
     if(!(char.negat)) df <- df[-df$vtest<0,,drop=FALSE]
     numposit<-nrow(df[df$vtest>0,])
     numnegat <- nrow(df[df$vtest<0,])

     if(numposit > max.posit) {
       df <- df[-((max.posit+1):numposit),,drop=FALSE]
       numposit <- max.posit
     }
     if(numnegat > max.negat) {
       df <- df[-((max.posit+1):(nrow(df)-max.negat)),,drop=FALSE]
     }

     colorXX <- c(rep(col.char.negat,nrow(df)- numposit),rep(col.char.posit,numposit)) 
     title <- names(x$CharWord[idoc]) 
icont <- icont+1

     pword[[icont]] <- ggplot(df) + geom_bar(aes(x=words,y=vtest),stat = "identity", color = col.lines, fill = colorXX)+ coord_flip() +
         ylab("vtest") + xlab("") + ggtitle(title) + theme
     }}
 suppressWarnings(marrangeGrob2(grobs=pword, nrow = numr, ncol = numc, top=top))
}
