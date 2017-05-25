#' @import ggplot2
#' @importFrom grDevices dev.new dev.interactive
#' @importFrom utils head write.table
#' @export
plot.TextData <- function (x, ndoc=25, nword=25, nseg=25, stop.word.tm=FALSE,
  stop.word.user=NULL, theme=theme_bw(), col.fill="grey", col.lines="black", text.size=12, ...)
{
if (!inherits(x, "TextData")) 
 stop("Object x should be TextData class")

idiom <- x$info$idiom[[1]]
theme$text$size <- text.size
words<-iFreq<-docnames<-rsegment<-frequency<-NULL


if(nseg>0){
if(!is.null(x$DocSeg)==TRUE) {
 titleseg <- paste0( nword," most frequent segments" )
 df <- data.frame(x$indexS$segOrderFreq[1:nseg,])
 colnames(df)[1] <- "rsegment"
 df$rsegment <- reorder(df$rsegment,df$frequency)
 row.has.na <- apply(df, 1, function(z){any(is.na(z))})
 df <- df[!row.has.na,]
 pseg <- ggplot(df)+ geom_bar(aes(x=rsegment, y=frequency),stat = "identity", color = col.lines, fill = col.fill) + coord_flip() + ylab("Repeated segment frequency") + xlab("") + ggtitle(titleseg) + theme
if (dev.interactive()) dev.new()
 print(pseg)
}}


df <- data.frame(x$indexW[,1], rownames(x$indexW))
colnames(df) <- c("iFreq", "words") 

# Remove stopwords
if(stop.word.tm==TRUE) {
   stop.word <- tm::stopwords(idiom)
   pos <- which(!df$words%in%stop.word)
   df <-  df[pos,]}

# Remove user words
if(!is.null(stop.word.user)) {
   pos <- which(!df$words%in%stop.word.user)
   df <-  df[pos,]}
nword <- min(nword,nrow(df))
if(nword>0){
df <- df[c(1:nword),]
title <- paste0( nword," most frequent words" )
# Present the results ordered by frequency
df$words <- reorder(rownames(df),df$iFreq)
pword <- ggplot(df)+ geom_bar(aes(x=words, y=iFreq),stat = "identity", color = col.lines, fill = col.fill)+ coord_flip() +
ylab("Word frequency") + xlab("") + ggtitle(title) + theme
if (dev.interactive()) dev.new()
print(pword)
}



# To build a data frame with three columns (rows, columns and frequencies) keeping compress structure
if(ndoc>0){
DT <- data.frame(x$DocTerm$i,x$DocTerm$j,x$DocTerm$v)
docnames <- rownames(x$DocTerm) 
colnames(DT) <- c("i","j","v")
iFreq <- sapply(split(DT, DT$i), function(z) sum(z$v))
df <- data.frame(docnames, iFreq)
# Presents the results ordered by frequency
df <- df[with(df, order(-iFreq,docnames)), ]
ndoc <- min(ndoc, nrow(df))
df <- df[c(1:ndoc),]
titledoc <- paste0(ndoc," documents with higher length" )
df$docnames <- reorder(df$docnames,df$iFreq)
pdoc <-ggplot2::ggplot(df)+ geom_bar(aes(x=docnames, y=iFreq),stat = "identity", color = col.lines, fill = col.fill)+ coord_flip() +
  ylab("Document length") + xlab("") + ggtitle(titledoc) + theme
if (dev.interactive()) dev.new()
print(pdoc)
}
}
