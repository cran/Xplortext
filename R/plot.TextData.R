#' @import ggplot2
#' @importFrom grDevices dev.new dev.interactive
#' @importFrom utils head write.table
#' @export
plot.TextData <- function (x, ndoc=25, nword=25, nseg=25, sel=NULL,
  stop.word.tm=FALSE, stop.word.user=NULL, theme=theme_bw(), title=NULL, xtitle=NULL,
  col.fill="grey", col.lines="black", text.size=12, ...)
{
if (!inherits(x, "TextData")) 
 stop("Object x should be TextData class")

idiom <- x$info$idiom[[1]]
theme$text$size <- text.size
words<-iFreq<-docnames<-rsegment<-frequency<-NULL

if(!is.null(sel)) {
if(sel=="seg") {ndoc<-0; nword<-0} 
if(sel=="word") {ndoc<-0; nseg<-0}
if(sel=="doc") {nword<-0; nseg<-0}
}

if(nseg>0){
if(!is.null(x$DocSeg)==TRUE) {
ifelse(is.null(title), ttitle <- paste0( nseg," most frequent segments" ), ttitle <- title)
ifelse(is.null(xtitle), txtitle <- "Repeated segment frequency" ,txtitle <- xtitle)

df <- data.frame(x$indexS$segOrderFreq[,])
 colnames(df)[1] <- "rsegment"
 df <- df[with(df, order(-frequency,rsegment)), ]
 row.has.na <- apply(df, 1, function(z){any(is.na(z))})
 df <- df[!row.has.na,]
 nseg <- min(nseg, nrow(df))
 df <- df[c(1:nseg),]
 rownames(df) <- c(1:nseg)
 df$rsegment <- reorder(df$rsegment,-as.numeric(rownames(df)))

 pseg <- ggplot(df)+ geom_bar(aes(x=rsegment, y=frequency),stat = "identity", color = col.lines, fill = col.fill) + coord_flip() +
 ylab(txtitle) + xlab("") + ggtitle(ttitle) + theme + theme(plot.title = element_text(hjust = 0.5))
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
ifelse(is.null(title), ttitle <- paste0( nword," most frequent words"), ttitle <- title)
ifelse(is.null(xtitle), txtitle <- "Word frequency", txtitle <- xtitle)

# Present the results ordered by frequency
 df <- df[with(df, order(-iFreq,words)), ]
 df <- df[c(1:nword),]
 rownames(df) <- c(1:nword)
 df$words <- reorder(df$words,-as.numeric(rownames(df)))
 pword <- ggplot(df)+ geom_bar(aes(x=words, y=iFreq),stat = "identity", color = col.lines, fill = col.fill)+ coord_flip() +
ylab(txtitle) + xlab("") + ggtitle(ttitle) + theme + theme(plot.title = element_text(hjust = 0.5))
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

ifelse(is.null(title), ttitle <- paste0( ndoc," documents with higher length"), ttitle <- title)
ifelse(is.null(xtitle), txtitle <- "Document length", txtitle <- xtitle)
df$docnames <- reorder(df$docnames,df$iFreq)
pdoc <-ggplot2::ggplot(df)+ geom_bar(aes(x=docnames, y=iFreq),stat = "identity", color = col.lines, fill = col.fill)+ coord_flip() +
  ylab(txtitle) + xlab("") + ggtitle(ttitle) + theme + theme(plot.title = element_text(hjust = 0.5))
if (dev.interactive()) dev.new()
print(pdoc)
}

theme_set(theme)
}
