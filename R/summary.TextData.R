###' @export
summary.TextData <- function(object,ndoc=10, nword=50, nseg=50, ordFreq = TRUE, file = NULL, sep=";",
      ...) 
{
  options(stringsAsFactors = FALSE)
 if (!inherits(object, "TextData")) 
 stop("Object should be TextData class")

if(is.null(ndoc)) ndoc<-0
if(is.null(nword)) nword<-0


if(ndoc=="ALL") ndoc <- nrow(object$summDoc)


if(nword=="ALL") nword <- object$info$Nword[[1]]
 nword <- min(nword, object$info$Nword[[1]])
 ndoc <- min(ndoc, nrow(object$summDoc))

 var.agg <- as.character(object$info$name.var.agg[[1]])

if(!is.null(object$DocSeg)) {
 if(nseg=="ALL") nseg <- ncol(object$DocSeg)
 if(is.null(nseg)) nseg<-0
 nseg <- min(nseg, ncol(object$DocSeg))
} else { nseg <-0}

if(!is.null(file)) sink(file)
cat("\nTextData summary\n\n")


if(is.null(file)) print(object$summGen)  else {
  out1 <- t(c("","Before", "After",sep="\t"))
  write.table(out1, quote=FALSE, sep = sep,col.names=FALSE, row.names=FALSE)
  write.table(object$summGen, quote=FALSE, sep = sep, col.names=FALSE)
}

if(var.agg=="") {
 name1 <- c("Occurrences", "DistinctWords", "PctLength", "Mean Length100")
 nameF <- c("DocName", name1, name1)
 name2 <- c("", rep("before",4), rep("after",4))
} else {
 nameF <- c("DocName", "Occurrences", "DistinctWords", "NumberDocs","PctLength", "MeanLength100", 
            "Occurrences", "DistinctWords", "PctLength", "MeanLength100")
 name2 <- c("", rep("before",5), rep("after",4))
}
A1 <- data.frame(nameF,name2)


# A3 7 para ordenar  
if(ordFreq==TRUE) {
  datDoc <- object$summDoc[with(object$summDoc, order(-object$summDoc$Occurrences.after)), ]
} else {datDoc <-object$summDoc  }

colnames(A1) <- NULL
A2 <- cbind(A1,t(datDoc[c(1:ndoc),]))


colnames(A2) <- NULL
A3 <- data.frame(t(A2)) 

rownames(A3) <- c(""," ",c(1:(nrow(A3)-2)))
colnames(A3) <- NULL


if(ndoc>0){
if(ndoc== object$info$Ndoc[[1]]) cat("\nStatistics for the documents\n")
 else cat("\nStatistics for the ", ndoc, " first docs\n")

  if(is.null(file))  print(A3[1:(ndoc+2),],na.print = "") else
   write.table(A3[1:(ndoc+2),], quote=FALSE, sep = sep)}

if(nword>0){
if (ordFreq) { 
if(nword== object$info$Nword[[1]]) cat("\nIndex of the words\n")
 else cat("\nIndex of the ", nword, " most frequent words\n")

   A1 <- cbind(format(as.data.frame(rownames(object$indexW[c(1:nword), ])), justify = "l"),
      object$indexW[c(1:nword), ])
      colnames(A1) <- c(format("Word",justify = "l"), "Frequency", "N.Documents")
 if(is.null(file))  print(A1) else {
   out1 <- t(c("","Word", "Frequency", "N.Documents",sep="\t"))
   write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
   write.table(A1, quote=FALSE, sep = sep, col.names=FALSE)}
 } else { 
     out1 <- t(c("","Word", "Frequency", "N.Documents",sep="\t"))
     Tabfq  <- object$indexW[(order(rownames(object$indexW)))[c(1:nword)],]

if(nword== object$info$Nword[[1]]) cat("\nIndex of the words in alphabetical order\n")
 else cat("\nIndex of the first ", nword, " words in alphabetical order\n")
     if(is.null(file)){  
    print(Tabfq) 
}else {
         cat("\nIndex of the words in alphabetical order\n")
         write.table(Tabfq, quote=FALSE, sep = sep, row.names=TRUE,col.names=FALSE)
    }
} 
}

if(nseg>0) { 
  cat("\nNumber of repeated segments\n")
  cat(" ", ncol(object$DocSeg),"\n") 
if (ordFreq) {
 df <- data.frame(object$indexS$segOrderFreq[,])
 df <- df[with(df, order(-frequency,segment)), ]
 row.has.na <- apply(df, 1, function(z){any(is.na(z))})
 A1 <- df[!row.has.na,]
 A1 <- A1[c(1:nseg),]
 nseg <- min(nseg, nrow(A1))
   colnames(A1) <- c(format("Segment",justify = "l"), "Frequency", "Long")
 if(nseg== ncol(object$DocSeg)) cat("\nIndex of the repeated segments\n")
 else cat("\nIndex of the ", nseg, " most frequent repeated segments\n")
  if(is.null(file))  print(A1) else {
     out1 <- t(c("Number","Segment", "Frequency", "Long",sep="\t"))
     write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
     write.table(object$indexS$segOrderFreq[1:nseg,], quote=FALSE, sep = sep, col.names=FALSE)}
} else {
  if(is.null(file)) {
 if(nseg== ncol(object$DocSeg)) cat("\nIndex of the repeated segments in alphabetical order\n")
  else cat("\nIndex of the first ", nseg, " repeated segments in alphabetical order\n")
    colnames(object$indexS$segOrderlist) <- c("Segment", "Frequency", "Long")
    print(object$indexS$segOrderlist[1:nseg,]) 
} else {
    cat("\nIndex of the repeated segments in alphabetical order\n")
    out1 <- t(c("Number","Segment", "Frequency", "Long",sep="\t"))
    write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
    write.table(object$indexS$segOrderlist[1:nseg,], quote=FALSE, sep = sep, col.names=FALSE)}
}} 


if(var.agg=="") {
 if(!is.null(object$context$quali))
  if(dim(object$context$quali)[[2]]>0)
{
  cat("\nSummary of the contextual categorical variables\n")
  print(summary(object$context$quali))}

if(!is.null(object$context$quanti))
if(ncol(object$context$quanti)>0){
 cat("\nSummary of the contextual quantitative variables\n")
  print(summary(object$context$quanti))}

} else {
  if(!is.null(object$context$quanti)){
  cat("\nSummary of the contextual quantitative variables\n")
 Tabfq <- cbind(object$context$quanti, object$summDoc[,4])
 sname <- colnames(object$context$quanti)
 colnames(Tabfq) <- c(sname,"Ndocs")
 print(Tabfq)}

if(!is.null(object$context$quali$qualitable))
{
  cat("\nSummary of the contextual categorical variables\n")
 cat("Supplementary tables in: $context$quali$qualitable\n")
}}
if(!is.null(file)) sink()
}


