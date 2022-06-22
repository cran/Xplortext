#' @export
print.TextData <- function (x, file = NULL, sep = ";", ...) 
{
  options(stringsAsFactors = FALSE)
  
    res.TextData <- x
  if (!inherits(res.TextData, "TextData")) 
        stop("non convenient data")

sink.reset <- function(){
    for(i in seq_len(sink.number())){sink()}}
sink.reset

    cat("*The results are available in the following objects:\n\n")
addindex <- 0 
if(is.null(x$DocSeg)) indexp <- 7 else indexp <- 9 
if(!is.null(x$SourceTerm)) addindex <- 1 
if(!is.null(x$FullDocTerm)) addindex <- addindex +1
index <- indexp+addindex
    res <- array("", c(index, 2), list(1:index, c("name", "description")))
    res[1, ] <- c("$summGen", "General summary")
    res[2, ] <- c("$summDoc", "Documents summary")
    res[3, ] <- c("$indexW", "Index of words")
    res[4, ] <- c("$DocTerm", "Documents by Words table (*1)")
    res[5, ] <- c("$context", "Contextual variables")
    res[6, ] <- c("$info", "Information about selection of words")
    res[7, ] <- c("$remov.docs", "Removed empty documents")
if(!is.null(x$DocSeg)){
     res[8, ] <- c("$indexS", "Index of segments")
     res[9, ] <- c("$DocSeg", "Documents by Segments table (*2)")}
if(!is.null(x$SourceTerm)) {
    res[indexp+1, ] <- c("$SourceTerm", "Sequency Documents by Words table (*3)")
}

    print(res[c(1:index),])
    print(paste("(*1) $DocTerm has a compressed format"))
    print(paste("Use print(as.matrix(x$DocTerm))"))
 
if(!is.null(x$DocSeg)) {
    print(paste("(*2) $DocSeg has a compressed format"))
    print(paste("Use print(as.matrix(x$DocSeg))"))
}

if(!is.null(x$SourceTerm)) {
    print(paste("(*3) $SourceTerm has a compressed format"))
    print(paste("Use print(as.matrix(x$SourceTerm))"))
}

if(!is.null(x$FullDocTerm)) {
    print(paste("(*4) $FullDocTerm has a compressed format"))
    print(paste("Use print(as.matrix(x$FullDocTerm))"))
}

    if (!is.null(file)) {
for(i in seq_len(sink.number())){sink(NULL)}
 sink(file)
 ndoc <- x$info$Ndoc[[1]]
 var.agg <- as.character(x$info$name.var.agg[[1]])

cat("\nGeneral summary (summGen)\n")
if(is.null(file)) print(x$summGen)  else {
  out1 <- t(c("","Before", "After",sep="\t"))
  write.table(out1, quote=FALSE, sep = sep,col.names=FALSE, row.names=FALSE)
  write.table(x$summGen, quote=FALSE, sep = sep, col.names=FALSE)}

cat("\nDocuments summary (summDoc)\n")  
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
colnames(A1) <- NULL
A2 <- cbind(A1,t(x$summDoc))
colnames(A2) <- NULL
A3 <- data.frame(t(A2)) 
colnames(A3) <- NULL
rownames(A3) <- c(""," ",c(1:(nrow(A3)-2)))
write.table(A3[1:(ndoc+2),], quote=FALSE, sep = sep)

cat("\nIndex of words (indexW)\n")
   cat("\nIndex of the most frequent words\n")
   A1 <- cbind(format(as.data.frame(rownames(x$indexW)), justify = "l"), x$indexW)
   colnames(A1) <- c(format("Word",justify = "l"), "Frequency", "N.Documents")
   out1 <- t(c("","Word", "Frequency", "N.Documents",sep="\t"))
   write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
   write.table(A1, quote=FALSE, sep = sep, col.names=FALSE)
   out1 <- t(c("","Word", "Frequency", "N.Documents",sep="\t"))
   Tabfq  <- x$indexW[order(rownames(x$indexW)),]
   cat("\nIndex of the words in alphabetical order\n")
   write.table(Tabfq, quote=FALSE, sep = sep, col.names=FALSE)

cat("\nDocuments by Words table (DocTerm)\n")
 out1 <- t(c("", colnames(x$DocTerm),sep="\t"))
 write.table(out1, quote=FALSE, sep = sep,col.names=FALSE, row.names=FALSE)
 write.table(as.matrix(x$DocTerm), quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE)

# Repeated segments
if(!is.null(x$DocSeg)) { 
 cat("\nNumber of repeated segments\n")
 cat(" ", ncol(x$DocSeg),"\n")
 cat("\nIndex of segments (indexS)\n")
 cat("\nSegments ordered by frequency (indexS$segOrderFreq)\n")
  A1 <- cbind(format(as.data.frame(rownames(x$indexS$segOrderFreq)), justify = "l"),
     x$indexS$segOrderFreq)
     colnames(A1) <- c("Number", format("Segment",justify = "l"), "Frequency", "Long")
     out1 <- t(c("Number","Segment", "Frequency", "Long",sep="\t"))
     write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
     write.table(x$indexS$segOrderFreq, quote=FALSE, sep = sep, col.names=FALSE)
 cat("\nSegments in alphabetical order (indexS$segOrderlist)\n")
 out1 <- t(c("Number","Segment", "Frequency", "Long",sep="\t"))
    write.table(out1, quote=FALSE, sep = sep, col.names=FALSE, row.names=FALSE)
    write.table(x$indexS$segOrderlist, quote=FALSE, sep = sep, col.names=FALSE)
 cat("\nDocuments by Repeated segments table (DocSeg)\n")
   out1 <- t(c("", colnames(x$DocSeg),sep="\t"))
  write.table(out1, quote=FALSE, sep = sep,col.names=FALSE, row.names=FALSE)
  write.table(as.matrix(x$DocSeg), quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE)
}

# Supplementary variables
if(var.agg=="") {
 if(!is.null(x$context$quali))
 if(length(x$context$quali)>0) {
  cat("\nSummary of contextual qualitative variables\n")
  print(summary(x$context$quali))}

if(!is.null(x$context$quanti))
if(ncol(x$context$quanti)>0){
 cat("\nSummary of contextual quantitative variables\n")
  print(summary(x$context$quanti))}
} else {
  if(!is.null(x$context$quanti)){
  cat("\nSummary of contextual quantitative variables\n")
 Tabfq <- cbind(x$context$quanti, x$summDoc[,4])
 sname <- colnames(x$context$quanti)
 colnames(Tabfq) <- c(sname,"Ndocs")
  out1 <- t(c("",colnames(Tabfq)))
  write.table(out1, quote=FALSE, sep = sep,col.names=FALSE, row.names=FALSE)
  write.table(Tabfq, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE)
}
if(!is.null(x$context$quali$qualivar))
 if(ncol(x$context$quali$qualivar)>0){
  cat("\nSummary of contextual qualitative variables\n")
  out1 <- t(c("",colnames(x$context$quali$qualitable)))
  write.table(out1, quote=FALSE, sep = sep,col.names=FALSE, row.names=FALSE)
  write.table(x$context$quali$qualitable, quote=FALSE, sep = sep, col.names=FALSE, row.names=TRUE)
 }
}

if(!is.null(x$SourceTerm)) {
  cat("\nFirst 5 documents in the source documents by words table\n")
  cat("For print all documents to a file use  write.table(as.matrix(x$SourceTerm),file=''...)\n")
  write.table(head(as.matrix(x$SourceTerm)))}

 cat("\nInformation about selected options in TextData\n")
 writeLines(unlist(lapply(x$info, paste, collapse=" ")))
  sink()	
}

if (!is.null(file)) {	
        	print(paste("All the results are in the file", file))
    	}
}
