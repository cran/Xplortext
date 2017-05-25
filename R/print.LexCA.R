#' @export
print.LexCA <- function (x, file = NULL, sep = ";", ...) 
{

sink.reset <- function(){
    for(i in seq_len(sink.number())){
        sink()}}
sink.reset
res.LexCA <- x

	if (!inherits(res.LexCA, "LexCA")) stop("x object should be LexCA class")
	cat("**Results for CA and Aggregate Lexical Table (LexCA)**\n")
    	cat("*The results are available in the following objects:\n\n")
       indice <- 12 
    	res <- array("", c(12, 2), list(1:12, c("name", "description")))
   	res[1, ] <- c("$eig", "Eigenvalues and % of variance")
   	res[2, ] <- c("$row", "CA results for the active documents/aggregate documents")
   	res[3, ] <- c("$col", "CA results for active words")
   	res[4, ] <- c("$row.sup", "CA results for the supplementary documents/aggregate documents")
   	res[5, ] <- c("$col.sup", "CA results for supplementary words")
   	res[6, ] <- c("$quanti.sup", "CA results for the supplementary continuous variables")
   	res[7, ] <- c("$quali.sup", "CA results for the supplementary categorical variables")
   	res[8, ] <- c("$meta", "list of Keydocs selected documents and Metakeys selected words")
   	res[9, ] <- c("$VCr", "Cramer's V coefficient")
   	res[10, ] <- c("Inertia", "total inertia")
   	res[11, ] <- c("segment", "CA results for repeated segments")
   	res[12, ] <- c("var.agg","name of the aggregate variable used")


    print(res[1:indice, ])
    	if (!is.null(file)) {
 sink(file)

  cat("\nEigenvalues and percentage of explained variance\n")
  out1 <- t(c("",colnames(x$eig)))
  write.table(out1, quote=FALSE, sep = sep,col.names=FALSE, row.names=FALSE)
  x$eig <- format(x$eig, digits=6, scientific = FALSE)
  write.table(x$eig, quote=FALSE, sep = sep, col.names=FALSE)

if(is.null(x$var.agg)) {
  msg1 <- " for active documents\n" 
  msg4 <- " for supplementary documents\n"
} else{  
  msg1 <- " for active aggregate documents\n"
  msg4 <- " for supplementary aggregate documents\n"
 }
 msg3 <- " for active words\n"
 cat("\n\nCA results", msg1)
 cat("\nCoordinates ", msg1)
 msg2 <- t(c("",colnames(x$row$coord)))
 write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
 write.table(x$row$coord, quote=FALSE, sep = sep, col.names=FALSE)

 cat("\nContributions ", msg1)
 write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
 write.table(x$row$contrib, quote=FALSE, sep = sep, col.names=FALSE)

 cat("\nSquared Cosinus ", msg1)
 write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
 write.table(x$row$contrib, quote=FALSE, sep = sep, col.names=FALSE)

 cat("\nInertia ", msg1)
 x$row$inertia <- data.frame(x$row$inertia)
 rownames(x$row$inertia) <- rownames(x$row$coord) 
 write.table(x$row$inertia, quote=FALSE, sep = sep,row.names=TRUE, col.names=FALSE)

 cat("\nCoordinates ", msg3)
 msg2 <- t(c("",colnames(x$col$coord)))
 write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
 write.table(x$col$coord, quote=FALSE, sep = sep, col.names=FALSE)

 cat("\nContributions ", msg3)
 write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
 write.table(x$col$contrib, quote=FALSE, sep = sep, col.names=FALSE)

 cat("\nSquared Cosinus ", msg3)
 write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
 write.table(x$col$contrib, quote=FALSE, sep = sep, col.names=FALSE)

 cat("\nInertia ", msg3)
 x$col$inertia <- data.frame(x$col$inertia)
 rownames(x$col$inertia) <- rownames(x$col$coord) 
 write.table(x$col$inertia, quote=FALSE, sep = sep,row.names=TRUE, col.names=FALSE)

# row.sup
if(!is.null(x$row.sup)) {
 cat("\nCoordinates of ", msg4)
 write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
 write.table(x$row.sup$coord, quote=FALSE, sep = sep, col.names=FALSE)

 cat("\nSquared Cosinus ", msg4)
 write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
 write.table(x$row.sup$cos2,quote=FALSE, sep = sep, col.names=FALSE)
}

# col.sup
if(!is.null(x$col.sup)) {
 ncsup <- nrow(x$col.sup$coord)

 if(!is.null(x$segment)) ncsup <- ncsup - nrow(x$segment$coord)
  cat("\nCoordinates of the supplementary words\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$col.sup$coord[(1:ncsup),], quote=FALSE, sep = sep, col.names=FALSE)

  cat("\nSquared Cosinus of the supplementary words\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$col.sup$cos2[1:ncsup,],quote=FALSE, sep = sep, col.names=FALSE)
}

# quanti.sup
if(!is.null(x$quanti.sup)){
  cat("\nCoordinates of the supplementary quantitative variables\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$quanti.sup$coord, quote=FALSE, sep = sep, col.names=FALSE)

  cat("\nSquared Cosinus of quantitative variables\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$quanti.sup$cos2,quote=FALSE, sep = sep, col.names=FALSE)
}

# quali.sup
if(is.null(x$var.agg)){
if(!is.null(x$quali.sup)){
  cat("\nCoordinates of the supplementary categories\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$quali.sup$coord, quote=FALSE, sep = sep, col.names=FALSE)

  cat("\nv-test of the supplementary categories\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$quali.sup$v.test,quote=FALSE, sep = sep, col.names=FALSE)

  cat("\nEta Squared of the supplementary categories\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$quali.sup$eta2,quote=FALSE, sep = sep, col.names=FALSE)
}}



# keys
nEig <- max(x$meta$Doc$Dim, x$meta$Word$Dim)

if(!is.null(x$meta)){
 cat("\nlist of Keydocs documents\n")
  msg5 <- t(c("","Docs", "Dim", "Coordinate"))
  write.table(msg5, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$meta$Doc, quote=FALSE, sep = sep,row.names=TRUE, col.names=FALSE)

cat("\n\nDocuments whose contribution is over ", x$meta$lmd, " times the average document contribution, axis by axis\n")
for(i in 1:nEig) {
 cat("\nDimension ", i, "+\n")
 cat(x$meta$Doc$Docs[which(x$meta$Doc$Dim==i & x$meta$Doc$Coordinates>0)])
 cat("\nDimension ", i, "-\n")
 cat(x$meta$Doc$Docs[which(x$meta$Doc$Dim==i & x$meta$Doc$Coordinates<0)],"\n")
}

dfD <- x$meta$Doc
dfD$Coordinates <- ifelse(dfD$Coordinates>0, as.numeric(dfD$Dim), -1* as.numeric(dfD$Dim))
TDsign <- table(dfD$Docs, dfD$Coordinates)
colnames(TDsign) <- paste("DIM_",colnames(TDsign))
TDsign <- addmargins(TDsign)
  cat("\n\nDocuments whose contribution is over ", x$meta$lmd, " times the average document contribution, axis by axis\n")
  write.table(t(c("",colnames(TDsign))), quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(TDsign, quote=FALSE, sep = sep,row.names=TRUE, col.names=FALSE)

  cat("\nlist of Metakeys selected words\n")
  msg6 <- t(c("","Words", "Dim", "Coordinate"))
  write.table(msg6, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$meta$Word, quote=FALSE, sep = sep,row.names=TRUE, col.names=FALSE)


cat("\n\nWords whose contribution is over ", x$meta$lmw , " times the average word contribution, axis by axis\n")
for(i in 1:nEig) {
 cat("\nDimension ", i, "+\n")
 cat(x$meta$Word$Words[which(x$meta$Word$Dim==i & x$meta$Word$Coordinates>0)])
 cat("\nDimension ", i, "-\n")
 cat(x$meta$Word$Words[which(x$meta$Word$Dim==i & x$meta$Word$Coordinates<0)],"\n")
}
dfW <- x$meta$Word
dfW$Coordinates <- ifelse(dfW$Coordinates>0, as.numeric(dfW$Dim), -1* as.numeric(dfW$Dim))
TDsign <- table(dfW$Words, dfW$Coordinates)
colnames(TDsign) <- paste("DIM_",colnames(TDsign))
TDsign <- addmargins(TDsign)
  cat("\n\nWords whose contribution is over ", x$meta$lmw, " times the average word contribution, axis by axis\n")
  write.table(t(c("",colnames(TDsign))), quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(TDsign, quote=FALSE, sep = sep,row.names=TRUE, col.names=FALSE)
}

# VCr
 cat("\n\nCramer's V coefficient\n")
 cat(x$VCr)

# Inertia
 cat("\n\nInertia\n")
 cat(x$Inertia)

# segment
if(!is.null(x$segment$coord)){
 cat("\nRepeated segments\n")
 cat("\nCoordinates of the repeated segments\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$segment$coord, quote=FALSE, sep = sep, col.names=FALSE)

  cat("\nSquared Cosinus of of the repeated segments\n")
  write.table(msg2, quote=FALSE, sep=sep, col.names=FALSE, row.names=FALSE)
  write.table(x$segment$cos2, quote=FALSE, sep = sep, col.names=FALSE)
}

 sink()}
       	print(paste("All the results are in file", file))
    }
