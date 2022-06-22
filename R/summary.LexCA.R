#' @export
summary.LexCA <- function (object, ncp=5, nb.dec = 3, ndoc=10, nword=10, nseg=10, 
  nsup=10, metaDocs=FALSE, metaWords=FALSE, file = NULL, ...) 
{
if (!inherits(object, "LexCA")) stop("non convenient object")
  options(stringsAsFactors = FALSE)
  
printE <- function(mat, file = "") {
            mat <- cbind(format(rownames(mat)), mat)
            mat <- rbind(colnames(mat), mat)
            mat2 <- cbind(format(mat[, 1], justify = "right"),format(mat[, 2], justify = "right"))
            for (k in 3:ncol(mat)) mat2 <- cbind(mat2, format(mat[,k], justify = "right"))
            mat2 <- cbind(mat2, "\n")
            for (i in 1:nrow(mat2)) cat(mat2[i, ])} 

printC <- function(obj, text1, nval, nDim, nIni) {
 cat(paste("\n",text1,"\n"))
 if(nrow(obj)> nval) cat(paste("Only the first ", (nval-nIni+1), "elements are shown\n"))
  print(round(obj[nIni:(nval+nIni-1),1:nDim,drop=FALSE],nb.dec), quote=FALSE)
}


printI <- function(obj, obj2, text1, nval) {
# cat(paste("\n",text1,"\n"))
 cat(paste("\n"))
 obj <- matrix(round(obj[1:nval],nb.dec),nrow=nval,ncol=1)
 rownames(obj) <- rownames(obj2)[1:nrow(obj)]; colnames(obj) <- "Inertia"
# if(nrow(obj2)> nval) cat(paste("Only the first ", nval, "elements are shown\n"))
 print(obj)
}

if(is.null(ndoc)) ndoc<-0
if(is.null(nword)) nword<-0
if(is.null(nseg)) nseg<-0
if(is.null(nsup)) nsup<-0
if(is.null(object$segment$coord)) nseg<-0



if(ndoc=="ALL") ndoc<- nrow(object$row$coord)
if(nword=="ALL") nword<- nrow(object$col$coord)
if(nseg=="ALL") nseg<- nrow(object$segment$coord)

ndoc <- min(nrow(object$row$coord),ndoc)
nword <- min(nrow(object$col$coord),nword)
nseg <- min(nseg,nrow(object$segment$coord))

# for(i in seq_len(sink.number())){sink(NULL)}
if(!is.null(file)) sink(file) else file=""

cat("Correspondence analysis summary\n")
cat("\nEigenvalues\n")

nEig <- min(ncp,ncol(object$row$coord))
eige <- format(t(round(object$eig[1:nEig, 1:3], nb.dec)), justify = "right")
   rownames(eige) <- c("Variance", "% of var.", "Cumulative % of var.")
    nEig <- min(nEig, nrow(object$eig))
    colnames(eige[, 1:nEig]) <- paste("Dim", 1:nEig, sep = ".")
    printE(t(eige[, 1:nEig]))
    cat(paste("\nCramer's V ", round(object$VCr,nb.dec) ,"   Inertia ", round(object$Inertia, nb.dec), "\n", sep=" "))

if(ndoc>0) {
cat("\n\nDOCUMENTS")
if(!is.null(object$var.agg)) cat("\nAll documents are aggregate documents\n")
printC(object$row$coord,"Coordinates", ndoc, nEig,1)
printC(object$row$contrib,"Contributions (by column total=100)", ndoc, nEig,1)
printC(object$row$cos2,"Square cosinus (by row total=1)", ndoc, nEig,1)
printI(object$row$inertia,object$row$coord, "", ndoc)}

if(metaDocs!=FALSE){
cat("\n\nDocuments whose contribution is over ", object$meta$lmd, " times the average document contribution\n")
for(i in 1:nEig) {
 cat("\nDimension ", i, "+\n")
 cat(object$meta$Doc$Docs[which(object$meta$Doc$Dim==i & object$meta$Doc$Coordinates>0)])
 cat("\nDimension ", i, "-\n")
 cat(object$meta$Doc$Docs[which(object$meta$Doc$Dim==i & object$meta$Doc$Coordinates<0)],"\n")
}}



if(nword>0) {
cat("\n\nWORDS\n")
printC(object$col$coord,"Coordinates", nword, nEig,1)
printC(object$col$contrib,"Contributions (by-column total=100)", nword, nEig,1)
printC(object$col$cos2,"Square cosinus (by-row total=1)", nword, nEig,1)
printI(object$col$inertia,object$col$coord, "", nword)}

if(metaWords!=FALSE){
cat("\n\nWords whose contribution is over ", object$meta$lmw , " times the average word contribution\n")
for(i in 1:nEig) {
 cat("\nDimension ", i, "+\n")
 cat(object$meta$Word$Words[which(object$meta$Word$Dim==i & object$meta$Word$Coordinates>0)])
 cat("\nDimension ", i, "-\n")
 cat(object$meta$Word$Words[which(object$meta$Word$Dim==i & object$meta$Word$Coordinates<0)],"\n")
}}


nSeg<- ifelse(is.null(object$segment$coord),0 ,nrow(object$segment$coord))
if(nseg>0) {
cat("\n\nREPEATED SEGMENTS\n")
nS <- min(nSeg,nseg)
   printC(object$segment$coord,"Coordinates", nS , nEig, 1)
   printC(object$segment$cos2,"Square cosinus", nS , nEig, 1)
 }



if(nsup!=0) {

if(!is.null(object$row.sup)) {
  if(nsup=="ALL") nDsup<-nrow(object$row.sup$coord) 
     else nDsup <- min(nrow(object$row.sup$coord),nsup)} else {nDsup<-0}

nWsup<- ifelse(is.null(object$col.sup$coord),0 ,nrow(object$col.sup$coord)-nSeg)
if(nWsup>0) {if(nsup!="ALL") nWsup<- min(nWsup,nsup)} 


if(!is.null(object$quali.sup)) {
  nQsup<- nrow(object$quali.sup$coord)
  if(nsup!="ALL") nQsup<- min(nQsup,nsup)}else {nQsup<-0}
if(!is.null(object$quanti.sup)) nQQsup<- nrow(object$quanti.sup$coord) else nQQsup<-0
if(nDsup + nWsup+ nQsup + nQQsup>0) {


if(nDsup>0){
 cat("\n\nSUPPLEMENTARY DOCUMENTS\n")
 printC(object$row.sup$coord,"Coordinates", nDsup, nEig,1)
 printC(object$row.sup$cos2,"Square cosinus", nDsup, nEig,1)}

if(nWsup>0) {
 cat("\n\nSUPPLEMENTARY WORDS\n")
 printC(object$col.sup$coord,"Coordinates", nWsup,nEig,1)
 printC(object$col.sup$cos2,"Square cosinus", nWsup, nEig,1)}


if(nQsup>0) {
 cat("\n\nSUPPLEMENTARY CATEGORIES\n")
 printC(object$quali.sup$coord,"Coordinates",nQsup, nEig,1)
 printC(object$quali.sup$cos2,"Square cosinus", 
  nQsup, nEig,1)
 if(!is.null(object$quali.sup$v.test)) {
 printC(object$quali.sup$v.test,"v.test of the supplementary categories", 
   nQsup, nEig,1)}}


if(nQQsup>0) {
cat("\n\nSUPPLEMENTARY QUANTITATIVE VARIABLES\n")
   cat("\nCoordinates\n")
   cat("\nEqual to the correlations of the variable with the axes\n")
   print(data.frame(round(object$quanti.sup$coord[,1:nEig,drop=FALSE],nb.dec),
      row.names=rownames(object$quanti.sup$coord)))
   cat("\nSquare cosinus\n")
   print(data.frame(round(object$quanti.sup$cos2[,1:nEig,drop=FALSE],nb.dec),
      row.names=rownames(object$quanti.sup$cos2)))}
}
}


if(file!=""){
  sink()
cat("\nAll the results are in file: ",file,"\n")
}
}
