\name{plot.LexCA}
\encoding{latin1}
\alias{plot.LexCA}

\title{Plot of LexCA objects}

\description{
Plots textual correspondence analysis (CA) graphs from a LexCA object. 
}

\usage{
\method{plot}{LexCA}(x, selDoc="ALL", selWord="ALL", selSeg=NULL, selDocSup=NULL,
  selWordSup=NULL, quanti.sup=NULL, quali.sup=NULL, maxDocs=20, eigen=FALSE, 
  title=NULL, axes=c(1,2), col.doc="blue", col.word="red", col.doc.sup="darkblue", 
  col.word.sup="darkred", col.quanti.sup = "blue", col.quali.sup="darkgreen", 
  col.seg="cyan4", col="grey", cex=1, xlim=NULL, ylim=NULL, shadowtext=FALSE,
  habillage="none", unselect=1, label="all", autoLab=c("auto", "yes", "no"), 
  plot.new=TRUE, ...)
}
\arguments{
  \item{x}{object of LexCA class}
  \item{selDoc}{vector with the active documents to plot (indexes, names or rules; see details; by default "ALL")}
  \item{selWord}{vector with the active words to plot (indexes, names or rules; see details; by default "ALL")}
  \item{selSeg}{vector with the supplementary repeated segments to plot (indexes, names or rules; see details; by default NULL)}
  \item{selDocSup}{vector with the supplementary documents to plot (indexes, names or rules; see details; by default NULL)}
  \item{selWordSup}{vector of the supplementary words to plot (indexes, names or rules; see details; by default NULL)} 
  \item{quanti.sup}{vector of the supplementary quantitative variables to plot (indexes, names or rules; see details; by default NULL)}
  \item{quali.sup}{vector with the supplementary categorical variables/categories to plot (indexes, names or rules; see details; by default NULL). The selected categories (through the variables or directly) are plotted}
  \item{maxDocs}{limit to the number of active documents in the lexical table when selecting the words to be plotted for being characteristic of the selected documents (by default 20)}
  \item{eigen}{if TRUE, the eigenvalues barplot is drawn (by default FALSE); no other elements can be simultaneously selected}
  \item{title}{title of the graph (by default NULL and the title is automatically assigned)}
  \item{axes}{length-2 vector indicating the axes considered in the graph (by default c(1,2))}
  \item{col.doc}{color for the point-documents (by default "blue")}
  \item{col.word}{color for the point-words (by default "red")}
  \item{col.doc.sup}{color for the supplementary point-documents (by default "darkblue")}
  \item{col.word.sup}{color for the supplementary point-words (by default "darkred")}
  \item{col.quanti.sup}{color for the quanti.sup variables (by default "blue")}
  \item{col.quali.sup}{color for the categorical supplementary point-categories, (by default "darkgreen")}
  \item{col.seg}{color for the supplementary point-repeated segments, (by default "cyan4")}
  \item{col}{color for the bars in the eigenvalues barplot (by default "grey")}
  \item{cex}{text and symbol size is scaled by cex, in relation to size 1 (by default 1)}
  \item{xlim}{range for 'x' values on the graph, defaulting to the finite values of 'x' range (by default NULL)}
  \item{ylim}{range for the 'y' values on the graph, defaulting to the the finite values of 'y' range (by default NULL)}
  \item{shadowtext}{if TRUE, shadow on the labels (rectangles are written under the labels which may lead to difficulties to modify the graph with another program) (by default FALSE)}
  \item{habillage}{index or name of the categorical variable used to differentiate the documents by colors given according to the category; by default "none")}
  \item{unselect}{either a value between 0 and 1 or a color. In the first case, transparency level of the unselected objects (if unselect=1 the transparency is total and the elements are not represented; if unselect=0 the elements are represented as usual but without any label); in the case of a color (e.g. unselect="grey60"), the non-selected points are given this color (by default 1)}
  \item{label}{a list of character for the variables which are labelled (by default NULL and all the drawn variables are labelled). You can label all the active variables by putting "var" and/or all the supplementary variables by putting "quanti.sup" and/or a list with the names of the variables which should be labelled. Value should be one of "all", "none", "row", "row.sup", "col", "col.sup", "quali.sup" or NULL.}
  \item{autoLab}{if autoLab="auto", autoLab turns to be equal to "yes" if there are less than 50 elements and equal to "no" otherwise; if "yes", the labels are moved, as little as possible, to avoid overlapping (time-consuming if many elements); if "no" the labels are placed quickly but may overlap}
  \item{plot.new}{if TRUE, a new graphical device is created (by default TRUE)}
  \item{\dots}{further arguments passed from other methods...}
}



\details{
The argument autoLab = "yes" is time-consuming if many overlapping labels. Furthermore, the visualization of the words cloud can result distorted because of the apparent greater dispersion of the words labels. An alternative would be reducing the character size of the words labels to reduce overlapping (e.g. cex=0.7).

selDoc, selWord, selSeg, selDocSup, selWordSup, quanti.sup and quali.sup allow for selecting all or part of the elements of the corresponding type, using either labels, indexes or rules.\cr

The syntax is the same for all types.

1. Using labels:
\preformatted{
selDoc = c("doc1","doc5"): only the documents with labels doc1 and doc5 are plotted.
quali.sup=c("varcateg1","category12"): only the categories (all of them) of 
   categorical variable labeled "varcateg1" and the category labeled "category12"
   are plotted.}

2.- Using indexes:
\preformatted{
selDoc = c(1:5): documents 1 to 5 are plotted.
quali.sup=c(1:5,7): categories 1 to 5 and 7 are plotted. The numbering of the
   categories have to be consulted in the LexCA numerical results.}

3.- Using rules:    
 Rules are based on the coordinates (coord), the contribution (contrib or meta; concerning only
 active elements) or the square cosine (cos2).\cr
 Somes examples are given hereafter:

\preformatted{
selDoc="coord 10": only the 10 documents with the highest coordinates, as globally
   computed on the 2 axes, are plotted.
selWord="contrib 10": the words with a contribution to the inertia, of any of 
   the 2 axes, over 10\% of the axis inertia are plotted.
selWord="meta 3": the words with a contribution over 3 times the average word 
   contribution on any of the two axes are plotted. Only active words or documents 
   can be selected.
selDocSup="cos2 .85": the supplementary documents with a cos2 over 0.85, as summed
   on the 2 axes, are plotted.
selWord="char 0.05": only the characteristic words of the documents selected in 
   SelDoc are plotted. The selection of the words follow the rationale used in 
   function LexChar using as limit for the p-value the value given, here.0.05.}

}



\references{
Husson F., Le S., Pages J. (2011). Exploratory Multivariate Analysis by Example Using R. Chapman & Hall/CRC.


}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

\seealso{\code{\link{LexCA}}, \code{\link{print.LexCA}}, \code{\link{summary.LexCA}}}

\examples{
data(open.question)
res.TD<-TextData(open.question,var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10,
        remov.number=TRUE, stop.word.tm=TRUE)
resCA <- LexCA(res.TD, graph=FALSE)
plot(resCA, selDoc="contrib 30", selWord="coord 20")
}


\keyword{plot}

