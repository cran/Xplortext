\name{plot.LexGalt}
\encoding{UTF-8}
\alias{plot.LexGalt}

\title{Plot LexGalt objects}

\description{
Plots Generalised Aggregate Lexical Tables (LexGalt) graphs from a LexGalt object
}

\usage{
\method{plot}{LexGalt}(x,type="QL", selDoc=NULL, selWord=NULL, selQualiVar=NULL,
  selQuantiVar=NULL, conf.ellip=FALSE, selWordEllip=NULL, selQualiVarEllip=NULL,
  selQuantiVarEllip=NULL, level.conf=0.95, eigen=FALSE, title = NULL, axes = c(1, 2),
  xlim = NULL, ylim = NULL, col.eig="grey", col.doc = "black", col.word = NULL,
  col.quali = "blue", col.quanti = "blue", col="grey", pch = 20, label = TRUE, 
  autoLab = c("auto", "yes", "no"), palette = NULL, unselect = 1, 
  selCov=FALSE, selGroup="ALL", partial=FALSE, plot.group=FALSE, 
  col.group=NULL, label.group=NULL, legend=TRUE, pos.legend="topleft", 
  new.plot = TRUE, cex=1,...)
}
\arguments{
  \item{x}{object of LexGalt class}
  \item{type}{results from a qualitative analysis (type="QL") or quantitative analysis (type="QN"); see details; by default Q)}
  \item{selDoc}{vector with the documents to plot (indexes, names or rules; see details; by default NULL)}
  \item{selWord}{vector with the words to plot (indexes, names or rules (indexes, names or rules; see details; by default NULL)}
  \item{selQualiVar}{vector with the categories of categorical variables to plot (indexes, names or rules; see details; by default NULL)}
  \item{selQuantiVar}{vector with the numerical variables to plot (indexes, names or rules; see details; by default NULL)}
  \item{conf.ellip}{to drawn confidence ellipses, by default FALSE}
  \item{selWordEllip}{vector with the words that defines which ellipses are drawn (indexes, names or rules; see details; by default NULL)}
  \item{selQualiVarEllip}{vector with the categories of categorical variables which ellipses are drawn (indexes, names or rules; see details; by default NULL)}
  \item{selQuantiVarEllip}{vector with the numerical variables which ellipses are drawn(indexes, names or rules; see details; by default NULL)}
  \item{level.conf}{level of confidence used to construct the ellipses; by default 0.95}
  \item{eigen}{if TRUE, the eigenvalues barplot is drawn (by default FALSE); other elements can be simultaneously selected}
  \item{title}{title of the graph (by default NULL and the title is automatically assigned)}
  \item{axes}{length-2 vector indicating the axes considered in the graph; by default c(1,2)}
  \item{xlim}{range for 'x' values on the graph, defaulting to the finite values of 'x' range (by default NULL)}
  \item{ylim}{range for the 'y' values on the graph, defaulting to the the finite values of 'y' range (by default NULL)}
  \item{col.eig}{value or vector with colors for the bars of eigenvalues (by default "grey")}
  \item{col.doc}{color for the point-documents(by default "black")}
  \item{col.word}{color for the point-words (by default NULL is darkred in simple analysis; see details)}
  \item{col.quali}{color for the categories of categorical variables (by default "blue")}
  \item{col.quanti}{color for the numerical variables (by default "blue")}
  \item{col}{color for the bars in the eigenvalues barplot (by default "grey")}
  \item{pch}{plotting character for coordinates, cf. \code{\link[graphics]{points}} function in the graphics package}
  \item{label}{a list of character for the elements which are labelled (by default TRUE and all the drawn elements are labelled).}
  \item{autoLab}{if autoLab="auto", autoLab turns to be equal to "yes" if there are less than 50 elements and equal to "no" otherwise; if "yes", the labels are moved, as little as possible, to avoid overlapping (time-consuming if many elements); if "no" the labels are placed quickly but may overlap}
  \item{palette}{the color palette used to draw the points. By default colors are chosen. If you want to define the colors : palette=c("black", "red", "blue"); or you can use: palette=rainbow(10), or in black and white for example: palette=gray(seq(0,.9,len=3))}
  \item{unselect}{may be either a value between 0 and 1 that gives the transparency of the unselected objects (if unselect=1 the transparceny is total and the elements are not drawn, if unselect=0 the elements are drawn as usual but without any label) or may be a color (for example unselect="grey60")}
  \item{selCov}{a boolean, if TRUE then data are scaled to unit variance (by default TRUE)}
  \item{selGroup}{vector with the groups to plot if multiple analysis was performed (indexes, names or rules; see details; by default NULL)}
  \item{partial}{if TRUE partial elements (results for the groups) are shown, if ALL results for the conjoint analysis are superimposed; by default FALSE}
  \item{plot.group}{draw a plot comparing the groups in multiple case (by default TRUE)} 
  \item{col.group}{color for the groups if multiple analysis was performed (by default NULL and they are selected from palette)}
  \item{label.group}{a vector containing the new name of the groups. If "BLANK" no labels with the group are added at the end of the drawn elements (by default, NULL and the name of each group is added)}
  \item{legend}{show the legend of labels of groups. See \code{\link[graphics]{legend}} from graphics package (by default TRUE}
  \item{pos.legend}{position of the legend of labels of groups. See \code{\link[graphics]{legend}} from graphics package (by default "topleft")}
  \item{new.plot}{if TRUE, a new graphical device is created (by default TRUE)}
  \item{cex}{text and symbol size is scaled by cex, in relation to size 1 (by default 1)}
  \item{\dots}{further arguments passed from other methods...}
}

\details{
The argument autoLab = "yes" is time-consuming if many overlapping labels. Furthermore, the visualization of the words cloud can result distorted because of the apparent greater dispersion of the words labels. An alternative would be reducing the character size of the words labels to reduce overlapping (e.g. cex=0.7).

selDoc, selWord, selQualiVar, selQuantiVar, selWordEllip, selQualiVarEllip,
    selQuantiVarEllip allow for selecting all or part of the elements of the corresponding type, using either labels, indexes or rules.\cr

The syntax is the same for all types.

1. Using labels:
\preformatted{
selDoc = c("doc1","doc5"): only the documents with labels doc1 and doc5 are plotted.
selQualiVar=c("category1","category2"): only the categories labeled category1 and
 category2 are plotted.}

2.- Using indexes:
\preformatted{
selDoc = c(1:5): documents 1 to 5 are plotted.
quali.sup=c(1:5,7): categories 1 to 5 and 7 are plotted. The numbering of the
   categories have to be consulted in the LexGalt numerical results.}

3.- Using rules:    
 Rules are based on the coordinates (coord), the contribution (contrib or meta) or the square cosine (cos2).\cr
 Somes examples are given hereafter:

\preformatted{
selDoc="coord 10": only the 10 documents with the highest coordinates, as globally
   computed on the 2 axes, are plotted.
selWord="contrib 10": the words with a contribution to the inertia, of any of 
   the 2 axes.
selWord="meta 3": the words with a contribution over 3 times the average word 
   contribution on any of the two axes are plotted.
selWord="cos2 .85": the words with a cos2 over 0.85, as summed
   on the 2 axes, are plotted.
 }
\preformatted{
col.word by default NULL is "darkred" for simple analysis, if it is null takes
the colors from col.group 
i.e. col.group=c("red","blue"). To select the colors for some words in object res, 
we can use:
str.col.words <- rep("darkred",nrow(res$MQL$word$coord))
str.col.words[which(rownames(res$MQL$word$coord) == "kids")] <- "red"
str.col.words[which(rownames(res$MQL$word$coord) == "friends")] <- "green"
str.col.words[which(rownames(res$MQL$word$coord) == "job")] <- "pink"
plot(res, selGroup=1, selWord=c("friends", "job", "kids", "at"),new.plot=FALSE, 
col.group=c("darkred","blue"), autoLab = "yes", col.word=str.col.words)
} 
 
}


\references{
 
Bécue-Bertaut M. and Pagès J. (2015). Correspondence analysis of textual data involving contextual information: CA-GALT on principal components. Advances in Data Analysis and Classification, vol.(9) 2: 125-142.

Bécue-Bertaut M., Pagès J. and Kostov B. (2014). Untangling the influence of several contextual variables on the respondents' lexical choices. A statistical approach. SORT - Statistics and Operations Research Transactions, vol.(38) 2: 285-302.

Kostov B. A. (2015). A principal component method to analyse disconnected frequency tables by means of contextual information. (Doctoral dissertation). Retrieved from \url{http://upcommons.upc.edu/handle/2117/95759}.

}

\author{
Belchin Kostov, Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, 
Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{LexGalt}}}

\examples{
\dontrun{
data(open.question)

res.TD<-TextData(open.question,var.text=c(9,10),  Fmin=10, Dmin=10,
 context.quali=c("Gender", "Age_Group", "Education"),
 remov.number=TRUE, stop.word.tm=TRUE)

res.LexGalt <- LexGalt(res.TD, graph=FALSE, nb.ellip =0)
plot(res.LexGalt, selQualiVar="ALL")
}
}

\keyword{plot}

