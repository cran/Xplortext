\name{plot.TextData}
\encoding{UTF-8}
\alias{plot.TextData}

\title{Plot TextData objects}

\description{
Draws the barcharts of the longest documents, most frequent words and segments from a TextData object.
}

\usage{
\method{plot}{TextData}(x, ndoc=25, nword=25, nseg=25, sel=NULL, stop.word.tm=FALSE, 
  stop.word.user=NULL, theme=theme_bw(), title=NULL, xtitle=NULL, col.fill="grey",
  col.lines="black", text.size=12, freq=NULL, vline=NULL,...) 
}
\arguments{
  \item{x}{object of TextData class}
  \item{ndoc}{number of documents in the barchart (by default 25)}
  \item{nword}{number of words in the barchart  (by default 25)}
  \item{nseg}{number of segments in the barchart (by default 25)}
  \item{sel}{type of barchart (doc, word or seg for documents, words or repeated segments) (by default NULL and all the barchart are draw)}
  \item{stop.word.tm}{the tm stopwords are not considered for the barchart (by default FALSE)}
  \item{stop.word.user}{the user's stopwords are not considered for the barchart (by default NULL)}
  \item{theme}{theme settings (see ggplot2 package; by default theme_bw())}
  \item{title}{title of the graph (by default NULL and the title is automatically assigned)}
  \item{xtitle}{x title of the graph (by default NULL and the x title is automatically assigned)}
  \item{col.fill}{background color for the barChart bars (by default grey)}
  \item{col.lines}{lines color for the barChart bars (by default black)}
  \item{text.size}{text font size (by default 12)}
  \item{freq}{add frequencies to word and document barplots, see details (by default NULL)}
  \item{vline}{if "YES" or TRUE add vertical line to barplot, see details (by default NULL)}
  \item{\dots}{further arguments passed to or from other methods...}
}

\details{
freq adds frequencies to barplot (by default NULL). If "YES" or TRUE displays the frequencies at the right of the bars at +5 position.
Numerical values display the frequencies at the right positions (positive values) or at the left (negative values).

vline adds a vertical line to barplot (by default NULL). If TRUE a vertical line is added at mean level. If "median" a vertical line is added at median level.
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{TextData}}, \code{\link{print.TextData}},  \code{\link{summary.TextData}}}

\examples{
# Non aggregate analysis

 data(open.question)
 res.TD<-TextData(open.question, var.text=c(9,10), remov.number=TRUE, Fmin=10, Dmin=10,  
 stop.word.tm=TRUE, context.quali=c("Gender","Age_Group","Education"), context.quanti=c("Age"))
 plot(res.TD)

# Aggregate analysis
 data(open.question)
 res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Gen_Age", remov.number=TRUE, 
 Fmin=10, Dmin=10, stop.word.tm=TRUE, context.quali=c("Gender","Age_Group","Education"), 
 context.quanti=c("Age"), segment=TRUE)
 plot(res.TD)

}

\keyword{plot}
