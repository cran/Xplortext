\name{plot.TextData}
\encoding{UTF-8}
\alias{plot.TextData}

\title{Plot TextData objects}

\description{
Draws the barcharts of the longest documents, most frequent words and segments from a TextData object.
}

\usage{
\method{plot}{TextData}(x, ndoc=25, nword=25, nseg=25, sel=NULL, ordFreq=TRUE, 
 stop.word.tm=FALSE, idiom="en", stop.word.user=NULL, theme=theme_bw(), title=NULL,
 xtitle=NULL, col.fill="grey", col.lines="black", text.size=12, freq=NULL, vline=NULL, 
 interact=FALSE, round.dec = 4,...) 
}
\arguments{
  \item{x}{object of TextData class}
  \item{ndoc}{number of documents in the barchart (by default 25)}
  \item{nword}{number of words in the barchart  (by default 25)}
  \item{nseg}{number of segments in the barchart (by default 25)}
  \item{sel}{type of barchart (doc, word or seg for documents, words or repeated segments) (by default NULL and all he graphs are drawn)), see details}
   \item{ordFreq}{if ordFreq=TRUE, glossaries of words and repeated segments, are drawn in frequency order; if ordFreq=FALSE, glossaries are drown in alphabetic order (by default TRUE)}
  \item{stop.word.tm}{if TRUE, the tm stopwords (if the words are selected in TextData object) are not considered for the barchart (by default FALSE)}
  \item{idiom}{declared idiom for the textual column(s) (by default English "en", see IETF language in package NLP)}
  \item{stop.word.user}{the user's stopwords (if the words are selected in TextData object) are not considered for the barchart (by default NULL)}
  \item{theme}{theme settings (see ggplot2 package; by default theme_bw())}
  \item{title}{title of the graph (by default NULL and the title is automatically assigned)}
  \item{xtitle}{x title of the graph (by default NULL and the x title is automatically assigned)}
  \item{col.fill}{background color for the barChart bars (by default grey)}
  \item{col.lines}{lines color for the barChart bars (by default black)}
  \item{text.size}{text font size (by default 12)}
  \item{freq}{add frequencies to word and document barplots, see details (by default NULL)}
  \item{vline}{if "YES" or TRUE add vertical line to barplot, see details (by default NULL)}
  \item{interact}{if FALSE a ggplot graph, if TRUE an interactive plotly graph, see details (by default FALSE)}
  \item{round.dec}{number of decimals (by default 4)}
  \item{\dots}{further arguments passed to or from other methods...}
}



\details{
freq adds frequencies to barplot (by default NULL). If "YES" or TRUE displays the frequencies at the right of the bars at +5 position.
Numerical values display the frequencies at the right positions (positive values) or at the left (negative values).

vline adds two vertical line to word and document barplot (by default NULL). If TRUE a first vertical row line is added at mean level computed from the selected items from TextData,
and a second vertical blue line with the frequency mean of words/documents selected to plot in plot.TextData. If row and blue lines are the same, only blue line is shown.
If vline is a number, a line is show with this value.

Barchart selected in sel argument (doc, word and/or repeated segments) is in ggplot format. 
Barchart is used with geom_bar function of ggplot package. If it is only one element in sel argument the plot can
be saved in ggplot format:  newobject <- plot(TextDataObject,sel="word")

Selection of docs, words or segments can be done by numbers sel=list(type="doc", select=c(1,2:4,6)) or names sel= list(type="doc", select=c("M31_55", "M>55")).

If interact, rank for words/docs/segments from TextData selection are shown.

}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Mónica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

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
