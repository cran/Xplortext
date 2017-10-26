\name{TextData}
\encoding{latin1}
\alias{TextData}

\title{Building textual and contextual tables (TextData)}

\description{
Creates a textual and contextual working-base (TextData format) from a source-base (data frame format).
}

\usage{
TextData(base, var.text=NULL, var.agg=NULL, context.quali=NULL, context.quanti= NULL,
 selDoc="ALL", lower=TRUE, remov.number=TRUE,lminword=1, Fmin=Dmin,Dmin=1, Fmax=Inf,
 stop.word.tm=FALSE, idiom="en", stop.word.user=NULL, segment=FALSE,
 sep.strong="\u005B()\u00BF?./:\u00A1!=+;{}-\u005D", seg.nfreq=10, seg.nfreq2=10,
 seg.nfreq3=10, graph=FALSE)
}
\arguments{
 \item{base}{source data frame with at least one textual column}
 \item{var.text}{vector with index(es) or name(s) of the selected textual column(s) (by default NULL)}
 \item{var.agg}{index or name of the aggregation categorical variable (by default NULL)}
 \item{context.quali}{vector with index(es) or name(s) of the selected categorical variable(s) (by default NULL)}
 \item{context.quanti}{vector with index(es) or name(s) of the selected quantitative variable(s) (by default NULL)}
 \item{selDoc}{vector with index(es) or name(s) of the selected source-documents (rows of the source-base) (by default "ALL")}
 \item{lower}{if TRUE, the corpus is converted into lowercase (by default TRUE)}
 \item{remov.number}{if TRUE, numbers are removed (by default TRUE)}  
 \item{lminword}{minimum length of a word to be selected (by default 1)}
 \item{Fmin}{minimum frequency of a word to be selected (by default Dmin)}
 \item{Dmin}{a word has to be used in at least Dmin source-documents to be selected (by default 1)}
 \item{Fmax}{maximum frequency of a word to be selected (by default Inf)}
 \item{stop.word.tm}{if TRUE, stoplist automatically provided in accordance with the idiom (by default FALSE)}
 \item{idiom}{declared idiom for the textual column(s) (by default English "en", see IETF language in package NLP)}
 \item{stop.word.user}{stoplist provided by the user}
 \item{segment}{if TRUE, the repeated segments are identified (by default FALSE)}
 \item{sep.strong}{string with the characters marking out the repeated segments (by default "[()¿?./:¡!=+;{}-]\")}
 \item{seg.nfreq}{minimum frequency of a more-than-three-words-long repeated segment (by default 10)}
 \item{seg.nfreq2}{minimum frequency of a two-words-long repeated segment (by default 10)}
 \item{seg.nfreq3}{minimum frequency of a three-words-long repeated segment (by default 10)}
 \item{graph}{if TRUE, documents, words and repeated segments barcharts are displayed; use plot.TextData to use more options (by default FALSE)}
}


\value{
A list including:
\item{summGen}{general summary}
\item{summDoc}{document summary}
\item{indexW}{index of words}
\item{DocTerm}{working lexical table (non-aggregate or aggregate table depending on var.agg value); 
working-documents by words table in slam package compressed format}
\item{context}{contextual variables if context.quali or context.quanti are non-NULL; the structure greatly differs 
in accordance with the nature of DocTerm table (non-aggregate/ aggregate), see details}
\item{info}{information about the selection of words}
\item{var.agg}{a one-column data frame with the values of the aggregation variable; NULL if non-aggregate analysis}
\item{SourceTerm}{in the case of DocTerm being an aggregate analysis, the source-documents by words table is kept in 
this data structure, in slam package compressed format}
\item{indexS}{working-documents by repeated-segments table, in slam package compressed format}
\item{remov.docs}{vector with the names of the removed empty source-documents}
}

\details{
Each row of the source-base is considered as a source-document.
TextData function builds the working-documents-by-words table, submitted to the analysis.

Information related to context.quanti and context.quali arguments:
  \enumerate{
    \item If numeric, contextual variables can be included in both vectors. The function TextData converts the numeric variable into factor to include it in context.quali vector. This possibility is interesting in some cases. For example, when treating open-ended questions, we can be interested in computing the correlation between the contextual variable "Age" and the axes and, at the same time, to draw the trajectory of the different values of "Age" (year by year) on the CA maps. 
    \item In the case of one or several columns with textual data not selected in vector var.text, if the argument context.quali is equal to "ALL", these columns will be considered as categorical variables. 
  }


Non-aggregate table versus aggregate table.

If var.agg=NULL: 
\enumerate{
 \item The work-documents are the non-empty-source-documents.
 \item DocTerm: non-aggregate lexical table with:
 \tabular{ll}{
 \tab as many rows as non-empty source-documents\cr
 \tab as many columns as words are selected.
  }
\item context$quali: data frame crossing the non-empty source-documents (rows) and the categorical contextual-variables (columns).
\item context$quanti: data frame crossing the non-empty source-documents (rows) and the quantitative contextual-variables (columns).
Both contextual tables can be juxtaposed row-wise to DocTerm table. 
}

If var.agg is NON-NULL:
\enumerate{
 \item The work-documents are aggregate-documents, issued from aggregating the source-documents depending on the categories of the aggregation variable; the aggregate-documents inherit the names of the corresponding categories.
 \item DocTerm is an aggregate table with:
 \tabular{ll}{
 \tab as many rows as as categories the aggregation variable has\cr
 \tab as many columns as words are selected.
  }
\item context$quali$qualitable: juxtaposes as many supplementary aggregate tables as categorical contextual variables. Each table has 
 \tabular{ll}{
 \tab as many rows as categories the contextual categorical variable has\cr
 \tab as many columns as selected words, i.e. as many columns as DocTerm has.
  }
\item context$quali$qualivar: names of categories of the supplementary categorical variables.
\item context$quanti: data frame crossing the working aggregate-documents (rows) and the quantitative contextual-variables (columns).
 The value for an active aggregate-document is the mean-value of the source-documents belonging to this aggregate-document.
}
}


\references{
Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.). <doi:10.1007/978-94-017-1525-6>.
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{print.TextData}}, \code{\link{summary.TextData}},  \code{\link{plot.TextData}}}

\examples{
# Non aggregate analysis
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), remov.number=TRUE, Fmin=10, Dmin=10,  
 stop.word.tm=TRUE, context.quali=c("Gender","Age_Group","Education"), context.quanti=c("Age"))

# Aggregate analysis and repeated segments
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Gen_Age", remov.number=TRUE, 
 Fmin=10, Dmin=10, stop.word.tm=TRUE, context.quali=c("Gender","Age_Group","Education"), 
 context.quanti=c("Age"), segment=TRUE)
}

\keyword{multivariate}

