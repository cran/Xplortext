\name{LexGalt}
\encoding{UTF-8}
\alias{LexGalt}

\title{Generalised Aggregate Lexical Table (LexGalt)}

\description{
Performs Dual Correspondence Analysis in Generalised Aggregate Lexical Table
}

\usage{
LexGalt(object, context="ALL", ncp=NULL, level.ventil = 0, scale=TRUE,
 conf.ellip =FALSE, nb.ellip = 100, graph=TRUE, axes = c(1, 2), label.group=NULL)
}

\arguments{
 \item{object}{object(s) of TextData class (see details)}
 \item{context}{column index(es) or name(s) of the contextual qualitative or quantitative variables among those selected in TextData function (by default "ALL")}
 \item{ncp}{number of dimensions kept to compute (by default NULL to indicate all dimensions)}
 \item{level.ventil}{a proportion corresponding to the level under which the category is ventilated; by default, 0 and no ventilation is done}
 \item{scale}{variables are are scaled to unit variance (standardized) (by default TRUE)}
 \item{conf.ellip}{if condifence ellipses are computed (by default FALSE)}
 \item{nb.ellip}{number of samples drawn to evaluate the stability of the points only if conf.ellip=TRUE (by default 100)}
 \item{graph}{if TRUE, all several graphs are displayed; use \code{\link{plot.LexGalt}} to obtain detailed graphs (by default TRUE)}
 \item{axes}{length-2 vector indicating the axes to plot (by default axes=c(1,2))}
 \item{label.group}{a vector containing the name of the groups (by default, NULL and the group are named GROUP.1, GROUP.2 and so on)}
}


\value{
Returns a list including SQL object for simple qualitative analysis, SQN object for simple quantitative analysis, MQL object for multiple qualitative analysis,
 MQN object for multiple quantitative analysis:
  \item{eig}{matrix with the eigenvalues, the percentages of inertia and the cumulative percentages of inertia}
  \item{doc}{a list of matrices containing all the results for the documents (coordinates, square cosine)}
  \item{word}{a list of matrices containing all the results for the words (coordinates, square cosine, contributions)}
  \item{quali.var}{a list of matrices containing all the results for the categorical variables (coordinates of each categories of each variables, square cosine)}
  \item{quanti.var}{a list of matrices containing all the results for the quantitative variables (coordinates, correlation between variables and axes, square cosine)}
  \item{ellip}{coordinates for confidence ellipses are drawn. In simple qualitative analysis (SQL) they are words and categories}
  \item{group}{coordinates, contributions and square cosine for each group in multiple analysis (MQL or MQN)} 
Returns the factor maps. The plots may be improved using the plot.LexGalt function. 
}

\details{
If multiple analysis, each TextData object must be created in advance (object1, object2, object3) and joined in a list:
LexGalt(list=c(object1,object2,object3),...)
Ellipses only are computed for simple analysis.
}

\references{

Bécue-Bertaut M. and Pagès J. (2015). Correspondence analysis of textual data involving contextual information: CA-GALT on principal components. Advances in Data Analysis and Classification, vol.(9) 2: 125-142.

Bécue-Bertaut M., Pagès J. and Kostov B. (2014). Untangling the influence of several contextual variables on the respondents' lexical choices. A statistical approach. SORT - Statistics and Operations Research Transactions, vol.(38) 2: 285-302.

Kostov B. A. (2015). A principal component method to analyse disconnected frequency tables by means of contextual information. (Doctoral dissertation). Retrieved from \url{http://upcommons.upc.edu/handle/2117/95759}.

}

\author{Belchin Kostov, Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{plot.LexGalt}}}

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
\keyword{multivariate}
