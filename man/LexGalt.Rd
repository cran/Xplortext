\name{LexGalt}
\encoding{UTF-8}
\alias{LexGalt}

\title{Correspondence Analysis on a Simple or Multiple Generalized Aggregate Lexical Table (LexGalt)}

\description{
Performs an extension of correspondence analysis on either a simple or a multiple generalized aggregated lexical table. In the case of a multiple table, a multiple factor analysis approach is used
}

\usage{
LexGalt(object, context="ALL", conf.ellip =FALSE, nb.ellip = 100, graph=TRUE, 
        axes = c(1, 2), label.group=NULL)
}

\arguments{
 \item{object}{object or list of objects (s) of TextData class (see details)}
 \item{context}{column index(es) or name(s) of the contextual variables (either qualitative or quantitative) used to build the generalized aggregated lexical table(s). These variables must have been previously selected in TextData function (by default "ALL")}
 \item{conf.ellip}{computing confidence ellipses (available only in the case of a simple table) (by default FALSE)}
 \item{nb.ellip}{number of samples drawn to evaluate the stability of the points (by default 100) only if conf.ellip= TRUE }
 \item{graph}{if TRUE, all several graphs are displayed; use \code{\link{plot.LexGalt}} to obtain detailed graphs (by default TRUE)}
 \item{axes}{length-2 vector indicating the axes to plot (by default axes=c(1,2))}
 \item{label.group}{In the case of analyzing a multiple generalized aggregated lexical table, vector containing the name of the groups (by default, NULL and the group are named GROUP.1, GROUP.2 and so on)}
}


\value{
Returns a list including an object named SQL if the simple qualitative analysis is performed, SQN for simple quantitative analysis, MQL for multiple qualitative analysis or
 MQN for multiple quantitative analysis (see details):
  \item{eig}{eigenvalues, percentages of inertia and cumulative percentages of inertia}
  \item{word}{the results for the words (coordinates, square cosine, contributions)}
  \item{quali.var}{results for the categorical variables (coordinates of each categories of each variables, square cosines)}
  \item{quanti.var}{results for the quantitative variables (coordinates, correlation between variables and axes, square cosines)}
  \item{ellip}{coordinates for confidence ellipses (words and categories) are drawn}
  \item{group}{in the case of multiple analysis, results for the groups (coordinates, contributions and square cosines) (MQL or MQN)} 
Returns the factor maps. The plots may be improved using the plot.LexGalt function. 
}

\details{


The default "context" argument is "ALL" and may contain qualitative and/or quantitative variables (names or indexes). If both types of variables are included, two independent LexGalt analyses are performed, 
saving the results for the qualitative analysis into an object named SQL (or MQL in the multiple case) and for the quantitative analysis into the SQN object (or MQN in the multiple case). 

In the multiple case, each TextData object must be created from as many executions of the function TextData as there are tables. They are joined in a list in the call to LexGalt function:

LexGalt(list(object1,object2,object3),...).

The variable names of each object in the list must be the same as the name of the variables selected in object1.

}

\references{

Bécue-Bertaut M. and Pagès J. (2015). Correspondence analysis of textual data involving contextual information: CA-GALT on principal components. Advances in Data Analysis and Classification, vol.(9) 2: 125-142. 
\doi{10.1007/s11634-014-0171-9}

Bécue-Bertaut M., Pagès J. and Kostov B. (2014). Untangling the influence of several contextual variables on the respondents' lexical choices. A statistical approach. SORT - Statistics and Operations Research Transactions, vol.(38) 2: 285-302.

Kostov B. A. (2015). A principal component method to analyse disconnected frequency tables by means of contextual information. (Doctoral dissertation). Retrieved from \url{http://upcommons.upc.edu/handle/2117/95759}.

Kostov, B., Bécue-Bertaut, M., & Husson, F. (2015). Correspondence Analysis on Generalised Aggregated Lexical Tables (CA-GALT) in the FactoMineR Package. 
The R Journal, Vol.7, Num.1, 109-117. \doi{10.32614/RJ-2015-010}
}

\author{Belchin Kostov, Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{plot.LexGalt}}}

\examples{

data(open.question)

res.TD<-TextData(open.question,var.text=c(9,10), Fmin=10, Dmin=10,
 context.quali=c("Gender", "Age_Group", "Education"),
 remov.number=TRUE, stop.word.tm=TRUE)

res.LexGalt <- LexGalt(res.TD, graph=FALSE, conf.ellip =FALSE)
plot(res.LexGalt, selQualiVar="ALL")

}

\keyword{multivariate}
