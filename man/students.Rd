\name{students}
\docType{data}
\alias{students}
\title{dataset students}
\description{
A dataframe with data from a survey of 270 students attending a statistics course at the
Department of Economics and Business of the University of Catania in the academic year 2011/2012. 
It contains the following variables:
\itemize{
    \item \code{GENDER} gender of the respondent;
    \item \code{HEIGHT} height of the respondent, measured in centimeters;
    \item \code{WEIGHT} weight of the respondent, measured in kilograms;
    \item \code{HEIGHT.F} height of respondent's father, measured in centimeters.
  }
}
\usage{data(students)}

\format{A dataset}
\source{\url{ http://www.economia.unict.it/punzo/}}
\references{
Ingrassia, S., Minotti, S. C., and Punzo, A. (2014). Model-based clustering via linear cluster-weighted models. \emph{Computational Statistics and Data Analysis}, \bold{71}, 159-182.

}
\seealso{
\code{\link{flexCWM-package}}, \code{\link{cwm}}
}
\keyword{datasets}
\examples{
data("students")
attach(students)
str(students)
fit2 <- cwm(WEIGHT ~ HEIGHT + HEIGHT.F , Xnorm = cbind(HEIGHT, HEIGHT.F), k = 2, 
  initialization = "kmeans", modelXnorm = "EEE")
summary(fit2, concomitant = TRUE)
plot(fit2)

}