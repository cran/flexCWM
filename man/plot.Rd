\name{plot.cwm}
\alias{plot.cwm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Uncertainty Plot for Generalized Linear Mixed CWMs
}
\description{
This plot shows the classification uncertainty in Generalized Linear Mixed CWMs
}
\usage{
 \method{plot}{cwm}(x, criteria="BIC", regr=TRUE, ctype =c("Xnorm", "Xbin", "Xpois", "Xmult"),
        which="all",...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{cwm}.
}
  \item{criteria}{an optional character string. Among the models tested, it selects the model that fits best according to an information criteria; supported values are: \code{"AIC", "AICc", "AICu", "AIC3", "AWE","BIC", "CAIC", "ICL"}}.
  \item{regr}{boolean, allows for bivariate regression plot.}
  \item{ctype}{a vector with concomitant variables types to plot.}
  \item{which}{a vector with columns number to plot, or "all" for all the columns}
  \item{\dots}{further arguments for \code{\link[graphics:plot]{plot}}.
}


}

\examples{
data(tourism)
Y <- tourism$overnights
X <- tourism$attendance
#res <- cwm(Y=Y,Xcont=X,k=1:4,seed=1)
#summary(res)
#plot(res)
}
