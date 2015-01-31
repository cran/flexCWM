\name{plot.cwm}
\alias{plot.cwm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot for CWMs
}
\description{
Plot method for cwm class objects.
}
\usage{
 \method{plot}{cwm}(x, regr = TRUE, ctype = c("Xnorm","Xbin","Xpois",
                 "Xmult"), which = NULL, criterion = "BIC", k = NULL,
                 modelXnorm = NULL, familyY = NULL, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{cwm}.
}
  \item{regr}{boolean, allows for bivariate regression plot.}
  \item{ctype}{a vector with concomitant variables types to plot.}
  \item{which}{a vector with columns number to plot, or "all" for all the columns}
    \item{criterion}{a string with the information criterion to consider; supported values are:\code{"AIC", "AICc", "AICu", "AIC3", "AWE", "BIC", "CAIC", "ICL".} Default value is \code{"BIC"}. 
}  
  \item{k}{an optional vector containing the numbers of mixture components to consider. If not specified, all the estimated models are considered.
  }
\item{modelXnorm}{an optional vector of character strings indicating the parsimonious models to consider for \code{Xnorm}. If not specified, all the estimated models are considered. 
  }
\item{familyY}{an optional vector of character strings indicating the conditional distribution of \eqn{Y} in each mixture component to consider. If not specified, all the estimated models are considered. 
  }
  \item{\dots}{further arguments for \code{\link[graphics:plot]{plot}}.
}


}

\examples{
data("students")
attach(students)
str(students)
fit2 <- cwm(WEIGHT ~ HEIGHT + HEIGHT.F , Xnorm = cbind(HEIGHT, HEIGHT.F), k = 2, 
  initialization = "kmeans", modelXnorm = "EEE")
summary(fit2, concomitant = TRUE)
plot(fit2)

}
