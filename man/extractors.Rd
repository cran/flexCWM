\name{Extractor functions}
\alias{getBestModel}
\alias{getPosterior}
\alias{getSize}
\alias{getCluster}
\alias{getPar}
\alias{getParPrior}
\alias{getParGLM}
\alias{getParConcomitant}
\alias{getParXnorm}
\alias{getParXbin}
\alias{getParXpois}
\alias{getParXmult}
\alias{getIC}
\alias{whichBest}
\alias{summary.cwm}
\alias{print.cwm}
\alias{print.cwm.IC}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extractors for \code{cwm} class objects.
}
\description{
These functions extract values from \code{cwm} class objects. 
}
\usage{

getBestModel(object, criterion = "BIC", k = NULL, modelXnorm = NULL, familyY = NULL)
getPosterior(object, ...)
getSize(object, ...)
getCluster(object, ...)
getParGLM(object, ...)
getParConcomitant(object, name = NULL, ...)
getPar(object, ...)
getParPrior(object, ...)
getParXnorm(object, ...)
getParXbin(object, ...)
getParXpois(object, ...)
getParXmult(object, ...)
getIC(object,criteria)
whichBest(object, criteria = NULL, k = NULL, modelXnorm = NULL, familyY = NULL)

\method{summary}{cwm}(object, criterion = "BIC",  concomitant = FALSE, 
  digits = getOption("digits")-2, \dots)
\method{print}{cwm}(x, \dots)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{object, x}{a class \code{cwm} object.
  }
    \item{criterion}{a string with the information criterion to consider; supported values are:\code{"AIC", "AICc", "AICu", "AIC3", "AWE", "BIC", "CAIC", "ICL".} Default value is \code{"BIC"}. 
}  
    \item{criteria}{a vector of strings  with the names of information criteria to consider. If \code{NULL} all the supported infromation criteria are considered.
}
  \item{k}{an optional vector containing the numbers of mixture components to consider. If not specified, all the estimated models are considered.
  }
\item{modelXnorm}{an optional vector of character strings indicating the parsimonious models to consider for \code{Xnorm}. If not specified, all the estimated models are considered. 
  }
\item{familyY}{an optional vector of character strings indicating the conditional distribution of \eqn{Y} in each mixture component to consider. If not specified, all the estimated models are considered. 
  }
  \item{name}{an optional vector of strings specifing the names of distribution families of concomitant variables; if \code{NULL}, parameters estimated for all concomitant variables are returned.}
    \item{concomitant}{When \code{TRUE}, concomitant variables parameters are displayed. Default is \code{FALSE}.
}
  \item{digits}{ 
        integer used for number formatting.
}
  \item{\dots}{   
  additional arguments to be passed to \code{getBestModel} (or to \code{whichBest} for the \code{print} method). 
}
  
}
\details{
When several models have been estimated, these functions consider the best model according to the information criterion in \code{criterion}, among the estimated models having a number of components among those in \code{k} an error distribution among those in \code{familyY} and a parsimonious model among those in \code{modelXnorm}.
\code{getIC} provides values for the information criteria in \code{criteria}.

The \code{getBestModel} method returns a \code{cwm} object containing the best model only, selected as described above.
}

\examples{
#res <- cwm(Y=Y,Xcont=X,k=1:4,seed=1)
#summary(res)
#plot(res)
}
