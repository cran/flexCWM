\name{gcwm}
\alias{gcwm}
\title{Flexible Cluster-Weighted Modeling}
\description{Run the Generalized Linear Gaussian  Cluster-Weighted Model }
\usage{
gcwm (Y, X, familyY="Gaussian", k=2, ic=c("BIC", "AIC", "ICL"), mY=1, 
  method="Nelder-Mead", initialization="random.soft", start.z=NULL, 
  iter.max=1000, threshold=1.0e-04, loglikplot=FALSE, seed=NULL)              
}
                 
%- maybe also 'usage' for other objects documented here.
\arguments{                      
   \item{Y}{
  numerical vector for the response variable.
  }
  \item{X}{
  matrix for the covariates. 
  }
  \item{ familyY}{
  the exponential distribution used for Y|x in each cluster; it can be: 
    \itemize{
      \item \code{"Gaussian"}
      \item \code{"Poisson"}
      \item \code{"Binomial"}
      \item \code{"Gamma"}
    }
  Default value is "Gaussian".
  }

  \item{k}{
  a vector containing the numbers of clusters to be tried. The one with the lowest information criterion is selected. Default value is 2.
  }
  \item{ic}{the information criteria by which the best model is selected when \code{lencgth(k)>1}. Possible values are
    \itemize{
      \item \code{"BIC"} 
      \item \code{"AIC"} 
      \item \code{"ICL"}
    }
  }
  \item{mY}{
  When \code{familyY="Binomial"}, it sets the sample size. Default value is 1 (Bernoulli distribution).
  } 

  \item{method}{
  optimization method used in the M-step of the EM algorithm (see \code{\link[stats:optim]{optim}}). Default value is \code{"Nelder-Mead"}.
  }
  \item{initialization}{
  initialization strategy for the EM-algorithm. It can be:
  \itemize{
    \item \code{"random.soft"}
    \item \code{"random.hard"}
    \item \code{"manual"}
  }
  Default value is \code{"random.soft"}.
}
  \item{start.z}{
  matrix of soft or hard classification: it is used only if \code{initialization="manual"}.
}
  \item{iter.max}{
   maximum number of iterations in the EM-algorithm. Default value is 200.
}
  \item{threshold}{
  threshold for Aitken acceleration procedure. Default value is 1.0e-04.
}
  \item{loglikplot}{
  if TRUE, the log-likelihood values against the iterations are plotted. Default value FALSE.
}
  \item{seed}{
 the seed for the random number generator, when random initializations are used; if \code{NULL}, current seed is not changed. Default value is \code{NULL}.
}

}
\value{
This function returns a list of values related to the model selected. It contains:

\item{Y}{response variable}
\item{X}{covariates}
\item{familyY}{exponential distribution used for Y|x in each cluster}
\item{p}{number of covariates}
\item{k}{number of groups}
\item{n}{sample size}
\item{npar}{number of parameters}
\item{mY}{sample size, used when \code{familyY="Binomial"}}
\item{prior}{weights for the mixture components}
\item{muX}{covariates means}
\item{VarX}{covariates variances}
\item{PX}{marginal distribution of X for each cluster}
\item{beta}{regression coefficients}
\item{muY}{mean of Y}
\item{dispY}{dispersion parameter of Y}
\item{VarFunY}{variance function of Y}
\item{VarY}{variance of Y}
\item{nuY}{when familyY="Gamma", the gamma distribution is parameterized according to muY and nuY (see McCullagh, P. and Nelder, J. 1989)}
\item{PY}{conditional distribution of Y|x for each cluster}
\item{iter.stop}{number of iterations performed in EM algorithm}
\item{z}{matrix of posterior probabilities}
\item{group}{classification vector}
\item{loglik}{final log-likelihood value}
\item{AIC}{}
\item{BIC}{}
\item{ICL}{}
\item{call}{an object of class \code{call}}
}

\references{
Ingrassia, S., Minotti, S. C., and Vittadini, G. (2012). Local statistical modeling via the cluster-weighted approach with elliptical distributions. \emph{Journal of Classification}, \bold{29}(3), 363-401.

Ingrassia, S., Minotti, S. C., Punzo, A., and Vittadini, G. (2012). Generalized linear Gaussian cluster-weighted
modeling. \emph{arXiv.org} e-print \bold{1211.1171}, available at: http://arxiv.org/abs/1211.1171.

McCullagh, P. and Nelder, J. (1989). Generalized Linear Models. Chapman & Hall, Boca Raton, 2nd edition
}
\author{
Incarbone G., Mazza A., Punzo A., Ingrassia S.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{                  
\code{\link{flexCWM-package}}, \code{\link{tourism}}
}
\examples{
data(tourism)
Y <- tourism$overnights
X <- tourism$attendance
res <- gcwm(Y=Y,X=X,k=1:4,seed=1)
plot(cbind(Y,X),col=res$best$group)


}

