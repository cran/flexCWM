\name{cwm}
\alias{cwm}
\alias{student.t}
\title{Fit for the CWM}
\description{Maximum likelihood fitting of the cluster-weighted model by the EM algorithm.}
\usage{
cwm(formulaY = NULL, familyY = gaussian, data, Xnorm = NULL, Xbin = NULL,
  Xpois = NULL, Xmult = NULL, modelXnorm = NULL, Xbtrials = NULL, k = 1:3, 
  initialization = c("random.soft", "random.hard", "kmeans", "mclust", "manual"), 
  start.z = NULL, seed = NULL, maxR = 1, iter.max = 1000, threshold = 1.0e-04, 
  eps = 1e-100, parallel = FALSE,   pwarning = FALSE)
}
                            
%- maybe also 'usage' for other objects documented here.
\arguments{                      
  \item{formulaY}{an optional object of class "\code{\link[stats:formula]{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.    
  }
    \item{familyY}{a description of the error distribution and link function to be used for the conditional distribution of \eqn{Y} in each mixture component. This can be a character string naming a \code{\link[stats:family]{family function}}, a family function or the result of a call to a family function. 
 The following family functions are supported:  
    \itemize{
      \item \code{binomial(link = "logit")} 
      \item \code{gaussian(link = "identity")} 
      \item \code{Gamma(link = "log")}
      \item \code{inverse.gaussian(link = "1/mu^2")} 
      \item \code{poisson(link = "log")}
      \item \code{student.t(link = "identity")} 
    }
  Default value is \code{gaussian(link = "identity")}.
  }
  \item{data}{
  an optional \code{\link[base:data.frame]{data.frame}}, \code{\link[base:list]{list}}, or \code{\link[base:environment]{environment}} with the variables needed to use \code{formulaY}.
  }
  \item{Xnorm, Xbin, Xpois, Xmult}{
  an optional matrix containing variables to be used for marginalization having normal, binomial, Poisson and multinomial distributions. 
  }
  \item{modelXnorm}{an optional vector of character strings indicating the parsimonious models to be fitted for variables in \code{Xnorm}. The default is  \code{c("E", "V")} for a single continuous covariate, and \code{c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")} for multivariate continuous covariates (see \code{\link[mixture:gpcm]{mixture:gpcm} for details}). 
  }
  \item{Xbtrials}{
   an optional vector containing the number of trials for each column in \code{Xbin}. If omitted, the maximum of each column in \code{Xbin} is used.
   }
  \item{k}{
  an optional vector containing the numbers of mixture components to be tried. Default value is \code{1:3}.
  }
  \item{initialization}{
  an optional character string. It sets the initialization strategy for the EM-algorithm. It can be:
  \itemize{
    \item \code{"random.soft"}
    \item \code{"random.hard"}
    \item \code{"kmeans"}
    \item \code{"mclust"} 
    \item \code{"manual"}
  }
   Default value is \code{"random.soft"}.
}
  \item{start.z}{
  matrix of soft or hard classification: it is used only if \code{initialization = "manual"}.Only models with the same number of mixture components as the columns of this matrix will be fit.
}
  \item{seed}{
 an optional scalar. It sets the seed for the random number generator, when random initializations are used; if \code{NULL}, current seed is not changed. Default value is \code{NULL}.
}

 \item{maxR}{
 number of initializations to be tried. Default value is 1.
  } 
  \item{iter.max}{
   an optional scalar. It sets the maximum number of iterations in the EM-algorithm. Default value is 200.
}
  \item{threshold}{
  an optional scalar. It sets the threshold for the Aitken acceleration procedure. Default value is 1.0e-04.
}  
  \item{eps}{
  an optional scalar. It sets the smallest value for eigenvalues of covariance matrices for \code{Xnorm}. Default value is 1e-100.
}
  \item{parallel}{
  When \code{TRUE}, the package \code{\link[parallel:parallel-package]{parallel}} is used for parallel computation. When several models are estimated, computational time is reduced. The number of cores to use may be set with the global option \code{cl.cores}; default value is detected using \code{\link[parallel:detectCores]{detectCores()}}.
} 
  \item{pwarning}{When \code{TRUE}, warnings are printed.}
}
\details{When \code{familyY = binomial}, the response variable must be a matrix with two columns, where the first column is the number of "successes" and the second column is the number of "failures".  
When several models have been estimated, methods \code{summary} and \code{print} consider the best model according to the information criterion in \code{criterion}, among the estimated models having a number of components among those in \code{k} an error distribution among those in \code{familyY} and a parsimonious model among those in \code{modelXnorm}.
}
\value{
This function returns a class \code{cwm} object, which is a list of values related to the model selected. It contains:

\item{call}{an object of class \code{call}.}
\item{formulaY}{an object of class \code{\link[stats:formula]{formula}} containing a symbolic description of the model fitted.}
\item{familyY}{the distribution used for the conditional distribution of \eqn{Y} in each mixture component.}
\item{data}{a \code{\link[base:data.frame]{data.frame}} with the variables needed to use \code{formulaY}.}
\item{concomitant}{a list containing \code{Xnorm}, \code{Xbin}, \code{Xpois}, \code{Xmult}.}
\item{\code{Xbtrials}}{number of trials used for \code{Xbin}.}
\item{models}{a list; each element is related to one of the models fitted. Each element is a list and contains:}
\itemize{
  \item{\code{posterior}} {posterior probabilities}
  \item{\code{iter}} {number of iterations performed in EM algorithm}
  \item{\code{k }} {number of (fitted) mixture components.}
  \item{\code{size}} {estimated size of the groups.}
  \item{\code{cluster }} {classification vector}
  \item{\code{loglik }} {final log-likelihood value}
  \item{\code{df }} {overall number of estimated parameters}
  \item{\code{prior }} {weights for the mixture components}
  \item{\code{IC }} {list containing values of the information criteria }
  \item{\code{converged}} {logical; \code{TRUE} if EM algorithm converged}
  \item{\code{GLModels}} {a list;  each element is related to a mixture component and contains:}
  \itemize{
    \item{\code{model }} {a "\code{\link[stats:glm]{glm}}" class object.}
    \item{\code{sigma }} {estimated local scale parameters of the conditional distribution of \eqn{Y}, when \code{familyY} is \code{gaussian} or \code{student.t}}
    \item{\code{t_df}} {estimated degrees of freedom of the t distribution, when \code{familyY} is \code{student.t}}
    \item{\code{nuY}} {estimated shape parameter, when \code{familyY} is \code{Gamma}. The gamma distribution is parameterized according to McCullagh & Nelder (1989, p. 30)}
  }
  \item{\code{concomitant} {a list with estimated concomitant variables parameters for each mixture component}}
  \itemize{
    \item{\code{normal.d, multinomial.d, poisson.d, binomial.d}} {marginal distribution of concomitant variables}
    \item{\code{normal.mu}} {mixture component means for \code{Xnorm}  }
    \item{\code{normal.Sigma}} {mixture component covariance matrices for \code{Xnorm}}
    \item{\code{normal.model}} {models fitted for \code{Xnorm} } 
    \item{\code{multinomial.probs}} {multinomial distribution probabilities for  \code{Xmult}} 
    \item{\code{poisson.lambda}} {lambda parameters for \code{Xpois} } 
    \item{\code{binomial.p}} {binomial probabilities for \code{Xbin} } 
  }
  }
}

\references{

Mazza, A., Ingrassia, S., and Punzo, A. (2018). {flexCWM}: A Flexible Framework for Cluster-Weighted Models.  \emph{Journal of Statistical Software}, \bold{86}(2), 1-30.

Ingrassia, S., Minotti, S. C., and Vittadini, G. (2012). Local Statistical Modeling via the Cluster-Weighted Approach with Elliptical Distributions. \emph{Journal of Classification}, \bold{29}(3), 363-401.

Ingrassia, S., Minotti, S. C., and Punzo, A. (2014). Model-based clustering via linear cluster-weighted models. \emph{Computational Statistics and Data Analysis}, \bold{71}, 159-182.

Ingrassia, S., Punzo, A., and Vittadini, G. (2015). The Generalized Linear Mixed Cluster-Weighted Model. \emph{Journal of Classification}, \bold{32}(forthcoming)

McCullagh, P. and Nelder, J. (1989). Generalized Linear Models. Chapman & Hall, Boca Raton, 2nd edition

Punzo, A. (2014). Flexible Mixture Modeling with the Polynomial Gaussian Cluster-Weighted Model. \emph{Statistical Modelling}, \bold{14}(3), 257-291.

%Punzo, A. and Ingrassia, S. (2015a). Clustering Bivariate Mixed-Type Data via the Cluster-Weighted Model. Computational Statistics, ..(..), ..-..

%Punzo, A. and Ingrassia, S. (2015b). Parsimonious Generalized Linear Gaussian Cluster-Weighted Models. In T. Minerva, I. Morlini, and F. Palumbo (eds.), \emph{Statistical Models for Data Analysis}, Studies in Classification, Data Analysis and Knowledge Organization, pp. ..-.. Springer International Publishing, Switzerland.
}
\author{
Mazza A., Punzo A., Ingrassia S.
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{                  
\code{\link{flexCWM-package}}
}
\examples{
## an exemple with artificial data
data("ExCWM")
attach(ExCWM)
str(ExCWM)

# mixtures of binomial distributions
resXbin <- cwm(Xbin = Xbin, k = 1:2, initialization = "kmeans")
getParXbin(resXbin)

# Mixtures of Poisson distributions
resXpois <- cwm(Xpois = Xpois, k = 1:2, initialization = "kmeans")
getParXpois(resXpois)

# parsimonious mixtures of multivariate normal distributions
resXnorm <- cwm(Xnorm = cbind(Xnorm1,Xnorm2), k = 1:2, initialization = "kmeans")
getParXnorm(resXnorm)

## an exemple with real data
data("students")
attach(students)
str(students)
# CWM
fit2 <- cwm(WEIGHT ~ HEIGHT + HEIGHT.F , Xnorm = cbind(HEIGHT, HEIGHT.F), 
  k = 2, initialization = "kmeans", modelXnorm = "EEE")
summary(fit2, concomitant = TRUE)
plot(fit2)


}

