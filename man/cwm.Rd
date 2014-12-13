\name{cwm}
\alias{cwm}
\alias{summary.cwm}
\alias{print.cwm}
\alias{modelget}
\alias{ICget}
\alias{bestmodel}
\alias{print.cwm.IC}
\title{Fit for the Generalized Linear Mixed CWM}
\description{Maximum likelihood fitting of the generalized linear mixed  cluster-weighted model by the EM algorithm.}
\usage{
cwm(formulaY, familyY=gaussian, data,link, Xnorm=NULL, modelXnorm=NULL, Xbin=NULL,
  Xbtrials=NULL, Xpois=NULL, Xmult=NULL, k=1:3, initialization=c("random.soft", 
  "random.hard", "kmeans", "mclust", "manual"), start.z=NULL, seed=NULL, maxR=1,
  iter.max=1000, threshold=1.0e-04, parallel=FALSE)

ICget(object, criteria)
bestmodel(object, criteria, k=NULL, modelXnorm=NULL)
modelget(object, criteria="BIC", k=NULL, modelXnorm=NULL)
\method{summary}{cwm}(object, criteria="BIC", k=NULL, modelXnorm=NULL, concomitant=FALSE, 
  digits = getOption("digits")-2, \dots)

}
                            
%- maybe also 'usage' for other objects documented here.
\arguments{                      
  \item{formulaY}{an optional object of class "\code{\link[stats:formula]{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.    
  }
    \item{ familyY}{
  the distribution used for \eqn{Y|x} in each mixture component; it can be: 
    \itemize{
      \item \code{"gaussian"} with default \code{"link=identity"}
      \item \code{"poisson"} with default \code{"link=log"} 
      \item \code{"binomial"} with default \code{"link=logit"} 
      \item \code{"Gamma"} with default \code{"link=log"} 
      \item \code{"inverse.gaussian"} with default \code{"link=log"} 
      \item \code{"t"} with default \code{"link=identity"} 
    }
  Default value is \code{"gaussian"}.
  }
  \item{data}{
  an optional \code{\link[base:data.frame]{data.frame}}, \code{\link[base:list]{list}}, or \code{\link[base:environment]{environment}} with the variables needed to use \code{formulaY}.
  }
    \item{link}{
  a specification for the model link function to be used. See \code{link} argument in \code{\link[stats:family]{family}}.   
  }
  \item{Xnorm, Xbin, Xpois, Xmult}{
  an optional matrix containing variables to be used for marginalization having normal, binomial, Poisson, multinomial distributions. 
  }
  \item{modelXnorm}{an optional vector of character strings indicating the parsimonious models to be fitted. The default is  \code{c("E", "V")} for a single continuous covariate, and \code{c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")} for multivariate continuous covariates. For details, see \code{\link[mixture:gpcm]{mixture:gpcm}}. 
  }
  \item{Xbtrials}{
   an optional vector containing the number of trials for each column in \code{Xbin}. If omitted, the maximum of each column in \code{Xbin} is chosen.
   }
  \item{k}{
  an optional vector containing the numbers of mixture components to be tried. Default value is 1:3.
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
  matrix of soft or hard classification: it is used only if \code{initialization="manual"}.
}
  \item{seed}{
 an optional scalar. It sets the seed for the random number generator, when random initializations are used; if \code{NULL}, current seed is not changed. Default value is \code{NULL}.
}

 \item{maxR}{
 number of initializations to be tried.Default value is 1.
  } 
  \item{iter.max}{
   an optional scalar. It sets the maximum number of iterations in the EM-algorithm. Default value is 200.
}
  \item{threshold}{
  an optional scalar. It sets the threshold for Aitken acceleration procedure. Default value is 1.0e-04.
}
  \item{parallel}{
  When \code{TRUE}, the package \code{\link[parallel:parallel-package]{parallel}} is used for parallel computation. When several models are estimated, computational time is reduced. The number of cores to use may be set with the global option \code{cl.cores}; default value is detected using \code{\link[parallel:detectCores]{detectCores()}}.
}
  
  \item{object}{a class \code{cwm} object.
}
  \item{concomitant}{When \code{TRUE}, concomitant variables parameters are displayed. Default is \code{FALSE}
}
  \item{digits}{ 
        integer used for number formatting.
}
  \item{criteria}{an optional character string. It sets the information criteria to consider; supported values are: \code{"AIC", "AICc", "AICu", "AIC3", "AWE", "BIC", "CAIC", "ICL".} Default value is \code{"BIC"}.
}
  \item{\dots}{   
  additional arguments affecting the summary produced. 
}
}
\details{When several models have been estimated, methods \code{bestmodel}, \code{summary} and \code{print} consider the models with the best information criteria in \code{criteria}, among those with  \code{k} groups and \code{modelXnorm} parsimonious model. If \code{criteria} is missing, the model with best BIC is returned.
The \code{modelget} method returns a \code{cwm} object containing the best model according to a single criterion in  .
}
\value{
This function returns a class \code{cwm} object, which is a list of values related to the model selected. It contains:

\item{call}{an object of class \code{call}.}
\item{formulaY}{an object of class \code{\link[stats:formula]{formula}} containing a symbolic description of the model fitted.}
\item{ familyY}{the distribution used for \eqn{Y|x} in each mixture component.}
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
  \item{\code{GLModels}} {a list;  each element is related to a mixture component. Each element is a list and contains:}
  \itemize{
    \item{\code{model }} {a "\code{\link[stats:glm]{glm}}" class object.}
    \item{\code{sigma }} {estimated local scale parameters of \eqn{Y|x}, when \code{familyY} is \code{gaussian} or \code{t}}
    \item{\code{t_df}} {estimated degrees of freedom of the t distribution, when \code{familyY} is \code{t}}
    \item{\code{nuY}} {estimated shape parameter, when \code{familyY} is \code{Gamma}. The gamma distribution is parameterized according to McCullagh, P. and Nelder, J. 1989, p. 30}
  }
  \item{\code{concomitant} {a list with estimated concomitant variables parameters for each mixture component}}
  \itemize{
    \item{\code{normal.dnorm, multinomial.dmulti, poisson.dpois, binomial.dbin}} {marginal distribution of concomitant variables}
    \item{\code{normal.mu}} {mixture component means for \code{Xnorm}  }
    \item{\code{normal.Sigma}} {mixture component covariance matrices for \code{Xnorm}}
    \item{\code{normal.model}} {models fitted for \code{Xnorm} } 
    \item{\code{multinomial.model}} {model fitted for \code{Xmult} } 
    \item{\code{multinomial.probs}} {multinomial distribution probabilities for  \code{Xmult}} 
    \item{\code{poisson.lambda}} {lambda parameters for \code{Xpois} } 
    \item{\code{binomial.p}} {binomial probability for \code{Xbin} } 
  }
  }
}

\references{

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
\code{\link{flexCWM-package}}, \code{\link{tourism}}
}
\examples{
data("students")
str(students)
attach(students)

# mixture of Gaussian distributions
res <- cwm(Xnorm=HEIGHT, k=1:3, initialization="kmeans")
summary(res)
plot(res)

# mixture of Gaussian regressions
res2 <- cwm(HEIGHT ~ HEIGHT.F, k=1:3, initialization="mclust")
summary(res2)
plot(res2)

%# CWM

%res3 <- cwm(WEIGHT ~ HEIGHT + HEIGHT.F, familyY="t", Xnorm=cbind(HEIGHT,HEIGHT.F), 
%            modelXnorm="VVV", k=1:2, initialization="mclust")
%summary(res3,concomitant=TRUE)
%plot(res3)
}

