############################################
## Classification matrix given the groups ##
############################################

class <- function(groups,k){
  
  n <- length(groups)
  z <- array(0,c(n,k),dimnames=list(1:n,paste("comp",1:k,sep="")))
  for(i in 1:n) 
    z[i,groups[i]] <- 1
  return(z)
  
}

########################################
## Converter for the categorical data ##
########################################

converter <- function(X,m=NULL){
  
  # X: matrix of integers (having the variables by column)
  # m: number of levels for each categorical variable (starting by 1)
  
  if(is.null(X))     stop('Hey, we need some data, please! X is null')
  if(!is.matrix(X))  stop('X needs to be in matrix form')
  #if(!is.numeric(X)) stop('X is required to be numeric')
  
  d <- ncol(X)
  n <- nrow(X)
  
  if(is.null(m)) m <- apply(X,2,max)
  
  binary <- NULL
  
  for(u in 1:d)
    binary[[u]] <- class(X[,u],m[u]) 
  
  return(
    list(
      m = m,
      binary = binary,
      d = d,
      n = n
    )
  )
  
}

#########################################
## Re-parameterized Gamma distribution ##
#########################################

dgam <- function(x,mu,nu,log=FALSE){
  
  dgamma(x, shape=nu, scale = mu/nu, log = log)
  
}

rgam <- function(x,mu,nu){
  
  rgamma(x, shape=nu, scale = mu/nu)
  
}

################################
## Weighted Covariance Matrix ##
################################

mstep <- function(x,wt){
  
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  else if (!is.matrix(x)) 
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x))) 
    stop("'x' must contain finite values only")
  p   <- ncol(x)
  n   <- nrow(x)
  mu  <- array(colSums(wt/sum(wt)*x),c(p),dimnames=list(paste("X.",1:p,sep="")))
  cov <- array(crossprod(sqrt(wt/sum(wt))*(x-matrix(rep(mu,n),n,p,byrow=TRUE))),c(p,p),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep="")))
  return(
    list(
      center = mu,
      cov = cov      
    )
  )
}

#########################
## Multinomial case: X ##
#########################

# MultX <- function(X,m=NULL,weights){ #,method="Nelder-Mead",start=FALSE,initial.values=NULL,maxit=1000,reltol=1e-15,trace=1
#   
#   # X: matrix of integers (having the variables by column)
#   # m: number of levels for each categorical variable (starting by 1)
#   # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
#   
#   temp <- converter(X=X,m=m)
#   
#   n      <- temp$n
#   d      <- temp$d
#   m      <- temp$m
#   binary <- temp$binary
#   
#   # binary matrix creation
#   
#   Xmod <- binary[[1]]
#   if(d>1){
#     for(j in 2:d)
#       Xmod <- cbind(Xmod,binary[[j]])
#   }  
#   
#   alpha.obs <- colMeans(m)
#   
#   alpha.obs  <- weighted.mean(x=X,w=weights)
#   
#   f <- function(par,X,weights){
#     
#     lambda <- par
#     
#     l <- -sum(weights*dpois(x=X, lambda = lambda, log = TRUE))
#     
#   }
#   
#   res <- optimize(f=f, interval=c(0,2*lambda.obs), X=X, weights=weights)
#   
#   loglik     <- -res$objective
#   lambda.hat<-  res$minimum
#   
#   return(
#     list(
#       loglik  = loglik,
#       lambdaX = lambda.hat
#     )
#   )  
#   
# }

########################
## Gaussian case: Y|x ## 
########################

GaussianY <- function(Y,X,weights,method="Nelder-Mead",start=FALSE,initial.values=NULL,maxit=1000,reltol=1e-15,trace=0){ 
  
  # Y: a numerical vector for the response variable
  # X: a matrix for the covariates
  # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
  
  p <- ncol(X) 
  
  if(start==FALSE){
    
    modelY      <- glm(Y ~ X,family=gaussian(),weights=weights,control=list(trace=FALSE,epsilon=1e-14)) 
    beta.obs    <- modelY$coefficients
    sigma.obs   <- sqrt(summary(modelY)$dispersion)
    
    starting.values <- c(beta.obs,log(sigma.obs))
    
  }
  
  f <- function(par,Y,X,weights){
    
    beta0 <- par[1]
    beta1 <- par[2:(p+1)]
    sigma <- exp(par[2+p])
    
    l <- -sum(weights*dnorm(x=Y, mean = beta0 + X %*% beta1, sd = sigma, log = TRUE))
    
  }
  
  if(start==FALSE)
    res <- optim(par=starting.values, fn=f, Y=Y, X=X, weights=weights, method=method, control = list(maxit=maxit,reltol=reltol,trace=trace))
  if(start==TRUE)
    res <- optim(par=initial.values, fn=f, Y=Y,X=X, weights=weights, method=method, control = list(maxit=maxit,reltol=reltol,trace=trace))
  
  loglik  <- -res$value
  est     <-  res$par
  
  beta0.hat <- est[1]
  beta1.hat <- est[2:(p+1)]
  sigma.hat <- exp(est[2+p])
  
  return(
    list(
      loglik = loglik,
      beta0  = beta0.hat,
      beta1  = beta1.hat,
      sigmaY = sigma.hat
    )
  )  
  
}

#######################
## Poisson case: Y|x ## 
#######################

PoissonY <- function(Y,X,weights,method="Nelder-Mead",start=FALSE,initial.values=NULL,maxit=1000,reltol=1e-15,trace=0){ 
  
  # Y: a numerical vector for the response variable
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
  
  p <- ncol(X) 
  
  if(start==FALSE){
    
    modelY      <- glm(Y ~ X,family=poisson(),weights=weights,control=list(trace=FALSE,epsilon=1e-14)) 
    beta.obs    <- modelY$coefficients
    
    starting.values <- c(beta.obs)
    
  }
  
  f <- function(par,Y,X,weights){
    
    beta0 <- par[1]
    beta1 <- par[2:(p+1)]
    
    l <- -sum(weights*dpois(x=Y, lambda = exp(beta0 + X %*% beta1), log = TRUE))
    
  }
  
  if(start==FALSE)
    res <- optim(par=starting.values, fn=f, Y=Y, X=X, weights=weights, method=method, control = list(maxit=maxit,reltol=reltol,trace=trace))
  if(start==TRUE)
    res <- optim(par=initial.values, fn=f, Y=Y, X=X, weights=weights, method=method, control = list(maxit=maxit,reltol=reltol,trace=trace))
  
  loglik  <- -res$value
  est     <-  res$par
  
  beta0.hat <- est[1]
  beta1.hat <- est[2:(p+1)]
  
  return(
    list(
      loglik = loglik,
      beta0  = beta0.hat,
      beta1  = beta1.hat
    )
  )  
  
}

########################
## Binomial case: Y|x ## 
########################

BinomialY <- function(Y,X,weights,m,method="Nelder-Mead",start=FALSE,initial.values=NULL,maxit=1000,reltol=1e-15,trace=0){ 
  
  # Y: a numerical vector for the response variable
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
 
  p <- ncol(X)
  
  if(start==FALSE){
    
    modelY   <- glm(cbind(Y,m-Y) ~ X,family=binomial(),weights=weights,control=list(trace=FALSE,epsilon=1e-14)) 
    beta.obs <- modelY$coefficients
    
    starting.values <- c(beta.obs)
    
  }
  
  f <- function(par,Y,X,weights,m){
    
    beta0 <- par[1]
    beta1 <- par[2:(p+1)]
    
    l <- -sum(weights*dbinom(x=Y, size = m, prob = exp(beta0 + X %*% beta1)/(1+exp(beta0 + X %*% beta1)), log = TRUE))
    
  }
  
  if(start==FALSE)
    res <- optim(par=starting.values, fn=f, Y=Y, X=X, weights=weights, m=m, method=method, control = list(maxit=maxit,reltol=reltol,trace=trace))
  if(start==TRUE)
    res <- optim(par=initial.values, fn=f, Y=Y, X=X, weights=weights, m=m, method=method, control = list(maxit=maxit,reltol=reltol,trace=trace))
  
  loglik  <- -res$value
  est     <-  res$par
  
  beta0.hat <- est[1]
  beta1.hat <- est[2:(p+1)]
  
  return(
    list(
      loglik = loglik,
      beta0  = beta0.hat,
      beta1  = beta1.hat
    )
  )  
  
}

#####################
## Gamma case: Y|x ## 
#####################

GammaY <- function(Y,X,weights,method="Nelder-Mead",start=FALSE,initial.values=NULL,maxit=1000,reltol=1e-15,trace=1){ 
  
  # Y: a numerical vector for the response variable
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
  
  p <- ncol(X)
  
  if(start==FALSE){
    
    modelY    <- glm(Y ~ X,family=Gamma(link="log"),weights=weights,control=list(trace=FALSE,epsilon=1e-14)) 
    beta.obs  <- modelY$coefficients
    nu.obs    <- 1/summary(modelY)$dispersion
    
    starting.values <- c(beta.obs,log(nu.obs))
    
  }
  
  f <- function(par,Y,X,weights){
    
    beta0 <- par[1]
    beta1 <- par[2:(p+1)]
    nu    <- exp(par[2+p])
    
    l <- -sum(weights*dgam(x=Y, mu = exp(beta0 + X %*% beta1), nu = nu, log = TRUE))
    
  }
  
  if(start==FALSE)
    res <- optim(par=starting.values, fn=f, Y=Y, X=X, weights=weights, method=method, control = list(maxit=maxit,reltol=reltol,trace=trace))
  if(start==TRUE)
    res <- optim(par=initial.values, fn=f, Y=Y, X=X, weights=weights, method=method, control = list(maxit=maxit,reltol=reltol,trace=trace))
  
  loglik  <- -res$value
  est     <-  res$par
  
  beta0.hat <- est[1]
  beta1.hat <- est[2:(p+1)]
  nu.hat    <- exp(est[2+p])
  
  return(
    list(
      loglik = loglik,
      beta0  = beta0.hat,
      beta1  = beta1.hat,
      nuY    = nu.hat
    )
  )  
  
}

#####################
## Generalized CWM ##
#####################

glmcwm2 <- function(
	Y,			                 # numerical vector for the response variable
  Xcont=NULL,			                 # (nxp) matrix for the continuous covariates
	Xcate=NULL,  		                 # (nxd) matrix for the categorical covariates
  m=NULL,                  # m: number of levels for each categorical variable in Xcate (starting by 1)
	familyY = "Gaussian",    # the exponential distribution used for Y|x 
	k=2,                     # number of groups
	mY=1,                         # mY=1 for Bernoulli; mY>1 for Binomial
	method="Nelder-Mead",         # optimization method
  initialization="random.soft", # initialization procedure: "random.soft", "random.hard", and "manual" 
	start.z=NULL,                 # (n x k)-matrix of soft or hard classification: it is used only if initialization="manual"		
	iter.max=200,                # maximum number of iterations in the EM-algorithm
	threshold=1.0e-04,            # stopping rule in the Aitken procedure
	loglikplot=FALSE,               # if TRUE, the log-likelihood values against the iterations are plotted
	seed = NULL
  )
{
  
  if(is.null(Xcont) & is.null(Xcate)) 
    stop('Hey, we need some data, please!')
  if(!is.numeric(Xcont) & !is.null(Xcont)) 
    stop('Xcont is required to be numeric')
  if(!is.numeric(Xcate) & !is.null(Xcate)) 
    stop('Xcate is required to be numeric')
  n <- length(Y)    # sample size
  if(n == 1)   
    stop('Hey, we need more than one observation, please!')
  #if(any(is.na(Xcont)) & !is.null(Xcont))  
    #stop('No NAs allowed in Xcont.')
  #if(any(is.na(Xcate)) & !is.null(Xcate))  
    #stop('No NAs allowed in Xcate.')
  if(is.null(k)) 
    stop('k is NULL')
  k <- as.integer(ceiling(k))
  if(!is.integer(k)) 
    stop('k is not a integer')
  if(any(k < 1)) 
    stop('k is not a positive integer')
  
  if(!is.null(Xcont))
    Xcont <- as.matrix(Xcont)
  if(!is.null(Xcate))
    Xcate <- as.matrix(Xcate)
  X <- cbind(Xcont,Xcate)
  p <- 0            # number of variables for the continuous part of X
  d <- 0            # number of variables for the categorical part of X
  if(!is.null(Xcont))
    p <- ncol(Xcont)  
  if(!is.null(Xcate))
    d <- ncol(Xcate)  
 
  # modification of Xcate in order to apply the multinomial model
  
  if(!is.null(Xcate)){
    temp <- converter(X=Xcate,m=m)
    Xmod <- temp$binary
    m    <- temp$m
  }
  
# parameters definition
  
prior   <- numeric(k) # weights

# Xcont

muX <- VarX <- invVarX <- NULL
if(!is.null(Xcont)){
  muX     <- array(0,c(p,k),dimnames=list(paste("X.",1:p,sep=""),paste("comp.",1:k,sep="")))
  VarX    <- array(0,c(p,p,k),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("comp.",1:k,sep="")))
  invVarX <- array(0,c(p,p,k),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("comp.",1:k,sep="")))
}
PXcont  <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))

# Xcate

alpha <- NULL
if(!is.null(Xcate)){
  for(u in 1:d)
    alpha[[u]] <- array(0,c(m[u],k),dimnames=list(paste("level.",1:m[u],sep=""),paste("comp.",1:k,sep="")))
}
PXcate <- array(1,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  
# Y|x

beta0   <- numeric(k)
beta1   <- array(0,c(p+d,k),dimnames=list(paste("X.",1:(p+d),sep=""),paste("comp.",1:k,sep="")))
dispY   <- numeric(k)   # Conditional Dispersion for each group
nuY     <- NULL         # only for the Gamma distribution   
VarFunY <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
VarY    <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
muY     <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
PY      <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))

# ------------------------- #
# posteriors initialization #
# ------------------------- #

# RANDOM INITIALIZATION

if(initialization=="kmeans"){
  
  D         <- cbind(Y,X)                # the complete set of data
  clusters  <- kmeans(x=D,centers=k)               # clusters on D
  z         <- class(clusters$cluster,k)
  
  } 

if(initialization=="random.soft"){

  if(!is.null(seed)) set.seed(seed)
	z  <- array(runif(n*k),c(n,k)) # soft posterior probabilities (no-normalized) (n x k) 
	z  <- z/rowSums(z)             # soft posterior probabilities (n x k)
  
  } 

if(initialization=="random.hard"){

  if(!is.null(seed)) set.seed(seed)
  z  <- t(rmultinom(n, size = 1, prob=rep(1/k,k)))  # hard posterior probabilities (n x k)
  
  } 

if(initialization=="manual"){ # z.start can be both soft and hard initialization
  
  z  <- start.z      # posterior probabilities (n x k) no-normalized

  } 
  
# ------------ #
# EM algorithm #
# ------------ #

# Preliminary definition of convergence criterions

check     <- 0
iteration <- 1
loglik    <- NULL
aloglik   <- NULL
aloglik   <- c(0,0)
a         <- NULL
a         <- c(0,0)

while(check<1){

	# ++++++ #
	# M-step #
	# ++++++ #

  # ------- #
  # Weights #
  # ------- #

	prior <- colMeans(z)
  nj    <- colSums(z)

  # ----- #
	# Xcont #
	# ----- #
  
	if(!is.null(Xcont)){
    for(j in 1:k){ 
	  temp          <- mstep(x=Xcont,wt=z[,j])
	  muX[,j]       <- temp$center
	  VarX[,,j]     <- temp$cov
	  invVarX[,,j]  <- solve(temp$cov)
	  PXcont[,j]    <- (2*pi)^(-p/2)*( ifelse(p>1,det(VarX[,,j]),VarX[,,j]) )^(-1/2)*exp(-1/2*mahalanobis(x=Xcont, center=muX[,j], cov=invVarX[,,j], inverted=TRUE))                                   
	  }
	}
    
	# ----- #
	# Xcate #
	# ----- #
	
  if(!is.null(Xcate)){
    for(u in 1:d)
      alpha[[u]] <- (t(Xmod[[u]]) %*% z)/matrix(rep(nj,m[u]),nrow=m[u],byrow=T)
    for(j in 1:k)
      for(i in 1:n){
        PXcate[i,j] <- 1
        for(u in 1:d)
          PXcate[i,j] <- PXcate[i,j]*dmultinom(Xmod[[u]][i,],size=1,prob=alpha[[u]][,j])
      }
  }
  
  # --- #
	# Y|x #
	# --- #
  
  for(h in 1:k){
	  if(familyY=="Gaussian"){
      modelY      <- GaussianY(Y,X,weights=z[,h],method=method) #,mustart=muY[,h]
      beta0[h]    <- modelY$beta0
      beta1[,h]   <- modelY$beta1
      muY[,h]     <- beta0[h] + X %*% beta1[,h] 
      dispY[h]    <- modelY$sigmaY^2
      VarFunY[,h] <- rep(1,n)     
      VarY[,h]    <- dispY[h]*VarFunY[,h]     
      PY[,h]      <- dnorm(Y,mean=muY[,h],sd=sqrt(VarY[,h]))              
	  }
	}
	for(h in 1:k){
	  if(familyY=="Poisson"){
	    modelY      <- PoissonY(Y,X,weights=z[,h],method=method) #,mustart=muY[,h]
	    beta0[h]    <- modelY$beta0
	    beta1[,h]   <- modelY$beta1
	    muY[,h]     <- exp(beta0[h] + X %*% beta1[,h])
	    dispY[h]    <- 1
	    VarFunY[,h] <- muY[,h]     
	    VarY[,h]    <- dispY[h]*VarFunY[,h]     
	    PY[,h]      <- dpois(Y,lambda=muY[,h])              
	  }
	}
	for(h in 1:k){
	  if(familyY=="Binomial"){
	    modelY      <- BinomialY(Y,X,weights=z[,h],m=mY,method=method) #,mustart=muY[,h]
	    beta0[h]    <- modelY$beta0
	    beta1[,h]   <- modelY$beta1
	    muY[,h]     <- exp(beta0[h] + X %*% beta1[,h])/(1+exp(beta0[h] + X %*% beta1[,h]))
	    dispY[h]    <- 1/mY
	    VarFunY[,h] <- muY[,h]*(1-muY[,h])     
	    VarY[,h]    <- dispY[h]*VarFunY[,h]     
	    PY[,h]      <- dbinom(Y,size=mY,prob=muY[,h])              
	  }
	}
	for(h in 1:k){
	  if(familyY=="Gamma"){
	    modelY      <- GammaY(Y,X,weights=z[,h],method=method) #,mustart=muY[,h]
	    beta0[h]    <- modelY$beta0
	    beta1[,h]   <- modelY$beta1
	    muY[,h]     <- exp(beta0[h] + X %*% beta1[,h])
	    dispY[h]    <- 1/modelY$nuY
	    VarFunY[,h] <- muY[,h]^2     
	    VarY[,h]    <- dispY[h]*VarFunY[,h]     
	    nuY[h]      <- modelY$nuY
      PY[,h]      <- dgam(Y, mu = muY[,h], nu = nuY[h])              
	  }
	}
  
  # ------------------------------------- # 
  # Global - Observed-data log-likelihood # 
  # ------------------------------------- #
  
	llvalues          <- log(rowSums(matrix(rep(prior,n),n,k,byrow=TRUE)*PY*PXcont*PXcate))
  loglik[iteration] <- sum(llvalues)
  
  # ----------------------------------------------- #
	# Aitkane's Acceleration-Based Stopping Criterion #
	# ----------------------------------------------- #
	
	if(iteration > 2 & k > 1){
	  a[iteration-1]      <- (loglik[iteration]-loglik[iteration-1])/(loglik[iteration-1]-loglik[iteration-2])
	  aloglik[iteration]  <- loglik[iteration-1]+(1/(1-a[iteration-1])*(loglik[iteration]-loglik[iteration-1]))
	  if(abs(aloglik[iteration]-loglik[iteration]) < threshold) 
	    check <- 1 
	}
	
	if(iteration==iter.max | k==1) check <- 1
  
	cat("*")
	iteration <- iteration + 1
  
	# ++++++ #
	# E-Step #
	# ++++++ #
	
	z.num  <- matrix(rep(prior,n),n,k,byrow=TRUE)*PY*PXcont*PXcate  # (n x k)
	z.den  <- rowSums(z.num)                                        # n-vector
	z      <- z.num/matrix(rep(z.den,k),ncol=k)                     # (n x k)
  
}

finalloglik <- loglik[iteration-1] 
	
# ----------------------------------------------------------------------- #
# The EM-algorithm is finished                                            #
# Check on the EM-monotonicity                                            #
# Plotting the values of the observed loglikelihood versus the iterations #
# ----------------------------------------------------------------------- #

cat("\n")
if(loglikplot==TRUE & k>1){

  par(mai=c(0.84,0.8,0.012,0.004))
  par(las = 3)
  par(cex.axis=0.7)
  par(cex.lab=1.2)
  plot(0:(iteration-2),loglik[1:(iteration-1)],type="l",axes = FALSE,xlab="iterations",ylab="log-likelihood",lwd=2)
  axis(1, at = 0:(iteration-2),labels = 0:(iteration-2)) 
  axis(2)
  box(col = "black")

}

# --------------------- #
# Classification Matrix #
# --------------------- #

group <- apply(z,1,which.max)  

##########################
## Information criteria ##
##########################
    
if(familyY=="Gaussian" | familyY=="Gamma")  
  npar <- (k-1) + p*k + k*(p*(p+1)/2) + k*(sum(m)-d) + (p+d+2)*k 
if(familyY=="Binomial" | familyY=="Poisson")  
  npar <- (k-1) + p*k + k*(p*(p+1)/2) + k*(sum(m)-d) + (p+d+1)*k
    
AIC       <- 2*finalloglik - npar*2
BIC       <- 2*finalloglik - npar*log(n)
AIC3      <- 2*finalloglik - npar*3  
AICc      <- AIC - (2*npar*(npar+1))/(n-npar-1)
AICu      <- AICc - n*log(n/(n-npar-1))
CAIC      <- 2*finalloglik - npar*(1+log(n))
AWE       <- 2*finalloglik - 2*npar*(3/2+log(n))   

z.const    <- (z<10^(-322))*10^(-322)+(z>10^(-322))*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
hard.z     <- (matrix(rep(apply(z,1,max),k),n,k,byrow=F)==z)*1
ECM        <- sum(hard.z*log(z.const))
ICL        <- BIC+ECM

beta <- rbind(beta0,beta1)
dimnames(beta)[[1]] <- paste("beta",0:(p+d),sep="")
dimnames(beta)[[2]] <- paste("comp",1:k,sep="")

result <- list(
  Y         = Y,            
  familyY   = familyY,       
  X         = X,            
  p         = p,
  d         = d,
  k         = k,            
  n         = n,
  m         = m,
  npar      = npar,
  alpha     = alpha,
  mY        = mY,           # mY=1 for Bernoulli; mY>1 for Binomial
  prior     = prior,
  muX       = muX,
  VarX      = VarX,
  PXcont    = PXcont,
  PXcate    = PXcate,
  beta      = beta,
  muY       = muY,
  dispY     = dispY,
  VarFunY   = VarFunY,
  VarY      = VarY,
  nuY       = nuY,
  PY        = PY,
  iter.stop = iteration,
  z         = z,
  group     = group,
  loglik    = finalloglik,
  AIC       = AIC,
  AIC3      = AIC3,
  AICc      = AICc,
  AICu      = AICu,
  AWE       = AWE,
  CAIC      = CAIC,
  BIC       = BIC,
  ICL       = ICL,      # alla McNicholas
  call      = match.call()
)

class(result) <- "glmcwm"
return(result)

alarm()

}

#####################
## model Selection ##
#####################

glmcwm <- function(
  Y,  		                 # numerical vector for the response variable
  Xcont=NULL,			                 # (nxp) matrix for the continuous covariates
  Xcate=NULL,  		                 # (nxd) matrix for the categorical covariates
  m=NULL,                  # m: number of levels for each categorical variable in Xcate (starting by 1)
  familyY = "Gaussian",    # the exponential distribution used for Y|x 
  k=2,                  # minimum number of groups
  mY=1,                         # mY=1 for Bernoulli; mY>1 for Binomial
  method="Nelder-Mead",         # optimization method
  initialization="random.soft", # initialization procedure: "random.soft", "random.hard", and "manua" 
  start.z=NULL,                 # (n x k)-matrix of soft or hard classification: it is used only if initialization="manual"		
  iter.max=1000,                # maximum number of iterations in the EM-algorithm
  threshold=1.0e-04,             # stopping rule in the Aitken rule
  loglikplot=FALSE,
  seed=NULL
)              
{
  call=match.call()
  
  n     <- length(Y)
  gridk <- k
  numk  <- length(gridk)
  
  IC <- array(0,c(numk,9),dimnames=list(paste(gridk,"groups",sep=" "),c("loglik","AIC","AICc","AICu","AIC3","AWE","BIC","CAIC","ICL")))

  bestic <- -1e100
  par <- list()
  for(i in 1:numk){
    par[[i]] <- glmcwm2(
      Y=Y,    	                 
      Xcont=Xcont,			                 
      Xcate=Xcate,  		                 
      m=m,                 
      familyY = familyY,     
      k=gridk[i],                            
      mY=mY,                         
      method=method,                
      initialization=initialization,  
      start.z=start.z,                 		
      iter.max=iter.max,                
      threshold=threshold,            
      loglikplot=FALSE,
      seed=seed
    )
  }
  ic <- c("AIC", "AICc", "AICu", "AIC3", "AWE", "BIC", "CAIC", "ICL")
  cr <- t(sapply(par,function(x) sapply(parse(text=paste0("x$",ic)),function(g) eval(g) )))
  best <- apply(cr,2,function(x) c(which.max(x),max(x)))
  colnames(best) <- ic
  rownames(best) <- c("n. of groups","value")
  
  cat("\n\n")
  cat("# ----------------------- #","\n")
  cat("# Model Selection Results #","\n")
  cat("# ----------------------- #","\n\n")
  print(best)
  
  
  bestpar <- apply(best,2,function(x) par[[x[1]]])
  
   return(
    list(best=best,bestpar=bestpar,par=par)
  )  
  
}