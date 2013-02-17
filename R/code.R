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

PoissonY <- function(Y,X,weights,method="Nelder-Mead",start=FALSE,initial.values=NULL,maxit=1000,reltol=1e-15,trace=1){ 
  
  # Y: a numerical vector for the response variable
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
  
  p <- ncol(X) 
  
  if(start==FALSE){
    
    modelY      <- glm(Y ~ X,family=poisson(),weights=weights,control=list(trace=TRUE,epsilon=1e-14)) 
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

BinomialY <- function(Y,X,weights,m,method="Nelder-Mead",start=FALSE,initial.values=NULL,maxit=1000,reltol=1e-15,trace=1){ 
  
  # Y: a numerical vector for the response variable
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the weights to use for elements of X.
 
  p <- ncol(X)
  
  if(start==FALSE){
    
    modelY   <- glm(cbind(Y,m-Y) ~ X,family=binomial(),weights=weights,control=list(trace=TRUE,epsilon=1e-14)) 
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
    
    modelY    <- glm(Y ~ X,family=Gamma(link="log"),weights=weights,control=list(trace=TRUE,epsilon=1e-14)) 
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



glgcwm <- function(
	Y,			                 # numerical vector for the response variable
  X,			                 # (nxp) matrix for the covariates
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
  
X <- as.matrix(X)
n <- length(Y)  # sample size
p <- ncol(X)    # number of variables for X
  
  if(is.null(X))     stop('Hey, we need some data, please! X is null')
  if(!is.matrix(X))  stop('X needs to be in matrix form')
  if(!is.numeric(X)) stop('X is required to be numeric')
  if(n == 1)   stop('nrow(X) is equal to 1')
  if(any(is.na(X)))  stop('No NAs allowed.')
  if(is.null(k)) stop('k is NULL')
  k <- as.integer(ceiling(k))
  if(!is.integer(k)) stop('k is not a integer')
  if(any(k < 1)) stop('k is not a positive integer')
  
prior   <- numeric(k) # weights

# X

muX     <- array(0,c(p,k),dimnames=list(paste("X.",1:p,sep=""),paste("comp.",1:k,sep="")))
VarX    <- array(0,c(p,p,k),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("comp.",1:k,sep="")))
invVarX <- array(0,c(p,p,k),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("comp.",1:k,sep="")))
PX      <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))

# Y|x

beta0   <- numeric(k)
beta1   <- array(0,c(p,k),dimnames=list(paste("X.",1:p,sep=""),paste("comp.",1:k,sep="")))
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

  if (!is.null(seed)) set.seed(seed)
	z  <- array(runif(n*k),c(n,k)) # soft posterior probabilities (no-normalized) (n x k) 
	z  <- z/rowSums(z)             # soft posterior probabilities (n x k)
  } 

if(initialization=="random.hard"){

  if (!is.null(seed)) set.seed(seed)
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

  # - #
	# X #
	# - #
  
	for(j in 1:k){ 
	  temp          <- mstep(x=X,wt=z[,j])
	  muX[,j]       <- temp$center
	  VarX[,,j]     <- temp$cov
	  invVarX[,,j]  <- solve(temp$cov)
	  PX[,j]        <- (2*pi)^(-p/2)*( ifelse(p>1,det(VarX[,,j]),VarX[,,j]) )^(-1/2)*exp(-1/2*mahalanobis(x=X, center=muX[,j], cov=invVarX[,,j], inverted=TRUE))                                   
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
	    modelY      <- PoissonY(Y,X,weights=z[,h]) #,mustart=muY[,h]
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
	    modelY      <- BinomialY(Y,X,weights=z[,h],m=mY) #,mustart=muY[,h]
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
  
	llvalues          <- log(rowSums(matrix(rep(prior,n),n,k,byrow=TRUE)*PY*PX)) 
	loglik[iteration] <- sum(llvalues)
  
  # ----------------------------------------------- #
	# Aitkane's Acceleration-Based Stopping Criterion #
	# ----------------------------------------------- #
	
	if(iteration>2 & k > 1){
	  a[iteration-1]      <- (loglik[iteration]-loglik[iteration-1])/(loglik[iteration-1]-loglik[iteration-2])
	  aloglik[iteration]  <- loglik[iteration-1]+(1/(1-a[iteration-1])*(loglik[iteration]-loglik[iteration-1]))
	  if(abs(aloglik[iteration]-loglik[iteration])<threshold) 
	    check <- 1 
	}
	
	if(iteration==iter.max | k==1) check <- 1
  
	cat("*")
	iteration <- iteration + 1
  
	# ++++++ #
	# E-Step #
	# ++++++ #
	
	z.num  <- matrix(rep(prior,n),n,k,byrow=TRUE)*PY*PX  # (n x k)
	z.den  <- rowSums(z.num)                             # n-vector
	z      <- z.num/matrix(rep(z.den,k),ncol=k)          # (n x k)
  
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

###########################
## Information criteria  ##
###########################
    
if(familyY=="Gaussian" | familyY=="Gamma")  
  npar <- (k-1) + (p+2)*k + p*k + mixture:::ncovpar(modelname="VVV", p=p, G=k)  
if(familyY=="Binomial" | familyY=="Poisson")  
  npar <- (k-1) + (p+1)*k + p*k + mixture:::ncovpar(modelname="VVV", p=p, G=k)  
    
AIC       <- 2*finalloglik - npar*2
BIC       <- 2*finalloglik - npar*log(n)

z.const    <- (z<10^(-322))*10^(-322)+(z>10^(-322))*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
hard.z     <- (matrix(rep(apply(z,1,max),k),n,k,byrow=F)==z)*1
ECM        <- sum(hard.z*log(z.const))
ICL        <- BIC+ECM

beta <- rbind(beta0,beta1)
dimnames(beta)[[1]] <- paste("beta",0:p,sep="")
dimnames(beta)[[2]] <- paste("comp",1:k,sep="")

result <- list(
  Y         = Y,            # Y: a numerical vector for the response variable
  X         = X,            # X: a numerical vector for the unique covariate
  familyY   = familyY,      # the exponential distribution used for Y|x 
  p         = p,            # number of covariates
  k         = k,            # number of groups
  n         = n,            # sample size
  npar      = npar,
  mY        = mY,           # mY=1 for Bernoulli; mY>1 for Binomial
  prior     = prior,
  muX       = muX,
  VarX      = VarX,
  PX        = PX,
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
  BIC       = BIC,
  ICL       = ICL,      # alla McNicholas
  call      = match.call()
)

class(result) <- "glgcwm"
return(result)

alarm()

}

dgam <- function(x,mu,nu,log=FALSE){
  
  dgamma(x, shape=nu, scale = mu/nu, log = log)
  
}

rgam <- function(x,mu,nu){
  
  rgamma(x, shape=nu, scale = mu/nu)
  
}



#####################
## model Selection ##
#####################

gcwm <- function(
  Y,  		                # numerical vector for the response variable
  X,                       # (nxp) matrix for the covariates
  familyY = "Gaussian",    # the exponential distribution used for Y|x 
  k=2,                  # minimum number of groups
  ic=c("BIC", "AIC", "ICL"),
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
  ic <- match.arg(ic)
  call=match.call()
  
  n <- length(X)
  gridk    <- k
  numk     <- length(gridk)
  
  IC <- array(0,c(numk,4),dimnames=list(paste(gridk,"groups",sep=" "),c("loglik","AIC","BIC","ICL")))

  bestic <- -1e100
  for(i in 1:numk){
    temp <- glgcwm(
      Y=Y,  		                 
      X=X,			                 
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
    if (eval(parse(text=paste0("temp$",ic))) > bestic) {
      besttemp <- temp
      bestic <- eval(parse(text=paste0("temp$",ic)))
    }
    IC[i,1] <- temp$loglik
    IC[i,2] <- temp$AIC
    IC[i,3] <- temp$BIC
    IC[i,4] <- temp$ICL
    
  }

  cat("\n\n")
  cat("# ----------------------- #","\n")
  cat("# Model Selection Results #","\n")
  cat("# ----------------------- #","\n\n")
  
  cat("Best ",ic," value of",eval(parse(text=paste0("besttemp$",ic))),"obtained for k =",besttemp$k,"group(s)","\n\n")
  if(loglikplot==TRUE & besttemp$k>1){
    temp <- glgcwm(
      Y=Y,    	                 
      X=X,			                 
      familyY = familyY,     
      k=besttemp$k,                            
      mY=mY,                         
      method=method,                
      initialization=initialization,  
      start.z=start.z,                 		
      iter.max=iter.max,                
      threshold=threshold,            
      loglikplot=loglikplot,
      seed=seed
    ) 
  }
  besttemp$call <- call
  return(
    c(besttemp,list(IC=IC))
  )  
  
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



