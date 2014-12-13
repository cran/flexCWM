ICget <- function(object,criteria){
  if(missing(criteria)) criteria <- .ICnames()
  else criteria <- match.arg(criteria, .ICnames(),several.ok =TRUE) 
  lc <- length(criteria)
  lm <- length(object$models)
  res <- matrix(NA,lm,lc,dimnames=list(1:length(object$models),criteria))
  k <- m <- vector(length =lm)
  for(i in 1:lm){
    res[i,] <- sapply(1:lc, FUN=function(j) object$models[[i]]$IC[[criteria[j]]])
    k[i] <- object$models[[i]]$k
    if(!is.null(object$models[[i]]$concomitant$normal.model)) m[i] <-  object$models[[i]]$concomitant$normal.model
  }
  attr(res,"k") <- k
  if(!is.null(object$models[[i]]$concomitant$normal.model))attr(res,"normal.model") <- m
  class(res) <- "cwm.IC"
  res
}
print.cwm.IC <- function(x, digits = max(3L, getOption("digits") - 2L),...){
  class(x) <-"matrix"
  d <- data.frame(x)
  if (!is.null(attributes(x)$normal.model)){
    d <- data.frame(model= attributes(x)$normal.model,d)
  }
  
  d <- data.frame(k = attributes(x)$k,d)
  
  print(d,digits=digits) 
}
bestmodel <- function(object, criteria, k=NULL, modelXnorm=NULL){
  if(missing(criteria)) criteria <- .ICnames()
  else criteria <- match.arg(criteria, .ICnames(),several.ok =TRUE) 
  a   <- ICget(object,criteria)
  w <- TRUE
  if(!is.null(k))  w <- w & attr(a,"k") %in% k
  if(!is.null(modelXnorm))  w <- w & attr(a,"normal.model") %in% modelXnorm
  a <-a[w,,drop=FALSE]
  best <- strtoi(rownames(a)[sapply(1:length(criteria), function (i) which(t(a[,i])==max(a[,i])))])
  names(best) <- criteria
  best
}
modelget <- function(object,criteria="BIC",k=NULL,modelXnorm=NULL){
  criteria <- match.arg(criteria,.ICnames())
  n <- bestmodel(object,criteria,k,modelXnorm)
  if (is.null(criteria)) stop("No model found with specified args.")
  if (length(n)>1) stop("More than one model matches the conditions specified.")
  foo <- object$models[[n]]
  object$models <- NULL
  object$models[[1]] <- foo
  invisible(object)
}
.ModelNames <- function (model){
  type <- switch(EXPR = as.character(model), E = "univariate, equal variance", 
                 V = "univariate, unequal variance", EII = "spherical, equal volume", 
                 VII = "spherical, varying volume", EEI = "diagonal, equal volume and shape", 
                 VEI = "diagonal, equal shape", EVI = "diagonal, equal volume, varying shape", 
                 VVI = "diagonal, varying volume and shape", EEE = "elliposidal, equal volume, shape and orientation", 
                 VEE = "elliposidal, equal shape and orientation", EVE = "elliposidal, equal volume and orientation", 
                 VVE = "ellipsoidal, equal orientation", EEV = "ellipsoidal, equal volume and shape", 
                 VEV = "ellipsoidal, equal shape", EVV = "elliposidal, equal volume", 
                 VVV = "ellipsoidal, varying volume, shape, and orientation", 
                 X = "univariate normal", XII = "spherical multivariate normal", 
                 XXI = "diagonal multivariate normal", XXX = "elliposidal multivariate normal", 
                 warning("invalid model"))
  return(list(model = model, type = type))
}
.ICnames <- function(x=NULL){
  v <- c("AIC", "AICc", "AICu", "AIC3", "AWE","BIC", "CAIC", "ICL")
  if (is.null(x)) return(v)
  else return (v[x])
}