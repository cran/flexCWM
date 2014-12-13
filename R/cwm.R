cwm <- function(
  formulaY=NULL,
  familyY= gaussian, # the exponential distribution used for Y|x 
  data=NULL,
  link,
  Xnorm=NULL,
  modelXnorm=NULL,            # models to be fitted in the EM phase of clustering
  Xbin=NULL,
  Xbtrials=NULL,
  Xpois=NULL, #list of strings with names of covariates with Poisson distribuntion 
  Xmult=NULL,
  k=1:3,                  # groups
  initialization=c("random.soft","random.hard","kmeans","mclust","manual"),  # initialization procedure
  start.z=NULL,                 # (n x k)-matrix of soft or hard classification: it is used only if initialization="manual"    
  seed=NULL,
  maxR=1,                  #number of initial random points to be tested
  iter.max=1000,                # maximum number of iterations in the EM-algorithm
  threshold=1.0e-04,             # stopping rule in the Aitken rule
  parallel=FALSE
)              
{
  # Preliminary checks and init ----------------------------------------------
  if(is.null(formulaY) & is.null(Xnorm) & is.null(Xbin) & is.null(Xpois) & is.null(Xmult)) stop("No data were entered.")
  
  initialization<- match.arg(initialization)
  # Get regression terms ----------------------------------------  
  if(!is.null(formulaY)) {
    if (class(familyY)=="character"){
      familyYname <- match.arg(familyY,c("gaussian","poisson","binomial","Gamma","t","inverse.gaussian"))
      if (missing(link)){
        familyY <- switch(familyYname,
                          "Gamma" = Gamma(link="log"),
                          "t" = gaussian(),
                          "inverse.gaussian" = inverse.gaussian(link="log"),
                          do.call(familyYname, list())
        )
      } 
      else{
        if (familyYname=="t") familyY <- do.call("gaussian",list(link=link))
        else familyY <- do.call(familyYname,list(link=link))  
      } 
    } 
    else {
      if (is.function(family)) familyYname <- familyY()$family
      else stop("Invalid value for argument 'familyY'")
    }
    if (is.null(data)){
      data <- model.frame(as.formula(formulaY))
    } else data <- model.frame(as.formula(formulaY),data)
    Y <- as.matrix(model.response(data))
    
    
    if(familyYname=="Binomial" & !is.factor(Y)) mY <- rowSums(Y)
    else mY <- rep(1,nrow(Y))
    
  } 
  else {
    Y <- data <- familyYname <-NULL
  }
  
  n <- unlist(sapply(list(Y,Xnorm,Xpois,Xbin,Xmult), function(x)nrow(cbind(x))))
  if (max(n) != min(n)) stop("Data length mismatch")
  n <- n[1]
  
  
  # Fix variables to use for marginal --------------------------
  colXn <- colXm <- m <-colXb <-colXp<-0
  if(!is.null(Xmult)){
    Xmult <- as.data.frame(Xmult)
    colXm <- ncol(Xmult)
    m <- sapply(1:ncol(Xmult),function(u) length(levels(Xmult[,u])))
    Xmod <- lapply(1:ncol(Xmult),function(u)mclust::unmap(Xmult[,u],levels(Xmult)[u]))
    names(Xmod) <- colnames(Xmult)
  } 
  if(!is.null(Xnorm)){
    Xnorm <- as.matrix(Xnorm)
    colXn <- ncol(Xnorm)
    if(is.null(colnames(Xnorm))) colnames(Xnorm) <- paste0("Xnorm",1:colXn)
  }
  if(!missing(Xpois) && !is.null(Xpois)){
    Xpois <- as.matrix(Xpois)
    colXp <- ncol(Xpois)
    if(is.null(colnames(Xpois))) colnames(Xpois) <- paste0("Xpois",1:colXp)
  } 
  if(!is.null(Xbin)){
    Xbin  <- as.matrix(Xbin)
    colXb <- ncol(Xbin)
    if(is.null(colnames(Xbin))) colnames(Xbin) <- paste0("Xbin",1:colXb)
    if(!is.null(Xbtrials)){
      Xbtrials <- as.vector(Xbtrials)
      if (length(Xbtrials) != colXb) {
        stop(paste0("If provided, Xbtrials has to be a vector of length equal to the number of columns of Xbin."))
      }
    } else Xbtrials <- sapply(1:colXb,function(i) max(Xbin[,i]))
  }
  # More checks ------------------------------------------
  if((initialization=="mclust" | initialization=="kmeans") &
       is.null(Xnorm) & is.null(Xpois) & is.null(Xbin) & is.null(data)){  
    stop(paste0("when initialization is '",initialization, "', numeric variables are needed."))}
  
  lm <- list(k=k)
  if(colXn>0){
    if(colXn==1){
      modelXnormNames <- c("E","V")
      eqmod <- data.frame(name=modelXnormNames,number=c(1,1))
    }
    else{
      modelXnormNames <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","EEV","VVE","VEV","EVV","VVV")
      eqmod <- data.frame(name=modelXnormNames, number=c(1,1,2,2,2,2,3,3,3,3,3,3,3,3))
      
    }  
    if(is.null(modelXnorm)) modelXnorm <- modelXnormNames 
    else modelXnorm <- match.arg(modelXnorm, modelXnormNames, several.ok = TRUE)
    lm$modelXnorm <- modelXnorm
  } else {
    modelXnorm <- eqmod <- NULL
  }
  
  k <- as.integer(ceiling(k))
  par  <- list()
  cc <- 0
  # Compute ----------------------------------------------
  mm <- expand.grid(lm)
  if (!is.null(modelXnorm) & 1 %in% k){
    mm <- merge(mm,eqmod,by.x="modelXnorm",by.y="name")
    mm1 <- mm[mm$k==1,]
    mm2 <- mm1[!duplicated(subset(mm1,select= -modelXnorm, drop=FALSE)),]
  
    if (nrow(mm1) > nrow(mm2)){
      cat("\nWith k=1, some models in modelXnorm are equivalent, so only the first is estimated.\n")
    }
    mm <- rbind(mm2,mm[mm$k!=1,])
  }
  mm <- mm[order(mm$k),,drop=FALSE]
 job <- function(i){
    cat(paste("\nEstimating model ",mm$modelXnorm[i]," with k=",mm$k[i]," ",sep=""))
    cwm2(formulaY=formulaY, familyYname=familyYname,
         data=data, Y=Y,
         Xnorm=Xnorm, Xmult=Xmult, Xbin=Xbin,Xpois=Xpois,
         Xbtrials=Xbtrials,n=n, m=m,colXn=colXn,colXp=colXp,colXb=colXb,colXm=colXm, Xmod=Xmod,
         k=mm$k[i], modelXnorm = mm$modelXnorm[i],           
         familyY = familyY, mY=mY, method="Nelder-Mead", initialization=initialization,  
         start.z=start.z, iter.max=iter.max, threshold=threshold, seed=seed, maxR=maxR)
  }
  if(parallel){
    cores <- getOption("cl.cores", detectCores())
    cat(paste("Using",cores,"cores\n"))
    cl <- makeCluster(cores)
    #clusterExport(cl,envir=environment())
    par <- parLapply(cl=cl,1:nrow(mm),function(i) job(i))
    stopCluster(cl)
  }
 else {
  par <- lapply(1:nrow(mm),function(i) job(i))
 }
  
  
  cat("\n")
  res <- structure(list(
    call=match.call(),
    formulaY=formulaY,
    familyY=familyYname,
    data=data,
    concomitant=list(Xnorm=as.df(Xnorm,"Xnorm"),Xbin=as.df(Xbin,"Xbin"),Xpois=as.df(Xpois,"Xpois"),Xmult=as.df(Xmult,"Xmult")),
    Xbtrials=Xbtrials,
    models=par
  ),class = "cwm")
  print(res)
  invisible(res)
}
as.df <- function(ma,type){
  if (!is.null(ma)){
    df <- data.frame(ma)
    for(i in 1:ncol(ma)) attr(df[[i]],"type") <- type
    df
  } else NULL
}
