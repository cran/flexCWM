plot.cwm <- function(x,criteria="BIC",regr=TRUE, ctype=c("Xnorm","Xbin","Xpois","Xmult"), which="all",...){
  ctype <- match.arg(ctype,several.ok =TRUE)
  criteria <- match.arg(criteria,.ICnames())
  x <- modelget(x,criteria)
  fp <- FALSE
  if (regr &  length(x$models[[1]]$GLModel)>0) {
    .plotregr(x,...)
    fp <- TRUE
  }
  .plotConc(x=x,ctype=ctype,which=which,fp=fp,...)
 # if (i<length(df)) readline("Press <Enter> to continue")
}  
.plotConc <- function(x,ctype,which,fp,main="",col=c(2:(x$models[[1]]$k+1)),lwd=1,lty=2,...){
  
  for (i in ctype){
    df <-x$concomitant[[i]]
    if (!is.null(df))
      for(j in if(which=="all") seq_len(ncol(df)) else which){
        if (fp) readline("Press <Enter> to continue")
        fp <- TRUE
        if (i=="Xnorm") .plotConcNorm(x, as.matrix(df[j]),j,xlab=colnames(df[j]),main=main,col=col,lwd=lwd,lty=lty,...)
        else .plotConcBar(x=x, df=as.matrix(df[j]),ind=j,dtype=i,xlab=colnames(df[j]),type="b",main=main,col=col,lwd=lwd,lty=lty,...) 
      }
    }
 }

.plotConcNorm <- function(x,df,ind,ylab="density",ylim,col,lwd,lty,type,...){
  data <- df
  k <- x$models[[1]]$k
  group <- x$models[[1]]$cluster
  m <- x$models[[1]]$concomitant$normal.mu[ind,]
  s <- x$models[[1]]$concomitant$normal.Sigma[ind,ind,]
  w <- x$models[[1]]$prior
  
  h <- hist(data, breaks=31,plot=FALSE)
  xx <- seq(min(h$breaks),max(h$breaks),l=1001)
  yy <- lapply(1:k, function(i) w[i]*dnorm(xx,m[i],sqrt(s[i])))
  yline <- rowSums(sapply(1:k, function(i) yy[[i]]))
  
  if(missing(ylim)) ylim <- c(0,max(density(data)$y,h$density,yline))*1.15
  
  plot(h,ylim =ylim,freq=F,ylab=ylab,...)
  
  if (k>1) lines(x= xx,y=yline,lwd=lwd+1)
  for(i in 1:k) lines(x= xx,y=yy[[i]],lty=lty,col=col[i],lwd=lwd)
  #d <- lapply(1:k, function(i) curve(w[i]*dnorm(x,m[i],s[i]),add=T,lty=lty,n=1001,col=col[i],lwd=lwd))
    
  for(i in 1:k) rug(data[which(group==i)],col=col[i])
}

.plotConcBar <- function(x,df,dtype,ind,xlab=colnames(df),ylab="probability",xlim=NULL,ylim,col,lwd,lty,type,...){
  k <- x$models[[1]]$k
  group <- x$models[[1]]$cluster
  w <- x$models[[1]]$prior
  if (dtype=="Xbin"){
    n <- x$Xbtrials[ind]
    data <- factor(df, levels=0:n)
    th <- matrix(nrow=k,ncol=n+1)
    for(i in 1:k) th[i,] <- w[i]* dbinom(0:n,n, prob=x$models[[1]]$concomitant$binomial.p[ind,i])
  }
  if (dtype=="Xpois"){
    n <- max(df)
    data <- factor(df, levels=0:n)
    th <- matrix(nrow=k,ncol=n+1)
    for(i in 1:k) th[i,] <- w[i]*dpois(0:n,x$models[[1]]$concomitant$poisson.lambda[ind,i])
  }
  if (dtype=="Xmult"){
    data <- df
    th <- w[i]*t(x$models[[1]]$concomitant$multinomial.prob[[ind]])
  }
  dt <- table(data)
  dt <- dt/sum(dt)
  if(missing(ylim)) ylim <- c(0,max(colSums(th),dt)*1.15)
  mp <- barplot(height=dt,xlim=xlim,ylim =ylim,xlab=xlab,ylab=ylab,...)
  lines(x=mp, y=colSums(th),col=1,type=type,lwd=lwd+1,...)
  if (k>1) for(i in 1:k) lines(x=mp, y=th[i,],col=col[i],type=type,lwd=lwd,...)

}
.plotregr <- function(x, quantiles = c(0.75, 0.95),uncertanty=FALSE, col=c(2:(x$models[[1]]$k+1)),cex = 1, pch=NULL,type,...) {
  if (uncertanty) {
    uncert <- 1 - apply(x$models[[1]]$posterior, 1, max)
    breaks <- quantile(uncert, probs = sort(quantiles))
    I <- rep(3,length(uncert))
    I[uncert < breaks[2]] <- 2
    I[uncert < breaks[1]] <- 1
    cex <- I*cex
  }

  colp  <- col[x$models[[1]]$cluster]
  pchp  <- pch[x$models[[1]]$cluster]
  
    if(ncol(as.matrix(x$data))==2){
      plot(y=x$data[,1],x=x$data[,2],cex=cex, col=colp, pch = pchp,...)
      xnew <- data.frame(seq(min(x$data[,2]), max(x$data[,2]), length.out = 100))
      names(xnew) <- names(x$data)[2]
      for (i in 1:x$models[[1]]$k){
        yhat <- predict(x$models[[1]]$GLModel[[i]]$model,newdata=xnew,type="response")
        lines(y=yhat, x=xnew[[1]], col=col[i],lwd=2)
        
        #yhat <- fitted(x$models[[1]]$GLModel[[i]]$model)[x$models[[1]]$cluster==i]
        #ex   <- x$data[x$models[[1]]$cluster==i,2]
        #lines(y=yhat, x=ex, col=col[i])
      }
    }
    if(ncol(as.matrix(x$data))>2) pairs(as.matrix(x$data),cex=cex, col=colp, pch = pch)

}