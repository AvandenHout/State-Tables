# Function to summarise a simulation study

# Ardo, UCL 2015

# Function to generate results:
SimSummary <- function(p.true, results, results.SE, digits=3, PLOT=FALSE,
                       col.box="grey",main=NULL, TRACE=FALSE, 
                       PLOT.D = FALSE, groups = NULL, DENSE=TRUE){
  
  # Checks:
  if(length(p.true)!=ncol(results)){
    stop("\n\nERROR. Number of true parameter values is not consistent
    with matrix with results.\n\n")
  }
  if(PLOT.D){ PLOT <- FALSE; TRACE <- FALSE}
  
  # Deal with NAs:
  if(any(is.na(results.SE))){
    S <- nrow(results)
    for(j in 1:ncol(results.SE)){
        if(any(is.na(results.SE[,j]))){
         ind <- which(is.na(results.SE[,j]))
         results.SE <- results.SE[-ind,]
         results    <- results[-ind,]
         }
     }
     cat("\n\nWARNING. ",S-nrow(results)," replication(s) are ignored because of
          NAs in the matrix with the SEs.\n\n")
  }
   
     
  # Prelim:
  L <- length(p.true)
  S <- nrow(results)
  cat("Number of simulation S =",S,".\n")
  cat("Summary simulations results:\n")
  
  # Bias etc:
  estim <- rep(NA,L); for(i in 1:L){estim[i]<-mean(results[,i])}
  bias  <- estim-p.true
  rel.bias  <- abs(bias/p.true) 
  RMSE <- rep(NaN,L)
  for(i in 1:L){RMSE[i] <- sqrt(mean((results[,i]-p.true[i])^2))}
  
  # For median:
  estim.mdn <- rep(NA,L)
  for(i in 1:L){estim.mdn[i] <- median(results[,i])}
  bias.mdn  <- estim.mdn-p.true
  rel.bias.mdn  <- abs(bias.mdn/p.true) 
  
  # Emperical SD and average SE:
  empSD <- rep(NaN,L)
  for(i in 1:L){ empSD[i] <- sqrt(var(results[,i])) }
  averageSE <- rep(NaN,L)
  for(i in 1:L){ averageSE[i] <- sqrt(mean(results.SE[,i]^2))}
      
  # ACPs:
  if(is.null(results.SE)){
   ACP <- rep(NA,L)
  }else{
   ACP <- rep(0,L)
   for(s in 1:S){
    for(i in 1:L){
     lb <- results[s,i]-1.96*results.SE[s,i]
     ub <- results[s,i]+1.96*results.SE[s,i]
     if(lb< p.true[i] & p.true[i]<ub){ACP[i] <- ACP[i]+1/S}
    }
   }
  }

  # Output to screen:
  if(!DENSE){
  print(round(cbind(true=p.true,mean=estim, bias=bias,rel.bias=rel.bias,
          RMSE=RMSE,empSD=empSD,aveSE=averageSE,
          ACP=ACP,bias.mdn=estim.mdn-p.true,
          rel.bias.mdn=rel.bias.mdn),digits))
  }else{
  print(round(cbind(true=p.true,mean=estim, bias=bias,rel.bias=rel.bias,
                      RMSE=RMSE),digits))
    
  }
 
 # Name labels for plots:
  names <- c(expression(theta[1]),expression(theta[2]),expression(theta[3]),
         expression(theta[4]),expression(theta[5]),expression(theta[6]),
         expression(theta[7]),expression(theta[8]),expression(theta[9]),
         expression(theta[10]),expression(theta[11]),expression(theta[12]))[1:L]
  names <- c("beta1","beta2") 
  
  # Set up plotting framwork:
  if(PLOT & TRACE){opar<-par(mfrow=c(1,2))}
  if(PLOT & !TRACE){
     if(is.null(groups)){ opar <- par(mfrow=c(1,1)) }
     if(!is.null(groups)){ 
       ngroups <- length(groups)
       opar <- par(mfrow=c(1,ngroups))
     }   
  }
  if(!PLOT & TRACE){opar<-par(mfrow=c(1,1))}
  # Boxplots:
  if(PLOT){
    boxplot.dta <- as.data.frame(results-matrix(p.true,S,L,byrow=TRUE))
    if(is.null(main)){
      main <- paste(paste("Estimates minus true value (",as.character(S))," replications)")
    }
    if(is.null(groups)){
      boxplot(boxplot.dta,col=col.box,names=names,main=main)
      abline(h=0,lwd=2,col=1,lty=2)
    }
    if(!is.null(groups)){
      for( i in 1:ngroups){
        index <- groups[[i]] 
        boxplot(boxplot.dta[,index],col=col.box,names=names[index],main=main)
        abline(h=0,lwd=2,col=1,lty=2)
      }
    }
  }
  
  # Trace plots:
  if(TRACE){
   eps <- max(c( (results[1,1]-p.true[1])/p.true[1], (results[1,L]-p.true[L])/p.true[L] ))
   plot(c(1,S),c(1-eps,L+eps),type="n",main="Trace plot relative bias",
     ylab="",xlab="Replication index",axes=FALSE,cex.lab=1.2)
   axis(1,at=c(1,S), label=c(1,S),tick=TRUE)
   axis(2,at=1:L, label=names,tick=FALSE, las = 1)
   for(i in 1:L){
     # Compute running bias:
     running.rel.bias <- rep(NA,S)
     for(j in 1:S){running.rel.bias[j] <- (mean(results[1:j,i])-p.true[i])/p.true[i]}
     # Plot running rel.bias:
     lines(1:S,i+running.rel.bias,col=1,lwd=3)
     abline(h=i,col="grey",lty=2,lwd=2)
   }
   }
   
   # Density plots:
   if(PLOT.D){
    opar <- par(mfrow=c(3,floor(L/3)+1))
    for(i in 1:L){
     plot(density(results[,i]),lwd=2,col=1,main="",xlab=names[i])
     abline(v= estim[i],col="black",lwd=4)
     abline(v= estim.mdn[i],col="green",lwd=2)
     abline(v= p.true[i],col="grey",lwd=3)     
     }
   }
  
   # Tidy up:
   opar <- par(mfrow=c(1,1))  
}  


