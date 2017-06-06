LRtest <- function(...,main=NULL,latex=FALSE){
  suc <- is.null(main)  #Succesive testing?
  args <- list(...)
  if(is.null(names(args)))names(args) <- sapply(substitute(f(...)),deparse)[-1]
  nparms <- sapply(args,function(x)length(x$fit$par) + ncol(x$data$X))
  lik <- sapply(args,function(x)sum(x$fit$objective))  
  casenames <- sub("env","Model",names(args))
  if(suc)minus2logQ <- c(NA,2*diff(lik)) else minus2logQ <- 2*(lik-main$fit$objective)
  if(suc)df <- c(NA, -diff(nparms)) else df <- length(main$par)-nparms
  p <- 1-pchisq(minus2logQ,df=df)
  m <- cbind(lik,minus2logQ,nparms,df,p)
  rownames(m) <- casenames
  if(latex)colnames(m)[1:2] <- c('$-\\log L$','$-2\\log Q$')
  else colnames(m)[1:2] <- c('-log L','-2log Q')
  m
}
