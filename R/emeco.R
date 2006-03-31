
###
### main function
###
ecoML <- function(formula, data = parent.frame(), N=NULL, supplement = NULL, 
                  theta.start = c(0,0,1,1,0), fix.rho = TRUE,
                  context = FALSE, sem = TRUE, epsilon=10^(-10),
                  maxit = 1000, loglik = TRUE, verbose= TRUE) { 

  
  ## getting X and Y
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- model.response(model.frame(tt, data=data))
  
  #n.var: total number of parameters involved in the estimation
  #n.par: number of nonstatic paramters need to estimate through EM 
  #       also need SEM
  #ndim: dimension of the multivariate normal distribution

  ndim<-2
  if (context) ndim<-3

  n.var<-2*ndim+ ndim*(ndim-1)/2
 
  n.par<-n.S<-n.var 
  if (context) {
      n.par<-n.S<-n.var-2
   }

  if (fix.rho) n.par<-n.par-1

  r12<-NULL
  if (fix.rho) 
     r12<-theta.start[n.par+1]

  flag<-as.integer(context)+2*as.integer(fix.rho)+2^2*as.integer(sem)

  ##checking data
  tmp <- checkdata(X, Y, supplement, ndim)
  bdd <- ecoBD(formula=formula, data=data)
  n <- tmp$n.samp+tmp$survey.samp+tmp$samp.X1+tmp$samp.X0
  inSample.length <- ndim*tmp$n.samp


  #if NCAR and the user did not provide a theta.start
  if (context && (length(theta.start)==5) ) 
    theta.start<-c(0,0,1,1,0,0,0)

  ## Fitting the model via EM  
  res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
            as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
            as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
            as.double(tmp$survey.data),
            as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
            as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
            as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
            as.integer(flag),as.integer(verbose),as.integer(loglik),
            optTheta=rep(-1.1,n.var), pdTheta=double(n.var),
            S=double(n.S+1),inSample=double(inSample.length),DMmatrix=double(n.par*n.par),
            itersUsed=as.integer(0),history=double((maxit+1)*(n.var+1)),
            PACKAGE="eco")

  ##record results from EM
  theta.em<-res$pdTheta
  theta.fisher<-param.trans(theta.em, transformation="Fisher")
  iters.em<-res$itersUsed
  mu.em <- matrix(rep(NA,iters.em*ndim),ncol=ndim)
  sigma.log.em <- matrix(rep(NA,iters.em*ndim),ncol=ndim)
  loglike.em <- as.double(rep(NA,iters.em))
  nrho<-length(theta.em)-2*ndim
  rho.fisher.em <- matrix(rep(NA,iters.em*nrho),ncol=nrho)
  for(i in 1:iters.em) {
    mu.em[i,1:ndim]=res$history[(i-1)*(n.var+1)+(1:ndim)]
    sigma.log.em[i,1:ndim]=res$history[(i-1)*(n.var+1)+ndim+(1:ndim)]
     if (nrho!=0)
    rho.fisher.em[i, 1:nrho]=res$history[(i-1)*(n.var+1)+2*ndim+(1:nrho)]
    loglike.em[i]=res$history[(i-1)*(n.var+1)+2*ndim+nrho+1]
  }

  ## In sample prediction of W
  W <- matrix(rep(NA,inSample.length),ncol=ndim)
  for (i in 1:tmp$n.samp)
    for (j in 1:ndim)
      W[i,j]=res$inSample[(i-1)*2+j]

  ## SEM step
  iters.sem<-0

  if (sem) 
  {
    DM <- matrix(rep(NA,n.par*n.par),ncol=n.par)

    res <- .C("cEMeco", as.double(tmp$d), as.double(theta.start),
              as.integer(tmp$n.samp),  as.integer(maxit), as.double(epsilon),
              as.integer(tmp$survey.yes), as.integer(tmp$survey.samp), 
              as.double(tmp$survey.data),
              as.integer(tmp$X1type), as.integer(tmp$samp.X1), as.double(tmp$X1.W1),
              as.integer(tmp$X0type), as.integer(tmp$samp.X0), as.double(tmp$X0.W2),
              as.double(bdd$Wmin[,1,1]), as.double(bdd$Wmax[,1,1]),
              as.integer(flag),as.integer(verbose),as.integer(loglik),
              res$pdTheta, pdTheta=double(n.var), S=double(n.S+1),
              inSample=double(inSample.length),DMmatrix=double(n.par*n.par),
              itersUsed=as.integer(0),history=double((maxit+1)*(n.var+1)),
              PACKAGE="eco")     
  
    iters.sem<-res$itersUsed
    for(i in 1:n.par)
      for(j in 1:n.par)
        DM[i,j]=res$DMmatrix[(i-1)*n.par+j]

   print(n.par)
   print(DM)

} 

    suff.stat<-res$S
  if (context && (!fix.rho))
      {
	 suff.stat<-rep(0,(n.var+1))
         suff.stat[1]<-mean(logit(c(X,supplement[,3])))
         suff.stat[2:3]<-res$S[1:2]
         suff.stat[4]<-mean((logit(c(X, supplement[,3])))^2)
         suff.stat[5:6]<-res$S[3:4]
         suff.stat[7:8]<-res$S[6:7]
         suff.stat[9]<-res$S[5]
         suff.stat[10]<-res$S[8]
      }
   




  ## output
  res.out<-list(call = mf, Y = Y, X = X, N = N, 
                fix.rho = fix.rho, context = context, sem=sem, epsilon=epsilon,
		theta.em=theta.em, r12=r12, 
               sigma.log = theta.fisher[(ndim+1):(2*ndim)], suff = res$S[1:n.par],
                loglik = res$S[n.par+1], iters.em = iters.em, 
                iters.sem = iters.sem, mu.em = mu.em,
                sigma.log.em = sigma.log.em,
                rho.fisher.em = rho.fisher.em, loglike.em = loglike.em,
                W = W)
  
  if (sem) {
    res.out$DM<-DM

  }

  class(res.out) <- "ecoML"
  return(res.out)
}
