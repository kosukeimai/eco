".First.lib" <- function(lib, pkg) 
  library.dynam("eco", pkg, lib)


eco <- function(n.gen, prior=TRUE, link=1, Y, X, nu0=4, tau0=1,
mu0=c(0,0), S0=diag(10,2), a.update=FALSE, burn.in=0, nth.draw=5,
predict=TRUE, a0=0.1, b0=0.1){ 	
                
if (a.update & !prior) stop("a.update=TRUE only when prior is TRUE!")
if (burn.in>=n.gen) stop("Gibbs draws less than number of burn-ins!") 
  n.samp <- length(Y)	 
  d <- cbind(X, Y)
  if (prior){	
    n.a<- floor((n.gen-burn.in)/nth.draw)
    n.par <- n.a * n.samp 
    n.w <- n.a * n.samp
    unit.a<-1
    unit.par<-n.samp
    unit.w<-n.samp	
   res <- .C("cDPeco", as.integer(n.gen), as.integer(link),  as.double(d), as.integer(n.samp), as.integer(nu0), as.double(tau0), as.double(mu0), as.double(S0), as.integer(a.update),
            as.integer(burn.in), as.integer(nth.draw), as.integer(predict),
            as.double(a0), as.double(b0),  pdSMu0=double(n.par),
            pdSMu1=double(n.par), pdSSig00=double(n.par),
            pdSSig01=double(n.par), pdSSig11=double(n.par),
            pdSW1=double(n.w), pdSW2=double(n.w), 
	    pdSWt1=double(n.w), pdSWt2=double(n.w), 
            pdSa=double(n.a), pdSn=integer(n.a))

  mu1.post<-matrix(res$pdSMu0, n.a, unit.par, byrow=T)
  mu2.post<-matrix(res$pdSMu1, n.a, unit.par, byrow=T)
  Sigma11.post<-matrix(res$pdSSig00, n.a, unit.par, byrow=T)
  Sigma12.post<-matrix(res$pdSSig01, n.a, unit.par, byrow=T)
  Sigma22.post<-matrix(res$pdSSig11, n.a, unit.par, byrow=T)
  w1.post<-matrix(res$pdSW1, n.a, unit.w, byrow=T)
  w2.post<-matrix(res$pdSW2, n.a, unit.w, byrow=T)
  w1.pred<-matrix(res$pdSWt1, n.a, unit.w, byrow=T)
  w2.pred<-matrix(res$pdSWt2, n.a, unit.w, byrow=T)
  a.post<-matrix(res$pdSa, n.a, unit.a, byrow=T)
  nstar<-matrix(res$pdSn, n.a, unit.a, byrow=T)
res.out <- list(model="Dirichlet Process Prior", ALPHA=a.update, burn.in=burn.in, X=X, Y=Y, nu0=nu0, tau0=tau0, mu0=mu0, S0=S0, a0=a0, b0=b0, mu1.post=mu1.post, mu2.post=mu2.post, Sigma11.post=Sigma11.post, Sigma12.post=Sigma12.post, Sigma22.post=Sigma22.post, w1.post=w1.post, w2.post=w2.post, w1.pred=w1.pred, w2.pred=w2.pred, a.post=a.post, nstar=nstar)  
  }	 
  else{
    n.a<- floor((n.gen-burn.in)/nth.draw)
    n.par <- n.a
    n.w <- n.a * n.samp
    unit.a<-1
    unit.par<-1
    unit.w<-n.samp	
    res <- .C("cBaseeco", as.integer(n.gen), as.integer(link),  as.double(d), as.integer(n.samp), as.integer(nu0), as.double(tau0), as.double(mu0), as.double(S0),as.integer(burn.in), as.integer(nth.draw), as.integer(predict),
            pdSMu0=double(n.par),
            pdSMu1=double(n.par), pdSSig00=double(n.par),
            pdSSig01=double(n.par), pdSSig11=double(n.par),
            pdSW1=double(n.w), pdSW2=double(n.w), 
	    pdSWt1=double(n.w), pdSWt2=double(n.w))
  mu1.post<-matrix(res$pdSMu0, n.a, unit.par, byrow=T)
  mu2.post<-matrix(res$pdSMu1, n.a, unit.par, byrow=T)
  Sigma11.post<-matrix(res$pdSSig00, n.a, unit.par, byrow=T)
  Sigma12.post<-matrix(res$pdSSig01, n.a, unit.par, byrow=T)
  Sigma22.post<-matrix(res$pdSSig11, n.a, unit.par, byrow=T)
  w1.post<-matrix(res$pdSW1, n.a, unit.w, byrow=T)
  w2.post<-matrix(res$pdSW2, n.a, unit.w, byrow=T)
  w1.pred<-matrix(res$pdSWt1, n.a, unit.w, byrow=T)
  w2.pred<-matrix(res$pdSWt2, n.a, unit.w, byrow=T)
res.out <- list(model="Normal prior", burn.in=burn.in, X=X, Y=Y, nu0=nu0, tau0=tau0, mu0=mu0, S0=S0, mu1.post=mu1.post, mu2.post=mu2.post, Sigma11.post=Sigma11.post, Sigma12.post=Sigma12.post, Sigma22.post=Sigma22.post, w1.post=w1.post, w2.post=w2.post, w1.pred=w1.pred, w2.pred=w2.pred)  


}

  class(res.out) <- "eco"
  return(res.out)

}


