


eco.em <- function(Y, X, data = parent.frame(),supplement=NULL, 
      theta.old=c(0,0,1,0,1), 
      convergence=0.0001, iteration.max=20,
      n.draws = 10, draw.max=200, printon=TRUE) {

  ## checking inputs
  if ((dim(supplement)[2] != 2) && (length(supplement)>0))
    stop("Error: use n by 2 matrix for survey data")
  
  call <- match.call()

  ff <- as.formula(paste(call$Y, "~ -1 +", call$X))
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(ff, data)
  Y <- model.response(model.frame(ff, data=data))

  ##survey data
  if (length(supplement) == 0) {
    survey.samp <- 0
    survey.data <- 0
    survey.yes<-0
  }
  else {
    survey.samp <- length(supplement[,1])
    survey.data <- as.matrix(supplement)
    survey.yes<-1
  }
  
  ind<-c(1:length(X))
  X1type<-0
  X0type<-0
  samp.X1<-0
  samp.X0<-0
  X1.W1<-0
  X0.W2<-0
  
  ##Xtype x=1
  X1.ind<-ind[along=(X==1)]
  if (length(X[X!=1])<length(X)){
    X1type<-1
    samp.X1<-length(X1.ind)
    X1.W1<-Y[X1.ind]
  }
  
  ##Xtype x=0
  X0.ind<-ind[along=(X==0)]
  if (length(X[X!=0])<length(X)){
    X0type<-1
    samp.X0<-length(X0.ind)
    X0.W2<-Y[X0.ind]
  }
  
  XX.ind<-setdiff(ind, union(X0.ind, X1.ind))
  X.use<-X[XX.ind]
  Y.use<-Y[XX.ind]

  order.old<-order(c(XX.ind, X0.ind, X1.ind))
  
  ## fitting the model
  n.samp <- length(Y.use)	 
  d <- cbind(X.use, Y.use)
  n.var<-5

  cdiff<-1
  i<-1

while ((cdiff > convergence) && (i<iteration.max))
{
  res <- .C("cEMeco", as.double(d), as.double(theta.old),
	      as.integer(n.samp),  as.integer(n.draws), 
              as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
   	      as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
   	      as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
	      pdTheta=double(n.var),
              PACKAGE="eco")

  temp<-res$pdTheta
  cdiff<-as.numeric(t(temp-theta.old)%*%(temp-theta.old))
  theta.old<-temp
  if (printon) {
  cat("\ni=  ", i, "  cdiff=", cdiff, "\n")
  cat(theta.old)
  }
  i<-i+1
  if (draw<=draw.max) draw<-draw+50
}

#  theta<-list(mu=c(theta.old[1], theta.old[2]), Sigma=matrix(theta.old[c(3,4,4,5)], 2,2))



cat("\n")
cat("Estimates based on EM", "\n")
cat("Mean:", "\n")
print(theta.old[1:2])
cat("\n")
cat("Covarianace Matrix:", "\n")
print(matrix(theta.old[c(3,4,4,5)],2,2))

res.out<-list(theta=theta.old)
  class(res.out) <- "eco"
  return(res.out)

}


eco.sem<-function(Y, X, data = parent.frame(),supplement=NULL, 
      theta.old=c(0,0,1,0,1), theta.em=NULL,
      R.convergence=0.001, iteration.max=50,
      n.draws = 10, draw.max=200, printon=TRUE) {

  ## checking inputs
  if ((dim(supplement)[2] != 2) && (length(supplement)>0))
    stop("Error: use n by 2 matrix for survey data")
  
  call <- match.call()

  ff <- as.formula(paste(call$Y, "~ -1 +", call$X))
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(ff, data)
  Y <- model.response(model.frame(ff, data=data))

  ##survey data
  if (length(supplement) == 0) {
    survey.samp <- 0
    survey.data <- 0
    survey.yes<-0
  }
  else {
    survey.samp <- length(supplement[,1])
    survey.data <- as.matrix(supplement)
    survey.yes<-1
  }
  
  ind<-c(1:length(X))
  X1type<-0
  X0type<-0
  samp.X1<-0
  samp.X0<-0
  X1.W1<-0
  X0.W2<-0
  
  ##Xtype x=1
  X1.ind<-ind[along=(X==1)]
  if (length(X[X!=1])<length(X)){
    X1type<-1
    samp.X1<-length(X1.ind)
    X1.W1<-Y[X1.ind]
  }
  
  ##Xtype x=0
  X0.ind<-ind[along=(X==0)]
  if (length(X[X!=0])<length(X)){
    X0type<-1
    samp.X0<-length(X0.ind)
    X0.W2<-Y[X0.ind]
  }
  
  XX.ind<-setdiff(ind, union(X0.ind, X1.ind))
  X.use<-X[XX.ind]
  Y.use<-Y[XX.ind]

  order.old<-order(c(XX.ind, X0.ind, X1.ind))
  
  ## fitting the model
  n.samp <- length(Y.use)	 
  d <- cbind(X.use, Y.use)
  n.var<-5

  R.t1<-matrix(0, n.var, n.var)
  R.t2<-matrix(0, n.var, n.var)
  rowdiff<-rep(1, n.var)

  k<-1

  Rconverge<-FALSE

while (!Rconverge && (k<iteration.max))
{
  res <- .C("cEMeco", as.double(d), as.double(theta.old),
	      as.integer(n.samp),  as.integer(n.draws), 
              as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
   	      as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
   	      as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
	      pdTheta=double(n.var),
              PACKAGE="eco")

  theta.t<-res$pdTheta
  #cat("\n theta.t \n")
  #print(theta.t)
  #theta.t.i<-rep(NA, n.var)

  Rconverge<-TRUE

  for (i in 1:n.var)
  {
    rowdiff.temp<-0
    if (rowdiff[i]>R.convergence) {
    Rconverge<-FALSE
    theta.t.i<-theta.em
    theta.t.i[i]<-theta.t[i]
   #cat("\n theta.t.i \n")
   #print(theta.t.i)

    temp<-.C("cEMeco", as.double(d), as.double(theta.t.i),
	      as.integer(n.samp),  as.integer(n.draws), 
              as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
   	      as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
   	      as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
	      pdTheta=double(n.var),
              PACKAGE="eco")$pdTheta

  for (j in 1:n.var)
  {
    R.t2[i,j]<-(temp[j]-theta.em[j])/(theta.t[i]-theta.em[i])
    rowdiff.temp<-abs(R.t2[i,j]-R.t1[i,j])+rowdiff.temp    
    }   
   rowdiff[i]<-rowdiff.temp
 }
}

  theta.old<-theta.t
  if (printon) {
  cat("\nk=  ", k, "  rowdiff=:\n")
  print(rowdiff)

  R.t1<-R.t2
  cat("\n R")
  print(R.t2)
  cat("\n")
  print(theta.old)
  }
  k<-k+1
  if (draw<=draw.max) draw<-draw+10
}



cat("\n")
cat("Estimates based on EM", "\n")
cat("Mean:", "\n")
print(theta.em[1:2])
cat("\n")
cat("Covarianace Matrix:", "\n")
print(matrix(theta.em[c(3,4,4,5)],2,2))
cat("DM matrix:", "\n")
print(R.t2)

  res.out<-list(theta=theta.em, DM=R.t2)
  class(res.out) <- "eco"
  return(res.out)

}
  
