logit<-function(p){
if ((p==0) | (p==1)) stop("probability equals one or zero in logit function")
g<-log(p/(1-p))
return(g)
}
