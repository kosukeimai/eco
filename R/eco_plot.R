tomoplot<-function(X,Y,truep=NA, truer=NA)
{
	plot(c(0,0), c(0,0), type="l", xlim=c(0,1), ylim=c(0,1), xlab=expression(W[1]), ylab=expression(W[2]))
	for (i in 1:length(X))
	{
	if ((X[i]<1/2) & (Y[i]>X[i]) & (Y[i]<(1-X[i])))
	{
		lines(c(0,1), c(Y[i]/(1-X[i]), (Y[i]-X[i])/(1-X[i])), lwd=0.5)
	}

	if ((Y[i]>1/2) & (X[i]<Y[i]) & ((1-Y[i])<X[i]))
	{
		lines(c((X[i]+Y[i]-1)/X[i], 1), c(1, (Y[i]-X[i])/(1-X[i])), lwd=0.5)
	}

	if ((X[i]>1/2) & (X[i]>Y[i]) & (Y[i]>(1-X[i])))
	{
		lines(c((X[i]+Y[i]-1)/X[i], Y[i]/X[i]), c(1,0), lwd=0.5)
	}

	if ((Y[i]<1/2) & (Y[i]<X[i]) & (X[i]<(1-Y[i])))
	{
		lines(c(0, Y[i]/X[i]), c(Y[i]/(1-X[i]), 0), lwd=0.5)
	}

	if ((!is.na(truep)) & (!is.na(truer)))
	{
		points(truep, truer, pch=20, cex=0.7, col="red")
	}
	}
}

library(MASS)
library(mvtnorm)
library(MCMCpack)
library(KernSmooth)

eco.plot<-function(figname="plot.ps", burn.in=10000, draw=10000, simdata=16, true.den=TRUE, true.W1=NA, true.W2=NA, x=NA, y=NA, z=NA){

plot.star<-FALSE
dpW1.post.vec<-scan("dpw1.txt")
dpW2.post.vec<-scan("dpw2.txt")
dpaW1.post.vec<-scan("dpaw1.txt")
dpaW2.post.vec<-scan("dpaw2.txt")
#dpMu1.post.vec<-scan("dpmu1.txt")
#dpMu2.post.vec<-scan("dpmu2.txt")

nW1.post.vec<-scan("nw1.txt")
nW2.post.vec<-scan("nw2.txt")


#dpWstar1.post.vec<-log(dpW1.post.vec/(1-dpW1.post.vec))
#dpWstar2.post.vec<-log(dpW2.post.vec/(1-dpW2.post.vec))
#nWstar1.post.vec<-log(nW1.post.vec/(1-nW1.post.vec))
#nWstar2.post.vec<-log(nW2.post.vec/(1-nW2.post.vec))

n<-length(dpW1.post.vec)
rsamp<-sample((burn.in+1):n, draw)

postscript(paste(figname))
#par(mfrow=c(2,2))
#par(mfcol=c(2,2), cex.lab=0.8, cex.axis=0.65, oma=c(2,4,4,2), mar=c(1.5,1.5,0.5,0.5), mgp=c(1.35,0.4,0),tcl=-0.2, usr=c(0,1,0,1), pin=c(1.03,1.03))
par(mfrow=c(2,2), cex.lab=0.8, cex.axis=0.65, oma=c(2,4,4,2), mar=c(1.5,1.5,0.5,0.5))
#W plots
   plot(c(0,1), c(0,0), xlim=c(0,1), ylim=c(0,4), xlab=NA, ylab=NA)
   if (true.den==F){
    polygon(bkde(dataXY[,4], range.x=c(0,1), truncate=FALSE), border=NA, col="grey")}
   else {
    polygon(true.W1, border=NA, col="grey")}
    lines(bkde(dpW1.post.vec[rsamp], range.x=c(0,1), truncate=FALSE), type='l', lty=)
    lines(bkde(dpaW1.post.vec[rsamp], range.x=c(0,1), truncate=FALSE), lty=5, lwd=0.7)
    lines(bkde(nW1.post.vec[rsamp], range.x=c(0,1), truncate=FALSE), lty=4, lwd=0.7)
    legend(0.45, 3, c("DP", "DPa", "Normal"), lty=c(1,5,4), cex=0.6, bty="n")  
  title("compare density of W1", cex.main=1)

   plot(c(0,1), c(0,0), xlim=c(0,1), ylim=c(0,4), xlab=NA, ylab=NA)
   if (true.den==F){
    polygon(bkde(dataXY[,5], range.x=c(0,1), truncate=FALSE), border=NA, col="grey")}
   else {
     polygon(true.W2, border=NA, col="grey")}
  lines(bkde(dpW2.post.vec[rsamp], range.x=c(0,1), truncate=FALSE), type='l', lty=1)
    lines(bkde(dpaW2.post.vec[rsamp], range.x=c(0,1), truncate=FALSE), lty=5, lwd=0.7)
    lines(bkde(nW2.post.vec[rsamp], range.x=c(0,1), truncate=FALSE), lty=4, lwd=0.7)

    legend(0.45, 3, c("DP", "DPa", "normal"), lty=c(1,5,4), cex=0.6, bty="n")  
  title("compare density of W2", cex.main=1)



#Wstar plots
if (plot.star==TRUE) {
  plot(bkde(dpWstar1.post.vec[rsamp],  truncate=FALSE), type='l', lty=1, ylim=c(0, 0.5), lwd=0.7, xlab=NA, ylab=NA)
    lines(bkde(nWstar1.post.vec[rsamp], truncate=FALSE), lty=4, lwd=0.7)
    lines(bkde(log(dataXY[,4]/(1-dataXY[,4])), truncate=FALSE), lty=2, lwd=0.7)
    legend(0.45, 3, c("DP", "normal", "truth"), lty=c(1,4,2), cex=0.6, bty="n")  
  title("compare density of Wstar1", cex.main=1)
  plot(bkde(dpWstar2.post.vec[rsamp], truncate=FALSE), type='l', lty=1,  ylim=c(0, 0.5), lwd=0.7, xlab=NA, ylab=NA)
    lines(bkde(nWstar2.post.vec[rsamp], truncate=FALSE), lty=4, lwd=0.7)
    lines(bkde(log(dataXY[,5]/(1-dataXY[,5])), truncate=FALSE), lty=2, lwd=0.7)
    legend(0.45, 3, c("DP", "normal", "truth"), lty=c(1,4,2), cex=0.6, bty="n")  
  title("compare density of Wstar2", cex.main=1)

mtext(paste("simulate data  ", simdata), side=3, outer=T, cex=1.2, line=0.5, at=0.5)
}

par(mfrow=c(1,1))
#contour plots


contour(x,y,z)

points(dpW1.post.vec[rsamp],dpW2.post.vec[rsamp], pch='.', cex=5)
title(main="W contour", sub=paste("simulate data  ", simdata))

contour(x,y,z)

points(dpaW1.post.vec[rsamp],dpaW2.post.vec[rsamp], pch='.', cex=5)
title(main="W contour-a vary", sub=paste("simulate data  ", simdata))

#contour(x,y,z)
#points(dpMu1.post.vec[rsamp],dpMu2.post.vec[rsamp], pch='.', cex=5)
#title(main="Mu contour", sub=paste("simulate data  ", simdata))

dev.off()
}



##bounds conditions

bounds<-function(X,Y) {
    temp<-(X+Y-1)/X
    temp<-as.matrix(temp)
    W1.min<-apply(temp,1,max,0)
    temp<-Y/X
    temp<-as.matrix(temp)
    W1.max<-apply(temp,1,min, 1)
    temp<-(Y-X)/(1-X)
    temp<-as.matrix(temp)
    W2.min<-apply(temp, 1, max, 0)
    temp<-Y/(1-X)
    temp<-as.matrix(temp)
    W2.max<-apply(temp, 1, min, 1)

    bounds<-as.data.frame(cbind(W1.min, W1.max, W2.min, W2.max))
    return(bounds)
}
