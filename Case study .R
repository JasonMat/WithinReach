#########################################################
##                  1. DATA IMPORT                     ##
#########################################################

# Importing data
load("/Users/Jason/Dropbox/2. Science/1. Papers/(2 In progress) Sessility to mobility/Code/Am Nat 3rd revision/HarbourSeal_data.rdata")


# Box from within which to learn about the environment
xst1<-ceiling(0.78*ndx)
Xst1<-floor(0.88*ndx)
yst1<-ceiling(0.11*ndy)
Yst1<-floor(0.95*ndy)

xst2<-ceiling(max(1,min(md3$Txrd)-0*sg))
Xst2<-floor(min(ndx, max(md3$Txrd)+0*sg))
yst2<-ceiling(max(1,min(md3$Tyrd)-0*sg))
Yst2<-floor(min(ndy, max(md3$Tyrd)+0*sg))
outrange<-0 # Decides if the environment will be sampled from the area of the telemetry data (0) or outside (1)

if (outrange==0)
  {
  xst<-xst2
  Xst<-Xst2
  yst<-yst2
  Yst<-Yst2
} else
{
  xst<-xst1
  Xst<-Xst1
  yst<-yst1
  Yst<-Yst1
}


lstxs<-as.matrix(expand.grid(xst:Xst,yst:Yst))
path<-cbind(ndx*(md3$xrd-Xmin)/(Xmax-Xmin),ndy*(md3$yrd-Ymin)/(Ymax-Ymin))
library(RColorBrewer)
cols <- brewer.pal(3, "BrBG")
pal <- colorRampPalette(cols)
par(mfrow=c(1,2))
image(x,y,u1/mask, col=pal(30))
contour(x,y,mask, levels=c(0.5),add=TRUE,drawlabels=FALSE)
lines(path)
polygon(c(xst,xst,Xst,Xst),c(yst,Yst,Yst,yst))
image(x,y,u2/mask, col=pal(10))
contour(x,y,mask, levels=c(0.5),add=TRUE,drawlabels=FALSE)
lines(path)
polygon(c(xst,xst,Xst,Xst),c(yst,Yst,Yst,yst))
par(mfrow=c(1,1))



#########################################################
##       2. STEP SELECTION BY G-SPACE SAMPLING         ##
#########################################################
library(survival)
library(extraDistr)
# Data-frame creation
n<-200 # Number of controls

dat<-data.frame("x"=md3$Txrd, "y"=md3$Tyrd, "match"=1:nc, "case"=rep(1,nc), "Depth"=md3$Tdepth, "Slib"=md3$Tslib)
for(i in 1:nc)
{
  di<-rrayleigh(n,sg) # Radius selection 
  th<-runif(n,0,2*pi)
  xo<-round(md3$Txrd0[i]+di*cos(th))
  yo<-round(md3$Tyrd0[i]+di*sin(th))
  xstand<-pmax(1,pmin(ndx,xo))
  ystand<-pmax(1,pmin(ndy,yo))
  dat<-rbind(dat, cbind("x"=xo,"y"=yo, "match"=rep(dat$match[i],n),"case"=rep(0,n),"Depth"=u1[cbind(xstand,ystand)],"Slib"=u2[cbind(xstand,ystand)]))
  if(i/1000==round(i/1000)) print(paste(round(100*i/nc),"% done"))
}



# Conditional logistic regression model
mod<-clogit(case~Depth+Slib+strata(match), data=dat)
mod

# Extraction of coefficients and se from model object
ma1<-summary(mod)$coefficients[1,1]
ma2<-summary(mod)$coefficients[2,1]
ea1<-summary(mod)$coefficients[1,3]
ea2<-summary(mod)$coefficients[2,3]

#########################################################
##       3. STEP SELECTION BY E-SPACE APPROXIMATION    ##
#########################################################
comp<-50 # Number of mixture components  
nbrk<-50 # Number of bins for environmental dimensions

### Spatial covariate 1: Bathymetry ###
# Correlogram for spatial layer
cort<-autocor(u1[xst:Xst,yst:Yst],its=50000,bins=100,ranp=0.5)
ro1<-pmax(0,cort$ro)
dx1<-c(cort$lag)

# Constructing the one-dimensional E-space
Espace1<-Mclust(c(u1[lstxs]), G=comp, modelNames="E") # Calculates mixture approximation
#Espace1<-Mclust(dat$Depth, G=comp, modelNames="E") # Calculates mixture approximation
wMi1<-Espace1$parameters$pro # Means of the mixture components
mus1<-Espace1$parameters$mean # Weights of the mixture components
si1<-sqrt(Espace1$parameters$variance$sigmasq)
xss1<-seq(min(u1), max(u1), length.out=nbrk) # Vector of values of environmental variable

### Spatial covariate 2: Sediment ###
# Correlogram for spatial layer
cort<-autocor(u2[xst:Xst,yst:Yst],its=50000,bins=100,ranp=0.5)
ro2<-pmax(0,c(cort$ro))
dx2<-c(cort$lag)

### Constructing the one-dimensional E-space
Espace2<-Mclust(c(u2[lstxs]), G=comp, modelNames="E") # Calculates mixture approximation
#Espace2<-Mclust(dat$Slib, G=comp, modelNames="E") # Calculates mixture approximation
wMi2<-Espace2$parameters$pro # Means of the mixture components
mus2<-Espace2$parameters$mean # Weights of the mixture components
si2<-sqrt(Espace2$parameters$variance$sigmasq)
xss2<-seq(min(u2), max(u2), length.out=nbrk) # Vector of values of environmental variable


# Visualisation 
par(mfrow=c(2,2))
plot(dx1,ro1, type="l", ylim=c(-1,1), xlim=c(0,50), main="Depth autocorrelation")
lines(dx1,ro1)

plot(dx2,ro2, type="l", ylim=c(-1,1), xlim=c(0,50), main="Sediment autocorrelation")
lines(dx2,ro2)

plot(xss1, dens(mus1, wMi1, si1, xss1), type="l", main="Depth marginal availability", xlab="E-space", ylab="Density")
h<-hist(u1[lstxs], plot=FALSE, breaks=nbrk)
points(h$mids, h$counts/sum(h$counts)/(h$mids[2]-h$mids[1])) # Adds histogram of true availability

plot(xss2, dens(mus2, wMi2, si2, xss2), type="l", main="Sediment marginal availability", xlab="E-space", ylab="Density")
h<-hist(u2[lstxs], plot=FALSE, breaks=nbrk)
points(h$mids, h$counts/sum(h$counts)/(h$mids[2]-h$mids[1])) # Adds histogram of true availability
par(mfrow=c(1,1))



if (outrange==0)
{
  ro1In<-ro1
  ro2In<-ro2
  lstxsIn<-lstxs
  dens1In<-dens(mus1, wMi1, si1, xss1)
  dens2In<-dens(mus2, wMi2, si2, xss2)
} else
{
  ro1Out<-ro1
  ro2Out<-ro2
  lstxsOut<-lstxs
  dens1Out<-dens(mus1, wMi1, si1, xss1)
  dens2Out<-dens(mus2, wMi2, si2, xss2)
}


# Fiting model object to connect correlations to initialisation parameters
roTohi1<-RoToHi(xss1, -0.01, .1, nbrk, mus1, si1, wMi1, plot=T, verbose=T)
roTohi2<-RoToHi(xss2, -0.1, 10, nbrk, mus2, si2, wMi2, plot=T, verbose=T)

psi1<-Psi(sg, dx1,ro1, roTohi1, mus1, wMi1)
psi2<-Psi(sg, dx2,ro2, roTohi2, mus2, wMi2)




# Conditional logistic regression model
# ga: Candidate parameter values. First row are 1st order coefs. Second row is 2nd order coefs.
# dat: Covariate data for detections
# wM: Matrix of mixture weights (columns) by variable (rows)
# Mus: Matrix of mixture means (columns) by variable (rows)
# Sis: Vector of mixture st devs 
fjm<-function(datt0,wM,Mus,Sis)
{
  J<-length(datt0[,1]) #Data points
  K<-length(datt0[1,]) #Number of covariates
  fjm<-rep(1,J)
  for(j in 1:J) # Length of data
  {
    for(i in 1:K)
    {
      fjm[j]<-fjm[j]*dens(Mus[i,], wM[i,], Sis[i], datt0[j,i])
    }
  }
  return(fjm)
}



lla<-function(ga)
{
  J<-length(datt[,1]) #Data points
  K<-length(datt[1,]) #Number of covariates
  w1<-exp(ga[1]*datt[,1]+ga[2]*datt[,2]) # Numerator of likelihood
  w2<-w1*0
  L<-length(Mus[1,]) # Number of mixture components
  
  theta<-matrix(NA,nrow=K, ncol=L )
  
  for(i in 1:K) #Covariate loop
  {
    #theta[i,]<-sqrt(1/(1-2*ga[2,i]*Sis[i]^2))*exp(-Mus[i,]^2/(2*Sis[i]^2)+(ga[1,i]*Sis[i]^2+Mus[i,])^2/(2*Sis[i]^2*(1-2*ga[2,i]*Sis[i]^2)))
    theta[i,]<-exp(-Mus[i,]^2/(2*Sis[i]^2)+(ga[i]*Sis[i]^2+Mus[i,])^2/(2*Sis[i]^2))
    }
  
  
  for(j in 1:J) # Length of data
  {
    fj<-1
    for(i in 1:K)
    {
      intsum<-dnorm(datt0[j,i],Mus[i,], Sis[i]) %*% psi[[i]] %*% t(theta)[,i]
      fj<-fj*intsum
    }
    w2[j]<-fj/fjmV[j]
  }
  
  
  return(-sum(log(w1))+sum(log(w2)))
}


wM<-rbind(wMi1,wMi2)
Mus<-rbind(mus1, mus2)
Sis<-c(si1,si2)
psi<- list(psi1,psi2)
datt0<-cbind(md3$Tdepth0,md3$Tslib0)
datt<-cbind(md3$Tdepth,md3$Tslib)
fjmV<-fjm(datt0,wM,Mus,Sis)
optOut<-optim(c(0,0), fn=lla)

if (outrange==0)
{
  a1_0<-optOut$par[1]
  a2_0<-optOut$par[2]
} else
{
  a1_1<-optOut$par[1]
  a2_1<-optOut$par[2]
}


#########################################################
##           4. Exploration of likelihood              ##
#########################################################
ll<-function(a1,a2,dat,nc)
{
  w1<-exp(a1*dat$Depth[1:nc]+a2*dat$Slib[1:nc])
  w2<-exp(a1*dat$Depth[(nc+1):length(dat$x)]+a2*dat$Slib[(nc+1):length(dat$x)])
  w2<-tapply(w2,dat$match[(nc+1):length(dat$x)],sum)
  return(sum(log(w1))-sum(log(w1+w2)))
}



bns<-50
prs1<-seq(ma1-6*ea1,ma1+6*ea1, length.out = bns)
prs2<-seq(ma2-6*ea2,ma2+6*ea2, length.out = bns)
ll1<-matrix(0, bns, bns)
ll2<-matrix(0, bns, bns)

datt0<-cbind(md3$Tdepth0, md3$Tslib0)
datt<-cbind(md3$Tdepth, md3$Tslib)
wM<-rbind(wMi1,wMi2)
Mus<-rbind(mus1, mus2)
Sis<-c(si1,si2)
psi<- list(psi1,psi2)
fjmV<-fjm(datt0,wM,Mus,Sis)

for(i in 1:length(prs1))
{
  for(j in 1:length(prs2))
    {
    ll1[i,j]<--ll(prs1[i],prs2[j],dat,nc)
    ll2[i,j]<-lla(c(prs1[i],prs2[j]))
  }
  print(paste(round(100*i/bns),"% completed"))
}


if (outrange==0)
{
  ll2_0<-ll2
} else
{
  ll2_1<-ll2
}

#########################################################
##     5. Visual comparison between two approaches     ##
#########################################################


par(mfrow=c(4,3))
labelsX<-round(seq(Xmin,Xmax,length.out=length(x)))
labelsY<-round(seq(Ymin,Ymax,length.out=length(y)))
image(labelsX,labelsY,u1/mask, col=pal(30), main="a. Depth", xlab="longitude",ylab="latitude")
contour(labelsX,labelsY,mask, levels=c(0.5),add=TRUE,drawlabels=FALSE, col="grey")

image(labelsX,labelsY,u2/mask, col=pal(30),  main="b. Sediment", xlab="longitude",ylab="latitude")
contour(labelsX,labelsY,mask, levels=c(0.5),add=TRUE,drawlabels=FALSE, col="grey")

plot(path, xlim=range(x), ylim=range(y), main="c. Telemetry", type="l", xlab="longitude",ylab="latitude", axes=F)
polygon(c(xst2,xst2,Xst2,Xst2),c(yst2,Yst2,Yst2,yst2), col="beige",border="grey")
polygon(c(xst1,xst1,Xst1,Xst1),c(yst1,Yst1,Yst1,yst1), col="turquoise",border="grey")
lines(path, col = rgb(red = 0.1, green = 0.1, blue = 0.1, alpha=0.5),lwd=1)
box()


h<-hist(dat$Depth, plot=FALSE, breaks=nbrk)
plot(h$mids, h$counts/sum(h$counts)/(h$mids[2]-h$mids[1]),type="l", main="d. Depth marginal availability", xlab="E-space", ylab="Density") # Adds histogram of true availability
lines(xss1, dens1In, col="brown")
lines(xss1, dens1Out, col="brown", lty=3)

h<-hist(dat$Slib, plot=FALSE, breaks=nbrk)
plot(h$mids, h$counts/sum(h$counts)/(h$mids[2]-h$mids[1]), type="l", main="e. Sediment marginal availability", xlab="E-space", ylab="Density") # Adds histogram of true availability
lines(xss2, dens2In, col="blue")
lines(xss2, dens2Out, col="blue", lty=3)

plot(dx1,ro1In, type="l", ylim=c(0,1), xlim=c(0,10), xlab="Spatial lag (grid cells)", ylab="ro", main="f. Autocorrelation")
lines(dx1,ro1In, col="brown")
lines(dx1,ro1Out, col="brown", lty=2)
lines(dx2,ro2In, col="blue")
lines(dx2,ro2Out, col="blue", lty=2)

h<-exp(mod$coefficients[1]*u1+mod$coefficients[2]*u2)
h<-log(h/sum(h))
h0<-exp(a1_0*u1+a2_0*u2)
h0<-log(h0/sum(h0))
h1<-exp(a1_1*u1+a2_1*u2)
h1<-log(h1/sum(h1))
hr<-c(min(h,h1,h0),max(h,h1,h0))

image(labelsX,labelsY,h/mask, col=pal(30), main="g. SSF by G-space sampling", xlab="longitude",ylab="latitude")
contour(labelsX,labelsY,mask, levels=c(0.5),add=TRUE,drawlabels=FALSE, col="grey")

image(labelsX,labelsY,h0/mask, col=pal(30), main="h. SSF by E-space approx. 1", xlab="longitude",ylab="latitude")
contour(labelsX,labelsY,mask, levels=c(0.5),add=TRUE,drawlabels=FALSE, col="grey")

image(labelsX,labelsY, h1/mask,  col=pal(30), main="i. SSF by E-space approx. 2", xlab="longitude",ylab="latitude")
contour(labelsX,labelsY,mask, levels=c(0.5),add=TRUE,drawlabels=FALSE, col="grey")

cols <- brewer.pal(3, "PuOr")
pal2 <- colorRampPalette(cols) 

image(prs1, prs2,ll1, main="j. Likelihood by G-space sampling", xlab="Depth",ylab="Sediment", col=pal2(10))
CI1<-contour(prs1,prs2, ll1,  add=TRUE, drawlabels=F)
CI1<-contour(prs1,prs2, ll1, levels=c(min(ll1)+1.92), add=TRUE, drawlabels=F, lwd=3, col="white")
abline(v=ma1, lwd=3, col="white")
abline(h=ma2, lwd=3, col="white")

image(prs1, prs2,ll2_0, main="k. Likelihood by E-space approx. 1", xlab="Depth",ylab="Sediment", col=pal2(10))
CI1<-contour(prs1,prs2, ll2_0,  add=TRUE, drawlabels=F)
CI1<-contour(prs1,prs2, ll2_0, levels=c(min(ll2_0)+1.92), add=TRUE, drawlabels=F, lwd=3, col="white")
abline(v=a1_0, lwd=3, col="white")
abline(h=a2_0, lwd=3, col="white")

image(prs1, prs2,ll2_1, main="l. Likelihood by E-space approx. 2", xlab="Depth",ylab="Sediment", col=pal2(10))
CI1<-contour(prs1,prs2, ll2_1,  add=TRUE, drawlabels=F)
CI1<-contour(prs1,prs2, ll2_1, levels=c(min(ll2_1)+1.92), add=TRUE, drawlabels=F, lwd=3, col="white")
abline(v=a1_1, lwd=3, col="white")
abline(h=a2_1, lwd=3, col="white")
par(mfrow=c(1,1))


