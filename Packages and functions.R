###################
## 1. FUNCTIONS  ##
###################

## Libraries #################################################
library(boot)
library(mgcv)
library(spatial)  
library(mclust)
library(mvtnorm)
library(pgirmess)
library(plotrix)
library(graphics)

## Main functions ############################################  

# FUNCTION: Re-normalises matrix to given marginal weights 
# Works by repeatedly normalinsing by rows and columns
# Argumemts are:
# w: The initial matrix of joint weights
# wMi: The vector of marginal weights
# iterations: The nuber of normalising iterations
# Five iterations are usually more than enough
renorm<-function(w, wMi, iterations=50)
{
  n<-length(wMi)
  for(i in 1:iterations)
  {
    wc<-matrix(wMi/colSums(w), n, n, byrow=T)
    w<-w*wc
    wr<-matrix(wMi/rowSums(w), n, n)
    w<-w*wr
  }
  return(w)
}

# FUNCTION: Generation of initial matrix of weights
# Creates an initial matrix of weights exponentially decline away from the diagonal
# The rate of decay with distance between the component means is 
# determined by the distance between any two means and 
# the decay parameter Hi
# Arguments are:
# mus: The vector of mixture component means
# Hi: The decay parameter. 
# This has default of 0 (no decay with distance, leading to uniformly weighted joint components)
init.joint.weights<-function(mus, Hi=0)
{
  
  n<-length(mus)
  comb<-expand.grid(mus,mus)
  dcomb<-comb[,1]-comb[,2]
  w<-matrix(exp(-Hi*dcomb^2), n, n)
  w<-w/sum(w)
  return(w)
}


# FUNCTION: Calculates full joint mixture distribution on a grid of particular resolution
# x: vector of grid tickpoints
# w: matrix of joint mixture weights
# mus: vector of means for marginal mixture
# si: standard deviation of mixture components
joint.mixture<-function(x, w, mus, si)
{
  bins<-length(x)
  joint.dens<-matrix(0, bins, bins)
  n<-length(mus)
  locs<-as.matrix(expand.grid(x,x))
  
  for(k1 in 1:n)
  {
    for(k2 in 1:n)
    {
      joint.dens<-joint.dens+w[k1,k2]*dmvnorm(locs, mean=c(mus[k1],mus[k2]), sigma=diag(c(si^2,si^2)), log=FALSE)
    }
  }
  joint.dens<-joint.dens/sum(joint.dens)
  return(joint.dens)
}

# FUNCTION: Fitting a relationship between ro and Hi
# Generates an empirical relationship between Hi, the 
# parameter that is used to initialise the correlation 
# structure in the weights matrix, and ro, the eventual
# correlarion generated in the joint distribution
# Arguments are:
# x: The vector of environmental values to be used
# Himin, Himax: The minimum/maximum values to be used for Hi
# Histps: The number of steps to be used through that range
# mus: The means of the marginal mixture components
# si: The standard deviation of the mixture components
# wMi: The weights of the mixture components
RoToHi<-function(x, Himin, Himax, Histps,  mus, si, wMi, plot=F, verbose=F)
{  
  if(Himin<0) print("WARNING: Initialisation value h for distance dependence is negative, this can cause unstable results")
  Hi<-seq(Himin,Himax, (Himax-Himin)/Histps)
  data<-cbind(Hi,"ro"=0*Hi)
  
  for(h in 1:length(Hi))
  {
    if(verbose==T) print(paste(h, " out of", length(Hi)))
    w<-init.joint.weights(mus,Hi[h])
    w<-renorm(w, wMi) 
    joint.dens<-joint.mixture(x,w,mus,si)
    data[h,2]<-corr(expand.grid(x,x), c(joint.dens))
  }
  
  
  data<-as.data.frame(data)
  plot(data$ro, data$Hi)
  RoToHi<-gam(Hi~s(ro, k=ceiling(Histps)), data, family=gaussian()) # Model approximating the relationship between decay parameter in initialisation of weights and resulting correlation structure
  
  if (plot==T)
  {
    plot(data$ro, data$Hi)
    lines(data$ro, fitted(RoToHi))
  }
  return(RoToHi)
}
# Himin<-0
# Himax<-1
# Histps<-30
# FUNCTION: Calculation of new weights for joint mixture
# These are the aggregate weights with distance integrated out
# Arguments include:
# sg: The standard deviation of the mobility constraint
# dx: The vector of distaces to be used for the approximation to the integral
# ro: The vector of environmental auto-correlatons observed at these distances
# roTohi: A model object linking initialisation parameters for the weight matrix with correlations

Psi<-function(sg, dx, ro, roTohi, mus, wMi)
{
  ni<-length(mus) # Number of mixture components
  Kg<-dx/sg^2*exp(-dx^2/(2*sg^2)) # Mobility constraint, a Rayleigh distribution
  His<-predict(roTohi, newdata=as.data.frame(ro)) # Generates initialisation values that correspond to desired correlation values
  Psi<-matrix(0, ni, ni) # Will store new, distance-integrated weights for joint mixture
  for(dI in 2:length(dx)) # Loops through the distance steps
  {
    Hi<-His[dI] # Picks appropriate initialisation parameter that will generate the correct correlation for this distance
    w<-init.joint.weights(mus,Hi) # Initialisation according to this parameter value
    w<-renorm(w, wMi) # Renormalisation of the joint weights
    Psi<-Psi+w*Kg[dI]*(dx[dI]-dx[dI-1]) # Combination of envirnmental correlation with mobility constraint at the current distance
    
  }
  Psi<-renorm(Psi, wMi) # Normalisation of the aggregate weights
  return(Psi)
}

## Helper functions ##########################################



# FUNCTION: Generates a random environmental layer in a 
# square dxd arena using a total of x seed points
# Arguments are:
# d: The dimensions of the square grid
# x: The number of seeds/kernels
# bw: The kernel bandwidth
environ<-function(d,x,bw,seeds=NULL)
  
{
  ar<-array(0, dim=c(d,d))
  # Places seeds in arena
  if(length(seeds)==0) seeds<-cbind(runif(x, min=1, max=d), runif(x, min=1, max=d))
  
  # Smooths seeds to create spatial autocorrelation 
  require(KernSmooth)
  sarx<-bkde2D(seeds, bandwidth = c(bw,bw), gridsize=c(d,d),
               range.x=list(c(1,d),c(1,d)))
  sarx2<-bkde2D(seeds, bandwidth = c(2*bw,2*bw), gridsize=c(d,d),
                range.x=list(c(1,d),c(1,d)))
  lay<-sarx$fhat+sarx2$fhat
  lay<-lay/sum(lay)
  return(lay)
}


# FUNCTION: Generates probability density of marginal mixture
# Arguments are:
# mus: The vector of means for the mixture
# wMi: The vector of weights for the mixture
# si: The standard deviation of mixture components
# x: The vector of positions at which the mixture likelihood is required
dens<-function(mus, wMi, si, x)
{
  ll<-rep(0,length(x))
  for(i in 1:length(x))
  {
    n<-length(mus)
    ll[i]<-wMi%*%dnorm(x[i],mus, si)
  }

  return(ll)
}



# FUNCTION: Visualises mixture distribution and its marginals
# Arguments are:
# x: The tickpoints of the marginal and joint grids
# joint: The matrix of  joint densities on the grid
# marginal: The vector of marginal densities on the grid
joint.mixture.plot<-function(x, joint, marginal)
{
  par(mar=c(.1,.1,.1,.1))
  layout(matrix(c(2,1,0,3), 2, 2), widths=c(6,1),heights=c(1,6), respect=TRUE)
  
  image(x,x, joint, xlab=NA, ylab=NA, axes=F)
  contour(x,x,joint, add=T, levels=seq(0,0.0001,0.00001))
  box()
  plot(x, marginal, type="l", xlab=NA, axes=F)
  box()
  plot(marginal[length(x):1], x[length(x):1], type="l", ylab=NA, axes=F)
  box()
  layout(matrix(c(1), 1,1 ))
  par(mar=c(5,4,4,2))
}

# FUNCTION: Calculates autocorrelation from a map of values in matrix form
# Arguments are:
# u: The matrix of values
# its: The number of pairs of cells to be used for each distance bin
# bins: The number of bins to be used for subdividing the distance axis
# ranp: The proportion of the full range to be used on the distance axis. The default of 0.5 goes from 0 to half the maximum distance.
autocor<-function(u,its=5000,bins=100,ranp=0.5)
{
  xl<-length(u[,1])
  yl<-length(u[1,])
  ran<-ceiling(ranp*sqrt((xl)^2+(yl)^2))
  binw<-ran/bins
  binmid<-seq(binw/2,ran-binw/2,by=binw)
  ro<-binmid*0

  for(i in 1:length(binmid))
  {
    dat<-matrix(0,nrow=its,ncol=2)
    x0<-runif(its,1,xl)
    y0<-runif(its,1,yl)
    th<-runif(its,0,2*pi)
    x1<-x0+binmid[i]*cos(th)
    y1<-y0+binmid[i]*sin(th)
    inside<-as.logical((x1<=xl)*(x1>=1)*(y1<=yl)*(y1>=1))
    z1<-u[cbind(round(x0[inside]),round(y0[inside]))]
    z2<-u[cbind(round(x1[inside]),round(y1[inside]))]
    ro[i]<-corr(cbind(z1,z2))
  }
  dt<-data.frame("lag"=binmid,"ro"=ro)
  return(dt)
}

