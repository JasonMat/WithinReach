################################
## 2. ILLUSTRATING SCENARIO   ## 
################################


### A two dimensional covariate layer
d<-300 # Dimension of square G-space

rng<-round(d/6):round(5/6*d) # Subsetting the range, to avoid edge effects
bw<-10# Bandwidth for generation of environmental layer (creates the autocorrelation in the environment)
uor<-3000*environ(d,0.01*d^2, bw)+3000*environ(d,0.01*d^2, bw)^0.01+3000*environ(d,0.01*d^2, bw)^0.001 # Generation of environmental layer in G-space
u<-uor[50:250, 50:250]

cort<-autocor(u,its=50000,bins=300,ranp=0.5)
dx<-c(cort$lag) # Spatial lag
ro<-pmax(0,cort$ro) # Corresponding autocorrelation

plot(dx,ro)

### Constructing the one-dimensional E-space
n<-50 # Number of mixture components
Espace<-Mclust(c(u), G=n, modelNames="E") # Calculates mixture approximation
wMi<-Espace$parameters$pro # Means of the mixture components
mus<-Espace$parameters$mean # Weights of the mixture components
si<-sqrt(Espace$parameters$variance$sigmasq) # Shared variance of the mixture components
nbrk<-200 # Number of bins in environmental dimension

### Graphical output
par(mfrow=c(1,3))
# a. Synthetic environmental variable in G-space
image(seq(50,d-50), seq(50,d-50), u, col=terrain.colors(30), main="a") # plotting layer
# b. Autocorrelation plot
plot(dx,ro, xlab="Spatial lag", ylab="Autocorrelation", main="b") 
# c. Approximation of habitat availability
x<-seq(min(u), max(u), length.out=nbrk) # Vector of values of environmental variable
plot(x, dens(mus, wMi, si, x), type="l", main="c", xlab="E-space", ylab="Density")
h<-hist(u, plot=FALSE, breaks=200)
points(h$mids, h$counts/sum(h$counts)/(h$mids[2]-h$mids[1])) # Adds histogram of true availability
par(mfrow=c(1,1))


###########################################
## 3. EXAMPLES OF USAGE OF THE FUNCTIONS ##
###########################################

# A simple illustration of the renorm() function
Hi<-10 # The parameter that introduces correlation in the joint distribution (zero implies no correlation)
w<-init.joint.weights(mus,Hi)
w1<-w # Copy of original weights
w<-renorm(w, wMi) 

# Diagnostic plots of the renormalisation 
par(mfrow=c(2,2))
image(1:n, 1:n, w1, xlab="weights", ylab="weights", main="a. Original weights")
image(1:n, 1:n, w, xlab="weights", ylab="weights", main="b. Normalised weights")
plot(x+(max(x)-min(x))/300, dens(mus, wMi, si, x), type="l", main="c. Marginal distributions", xlab="E-space", ylab="Density")
lines(x, dens(mus, colSums(w), si, x), col="red")
plot(c(wMi,wMi),c((colSums(w)- wMi)/wMi, (rowSums(w)-wMi)/wMi), main="d. Column & row sum residuals", xlab="True marginal weights", ylab="Standardised residuals", ylim=c(-1,1))
par(mfrow=c(1,1))


# Visualisation of full joint distribution
joint.dens<-joint.mixture(x,w,mus,si)
densi<-dens(mus, wMi, si, x)
densi<-densi/sum(densi)
joint.mixture.plot(x, joint.dens,densi)

#####################################
## 4. CALCULATION OF PSI WEIGHTS   ##
#####################################

# Fiting model object to connect correlations to initialisation parameters
roTohi<-RoToHi(x, -0.3, 2, 50, mus, si, wMi, plot=T, verbose=T)

# Calculation of new weights for joint mixture
sg<-10
psi<-Psi(sg, dx, ro, roTohi, mus, wMi)

# Visualisation of full joint distribution using aggregate weights 
# (i.e. distance independent joint distribution)
joint.dens<-joint.mixture(x,psi,mus,si)
densi<-dens(mus, wMi, si, x)
densi<-densi/sum(densi)
joint.mixture.plot(x, joint.dens,densi)
image(joint.dens/densi)

#####################################
## 5.       VALIDATION DATA        ##
#####################################
ssize<-50
joint.dens.real<-matrix(0, nbrk, nbrk)
for(xr in 50:(d-50))
{
  for(yr in 50:(d-50))
  {
    xs<-pmin(d, pmax(1, round(rnorm(ssize, xr, sg))))
    ys<-pmin(d, pmax(1, round(rnorm(ssize, yr, sg))))
    chi<-which.min(abs(x-uor[xr,yr]))
    for( i in 1:ssize)
    {
      ex<-which.min(abs(x-uor[xs[i],ys[i]]))
      joint.dens.real[chi,ex]<-joint.dens.real[chi,ex]+1
    }
    
  }
}
joint.dens.real<-joint.dens.real/sum(joint.dens.real)
densi.real<-rowSums(joint.dens.real)
densi.real<-densi.real/sum(densi.real)
#dev.new()
joint.mixture.plot(x, joint.dens.real,densi.real)
image(joint.dens.real/matrix(densi.real, 200, 200))

par(mfrow=c(1,2))
image(joint.dens, zlim=range(joint.dens), main="High mobility - Model")
contour(joint.dens, add=T, levels=seq(0,0.0001,0.00001))
image(joint.dens.real, zlim=range(joint.dens), main="High mobility - Realisation")
contour(joint.dens.real, add=T, levels=seq(0,0.0001,0.00001))
par(mfrow=c(1,1))




