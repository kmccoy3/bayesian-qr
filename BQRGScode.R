#####################################################################
##### ---- Bayesian Quantile Regression Gibbs Sampler Code ---- #####
##### Note: Algorithm based on "Bayesian quantile regression    #####
##### for longitudinal data models" by Luo et al. (2011)        #####
#####################################################################

library(MASS)
library(GIGrvg) # different rgig function
library(GeneralizedHyperbolic) # different rgig function
library(invgamma)

#set.seed(10152)

### Parameter and Data Initialization
tau <- 0.25
I<-50   # no. of subjects
J<-25    # no. of longitudinal obs. per subject
xx <- matrix(runif(I*J,0,1), nrow=I, ncol=J, byrow=TRUE)   # all the x_ij's
betas <- c(1,5)
alphas <- mvrnorm(n=I,mu=c(0,0),Sigma=diag(c(1,1)))

y.data <- matrix(NA, nrow=I, ncol=J)
for (i in 1:I){
  for (j in 1:J){
    err <- rnorm(n=1,mean=0,sd=sqrt(0.1)) #!!! ERROR DISTRIBUTION
    y.help <- betas[1] + betas[2]*xx[i,j] + alphas[i,1] + alphas[i,2]*xx[i,j] + err
    y.data[i,j] <- y.help
  }
}

### ----------------------------
### Gibbs sampler algorithm:
### ----------------------------
S <- 25000
num.rand.params <- 2
c0<-0.01; d0<-0.01
k0<-0.01; w0<-0.01
k1 <- (1-(2*tau))/(tau*(1-tau))
k2 <- 2/(tau*(1-tau))
B0 <- 100*diag(c(1,1))
b0.vec <- c(0,0)

### 1: Generate initial values sigma.0, alpha.i.0, beta.0
betas.P <- matrix(ncol=2,nrow=0)    # beta0, beta1
alpha0.P <- matrix(ncol=I,nrow=0)   # alpha0.i
alpha1.P <- matrix(ncol=I,nrow=0)   # alpha1.i
sigma.P <- matrix(ncol=1,nrow=0)    # sigma
phi2.P <- matrix(ncol=1,nrow=0)     # phi2
ee.P <- matrix(ncol=(I*J),nrow=0)   # e.ij

# current values
betas.c <- mvrnorm(n=1,mu=c(0,3),Sigma=diag(c(1,1)))
alpha0.c <- mvrnorm(n=1,mu=rep(0,I),Sigma=diag(nrow=I))
alpha1.c <- mvrnorm(n=1,mu=rep(0,I),Sigma=diag(nrow=I))
sigma.c <- 0.5
phi2.c <- 0.5
ee.c <- matrix(rexp(n=I*J,rate=phi2.c), nrow=I, ncol=J, byrow=TRUE)

# save post. samples
betas.P <- rbind(betas.P, betas.c)
alpha0.P <- rbind(alpha0.P, alpha0.c)
alpha1.P <- rbind(alpha1.P, alpha1.c)
sigma.P <- rbind(sigma.P, sigma.c)
phi2.P <- rbind(phi2.P, phi2.c)
ee.P <- ee.c

ee11 <- c()
ee114 <- c()

for (s in 2:S){
  if (s%%1000==0){print(paste0(s, " samples drawn"))}
  
  ### 2: Generate eij
  gamma2.ij <- ((k1^2)/(k2*sigma.c))+(2/sigma.c)
  for (i in 1:I){
    for (j in 1:J){
      delta2.ij <- (y.data[i,j] - betas.c[1] - betas.c[2]*xx[i,j] - alpha0.c[i] - alpha1.c[i]*xx[i,j])^2/(k2*sigma.c)
      #ee.c[i,j] <- rgig(1, chi=sqrt(delta2.ij), psi=sqrt(gamma2.ij), lambda=0.5) # Update current val TRY A DIFFERENT GIG Random GENERATOR
      ee.c[i,j] <- GIGrvg::rgig(n=1, lambda=0.5, chi=delta2.ij, psi=gamma2.ij)
      if (i==1 & j==1){
        ee11 <- append(ee11, ee.c[i,j])
      }
      if (i==11 & j==4){
        ee114 <- append(ee114, ee.c[i,j])
      }
    }
  }
  ### 3: Generate sigma
  nu <- ((I*J)/2)+c0
  omega.help <-0
  for (i in 1:I){
    for (j in 1:J){
      val <- (y.data[i,j] - betas.c[1] - betas.c[2]*xx[i,j] - alpha0.c[i] - alpha1.c[i]*xx[i,j] - (k1*ee.c[i,j]))^2
      val <- val/ee.c[i,j]
      omega.help <- omega.help+val
    }
  }
  omega.help <- omega.help/(2*k2)
  omega <- d0 + omega.help
  #sigma.c <- 1/rgamma(1,nu,omega) # Update current val
  sigma.c <- rinvgamma(1,nu,omega)
  
  ### 4: Generate phi2
  arg1 <- (I*num.rand.params/2)+k0
  sum.help <- 0
  for (i in 1:I){
    sum.help <- sum.help + alpha0.c[i]^2 + alpha1.c[i]^2
  }
  arg2 <- w0 + 0.5*sum.help
  #phi2.c <- 1/rgamma(1,arg1,arg2) # Update current val
  phi2.c <- rinvgamma(1,arg1,arg2)
  
  ### 5 Generate alpha.i
  for (i in 1:I){
    a0.help<-0; a1.help<-0; Amat.help<-0
    for (j in 1:J){
      a0.help <- a0.help + ((1)*(y.data[i,j] - betas.c[1] - (betas.c[2]*xx[i,j]) - (k1*ee.c[i,j]))/ee.c[i,j])
      a1.help <- a1.help + ((xx[i,j])*(y.data[i,j] - betas.c[1] - betas.c[2]*xx[i,j] - (k1*ee.c[i,j]))/ee.c[i,j])
      mat.help <- matrix( c(1,xx[i,j],xx[i,j],xx[i,j]^2), ncol=num.rand.params, nrow=num.rand.params, byrow=TRUE)
      Amat.help <- Amat.help + (mat.help/ee.c[i,j])
    }
    a0.help <- a0.help/(k2*sigma.c)
    a1.help <- a1.help/(k2*sigma.c)
    Amat.inv <- (Amat.help/(k2*sigma.c)) + (diag(num.rand.params)/phi2.c)
    
    alpha.i.c <- mvrnorm(n=1, solve(Amat.inv)%*%matrix(c(a0.help,a1.help),ncol=1), solve(Amat.inv))
    alpha0.c[i] <- alpha.i.c[1] # Update current val
    alpha1.c[i] <- alpha.i.c[2] # Update current val
  }
  
  ### 6 Generate beta
  b0.help<-0; b1.help<-0; Bmat.help<-0
  for (i in 1:I){
    for (j in 1:J){
      b0.help <- b0.help + ((1)*(y.data[i,j] - alpha0.c[i] - alpha1.c[i]*xx[i,j] - (k1*ee.c[i,j]))/ee.c[i,j])
      b1.help <- b1.help + ((xx[i,j])*(y.data[i,j] - alpha0.c[i] - alpha1.c[i]*xx[i,j] - (k1*ee.c[i,j]))/ee.c[i,j])
      mat.help <- matrix( c(1,xx[i,j],xx[i,j],xx[i,j]^2), ncol=num.rand.params, nrow=num.rand.params, byrow=TRUE)
      Bmat.help <- Bmat.help + (mat.help/ee.c[i,j])
    }
  }
  Bmat.inv <- (Bmat.help/(k2*sigma.c)) + solve(B0)
  b0.help <- b0.help/(k2*sigma.c) + solve(B0)[1,]%*%b0.vec
  b1.help <- b1.help/(k2*sigma.c) + solve(B0)[2,]%*%b0.vec # b0.vec refers to the prior means
  betas.c <- mvrnorm(n=1, solve(Bmat.inv)%*%matrix(c(b0.help,b1.help),ncol=1), solve(Bmat.inv))
  
  ### 7 Save all current samples (except the e.ij's)
  betas.P <- rbind(betas.P, betas.c)
  alpha0.P <- rbind(alpha0.P, alpha0.c)
  alpha1.P <- rbind(alpha1.P, alpha1.c)
  sigma.P <- rbind(sigma.P, sigma.c)
  phi2.P <- rbind(phi2.P, phi2.c)
  ee.P <- ee.c
}

plot(betas.P[,1], type='l')
plot(alpha0.P[,6], type='l')
plot(sigma.P[,1], type='l')
plot(phi2.P[,1], type='l')
plot(ee114, type="l")

hist(betas.P[5000:S,1], breaks='scott')
hist(betas.P[5000:S,2], breaks='scott')
mean(betas.P[5000:S,1])
mean(betas.P[5000:S,2])

acf(betas.P[,2])
