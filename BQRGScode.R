library(MASS)
library(invgamma)

#set.seed(10107)

### Parameter and Data Initialization
tau <- 0.75
I<-20   # no. of subjects
J<-5    # no. of longitudinal obs. per subject
xij <- matrix(runif(I*J,0,1), nrow=I, ncol=J)
betas <- c(1,5)
alphas <- mvrnorm(n=I,mu=c(0,0),Sigma=diag(c(1,1)))

y.data <- matrix(NA, nrow=I, ncol=J)
for (i in 1:I){
  for (j in 1:J){
    err <- rnorm(n=1,mean=0,sd=sqrt(0.1)) #!!! ERROR DISTRIBUTION
    y.help <- betas[1] + betas[2]*xij[i,j] + alphas[i,1] + alphas[i,2]*xij[i,j] + err
    y.data[i,j] <- y.help
  }
}

### ----------------------------
### Metropolis-Hastings sampling
### ----------------------------
S <- 40000
no.params <- 2
sig2.beta <- 0.005        # Used to control acceptance rate of betas
sig2.alpha <- 0.05         # Used to control acceptance rate of alphas
B0 <- 100*diag(c(1,1))
c0<-0.01; d0<-0.01
k0<-0.01; w0<-0.01
alphas.accept<-0; betas.accept<-0

# 1. Initialize first parameter values
PARAMS <- matrix(ncol=4,nrow=0)
ALPHAS.0 <- matrix(ncol=I,nrow=0)
ALPHAS.1 <- matrix(ncol=I,nrow=0)

betas.1 <- mvrnorm(n=1,mu=c(0,3),Sigma=diag(c(1,1)))    # Indices: 1-2
sigma.1 <- 0.5                                          # Indices: 3
phi2.1 <- 0.5                                           # Indices: 4
PARAMS <- rbind(PARAMS, c(betas.1,sigma.1,phi2.1))
ALPHAS.0 <- rbind(ALPHAS.0, mvrnorm(n=1,mu=rep(0,I),Sigma=diag(nrow=I)))
ALPHAS.1 <- rbind(ALPHAS.1, mvrnorm(n=1,mu=rep(0,I),Sigma=diag(nrow=I)))

### --- 2. Define useful functions
double.summer <- function(betas.in, alphas0.in, alphas1.in){
  val <- 0
  for (i in 1:I){
    for (j in 1:J){
      val.help <- y.data[i,j] - 
        betas.in[1] - betas.in[2]*xij[i,j] - 
        alphas0.in[i] - alphas1.in[i]*xij[i,j]
      val.help <- val.help*(tau-ifelse(val.help<=0,1,0))
      val <- val + val.help
    }
  }
  return(val)
}

single.summer <- function(betas.in, alphas.in, i.in){
  val <- 0
  for (j in 1:J){
    val.help <- y.data[i.in,j] - betas.in[1] - betas.in[2]*xij[i.in,j] - alphas.in[1] - alphas.in[2]*xij[i.in,j]
    val.help <- val.help*(tau-ifelse(val.help<=0,1,0))
    val <- val + val.help
  }
  return(val)
}


### --- 3. MH Sampling
for (s in 2:S){
  if (s%%1000==0){print(paste0("At step:", s))}
  # 3.a Update betas
  betas.p <- mvrnorm(n=1, mu=PARAMS[s-1,1:2], Sigma=(sig2.beta*diag(c(1,1))))
  betas.c <- PARAMS[s-1,1:2]
  
  r <- exp( (-t(betas.p)%*%solve(B0)%*%betas.p/2)-(double.summer(betas.p,ALPHAS.0[s-1,],ALPHAS.1[s-1,])/PARAMS[s-1,3]) ) /
    exp( (-t(betas.c)%*%solve(B0)%*%betas.c/2)-(double.summer(betas.c,ALPHAS.0[s-1,],ALPHAS.1[s-1,])/PARAMS[s-1,3]) )
  
  print(paste0("beta r:", r))
  if (runif(1)<r) {betas.c <- betas.p; betas.accept<-betas.accept+1}
  
  # 3.b Update alphas
  alphas0.c <- c()
  alphas1.c <- c()
  for (i in 1:I){
    alphas.p <- mvrnorm(n=1, mu=c(ALPHAS.0[s-1,i],ALPHAS.1[s-1,i]), Sigma=(sig2.alpha*diag(c(1,1))))
    alphas.c <- c(ALPHAS.0[s-1,i],ALPHAS.1[s-1,i])
    
    r <- exp( (-t(alphas.p)%*%alphas.p/(2*PARAMS[s-1,4])) - (single.summer(betas.c,alphas.p,i)/PARAMS[s-1,3])) /
      exp( (-t(alphas.c)%*%alphas.c/(2*PARAMS[s-1,4])) - (single.summer(betas.c,alphas.c,i)/PARAMS[s-1,3]))
    
    print(paste0("alpha r:", r))
    
    if (runif(1)<r){
      alphas0.c <- append(alphas0.c,alphas.p[1])
      alphas1.c <- append(alphas1.c,alphas.p[2])
      alphas.accept <- alphas.accept+1}
    else{
      alphas0.c <- append(alphas0.c,alphas.c[1])
      alphas1.c <- append(alphas1.c,alphas.c[2])}
  }
  
  # 3.c Update sigma
  arg1 <- I*J+c0
  arg2 <- d0 + double.summer(betas.c, alphas0.c, alphas1.c)
  sigma.c <- rinvgamma(1,arg1,arg2)
  
  # 3.d Update phi2
  arg1 <- (I*no.params/2) + k0
  arg2 <- w0 + 0.5*(sum(alphas0.c^2)+sum(alphas1.c^2))
  phi2.c <- rinvgamma(1,arg1,arg2)
  
  # 3.d Update holding variables with updated samples
  PARAMS <- rbind(PARAMS, c(betas.c,sigma.c,phi2.c))
  ALPHAS.0 <- rbind(ALPHAS.0, alphas0.c)
  ALPHAS.1 <- rbind(ALPHAS.1, alphas1.c)
}

plot(PARAMS[1:S,1],type='l')
plot(PARAMS[1:S,2],type='l')
betas.accept/S
alphas.accept/S/I

acf(PARAMS[1:S,1], lag.max=100)
acf(PARAMS[seq(1,S,20),1], lag.max=100)
hist(PARAMS[1:S,2], breaks="scott")

hist(PARAMS[seq(5000,S,10),1], breaks="scott")
hist(PARAMS[seq(5000,S,10),2], breaks="scott")

mean(PARAMS[5000:S,1])
mean(PARAMS[5000:S,2])

mean(PARAMS[seq(5000,S,10),1])
mean(PARAMS[seq(5000,S,10),2])
