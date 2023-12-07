#### Ordinary Bayesian Quantile Regression

library(bayesQR)

#set.seed(10108)

### Parameter and Data Initialization
tau <- 0.75
I<-50   # no. of subjects
J<-25    # no. of longitudinal obs. per subject


### ----------------------------
### Ordinary Bayesian Quantile Regression
### ----------------------------
N <- 4
S <- 40000
betas <- c(1,5)

beta0.means <- c()
beta1.means <- c()
for (nn in 1:N){
  print(paste0("At step: ", nn))
  
  xx <- matrix(runif(I*J,0,1), nrow=I, ncol=J, byrow=TRUE)   # all the x_ij's
  betas <- c(1,5)
  alphas <- mvrnorm(n=I,mu=c(0,0),Sigma=diag(c(1,1)))
  y.data <- matrix(NA, nrow=I, ncol=J)
  for (i in 1:I){
    for (j in 1:J){
      err <- rnorm(n=1,mean=0,sd=sqrt(0.1)) #!!! ERROR DISTRIBUTION
      #err <- rt(n=1,df=3) #!!! ERROR DISTRIBUTION
      y.help <- betas[1] + betas[2]*xx[i,j] + alphas[i,1] + alphas[i,2]*xx[i,j] + err
      y.data[i,j] <- y.help
    }
  }
  xx <- as.vector(xx)
  y.data <- as.vector(y.data)
  
  bqr <- bayesQR(y.data~xx, quantile=c(tau), ndraw=S, )
  sum <- summary(bqr, burnin=5000)
  
  beta0.means <- append(beta0.means,sum[[1]]$betadraw[1,1])
  beta1.means <- append(beta1.means,sum[[1]]$betadraw[2,1])
}

mean(beta0.means)
mean(beta1.means)

hist(beta0.means)
hist(beta1.means)

write.table(cbind(beta0.means, beta1.means), file="BQR_0.25_t.txt", row.names=FALSE)

  
  
  
  
