
###############################################################################

# import necessary libraries
library(quantreg)
library(Brq)
library(bayesQR)
library(ggplot2)
library(latex2exp)
library(bayesreg)


###############################################################################

# Normal Errors


set.seed(2023)

mu <- 5.0
n <- 100 

x <- seq(from = 0, to = 10, length.out = n)

x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + rnorm(n, 0, 20)
df <- data.frame(x)
df$y <- y
df$y_true <- y_true



# geom_line(aes(x=x, y=y_true, color="red"), show.legend=FALSE)

ggplot(df, aes(x=x, y=y)) + geom_point()  + 
  ggtitle("Data Generation") + xlab("X") + 
  ylab("Y")

ggsave('data_gen_pre.pdf')

SNR <- mean(y)/20






pdf('FRQ.pdf')

set.seed(2023)

# Generate Data

mu <- 5.0
n <- 100 

x <- seq(from = 0, to = 10, length.out = n)

x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + rnorm(n, 0, 20)
df <- data.frame(x)
df$y <- y

plot(x,y, xlab="X", ylab="Y", main="Frequentist Quantile Regression Fit", col=rgb(red = 0, green = 0, blue = 0))


# Plot Data

cols <- c("red", "orange", "green", "orange", "red") 
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)

offsets <- 20*c(-1.65, -0.68, 0, 0.68, 1.65)

for (i in 1:5){
  lines(x, y_true+offsets[i], col="blue")
  
  
  tau <- taus[i]
  model <- rq(formula = y ~ x + x2 + sin(x) + 1, tau=tau, data=df)
  y_pred <-  model$fitted.values
  lines(x, y_pred, lty=2, col=cols[i])
  
  
}

legend(x=-0.30,y=max(y)+0.5,legend=c(.05,.25,.50,.75,.95, "truth"),lty=c(2, 2, 2, 2, 2, 1),
  lwd=c(1,1,1,1,1),col=c("red", "orange", "green", "orange", "red", "blue"),title="Quantile")

dev.off()





## Bayesian Methods





# pdf('BRQ.pdf')

set.seed(2023)

# Generate Data

mu <- 5.0
n <- 100

x <- seq(from = 0, to = 10, length.out = n)
x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + rnorm(n, 0, 20)
df <- data.frame(x)
df$y <- y


# Plot Data

plot(x, y, xlab="X", ylab="Y", main="Bayesian Quantile Regression Fit", col=rgb(red = 0, green = 0, blue = 0))

cols <- c("red", "orange", "green", "orange", "red") 
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)

offsets <- 20*c(-1.65, -0.68, 0, 0.68, 1.65)

for (i in 1:5){
  lines(x, y_true+offsets[i], col="blue")
  
  
  p <- taus[i]
  
  model = Brq(y~x + x2 + sin(x) + 1, tau=p, runs=1500, burn=500)
  y_pred <-  model$fitted.values
  lines(x, y_pred, lty=2, col=cols[i])
  
}


legend(x=-0.30,y=max(y)+0.5,legend=c(.05,.25,.50,.75,.95, "truth"),lty=c(2, 2, 2, 2, 2, 1),
  lwd=c(1,1,1,1,1),col=c("red", "orange", "green", "orange", "red", "blue"),title="Quantile")

# dev.off()






# Chi Squared Errors





set.seed(2023)

mu <- 5.0
n <- 100 

x <- seq(from = 0, to = 10, length.out = n)

x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + 10*rchisq(n, 4)
df <- data.frame(x)
df$y <- y
df$y_true <- y_true





ggplot(df, aes(x=x, y=y)) + geom_point()  + geom_line(aes(x=x, y=y_true, color="red"), show.legend=FALSE) +
  ggtitle("Data Generation") + xlab("X") + 
  ylab("Y")

ggsave('data_gen_pre2.pdf')

SNR <- mean(y)/20






pdf('FRQ2.pdf')

set.seed(2023)

# Generate Data

mu <- 5.0
n <- 100 

x <- seq(from = 0, to = 10, length.out = n)

x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + 10*rchisq(n, 4)
df <- data.frame(x)
df$y <- y

plot(x,y, xlab="X", ylab="Y", main="Frequentist Quantile Regression Fit", col=rgb(red = 0, green = 0, blue = 0))


# Plot Data

cols <- c("red", "orange", "green", "orange", "red") 
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)

offsets <- 10*qchisq(taus, 4)

for (i in 1:5){
  lines(x, y_true+offsets[i], col="blue")
  
  
  tau <- taus[i]
  model <- rq(formula = y ~ x + x2 + sin(x) + 1, tau=tau, data=df)
  y_pred <-  model$fitted.values
  lines(x, y_pred, lty=2, col=cols[i])
  
  
}

legend(x=-0.30,y=max(y)+0.5,legend=c(.05,.25,.50,.75,.95, "truth"),lty=c(2, 2, 2, 2, 2, 1),
  lwd=c(1,1,1,1,1),col=c("red", "orange", "green", "orange", "red", "blue"),title="Quantile")

dev.off()





## Bayesian Methods





pdf('BRQ2.pdf')

set.seed(2023)

# Generate Data

mu <- 5.0
n <- 100

x <- seq(from = 0, to = 10, length.out = n)
x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + 10*rchisq(n,4)
df <- data.frame(x)
df$y <- y


# Plot Data

plot(x, y, xlab="X", ylab="Y", main="Bayesian Quantile Regression Fit", col=rgb(red = 0, green = 0, blue = 0))

cols <- c("red", "orange", "green", "orange", "red") 
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)

offsets <- 10*qchisq(taus, 4)

for (i in 1:5){
  lines(x, y_true+offsets[i], col="blue")
  
  
  p <- taus[i]
  
  model = Brq(y~x + x2 + sin(x) + 1, tau=p, runs=1500, burn=500)
  y_pred <-  model$fitted.values
  lines(x, y_pred, lty=2, col=cols[i])
  
}


legend(x=-0.30,y=max(y)+0.5,legend=c(.05,.25,.50,.75,.95, "truth"),lty=c(2, 2, 2, 2, 2, 1),
  lwd=c(1,1,1,1,1),col=c("red", "orange", "green", "orange", "red", "blue"),title="Quantile")

dev.off()




# Credible Intervals




set.seed(2023)

# Generate Data

mu <- 5.0
n <- 100

x <- seq(from = 0, to = 10, length.out = n)
x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + rnorm(n, 0, 20)
df <- data.frame(x)
df$y <- y


# Plot Data

plot(x, y, xlab="x", main="Scatterplot and Quantile Regression Fit", col=rgb(red = 0, green = 0, blue = 0))

cols <- c("red", "red") 
taus <- c(0.05, 0.95)

offsets <- 20*c(-1.65, 1.65)

func <- function(coeffs, x){
  return(coeffs %*% x)
  
}

X <- cbind(rep(1, length(x)), x, x2, sin(x))

for (i in 1:2){
  lines(x, y_true+offsets[i], col="blue")
  
  
  p <- taus[i]
  x2 <- x^2
  model = Brq(y~x + x2 + sin(x) + 1, tau=p, runs=1500, burn=500)
  y_pred <-  model$fitted.values
  lines(x, y_pred, lty=2, col=cols[i])
  
  
  
  ys <- X %*% t(model$beta) 
  
  uppers <- c()
  lowers <- c()
  
  for (j in 1:100){
    upper <- quantile(ys[j,], probs=c(0.95))
    lower <- quantile(ys[j,], probs=c(0.05))
    
    uppers <- c(uppers, upper)
    lowers <- c(lowers, lower)
  }
  
  lines(x, uppers)
  lines(x, lowers)

  
  
}








set.seed(2023)

# Generate Data

mu <- 5.0
n <- 100

x <- seq(from = 0, to = 10, length.out = n)
x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + 10*rchisq(n, 4)
df <- data.frame(x)
df$y <- y


# Plot Data

plot(x, y, xlab="x", main="Scatterplot and Quantile Regression Fit", col=rgb(red = 0, green = 0, blue = 0))

cols <- c("red", "red") 
taus <- c(0.05, 0.95)

offsets <- 10*qchisq(taus, 4)

func <- function(coeffs, x){
  return(coeffs %*% x)
  
}

X <- cbind(rep(1, length(x)), x, x2, sin(x))

for (i in 1:2){
  lines(x, y_true+offsets[i], col="blue")
  
  
  p <- taus[i]
  x2 <- x^2
  model = Brq(y~x + x2 + sin(x) + 1, tau=p, runs=1500, burn=500)
  y_pred <-  model$fitted.values
  lines(x, y_pred, lty=2, col=cols[i])
  
  
  
  ys <- X %*% t(model$beta) 
  
  uppers <- c()
  lowers <- c()
  
  for (j in 1:100){
    upper <- quantile(ys[j,], probs=c(0.95))
    lower <- quantile(ys[j,], probs=c(0.05))
    
    uppers <- c(uppers, upper)
    lowers <- c(lowers, lower)
  }
  
  lines(x, uppers)
  lines(x, lowers)

  
  
}






set.seed(2023)

# Generate Data

mu <- 5.0
n <- 100

x <- seq(from = 0, to = 10, length.out = n)
x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + rnorm(n, 0, 20)
df <- data.frame(x)
df$y <- y


# Plot Data

plot(x, y, xlab="x", main="Scatterplot and Quantile Regression Fit", col=rgb(red = 0, green = 0, blue = 0))

cols <- c("red", "orange", "green", "orange", "red") 
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)

offsets <- 20*c(-1.65, -0.68, 0, 0.68, 1.65)

pr <- prior(y~x + x2 + sin(x) + 1, data=df, beta0=rep(1,4), V0=diag(4), shape0=1, scale0=1)

str(pr)


out = bayesQR(y~x + x2 + sin(x) + 1, data=df, taus, ndraw=1500, prior=pr)

sum <- summary(out, burnin=50)
for (i in 1:length(sum)){
  beta0 <- sum[[i]]$betadraw[1,1]
  beta1 <- sum[[i]]$betadraw[2,1]
  beta2 <- sum[[i]]$betadraw[3,1]
  beta3 <- sum[[i]]$betadraw[4,1]
  
lines(x, beta0+x*beta1+beta2*x*x+beta3*sin(x),col=cols[i])
}

hist(out[[1]][["betadraw"]])
hist(out[[2]][["betadraw"]])
hist(out[[3]][["betadraw"]])
hist(out[[4]][["betadraw"]])



# Sensitivity Analysis



set.seed(2023)

# Generate Data

mu <- 5.0
n <- 100

x <- seq(from = 0, to = 10, length.out = n)
x2 <- x^2

y_true <- mu + x + x2 + 20*sin(x)

y <- mu + x + x2 + 20*sin(x) + rnorm(n, 0, 20)
df_og <- data.frame(x, y)

df <- data.frame(val=character(), x=numeric(), y=numeric(), ys1=numeric(), ys2=numeric())


# Plot Data

cols <- c("red", "orange", "green", "orange", "red") 
cols2 <- c("red", "green", "blue")
taus <- c(0.05, 0.25, 0.5, 0.75, 0.95)

df_hists <- data.frame()

vals <- c(1, 10, 100)
names <- c("1", "10", "100")

for (i in 1:length(vals)){
  val <- vals[i]
  pr <- prior(y~x + x2 + sin(x) + 1, data=df_og, beta0=rep(1,4), V0=val*diag(4))
  out = bayesQR(y~x + x2 + sin(x) + 1, data=df_og, taus, ndraw=1500, prior=pr, normal.approx=TRUE)
  
  new_df <- data.frame(val=names[i], beta0=out[[1]][["betadraw"]][,1], beta1=out[[1]][["betadraw"]][,2], beta2=out[[1]][["betadraw"]][,3], beta3=out[[1]][["betadraw"]][,4])
  
  df_hists <- rbind(df_hists, new_df)
  
  
  

  sum <- summary(out, burnin=50)

  beta0 <- sum[[1]]$betadraw[1,1]
  beta1 <- sum[[1]]$betadraw[2,1]
  beta2 <- sum[[1]]$betadraw[3,1]
  beta3 <- sum[[1]]$betadraw[4,1]
  ys1 <- beta0+x*beta1+beta2*x*x+beta3*sin(x)
  
  beta0 <- sum[[5]]$betadraw[1,1]
  beta1 <- sum[[5]]$betadraw[2,1]
  beta2 <- sum[[5]]$betadraw[3,1]
  beta3 <- sum[[5]]$betadraw[4,1]
  ys2 <- beta0+x*beta1+beta2*x*x+beta3*sin(x)
  
  new_df <- data.frame(val=names[i], x, y, ys1, ys2)
  
  df <- rbind(df, new_df)
  


}






ggplot(df, aes(x=x, y=y, fill=val)) + geom_point(show.legend=FALSE) + 
  geom_line(aes(x=x, y=ys2, fill=val, color=val)) + 
  geom_line(aes(x=x, y=ys1, fill=val, color=val)) + 
  xlab("X") + 
  ylab("Y") + labs(color=TeX('$\\sigma^2$'))


ggsave('sens_analysis_1.pdf')










require(gridExtra)
require(latex2exp)

plot0 <- ggplot(df_hists, aes(x=beta0, fill=val, after_stat(density))) +
  geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + 
  ggtitle(TeX('$\\beta_{0}$ Posterior Distribution')) + xlab(TeX("$\\beta_0$")) +
  ylab("density") + labs(fill=TeX('$\\sigma^2$'))

plot1 <- ggplot(df_hists, aes(x=beta1, fill=val, after_stat(density))) +
  geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + 
  ggtitle(TeX('$\\beta_{1}$ Posterior Distribution')) + xlab(TeX("$\\beta_1$")) + 
  ylab("density") + labs(fill=TeX('$\\sigma^2$'))

plot2 <- ggplot(df_hists, aes(x=beta2, fill=val, after_stat(density))) +
  geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + 
  ggtitle(TeX('$\\beta_{2}$ Posterior Distribution')) + xlab(TeX("$\\beta_2$")) + 
  ylab("density") + labs(fill=TeX('$\\sigma^2$'))

plot3 <- ggplot(df_hists, aes(x=beta3, fill=val, after_stat(density))) +
  geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + 
  ggtitle(TeX('$\\beta_{3}$ Posterior Distribution')) + xlab(TeX("$\\beta_3$")) + 
  ylab("density") + labs(fill=TeX('$\\sigma^2$'))


whole_plot <- grid.arrange(plot0, plot1, plot2, plot3, ncol=2)

ggsave('sens_analysis.pdf', whole_plot)




# Convergence Properties

  
df2 <- df_hists[df_hists$val==1,]

par(mfrow = c(2, 2))

acf(df2$beta0)
acf(df2$beta1)
acf(df2$beta2)
acf(df2$beta3)










