library(lme4)

InsampleTown <- subset(Insample, TypeTown==1)

# The negative binomial frequency model: N
N <- InsampleTown$Freq
freq_lik <- function(parm) {
  r    <- parm[1]
  beta <- parm[2]
  lik  <- -sum(dnbinom(N, size=r, prob=1/(1+beta), log=TRUE))
  return(lik)
}
init.parm.freq <- c(mean(N)/(var(N)/mean(N)-1), var(N)/mean(N)-1) # initial estimates by method of moments
freq_mod       <- optim(par=init.parm.freq, fn=freq_lik) # Maximum likelihood estimation for the frequency model

r.est    <- freq_mod$par[1]
beta.est <- freq_mod$par[2]

dlomax <- function(x, shape=1, scale=1, log=FALSE) {
  if (log==FALSE) result <- shape*scale^shape / (scale + x)^(shape+1)
  if (log==TRUE ) result <- log(shape)+shape*log(scale)-(shape+1)*log(scale + x)
  return(result) }

rlomax <- function(n, shape=1, scale=1) {
  u <- runif(n)
  result <- scale*((1 - u)^(-1/shape)-1)
  return(result) }

# The Pareto severity model: Xbar
Xbar <- InsampleTown$yAvg[which(InsampleTown$Freq>0)]
sev_lik <- function(parm) {
  alpha <-  parm[1]
  theta <-  parm[2]
  lik   <- -sum(dlomax(Xbar, shape=alpha, scale=theta, log=TRUE))
  return(lik)
}
init.parm.sev <- c( 2/(1-mean(Xbar)^2/var(Xbar))  , mean(Xbar)*(2/(1-mean(Xbar)^2/var(Xbar))-1) ) # initial estimates by method of moments
sev_mod       <- optim(par=init.parm.sev, fn=sev_lik, method="L-BFGS-B") # Maximum likelihood estimation for the severity model
alpha.est     <- sev_mod$par[1]
theta.est     <- sev_mod$par[2]


