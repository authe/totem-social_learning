library(FKF)
library(zoo)

K <- 100  # number of simulations

# set parameters

T <- 150             # number of submission dates
N <- 20               # number of submitters
rho <- 0.8            # AR1 coefficient for underlying: theta_t = rho theta_{t-1} + e_t
theta.bar <- 0.2      # LR mean of theta_t
sig2.e <- 0.003       #(1-rho^2)*(0.1)^2         # variance of fundamental shock e_t ~ N(0, sig2.e)
sig2.u <- (0.025)^2    # variance of shock to public signal S_t = a + b theta_t + u_t with u_t ~ N(0,sig2.u)
sig2.v <- (0.02)^2          
sig2.z <- 0     #(0.005)^2    # variance of measurement error
a <- 0.1              # public signal is S_t = a + b theta_t + u_t
b <- 0.9

paras.sim <- c(T=T, N=N, rho=rho, theta.bar=theta.bar, sig2.u=sig2.u, sig2.v=sig2.v, sig2.e=sig2.e, a=a, b=b, sig2.z=sig2.z)

# derived parameters for steady state Kalman filter
# steady state variance of theta_t
sig2 <- (2*b^2*sig2.e^2 + rho^2*sig2.e*sig2.u + 2*b^2*sig2.e*sig2.v + b^2*rho^2*sig2.e*sig2.v - rho^2*sig2.u*sig2.v + rho^4*sig2.u*sig2.v + rho^2*sqrt(sig2.e^2*sig2.u^2 + 2*b^2*sig2.e^2*sig2.u*sig2.v + 2*sig2.e*sig2.u^2*sig2.v + 2*rho^2*sig2.e*sig2.u^2*sig2.v + b^4*sig2.e^2*sig2.v^2 + 2*b^2*sig2.e*sig2.u*sig2.v^2 + 2*b^2*rho^2*sig2.e*sig2.u*sig2.v^2 + sig2.u^2*sig2.v^2 - 2*rho^2*sig2.u^2*sig2.v^2 + rho^4*sig2.u^2*sig2.v^2))/(2*(b^2*sig2.e + rho^2*sig2.u + b^2*sig2.v))
# steady state Kalman gains
k1 <- (b*rho*(sig2 - sig2.e)*sig2.v)/(rho^2*sig2.u*(sig2 + sig2.v) - b^2*(-sig2 + sig2.e)*(sig2.e + sig2.v))
k2 <- (b^2*sig2.e*(-sig2 + sig2.e) - rho^2*sig2*sig2.u)/(-(rho^2*sig2.u*(sig2 + sig2.v)) + b^2*(-sig2 + sig2.e)*(sig2.e + sig2.v))
k <- k2 + (b/rho)*k1

# initial conditions for optimization
paras0 = c(trans.rho=1, theta.bar=0.15, log.sig2.u = -5, log.sig2.v=-5, log.sig2.e=-5, a=0, b=1)


# define model

# create state space model from parameters

SP.model <- function(trans.rho, theta.bar, log.sig2.u, log.sig2.v, log.sig2.e, a, b, N){
  
  rho <- trans.rho^2/(1+trans.rho^2)
  sig2.e <- exp(log.sig2.e)
  sig2.u <- exp(log.sig2.u)
  sig2.v <- exp(log.sig2.v)
  sig2.z <- 0   # no measurement error
  
  sig2 <- (2*b^2*sig2.e^2 + rho^2*sig2.e*sig2.u + 2*b^2*sig2.e*sig2.v + b^2*rho^2*sig2.e*sig2.v - rho^2*sig2.u*sig2.v + rho^4*sig2.u*sig2.v + rho^2*sqrt(sig2.e^2*sig2.u^2 + 2*b^2*sig2.e^2*sig2.u*sig2.v + 2*sig2.e*sig2.u^2*sig2.v + 2*rho^2*sig2.e*sig2.u^2*sig2.v + b^4*sig2.e^2*sig2.v^2 + 2*b^2*sig2.e*sig2.u*sig2.v^2 + 2*b^2*rho^2*sig2.e*sig2.u*sig2.v^2 + sig2.u^2*sig2.v^2 - 2*rho^2*sig2.u^2*sig2.v^2 + rho^4*sig2.u^2*sig2.v^2))/(2*(b^2*sig2.e + rho^2*sig2.u + b^2*sig2.v))
  k1 <- (b*rho*(sig2 - sig2.e)*sig2.v)/(rho^2*sig2.u*(sig2 + sig2.v) - b^2*(-sig2 + sig2.e)*(sig2.e + sig2.v))
  k2 <- (b^2*sig2.e*(-sig2 + sig2.e) - rho^2*sig2*sig2.u)/(-(rho^2*sig2.u*(sig2 + sig2.v)) + b^2*(-sig2 + sig2.e)*(sig2.e + sig2.v))
  k <- k2 + (b/rho)*k1
  
  # transition equation: alpha_t = dt + Tt alpha_{t-1} + Ht eta_t
  # drift
  dt <- matrix((1-rho) * theta.bar, nrow = N+1)
  dt <- rbind(dt, 0)
  dt <- rbind(dt, 0)
  
  # transition matrix
  Tt <- (1-k) * rho * diag(N)
  Tt <- cbind( rep(k * rho, N), Tt)
  Tt <- rbind(c(rho, rep(0,N)), Tt)
  Tt <- rbind(Tt, rep(0,N+1))
  Tt <- rbind(Tt, rep(0,N+1))
  Tt <- cbind(Tt, rep(0,N+3))
  Tt <- cbind(Tt, rep(0,N+3))
  
  # shock variance matrix HHt = E(Ht eta eta' Ht') = Ht Ht'
  Ht <- k2* sqrt(sig2.v) * diag(N)
  Ht <- cbind( rep( k1 * sqrt(sig2.u), N), Ht)
  Ht <- cbind( rep( k2*sqrt(sig2.e), N), Ht)
  Ht <- rbind(c(sqrt(sig2.e), rep(0,N+1)), Ht)
  Ht <- rbind(Ht, c(sqrt(sig2.e),0,rep(0,N)))
  Ht <- rbind(Ht, c(0,sqrt(sig2.u),rep(0,N)))
  
  HHt <- Ht %*% t(Ht)
  
  # observation equation: y_t = ct + Zt alpha_t + Gt epsilon_t
  # drift
  ct <- matrix(0, nrow = N)
  ct <- rbind(ct,a - b*(1-rho)*theta.bar /rho)
  
  # observation matrix
  Zt <- cbind(rep(0, N), diag(N))
  Zt <- cbind(Zt, rep(0,N))
  Zt <- cbind(Zt, rep(0,N))
  Zt <- rbind(Zt, c(b/rho,rep(0,N),-b/rho,1))
  
  # shock variance matrix GGt = E(Gt e e' Gt') = Gt Gt'
  Gt = diag(sqrt(sig2.z), N)
  Gt = rbind(Gt,rep(0,N))
  GGt <- Gt %*% t(Gt)
  #GGt <- diag(0,N+1)
  
  return(list(dt = dt, Tt = Tt, HHt = HHt, ct = ct, Zt = Zt, GGt = GGt))
} 

# objective function passed to 'optim'
objective <- function(paras, yt) {
  N <- dim(yt)[1] - 1
  sp <- SP.model(paras[1], paras[2], paras[3], paras[4], paras[5], paras[6], paras[7], N)
  a0 <- c(rep(0.15,N+1),0,0)
  P0 <- diag(0.001, N+3)
  ans <- fkf(a0 = a0, P0 = P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
  return(-ans$logLik)
}

# run estimation for K simulated datasets

results <- list()

for (s in 1:K){
  
  set.seed(s)
  print(s)
  
  # simulate data
  
  # AR(1) process for state
  #theta <- arima.sim(model = list(ar = rho), n = T, innov = rnorm(T) * sqrt(sig2.e), n.start=1, start.innov=theta.hat.0)
  eps = sqrt(sig2.e) * rnorm(T+1)
  theta <- theta.bar 
  aux <- theta
  for (t in 2:(T+1)){
    aux <- (1-rho)*theta.bar + rho*aux + eps[t]
    theta <- c(theta,aux)
  }
  
  # draw initial beliefs mu0
  mu0 <- theta.bar + sqrt(sig2) * rnorm(N)
  
  # create time series of beliefs for N participants
  mu <- mu0
  pub.shocks <- c()
  
  # recursion for submitters beliefs
  # mu_{i,t} = (1-rho) theta.bar + (1-k) rho mu_{i,t-1} + k rho theta_{t-1} + k2 v_{i,t} + k1 u_t + k2 e_t  
  
  aux <- mu0
  for (t in 2:(T+1)){
    shock <- sqrt(sig2.u) * rnorm(1)
    aux <- (1-rho) * theta.bar + (1-k) * rho * aux + k * rho * theta[t-1] + k2 * sqrt(sig2.v) * rnorm(N) + k1 * shock + k2 * eps[t]
    mu <- rbind(mu,aux)
    pub.shocks <- c(pub.shocks,shock)
  }
  
  # add measurement error
  for (t in 1:(T+1)){
    mu[t,] = mu[t,] + sqrt(sig2.z) * rnorm(N)
  }
  
  sim.beliefs <- t(mu[2:(T+1),])
  colnames(sim.beliefs) <- 1:T
  signal <- a + b * theta[1:T] + pub.shocks
  
  data.est <- rbind(sim.beliefs, signal)
  
  # estimation
  
  fit <- optim(paras0, objective, yt = data.est, method="BFGS", hessian = TRUE, control = list(maxit = 10000))
  
  results[[s]] <- fit
  #estimates <- as.numeric(fit$par)
  #estimates <- c(estimates[1]^2/(1 + estimates[1]^2), estimates[2], exp(estimates[3:5]), estimates[6:7], fit$convergence)
  
  #results <- rbind(results, estimates)
  #print(rbind(c(paras.sim[3:9], 0), colMeans(results)))
  
  #points(log(estimates[3]), log(estimates[4]))
  #plot(data.frame(results))
  
  
}

save(results, paras.sim, file="results_010717.Rda")
