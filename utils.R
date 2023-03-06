
## Post hoc utilities
postproc <- function(env, ycols=8:10) {
  env$historical_data <- env$historical_data[seq(env$hdx),]
  I <- env$historical_data$item_no
  J <- env$historical_data$supplier
  idx <- ((I-1L)*env$nsupp+(J-1))*env$tau+env$historical_data$order_on
  env$historical_data[,ycols[1]] <- env$y[idx,1]
  env$historical_data[,ycols[2]] <- env$y[idx,2]
  env$historical_data[,ycols[3]] <- env$y[idx,3]
  ## Remove rows where outcome is NA
  if (length(rem <- which(is.na(env$historical_data$supplier)))) {
    env$historical_data <- env$historical_data[-rem,]
  }
}

# Average cost computation
avg_cost <- function(env, burnin=365L) {
  idx <- historical_data$order_on > burnin
  sum(historical_data$quantity[idx]*(as.matrix(historical_data[idx,8:10]) %*% env$w), na.rm = TRUE) / 
    sum(historical_data$quantity[idx])
}

## General utilities
# Softmax function
softmax <- function(u) exp(u)/sum(exp(u))

# Convert 2D index to 1D index
get.idx <- function(x,N,M) {
  if (x>N*M) return(-1)
  n <- c(min(N,M),max(N,M))[(N<=M)+1L]
  c(x%/%n,x%%n)+1L
}

## Modified Lewis algorithm
rEventTimes <- function(tau, lambda.dom, threshold) {
  te <- vector(mode = "numeric", length = 1e4)
  t <- j <- 1L
  while(1) {
    v <- runif(1L)
    t <- t - (log(v) / lambda.dom[floor(t)])
    if (t > tau) break
    if (runif(1L) < threshold[floor(t)]) {
      te[j] <- t
      j <- j + 1L
    }
  }
  te[1L:(j-1L)]
}

rEventTimes.test <- function(env) {
  te <- vector(mode = "list", length = env$n)
  for (i in seq(env$n)) {
    idx     <- block.tau(env, i)
    te[[i]] <- rEventTimes(env$tau, env$lambda.dom[idx], env$threshold[idx])
  }
  te
}

## Simulating the next event times over a rollout period of 3 days in the future
## This module might depend on the configuration of the intensity function.
## The idea that each requisition event schedules the next will not work nicely
## because every decision day the information changes. Therefore, it is better to
## schedule the arrivals over a rollout period in the future on a decision day.
rNextEventTime <- function(currday, lambda.dom.t, threshold.t, rollout=3.0) {
  te <- c()
  t  <- currday
  while(1) {
    v <- runif(1L)
    t <- t - (log(v) / lambda.dom.t)
    if (t > currday+rollout-1e-6) break
    if (runif(1L) < threshold.t) te <- c(te,t)
  }
  te
}

## Generate the event history process with rNextEventTime to compare with rEventTimes
rNextEventTime.test <- function(env) {
  n   <- env$n
  tau <- env$tau
  te <- vector(mode = "list", length = n)
  idx <- apply(t(seq(n)), 2, function(i) block.tau(env, i))
  for (t in seq(tau)) {
    for (i in seq(n)) {
      if (!is.null(e.t <- rNextEventTime(t, env$lambda[idx[t,i]], env$threshold[idx[t,i]]))) {
        te[[i]] <- c(te[[i]], e.t[e.t<=tau])
      }
    }
  }
  te
}

## Weibull function
weibull <- function(t, shape, scale) shape * scale^(-shape) * (t / scale)^(shape-1)

# Get the standard coefficients of a harmonic
harmonic.coef <- function(env) {
  phi <- env$phi
  n   <- length(phi)
  c(sin(phi), cos(phi))[c(t(matrix(c(c(1:n),n+c(1:n)), ncol = 2)))]
}

# Generate harmonic feature with period Ts over vectors (1,...,t) and phi=(phi.1,...,phi.n)
harmonics <- function(env, ...) {
  args <- list(...)
  if (is.null(Ts <- args$Ts)) {
    Ts <- env$Ts
  }
  if (is.null(t <- args$t)) {
    t <- seq(env$tau)
  }
  if (is.null(phi <- args$phi)) {
    phi <- env$phi
  }
  # Ts  <- env$Ts
  m   <- length(t)
  n   <- length(phi)
  cbind(
    crossprod(rbind(cos(2*pi*t/Ts)), rbind(sin(phi))),
    crossprod(rbind(sin(2*pi*t/Ts)), rbind(cos(phi)))
  )[,c(t(matrix(c(c(1:n),n+c(1:n)), ncol = 2)))]
}

# Smooth shock feature generator
# tau: simulation horizon
# Ts:  period of interest (tau - burnin)
# loc: temporal location of the shock
# spike.tip:  magnitude of the shock
# spike.band: bandwidth of the shock
smooth.shock <- function(env, ...) {
  args <- list(...)
  spike.loc <- env$spike.loc
  spike.tip <- env$spike.tip
  spike.bnd <- env$spike.bnd
  span <- with(env, (spike.loc+c(1:tau))/Ts) + 1e-32
  if (is.null(args$bound)) {
    return(spike.tip*sin(spike.bnd*pi*span)/(spike.bnd*pi*span))
  } else {
    return(apply(t(spike.tip/(spike.bnd*pi*abs(span))), 2, function(x) min(x, spike.tip)))
  }
}

## Unobserved heterogeneity
draw_ind_gamma_rvs <- function(n, shape = 2.0, scale = 1.0/2.0, seed=NULL) {
  if (is.null(seed)) {
    seed <- sample.int(12345,size=n,replace=FALSE)
  }
  zeta <- rep(NA,n)
  for (i in seq_along(seed)) {
    set.seed(seed[i])
    zeta[i] <- rgamma(1, shape = 2.0, scale = 1.0/2.0)
  }
  zeta
}

## Access block of fixed size tau for individual i
block.tau <- function(env, i) {
  (i-1L) * env$tau + seq(env$tau)
}

## Custom style for ggplot
ggstyle <- function(legend.position="top", axis.text.x.angle=0, size=12L) {
  theme_minimal() +
    theme(
      strip.text.x     = element_text(family="serif", size=size, colour="black", angle=0),
      axis.title.x     = element_text(family="serif", size=size),
      axis.title.y     = element_text(family="serif", size=size),
      axis.text.y      = element_text(family="serif", size=size),
      axis.text.x      = element_text(family="serif", angle=axis.text.x.angle, hjust=1, size=size),
      legend.position  = legend.position,
      legend.text      = element_text(family="serif", size=size),
      legend.title     = element_text(family="serif", size=size),
      plot.title       = element_text(size = size+4L, family="serif", hjust = 0, colour = "black", face = "bold"),
      plot.subtitle    = element_text(size = size+2L, family="serif", hjust = 0, face = "italic"),
      strip.background = element_rect(colour = "steelblue", fill = "white"),
      panel.border     = element_rect(linetype = "solid", fill = NA)
    )
}
