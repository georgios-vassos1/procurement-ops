## Configuration of exogenous effects
Xb2_generator <- function(env) {
  crng <- c(20,500) # cost range
  # Non-linear features:
  #   e01: number of items in the requisition
  #   e02: cost of requisition
  env$e01=I(-seq(5)^2)
  env$e02=c(rep(0,crng[1]),-5.0*seq(crng[1],crng[2])/(crng[2]),rep(-5.0,10L*crng[2]))^3
}

## Sample requisition
rA1 <- function(conf, ...) {
  # Nonlinear features for the decision-to-stop model
  e01      <- conf[["e01"]]
  e02      <- conf[["e02"]]
  # Parameter vector for the endogenous part of the decision-to-stop model
  beta0    <- conf[["beta0"]]
  # Supplier-independent cost estimate of each product
  cost.est <- conf[["cost.est"]]
  # Initialization
  n   <- conf[["nprod"]]
  tau <- conf[["tau"]]
  A   <- seq(n)
  As  <- A
  s   <- k <- c()
  a0  <- 1
  u1  <- NA
  # Execution
  while (a0 && length(s) < n) {
    if (length(s)==n-1L) {
      a <- As
    } else {
      u1  <- conf[["u1.t"]][As] # + endogenous part
      p1  <- (exp(u1) / sum(exp(u1)))
      a   <- sample(As, 1, prob = p1)
    }
    s        <- c(s,a)
    lambda.t <- conf[["lambda.t"]][a] # + endogenous part
    k        <- c(k, max(1L, rpois(1L, lambda = lambda.t)))
    As       <- A[-s]
    cdx      <- round(sum(cost.est[s]*k))
    u0       <- plogis(c(conf[["u0.t"]], e01[length(s)], e02[cdx]) %*% beta0)
    a0       <- rbinom(1, 1, u0)
  }
  list(
    items      = s,
    quantities = k
  )
}

demand_test <- function(env) {

  env$tau <- 365L
  # Number of individuals
  env$n     <- 2L
  env$nprod <- 5L

  demand <- Demand()
  demand$event_generator <- function(lambda.t) rexp(1L, 1.0 / lambda.t)
  demand$requisition_genenerator <- rA1

  ## Configuration
  Xb2_generator(env)
  config.1 <- list(
    tau     =tau,
    nprod   =nprod,
    cost.est=c(100, 50, 50, 10, 10),
    e01     =e01,
    e02     =e02,
    beta0   =c(1.0, 1.0, 5.0),
    u0.t    =0.0,
    u1.t    =0.0,
    lambda.t=0.0
  )

  i <- 1L  # Individual i
  event_times <- c()
  # One lambda per ship
  lambda      <- c(5L, 15L)
  t <- 0.9 # Starting time
  reqsize <- rep(0L, nprod)
  total_q <- rep(0L, nprod)
  while (1L) {
    t <- t + demand$event_generator(lambda[i])
    if (t > 365L) {
      break
    }
    event_times <- c(event_times, t)
    config.1[["u0.t"]]     <- rnorm(1L, 5.0) # Baseline we prefer more items
    config.1[["u1.t"]]     <- rnorm(nprod)
    config.1[["lambda.t"]] <- c(0.1, 0.5, 0.5, 2.0, 2.0)
    reqq <- demand$requisition_genenerator(config.1)
    reqsize[length(reqq$items)] <- reqsize[length(reqq$items)] + 1L
    total_q[reqq$items] <- total_q[reqq$items] + reqq$quantities
  }

  list(reqsize=reqsize, total_q=total_q)
}

