rm(list=ls())

source("~/procops/utils.R")
source("~/procops/interface.R")
source("~/procops/demand.R")
source("~/procops/ops.R")

# Simulation example object needed for initialization of parameters
sim <- Simulation(tau = 2L*365L, nprod = 3L, n = 50L)

# Supplier selection parameters
betas <- matrix(
  c(
    # Product 1 - Supplier 1
    c(75.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, -0.04, 0.5, 10.0),
    c(50.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, -0.00, 0.0, 10.0),
    c(50.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, -0.00, 0.0, 10.0),
    # Product 1 - Supplier 2
    c(50.0, 0.6, 0.0, 0.0, 0.1, 0.0, 0.0, -0.05, 1.0, -10.0),
    c(25.0, 0.0, 0.6, 0.0, 0.0, 0.1, 0.0, -0.01, 0.0, -10.0),
    c(25.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.1, -0.01, 0.0, -10.0),
    # Product 2 - Supplier 1
    c(25.0, 0.7, 0.0, 0.0, 0.15, 0.0, 0.0, -0.07, 1.0, 10.0),
    c(15.0, 0.0, 0.7, 0.0, 0.00, 0.0, 0.0, -0.05, 0.0, 10.0),
    c(15.0, 0.0, 0.0, 0.7, 0.00, 0.0, 0.1, -0.05, 0.0, -0.0),
    # Product 2 - Supplier 2
    c(40.0, 0.5, 0.0, 0.0, 0.2, 0.0, 0.0, -0.01, 0.5, 0.0),
    c(30.0, 0.0, 0.5, 0.0, 0.0, 0.2, 0.0, -0.00, 0.0, 0.0),
    c(30.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.2, -0.00, 0.0, 0.0),
    # Product 3 - Supplier 1
    c(0.1, 0.3, 0.0, 0.0, 0.1, 0.0, 0.0, 0.5, 0.0, 1.0),
    c(0.1, 0.0, 0.3, 0.0, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0),
    c(0.1, 0.0, 0.0, 0.3, 0.0, 0.0, 0.1, 0.5, 0.0, 1.0),
    # Product 3 - Supplier 2
    c(0.1, 0.3, 0.0, 0.0, 0.1, 0.0, 0.0, 0.5, 0.0, 1.0),
    c(0.1, 0.0, 0.3, 0.0, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0),
    c(0.1, 0.0, 0.0, 0.3, 0.0, 0.0, 0.1, 0.5, 0.0, 1.0)),
    # # Product 4 - Supplier 1
    # c(0.1, 0.3, 0.0, 0.0, 0.1, 0.0, 0.0, 0.5, 0.0, 1.0),
    # c(0.1, 0.0, 0.3, 0.0, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0),
    # c(0.1, 0.0, 0.0, 0.3, 0.0, 0.0, 0.1, 0.5, 0.0, 1.0),
    # # Product 4 - Supplier 2
    # c(0.1, 0.3, 0.0, 0.0, 0.1, 0.0, 0.0, 0.5, 0.0, 1.0),
    # c(0.1, 0.0, 0.3, 0.0, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0),
    # c(0.1, 0.0, 0.0, 0.3, 0.0, 0.0, 0.1, 0.5, 0.0, 1.0)),
  nrow = sim$nprod*sim$nsupp*sim$n.y, byrow = TRUE)

# Time-to-request model parameters
lambda <- rep(90.0, sim$n) # One lambda per ship

# Requisitions sampling model parameters
config <- list(
  tau     =sim$tau,
  nprod   =sim$nprod,
  cost.est=c(100, 50, 50, 10, 10),
  e01     =sim$e01,
  e02     =sim$e02,
  beta0   =c(0.1, 1.0, 2.0),
  u0.t=rnorm(1L, 5.0),
  u1.t=rnorm(sim$nprod),
  lambda.t=c(0.1, 0.5, 0.5, 2.0, 2.0)
)

# Supplier selection policies
policy <- list(
  "supplier.1"=function(...) {
    args <- list(...)
    k    <- args[["k"]]
    n.y  <- args[["n.y"]]
    matrix(rep(c(500.0, 0.0)[k], n.y), nrow=1L)
  },
  "supplier.2"=function(...) {
    args <- list(...)
    k    <- args[["k"]]
    n.y  <- args[["n.y"]]
    matrix(rep(c(0.0, 500.0)[k], n.y), nrow=1L)
  },
  "random"    =function(...) {
    args <- list(...)
    n.y <- args[["n.y"]]
    matrix(rep(1.0, n.y), nrow=1L, ncol=n.y)
  },
  "bandits"   =utility.func,
  "static"    =utility.func
)

# Multi-threaded Monte Carlo simulations
library(foreach)
library(doParallel)

# setup parallel backend to use many processors
cores <- detectCores()
cl    <- makeCluster(cores[1L]-4L) # not to overload your computer
registerDoParallel(cl)

joint.res <- foreach (i=seq(10L), .combine=rbind, .multicombine=TRUE) %dopar% {
  thread.res <- run_sim(policy[["static"]], betas=betas, lambda=lambda, config=config) # calling a function
  thread.res # Equivalent to joint.res = cbind(joint.res, thread.res)
}
# stop cluster
stopCluster(cl)

