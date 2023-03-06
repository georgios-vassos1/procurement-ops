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

dir.create("~/procops/dumps")
policies <- names(policy)

for (idx in policies) {
  joint.res <- foreach (i=seq(1000L), .combine=rbind, .multicombine=TRUE) %dopar% {
    thread.res <- run_sim(policy[[idx]], idx=="bandits", betas=betas, lambda=lambda, config=config) # calling a function
    thread.res # Equivalent to joint.res = cbind(joint.res, thread.res)
  }
  saveRDS(joint.res[,1], paste0("~/procops/dumps/tempcost-", idx))
  saveRDS(joint.res[,2], paste0("~/procops/dumps/regret-", idx))
}

# stop cluster
stopCluster(cl)

## Results (processing output dumps)
library(ggplot2)

process_dump_violin <- function(file, policy, ...) {
  stopifnot(!is.null(policy))

  do.call(rbind, readRDS(file)) |>
    tibble::rownames_to_column() |>
    dplyr::transmute(
      exp=factor(purrr::map_chr(stringr::str_split(rowname, stringr::fixed(".")), 2L)),
      t    =t,
      r.obs=r.obs,
      r.hat=r.hat
    ) |>
    dplyr::group_by(exp) |>
    dplyr::transmute(
      t     =t,
      Actual=cumsum(r.obs),
      Bandit=cumsum(r.hat)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(exp) |>
    dplyr::summarise(
      Actual=dplyr::last(Actual),
      Bandit=dplyr::last(Bandit),
      policy=policy
    ) |>
    dplyr::ungroup()
}

# Data dumps with results from experiments
dumps <- list(
  Supplier.1="~/procops/dumps/regret-supplier.1",
  Supplier.2="~/procops/dumps/regret-supplier.2",
  Random    ="~/procops/dumps/regret-random",
  Utility   ="~/procops/dumps/regret-static",
  Bandit    ="~/procops/dumps/regret-bandits"
)

# Process dumps into a data frame
lapply(names(dumps), function(x) process_dump_violin(dumps[[x]], x)) -> L
do.call(rbind, L) -> df

## Box plot
df |>
  ggplot(aes(x=policy, y=(Actual))) +
  # geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, fill="white") +
  xlab("") + ylab("Total regret") +
  # ggtitle("1000 Monte Carlo simulations") +
  scale_color_grey() +
  scale_fill_grey() +
  ggstyle(legend="None", axis.text.x.angle=45, size=22L)

## Spaghetti plot
do.call(rbind, readRDS("~/procops/dumps/regret-static")) |>
  tibble::rownames_to_column() |>
  dplyr::transmute(
    exp=factor(purrr::map_chr(stringr::str_split(rowname, stringr::fixed(".")), 2L)),
    t    =t,
    r.obs=r.obs,
    r.hat=r.hat
  ) |>
  dplyr::group_by(exp) |>
  dplyr::transmute(
    t         =t,
    Utility   =cumsum(r.obs),
    Bandit    =cumsum(r.hat)
  ) |>
  dplyr::ungroup() |>
  reshape2::melt(id.vars=c("exp","t"), variable.name="model") |>
  ggplot(aes(x=t, y=value, group=interaction(exp, model), color=model)) +
  geom_line(size=0.1, alpha=0.4) +
  geom_smooth(aes(group=model, color=model), colour=rep(gray.colors(2L, 0.1, 0.6), each=80L), se = TRUE) +
  scale_color_grey() +
  # facet_wrap(~regret, scales = "free") +
  # ggtitle("1000 Monte Carlo simulations") +
  xlab("day") + ylab("Cumulative regret") +
  ggstyle(legend="bottom", size=22L) +
  guides(colour = guide_legend(override.aes = list(size=1.0, alpha=1.0)))
