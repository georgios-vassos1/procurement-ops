
Demand <- function(...) {
  args   <- list(...)
  demand <- new.env(hash = FALSE, parent = emptyenv())
  demand$event_generator <- function(...) {
    return(waiting_time <- NULL)
  }
  demand$requisition_genenerator <- function(...) {
    return(list(items=NULL, quantities=NULL))
  }
  class(demand) <- "Demand"
  return(demand)
}

Decision <- function(...) {
  args   <- list(...)
  decision <- new.env(hash = FALSE, parent = emptyenv())
  decision$place_an_order     <- NULL
  decision$supplier_selection <- NULL
  class(decision) <- "Decision"
  return(decision)
}

Outcome <- function(...) {
  args   <- list(...)
  outcome <- new.env(hash = FALSE, parent = emptyenv())
  outcome$generation <- NULL
  class(outcome) <- "Outcome"
  return(outcome)
}

Simulation <- function(tau=365L, n=2L, nprod=5L, nsupp=2L, n.y=3L, w=c(0.5, 0.25, 0.25), ...) {
  args   <- list(...)
  Simulation <- new.env(hash = FALSE, parent = emptyenv())
  Simulation$tau                    <- tau    # termination of the simulation
  Simulation$n                      <- n      # Number of operational sites
  Simulation$w                      <- w      # Business weight for reward calculation
  Simulation$nprod                  <- nprod  # Number of products
  Simulation$nsupp                  <- nsupp  # Number of suppliers
  Simulation$n.y                    <- n.y    # Number of outcomes
  Simulation$sim_time               <- 0.0    # Simulation clock
  Simulation$next_event_type        <- NA
  Simulation$unresolved_requests    <- list()
  Simulation$next_event             <- rep(NA, 4L) # Event list
  class(Simulation) <- "Simulation"
  return(Simulation)
}

