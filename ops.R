
## Discrete-event simulation
initialization <- function(env, B, ...) {
  args <- list(...)
  # Auxiliary effects generation
  Xb2_generator(env)
  env$unresolved_requests <- list()
  # Counter: current number of requisitions for each individual
  env$reqno <- 0L
  env$hdx   <- 0L  # Current index in historical data container
  env$nlit  <- c() # Number of line items in requisitions
  env$historical_data <- data.frame(
    event_time = rep(NA,1e4), ship_id  = NA, req_no = NA, item_no   = NA, quantity = NA,
    supplier   = NA,          order_on = NA, cost   = NA, lead_time = NA, quality = NA)
  env$product  <- NULL
  env$doptim   <- NULL
  env$dobs     <- NULL
  env$dhat     <- NULL
  env$u.obs    <- NULL
  env$u.hat    <- NULL
  env$u.opt    <- NULL
  env$tdx      <- NULL
  env$y        <- matrix(NA, nrow = env$nprod*env$nsupp*env$tau, ncol = env$n.y)
  env$nylags   <- ifelse(!is.null(args[["nylags"]]), args[["nylags"]], 2L)
  env$ylags    <- vector(mode = "numeric", length = env$nprod*env$nsupp*env$nylags*env$n.y)
  env$nalloclags <- ifelse(!is.null(args[["nalloclags"]]), args[["nalloclags"]], 90L)
  env$allocbuffr <- vector(mode = "integer", length = env$nprod*env$nsupp*env$nalloclags)
  env$allocation <- vector(mode = "integer", length = env$nprod*env$nsupp)
  env$curr     <- rep(1L, env$n)
  # Context variables
  env$clen     <- 10L      # Number of context variables
  env$rmcols   <- c(2L:9L) # Introduce unobserved context by removing context variables
  env$p        <- env$clen-length(env$rmcols) # Observed context indices
  # Thompson sampling
  env$mu <- matrix(0.0001, nrow = (env$nprod*env$nsupp*env$p), ncol = env$n.y)
  env$S  <- diag(env$nprod*env$nsupp*env$p) * 70.0
  # Start at the outset of day 1.0
  env$sim_time <- 0.0
  env$next_req_time <- vector(mode = "list", length = env$n)
  for (i in seq(env$n)) {
    env$next_req_time[[i]] <- B$event_generator(args[["lambda.t"]][i])
  }
  # Store the next ship id
  env$next_ship_id <- which.min(plyr::laply(env$next_req_time, function(x) { ifelse(!is.null(x), min(x), Inf) }))
  # Initialize the event list
  if (is.null(x <- env$next_req_time[[env$next_ship_id]])) {
    env$next_event <- c(Inf, 2.0, env$tau, Inf)
  } else {
    env$next_event <- c(min(x), 2.0, env$tau, Inf)
  }
}

timing <- function(env) {
  idx                 <- which.min(env$next_event)
  env$next_event_type <- idx
  env$sim_time        <- env$next_event[idx]
}

requisition_generation <- function(env, B, ...) {
  args <- list(...)
  env$reqno <- env$reqno + 1L
  # Ship index
  i <- env$next_ship_id
  # Simulation time
  t <- env$sim_time
  # Sample a requisition
  requisition <- B$requisition_genenerator(args[["config"]])
  line_items  <- requisition$items
  n_items     <- length(line_items)
  # Requisition index for ship
  env$nlit    <- c(env$nlit, n_items)
  j <- env$curr[i]
  store_data(env, i, j, requisition)
  # Update historical counters and data records
  env$curr[i] <- j+1L
  # Update list of unresolved requests
  env$unresolved_requests[[as.character(env$reqno)]] <- list(
    items=line_items, quantities=requisition$quantities,
    status=rep(FALSE,n_items), gen_at=t)
  # Update the event list
  n <- env$n
  env$next_req_time[[i]] <- sort(c(env$next_req_time[[i]], env$sim_time + B$event_generator(args[["lambda.t"]][i])))
  if (j == length(env$next_req_time[[i]])) {
    env$next_event[1] <- Inf
  } else {
    e.t <- apply(t(seq(n)), 2L, function(x) {
      L <- env$next_req_time[[x]]
      ifelse(!is.null(L), ifelse(is.na(L[env$curr[x]]), Inf, (L[env$curr[x]])), Inf) }
    )
    # Get next ship index
    i <- which.min(e.t)
    env$next_ship_id <- i
    # Update the event list
    env$next_event[1] <- e.t[i]
  }
}

decision_point <- function(env, B, D, Y, bandit=FALSE, ...) {
  args    <- list(...)
  betas   <- args[["betas"]]
  if (is.null(rls.lambda <- args[["rls.lambda"]])) {
    rls.lambda <- 1.0
  }
  tau     <- env$tau
  n       <- env$n
  npd     <- env$nprod
  nsp     <- env$nsupp
  n.y     <- env$n.y
  nylags  <- env$nylags
  nalloclags <- env$nalloclags
  t       <- as.integer(env$sim_time)
  ctr.t   <- rep(0L, npd*nsp)
  Jxt     <- rep(FALSE, npd*nsp) # Visited combinations
  ## Store the context for the true reward model
  rmcols  <- env$rmcols
  # rmdx    <- c(outer(rmcols, seq(npd*nsp), '*')) # Exclude context variables
  # clen    <- ncol(betas)
  clen    <- env$clen
  stopifnot(clen == ncol(betas))
  env$context <- vector(mode = "numeric", length = npd*nsp*clen)
  ## Thompson sampling
  y.0     <- matrix(NA, nrow = nsp, ncol = n.y)
  p       <- env$p
  for (ydx in seq(npd*nsp)-1L) {
    # Simulate outcome for all combinations
    combo <- get.idx(ydx,npd,nsp)
    i   <- combo[1L]
    j   <- combo[2L]
    adx <- c(outer(seq(nalloclags), ((seq(npd)-1)*nsp + (j-1))*nalloclags, '+')) # allocation buffer access
    ylx <- ydx*nylags*n.y ## lag access
    # print(all(ydx == ((i-1L)*nsp+(j-1L))))
    env$context[ydx*clen+seq(clen)] <- c(
      1.0,
      env$ylags[ylx+seq(nylags*n.y)],
      sum(env$allocbuffr[adx]),
      sqrt(sum(env$allocation[(j-1L) + seq(1L, npd*nsp, nsp)]) / npd),
      sin((2.0*pi*t/365.0)+(pi/6.0)))
  }
  # Copy parameters to avoid cheating by updating with every line item
  # instead of updating once at the end of the day (minor impact)
  mu.t <- env$mu
  S.t  <- env$S
  # Loop the list of unresorved requisitions
  for (l in names(env$unresolved_requests)) {
    j_   <- as.integer(l)-1L
    idx  <- sum(env$nlit[1:j_]) * (j_>0) * 1.0
    u0.t <- D$place_an_order() # Propensity-to-order model
    for (m in which(!env$unresolved_requests[[l]]$status)) {
      if (rbinom(1L, 1L, prob=u0.t)) {
        i    <- env$unresolved_requests[[l]]$items[m]
        # u1.t <- D$supplier_selection(nsupp=nsp) # Utility model

        ydx <- ((i-1L)*nsp+(seq(nsp)-1L))
        u1.t   <- vector(mode = "numeric", length = nsp)
        u.true <- vector(mode = "numeric", length = nsp)
        u.r    <- vector(mode = "numeric", length = nsp)
        # Sample parameter vector
        idu    <- (i-1L)*nsp*p+seq(nsp*p)
        beta.r <- matrix(mvtnorm::rmvnorm(1, mean = c(mu.t[idu,]), (diag(1.0, n.y) %x% S.t[idu,idu])), ncol = n.y) # Slow line
        for (k in seq(nsp)) {
          u1.t[k] <- - c(D$supplier_utility(
            n.y  =n.y,
            beta =betas[ydx[k]*n.y+seq(n.y),-rmcols,drop=FALSE],
            x    =env$context[ydx[k]*clen+seq(clen)][-rmcols],
            Sigma=diag(0L, n.y),
            k    =k
          ) %*% env$w)
          y.0[k,] <- Y$generation(
            n.y  =n.y,
            beta =betas[ydx[k]*n.y+seq(n.y),],
            x    =env$context[ydx[k]*clen+seq(clen)],
            Sigma=diag(0L, n.y)
          )
          u.true[k] <- - c(y.0[k,] %*% env$w)
          u.r[k] <- - c(D$bandit_utility(
            n.y  =n.y,
            beta =t(beta.r[(k-1L)*p+seq(p),,drop=FALSE]),
            x    =env$context[ydx[k]*clen+seq(clen)][-rmcols],
            Sigma=diag(0L, n.y),
            k    =k
          ) %*% env$w)
        }
        # Supplier selection
        if (t > (tau %/% 2L)) {
          if (bandit) {
            # Bandit takes decision
            d.jk  <- which.max(u.r)
          } else {
            # Stochastic policy
            p1.t  <- softmax(u1.t)
            d.jk  <- sample(seq(nsp), 1L, prob=p1.t)
          }
          # Compute optimal decision
          d.opt <- which.max(u.true)
          # Store decisions
          env$tdx     <- c(env$tdx, t)
          env$product <- c(env$product, i)
          env$doptim  <- c(env$doptim, d.opt)
          env$dobs    <- c(env$dobs, d.jk)
          env$u.obs   <- c(env$u.obs, u.true[d.jk])
          env$u.opt   <- c(env$u.opt, u.true[d.opt])
          # Thompson sampling
          d.hat       <- which.max(u.r)
          env$dhat    <- c(env$dhat, d.hat)
          env$u.hat   <- c(env$u.hat, u.true[d.hat])
        } else {
          # Random policy
          d.jk <- sample(seq(nsp), 1L)
        }
        # Update counters
        ctx        <- (i-1)*nsp + d.jk
        ctr.t[ctx] <- ctr.t[ctx] + env$unresolved_requests[[l]]$quantities[m]
        Jxt[ctx]   <- TRUE
        adx <- (i-1)*nsp + d.jk
        env$allocation[adx] <- env$allocation[adx] + env$unresolved_requests[[l]]$quantities[m]
        # Update records
        env$historical_data[idx+m,6:7] <- c(d.jk, t)
        env$unresolved_requests[[l]]$status[m] <- TRUE
        ## RLS update
        id.b <- ((i-1)*nsp+(d.jk-1))*p + seq(p) # supp.bx(i,d.jk,n.supp,p) 
        # rmdx <- c(outer(rmcols, seq(0L, clen*(nsp-1L), clen), '+')) # supp.bx(i,d.jk,n.supp,N,1)+t-1
        # X    <- t(env$context[(i-1L)*nsp*clen+seq(nsp*clen)][-rmdx])
        X    <- env$context[ydx[d.jk]*clen+seq(clen)][-rmcols]
        w.t  <- 1.0 # IPW
        env$S[id.b,id.b] <- (env$S[id.b,id.b]-(w.t*env$S[id.b,id.b]%*%X%*%X%*%env$S[id.b,id.b])/
                               as.numeric(1.0+w.t*X%*%env$S[id.b,id.b]%*%X))/rls.lambda
        eps              <- y.0[d.jk,,drop=FALSE]-X%*%env$mu[id.b,]
        env$mu[id.b,]    <- env$mu[id.b,]+(w.t)*env$S[id.b,id.b]%*%X%*%as.numeric(eps)
      }
    }
    if (all(env$unresolved_requests[[l]]$status)) {
      env$unresolved_requests[[l]] <- NULL
    }
  }
  ## Compute the outcome variable at time t
  #  At a single time point t, a given supplier will not provide two different prices for the same product,
  #  given that the same product can be part of two different requisitions.
  for (ydx in seq(npd*nsp)-1L) {
    # Simulate outcome for all combinations
    combo <- get.idx(ydx,npd,nsp)
    i   <- combo[1L]
    j   <- combo[2L]
    adx <- c(outer(seq(nalloclags), ((seq(npd)-1)*nsp + (j-1))*nalloclags, '+')) # allocation buffer access
    ylx <- ydx*nylags*n.y ## lag access
    # y.t <- Y$generation(n.y=n.y)
    # print(all(ydx == ((i-1L)*nsp+(j-1L))))
    y.t <- Y$generation(
      beta =betas[ydx*n.y+seq(n.y),],
      x    =env$context[ydx*clen+seq(clen)],
      Sigma=args[["Sigma"]]
    )
    env$y[((i-1)*nsp+(j-1))*tau+t,] <- y.t
    # Observed combinations
    if (Jxt[ydx+1L]) {
      # Update data containers
      env$ylags[ylx+seq(nylags*n.y)] <- c(y.t, env$ylags[ylx+seq((nylags-1)*n.y)])
    }
    ## Update the allocation buffer every day
    # env$allocbuffr[adx] <- rdeque(env$allocbuffr[adx], ctr.t[ydx+1])
    env$allocbuffr[adx] <- c(env$allocbuffr[adx][-1L], ctr.t[ydx+1])
  }
  # Schedule next event times
  e.t <- apply(t(seq(n)), 2L, function(x) {
    L <- env$next_req_time[[x]]
    ifelse(!is.null(L), ifelse(is.na(L[env$curr[x]]), Inf, (L[env$curr[x]])), Inf) }
  )
  # Get next ship index
  i <- which.min(e.t)
  env$next_ship_id <- i
  # Update the event list
  env$next_event[1] <- e.t[i]
  env$next_event[2] <- env$next_event[2] + 1.0
}

supplier_evaluation <- function(env) {
  NULL
}

store_data <- function(env, i, j, requisition) {
  nreq <- env$reqno
  idx  <- env$hdx+seq(env$nlit[nreq])
  env$historical_data[idx,1:5] <- cbind(
    env$next_req_time[[i]][j], i, nreq, requisition$items, requisition$quantities)
  env$hdx <- env$hdx + env$nlit[nreq]
}

utility.func <- function(...) {
  args  <- list(...)
  beta  <- args[["beta"]]
  x     <- args[["x"]]
  Sigma <- args[["Sigma"]]
  c(mvtnorm::rmvnorm(1L, c(beta %*% x), Sigma))
}

run_sim <- function(supplier_utility, bandit, ...) {
  args <- list(...)

  sim <- Simulation(tau = 2L*365L, nprod = 3L, n = 50L)
  ## Demand configuration
  B   <- Demand()
  D   <- Decision()
  Y   <- Outcome()

  B$event_generator         <- function(lambda.t) rexp(1L, 1.0 / lambda.t)
  B$requisition_genenerator <- rA1

  D$place_an_order   <- function(...) runif(1L, 0.1, 0.9)
  D$bandit_utility   <- utility.func
  D$supplier_utility <- supplier_utility

  Y$generation <- function(...) {
    args  <- list(...)
    beta  <- args[["beta"]]
    x     <- args[["x"]]
    Sigma <- args[["Sigma"]]
    c(mvtnorm::rmvnorm(1L, c(beta %*% x), Sigma))
  }

  betas  <- args[["betas"]]
  lambda <- args[["lambda"]]
  config <- args[["config"]]

  initialization(sim, B, lambda.t=lambda)
  while (1) {
    # Determine the next event
    timing(sim)
    # Invoke the appropriate event routine
    j <- sim$next_event_type
    switch (j,
      requisition_generation(sim, B, lambda.t=lambda, config=config),
      decision_point(sim, B, D, Y, bandit, betas=betas, Sigma=diag(10L,sim$n.y), rls.lambda=0.98),
      break,
      supplier_evaluation(sim)
    )
  }

  cost.opt <- sum(-sim$u.opt)
  cost.obs <- sum(-sim$u.obs)
  cost.hat <- sum(-sim$u.hat)

  data.frame(t=sim$tdx-365L, r.obs=sim$u.opt-sim$u.obs, r.hat=sim$u.opt-sim$u.hat) |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      r.obs=sum(r.obs),
      r.hat=sum(r.hat)
    ) |>
    dplyr::ungroup() |>
    as.data.frame() -> df

  list(c(cost.opt,cost.obs,cost.hat), df)
}

