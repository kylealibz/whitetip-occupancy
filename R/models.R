make_reefs <- function(n, visits, habitat_effect) {
  habitat <- rnorm(n)
  psi <- plogis(-0.3 + habitat_effect * habitat)
  list(habitat = habitat, occupied = rbinom(n, 1, psi),
       visibility_noise = matrix(rnorm(n * visits), n),
       detection_draw = matrix(runif(n * visits), n))
}

observe_reefs <- function(reefs, visibility_slope) {
  visibility <- reefs$visibility_noise + visibility_slope * reefs$habitat
  detection <- plogis(-0.9 + 0.8 * visibility)
  list(visibility = visibility,
       y = 1L * (reefs$occupied == 1L & reefs$detection_draw < detection))
}

# Marginalize over the two possible occupancy states. A site with a
# detection can only be occupied. A zero history allows either state.
fit_mle <- function(habitat, y, visibility) {
  detected <- rowSums(y) > 0
  nll <- function(par) {
    psi <- plogis(par[1] + par[2] * habitat)
    eta <- par[3] + par[4] * visibility
    log_history <- rowSums(ifelse(y == 1, plogis(eta, log.p = TRUE),
                                   plogis(-eta, log.p = TRUE)))
    log_lik <- log(psi) + log_history
    log_lik[!detected] <- log1p(psi[!detected] * expm1(log_history[!detected]))
    -sum(log_lik)
  }
  fit <- optim(rep(0, 4), nll, method = "BFGS", hessian = TRUE,
               control = list(reltol = 1e-10))
  if (fit$convergence != 0L || max(abs(fit$par)) > 15) {
    stop("Optimizer did not find a stable interior solution")
  }
  covariance <- chol2inv(chol(fit$hessian))
  wald(fit$par[2], sqrt(covariance[2, 2]))
}

wald <- function(estimate, se) {
  c(estimate = estimate, lower = estimate - qnorm(0.975) * se,
    upper = estimate + qnorm(0.975) * se, rhat = NA_real_)
}

fit_pgocc <- function(habitat, y, visibility) {
  fit <- spOccupancy::PGOcc(
    occ.formula = ~ habitat, det.formula = ~ visibility,
    data = list(y = y, occ.covs = data.frame(habitat = habitat),
                det.covs = list(visibility = visibility)),
    priors = list(beta.normal = list(mean = 0, var = 6.25),
                  alpha.normal = list(mean = 0, var = 6.25)),
    n.samples = 1200, n.burn = 400, n.chains = 4, verbose = FALSE
  )
  draws <- fit$beta.samples[, 2]
  limits <- quantile(draws, c(0.025, 0.975), names = FALSE)
  c(estimate = mean(draws), lower = limits[1], upper = limits[2],
    rhat = unname(fit$rhat$beta[2]))
}

evaluate_models <- function(habitat, y, visibility, engine) {
  detected <- rowSums(y) > 0
  capture_fit <- function(label, fit) {
    warnings <- character()
    values <- tryCatch(withCallingHandlers(fit(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      warnings <<- c(warnings, conditionMessage(e))
      c(estimate = NA_real_, lower = NA_real_, upper = NA_real_, rhat = NA_real_)
    })
    data.frame(model = label, estimate = values[["estimate"]],
               lower = values[["lower"]], upper = values[["upper"]],
               rhat = values[["rhat"]], diagnostic = paste(warnings, collapse = " | "))
  }
  occupancy <- if (engine == "pgocc") fit_pgocc else fit_mle
  rbind(capture_fit("Occupancy", function() occupancy(habitat, y, visibility)),
        capture_fit("Ever detected", function() {
          fit <- glm(detected ~ habitat, family = binomial())
          if (!fit$converged || max(abs(coef(fit))) > 15) {
            stop("Logistic regression did not find a stable interior solution")
          }
          wald(unname(coef(fit)[2]), sqrt(vcov(fit)[2, 2]))
        }))
}

