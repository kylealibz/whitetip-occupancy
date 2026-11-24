# Run from the project directory: Rscript run.R [mle|pgocc] [replicates] [workers]
args <- commandArgs(trailingOnly = TRUE)
engine <- if (length(args)) args[1] else "pgocc"
stopifnot(engine %in% c("mle", "pgocc"))
reps <- if (length(args) >= 2L) as.integer(args[2]) else 100L
stopifnot(!is.na(reps), reps >= 1L)
if (dir.exists(".r-library")) .libPaths(c(".r-library", .libPaths()))
if (engine == "pgocc" && !requireNamespace("spOccupancy", quietly = TRUE)) {
  stop("Install spOccupancy with install.packages('spOccupancy') before using pgocc.")
}
source("R/models.R")
dir.create(file.path("results", engine), recursive = TRUE, showWarnings = FALSE)
out <- file.path("results", engine)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
set.seed(29117)
seeds <- matrix(sample.int(.Machine$integer.max, 7L * reps), nrow = reps)
efforts <- c(2L, 4L, 6L)
truth <- 0.9
scenarios <- c(Independent = 0, Confounded = -0.8)
run_replicate <- function(replicate) {
  rows <- vector("list", length(scenarios) * length(efforts))
  index <- 0L
  set.seed(seeds[replicate, 1])
  latent <- make_reefs(150L, max(efforts), truth)
  for (scenario in names(scenarios)) {
    survey <- observe_reefs(latent, scenarios[[scenario]])
    for (k in seq_along(efforts)) {
      index <- index + 1L
      days <- seq_len(efforts[k])
      set.seed(seeds[replicate, 1L + k + 3L * (scenario == "Confounded")])
      fits <- evaluate_models(latent$habitat, survey$y[, days],
                              survey$visibility[, days], engine)
      rows[[index]] <- data.frame(replicate, scenario, deployments = efforts[k],
                                  truth, fits)
    }
  }
  do.call(rbind, rows)
}

workers <- if (length(args) >= 3L) as.integer(args[3]) else 1L
stopifnot(!is.na(workers), workers >= 1L)
if (workers > 1L) {
  cluster <- parallel::makeCluster(min(reps, workers))
  results <- tryCatch({
    parallel::clusterCall(cluster, function(paths) {
      .libPaths(paths)
      source("R/models.R", local = .GlobalEnv)
      NULL
    }, .libPaths())
    parallel::clusterExport(cluster, c("engine", "seeds", "efforts", "truth",
                                      "scenarios", "run_replicate"))
    parallel::parLapplyLB(cluster, seq_len(reps), run_replicate)
  }, finally = parallel::stopCluster(cluster))
} else {
  results <- lapply(seq_len(reps), function(i) {
    result <- run_replicate(i)
    if (i %% 10L == 0L || i == reps) message(i, "/", reps)
    result
  })
}
results <- do.call(rbind, results)
if (!any(is.finite(results$estimate[results$model == "Occupancy"]))) {
  stop("All occupancy fits failed: ", results$diagnostic[results$model == "Occupancy"][1])
}
write.csv(results, file.path(out, "replicates.csv"), row.names = FALSE)
source("R/report.R")
report_results(results, out)



