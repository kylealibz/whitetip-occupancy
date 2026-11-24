report_results <- function(results, out) {
  groups <- split(results, list(results$scenario, results$deployments, results$model))
  metrics <- do.call(rbind, lapply(groups, function(group) {
    ok <- is.finite(group$estimate) & is.finite(group$lower) & is.finite(group$upper)
    error <- group$estimate[ok] - group$truth[ok]
    covered <- group$lower[ok] <= group$truth[ok] & group$upper[ok] >= group$truth[ok]
    coverage <- if (any(ok)) mean(covered) else NA_real_
    data.frame(scenario = group$scenario[1], deployments = group$deployments[1],
               model = group$model[1], attempted = nrow(group), usable = sum(ok),
               warnings = sum(nzchar(group$diagnostic) & ok),
               rhat_over_1_05 = if (any(is.finite(group$rhat))) {
                 sum(group$rhat > 1.05, na.rm = TRUE)
               } else NA_integer_,
               bias = mean(error), rmse = sqrt(mean(error^2)), coverage = coverage,
               coverage_mcse = if (any(ok)) sqrt(coverage * (1 - coverage) / sum(ok)) else NA_real_)
  }))
  write.csv(metrics, file.path(out, "summary.csv"), row.names = FALSE)

  png(file.path(out, "habitat-recovery.png"), width = 1500, height = 730, res = 135)
  on.exit(dev.off())
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  limits <- range(c(results$estimate, results$truth), finite = TRUE)
  for (scenario in unique(results$scenario)) {
    subset <- results[results$scenario == scenario, ]
    levels <- expand.grid(model = c("Occupancy", "Ever detected"),
                          deployments = sort(unique(subset$deployments)))
    groups <- lapply(seq_len(nrow(levels)), function(i) {
      values <- subset$estimate[subset$model == levels$model[i] &
                                subset$deployments == levels$deployments[i]]
      values[is.finite(values)]
    })
    boxplot(groups, at = c(1, 2, 4, 5, 7, 8), xaxt = "n", ylim = limits,
            col = c("#80BCD0", "#E9AD6C"), main = scenario,
            xlab = "Deployments per reef", ylab = "Estimated habitat coefficient")
    axis(1, at = c(1.5, 4.5, 7.5), labels = sort(unique(subset$deployments)))
    abline(h = subset$truth[1], lty = 2, lwd = 2)
    legend("bottomright", c("Occupancy", "Ever detected", "True habitat effect"),
           col = c("#80BCD0", "#E9AD6C", "black"), lwd = c(5, 5, 2),
           lty = c(1, 1, 2), bty = "n", cex = 0.8)
  }
  invisible(metrics)
}
