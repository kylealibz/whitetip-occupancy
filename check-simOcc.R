# Fit Doser et al.'s single-species model to a small simOcc() example.
# The main study constructs habitat-dependent visibility explicitly.
if (dir.exists(".r-library")) .libPaths(c(".r-library", .libPaths()))
source("R/models.R")
set.seed(753)
sim <- spOccupancy::simOcc(J.x = 15, J.y = 10, n.rep = rep(4L, 150),
                           beta = c(-0.3, 0.9), alpha = c(-0.9, 0.8))
stopifnot(all(dim(sim$y) == c(150L, 4L)), all(sim$y[sim$z == 0, ] == 0))
fits <- evaluate_models(sim$X[, 2], sim$y, sim$X.p[, , 2], "pgocc")
stopifnot(all(is.finite(fits$estimate)), is.finite(fits$rhat[1]))
print(data.frame(truth = 0.9, fits), row.names = FALSE)
