##' Fitting sets

# define all time points to use
# first pulseR run
allSets <- list(
  c(0, 1, 2, 4, 6, 8, 16)
)

# define subsets of time points to use
# uses first pulseR run with allSets as initial condition
timeSets <- list(
  c(0, 1, 2),
  c(0, 2, 4),
  c(0, 4, 6),
  c(0, 6, 8),
  c(0, 8, 16)
)

# timeSets <- list(
#   c(0, 1, 2, 4, 6),
#   c(0, 2, 4, 6, 8),
#   c(0, 4, 6, 8, 16)
# )
