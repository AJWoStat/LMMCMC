custom.prior <- function(family = "Custom", log.complexity.prior) {
  structure(list(family = family, hyper.parameters = log.complexity.prior), class = "prior")
}
