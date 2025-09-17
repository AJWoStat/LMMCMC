custom.prior <- function(family = "Custom", log.complexity.prior) {
  structure(list(family = "MatryoshkaDoll", hyper.parameters = log.complexity.prior), class = "prior")
}
