md.prior <- function(eta = 1) {
  if(!is.numeric(eta) || eta[1]<=0) eta=1
  structure(list(family = "MatryoshkaDoll", hyper.parameters = eta), class = "prior")
}
