setClass("coef.prior", slots = c(family = "character", hyper.parameters = "list", eff = "logical"))

zellner_siow.prior = function(scale = 1, eff = TRUE){
  if (is.na(scale) || !is.numeric(scale) || scale<=0) stop("scale must be > 0.")
  if (is.na(eff) || !is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("eff must be TRUE or FALSE.")
  structure(
    list(
      family = "zellner-siow",
      hyper.parameters = list(scale = scale),
      eff = eff
    ),
    class = "coef.prior"
  )
}

hyper_g.prior = function(scale = 1, shape = 3, eff = TRUE){
  if (is.na(scale) || !is.numeric(scale) || scale<=0) stop("scale must be > 0.")
  if (is.na(shape) || !is.numeric(shape) || shape<=0) stop("shape must be > 2.")
  if (is.na(eff) || !is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("eff must be TRUE or FALSE.")
  structure(
    list(
      family = "hyper-g",
      hyper.parameters = list(scale = scale, shape = shape),
      eff = eff
    ),
    class = "coef.prior"
  )
}

intrinsic.prior = function(scale = 1, eff = TRUE){
  if (is.na(scale) || !is.numeric(scale) || scale<=0) stop("scale must be > 0.")
  if (is.na(eff) || !is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("eff must be TRUE or FALSE.")
  structure(
    list(
      family = "intrinsic",
      hyper.parameters = list(scale = scale),
      eff = eff
    ),
    class = "coef.prior"
  )
}

inverse_gamma.prior = function(scale = 1, shape = 0.5, eff = TRUE){
  if (is.na(scale) || !is.numeric(scale) || scale<=0) stop("scale must be > 0.")
  if (is.na(shape) || !is.numeric(shape) || shape<=0) stop("shape must be > 0.")
  if (is.na(eff) || !is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("eff must be TRUE or FALSE.")
  structure(
    list(
      family = "inverse_gamma",
      hyper.parameters = list(scale = scale, shape = shape),
      eff = eff
    ),
    class = "coef.prior"
  )
}

beta_prime.prior = function(scale = 1, shape_0 = 1, shape_1 = 0.5, eff = TRUE){
  if (is.na(scale) || !is.numeric(scale) || scale<=0) stop("scale must be > 0.")
  if (is.na(shape_0) || !is.numeric(shape_0) || shape_0<=0) stop("shape_0 must be > 0.")
  if (is.na(shape_1) || !is.numeric(shape_1) || shape_1<=0) stop("shape_1 must be > 0.")
  if (is.na(eff) || !is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("eff must be TRUE or FALSE.")
  structure(
    list(
      family = "beta-prime",
      hyper.parameters = list(scale = scale, shape_0 = shape_0, shape_1 = shape_1),
      eff = eff
    ),
    class = "coef.prior"
  )
}

inverse_beta.prior = function(scale = 1, shape_0 = 0.5, shape_1 = 0.5, eff = TRUE){
  if (is.na(scale) || !is.numeric(scale) || scale<=0) stop("scale must be > 0.")
  if (is.na(shape_0) || !is.numeric(shape_0) || shape_0<=0) stop("shape_0 must be > 0.")
  if (is.na(shape_1) || !is.numeric(shape_1) || shape_1<=0) stop("shape_1 must be > 0.")
  if (is.na(eff) || !is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("eff must be TRUE or FALSE.")
  structure(
    list(
      family = "inverse-beta",
      hyper.parameters = list(scale = scale, shape_0 = shape_0, shape_1 = shape_1),
      eff = eff
    ),
    class = "coef.prior"
  )
}

