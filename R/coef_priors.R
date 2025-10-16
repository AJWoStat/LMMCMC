setClass("coef.prior", slots = c(family = "character", hyper.parameters = "list", eff = "logical"))

zellner_siow.prior = function(scale = 1.00, eff = TRUE){
  if (scale<=0 || !is.numeric(scale)) stop("zellner_siow.prior: scale must be > 0.")
  if (!is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("zellner_siow.prior: eff must be TRUE or FALSE.")
  structure(
    list(
      family = "zellner-siow",
      hyper.parameters = list(scale = scale),
      eff = eff
    ),
    class = "coef.prior"
  )
}

hyper_g.prior = function(scale = 1.00, shape = 3.00, eff = TRUE){
  if (scale<=0 || !is.numeric(scale)) stop("hyper_g.prior: scale must be > 0.")
  if (shape<=2 || !is.numeric(shape)) stop("hyper_g.prior: shape must be > 2.")
  if (!is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("hyper_g.prior: eff must be TRUE or FALSE.")
  structure(
    list(
      family = "hyper-g",
      hyper.parameters = list(scale = scale, shape = shape),
      eff = eff
    ),
    class = "coef.prior"
  )
}

intrinsic.prior = function(scale = 1.00, eff = TRUE){
  if (scale<=0 || !is.numeric(scale)) stop("intrinsic.prior: scale must be > 0.")
  if (!is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("intrinsic.prior: eff must be TRUE or FALSE.")
  structure(
    list(
      family = "intrinsic",
      hyper.parameters = list(scale = scale),
      eff = eff
    ),
    class = "coef.prior"
  )
}

gamma.prior = function(scale = 1.00, shape = 0.5, eff = TRUE){
  if (scale<=0 || !is.numeric(scale)) stop("gamma.prior: scale must be > 0.")
  if (shape<=0 || !is.numeric(shape)) stop("gamma.prior: shape must be > 0.")
  if (!is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("gamma.prior: eff must be TRUE or FALSE.")
  structure(
    list(
      family = "gamma",
      hyper.parameters = list(scale = scale, shape = shape),
      eff = eff
    ),
    class = "coef.prior"
  )
}

beta_prime.prior = function(scale = 1.00, shape_0 = 0.5, shape_inf = 1.00, eff = TRUE){
  if (scale<=0 || !is.numeric(scale)) stop("beta_prime.prior: scale must be > 0.")
  if (shape_0<=0 || !is.numeric(shape_0)) stop("beta_prime.prior: shape_0 must be > 0.")
  if (shape_inf<=0 || !is.numeric(shape_inf)) stop("beta_prime.prior: shape_inf must be > 0.")
  if (!is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("beta_prime.prior: eff must be TRUE or FALSE.")
  structure(
    list(
      family = "beta-prime",
      hyper.parameters = list(scale = scale, shape_0 = shape_0, shape_inf = shape_inf),
      eff = eff
    ),
    class = "coef.prior"
  )
}

scaled_beta.prior = function(scale = 1.00, shape_0 = 0.5, shape_1 = 0.5, eff = TRUE){
  if (scale<=0 || !is.numeric(scale)) stop("scaled_beta.prior: scale must be > 0.")
  if (shape_0<=0 || !is.numeric(shape_0)) stop("scaled_beta.prior: shape_0 must be > 0.")
  if (shape_1<=0 || !is.numeric(shape_1)) stop("scaled_beta.prior: shape_1 must be > 0.")
  if (!is.logical(eff) || (eff != TRUE && eff != FALSE)) stop("scaled_beta.prior: eff must be TRUE or FALSE.")
  structure(
    list(
      family = "scaled-beta",
      hyper.parameters = list(scale = scale, shape_0 = shape_0, shape_1 = shape_1),
      eff = eff
    ),
    class = "coef.prior"
  )
}

