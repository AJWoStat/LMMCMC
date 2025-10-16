setClass("model.prior", slots = c(family = "character", hyper.parameters = "list", trunc = "integer"))

matryoshka_doll.prior = function(prob_ratio = 1, type = 'descendants', trunc=NA){
  # truncated Poisson is given by type = "children" and Poisson rate is 1/prob_ratio
  if (is.na(prob_ratio) || !is.numeric(prob_ratio) || prob_ratio<=0) stop("prob_nesting_ratio must be > 0.")
  if (is.na(type) || !is.character(type) || (type != "descendants" && type != "children")) stop('type must be "descendants" or "children".')
  if (is.na(trunc)) trunc = as.integer(trunc)
  if (!is.integer(trunc)){
    remainder = trunc - floor(trunc)
    if(remainder < .Machine$double.eps){
      trunc = as.integer(round(trunc))
    }else{
      stop("trunc must be an integer.")
    }
  }
  if (trunc<=0 && !is.na(trunc)) stop("trunc must be >0 or NA.")
  
  structure(
    list(
      family = "matryoshka-doll",
      hyper.parameters = list(prob_ratio = prob_ratio, type = type),
      trunc = trunc
    ),
    class = "model.prior"
  )
}

binomial.prior = function(prob=0.5, trunc=NA){
  if (is.na(prob) || !is.numeric(prob) || prob<=0 || prob>=1) stop("prob must be between 0 and 1.")
  if (is.na(trunc)) trunc = as.integer(trunc)
  if (!is.integer(trunc)){
    remainder = trunc - floor(trunc)
    if(remainder < .Machine$double.eps){
      trunc = as.integer(round(trunc))
    }else{
      stop("trunc must be an integer.")
    }
  }
  if (trunc<=0 && !is.na(trunc)) stop("trunc must be >0 or NA.")
  
  structure(
    list(
      family = "binomial",
      hyper.parameters = list(prob=prob),
      trunc = trunc
    ),
    class = "model.prior"
  )
}

binomial_complexity.prior = function(intercept = 0, slope = 1, power = 1, trunc=NA){
  if (is.na(intercept) || !is.numeric(intercept) || intercept<0) stop("intercept must be >= 0.")
  if (is.na(slope) || !is.numeric(slope) || slope<1) stop("slope must be >= 1.")
  if (is.na(power) || !is.numeric(power) || power<=0) stop("power must be > 0.")
  if (is.na(trunc)) trunc = as.integer(trunc)
  if (!is.integer(trunc)){
    remainder = trunc - floor(trunc)
    if(remainder < .Machine$double.eps){
      trunc = as.integer(round(trunc))
    }else{
      stop("trunc must be an integer.")
    }
  }
  if (trunc<=0 && !is.na(trunc)) stop("trunc must be >0 or NA.")
  
  structure(
    list(
      family = "binomial-complexity",
      hyper.parameters = list(intercept = intercept, slope = slope, power = power),
      trunc = trunc
    ),
    class = "model.prior"
  )
}

beta_binomial.prior = function(alpha = 1, beta = 1, trunc=NA){
  if (is.na(alpha) || !is.numeric(alpha) || alpha<=0) stop("alpha must be > 0.")
  if (is.na(beta) || !is.numeric(beta) || beta<=0) stop("beta must be > 0.")
  if (is.na(trunc)) trunc = as.integer(trunc)
  if (!is.integer(trunc)){
    remainder = trunc - floor(trunc)
    if(remainder < .Machine$double.eps){
      trunc = as.integer(round(trunc))
    }else{
      stop("trunc must be an integer.")
    }
  }
  if (trunc<=0 && !is.na(trunc)) stop("trunc must be >0 or NA.")
  
  structure(
    list(
      family = "beta-binomial",
      hyper.parameters = list(alpha = alpha, beta = beta),
      trunc = trunc
    ),
    class = "model.prior"
  )
}

beta_binomial_complexity.prior = function(alpha = 1, intercept = 0, slope = 1, power = 1, trunc=NA){
  if (is.na(alpha) || !is.numeric(alpha) || alpha<=0) stop("alpha must be > 0.")
  if (is.na(intercept) || !is.numeric(intercept) || intercept<0) stop("intercept must be >= 0.")
  if (is.na(slope) || !is.numeric(slope) || slope<1) stop("slope must be >= 1.")
  if (is.na(power) || !is.numeric(power) || power<=0) stop("power must be > 0.")
  if (is.na(trunc)) trunc = as.integer(trunc)
  if (!is.integer(trunc)){
    remainder = trunc - floor(trunc)
    if(remainder < .Machine$double.eps){
      trunc = as.integer(round(trunc))
    }else{
      stop("trunc must be an integer.")
    }
  }
  if (trunc<=0 && !is.na(trunc)) stop("trunc must be >0 or NA.")
  
  structure(
    list(
      family = "beta-binomial-complexity",
      hyper.parameters = list(alpha = alpha, intercept = intercept, slope = slope, power = power),
      trunc = trunc
    ),
    class = "model.prior"
  )
  
}

negative_binomial.prior = function(size = 1, prob = 0.5, trunc = NA){
  # P(X=x)=gamma(shape+x)/(gamma(shape)*factorial(x))*prob^shape*(1-prob)^x
  if (is.na(size) || !is.numeric(size) || size<=0) stop("size must be > 0.")
  if (is.na(prob) || !is.numeric(prob) || prob<=0 || prob>=1) stop("prob must be between 0 and 1.")
  if (is.na(trunc)) trunc = as.integer(trunc)
  if (!is.integer(trunc)){
    remainder = trunc - floor(trunc)
    if(remainder < .Machine$double.eps){
      trunc = as.integer(round(trunc))
    }else{
      stop("trunc must be an integer.")
    }
  }
  if (trunc<=0 && !is.na(trunc)) stop("trunc must be >0 or NA.")
  
  structure(
    list(
      family = "negative-binomial",
      hyper.parameters = list(size = size, prob = prob),
      trunc = trunc
    ),
    class = "model.prior"
  )
}

negative_binomial_complexity.prior = function(size = 1, intercept = 0, slope = 1, power = 1, trunc=NA){
  # negative binomial with prob = 1 - 1/(intercept + slope*p^power)
  if (is.na(size) || !is.numeric(size) || size<=0) stop("size must be > 0.")
  if (is.na(intercept) || !is.numeric(intercept) || intercept<0) stop("intercept must be >= 0.")
  if (is.na(slope) || !is.numeric(slope) || slope<1) stop("slope must be >= 1.")
  if (is.na(power) || !is.numeric(power) || power<=0) stop("power must be > 0.")
  if (is.na(trunc)) trunc = as.integer(trunc)
  if (!is.integer(trunc)){
    remainder = trunc - floor(trunc)
    if(remainder < .Machine$double.eps){
      trunc = as.integer(round(trunc))
    }else{
      stop("trunc must be an integer.")
    }
  }
  if (trunc<=0 && !is.na(trunc)) stop("trunc must be >0 or NA.")
  
  structure(
    list(
      family = "negative-binomial-complexity",
      hyper.parameters = list(size = size, intercept = intercept, slope = slope, power = power),
      trunc = trunc
    ),
    class = "model.prior"
  )
}

custom.prior <- function(log.complexity.prior, family = "custom", trunc = NA) {
  if(is.na(family) || !is.character(family)) stop("family must be a string.")
  warning("User-defined prior on model complexity. Minimal checks performed.")
  if(!is.numeric(log.complexity.prior)) stop("log.complexity.prior must be numeric.")
  if(any(is.na(log.complexity.prior))) stop("NA not allowed in log.complexity.prior.")
  if(any(log.complexity.prior==Inf)) stop("Inf not allowed in log.complexity.prior.")
  m = max(log.complexity.prior)
  log_sum_exp = log(sum(exp(log.complexity.prior - m))) + m
  if(abs(log_sum_exp)>.Machine$double.eps){
    warning("log(sum(exp(log.complexity.prior))) != 0. Prior has been renormalized.")
    log.complexity.prior = log.complexity.prior - log_sum_exp
  }
  structure(
    list(
      family = family, 
      hyper.parameters = list(log.complexity.prior = log.complexity.prior),
      trunc = trunc
    ), 
    class = "prior"
  )
}

