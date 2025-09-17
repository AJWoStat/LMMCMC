lm.mcmc = function(y, X, model.prior = md.prior(1), burnin.iterations=1e4, mcmc.iterations=1e6, thin=1){

  weights = rep(1, length(y))
  base_model_indices = c()
  include_coef = 0
  include_vcov = 0
  warn_curr = options("warn")$warn
  options(warn=-1)
  aa = leaps::regsubsets(X,y, method="forward", nbest=1, nested = TRUE, nvmax=min(c(length(y)-1, dim(X)[2], 20)))
  options(warn = warn_curr)
  bb = summary(aa)
  cc = which.max(bb$bic==min(bb$bic))
  start_model_indices = which(bb$which[cc,-1])
  proposal_probs = c(0.45, 0.45, 0.1)
  trunc = ncol(X)
  hash_table_length = 1e6
  linked_list_length = 5
  max_load_factor = 0.75
  if(model.prior$family=="Beta-Binomial"){
    model.priorfamily = "Beta-Binomial"
    model.priorpar = model.prior$hyper.parameters[1:2]
    trunc = -1
  }else if(model.prior$family=="Bernoulli"){
    model.priorfamily = "Bernoulli"
    model.priorpar = model.prior$hyper.parameters[1]
    trunc = -1
  }else if(model.prior$family=="MatryoshkaDoll"){
    model.priorfamily = "MatryoshkaDoll"
    model.priorpar = model.prior$hyper.parameters[1]
    trunc = -1
  }else if(model.prior$family=="Trunc-Beta-Binomial"){
    model.priorfamily = "Trunc-Beta-Binomial"
    model.priorpar = model.prior$hyper.parameters[1:2]
    trunc = model.prior$hyper.parameters[3]
  }else if(model.prior$family=="Trunc-Poisson"){
    model.priorfamily = "Trunc-Poisson"
    model.priorpar = model.prior$hyper.parameters[1]
    trunc = model.prior$hyper.parameters[2]
  }else if(model.prior$family=="Trunc-Power-Prior"){
    model.priorfamily = "Trunc-Power-Prior"
    model.priorpar = model.prior$hyper.parameters[1]
    trunc = model.prior$hyper.parameters[2]
  }else if(model.prior$family=="Uniform"){
    model.priorfamily = "Uniform"
    model.priorpar = 0.5
    trunc = -1
  }else{
    model.priorfamily = model.prior$family
    model.priorpar = model.prior$hyper.parameters
    trunc = -1
  }

  result = .Call(
    "lm_mcmc_function", X, y, weights, as.integer(base_model_indices),
    model.priorfamily, model.priorpar, as.integer(trunc),
    as.integer(include_coef), as.integer(include_vcov),
    as.integer(start_model_indices),
    as.integer(burnin.iterations), as.integer(thin), as.integer(mcmc.iterations), proposal_probs,
    as.integer(hash_table_length), as.integer(linked_list_length), max_load_factor,
    PACKAGE = "LMMCMC"
  )
  out = list()
  out$which = result$which
  out$logmarg = result$logmarg
  out$priorprobs = exp(result$logprior)
  out$postprobs.RN = exp(result$logpost_RN)
  out$postprobs = result$post_sampling
  out$size = result$size
  out$rank = result$rank
  out$R2 = result$Rsq
  out$mse = result$residual_sd^2
  out$freq = result$mcmc_count
  out$mcmc_id = result$mcmc_id+1
  out$mcmc_draws = result$mcmc_draws+1

  return(out)
}
