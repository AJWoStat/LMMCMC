lm.mcmc = function(
    y, X, base.model.indices = c(),
    model.prior = matryoshka_doll.prior(), coef.prior = hyper_g.prior(), 
    burnin.iterations=1e4, mcmc.iterations=1e6, thin=1, 
    start.model.indices = NA, max.size.start = 20
){
  weights = rep(1, length(y))
  base_model_indices = as.integer(base.model.indices)
  base_model_indices = base_model_indices[base_model_indices<=ncol(X)]
  base_model_indices = base_model_indices[base_model_indices>0]
  include_coef = 0
  include_vcov = 0
  if(any(is.na(start.model.indices))){
    warn_curr = options("warn")$warn
    options(warn=-1)
    aa = leaps::regsubsets(X,y, method="forward", force.in = base_model_indices,
                           nbest=1, nested = TRUE, nvmax=min(c(length(y)-1, dim(X)[2], max.size.start)))
    options(warn = warn_curr)
    bb = summary(aa)
    cc = which.max(bb$bic==min(bb$bic))
    start_model_indices = which(bb$which[cc,-1])
  }else{
    start_model_indices = sort(unique(c(start.model.indices, base_model_indices)))
  }
  proposal_probs = c(0.45, 0.45, 0.1)
  hash_table_length = 1e6
  linked_list_length = 5
  max_load_factor = 0.75
  if(model.prior$family=="matryoshka-doll") model.prior$hyper.parameters$type = ifelse(model.prior$hyper.parameters$type == "descendants", 1, 0)
  model.priorfamily = model.prior$family
  model.priorpar = unlist(model.prior$hyper.parameters)
  trunc = ifelse(
    !is.na(model.prior$trunc) && model.prior$trunc<(ncol(X) - length(base_model_indices)),
    model.prior$trunc, 
    ncol(X) - length(base_model_indices)
  )
  coef.priorfamily = coef.prior$family
  coef.priorpar = unlist(coef.prior$hyper.parameters)
  eff = ifelse(coef.prior$eff, 1, 0)

  result = .Call(
    "lm_mcmc_function", X, y, weights, as.integer(base_model_indices),
    model.priorfamily, model.priorpar, as.integer(trunc),
    coef.priorfamily, coef.priorpar, as.integer(eff),
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
