\name{beta_binomial.prior}
\alias{beta_binomial.prior}
\title{Beta-Binomial Model Space Prior}
\usage{
  beta_binomial.prior(alpha = 1, beta = 1, trunc=NA)
}
\arguments{
  \item{alpha}{The shape of the Beta mixing distribution at 0. Must be positive.}
  \item{beta}{The shape of the Beta mixing distribution at 1. Must be positive.}
  \item{trunc}{Models with more than \code{trunc} predictors are assigned 0 prior probability. 
  Must be a positive integer or \code{NA}.}
}
\value{
  Returns an object of class "model.prior", which has the following slots.
  \item{family}{Value is the string \code{"beta-binomial"}.}
  \item{hyper.parameters}{Value is a list with components \code{alpha} and \code{beta}.}
  \item{trunc}{Value is the user supplied value of \code{trunc}.}
}
\description{
  Creates an object representing the Beta-Binomial prior distribution on regresion models.
}
\details{
  The Beta-Binomial prior treats the predictor inclusion indicators \eqn{\gamma_j} as \eqn{iid}
  Bernoulli\eqn{(\phi)} random variables and then gives \eqn{\phi} a Beta\eqn{(\alpha,\,\beta)}
  distribution. The probability of models with complexity \eqn{k} when there are \eqn{p} predictors is
  
  \deqn{
    P(|M| = k) = \dbinom{p}{k}\times \dfrac{B(\alpha+k,\,\beta + p-k)}{B(\alpha,\,\beta)}
  } 
  where
  \deqn{
    B(a, b) = \dfrac{\Gamma(a)\,\Gamma(b)}{\Gamma(a+b)}
  }
  is the Beta function. 
  
  If \code{trunc} \eqn{<p}, then models with \eqn{k>} \code{trunc} predictors get 
  prior probability 0 and the prior is renormalized. 
  If \code{trunc} \eqn{=} \code{NA}, then \code{trunc} is set to \eqn{p}.
}
\seealso{
  \code{\link{lm.mcmc}()}, \cr
  \code{\link{model.prior-class}},\cr
  \code{\link{matryoshka_doll.prior}()}, \cr
  \code{\link{binomial.prior}()},\cr 
  \code{\link{binomial_complexity.prior}()},\cr 
  \code{\link{beta_binomial_complexity.prior}()}, \cr
  \code{\link{negative_binomial.prior}()}, \cr\cr
  \code{\link{negative_binomial_complexity.prior}()}, \cr
  \code{\link{custom.prior}()}.
}
\author{
  Andrew Womack
}
