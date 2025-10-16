\name{binomial.prior}
\alias{binomial.prior}
\title{Binomial Model Space Prior}
\usage{
  binomial.prior(prob=0.5, trunc=NA)
}
\arguments{
  \item{prob}{The probability of success for the Binomial distribution. Must be between 0 and 1.}
  \item{trunc}{Models with more than \code{trunc} predictors are assigned 0 prior probability. 
  Must be a positive integer or \code{NA}.}
}
\value{
  Returns an object of class "model.prior", which has the following slots.
  \item{family}{Value is the string \code{"binomial"}.}
  \item{hyper.parameters}{Value is a list with component \code{prob}.}
  \item{trunc}{Value is the user supplied value of \code{trunc}.}
}
\description{
  Creates an object representing the Binomial prior distribution on regresion models.
}
\details{
  The Binomial prior treats the predictor inclusion indicators \eqn{\gamma_j} as \eqn{iid}
  Bernoulli\eqn{(\code{prob})} random variables. The probability of models with complexity \eqn{k} when there are \eqn{p} predictors is
  
  \deqn{
    P(|M| = k) = \dbinom{p}{k}\times \code{prob}^k(1-\code{prob})^{p-k}.
  } 
  
  If \code{trunc} \eqn{<p}, then models with \eqn{k>} \code{trunc} predictors get 
  prior probability 0 and the prior is renormalized. 
  If \code{trunc} \eqn{=} \code{NA}, then \code{trunc} is set to \eqn{p}.
}
\seealso{
  \code{\link{lm.mcmc}()},\cr
  \code{\link{model.prior-class}},\cr
  \code{\link{matryoshka_doll.prior}()},\cr
  \code{\link{binomial_complexity.prior}()},\cr
  \code{\link{beta_binomial.prior}()},\cr
  \code{\link{beta_binomial_complexity.prior}()},\cr
  \code{\link{negative_binomial.prior}()},\cr
  \code{\link{negative_binomial_complexity.prior}()},\cr
  \code{\link{custom.prior}()}.
}
\author{
  Andrew Womack
}
