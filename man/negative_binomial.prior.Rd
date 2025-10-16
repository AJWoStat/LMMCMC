\name{negative_binomial.prior}
\alias{negative_binomial.prior}
\title{Negative-Binomial Model Space Prior}
\usage{
  negative_binomial.prior(size = 1, prob = 0.5, trunc = NA)
}
\arguments{
  \item{size}{The number of successes defining the stopping rule for the Negative-Binomial distribution. Must be between positive.}
  \item{prob}{The probability of success for each trial for the Negative-Binomial distribution. Must be between 0 and 1.}
  \item{trunc}{Models with more than \code{trunc} predictors are assigned 0 prior probability. 
  Must be a positive integer or \code{NA}.}
}
\value{
  Returns an object of class "model.prior", which has the following slots.
  \item{family}{Value is the string \code{"negative-binomial"}.}
  \item{hyper.parameters}{Value is a list with components \code{size} and \code{prob}.}
  \item{trunc}{Value is the user supplied value of \code{trunc}.}
}
\description{
  Creates an object representing the Negative-Binomial prior distribution on regresion models.
}
\details{
  The Negative-Binomial prior is defined on the non-negative integers and is thus always truncated.
  The prior probability for models with complexity \eqn{k} when there are \eqn{p} predictors is
  
  \deqn{
    P(|M| = k) \propto \dfrac{\Gamma(\code{size} + k)}{k!}\times (1-\code{prob})^{k}.
  } 
  
  If \code{trunc} \eqn{<p}, then models with \eqn{k>} \code{trunc} predictors get 
  prior probability 0 and the prior is renormalized. 
  If \code{trunc} \eqn{=} \code{NA}, then \code{trunc} is set to \eqn{p}.
}
\seealso{
  \code{\link{lm.mcmc}()}, \cr
  \code{\link{model.prior-class}},\cr
  \code{\link{matryoshka_doll.prior}()}, \cr
  \code{\link{binomial.prior}()}, \cr
  \code{\link{binomial_complexity.prior}()}, \cr
  \code{\link{beta_binomial.prior}()}, \cr
  \code{\link{beta_binomial_complexity.prior}()}, \cr
  \code{\link{negative_binomial_complexity.prior}()}, \cr
  \code{\link{custom.prior}()}.
}
\author{
  Andrew Womack
}
