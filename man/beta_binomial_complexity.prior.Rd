\name{beta_binomial_complexity.prior}
\alias{beta_binomial_complexity.prior}
\title{Beta-Binomial with Complexity Penalization Model Space Prior}
\usage{
  beta_binomial_complexity.prior(alpha = 1, intercept = 0, slope = 1, power = 1, trunc=NA)
}
\arguments{
  \item{alpha}{The shape of the Beta mixing distribution at 0. Must be positive.}
  \item{intercept}{The intercept of the shape of the Beta mixing distribution at 1. Must be non-negative.}
  \item{slope}{The slope of the shape of the Beta mixing distribution at 1. Must be >=1.}
  \item{power}{The power on \eqn{p} (the number of predictors) of the shape of the Beta mixing distribution at 1. Must be positive.}
  \item{trunc}{Models with more than \code{trunc} predictors are assigned 0 prior probability. 
  Must be a positive integer or \code{NA}.}
}
\value{
  Returns an object of class "model.prior", which has the following slots.
  \item{family}{Value is the string \code{"beta-binomial-complexity"}.}
  \item{hyper.parameters}{Value is a list with components \code{alpha}, \code{intercept}, \code{slope}, and \code{power}.}
  \item{trunc}{Value is the user supplied value of \code{trunc}.}
}
\description{
  Creates an object representing the Beta-Binomial (with complexity penalization) prior distribution on regresion models.
}
\details{
  The prior on the model space is a Beta-Binomial\eqn{(\alpha, \beta)}prior (see \code{\link{beta_binomial.prior}()}) where 
  the shape parameter \eqn{\beta} depends on the total number of predictors \eqn{p}. 
  This shape parameter is defined by
  
  \deqn{\beta = \code{intercept} + \code{slope}*p^{\code{power}}}
  
  If \code{trunc} \eqn{<p}, then models with \eqn{k>} \code{trunc} predictors get 
  prior probability 0 and the prior is renormalized. 
  If \code{trunc} \eqn{=} \code{NA}, then \code{trunc} is set to \eqn{p}.
}
\seealso{
  \code{\link{lm.mcmc}()},\cr
  \code{\link{model.prior-class}},\cr
  \code{\link{matryoshka_doll.prior}()}, \cr
  \code{\link{binomial.prior}()}, \cr
  \code{\link{binomial_complexity.prior}()}, \cr
  \code{\link{beta_binomial.prior}()},\cr 
  \code{\link{negative_binomial.prior}()},\cr 
  \code{\link{negative_binomial_complexity.prior}()}, \cr
  \code{\link{custom.prior}()}.
}
\author{
  Andrew Womack
}
