\name{binomial_complexity.prior}
\alias{binomial_complexity.prior}
\title{Binomial with Complexity Penalization Model Space Prior}
\usage{
  binomial_complexity.prior(intercept = 0, slope = 1, power = 1, trunc=NA)
}
\arguments{
  \item{intercept}{The intercept of the transformed probability of success. Must be non-negative.}
  \item{slope}{The slope of the transformed probability of success. Must be at >=1.}
  \item{power}{The power on \eqn{p} (the number of predictors) of the transformed probability of success. Must be positive.}
  \item{trunc}{Models with more than \code{trunc} predictors are assigned 0 prior probability. 
  Must be a positive integer or \code{NA}.}
}
\value{
  Returns an object of class "model.prior", which has the following slots.
  \item{family}{Value is the string \code{"binomial-complexity"}.}
  \item{hyper.parameters}{Value is a list with components \code{intercept}, \code{slope}, and \code{power}.}
  \item{trunc}{Value is the user supplied value of \code{trunc}.}
}
\description{
  Creates an object representing the Binomial (with complexity penalization) prior distribution on regresion models.
}
\details{
  The prior on the model space is a Binomial\eqn{(\code{prob})} prior (see \code{\link{binomial.prior}()}) where 
  the parameter \code{prob} depends on the total number of predictors \eqn{p}. 
  This parameter is defined by
  
  \deqn{\dfrac{1}{\code{prob}} = \code{intercept} + \code{slope}*p^{\code{power}}}
  
  If \code{trunc} \eqn{<p}, then models with \eqn{k>} \code{trunc} predictors get 
  prior probability 0 and the prior is renormalized. 
  If \code{trunc} \eqn{=} \code{NA}, then \code{trunc} is set to \eqn{p}.
}
\seealso{
  \code{\link{lm.mcmc}()},\cr
  \code{\link{model.prior-class}},\cr
  \code{\link{matryoshka_doll.prior}()},\cr
  \code{\link{binomial.prior}()},\cr
  \code{\link{beta_binomial.prior}()},\cr
  \code{\link{beta_binomial_complexity.prior}()},\cr
  \code{\link{negative_binomial.prior}()},\cr
  \code{\link{negative_binomial_complexity.prior}()},\cr
  \code{\link{custom.prior}()}.
}
\author{
  Andrew Womack
}
