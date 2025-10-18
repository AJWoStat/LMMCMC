\name{custom.prior}
\alias{custom.prior}
\title{Custom Model Space Prior}
\usage{
  custom.prior(log.complexity.prior, family = "custom", trunc = NA)
}
\arguments{
  \item{log.complexity.prior}{Vector of logarithms of probabilities on model complexity. Should have length p+1 where p is the number of predictors in the full model. Index k+1 is the logarithm of the probability that k predictors are included in the model.}
  \item{family}{A string representing the name of the prior.}
  \item{trunc}{Models with more than \code{trunc} predictors are assigned 0 prior probability. 
  Must be a positive integer or \code{NA}.}
}
\value{
  Returns an object of class "prior" with the following components.
  \item{family}{Value is user supplied.}
  \item{hyper.parameters}{Value is a list with component \code{log.complexity.prior}.}
}
\description{
  Creates an object representing a custom prior distribution on models. User beware: minimal checks are performed on the input.
}
\seealso{
  \code{\link{lm.mcmc}()}, \cr
  \code{\link{model.prior-class}},\cr
  \code{\link{matryoshka_doll.prior}()},\cr
  \code{\link{binomial.prior}()}, \cr
  \code{\link{binomial_complexity.prior}()}, \cr
  \code{\link{beta_binomial.prior}()}, \cr
  \code{\link{beta_binomial_complexity.prior}()}, \cr
  \code{\link{negative_binomial.prior}()}, \cr
  \code{\link{negative_binomial_complexity.prior}()}.
}
\author{
  Andrew Womack
}
\examples{
  s = 3 #power
  p = 123 #number of predictors
  log_prior = -s*log(c(0:p) + 1)
  custom.prior(log.complexity.prior = log_prior, family = "zeta-3")
}
