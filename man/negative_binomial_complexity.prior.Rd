\name{negative_binomial_complexity.prior}
\alias{negative_binomial_complexity.prior}
\title{Negative-Binomial with Complexity Penalization Model Space Prior}
\usage{
  negative_binomial_complexity.prior(size = 1, intercept = 0, slope = 1, power = 1, trunc=NA)
}
\arguments{
  \item{size}{The number of successes defining the stopping rule for the Negative-Binomial distribution. Must be between positive.}
  \item{intercept}{The intercept of \code{1/(1-prob)} for the Negative-Binomial distribution. Must be non-negative.}
  \item{slope}{The slope of \code{1/(1-prob)} for the Negative-Binomial distribution. Must be >=1.}
  \item{power}{The power on \eqn{p} (the number of predictors) \code{1/(1-prob)} for the Negative-Binomial distribution. Must be positive.}
  \item{trunc}{Models with more than \code{trunc} predictors are assigned 0 prior probability. 
  Must be a positive integer or \code{NA}.}
}
\value{
  Returns an object of class "model.prior", which has the following slots.
  \item{family}{Value is the string \code{"negative-binomial-complexity"}.}
  \item{hyper.parameters}{Value is a list with components \code{size}, \code{intercept}, \code{slope}, and \code{power}.}
  \item{trunc}{Value is the user supplied value of \code{trunc}.}
}
\description{
  Creates an object representing the Negative-Binomial (with complexity penalization) prior distribution on regresion models.
}
\details{
  The prior on the model space is a Negative-Binomial\code{(size, prob)} prior (see (\code{\link{negative_binomial.prior}()})) where 
  the probability parameter \code{prob} depends on the total number of predictors \eqn{p}. 
  The defining equation for \code{prob}
  
  \deqn{
    \dfrac{1}{1-\code{prob}}= \code{intercept} + \code{slope}*p^{\code{power}}.
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
  \code{\link{negative_binomial.prior}()}, \cr
  \code{\link{custom.prior}()}.
}
\author{
  Andrew Womack
}
