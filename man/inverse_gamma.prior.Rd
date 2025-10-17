\name{inverse_gamma.prior}
\alias{inverse_gamma.prior}
\title{Inverse-Gamma Regression Coefficient Prior}
\usage{
  inverse_gamma.prior(scale = 1, shape = 0.5, eff=TRUE)
}
\arguments{
  \item{scale}{Scaling on \eqn{g}. Must be positive.}
  \item{shape}{Shape parameter for prior on \eqn{g}. Must be positive.}
  \item{eff}{Logical for whether the prior precision for a model is scaled by its effective number of parameters. See description in \code{\link{coef.prior-class}}.}
}
\value{
  Returns an object of class "coef.prior", which has the following slots.
  \item{family}{Value is the string \code{"inverse_gamma"}.}
  \item{hyper.parameters}{Value is a list with component \code{scale} and \code{shape}.}
  \item{eff}{Value is the user supplied value of \code{eff}.}
}
\description{
  Creates an object representing the Inverse-Gamma prior distribution on regresion coefficients.
}
\details{
  The prior distribution on 
  \deqn{
    w=\dfrac{\code{scale}}{g}
    \quad\text{is}\quad
    w\sim\text{Gamma}(\code{shape}, 1).
  }
}
\seealso{
  \code{\link{lm.mcmc}()},\cr
  \code{\link{coef.prior-class}},\cr
  \code{\link{zellner_siow.prior}()},\cr
  \code{\link{hyper_g.prior}()},\cr
  \code{\link{intrinsic.prior}()},\cr
  \code{\link{beta_prime.prior}()},\cr
  \code{\link{inverse_beta.prior}()}.
}
\author{
  Andrew Womack
}
