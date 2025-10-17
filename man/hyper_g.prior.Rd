\name{hyper_g.prior}
\alias{hyper_g.prior}
\title{Hyper_g Regression Coefficient Prior}
\usage{
  hyper_g.prior(scale = 1, shape = 3, eff=TRUE)
}
\arguments{
  \item{scale}{Scaling on \eqn{g}. Must be positive.}
  \item{shape}{Shape parameter for prior on \eqn{g}. Must be >2.}
  \item{eff}{Logical for whether the prior precision for a model is scaled by its effective number of parameters. See description in \code{\link{coef.prior-class}}.}
}
\value{
  Returns an object of class "coef.prior", which has the following slots.
  \item{family}{Value is the string \code{"hyper-g"}.}
  \item{hyper.parameters}{Value is a list with component \code{scale} and \code{shape}.}
  \item{eff}{Value is the user supplied value of \code{eff}.}
}
\description{
  Creates an object representing the Hyper-g prior distribution on regresion coefficients.
}
\details{
  The prior distribution on 
  \deqn{
  \zeta=\dfrac{g}{\code{scale}+g}
  \quad\text{is}\quad
  \zeta\sim\text{Beta}\left(1,\dfrac{\code{shape}}{2} - 1\right).
  } 
}
\seealso{
  \code{\link{lm.mcmc}()},\cr
  \code{\link{coef.prior-class}},\cr
  \code{\link{zellner_siow.prior}()},\cr
  \code{\link{intrinsic.prior}()},\cr
  \code{\link{inverse_gamma.prior}()},\cr
  \code{\link{beta_prime.prior}()},\cr
  \code{\link{inverse_beta.prior}()}.
}
\author{
  Andrew Womack
}
