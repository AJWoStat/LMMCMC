\name{inverse_beta.prior}
\alias{inverse_beta.prior}
\title{Inverse-Beta Regression Coefficient Prior}
\usage{
  inverse_beta.prior(scale = 1, shape_0 = 0.5, shape_1 = 0.5, eff=TRUE)
}
\arguments{
  \item{scale}{Scaling on \eqn{g}. Must be positive.}
  \item{shape_0}{Shape parameter at 0 for prior on \eqn{\frac{\code{scale}}{g}}. Must be positive.}
  \item{shape_1}{Shape parameter at 1 for prior on \eqn{\frac{\code{scale}}{g}}. Must be positive.}
  \item{eff}{Logical for whether the prior precision for a model is scaled by its effective number of parameters. See description in \code{\link{coef.prior-class}}.}
}
\value{
  Returns an object of class "coef.prior", which has the following slots.
  \item{family}{Value is the string \code{"inverse-beta"}.}
  \item{hyper.parameters}{Value is a list with component \code{scale}, \code{shape\_0}, and \code{shape\_1}.}
  \item{eff}{Value is the user supplied value of \code{eff}.}
}
\description{
  Creates an object representing the Inverse-Beta prior distribution on regresion coefficients.
}
\details{
  The prior distribution on 
  \deqn{
    w=\dfrac{\code{scale}}{g}
    \quad\text{is}\quad
    w\sim\text{Beta}\left(\code{shape\_0},\code{shape\_1}\right).
  }
}
\seealso{
  \code{\link{lm.mcmc}()},\cr
  \code{\link{coef.prior-class}},\cr
  \code{\link{zellner_siow.prior}()},\cr
  \code{\link{hyper_g.prior}()},\cr
  \code{\link{intrinsic.prior}()},\cr
  \code{\link{inverse_gamma.prior}()},\cr
  \code{\link{beta_prime.prior}()}.
}
\author{
  Andrew Womack
}
