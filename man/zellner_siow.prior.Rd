\name{zellner_siow.prior}
\alias{zellner_siow.prior}
\title{Zellner-Siow Regression Coefficient Prior}
\usage{
  zellner_siow.prior(scale = 1, eff=TRUE)
}
\arguments{
  \item{scale}{Scaling on \eqn{g}. Must be positive.}
  \item{eff}{Logical for whether the prior precision for a model is scaled by its effective number of parameters. See description in \code{\link{coef.prior-class}}.}
}
\value{
  Returns an object of class "coef.prior", which has the following slots.
  \item{family}{Value is the string \code{"zellner-siow"}.}
  \item{hyper.parameters}{Value is a list with component \code{scale}.}
  \item{eff}{Value is the user supplied value of \code{eff}.}
}
\description{
  Creates an object representing the Zellner-Siow prior distribution on regresion coefficients.
}
\details{
  The prior distribution on 
  \deqn{
    w=\dfrac{\code{scale}}{g}
    \quad\text{is}\quad
    w\sim\chi^2_1.
  }
}
\seealso{
  \code{\link{lm.mcmc}()},\cr
  \code{\link{coef.prior-class}},\cr
  \code{\link{hyper_g.prior}()},\cr
  \code{\link{intrinsic.prior}()},\cr
  \code{\link{inverse_gamma.prior}()},\cr
  \code{\link{beta_prime.prior}()},\cr
  \code{\link{inverse_beta.prior}()}.
}
\author{
  Andrew Womack
}
