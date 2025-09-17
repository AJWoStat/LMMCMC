\name{custom.prior}
\alias{custom.prior}
\title{Custom model space prior}
\usage{
  custom.prior(family, log.complexity.prior)
}
\arguments{
  \item{family}{A string representing the name of the prior.}
  \item{log.complexity.prior}{Vector of logarithms of probabilities on model complexity. Should have length p+1 where p is the number of predictors in the full model. Index k+1 is the logarithm of the probability that k predictors are included in the model.}
}
\value{
  Returns an object of class "prior" with the following components.
  \item{family}{Value is user supplied.}
  \item{hyper.parameters}{Value is the user supplied logarithm of probabilities on model complexity.}
}
\description{
  Creates an object representing a custom prior distribution on models in the style BAS. User beware: no checks are performed on the input.
}
\seealso{
  \code{\link{lm.mcmc}()}, \code{\link{md.prior}()}, \code{\link[BAS]{uniform}()},
  \code{\link[BAS]{Bernoulli}()},
  \code{\link[BAS]{beta.binomial}()},
  \code{\link[BAS]{tr.poisson}()},
  \code{\link[BAS]{tr.beta.binomial}()},
  \code{\link[BAS]{tr.power.prior}()}.
}
\author{
  Andrew Womack
}
\examples{
  custom.prior(family="tr.geom", log.complexity.prior = dgeom(0:10, 0.5, log=TRUE) - pgeom(10, 0.5, log.p=TRUE))
}
