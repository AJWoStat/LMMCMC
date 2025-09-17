\name{md.prior}
\alias{md.prior}
\title{Matryoshka Doll prior}
\usage{
  md.prior(eta=1)
}
\arguments{
  \item{eta}{The odds of a model against the models that nest it. Must be positive, forced to be 1 if a negative number is given.}
}
\value{
  Returns an object of class "prior", with the following components.
  \item{family}{Value is \code{"MatryoshkaDoll"}.}
  \item{hyper.parameters}{Value is \code{eta}.}
}
\description{
  Creates an object representing the Matryoshka Doll prior distribution on models in the style BAS.
}
\details{
  The Matryoshka Doll prior fixes the odds of a model against the models that nest it.
  This achieves a Bayesian version of multiplicity correction by treating each model as
  a local null hypothesis with set of alternatives that are the models that nest it.
}
\seealso{
  \code{\link{lm.mcmc}()}, \code{\link{custom.prior}()}, \code{\link[BAS]{uniform}()},
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
  md.prior(1)
}
