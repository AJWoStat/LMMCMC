\name{lm.mcmc}
\alias{lm.mcmc}
\title{MCMC for Gaussian linear models}
\description{
Proposal kernels include single random inclusion indicator change,
one in and one out inclusion indicator changes, and
random inclusion indicator vector replacement.
These are called with proportions 0.45, 0.45, and 0.1, respectively.
Current implementation uses a g/n-prior with g=1.
Starting model is the minimal forward stepwise BIC model (with a maximum of up to 20 predictors).
}
\usage{
lm.mcmc(y, X, model.prior = md.prior(1), burnin.iterations=1e4, mcmc.iterations=1e6, thin=1)
}
\arguments{
\item{y}{Response vector.}
\item{X}{Design matrix. A column for the intercept is added.}
\item{model.prior}{An object of class "prior" following the construction from BAS, see Details.}
\item{burnin.iterations}{Number of burn-in iterations.}
\item{mcmc.iterations}{Number of MCMC iterations to return.}
\item{thin}{Number (plus 1) of iterations to skip between MCMC recording.}
}
\value{
A list in the style of BAS. All but \code{mcmc.draws} share a common index.

\item{which}{(list) Model vector of indices for included predictors. Includes a 0 for the intercept.}
\item{logmarg}{(vector) Model logarithm of the Bayes factor to the intercept only model.}
\item{priorprobs}{(vector) Model prior probabilities.}
\item{postprobs.RN}{(vector) Model renormalized posterior probability.}
\item{postprobs}{(vector) Model MCMC frequency probability.}
\item{size}{(vector) Model number of predictors (including the intercept).}
\item{rank}{(vector) Model design matrix rank.}
\item{R2}{(vector) Model $R^2$.}
\item{mse}{(vector) Model response variance.}
\item{freq}{(vector) Model MCMC frequency.}
\item{mcmc.id}{(vector) Model ID in the sequence of MCMC draws. Unsampled models get the common value of 0.}
\item{mcmc.draws}{(vector) Sequence of integers representing the models sampled. Each element is a non-zero element of mcmc.id}
}

\details{
  Priors that work as \code{model.prior} input to \code{lm.mcmc}:

  From this package: \code{\link{md.prior}()} and those built using \code{\link{custom.prior}()}.

  Imported from BAS:
  \code{\link[BAS]{uniform}()},
  \code{\link[BAS]{Bernoulli}()},
  \code{\link[BAS]{beta.binomial}()},
  \code{\link[BAS]{tr.poisson}()},
  \code{\link[BAS]{tr.beta.binomial}()},
  \code{\link[BAS]{tr.power.prior}()}.


}

\seealso{
  \code{\link{md.prior}()}, \code{\link{custom.prior}()}, \code{\link[BAS]{uniform}()},
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
n = 100
p = 5
X = matrix(rnorm(n*p), n, p)
y = rnorm(n) + 1 + X\%*\%rep(c(1,0), times = c(2, p-2))
out = lm.mcmc(y, X, model.prior = md.prior(1), burnin.iterations=1e3, mcmc.iterations=1e4, thin=1)
}
