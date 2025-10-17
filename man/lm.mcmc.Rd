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
lm.mcmc(
  y, X, 
  model.prior = matryoshka_doll.prior(), 
  coef.prior = hyper_g.prior(),
  burnin.iterations=1e4, mcmc.iterations=1e6, thin=1
)
}
\arguments{
\item{y}{Response vector.}
\item{X}{Design matrix. A column for the intercept is added.}
\item{model.prior}{An object of class "model.prior", see Details.}
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
  Priors that work as \code{model.prior} input to \code{lm.mcmc} must be constructed from
  package-supplied functions that make objects of the class \code{model.prior} (see \code{\link{model.prior-class}}).
  Custom user-supplied vectors for \eqn{\log(P(|M|=k))} can be converted to this class using \code{\link{custom.prior}()}.
  
  The "number of predictors" \eqn{p} used in the descriptions in the help files for those model priors is \code{ncol(X)}.
}

\seealso{
  Model space priors:\cr
  \code{\link{model.prior-class}},\cr
  \code{\link{matryoshka_doll.prior}()}, \cr
  \code{\link{binomial.prior}()},\cr
  \code{\link{binomial_complexity.prior}()}, \cr
  \code{\link{beta_binomial.prior}()}, \cr
  \code{\link{beta_binomial_complexity.prior}()}, \cr
  \code{\link{negative_binomial.prior}()}, \cr
  \code{\link{negative_binomial_complexity.prior}()}, \cr
  \code{\link{custom.prior}()}.
  \cr\cr
  Regression coefficient priors:\cr
  \code{\link{coef.prior-class}},\cr
  \code{\link{g.prior}()},\cr
  \code{\link{zellner_siow.prior}()},\cr
  \code{\link{hyper_g.prior}()},\cr
  \code{\link{intrinsic.prior}()},\cr
  \code{\link{inverse_gamma.prior}()},\cr
  \code{\link{beta_prime.prior}()}\cr,
  \code{\link{inverse_beta.prior}()}.
  
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
