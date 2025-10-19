\name{model.prior-class}
\alias{model.prior-class}
\docType{class}
\title{
Class \code{model.prior} for Model Space Priors for Regression Models
}
\description{
A class with slots for \code{family}, \code{hyper.parameters}, and \code{trunc}. 
Objects of this class are created by the functions listed in \strong{See Also}.
These objects are interpreted by \code{\link{lm.mcmc}()} to compute model space priors.
There are no methods for the class.
}

\details{
  The \code{family} slot is a string. 
  
  The \code{hyper.parameters} slot is a list.
  
  The \code{trunc} slot is either \code{NA} or a positive integer. 
  If there are \eqn{p} predictors 
  and \code{trunc} \eqn{<p}, then models with \eqn{k>} \code{trunc} predictors get 
  prior probability 0 and the prior is renormalized. 
  If \code{trunc} \eqn{=} \code{NA}, then \code{trunc} is set to \eqn{p}.
  
  All of the priors specify the probability of a model \eqn{M} with \eqn{k} predictors as
  
  \deqn{
    P(M) = P(|M| = k)\times \binom{p}{k}^{-1}
  }
  
  where \eqn{p} is the total number of predictors under consideration, \eqn{M} is a model, and \eqn{|M|} is the number of predictors in model \eqn{M} (the so-called \emph{complexity} of model \eqn{M}). 
  
  The Binomial, Beta-Binomial, and (truncated) Negative-Binomial families have complexity penalization modifications 
  which explicitly use \eqn{p} in the definition of the prior on model complexity. 
  
  The Matryoshka-Doll family focuses on universal multiplicity correction and, as a consequence, implicitly penalizes model complexity.
  
  Some priors that arise as special cases are:
  \itemize{
    \item{The Uniform prior can be obtained using \code{\link{binomial.prior}()} with \code{prob=0.5}.}
    \item{The (truncated) Poisson\code{(lambda)} prior on model complexity can be obtained using \code{\link{matryoshka_doll.prior}()} 
    with arguments \code{prob_ratio = 1/lambda} and \code{type="children"}.}
    \item{The (truncated) Geometric\code{(prob)} prior on model complexity can be obtained using 
    \code{\link{negative_binomial.prior}()} with argument \code{size=1}.}
    \item{Custom sequences for log-probabilities on model complexity can be supplied by using
    \code{\link{custom.prior}()}.}
  }
}

\author{
  Andrew Womack
}

\seealso{
  \code{\link{lm.mcmc}()}, \cr
  \code{\link{matryoshka_doll.prior}()}, \cr
  \code{\link{binomial.prior}()},\cr
  \code{\link{binomial_complexity.prior}()}, \cr
  \code{\link{beta_binomial.prior}()}, \cr
  \code{\link{beta_binomial_complexity.prior}()}, \cr
  \code{\link{negative_binomial.prior}()}, \cr
  \code{\link{negative_binomial_complexity.prior}()}, \cr
  \code{\link{custom.prior}()}.
}
