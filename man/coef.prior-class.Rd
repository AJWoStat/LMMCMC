\name{coef.prior-class}
\alias{coef.prior-class}
\docType{class}
\title{
Class \code{coef.prior} for Regression Coefficients
}
\description{
A class with slots for \code{family}, \code{hyper.parameters}, and \code{eff}. 
Objects of this class are created by the functions listed in \strong{See Also}.
These objects are interpreted by \code{\link{lm.mcmc}()} to compute regression coefficient priors.
There are no methods for the class.
}

\details{
  The \code{family} slot is a string. 
  
  The \code{hyper.parameters} slot is a list.
  
  The \code{eff} slot is a logical. 
  
  
  For a model \eqn{M} with design matrix \eqn{(\boldsymbol{1}\ \mathbf{X}_M)} where \eqn{\mathbf{X}_M^\prime \boldsymbol{1} = \boldsymbol{0}}, the response \eqn{\mathbf{y}=(y_1, \ldots, y_n)} has model
  \deqn{
    \mathbf{y} = \alpha\boldsymbol{1} + \mathbf{X}_M\boldsymbol{\beta} + \frac{1}{\sqrt{\phi}}\boldsymbol{\epsilon}
  }
  where \eqn{\boldsymbol{\epsilon}\sim \mathcal{N}_n(\mathbf{0},\mathbf{I})} and \eqn{\phi} is the precision parameter. The prior for \eqn{(\alpha, \boldsymbol{\beta}, \phi)} is specified as a 
  null model based mixture of \eqn{g}-over-\eqn{n} priors. Specifically \eqn{(\alpha,\phi)\sim \phi^{-1}} and the prior for \eqn{\boldsymbol{\beta}} is
  \deqn{
    \boldsymbol{\beta}|g, \alpha, \phi, \mathbf{X}_M \sim \mathcal{N}_{|M|}\left(\mathbf{0}, \frac{gn}{\phi q_M}\left(\mathbf{X}_M^\prime \mathbf{X}_M\right)^{-1}\right).
  }
  The precision multiplier \eqn{q_M} is given by \eqn{q_M=\code{ifelse}(\code{eff}, |M|, 1)}.
  Setting \code{eff=FALSE} produces scaled mixtures of the unit-information prior.
  Setting \eqn{q_M = |M|} arises when 
  estimating the prior precision matrix when computing the Intrinsic Prior using minimal training samples.
  If \eqn{[n]=\{1,\ldots,n\}} and \eqn{\mathbf{Z}_S} for \eqn{S\subseteq [n]} is the design matrix whose rows
  correspond to those from \eqn{\mathbf{X}_M} whose indices are in \eqn{S}, then
  \deqn{
    \frac{q_M}{n}\,\mathbf{X}_M^\prime\mathbf{X}_M = \binom{n}{|M|}^{-1}
    \sum_{\substack{S\subset[n]\\|S|=|M|}} \mathbf{Z}_S^\prime \mathbf{Z}_S.
  }
  
  
  Three well-established scaled mixtures of \eqn{g}-priors are the Zellner-Siow, Hyper-\eqn{g}, and Intrinsic Priors. We include these priors on \eqn{g} and extend them.
  The Zellner-Siow Prior is extended by setting a Inverse-Gamma prior on \eqn{g}. 
  The Hyper-\eqn{g} Prior is extended by setting a Beta-Prime distribution on \eqn{g}. 
  The Intrinsic Prior is extended by setting a Inverse-Beta distribution on \eqn{g} (that is \eqn{w=1/g} follows a Beta distribution).
  
  All of the priors include a \code{scale} parameter as the first component of \code{hyperparameters}.
  Setting \code{eff=FALSE} and \code{scale=1/n} produces scaled mixtures of \eqn{g}-priors.
}

\author{
  Andrew Womack
}

\seealso{
  \code{\link{lm.mcmc}()}, \cr
  \code{\link{g.prior}()},\cr
  \code{\link{zellner_siow.prior}()},\cr
  \code{\link{hyper_g.prior}()},\cr
  \code{\link{intrinsic.prior}()},\cr
  \code{\link{inverse_gamma.prior}()},\cr
  \code{\link{beta_prime.prior}()},\cr
  \code{\link{inverse_beta.prior}()}.
  
}
