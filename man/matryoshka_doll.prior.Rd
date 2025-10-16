\name{matryoshka_doll.prior}
\alias{matryoshka_doll.prior}
\title{Matryoshka Doll Model Space Prior}
\usage{
  matryoshka_doll.prior(prob_ratio = 1, type = "descendants", trunc = NA)
}
\arguments{
  \item{prob_ratio}{The ratio of the probability of a model versus a set of models that nest it. Must be positive.}
  \item{type}{Description of the set of models used in the definition: "descendants" for all nesting model and "children" for nesting models that only contain one additional covariate.}
  \item{trunc}{Models with more than \code{trunc} predictors are assigned 0 prior probability. 
  Must be a positive integer or \code{NA}.}
}
\value{
  Returns an object of class "model.prior", with the following components.
  \item{family}{Value is the string \code{"matryoshka-doll"}.}
  \item{hyper.parameters}{Value is a list with components \code{prob_ratio} and \code{type}.}
  \item{trunc}{Value is the user supplied value of \code{trunc}.}
}
\description{
  Creates an object representing the Matryoshka Doll prior distribution on regresion models.
}
\details{
  The Matryoshka Doll prior fixes the ratio of probabilities of a model against a set models that nest it. 
  This nesting set is either all nesting models or only those nesting models that contain a single additional parameter.
  The Matryoshka Doll achieves a Bayesian version of multiplicity correction by treating each model as
  a local null hypothesis with set of alternatives that are the specified set of models that nest it.
  The prior ensures that the probabilities of each null and its corresponding set of alternatives are comparable. 
  
  If  \code{type = "descendants"}  and \eqn{\eta=} \code{prob_ratio} then
  
  \deqn{
    P(M) = \eta \times \sum_{M^\prime \supsetneq M} P(M^\prime)
  }
  
  When the number of predictors increases to \code{Inf}, then the prior using \code{type=descendants} converges to 
  Poisson distributions on model complexity with rate given by \code{log(1+1/prob_ratio)}. 
  
  If \code{type = "children"} and \eqn{\omega=} \code{prob_ratio} then
  
  \deqn{
    P(M) = \omega \times \sum_{\substack{M^\prime \supsetneq M\\|M^\prime| = |M|+1}} P(M^\prime).
  }
  
  The prior using \code{type="children"} is a truncated Poisson distribution on model complexity with rate given by \code{1/prob_ratio}.
  
  
 
  If \code{trunc} \eqn{<p}, then models with \eqn{k>} \code{trunc} predictors get 
  prior probability 0 and the prior is renormalized. 
  If \code{trunc} \eqn{=} \code{NA}, then \code{trunc} is set to \eqn{p}. 
}
\seealso{
  \code{\link{lm.mcmc}()}, \cr
  \code{\link{model.prior-class}},\cr
  \code{\link{binomial.prior}()}, \cr
  \code{\link{binomial_complexity.prior}()}, \cr
  \code{\link{beta_binomial.prior}()}, \cr
  \code{\link{beta_binomial_complexity.prior}()}, \cr
  \code{\link{negative_binomial.prior}()}, \cr
  \code{\link{negative_binomial_complexity.prior}()}, \cr
  \code{\link{custom.prior}()}.
}
\author{
  Andrew Womack
}
