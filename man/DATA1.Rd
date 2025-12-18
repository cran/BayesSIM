\name{DATA1}
\alias{DATA1}
\docType{data}
\title{Simulated Single–Index Data (n = 200, p = 4)}
\usage{data(DATA1)}
\format{
\describe{
  \item{X}{Numeric matrix of size 200 \eqn{\times} 4, entries i.i.d. \eqn{Unif(-1, 1)}}.
  \item{y}{Numeric vector of length 200.}
}
}
\description{
Synthetic data from a single–index model \eqn{y = f(X'\theta) + \varepsilon}
with \eqn{f(u) = u^2 \exp(u)} and \eqn{\varepsilon \sim N(0,\sigma^2)}.
The index vector is \eqn{\theta = (2,1,1,1) / \| (2,1,1,1) \|_2} and
\eqn{\sigma = 0.5}.
}
\examples{
data(DATA1)
str(DATA1)
}
\keyword{datasets}
