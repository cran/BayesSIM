\name{concrete}
\alias{concrete}
\docType{data}
\title{UCI Concrete Compressive Strength (n = 1030, p = 8)}
\usage{data(concrete)}
\format{
Numeric data.frame of size 1030 \eqn{\times} 8.
\describe{
\item{cement}{Numeric. Cement content (kg/m\eqn{^3}).}
  \item{blast_furnace_slag}{Numeric. Blast furnace slag (kg/m\eqn{^3}).}
  \item{fly_ash}{Numeric. Fly ash (kg/m\eqn{^3}).}
  \item{water}{Numeric. Mixing water (kg/m\eqn{^3}).}
  \item{superplasticizer}{Numeric. Superplasticizer (kg/m\eqn{^3}).}
  \item{coarse_aggreate}{Numeric. Coarse aggregate (kg/m\eqn{^3}).}
  \item{fine_aggregate}{Numeric. Fine aggregate (kg/m\eqn{^3}).}
  \item{age}{Numeric. Curing age (days; typically 1--365).}
  \item{strength}{Numeric. Concrete compressive strength (MPa).}
}
}
\description{
Concrete compressive strength dataset from the UCI Machine Learning Repository. No missing variables and there are 8 quantitative inputs and 1 quantitative output.
}
\details{
\strong{Source}
Concrete Compressive Strength in UCI Machine Learning Repository. This data is integrated by experimental data from 17 different sources to check the realiability of the strength. This dataset compiles experimental concrete mixes from multiple studies and is used to predict compressive strength and quantify how mixture ingredients and curing age affect that strength.

\strong{Variables.}
\itemize{
  \item Cement, Blast Furnace Slag, Fly Ash, Water, Superplasticizer, Coarse Aggregate, Fine Aggregate: quantities in kg per m\eqn{^3} of mixture.
  \item Age: curing time in days (1--365).
  \item Target(strength): compressive strength in MPa.
}

}
\references{
Yeh, I. (1998). Concrete Compressive Strength [Dataset]. UCI Machine Learning Repository. https://doi.org/10.24432/C5PK67.

Yeh, I. (1998). Modeling of strength of high-performance concrete using artificial neural networks. Cement and Concrete research, 28(12), 1797-1808.
}
\examples{
data(concrete)
str(concrete)
plot(density(concrete$strength), main = "Concrete compressive strength (MPa)")
}
\keyword{datasets}
