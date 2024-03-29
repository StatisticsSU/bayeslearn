% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{rDirichlet}
\alias{rDirichlet}
\title{Simulates from the Dirichlet distribution.}
\usage{
rDirichlet(n, alpha)
}
\arguments{
\item{n}{number of draws}

\item{alpha}{vector with hyperparameters}
}
\value{
n-by-p matrix where each row is a draw.
}
\description{
Simulates from the Dirichlet distribution.
}
\examples{
library(bayeslearn)
draws = rDirichlet(n = 1000, alpha = c(3,5,10))
colMeans(draws) # Close to theoretical c(3,5,10)/18 = 0.167 0.278 0.556
}
