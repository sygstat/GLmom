\name{magev.rlplot}
\alias{magev.rlplot}
\title{Return Level Plot for MAGEV}
\description{
Displays fitted return levels with 95\% confidence intervals against
return period on a log scale.
}
\usage{
magev.rlplot(par = NULL, se.vec = NULL, data = NULL)
}
\arguments{
  \item{par}{A numeric vector of GEV parameters (mu, sigma, xi) in Hosking style.
    Typically from \code{zx$surr$par} where \code{zx} is the output of \code{\link{ma.gev}}.}
  \item{se.vec}{A numeric vector of standard errors for the quantile estimates
    corresponding to the plotting positions. Typically from \code{zx$adj.se.ma}.}
  \item{data}{A numeric vector of observed data (annual maxima).}
}
\details{
The plot shows:
\itemize{
  \item Black line: Fitted return level curve
  \item Blue lines: 95\% confidence interval (mean +/- 1.96 * SE)
  \item Black points: Observed data at empirical return periods
}
The x-axis (return period) is on a log scale, ranging from 0.1 to 900 years.
}
\value{
NULL. The function produces a plot as a side effect.
}
\references{
Shin, Y., Shin, Y., & Park, J.-S. (2025). Model averaging with mixed criteria
for estimating high quantiles of extreme values.
arXiv preprint arXiv:2505.21417. \doi{10.48550/arXiv.2505.21417}
}
\author{
Yonggwan Shin, Yire Shin, Jeong-Soo Park
}
\examples{
\donttest{
data(streamflow)
ff <- c(seq(0.01, 0.09, by = 0.01), seq(0.1, 0.9, by = 0.1),
        0.93, 0.95, 0.98, 0.99, 0.995, 0.999)
zx <- ma.gev(streamflow$r1, quant = ff, weight = 'like1',
             numk = 9, varcom = TRUE)
magev.rlplot(par = zx$surr$par, se.vec = zx$adj.se.ma, data = streamflow$r1)
}
}
