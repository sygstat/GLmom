\name{magev.qqplot}
\alias{magev.qqplot}
\title{Q-Q Diagnostic Plot for MAGEV}
\description{
Creates a 2x2 panel of Q-Q plots comparing observed vs. fitted quantiles
for different estimation methods: MLE, LME, surrogate (MA), and REMLE.
}
\usage{
magev.qqplot(data = NULL, zx = NULL)
}
\arguments{
  \item{data}{A numeric vector of observed data.}
  \item{zx}{A list object returned by \code{\link{ma.gev}} with \code{remle = TRUE}.}
}
\details{
The function creates four Q-Q plots:
\itemize{
  \item Upper left: MLE (Maximum Likelihood Estimation)
  \item Upper right: LME (L-moment Estimation)
  \item Lower left: Surrogate MA (Model Averaging surrogate)
  \item Lower right: REMLE (Restricted MLE, if available)
}
Points close to the 45-degree diagonal line indicate good model fit.
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
qq <- c(seq(0.01, 0.99, by = 0.01), 0.995, 0.999)
zx <- ma.gev(streamflow$r1, quant = qq, weight = 'like1',
             numk = 9, varcom = FALSE, remle = TRUE)
magev.qqplot(data = streamflow$r1, zx = zx)
}
}
