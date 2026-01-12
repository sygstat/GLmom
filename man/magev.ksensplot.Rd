\name{magev.ksensplot}
\alias{magev.ksensplot}
\title{K Sensitivity Plot for MAGEV}
\description{
Plots return level estimates, standard errors, and first-order differences
across different numbers of candidate submodels K. This helps identify
stable regions where estimates converge and select an optimal K value.
}
\usage{
magev.ksensplot(data = NULL, q.cut = 0.6, mink = 4, maxk = 20,
                quant = c(0.99, 0.995))
}
\arguments{
  \item{data}{A numeric vector of data to be fitted (e.g., annual maxima).}
  \item{q.cut}{Quantile cutoff for determining stability (default 0.6).}
  \item{mink}{Minimum number of candidate submodels to test (default 4).}
  \item{maxk}{Maximum number of candidate submodels to test (default 20).}
  \item{quant}{The probabilities for high quantile estimation. Default is c(0.99, 0.995).}
}
\details{
The function computes MAGEV estimates for K ranging from \code{mink} to
\code{maxk}. For each K, it calculates return level estimates (black points),
normalized standard errors (blue line), and first-order differences (red triangles).

The optimal K is selected as the smallest value where both the normalized
standard error and first-order difference are below their respective
\code{q.cut} quantile cutoffs. The selected K is indicated by a purple
vertical line.
}
\value{
The optimal K value (integer) selected by the algorithm.
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
optimal_k <- magev.ksensplot(streamflow$r1)
print(optimal_k)
}
}
