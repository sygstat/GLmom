\name{bangkok}
\alias{bangkok}
\docType{data}
\title{Bangkok Maximum Rainfall Data}
\description{
Annual maximum daily rainfall data from Bangkok, Thailand.
This dataset is used for demonstrating model averaging methods
for high quantile estimation in extreme value analysis.
}
\usage{data(bangkok)}
\format{
A data frame containing annual maximum daily rainfall values.
}
\source{
Thai Meteorological Department (TMD; \url{https://www.tmd.go.th})
}
\references{
Shin, Y., Shin, Y. & Park, J.-S. (2025). Model averaging with mixed criteria
for estimating high quantiles of extreme values.
arXiv preprint arXiv:2505.21417. \doi{10.48550/arXiv.2505.21417}
}
\examples{
data(bangkok)
head(bangkok)

# Estimate high quantiles using model averaging
\donttest{
result <- ma.gev(bangkok[,1], quant = c(0.99, 0.995))
print(result$zp.ma)
}
}
\keyword{datasets}
