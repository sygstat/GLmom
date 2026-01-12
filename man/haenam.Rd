\name{haenam}
\alias{haenam}
\docType{data}
\title{Haenam Maximum Rainfall Data}
\description{
Annual maximum daily rainfall data from Haenam, South Korea.
This dataset is used for demonstrating model averaging methods
for high quantile estimation in extreme value analysis.
}
\usage{data(haenam)}
\format{
A data frame containing annual maximum daily rainfall values.
}
\source{
Korea Meteorological Administration (KMA; \url{https://www.kma.go.kr})
}
\references{
Shin, Y., Shin, Y. & Park, J.-S. (2025). Model averaging with mixed criteria
for estimating high quantiles of extreme values.
arXiv preprint arXiv:2505.21417. \doi{10.48550/arXiv.2505.21417}
}
\examples{
data(haenam)
head(haenam)

# Estimate high quantiles using model averaging
\donttest{
result <- ma.gev(haenam[,1], quant = c(0.98, 0.99, 0.995))
print(result$zp.ma)
}
}
\keyword{datasets}
