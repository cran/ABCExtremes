\name{summarize.stat}
\alias{summarize.stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summarize.stat
}
\description{
This function produces the summary statistic used in ABC rejection sampling.  It allows for spatial extremes data without unit Frechet margins.  It calls the >summary() function.
}
\usage{
summarize.stat(data, assign, frechet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Spatial extremes data in matrix form.  Each row is an independent replicate, and each column is one location.  The order of the columns must match the order of the rows in coord.  Need not have unit Frechet margins.
}
  \item{assign}{
The vector assignments of the triplets into one of the k groups.  This comes from the >cluster() function.
}
  \item{frechet}{
A logical indicator if the data already has unit Frechet margins or not.  If not, data are first transformed to unit Frechet.
}
}
\value{
Produces a k dimensional vector whose elements are the average extremal tripletwise coefficients in each group.
}
\examples{
## Simulate 100 years of a max-stable process at 25 randomly placed locations on a 10 by 10 square.
## Notice that the true parameter is (range,smooth) = (1,3)
library(SpatialExtremes)
n=10
a=runif(n,0,10)
b=runif(n,0,10)
coord=cbind(a,b)
yr=100
data=rmaxstab(yr, coord, cov.mod = "whitmat", nugget=1, range=3, smooth=1, grid=FALSE)
plot(coord)

## cluster the triplets into k groups.  The command will write the number of clusters that change 
## assignment each iteration, and the method stops when there are 0 changes.
k=5
assign=cluster(coord, k, "kmeans")

## Form an summary statistic for the observed data
sumstat=summarize.stat(data, assign, frechet="F")
sumstat
}
