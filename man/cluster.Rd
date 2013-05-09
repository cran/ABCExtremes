\name{cluster}
\alias{cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cluster
}
\description{
This function geometrically clusters the (n choose 3) triplets into k groups using either Ward's method, or the k-means++ algorithm}
\usage{
cluster(coord, k, method = c("Ward", "kmeans"))
}
\arguments{
  \item{coord}{
The two dimensional coordinates of the observed data.  Rows are locations, columns are coordinates.
}
  \item{k}{
The desired number of clusters k (same as the dimension of the summary statistic)
}
  \item{method}{
Choose either "Ward" or "kmeans" for Ward's and kmeans clustering, respectively.  K-means++ may run faster for 
datasets with a large number of locations (n).
}
}
\value{
Returns a vector of dimention (n choose 3), where n is the number of locations of observed data.  
Each element is the cluster assignment from 1 to k.
}

\examples{

library(SpatialExtremes)
n=10
a=runif(n,0,10)
b=runif(n,0,10)
coord=cbind(a,b)
k=5
yr=30
data=rmaxstab(yr, coord, cov.mod = "whitmat", nugget=1, range=3, smooth=1, grid=FALSE)

plot(coord)

out=cluster(coord, k, "kmeans") 
out
}

