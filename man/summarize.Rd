\name{summarize}
\alias{summarize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summarize
}
\description{
This function produces the summary statistic from unit Frechet data and assigned groups.  It is called by >summary.stat(), which allows for transforming the margins.  Use >summary.stat() to produce a summary statistic.
}
\usage{
summarize(datat, assign)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datat}{
A matrix of spatial extremes data transformed to unit Frechet margins.  Rows are replicates, columns are locations.  See the SpatialExtremes package for more information on transforming all margins to unit Frechet.
}
  \item{assign}{
A vector of assignments of each of the (n choose 3) triplets into one of the k groups.
}
}
\value{
A vector of dimension k, where each element is the average extremal tripletwise coefficient in the kth cluster.
}

