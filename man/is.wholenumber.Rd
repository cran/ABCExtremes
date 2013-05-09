\name{is.wholenumber}
\alias{is.wholenumber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
is.wholenumber
}
\description{
An internal function to check if a variable is a positive integer}
\usage{
is.wholenumber(x, tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{x}{
a number
}
  \item{tol}{
tolerance used to check if whole.
}
}
\details{
(There should be no reason to run or modify this function)
}
\examples{
is.wholenumber = function (x, tol = .Machine$double.eps^0.5) 
abs(x - round(x)) < tol
}