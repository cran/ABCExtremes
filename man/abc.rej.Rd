\name{abc.rej}
\alias{abc.rej}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
abc.rej
}
\description{
This function facilitates the ABC rejection sampling.  The intention is for the user to directly modify this
code to suit his or her purposes.  As presently written, this function will use a uniform(0,10) prior for
the range parameter, the weighted summary statistic from the Erhardt and Smith 2012 CSDA paper, and all 
draws will be independent.
}
\usage{
abc.rej(sum.stat, assign, sims, coord, yr, cov=c("whitmat", "powexp","cauchy")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sum.stat}{
summary statistic as a vector, taken from >summary.stat()
}
  \item{assign}{
The group assignments of each individual triplet.  Output from the >cluster() function
}
  \item{sims}{
The number of abc rejection draws.
}
  \item{coord}{
The 2 dimensional coordinates of the observed data.  Rows are individual locations, columns are coordinates.
}
  \item{yr}{
The number of independent replicates of the observed data (blocks, years, etc.)
}
  \item{cov}{
The covriance model underlying the max-stable process.  Can be "gauss", "whitmat", "powexp", etc. 
(see SpatialExtremes package for more details)
}
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

sum.stat=summarize.stat(data, assign, frechet="F")

## Run ABC rejection sampling

out=abc.rej(sum.stat, assign, sims=10, coord, yr=100, cov="whitmat")
out

## To run more simulations, consider running this command across many nodes of a research computing 
## cluster.  To run ABC for more than one parameter, with different priors, different summary 
## statistics, or any other changes, directly modify the function >abc.rej(), shown below.


####################################################################
## Here is the function >abc.rej() which should be modified directly 
####################################################################

abc.rej=function(sum.stat, assign, sims, coord, yr, cov=c("whitmat", "powexp","cauchy"))  {
  library(SpatialExtremes)
  library(combinat)
  if(missing(sum.stat) || missing(assign) || missing(sims) || missing(coord) || missing(cov)) {
    stop("Need to specify summary, assign, sims, coord, and cov")} else {}
  n=dim(coord)[1]
  library(combinat)
  C=combn(seq(1:n), 3)
  phi=sum.stat[[1]]
  len=sum.stat[[2]]
  wt=sqrt(-1*len*(len<0)+(len*(len>0)))

  running=array()
  for (i in 1:sims) {
    if(floor(i/100)-floor((i-1)/100)>0) {print(i)} else{}
    rastar=runif(1,0,10)
##  smstar=runif(1,0,10)
    smstar=1
    datastar=rmaxstab(yr, coord, cov.mod = cov, nugget=1, range=rastar, smooth=smstar, grid=FALSE)

    datastart=array()
    for (t in 1:n) {
      hold=gev2frech(datastar[,t], loc=gevmle(datastar[,t])[1], scale=gevmle(datastar[,t])[2], 
        shape=gevmle(datastar[,t])[3])
      datastart=cbind(datastart, hold)
    }
    datastart=datastart[,2:(n+1)]

    phistart=array()
    for (i in 1:dim(C)[2]) {
      phistart[i]=yr/sum(pmin(1/datastart[,C[,i][1]], 1/datastart[,C[,i][2]], 
        1/datastart[,C[,i][3]]))
    }

    phistar=array()
    for (q in 1:length(sum.stat[[1]])) {
      phistar[q]=mean(phistart[which(assign==q)])
    }
  dist=sum(wt*abs(phistar-phi))
  running=rbind(running,cbind(rastar, smstar, dist  ))
  } ## closes loop on number of simulations)
  out=running[-1,]
}
}
