###################################
## Stuff the user would need to run

## library(SpatialExtremes)
## n=25
## a=runif(n,0,10)
## b=runif(n,0,10)
## coord=cbind(a,b)
## k=20
## yr=30
## data=rmaxstab(yr, coord, cov.mod = "whitmat", nugget=1, range=3, smooth=1, grid=FALSE)

#######################################
## Functions begin here:

print("The ABC-REJ algorithm is a highly coputationally expensive method.")
print("This is particularly true for spatial extremes data.")
print("It is extremely likely you will need a high number of simulations,")
print("which is only feasible when abc.rej() is run in parallel on a research computing cluster.")


## 1. This is simply used to ensure k (number of groups) is a positive integer
is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

## 2. This function clusters the triplets into k groups using wither Ward's method, or the k-means++ algorithm
cluster = function(coord, k, method=c("Ward","kmeans")) {
  library(combinat)
  if(missing(method)) method="Ward"
  if(dim(coord)[2]!=2) stop("'coord' must be in 2 dimensions")
  if((is.wholenumber(k)=="FALSE") || (k<0)) stop("'k' must be a positive integer")

  ## enumerates all triplets
  n=dim(coord)[1]
  C=combn(seq(1:n), 3)
  D=matrix(0,length(C[1,]),3)
  for (i in 1:length(C[1,])) {
    x=C[1,i]
    y=C[2,i] 
    z=C[3,i]
    D[i,]=c(dist(rbind(coord[x,],coord[y,],coord[z,])))
  }

  ## sorts the three distances to go from smallest to largest for each row.
  temp=sort(D[1,])
  for (t in 2:length(D[,1])) {
    temp=rbind(temp, sort(D[t,]))
  }

  ## Ward's method, if used:
  if(method=="Ward") {
    W=dist(temp,"manhattan")
    ward=hclust(W, method="ward")
    assign=cutree(ward,k)
  }

  if(method=="kmeans") {

    sam=array()
    objs=length(temp[,1])
    sam[1]=sample(objs,1,replace=FALSE)

    perm2=function(x,y) {
      return(sum(abs(x-y)))
    }

    apply2 = function(x) {
      if (length(x)>3) {temp=apply(x,2,mean)} else {temp=x}
      return(temp)
    }

    ## This first bit initializes the k group centroids, trying to place across the space 
    groups=rep(1,objs)
    for (t in 2:k) {
      sampler=array()
      for (i in 1:objs) {
        sampler[i]=perm2(temp[i,],temp[sam[t-1],])
      }
      sampler[sam]=.0000000001
      pr=sampler^2/sum(sampler^2)
      sam[t]=sample(objs,1,replace=FALSE, prob=pr)
    } # closes t loop
    groups=temp[sam,]

    ## Step 1: Initial assignment.  For each triplet, compute distance to all centroids.  Assign by minimum
    assign=array()
    for (i in 1:objs) {
      ## print(i)
      hold=array()
      for (j in 1:k) {
        hold[j]=perm2(temp[i,],groups[j,])
      }
      temper=min(hold)
      assign[i]=which(hold == temper)
    }

    ## Step 2: Update centroids (always works)
    hold=apply2(temp[which(assign==1),])
    for (i in 2:k) {
      hold=rbind(hold,apply2(temp[which(assign==i),]))
    }
    groups=hold

    len=1
    iter=1
    while(len>0) {
      assignold=assign
      assign=array()
      for (i in 1:objs) {
        hold=array()
        for (j in 1:k) {
          hold[j]=perm2(temp[i,],groups[j,])
        }
        temper=min(hold)
        assign[i]=which(hold == temper)
      }

      len=length(which((assign==assignold)==FALSE))
      print(len)

      ## if there is an empty set, randomly re-assign groups and assign vector
      product=array()
      for (z in 1:length(k)) {
        product[z]=length(temp[which(assign==z),])
      }
      if (prod(product)<1) {
        for (z in 1:length(k)) {
          if(length(temp[which(assign==z),])<1) {groups[z]=sample(objs,1,replace=FALSE)}
        }
        assign=array()
        for (i in 1:objs) {
          hold=array()
          for (j in 1:k) {
            hold[j]=perm2(temp[i,],groups[j,])
          }
          temper=min(hold)
          assign[i]=which(hold == temper)
        }
      print("Fixed Empty Set")
      }
      else {}
      hold=apply2(temp[which(assign==1),])
      for (i in 2:k) {
        hold=rbind(hold,apply2(temp[which(assign==i),]))
      }
      groups=hold
      iter=iter+1
    } #closes while loop
  }
  return(assign)
}

## 3. This function produces the summary statistic from data and assigned groups.
summarize=function(datat, assign) {
  k=max(assign)
  library(combinat)
  yr=dim(data)[1]
  C=combn(seq(1:dim(data)[2]), 3)
  phidata=array()
  for (i in 1:dim(C)[2]) {
    phidata[i]=yr/sum(pmin(1/datat[,C[,i][1]], 1/datat[,C[,i][2]], 1/datat[,C[,i][3]]))
  }

  phi=array()
  for (i in 1:k) {
    phi[i]=mean(phidata[which(assign==i)])
  }

  len=array()
  for (i in 1:k) {
    len[i]=length(which(assign==i))
  }

  listout=list()
  listout[[1]]=phi
  listout[[2]]=len
  return(listout)
}

## 4. This function allows specification of unit-Frechet margins, otherwise it transforms the margins.
summarize.stat = function(data, assign, frechet) {
  n=dim(data)[2]
  if(missing(data)) stop("Data matrix is required")
  if(missing(assign)==F) {
    if(missing(frechet) || frechet!="TRUE") {
      datat=array()
      for (t in 1:n) {
        hold=gev2frech(data[,t], loc=gevmle(data[,t])[1], scale=gevmle(data[,t])[2], shape=gevmle(data[,t])[3])
        datat=cbind(datat, hold)
      }
      datat=datat[,2:(n+1)]
      out=summarize(datat, assign)
    } ## closes 'assign is defined, but data is not Frechet'
    else {
      out=summarize(data, assign)
    } ## closes 'assign is defined, data is Frechet'
  } 
  return(out)
}

## 5. Here's 

abc.rej=function(sum.stat, assign, sims, coord, yr, cov=c("whitmat", "powexp","cauchy"))  {
  library(SpatialExtremes)
  library(combinat)
  if(missing(sum.stat) || missing(assign) || missing(sims) || missing(coord) || missing(cov)) stop("Need to specify summary, assign, sims, coord, and cov") 
  n=dim(coord)[1]
  C=combn(seq(1:n), 3)
  phi=sum.stat[[1]]
  len=sum.stat[[2]]
  wt=sqrt(-1*len*(len<0)+(len*(len>0)))

  running=array()
  for (i in 1:sims) {
    if(floor(i/100)-floor((i-1)/100)>0) {print(i)} else{}
    rastar=runif(1,0,10)
 ## smstar=runif(1,0,5)
    smstar=1
    datastar=rmaxstab(yr, coord, cov.mod = cov, nugget=1, range=rastar, smooth=smstar, grid=FALSE)

    datastart=array()
    for (t in 1:n) {
      hold=gev2frech(datastar[,t], loc=gevmle(datastar[,t])[1], scale=gevmle(datastar[,t])[2], shape=gevmle(datastar[,t])[3])
      datastart=cbind(datastart, hold)
    }
    datastart=datastart[,2:(n+1)]

    phistart=array()
    for (i in 1:dim(C)[2]) {
      phistart[i]=yr/sum(pmin(1/datastart[,C[,i][1]], 1/datastart[,C[,i][2]], 1/datastart[,C[,i][3]]))
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





  




