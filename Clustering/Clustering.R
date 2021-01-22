#Fatima AlSaadeh kmeans assignment
#Kmean multivariate normal distribution with 3 clusters
#Make sure to have kmeans.csv in the directory

iterations <- 20
k <- 3

#First step of assigning each point to the cluster
#with the least squared distance between the mean and the point
my.kmeans <- function(X, k) {
  #take the first mean of a sample of the given points
  means <- X[sample(nrow(X), k),]
  #Initialize an old mean to the current to keep track of the mean values of the samples
  old_means <- means
  
  #Initialize the clusters vector
  cluster <- NULL
  
  #start 20 iterations
  while (iterations > 0) {
    #looping over observations k times
    for (j in k) {
      for (i in 1:nrow(X))
      {
        #initialize minimum distance to large number
        min_dist = 10e9
        #looping over means
        for (m in 1:nrow(means))
        {
          #calculate distance from the point to the mean
          d_from_mean = sum((means[m,] - X[i,]) ^ 2)
          
          #check if the calculated mean is the closest mean to the point
          if (d_from_mean <= min_dist)
          {
            #assign this mean to the point "nearest mean"
            cluster[i] = m
            #change the minimum distance
            min_dist = d_from_mean
          }
        }
      }
    }
    #the points coordinates are the mean of the calculated in the cluster
    means <-
      t(sapply(1:k, function(c)
        apply(X[cluster == c,], 2, mean)))
    
    #update the  iterator index, and old means
    iterations <- iterations - 1
    old_means <- means
  }
  return(list(means = means, cluster = cluster))
}

iterations <- 20
#start kmeans operation on multivariate normal distribution points
my.multivariate.kmeans <- function(X, k) {
  #calculate nearest mean vector for the 3 clusters on the given points
  km.initial <- my.kmeans(X, k)
  meanv <- km.initial$means
  
  #keep tracking of the previous iteration old means
  meanv_old <- meanv
  
  #inititalize the covariance matrix using the initialized points clusters
  cls <-
    sapply(1:k, function(i)
      length(which(km.initial$cluster == i)))
  cls <- cls / sum(cls)
  cov <- array(dim = c(ncol(X), ncol(X), k))
  
  #calculate the first iteration covariance matrix the summation E((x − meanx)(y − meany))
  for (i in 1:ncol(X))
    for (j in 1:ncol(X))
      for (c in 1:k)
        cov[i, j, c] <- 1 / nrow(X) *
    sum((X[km.initial$cluster == c, i] - meanv[c, i]) *
          (X[km.initial$cluster == c, j] - meanv[c, j]))
  
  while (iterations > 0) {
    # calculate the cluster probability from the mean and cov matrix
    multivariate.p <-
      sapply(1:k, function(c)
        pdf.multivariate(X, meanv[c, ], cov[, , c]))
    covi <-
      t(cls * t(multivariate.p)) / rowSums(t(cls * t(multivariate.p)))
    
    # use the calculated culsters for reoptimizing the covariance and mean metrices
    covsum <- colSums(covi)
    cls <- covsum / sum(covsum)
    meanv <-
      t(sapply(1:k, function(c)
        1 / covsum[c] * colSums(covi[, c] * X)))
    
    #calculate the covariance matrix the summation E((x − meanx)(y − meany))
    for (i in 1:ncol(X))
      for (j in 1:ncol(X))
        for (c in 1:k)
          cov[i, j, c] <- 1 / covsum[c] *
      sum(covi[, c] * (X[, i] - meanv[c, i]) *
            covi[, c] * (X[, j] - meanv[c, j]))
    
    #update the stop flag, iterator index, and old mean vector
    iterations <- iterations - 1
    meanv_old <- meanv
  }
  return(list(
    cov = covi,
    m = meanv,
    cluster = apply(covi, 1, which.max)
  ))
}

# calculate covariance matrix for the multivatiate distribution
cov.multivariate <- function(sdv) {
  E <- eigen(sdv)
  Lambda.inv <- diag(E$values ^ -1)
  Q <- E$vectors
  return(Q %*% Lambda.inv %*% t(Q))
}

#multivariate distribution probability density function for xn
pdf.multivariate.n <- function(xn, meanv, sdv) {
  param = 1 / sqrt((2 * pi) ^ length(xn) * det(sdv))
  matrixtranspose = t(xn - meanv)
  matrixmult = matrixtranspose %*% cov.multivariate(sdv) %*% (xn - meanv)
  param * exp(-(1 / 2) * matrixmult)
}

#multivariate distribution probability density function for X
pdf.multivariate <- function(X, meanv, sdv) {
  apply(X, 1, function(xn)
    pdf.multivariate.n(as.numeric(xn), meanv, sdv))
}

#read the data from the file "The file should be in the same directory"
if (!file.exists("./kmeans.csv")) {
  print("Make sure to have kmeans.csv in the directory")
} else {
  datav <- read.csv(file = "./kmeans.csv",
                    header = TRUE,
                    sep = ",")
  
  #prepare the observation vectors
  obs <- data.frame(
    a = datav$a,
    b = datav$b,
    c = datav$c,
    d = datav$d,
    e = datav$e
  )
  
  #pairs the observations vector before the kmeans process
  pairs(
    obs,
    pch = 21,
    lower.panel = NULL,
    col = c('black', 'yellow', 'green', 'red', 'blue')
  )
  
  #run kmeans multivariate normal distribution
  print("Kmean algorithm is running")
  multivariatekmeans <- my.multivariate.kmeans(datav, 3)
  print("The 3 groups covariance matrix")
  print(multivariatekmeans$cov)
  print("The 3 groups mean matrix")
  print(multivariatekmeans$m)
  print("Done!")
  #pairs the observations after the kmeans algorithm with cluster colors
  pairs(datav, lower.panel = NULL, col = multivariatekmeans$cluster)
}