#Fatima AlSaadeh
#FDR assignment 8 following the provided pseudo code
FDR <- function(p, Q) {
  #sort the p vector
  sortedp <- sort(p)
  #number of tests we have
  m <- length(p)
  #Plot sorted Pvalues(smallest to largest) vs line Q*c(1:m)/m
  q <- Q*c(1:m)/m
  

  
  #Find P* = P value <line
  p_star <-(sortedp<q)
  
  #Find the largest p* in p
  pmax<-max(sortedp[p_star])
  
  #Every P<=P* is “interesting”
  intrestingp <- sortedp<=pmax
  
  print(paste("We have ",length(sortedp[intrestingp])," interesting Points: "))
  print(sortedp[intrestingp])
  
  plot(
    c(1:m),
    sortedp,
    main = "RDF",
    ylab = "tests",
    xlab = "number of tests",
    type = "l"
  )
  legend("topleft",
         c("Q*c(1:m)/m", "sorted tests", "Interesting tests"),
         fill = c("red", "green", "blue")
  )
  lines(q, type = "l", col = "red")
  par(new = TRUE)
  points(sortedp, col = "green")
  points(c(1:m)[intrestingp], sortedp[intrestingp], col="blue")
}

v<-c(1e-5*runif(100), runif(900))
FDR(v,0.05)

