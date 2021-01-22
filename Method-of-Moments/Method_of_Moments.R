#bernoulli function taking data vector and return estimated param phat
bernoulli<-function(data) {
  phat=mean(data)
  phat
}
#binomial function taking data vector and return estimated param phat, nhat
binomial<-function(data) {
  svar=var(data)
  smean=mean(data)
  nhat=(smean^2)/(smean-svar)
  phat=smean/nhat
  c(nhat, phat)
}
#geometric function taking data vector and return estimated param phat
geometric<-function(data) {
  smean=mean(data)
  print(smean)
  phat=1/smean
  phat
}

#poisson function taking data vector and return estimated param lambdahat
poisson<-function(data) {
  smean=mean(data)
  lambdahat=smean
  lambdahat
}

#uniform function taking data vector and return estimated param alphahat, betahat
uniform<-function(data) {
  smean=mean(data)
  ssd=sd(data)
  alphahat=smean-(sqrt(3)*ssd)
  betahat=smean+(sqrt(3)*ssd)
  c(alphahat,betahat)
}

#normal function taking data vector and return estimated param sample mean, variance
normal<-function(data) {
  smean=mean(data)
  svar=var(data)
  c(smean,svar)
}
exponential<-function(data) {
  smean=mean(data)
  lambdahat=1/smean
  lambdahat
}
gamma<-function(data) {
  smean=mean(data)
  svar=var(data)
  alphahat=(smean)^2/svar
  betahat=svar/smean
  c(alphahat,betahat)
}
beta<-function(data) {
  smean=mean(data)
  svar=var(data)
  param=(((smean*(1-smean))/(svar))-1)
  betahat=(1-smean)*param
  alphahat=smean*param
  c(alphahat,betahat)
}
chisquare<-function(data) {
  smean=mean(data)
  betahat=smean
  betahat
}
mmestimator<-function(distribution, data){
  switch(distribution,
         "bernoulli" = bernoulli(data),
         "binomial"  = binomial(data),
         "geometric" = geometric(data),
         "poisson"   = poisson(data),
         "uniform"   = uniform(data),
         "normal"    = normal(data),
         "exponential" = exponential(data),
         "gamma"     = gamma(data),
         "beta"      = beta(data),
         "chisquare" = chisquare(data),
         print("Distribution doesn't exist") 
  )
}

#testing the functions area 

# vector for plots 
x <- seq(0, 7, length=100)

# generate uniform data using -1, 1 and compare the results
r1 <- runif(100, -1, 1)
print(mmestimator("uniform", r1))


# generate binomial data and compare the results
r2 <- rbinom(8,150,.4)
print(mmestimator("binomial", r2))

# generate geometric data and compare the results
r3 <- rgeom(1000, 1/2)
print(mmestimator("geometric", r3))

# generate poisson data and compare the results
r4 <- rpois(100, 10)
print(mmestimator("poisson", r4))

# generate normal data and compare the results
r6 <- rnorm(100, 100, sd=15) #rnorm takes the sd so the variance will be sd^2
print(mmestimator("normal", r6))

# generate exponential data and compare the results
r7 <- rexp(100, 1/3)
print(mmestimator("exponential", r7))

# generate gamma data and compare the results
r8 <- rgamma (n =100 , shape =9 , scale =.5)
print(mmestimator("gamma", r8))

# generate chisquare data and compare the results
r9 <- rchisq (n =100 , 2)
l <- mmestimator("chisquare", r9)
print(l)

# generate beta data and compare the results, and plot the density for both to show difference
r10 <- rbeta (n =100 , 2 , 5)
l <- mmestimator("beta", r10)
print(l)
hx <- dnorm(x)
colors <- c("red", "blue")
labels <- c("df=original", "df=estimated")
plot(x, hx, type="l", lty=2, xlab="x value",
     ylab="Density", main="Comparison of t Distributions")
lines(x, dbeta (x , 2 , 5), lwd=2, col=colors[1])
lines(x, dbeta (x , l[1] , l[2]), lwd=2, col=colors[2])
legend("topright", inset=.05, title="Distributions",
       labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)