#Fatima AlSaadeh
#Maximum Likelihood Estimator
#Edited remove an error and added a full ktest function for gamma distribution

#ktest for gamma distribution
ktest_gamma <- function(repeatval = 1000, x) {
  mlegamma <- mlestimator("gamma", x)
  n <- length(x)
  p0 <- ks.test(x, "pgamma", mlegamma[1], mlegamma[2], exact = FALSE)$statistic
  p1 <- NULL
  for (i in 1:repeatval) {
    x_star <- rgamma (n, shape =mlegamma[1] , scale =mlegamma[2])
    mlegamma_star <- mlestimator("gamma", x_star)
    p_star <- ks.test(x, "pgamma", mlegamma_star[1], mlegamma_star[2], exact = FALSE)$statistic
    p1 <- c(p1, p_star)
  }
  p_value <- sum(p1 > p0) /repeatval
  return(p_value)
}
#bernoulli function taking data vector and return estimated param phat
bernoulli_mle <- function(data) {
  phat = mean(data)
  phat
}
#binomial function taking data vector and return estimated param phat, nhat
binomial_mle <- function(data) {
  svar = var(data)
  smean = mean(data)
  nhat = (smean ^ 2) / (smean - svar)
  phat = smean / nhat
  c(nhat, phat)
}
#geometric function taking data vector and return estimated param phat
geometric_mle <- function(data) {
  smean = mean(data)
  phat = 1 / smean
  phat
}
#poisson function taking data vector and return estimated param lambdahat
poisson_mle <- function(data) {
  smean = mean(data)
  lambdahat = smean
  lambdahat
}

#uniform function taking data vector and return estimated param alphahat, betahat
uniform_mle <- function(data) {
  theta = max(data)
  theta
}
#normal function taking data vector and return estimated param sample mean, variance
normal_mle <- function(data) {
  smean = mean(data)
  svar = var(data)
  c(smean, svar)
}

exponential_mle <- function(data) {
  smean = mean(data)
  lambdahat = 1 / smean
  lambdahat
}

beta_mme<-function(data) {
  smean=mean(data)
  svar=var(data)
  param=(((smean*(1-smean))/(svar))-1)
  betahat=(1-smean)*param
  alphahat=smean*param
  c(alphahat,betahat)
}

beta_mle <- function(data) {
  n <- length(data)
  # alpha and beta hat params initialized with method of moment params
  mme <- beta_mme(data)
  alpha_init <- mme[1]
  beta_init  <- mme[2]
  
  #first step estimation of alpha
  deralpha <- (sum( log(data) ) / n) - digamma(alpha_init) + digamma(alpha_init+beta_init)
  #first step estimation of beta
  derbeta <-  (sum( log(1-data) ) / n) - digamma(beta_init) + digamma(alpha_init+beta_init)
  
  #calculate next alpha
  alphahat = alpha_init - deralpha
  betahat = beta_init - derbeta
  c(alphahat, betahat)
}

chisquare_mle <- function(data) {
  smean = mean(data)
  betahat = smean
  betahat
}
gamma_mme <- function(data) {
  smean = mean(data)
  svar = var(data)
  alphahat = (smean) ^ 2 / svar
  betahat = svar / smean
  c(alphahat, betahat)
}
gamma_mle <- function(data) {
  n <- length(data)
  mean_data <- mean(data)
  
  alpha_init <- gamma_mme(data)[1]
  beta_init  <- gamma_mme(data)[2]
  
  # alpha and beta hat params initialized with method of moment params
  alpha_hat <- alpha_init
  beta_hat  <-  beta_init
  
  #first step estimation of alpha
  firststep <- n * log(alpha_init / mean_data) - n * digamma(alpha_init) + sum(log(data))
  
  #calculate next alpha
  alpha_next <- alpha_init - firststep 
  
  alphahat <- alpha_next
  betahat <- mean_data / alpha_next
  
  c(alphahat, betahat)
}
mlestimator <- function(distribution, data) {
  switch(
    distribution,
    "bernoulli" = bernoulli_mle(data),
    "binomial"  = binomial_mle(data),
    "geometric" = geometric_mle(data),
    "poisson"   = poisson_mle(data),
    "uniform"   = uniform_mle(data),
    "normal"    = normal_mle(data),
    "exponential" = exponential_mle(data),
    "gamma"     = gamma_mle(data),
    "beta"      = beta_mle(data),
    "chisquare" = chisquare_mle(data),
    print("Distribution doesn't exist")
  )
}

#testing the functions area

# vector for plots
x <- seq(0, 7, length = 100)

# generate uniform data using -1, 1 and compare the results
r1 <- runif(100,-1, 1)
print("Uniform MLE:")
print(mlestimator("uniform", r1))


# generate binomial data and compare the results
r2 <- rbinom(8, 150, .4)
print("binomial MLE:")
print(mlestimator("binomial", r2))

# generate geometric data and compare the results
r3 <- rgeom(1000, 1 / 2)
print("geometric MLE:")
print(mlestimator("geometric", r3))

# generate poisson data and compare the results
r4 <- rpois(100, 10)
print("poisson MLE:")
print(mlestimator("poisson", r4))

# generate normal data and compare the results
r6 <-
  rnorm(100, 100, sd = 15) #rnorm takes the sd so the variance will be sd^2
print("normal MLE:")
print(mlestimator("normal", r6))

# generate exponential data and compare the results
r7 <- rexp(100, 1 / 3)
print("exponential MLE:")
print(mlestimator("exponential", r7))

# generate gamma data and compare the results
r8 <- rgamma (n =100 , shape =9 , scale =.5)
print("gamma MLE:")
mlegamma <- mlestimator("gamma", r8)
a=mlegamma[1]
b=mlegamma[2]
print(mlegamma)
print("Gamma Ktest:")
print(ktest_gamma(1000, r8))

# generate chisquare data and compare the results
r9 <- rchisq (n = 100 , 2)
print("chisquare MLE:")
l <- mlestimator("chisquare", r9)
print(l)

# generate beta data and compare the results, and plot the density for both to show difference
r10 <- rbeta (n = 100 , 2 , 5)
print("beta MLE:")
l <- mlestimator("beta", r10)
print(l)
hx <- dnorm(x)
colors <- c("red", "blue")
labels <- c("df=original", "df=estimated")
plot(
  x,
  hx,
  type = "l",
  lty = 2,
  xlab = "x value",
  ylab = "Density",
  main = "Comparison of t Distributions"
)
lines(x, dbeta (x , 2 , 5), lwd = 2, col = colors[1])
lines(x, dbeta (x , l[1] , l[2]), lwd = 2, col = colors[2])
legend(
  "topright",
  inset = .05,
  title = "Distributions",
  labels,
  lwd = 2,
  lty = c(1, 1, 1, 1, 2),
  col = colors
)
