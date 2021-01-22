#bernoulli function taking data vector and return estimated param phat
bernoulli <- function(data) {
  phat = mean(data)
  phat
}

#bernoulli maximum likelihood estimator function taking data vector and return estimated param phat
bernoulli_mle <- function(data) {
  phat = mean(data)
  phat
}

#likelihood bernoulli
#prior beta
#posterior beta
bernoulli_bayes <- function(data,
                            priora = 1,
                            priorb = 1) {
  s = sum(data)
  n = binomial(data)[2]
  print("prior beta, posterior beta with the params: ")
  posta = priora + s
  postb = priorb + n - s
  list(posta = posta, postb = postb)
}

#binomial function taking data vector and return estimated param phat, nhat
binomial <- function(data) {
  svar = var(data)
  smean = mean(data)
  nhat = (smean ^ 2) / (smean - svar)
  phat = smean / nhat
  c(nhat, phat)
}

#binomial maximum likelihood estimator function taking data vector and return estimated param phat, nhat
binomial_mle <- function(data) {
  svar = var(data)
  smean = mean(data)
  nhat = (smean ^ 2) / (smean - svar)
  phat = smean / nhat
  c(nhat, phat)
}

#likelihood binomial
#prior beta
#posterior beta
binomial_bayes <- function(data,
                           priora = 1,
                           priorb = 1) {
  s = sum(data)
  n = binomial(data)[2]
  print("prior beta, posterior beta with the params: ")
  posta = priora + s
  postb = priorb + n - s
  print(list(posta = posta, postb = postb))
}

#geometric function taking data vector and return estimated param phat
geometric <- function(data) {
  smean = mean(data)
  phat = 1 / smean
  phat
}

#geometric mle function taking data vector and return estimated param phat
geometric_mle <- function(data) {
  smean = mean(data)
  phat = 1 / smean
  phat
}

#likelihood geometric
#prior beta
#posterior beta
geometric_bayes <- function(data,
                            priora = 1,
                            priorb = 1) {
  phat = geometric(data)
  s = sum(data)
  print("prior beta, posterior beta with the params: ")
  postalpha = priora + s
  postbeta  = priorb + 1
  list(postalpha = postalpha, postbeta = postbeta)
}

#poisson function taking data vector and return estimated param lambdahat
poisson <- function(data) {
  smean = mean(data)
  lambdahat = smean
  lambdahat
}
#poisson mle function taking data vector and return estimated param lambdahat
poisson_mle <- function(data) {
  smean = mean(data)
  lambdahat = smean
  lambdahat
}

#likelihood poisson
#prior gamma
#posterior negative binomial
poisson_bayes <- function(data,
                          priora = 1,
                          priorb = 1) {
  s = sum(data)
  lambda = poisson(data)
  n = length(data)
  print("prior gamma, posterior negative binomial with the params: ")
  postalpha = priora + s
  postbeta  = priorb + n
  
  list(postalpha = postalpha, postbeta = postbeta)
}

#uniform function taking data vector and return estimated param alphahat, betahat
uniform <- function(data) {
  smean = mean(data)
  ssd = sd(data)
  alphahat = smean - (sqrt(3) * ssd)
  betahat = smean + (sqrt(3) * ssd)
  c(alphahat, betahat)
}

#uniform function taking data vector and return estimated param alphahat, betahat
uniform_mle <- function(data) {
  theta = max(data)
  theta
}

uniform_bayes <- function(data,
                          priora = 2.5,
                          priorb = 4) {
  n = length(data)
  maxdata = max(data)
  
  print("prior pareto, posterior pareto with the params: ")
  
  posta = n + priora
  postb = max(maxdata, priorb)
  c(posta, postb)
}

#normal function taking data vector and return estimated param sample mean, variance
normal <- function(data) {
  smean = mean(data)
  svar = var(data)
  c(smean, svar)
}

#normal function taking data vector and return estimated param sample mean, variance
normal_mle <- function(data) {
  smean = mean(data)
  svar = var(data)
  c(smean, svar)
}


#likelihood normal
#prior normal
#posterior normal
normal_bayes <- function(data,
                         priormean = 4,
                         priorvar = 2) {
  mme <- normal(data)
  lmean = mme[1]
  lvar = mme[2]
  print("prior normal, posterior normal with the params: ")
  n = length(data)
  a = 1 / (priorvar ^ 2)
  b = n / (lvar ^ 2)
  
  postmean = ((a * priormean) + (b * lmean)) / (a + b)
  postvar = 1 / (a + b)
  
  list(postmean = postmean, postvar = postvar)
}

exponential <- function(data) {
  smean = mean(data)
  lambdahat = 1 / smean
  lambdahat
}

#exponential mle function taking data vector and return estimated param lambdahat
exponential_mle <- function(data) {
  smean = mean(data)
  lambdahat = 1 / smean
  lambdahat
}

#likelihood exponential
#prior gamma
#posterior Lomax distribution
exponential_bayes <- function(data,
                              priora = 1,
                              priorb = 1) {
  lambdahat = exponential(data)
  n = length(data)
  print("prior gamma, posterior Lomax distribution with the params: ")
  #posterior Lomax distribution
  alphahat = priora + n
  betahat =  +sum(data)
  
  list(alphahat = alphahat, betahat = betahat)
}

gamma <- function(data) {
  smean = mean(data)
  svar = var(data)
  alphahat = (smean) ^ 2 / svar
  betahat = svar / smean
  c(alphahat, betahat)
}

#gamma mle function taking data vector and return estimated params
gamma_mle <- function(data) {
  n <- length(data)
  mean_data <- mean(data)
  
  mme <- gamma(data)
  alpha_init <- mme[1]
  beta_init  <- mme[2]
  
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

#likelihood gamma
#prior gamma
#posterior compound gamma distribution
gamma_bayes <- function(data) {
  mme <- gamma(data)
  laplha = mme[1]
  lbeta = mme[2]
  n = length(data)
  
  #prior gamma
  a = 1
  b = 1
  print("prior gamma, posterior compound gamma distribution with the params: ")
  #posterior compound gamma distribution
  alphahat = a + (n * laplha)
  betahat =  b + sum(data)
  
  list(alphahat = alphahat, betahat = betahat)
}

beta <- function(data) {
  smean = mean(data)
  svar = var(data)
  param = (((smean * (1 - smean)) / (svar)) - 1)
  betahat = (1 - smean) * param
  alphahat = smean * param
  c(alphahat, betahat)
}

#beta mle function taking data vector and return estimated params
beta_mle <- function(data) {
  n <- length(data)
  # alpha and beta hat params initialized with method of moment params
  mme <- beta(data)
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

beta_bayes <- function(data) {
  "No formal posterior"
}

chisquare <- function(data) {
  smean = mean(data)
  betahat = smean
  betahat
}

#chisquare mle function taking data vector and return estimated params
chisquare_mle <- function(data) {
  smean = mean(data)
  betahat = smean
  betahat
}

# special case of gamma with beta=n/2, alpha=2
#likelihood gamma
#prior gamma
#posterior compound gamma distribution
chisquare_bayes <- function(data) {
  betahat = chisquare(data)
  n = length(data)
  
  #prior gamma
  a = 1
  b = 1
  print("prior gamma, posterior compound gamma distribution with the params: ")
  #posterior compound gamma distribution
  alphahat = a + (n * 2)
  betahat =  b + sum(data)
  
  list(alphahat = alphahat, betahat = betahat)
}

mmestimator <- function(distribution, data) {
  switch(
    distribution,
    "bernoulli" = {
      print("Bernoulli Distribution: ")
      print("mme: ")
      print(bernoulli(data))
      print("mle: ")
      print(bernoulli_mle(data))
      
      bernoulli_bayes(data)
    },
    "binomial"  = {
      print("Binomial Distribution: ")
      print("mme: ")
      print(binomial(data))
      print("mle: ")
      print(binomial_mle(data))
      binomial_bayes(data)
    },
    "geometric" = {
      print("Geometric Distribution:")
      print("mme: ")
      print(geometric(data))
      print("mle: ")
      print(geometric_mle(data))
      geometric_bayes(data)
    } ,
    "poisson"   = {
      print("Poisson Distribution:")
      print("mme: ")
      print(poisson(data))
      print("mle: ")
      print(poisson_mle(data))
      poisson_bayes(data)
    },
    "uniform"   = {
      print("Uniform Distribution:")
      print("mme: ")
      print(uniform(data))
      print("mle: ")
      print(uniform_mle(data))
      uniform_bayes(data)
    },
    "normal"    = {
      print("Normal Distribution:")
      print("mme: ")
      print(normal(data))
      print("mle: ")
      print(normal_mle(data))
      normal_bayes(data)
    },
    "exponential" = {
      print("Exponential Distribution:")
      print("mme: ")
      print(exponential(data))
      print("mle: ")
      print(exponential_mle(data))
      exponential_bayes(data)
    },
    "gamma"     = {
      print("Gamma Distribution:")
      print("mme: ")
      print(gamma(data))
      print("mle: ")
      print(gamma_mle(data))
      gamma_bayes(data)
    },
    "beta"      = {
      print("Beta Distribution:")
      print("mme: ")
      print(beta(data))
      print("mle: ")
      print(beta_mle(data))
      beta_bayes(data)
    },
    "chisquare" = {
      print("Chisquare Distribution:")
      print("mme: ")
      print(chisquare(data))
      print("mle: ")
      print(chisquare_mle(data))
      chisquare_bayes(data)
    },
    print("Distribution doesn't exist")
  )
}

#testing the functions area

# vector for plots
x <- seq(0, 7, length = 100)

# generate uniform data using -1, 1 and compare the results
r1 <- runif(100,-1, 1)
print(mmestimator("uniform", r1))
print("_________________________________________")


# generate binomial data and compare the results
r2 <- rbinom(8, 150, .4)
params <- mmestimator("binomial", r2)
print("_________________________________________")

# generate geometric data and compare the results
r3 <- rgeom(1000, 1 / 2)
print(mmestimator("geometric", r3))
print("_________________________________________")

# generate poisson data and compare the results
r4 <- rpois(100, 10)
print(mmestimator("poisson", r4))
print("_________________________________________")

# generate normal data and compare the results
r6 <-
  rnorm(100, 100, sd = 15) #rnorm takes the sd so the variance will be sd^2
print(mmestimator("normal", r6))
print("_________________________________________")

# generate exponential data and compare the results
r7 <- rexp(100, 1 / 3)
print(mmestimator("exponential", r7))
print("_________________________________________")

# generate gamma data and compare the results
r8 <- rgamma (n = 100 , shape = 9 , scale = .5)
print(mmestimator("gamma", r8))
print("_________________________________________")

# generate chisquare data and compare the results
r9 <- rchisq (n = 100 , 2)
l <- mmestimator("chisquare", r9)
print(l)
print("_________________________________________")

# generate beta data and compare the results, and plot the density for both to show difference
r10 <- rbeta (n = 100 , 2 , 5)
print(mmestimator("beta", r10))
print("_________________________________________")

# hx <- dnorm(x)
# colors <- c("red", "blue")
# labels <- c("df=original", "df=estimated")
# plot(x, hx, type="l", lty=2, xlab="x value",
#      ylab="Density", main="Comparison of t Distributions")
# lines(x, dbeta (x , 2 , 5), lwd=2, col=colors[1])
# lines(x, dbeta (x , l[1] , l[2]), lwd=2, col=colors[2])
# legend("topright", inset=.05, title="Distributions",
#        labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)


#____________________________________________________
