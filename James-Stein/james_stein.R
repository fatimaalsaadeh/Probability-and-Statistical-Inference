library(LaplacesDemon)

#Fatima AlSaadeh(fya7)
#James-stein estimation

JSE <- function(k, n) {
  
  #define mean vector of length k and initialize by 0
  mu = c(1:k) * 0
  #get a multivariate normal matrix of k columns and n rows
  x <- rmvn(n, mu, diag(k))
  
  #find the mean of each column in the matrix
  thetai_mle = colMeans(x)
  
  #find the mle squared error (risk) sum((x-theta)^2)/n
  mle = sum((mu - thetai_mle) ^ 2) / n
  
  #find the James-stein risk 
  js <- ((1 - ((k - 2) / (n * (
    sum(thetai_mle) ^ 2
  )))))
  #get the positive values or replace by 0
  js <-  ifelse(js < 0, 0, js)
  #complete the James-stein risk calculation
  thetai_js <-  js * thetai_mle
  
  #find the jse squared error (risk) sum((thetaj-theta)^2)/n
  risk_js <- sum((mu - thetai_js) ^ 2) / n
  list(mle = mle, risk_js = risk_js)
}


#run for k=100, n=1000
n = 1000
js_risk <- NULL
mle_risk = NULL
for (k in 1:100) {
  est = JSE(k, n)
  mle_risk = c(mle_risk, est$mle)
  js_risk = c(js_risk, est$risk_js)
}
plot(
  c(1:k),
  mle_risk,
  main = "James-stein vs Mle risk estimation",
  ylab = "risk",
  xlab = "k",
  type = "l"
)
legend("topleft",
       c("MLE", "JSE"),
       fill = c("red", "green"))
lines(mle_risk, type = "l", col = "red")
par(new = TRUE)
lines(js_risk, type = "l", col = "green")
