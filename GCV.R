# clear variables
rm(list = ls(all = TRUE))

# install and load packages
libraries = c("locpol", "dplyr")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x) })
lapply(libraries, library, quietly = TRUE, character.only = TRUE)


# specify the model
set.seed(20130318)
n = 100  # number of observations
n_MC = 1000 # number of Monte-Carlo iterations

sim <- function() {
  x = runif(n)  # uniform sample
  sigma = sqrt(0.1)
  eps = rnorm(n) * sigma  # error
  f = (sin(2*pi*(x^3)))^3 %>% unname() # true line
  y = f + eps
  list(x = x, y = y, f = f)
}  

sim_one = sim()

xyf = cbind(sim_one$x, sim_one$y, sim_one$f)
xyf = xyf[order(xyf[, 1]), ]
x = xyf[, 1] 
y = xyf[, 2]
f = xyf[, 3]

# bandwidth grid
n_h = 25
h = seq(0.03, by = (0.15-0.03)/(n_h-1), length = n_h)

# Nadaraya-Watson estimator
fNW <- function(x, X, Y, h, K = dnorm) {
  Kx <- sapply(X, function(Xi) K((x - Xi) / h) / h)
  if (is.vector(Kx)) Kx = matrix(Kx, nrow = 1)
  W <- Kx / rowSums(Kx) 
  drop(W %*% Y)
}  
  
# compute empirical, LOO CV and GCV errors, variance and bias^2 
# for given data sample and a range of bandwidths
L_CV = matrix(0, n_h, 1)
L_GCV = L_CV
L_emp = L_CV
bias = L_CV
var = L_CV

for (k in 1:n_h) {
  # Nadaraya–Watson with Gaussian kernel 
  fh = fNW(x = x, X = x, Y = y, h = h[k])
  # Monte-Carlo estimates of true bias and variance
  biases = matrix(0, n_MC, 1)
  vars = biases
  for (i in 1:n_MC) {
    sim_MC = sim()
    fh_MC = fNW(x = sim_MC$x, X = sim_MC$x, Y = sim_MC$y, h = h[k])
    biases[i] = sum(sim_MC$y-fh_MC)/n
    vars[i] = sum((fh_MC-sum(fh_MC)/n)^2)/n
    }
  bias[k] = sum(biases)/n_MC 
  var[k] = sum(vars)/n_MC
  # empirical error
  L_emp[k] = sum((y-fh)^2)/n 
  # LOO CV
  fh_cv = sapply(1:n, function(i) 
    fNW(x = x[i], X = x[-i], Y = y[-i], h = h[k]))
  L_CV[k] = sum((y-fh_cv)^2)/n 
  # GCV
  tr_est = dnorm(0)/h[k]
  L_GCV[k] = 1/(1-tr_est/n)^2 * L_emp[k]
}


# plot
png("errors.png", width = 900, height = 900, bg = "transparent")
plot(h, L_CV, type = "l", lwd = 3, lty = 2, col = "black", xlab = "Bandwidth h", 
     ylab = "", cex.lab = 2, cex.axis = 2, ylim = c(min(L_emp), max(L_CV)))
title("Choosing optimal bandwidth", cex.main = 2)
lines(h, L_emp, lwd = 3, lty = 1, col = "blue3")
abline(h = sigma^2, col = "green", lwd = 1, lty = 1)
lines(h, L_GCV, lwd = 3, lty = 2, col = "red3")
legend("bottomright", c("L_emp", "L_CV", "L_GCV", "sigma^2"), lty = c(1, 2, 2, 1), 
       lwd = c(3, 3, 3, 1), col = c("blue3", "black", "red3", "green"), cex = 1.5)
dev.off()

png("bv.png", width = 900, height = 900, bg = "transparent")
par(mar = c(5, 4, 4, 4) + 0.3)
plot(h, var, type = "l", lwd = 3, col = "red3", xlab = "Bandwidth h", ylab = "Variance", 
     cex.lab = 2, cex.axis = 2, ylim = c(0, max(var)))
par(new = TRUE)
plot(h, bias^2, type = "l", lwd = 3, axes = FALSE, col = "blue3", ylab = "", xlab = "")
axis(side = 4, at = pretty(range(bias^2)), cex.axis = 2)
mtext("Bias^2", side = 4, line = 3, cex = 2)
title("Bias-Variance Tradeoff", cex.main = 2)
dev.off()

# choose optimal h acc CV
hopt_cv = h[which(L_CV == min(L_CV))]
nw_opt_cv = ksmooth(x, y, kernel = "normal", bandwidth = hopt_cv)

png("cv.png", width = 900, height = 900, bg = "transparent")
plot(x, f, type = "l",col = "blue3", lwd = 3, ylab = "", 
     xlab = "x", cex.lab = 2, cex.axis = 2, ylim = range(f))
title("Simulated Data Estimated with CV", cex.main = 2)
points(x, y, pch = 19, col = "red3", cex = 0.7)
lines(nw_opt_cv, lwd = 3)
dev.off()

# choose optimal h acc GCV
hopt_gcv = h[which(L_GCV == min(L_GCV))]
nw_opt_gcv = ksmooth(x, y, kernel = "normal", bandwidth = hopt_gcv)

png("gcv.png", width = 900, height = 900, bg = "transparent")
plot(x, f, type = "l", col = "blue3", lwd = 3, ylab = "", 
     xlab = "x", cex.lab = 2, cex.axis = 2, ylim = range(f))
title("Simulated Data Estimated with GCV", cex.main = 2)
points(x, y, pch = 19, col = "red3", cex = 0.7)
lines(nw_opt_gcv, lwd = 3)
dev.off()
