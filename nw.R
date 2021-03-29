# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# install and load packages
libraries = c("locpol", "dplyr")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)


# specify the model
set.seed(20130318)
n = 100  # number of observations
x = runif(n)  # uniform sample
sigma = sqrt(0.1)
eps = rnorm(n) * sigma  # error
f = (sin(2*pi*(x^3)))^3 %>% unname() # true line
y = f + eps

# NW from scratch
mNW <- function(x, X, Y, h = bandwidth, K = dnorm) {
  Kx <- sapply(X, function(Xi) K((x - Xi) / h) / h)
  W <- Kx / rowSums(Kx) 
  drop(W %*% Y)
}

# custom
xGrid <- seq(0, 1, l=100)

# the bandwidth in the native implementation corresponds to the double the interquartile distance (see ksmooth docs)
# the bandwidth in the manual implementation corresponds to the std. deviation of the kernel
bandwidth = 0.04
fct = 2 * (qnorm(0.75) - qnorm(0.25))

plot(x, y, col = "red", lwd = 3, ylab = "", 
     xlab = "x", cex.lab = 2, cex.axis = 2, ylim = range(f))
lines(ksmooth(x, y, bandwidth = fct * bandwidth, kernel = "normal", x.points = xGrid), lwd = 3, col= "blue")

mNW_results = mNW(xGrid, x, y)
lines(xGrid, mNW_results)
