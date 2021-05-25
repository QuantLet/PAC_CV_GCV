[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **PAC_CV_GCV** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet: PAC_CV_GCV

Published in: Metis

Description: 'Fitting the bandwidth of Nadaraya-Watson estimator to simulated data on an interval using cross-validation and generalised cross-validation in R and Python'

Keywords: smoothing, Nadaraya-Watson, cross-validation, generalised cross-validation, empirical error

Author: Anna Shchekina, Lili Matic, Wolfgang Karl Härdle

See also: QID-2142-SPMsimulmase

Submitted: 2021-05-17

```

![Picture1](bv.png)

![Picture2](cv.png)

![Picture3](errors.png)

![Picture4](gcv.png)

### R Code
```r

# clear variables
rm(list = ls(all = TRUE))

# install and load packages
libraries = c("locpol", "dplyr")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x) })
lapply(libraries, library, quietly = TRUE, character.only = TRUE)


# specify the model
set.seed(201053)
n = 100  # number of observations
n_MC = 10 # number of Monte-Carlo iterations
sigma2 = 0.3
sigma = sqrt(sigma2) 
x = runif(n) %>% sort() # uniform sample
f = (sin(2*pi*(x^3)))^3 %>% unname() # true line

y_sim <- function() {
  eps = rnorm(n) * sigma  # error
  y = f + eps
  return(y)
}  

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
bias_MC = L_CV
m1_MC = L_CV
m2_MC = L_CV

y = y_sim()

for (k in 1:n_h) {
  # Nadaraya–Watson with Gaussian kernel 
  fh = fNW(x = x, X = x, Y = y, h = h[k])
  # empirical error
  L_emp[k] = mean((y - fh)^2) 
  # LOO CV
  fh_cv = sapply(1:n, function(i) 
    fNW(x = x[i], X = x[-i], Y = y[-i], h = h[k]))
  L_CV[k] = mean((y - fh_cv)^2)
  # GCV
  tr_est = dnorm(0)/h[k]
  L_GCV[k] = 1/(1 - tr_est/n)^2 * L_emp[k]
}

# Monte-Carlo estimates of true bias and variance
for (j in 1:n_MC) {
  y_MC = y_sim()
  for (k in 1:n_h) {
    fh_MC = fNW(x = x, X = x, Y = y_MC, h = h[k])
    bias_MC[k] = bias_MC[k] + mean(fh_MC - f)
    m1_MC[k] = m1_MC[k] + mean(fh_MC)
    m2_MC[k] = m2_MC[k] + mean(fh_MC^2)
  }
}
bias_MC = bias_MC/n_MC
var_MC = m2_MC/n_MC - (m1_MC/n_MC)^2


# plot
png("errors.png", width = 900, height = 900, bg = "transparent")
plot(h, L_CV, type = "l", lwd = 3, lty = 2, col = "black", xlab = "Bandwidth h", 
     ylab = "", cex.lab = 2, cex.axis = 2, ylim = c(min(L_emp), max(L_CV)))
title("Choosing optimal bandwidth", cex.main = 2)
lines(h, L_emp, lwd = 3, lty = 1, col = "blue3")
abline(h = sigma2, col = "green", lwd = 1, lty = 1)
lines(h, L_GCV, lwd = 3, lty = 2, col = "red3")
legend("bottomright", c("L_emp", "L_CV", "L_GCV", "sigma^2"), lty = c(1, 2, 2, 1), 
       lwd = c(3, 3, 3, 1), col = c("blue3", "black", "red3", "green"), cex = 1.5)
dev.off()

png("bv.png", width = 900, height = 900, bg = "transparent")
par(mar = c(5, 4, 4, 4) + 0.3)
plot(h, var_MC, type = "l", lwd = 3, col = "red3", xlab = "Bandwidth h", ylab = "Variance", 
     cex.lab = 2, cex.axis = 2, ylim = c(0, max(var_MC)))
par(new = TRUE)
plot(h, bias_MC^2, type = "l", lwd = 3, axes = FALSE, col = "blue3", ylab = "", xlab = "")
axis(side = 4, at = pretty(range(bias_MC^2)), cex.axis = 2)
mtext("Bias^2", side = 4, line = 3, cex = 2)
title("Bias-Variance Tradeoff", cex.main = 2)
dev.off()

# choose optimal h acc CV
h_cv = h[which(L_CV == min(L_CV))]
f_cv = fNW(x = x, X = x, Y = y, h = h_cv)

png("cv.png", width = 900, height = 900, bg = "transparent")
plot(x, f, type = "l",col = "blue3", lwd = 3, ylab = "", 
     xlab = "x", cex.lab = 2, cex.axis = 2, ylim = range(y))
title("Simulated Data Estimated with CV", cex.main = 2)
points(x, y, pch = 19, col = "red3", cex = 0.7)
lines(x, f_cv, lwd = 3)
dev.off()

# choose optimal h acc GCV
h_gcv = h[which(L_GCV == min(L_GCV))]
f_gcv = fNW(x = x, X = x, Y = y, h = h_gcv)

png("gcv.png", width = 900, height = 900, bg = "transparent")
plot(x, f, type = "l", col = "blue3", lwd = 3, ylab = "", 
     xlab = "x", cex.lab = 2, cex.axis = 2, ylim = range(y))
title("Simulated Data Estimated with GCV", cex.main = 2)
points(x, y, pch = 19, col = "red3", cex = 0.7)
lines(x, f_gcv, lwd = 3)
dev.off()

```

automatically created on 2021-05-18

### IPYNB Code
```ipynb

{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "\n",
    "# PAC learning : Bias-variance trade-off,  $L_{emp}$, $L_{CV}$ and $L_{GCV}$ \n",
    "\n",
    "## Load packages\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import scipy as sp\n",
    "import random as r\n",
    "from scipy.stats import norm\n",
    "from sklearn.neighbors import KernelDensity\n",
    "import matplotlib.pyplot as plt\n",
    "import os"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Data generation process"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "\n",
    "# Fix some general setup variables\n",
    "\n",
    "r.seed(201053)\n",
    "\n",
    "# sample size n \n",
    "\n",
    "n = 100\n",
    "\n",
    "# Monte Carlo iterations J\n",
    "\n",
    "J = 10\n",
    "\n",
    " # number  of repetitions for the CV exercise J_cv \n",
    "\n",
    "J_cv = 1      \n",
    "\n",
    "# noise variance\n",
    "\n",
    "sig2 = 0.3\n",
    "\n",
    "# Define colors\n",
    "\n",
    "col_pink = 'hotpink'\n",
    "col_blue = 'b'\n",
    "col_red  = 'indianred'\n",
    "col_navy = 'navy'\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 1. Generate  X$\\sim U [0,1]$"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [],
   "source": [
    "x = sp.random.uniform(0,1,n)\n",
    "x.sort()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 2. Generate $y_i = \\sin^3(2\\pi x_i^3) + \\varepsilon_i,  i\\in [1,n]$, where $\\varepsilon_i \\sim N(0,\\sigma_\\varepsilon^2)$ and $\\sigma^2_\\varepsilon = 0.3$\n",
    "\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "2a) Generate $\\varepsilon_i \\sim N(0,\\sigma_\\varepsilon^2)$ with $\\sigma^2_\\varepsilon $ to be determined"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Only for illustration - regenerated in next block!\n",
    "\n",
    "sigma = np.sqrt(sig2)\n",
    "eps = np.random.normal(0, sigma,  n)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "2b) Generate $y_i = \\sin^3(2\\pi x_i^3) + \\varepsilon_i$ for $i\\in 1:n$"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "def f(x):\n",
    "    return (np.sin(2*np.pi*(x**3)))**3\n",
    "\n",
    "def gendata(f, sigma, n):\n",
    "    return ( f + np.random.normal(0, sigma, n) ) \n",
    "            \n",
    "f_x = f(x)\n",
    "y = gendata(f_x, sigma, n)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Data plot"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[<matplotlib.lines.Line2D at 0x7f84940c9bd0>]"
      ]
     },
     "execution_count": 7,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "fig, ax = plt.subplots(1, 1)\n",
    "ax.plot(x, y, 'o', color=col_red)\n",
    "ax.plot(x,f_x,'-g',color=col_blue, lw=3)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Kernel regression smoothing with NW kernel estimator"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Smoothing with Gaussian Kernel\n",
    "\n",
    "In this task we choose $K_h$ to be a Gaussian Kernel. The Naradaya Watson (NW)\n",
    "estimator is given by\n",
    "\\begin{align}\n",
    "\\hat{f}_h(t):=\\sum_{i=1}^n\\frac{K_h(t-x_i)}{\\sum_{i=1}^nK_h(t-x_i)}y_i\n",
    "\\end{align}\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Naradaya-Watson at points t, a mesh of grid points\n",
    "def nw_inner(t, x, y, h):\n",
    "    K_h = norm.pdf( (t-x)/h )/h\n",
    "    den = sum(K_h)\n",
    "    w   = K_h/den\n",
    "    nw  = np.dot(w, y) \n",
    "    return nw\n",
    "\n",
    "# f_h estimation \n",
    "def fh(ts, x, y, h):\n",
    "    return [nw_inner(t, x ,y ,h) for t in ts]\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [],
   "source": [
    "# trial bandwidth for Gaussian kernel, note that in terms of Uni (box car) kernel \n",
    "# this has to be multiplied with 1.35/0.77, i.e. we for h=0.01 are averaging over the window (t - 0.0175, t + 0.0175)\n",
    "# LM change the NW smoother !!!! \n",
    "\n",
    "h = 0.02\n",
    "\n",
    "# NW estimation\n",
    "\n",
    "f_h = fh(x, x, y, h) "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[<matplotlib.lines.Line2D at 0x7f848fb87050>]"
      ]
     },
     "execution_count": 10,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "fig, ax = plt.subplots(1, 1)\n",
    "ax.plot(x, y, 'o', color= col_red)\n",
    "ax.plot(x,f_h,'-g',color= col_blue, lw=3)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Bias-variance decomposition\n",
    "\n",
    "\\begin{equation}\n",
    "\\begin{aligned}\n",
    "\\mathrm{E}\\left[(y-\\hat{f})^{2}\\right] &=\\mathrm{E}\\left[y^{2}+\\hat{f}^{2}-2 y \\hat{f}\\right]=\\mathrm{E}\\left[y^{2}\\right]+\\mathrm{E}\\left[\\hat{f}^{2}\\right]-\\mathrm{E}[2 y \\hat{f}] \\\\\n",
    "&=\\operatorname{Var}[y]+\\mathrm{E}[y]^{2}+\\operatorname{Var}[\\hat{f}]+\\mathrm{E}[\\hat{f}]^{2}-2 f \\mathrm{E}[\\hat{f}]=\\operatorname{Var}[y]+\\operatorname{Var}[\\hat{f}]+(f-\\mathrm{E}[\\hat{f}])^{2} \\\\\n",
    "&=\\operatorname{Var}[y]+\\operatorname{Var}[\\hat{f}]+\\mathrm{E}[f-\\hat{f}]^{2}=\\sigma^{2}+\\operatorname{Var}[\\hat{f}]+\\operatorname{Bias}[\\hat{f}]^{2}\n",
    "\\end{aligned}\n",
    "\\end{equation}\n",
    "\n",
    "Source: Wikipedia\n",
    "\n",
    "Estimation for some bandwidths e.g. $h \\in [0.01, 0.02, \\dots, 0.1]$"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [],
   "source": [
    "h = list( np.linspace(0.01, 0.3, 25) )\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Find $\\hat{f}$ \n",
    "Kernel regression for all bandwidths\n",
    "\n",
    "Bias: $\\mathrm{E}[\\hat{f} - f] =  \\mathrm{E}[\\hat{f}]  - f $\n",
    "\n",
    "Variance: $E\\left[\\hat{f}-\\mathrm{E} \\hat{f}\\right]^{2}$"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Repeat data generation J times\n",
    "\n",
    "f_h = np.zeros( (n, len(h)) )\n",
    "E_fh =  np.zeros( (n, len(h)) )\n",
    "bias = np.zeros( (n, len(h)) )\n",
    "variance = np.zeros( (n, len(h)) )\n",
    "CV = np.zeros( len(h) )\n",
    "CV_all = np.zeros( (J_cv, len(h)) )   # collects the CVs over simulations\n",
    "L_emp = np.zeros( (J_cv, len(h)) )    # collects the emp error over simulations\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "metadata": {},
   "outputs": [],
   "source": [
    "#  Computation of bias\n",
    "\n",
    "r.seed(201053)      # fix seed to have bias and variance calculated on the same random set\n",
    "\n",
    "for j in range(1, J):\n",
    "    y = gendata(f_x, sigma, n)\n",
    "    for i in range(len(h)):\n",
    "        f_h[:,i] = fh(x, x, y, h[i])\n",
    "        E_fh[:,i] = E_fh[:,i] + f_h[:,i] \n",
    "    # end h loop\n",
    "        \n",
    "# end J loop\n",
    "\n",
    "E_fh = E_fh/J\n",
    "\n",
    "# now calc bias\n",
    "for i in range(len(h)):\n",
    "    bias[:,i] = E_fh[:,i] - f_x\n",
    "    # end h loop\n",
    "\n",
    "bias2 = bias**2\n",
    "Int_bias_h = np.sum(bias2, axis = 0)/n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[<matplotlib.lines.Line2D at 0x7f848fbfab90>]"
      ]
     },
     "execution_count": 14,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAXoAAAD4CAYAAADiry33AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADh0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uMy4xLjEsIGh0dHA6Ly9tYXRwbG90bGliLm9yZy8QZhcZAAAgAElEQVR4nO3de3xU1bn/8c9DEERRixCrAgoKWuHUWhzxjhdE0VbAiop3kYoepdZ66q21vyq1x3utVaviXVtEpGKpokjBu2gJKCJysAGtRCqNQLGCAoHn98eaOJNhIDvJZHZm5vt+vfJy9lp7T57l6DMra6+9lrk7IiJSvFrFHYCIiDQvJXoRkSKnRC8iUuSU6EVEipwSvYhIkWsddwCZOnXq5N26dYs7DBGRgjJr1qzP3L08W12LS/TdunWjoqIi7jBERAqKmf1jU3WRhm7MbKCZLTCzSjO7Mkt9PzObbWY1ZjY0o24XM3vBzOab2ftm1q2hDRARkcarN9GbWRlwF3As0As41cx6ZZz2MXAOMDbLWzwK3OzuewF9gX81JWAREWmYKEM3fYFKd18EYGbjgMHA+7UnuPtHyboN6RcmvxBau/vU5Hlf5CZsERGJKsrQTWdgcdpxVbIsij2Af5vZU2b2tpndnPwLoQ4zG2lmFWZWUV1dHfGtRUQkiiiJ3rKURV0gpzVwKPBTYD9gN8IQT903cx/j7gl3T5SXZ71pLCIijRQl0VcBXdOOuwBLIr5/FfC2uy9y9xrgaaBPw0IUEZGmiJLoZwI9zay7mbUBhgGTIr7/TKCDmdV2048kbWxfRESC99+Hn/0MVq/O/XvXm+iTPfFRwBRgPjDe3eeZ2WgzGwRgZvuZWRVwEnCvmc1LXrueMGwzzczmEoaB7st9M0RECs+KFXDPPbD//tC7N1x/PUycmPvfYy1tPfpEIuF6YEpEitX69fDXv8LDD4ekvmZN3fqjjoKpUxv+vmY2y90T2epa3JOxIiLFaMECeOQRePRR+OSTjevbtIFBg2D48Nz/biV6EZFmsnIljB8PDz0EM2ZkP2fffeGcc+DUU6Fjx+aJQ4leRCSH3OHll+H+++Gpp+DLLzc+p7wczjgjJPi9927+mJToRURy4Isv4LHH4K67YN68jetbt4bvfz8k9+OOgy22yF9sSvQiIk2wYAH8/vfh5urnn29cv/feYdz9tNNghx3yHh6gRC8i0mDr18Ozz4be+wsvbFy/9dZw5pkwciTssw9YtvUF8kiJXkQkomXL4IEHQg/+H1lWf99jD7joIjj7bNhuu/zHtylK9CIi9Zg9G+68Ex5/HL76qm6dGRx/PIwaBf37Q6sWuEGrEr2ISBbuMGUK3HBDmEWTafvt4bzz4IILoKXvfqpELyKSpqYGJkwICX7OnI3r+/SBH/0ITjkF2rXLf3yNoUQvIkIYknn4Ybj5Zli0qG5d69Zw8skhwe+/f/w3VxtKiV5EStrKlXD33fDb38LSpXXrttoqDM9ceinssks88eWCEr2IlKRPPw3J/e67N57/3qEDXHxxuMHaqVM88eWSEr2IlJSFC+GWW8L6M5krR3buDP/zP6EX3759PPE1ByV6ESkJ//gHjB4dVpBcv75u3Z57whVXwOmnh1Uki40SvYgUtU8/hV//GsaMgbVr69b17QtXXgmDB7fM+e+5okQvIkVp2TK46Sa4446NV5Ds3x9+/nM4/PDCm0HTGJG+w8xsoJktMLNKM7syS30/M5ttZjVmNjRL/bZm9omZ3ZmLoEVENuXzz8MQzW67hUSfnuQPOghefDHs8HTEEaWR5CFCj97MyoC7gAFAFTDTzCa5e/om3x8D5xD2h83mV0CWZ8tERHLjyy/DImM33BB68+n22ScM3xx7bOkk93RRhm76ApXuvgjAzMYBg4GvE727f5Ss25B5sZntC3wTeB7Iup+hiEhjrV0bNvm47jr45z/r1n3rW6F3f+KJxT0GX58oTe8MLE47rkqW1cvMWgG3ApfVc95IM6sws4rq6uooby0iJW7DhrDRx557hhUj05N8t27hKdf33oOTTirtJA/REn22P3Q84vtfCEx298WbO8ndx7h7wt0T5eXlEd9aRErVG2/AAQfAWWfBRx+lynfaKSwhvGBBWCq4rCy2EFuUKEM3VUDXtOMuwJKI738gcKiZXQi0B9qY2RfuvtENXRGR+ixeHOa7P/543fKOHeGqq+DCCwtnobF8ipLoZwI9zaw78AkwDDgtypu7++m1r83sHCChJC8iDbV6dZhBkzmLpm3b8CTrFVfAttvGF19LV+/QjbvXAKOAKcB8YLy7zzOz0WY2CMDM9jOzKuAk4F4zy7I1rohIw7jD2LFhHP7aa+sm+ZNOgvnzw2waJfnNM/eow+35kUgkvKKiIu4wRCRmM2fCj38MM2bULf/ud8NiZP36xRNXS2Vms9w968zGEr8XLSItzZIl4UZq3751k/wOO4RplDNnKsk3lJZAEJEW4auv4NZb4frrYdWqVHmbNnDJJWHJAg3RNI4SvYjEbtq0sPdqZWXd8hNOCDs+7b57PHEVCw3diEhsqqvDXPijjqqb5PfeOyT/p55Sks8FJXoRyTv38OTqXnuFp1trbbddeOBp9mw48sjYwis6GroRkbz64IMwTPPii3XLTzklzKbZccd44ipm6tGLSF6sXQu/+lUYlklP8rvuCpMnw7hxSvLNRT16EWl2r70GI0eGB5xqlZXBT34C11wDW28dW2glQYleRJrNihVheYL77qtbnkiEsn32iSeuUqOhGxHJOXd44olwszU9ybdvD7ffDm++qSSfT+rRi0hOVVfD+efDxIl1ywcPDvu3du2a/TppPkr0IpIzTz8dxuLT9w/aeWe4887w8JPEQ0M3ItJkK1fCOeeEZJ6e5M8/P9yAVZKPl3r0ItIk06bB8OFhU5BaO+8MDzwAAwfGF5ekqEcvIo2yejVcfHFYviA9yZ92WtirVUm+5VCPXkQa7K23who1H3yQKuvYEe65B4YOjS8uyS5Sj97MBprZAjOrNLONtgI0s35mNtvMasxsaFr5PmY2w8zmmdm7ZnZKLoMXkfxauxauvhoOOqhukv/+90MvXkm+Zaq3R29mZcBdwADCRuEzzWySu7+fdtrHwDnATzMuXw2c5e5/N7OdgVlmNsXd/52T6EUkb+bODb34d95JlW2zTVifZvhwMIsvNtm8KEM3fYFKd18EYGbjgMHA14ne3T9K1m1Iv9DdP0h7vcTM/gWUA0r0IgViw4awIcjVV4cefa3DDgsrUHbrFldkElWUoZvOQNqtFqqSZQ1iZn2BNsDCLHUjzazCzCqq0+dmiUisli6FY4+Fyy9PJfktt4TbboPp05XkC0WURJ/tD7IG7ShuZjsBjwHD3X1DZr27j3H3hLsnysvLG/LWItJMpk0LyxS88EKqLJEIa8Vfcgm00py9ghHlo6oC0h9a7gIsifoLzGxb4Fngand/s2HhiUi+1dSEYZoBA+DTT0OZGVx1FbzxRli/RgpLlDH6mUBPM+sOfAIMA06L8uZm1gaYCDzq7k82OkoRyYvFi+HUU+H111NlO+wQdoE6+uj44pKmqbdH7+41wChgCjAfGO/u88xstJkNAjCz/cysCjgJuNfM5iUvPxnoB5xjZu8kf7RmnUgLNGlSGKpJT/L9+8OcOUryhc7cGzTc3uwSiYRXVFTEHYZIyVizJqwZf/vtqbKyMhg9OpSXlcUXm0RnZrPcPZGtTk/GipSwysqwV+vs2amyrl1h7Fg45JD44pLc0n1zkRL1+OPQp0/dJD9oUHggSkm+uCjRi5SYVatgxIiw+Nh//hPKttgiPOH69NOw/fbxxie5p6EbkRJSWQk/+EFYzqBWjx4wbhzsu298cUnzUo9epEQ880x44Ck9yZ92GsyapSRf7JToRYrchg1wzTVw/PFhJyiAtm3Dpt1/+ANsu22s4UkeaOhGpIitWAFnnAGTJ6fKunaFp54KvXspDUr0IkXq3XfDXq2LFqXK+vcPs220pFRp0dCNSBEaOxYOOKBukr/8cnj+eSX5UqQevUgRWbcOLrus7lOu7duHdeNPPDG2sCRmSvQiReLTT+Hkk+HVV1Nle+4JEydqxclSp6EbkSIwY0aYIpme5E84Af72NyV5UaIXKWjucPfdYVu/JcldIlq1guuvhz/9SVMnJdDQjUiBWrsWLroI7r8/Vbb99uEp1wED4otLWh4lepEC9Nln4ebqK6+kyvr0Cb147eMqmTR0I1Jg3nsP9tuvbpI/4wx47TUleckuUqI3s4FmtsDMKs3syiz1/cxstpnVmNnQjLqzzezvyZ+zcxW4SCl65hk48ED46KNwbAY33ACPPgrt2sUamrRg9Q7dmFkZcBcwgLBR+Ewzm+Tu76ed9jFwDvDTjGu3B34JJAAHZiWvXZGb8EVKgzvcckvY8al2U7j27eGPfwxryItsTpQx+r5ApbsvAjCzccBg4OtE7+4fJes2ZFx7DDDV3Zcn66cCA4HHmxy5SIlYswZGjgy99lrduoU9Xr/97djCkgISJdF3BhanHVcB+0d8/2zXdo54rUjJW7o0zIefMSNVdsghYVEyLWUgUUUZo7csZVF3FI90rZmNNLMKM6uorq6O+NYixe2dd8JN1/Qkf+65MG2akrw0TJREXwV0TTvuAiyJ+P6RrnX3Me6ecPdEuf4LFuGpp+Dgg2Fx8u/hVq3gN78Jc+bbtIk3Nik8URL9TKCnmXU3szbAMGBSxPefAhxtZh3MrANwdLJMRLJwh1//OsyRX706lG27bZht85OfhFk2Ig1V7xi9u9eY2ShCgi4DHnT3eWY2Gqhw90lmth8wEegAHG9m17p7b3dfbma/InxZAIyuvTErInV99VUYmnk8barC7ruHm669esUXlxQ+c4863J4fiUTCKyoq4g5DJK+qq2HIEHjjjVTZEUfAk09Cx47xxSWFw8xmuXvWfcP0ZKxIzBYsCA9BpSf588+HKVOU5CU3lOhFYvTKKyHJL1wYjs3CTde774Yttog3NikeWtRMJCZ/+EMYk1+3Lhy3axe2ABwyJN64pPioRy+SZ+5w7bVw5pmpJP/Nb8LLLyvJS/NQj14kj9asgfPOg8ceS5X17g3PPgu77hpfXFLclOhF8mT5cvjBD0LPvdZRR8GECbDddvHFJcVPQzciebBwIRx0UN0kP2IETJ6sJC/NT4lepJm98QYccECYRlnr+uvhvvs0s0byQ0M3Is1o/Hg466wwNg/Qtm1Ybvjkk+ONS0qLevQizcAdbrwRTjklleQ7dYLp05XkJf/UoxfJsZoauPji8NBTrT33DDNrdt89vrikdCnRi+TQqlVw6qnwl7+kyvr1g4kTYfvt44tLSpuGbkRyZOnSsBBZepI/9VR44QUleYmXEr1IDnzwQVizZubMVNkVV4RlDtq2jS8uEdDQjUiTvfEGDBoEy5aF41at4I474MIL441LpJYSvUgTPPUUnH562DQEwsJk48aFxC/SUmjoRqSRbr8dhg5NJfnycnjpJSV5aXkiJXozG2hmC8ys0syuzFLf1syeSNa/ZWbdkuVbmNkjZjbXzOab2VW5DV8k/zZsgEsvhUsuCfPlAXr2hBkzoG/feGMTyabeRG9mZcBdwLFAL+BUM8vcwXIEsMLdewC3ATcmy08C2rr7t4F9gfNrvwRECtFXX8GwYXDbbamy2t2hNEdeWqooPfq+QKW7L3L3tcA4YHDGOYOBR5KvJwD9zcwAB7Y2s9ZAO2At8HlOIhfJs2XLwmqTTz6ZKjvhBJg2LTz1KtJSRUn0nYHFacdVybKs57h7DbAS6EhI+quAfwIfA7e4+/LMX2BmI82swswqqqurG9wIkeb24Ydw8MHw+uupsosvDkm/Xbv44hKJIkqityxlHvGcvsB6YGegO/A/ZrbbRie6j3H3hLsnysvLI4Qkkj8VFWF4Jn31yVtvhd/+FsrK4otLJKooib4K6Jp23AVYsqlzksM02wHLgdOA5919nbv/C3gdSDQ1aJF8mTwZDj88PPUK4eGn8ePDzVjL1r0RaYGiJPqZQE8z625mbYBhwKSMcyYBZydfDwWmu7sThmuOtGBr4ADg/3ITukjzuv/+MFVy1apw3KEDTJ0KJ50Ub1wiDVVvok+OuY8CpgDzgfHuPs/MRptZ7YzhB4COZlYJXArUTsG8C2gPvEf4wnjI3d/NcRtEcsod/t//C3u7rl8fynbdNYzPH3povLGJNIa5Zw63xyuRSHhFRUXcYUiJWrcuJPhHHkmV9ekDzzwDO+0UX1wi9TGzWe6edWhcSyCIJH3+eXjSderUVNnAgWFMfptt4otLpKm0BIIIsGRJWDc+Pcmfey5MmqQkL4VPiV5K3rx5YfPuOXNSZddcE27GavNuKQYaupGS9tJLMGQIrFwZjsvK4L77YPjwWMMSySkleilZ48bB2WfD2rXhuH17mDABjjkm3rhEck1DN1Jy3OGWW8I2f7VJfscd4ZVXlOSlOKlHLyWlpgZ+/GP4/e9TZXvtBc89F+bKixQjJXopGV98EZYYfvbZVNmhh8LTT2vzbiluGrqRklA7fTI9yQ8bBi+8oCQvxU+JXore3Llh+uTbb6fKrroK/vhH2HLL+OISyRcN3UhRmzoVTjwR/vOfcFxWBvfcAz/8YbxxieSTevRStB58EI47LpXkt9kmDN0oyUupUaKXouMOV18NI0aEWTYAXbrAa69p+qSUJg3dSFFZsyasUTN2bKpsn33C6pOdMzfAFCkRSvRSNJYvD8sZvPpqquy448ITsFqYTEqZhm6kKCxcCAcdVDfJX3AB/PnPSvIikRK9mQ00swVmVmlmV2apb2tmTyTr3zKzbml1e5vZDDObZ2ZzzUwT2iSn3nxz4827b7opPP3aWn+zitSf6M2sjLAl4LFAL+BUM+uVcdoIYIW79wBuA25MXtsa+ANwgbv3Bg4H1uUseil5Y8eGzburq8Nx7ebdl12mzbtFakXp0fcFKt19kbuvBcYBgzPOGQzUbr42AehvZgYcDbzr7nMA3H2Zu6/PTehSyjZsgF/8Ak4/PdyABejUCaZP1+bdIpmiJPrOwOK046pkWdZzkpuJrwQ6AnsAbmZTzGy2mV2e7ReY2UgzqzCziurarpnIJqxaBSefDNddlyrba68whHPQQfHFJdJSRUn02f4AztxRfFPntAYOAU5P/vMEM+u/0YnuY9w94e6J8vLyCCFJqaqqCmvW/OlPqbKBA2HGDNh99/jiEmnJoiT6KqBr2nEXYMmmzkmOy28HLE+Wv+zun7n7amAy0KepQUtp+tvfoG9fmD07VXbJJfCXv8B228UXl0hLFyXRzwR6mll3M2sDDAMmZZwzCTg7+XooMN3dHZgC7G1mWyW/AA4D3s9N6FJKnngCDjsM/vnPcNy6Ndx7L9x2m2bWiNSn3v9F3L3GzEYRknYZ8KC7zzOz0UCFu08CHgAeM7NKQk9+WPLaFWb2G8KXhQOT3f3ZrL9IJIsNG2D0aLj22lRZhw5h6OaII+KLS6SQWOh4txyJRMIrKiriDkNagNWrwybd48enyvbcMwzV9OwZX1wiLZGZzXL3RLY6/dErLdInn4TlDNK/848+OgzhfOMb8cUlUoi0BIK0OLNmhZuu6Ul+1KiwxLCSvEjDKdFLizJuXNjHdUlyXldZWVjK4I47dNNVpLH0v460COvWwRVXhFk0tb7xDXjySTjqqPjiEikGSvQSu6VL4ZRT4OWXU2V77BFuuu6xR3xxiRQLDd1IrN58E/bdt26SHzw4PBylJC+SG0r0Egv38MBTv35hhg2E1Savuw6eekpPuorkkoZuJO+++gouuihs3l2rQwd4/HHt6SrSHJToJa8+/hhOPLHu1Ml99gm9+O7d44tLpJhp6EbyZtq0MB6fnuTPPBNef11JXqQ5KdFLs3MPW/sdfTR89lkoa90a7rwTHnkEttoq3vhEip2GbqRZ/ec/Yb2a9PXjd9wRJkyAgw+OLy6RUqJEL83mvffCTlDz56fKDj44PAS1007xxSVSajR0IznnDnffDfvtVzfJ/+hHYU9XJXmR/FKPXnJq+XL44Q9h4sRUWbt2Yc78mWfGF5dIKVOil5x59VU4/XRYnLaV/N57h4XK9torvrhESp2GbqTJ1q8Pu0AdfnjdJH/RRfDWW0ryInGLlOjNbKCZLTCzSjO7Mkt9WzN7Iln/lpl1y6jfxcy+MLOf5iZsaSmqquDII+GXvwzb/kF4ynXixDB9csst441PRCIkejMrA+4CjgV6AaeaWa+M00YAK9y9B3AbcGNG/W3Ac00PV1qSP/8ZvvMdeOWVVFm/fjBnTtgdSkRahig9+r5Apbsvcve1wDhgcMY5g4FHkq8nAP3NzADMbAiwCJiXm5Albl99FWbQDBkSbr4CtGoF11wTZtV07RpreCKSIcrN2M5A2sgrVcD+mzrH3WvMbCXQ0cy+BK4ABgCbHLYxs5HASIBddtklcvCSf/Pnw7Bh8O67qbIuXWDs2LAzlIi0PFF69JalzCOecy1wm7t/sblf4O5j3D3h7ony8vIIIUm+ucN990EiUTfJDxkShmqU5EVarig9+iog/Y/xLsCSTZxTZWatge2A5YSe/1Azuwn4BrDBzL5y9zubHLnkzeLFcN55MGVKqqxtW/jNb+C//zusIy8iLVeURD8T6Glm3YFPgGHAaRnnTALOBmYAQ4Hp7u7A1/08M7sG+EJJvnC4hzXjL70UPv88Vb7XXvDEE/Dtb8cXm4hEV+/QjbvXAKOAKcB8YLy7zzOz0WY2KHnaA4Qx+UrgUmCjKZhSWBYvhmOPDU+51iZ5M/jJT8Iyw0ryIoXDQse75UgkEl6RvmC55NWmevE9e4byQw6JLzYR2TQzm+XuiWx1ejJWvlZVBccdl70X/847SvIihUpr3Qju8NBDIaGn9+J79AjlSvAihU09+hJXVQXf+x6MGFG3F3/JJWHapJK8SOFTj75EucPDD4de/MqVqfIePcJYvObFixQP9ehL0N//HvZvPffcVJI3gx//WA8/iRQj9ehLyJo1cOON8L//G17X2n33MBavBC9SnJToS8SLL4anWBcsSJW1agUXXwzXXQdbbx1fbCLSvJToi1x1Nfz0p/Doo3XLE4mwvV+fPvHEJSL5ozH6IrVhAzzwAHzrW3WT/DbbwB13wJtvKsmLlAr16IvQvHlwwQXw2mt1y086CX77W9h553jiEpF4qEdfRFavhquugn32qZvku3WDZ5+F8eOV5EVKkXr0ReK558Jm3B9+mCpr3TqMz//iF7DVVvHFJiLxUqIvcIsXh4ee/vSnuuUHHwz33AP/9V/xxCUiLYeGbgrUunVw881hbfj0JN+hQ9gJ6pVXlORFJFCPvgC98gpceGG46ZrurLNC8t9hh3jiEpGWSYm+gCxdCpdfvvGc+N694fe/h3794olLRFq2SEM3ZjbQzBaYWaWZbbR7lJm1NbMnkvVvmVm3ZPkAM5tlZnOT/zwyt+GXhvXr4e67N54Tv/XWcMst8PbbSvIismn19ujNrAy4CxhA2AR8pplNcvf3004bAaxw9x5mNgy4ETgF+Aw43t2XmNl/EbYj7JzrRhSzmTPDME3mpltDh8Jtt0GXLvHEJSKFI0qPvi9Q6e6L3H0tMA4YnHHOYOCR5OsJQH8zM3d/292XJMvnAVuaWdtcBF7sVqwICX7//esm+R494Pnn4cknleRFJJooib4zsDjtuIqNe+Vfn5PcTHwl0DHjnBOBt919TUY5ZjbSzCrMrKK6ujpq7EXJPQzP7LlnGK6p3dK3bVu49lqYOxeOOSbeGEWksES5GWtZyjJ3FN/sOWbWmzCcc3S2X+DuY4AxEDYHjxBTUZo/P6ww+fLLdcsHDoQ77wzLCYuINFSUHn0V0DXtuAuwZFPnmFlrYDtgefK4CzAROMvdFzY14GK0ejX87Gfwne/UTfJduoQ58pMnK8mLSONFSfQzgZ5m1t3M2gDDgEkZ50wCzk6+HgpMd3c3s28AzwJXufvruQq6mDz7bJgeef314SEogLKysHTB/Pnwgx+E3Z9ERBqr3kSfHHMfRZgxMx8Y7+7zzGy0mQ1KnvYA0NHMKoFLgdopmKOAHsAvzOyd5I8e5yEsXXDiifD978NHH6XKDz44TJe8+WZo3z628ESkiJh7yxoSTyQSXpE5l7CIrFsHv/sd/PKXsGpVqnz77eGmm2D48LDzk4hIQ5jZLHdPZKvTk7F5NGNGWCf+3Xfrlg8fHpJ8p07xxCUixU19xzxYvhxGjoSDDqqb5Hv3DuvWPPigkryINB8l+mZUUwNjxoQ58ffdlyrfaiu48cYwFn/oofHFJyKlQUM3zcA9PL162WUbrzA5aFAYo99113hiE5HSo0SfY3PmhKmRf/1r3fKuXcOm3IMzF48QEWlmGrrJkU8+gXPPhe9+t26Sb98errsO/u//lORFJB7q0TfRF1+EGTO33AJffpkqb9UKzjsvrE/zzW/GF5+IiBJ9I9XUwEMPhY23ly6tW/e974Xk36tXPLGJiKRTom+gzd1o/c534NZboX//eGITEclGiT6impqwLs3vfgfTp9et69wZfv1rOOOMsE6NiEhLokRfj48/hgceCD+ffFK3buut4cor4dJLw9x4EZGWSIk+i/Xrw9LA994Lzz0HGzbUra+90XrNNbDjjrGEKCISmRJ9mqqq0HO///7wOtMOO4QplOedB7vtlv/4REQao+QT/fr1MGVK6L0/88zGvXeAo44Ka9UMHgxt2uQ/RhGRpijJRP/vf8Orr4bdnJ58MozDZyovD6tKnnde2JBbRKRQlUSiX7YsrBL58svhZ86c1KbbmY44As4/H4YMCRtyi4gUukiJ3swGArcDZcD97n5DRn1b4FFgX2AZcIq7f5SsuwoYAawHLnb3KTmLfhP+9a+Q2F96KST2997b/PmdOsE554Te+x57NHd0IiL5VW+iN7My4C5gAGET8JlmNsnd3087bQSwwt17mNkw4EbgFDPrRdhjtjewM/BXM9vD3dfnshFr1sDEiake+/z5mz+/VSvo0wcOOwwOPxwGDFDvXUSKV5QefV+g0t0XAZjZOGAwkJ7oBwPXJF9PAO40M0uWj3P3NcCHyT1l+2WmFk4AAAVWSURBVAIzchN+YBZmw6SvNZOudWtIJEJiP+ywsC/rttvmMgIRkZYrSqLvDCxOO64C9t/UOe5eY2YrgY7J8jczru3c6Gg3oU0bOPDA1BOrW2wB+++fSuwHHqiNtkWkdEVJ9JalLPNW5qbOiXItZjYSGAmwyy67RAhpY8OHh92aDjsMDjgA2rVr1NuIiBSdKIm+CuiadtwFWLKJc6rMrDWwHbA84rW4+xhgDEAikdjEfJjNO+OMxlwlIlL8omw8MhPoaWbdzawN4ebqpIxzJgFnJ18PBaa7uyfLh5lZWzPrDvQE/pab0EVEJIp6e/TJMfdRwBTC9MoH3X2emY0GKtx9EvAA8FjyZutywpcByfPGE27c1gAX5XrGjYiIbJ75pp4cikkikfCKioq4wxARKShmNsvdE9nqtGesiEiRU6IXESlySvQiIkVOiV5EpMi1uJuxZlYN/COjuBPwWQzhNKdia1OxtQeKr03F1h4ovjY1pT27unt5tooWl+izMbOKTd1NLlTF1qZiaw8UX5uKrT1QfG1qrvZo6EZEpMgp0YuIFLlCSfRj4g6gGRRbm4qtPVB8bSq29kDxtalZ2lMQY/QiItJ4hdKjFxGRRlKiFxEpcrEnejMbaGYLzKzSzK7MUt/WzJ5I1r9lZt3S6q5Kli8ws2PyGfemNLY9ZtbNzL40s3eSP/fkO/ZNidCmfmY228xqzGxoRt3ZZvb35M/ZmdfGoYntWZ/2GWUu1x2bCG261MzeN7N3zWyame2aVleIn9Hm2lOon9EFZjY3GfdryT23a+ualuvcPbYfwrLHC4HdgDbAHKBXxjkXAvckXw8Dnki+7pU8vy3QPfk+ZQXcnm7Ae3HG34Q2dQP2Bh4FhqaVbw8sSv6zQ/J1h0JtT7Lui7g/k0a26Qhgq+Tr/077765QP6Os7Snwz2jbtNeDgOeTr5uc6+Lu0X+98bi7rwVqNx5PNxh4JPl6AtA/c+Nxd/8QqN14PE5NaU9LVW+b3P0jd38X2JBx7THAVHdf7u4rgKnAwHwEvRlNaU9LFaVNL7r76uThm4Td3qBwP6NNtaelitKmz9MOtya17WqTc13ciT7bxuOZm4fX2XgcSN94vL5r860p7QHobmZvm9nLZnZocwcbUVP+PRfqZ7Q5W5pZhZm9aWZDchtaozW0TSOA5xp5bT40pT1QwJ+RmV1kZguBm4CLG3Lt5kTZM7Y5NfvG43nWlPb8E9jF3ZeZ2b7A02bWO+NbPg5N+fdcqJ/R5uzi7kvMbDdgupnNdfeFOYqtsSK3yczOABLAYQ29No+a0h4o4M/I3e8C7jKz04CrCVu0NvkzirtH35CNx7FGbDyeZ41uT/LPsmUA7j6LMA63R7NHXL+m/Hsu1M9ok9x9SfKfi4CXgO/mMrhGitQmMzsK+DkwyN3XNOTaPGtKewr6M0ozDqj9a6Tpn1HMNyhaE27+dCd1g6J3xjkXUffm5fjk697UvUGxiPhvxjalPeW18RNu2HwCbB9ne6K2Ke3ch9n4ZuyHhJt8HZKvY21TE9vTAWibfN0J+DsZN9RaapsIyW4h0DOjvCA/o820p5A/o55pr48n7Mmdk1wXa+OTjTgO+CD5of08WTaa8C0NsCXwJOEGxN+A3dKu/XnyugXAsXG3pSntAU4E5iU/0NnA8XG3pQFt2o/Q61gFLAPmpV17brKtlcDwuNvSlPYABwFzk5/RXGBE3G1pQJv+CiwF3kn+TCrwzyhrewr8M7o9mQPeAV4k7YugqblOSyCIiBS5uMfoRUSkmSnRi4gUOSV6EZEip0QvIlLklOhFRIqcEr2ISJFTohcRKXL/H02cjtrsHjxoAAAAAElFTkSuQmCC\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.plot(h, Int_bias_h, color=col_blue, lw = 3)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "metadata": {
    "scrolled": true
   },
   "outputs": [],
   "source": [
    "r.seed(201053)\n",
    "for j in range(1,J):\n",
    "    y = gendata(f_x, sigma, n)\n",
    "    for i in range(len(h)):\n",
    "        f_h[:,i] = fh(x, x, y, h[i])\n",
    "        variance[:,i] = variance[:,i] +  (f_h[:,i]  - E_fh[:,i])**2\n",
    "    # end h loop\n",
    "        \n",
    "# end J loop  \n",
    "\n",
    "variance = variance/J\n",
    "\n",
    "Int_var_h = np.sum(variance, axis = 0)/n "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[<matplotlib.lines.Line2D at 0x7f84944aa750>]"
      ]
     },
     "execution_count": 16,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.plot(h, Int_var_h, color=col_red, lw = 3)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Plot bias-variance-tradeoff"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 2 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "fig, ax1 = plt.subplots()\n",
    "ax1.set_xlabel('bandwidth h')\n",
    "ax1.set_ylabel('squared bias', color=col_blue)\n",
    "ax1.plot(h, Int_bias_h, color=col_blue, lw = 3)\n",
    "ax1.tick_params(axis='y', labelcolor=col_blue)\n",
    "\n",
    "ax2 = ax1.twinx()  \n",
    "\n",
    "col_red = 'tab:red'\n",
    "ax2.set_ylabel('variance', color=col_red)  \n",
    "ax2.plot(h, Int_var_h, color=col_red, lw = 3)\n",
    "ax2.tick_params(axis='y', labelcolor=col_red)\n",
    "\n",
    "fig.tight_layout()  \n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[<matplotlib.lines.Line2D at 0x7fccf3637e20>]"
      ]
     },
     "execution_count": 17,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAXoAAAD4CAYAAADiry33AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADh0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uMy4yLjIsIGh0dHA6Ly9tYXRwbG90bGliLm9yZy+WH4yJAAAgAElEQVR4nO3deZyNdf/H8ddnxthJ4ZYsUYlU1km2ZF+zlbKEkpLuumm7S6nfXXfa63ZrtZTclhBFKlKkIktmIpJkUjGSJfs62/f3xxnHzBjmMMs155z38/GYhznfazmf7xzerrmu7/W9zDmHiIiErgivCxARkdyloBcRCXEKehGREKegFxEJcQp6EZEQV8DrAjJTpkwZV6VKFa/LEBEJGrGxsbucc2UzW5Yvg75KlSrExMR4XYaISNAws99PtUynbkREQpyCXkQkxCnoRURCnIJeRCTEKehFREKcgl5EJMQp6EVE8oHff9/LM88s5vff9+b4vvPlOHoRkXCwd+9RZs78kUmT1vD1175h8Ckpjscea5aj76OgFxHJQwkJycybt5HJk9fy0UcbOHYsOd3yyZPXMHz4NZhZjr2ngl5EJJc551i+PJ7Jk9cwffo6/vrryEnrREQYbdteTL9+tUhJcURGKuhFRPK9uLjdTJmyhsmT1xIXtzvTderVK0+/frXo1esKzj+/eK7UoaAXEclBBw8mMG3aD4wfv4ply+IzXady5XO4+eYr6du3FjVrZjoPWY5S0IuI5IA1a7YzZkwMkyat4cCBhJOWlyxZiBtvrEm/frW45poLiYjIuVMzWQko6M2sPTAKiATecs49l2F5DeAdoB4w3Dn3UpplpYC3gCsAB9zmnFuWM+WLiHjnyJFEZs78kdGjY1m6dMtJywsUiKBDh0vo168W1113KUWKRHlQZQBBb2aRwOtAGyAeWGlmc5xzP6ZZbTcwBOiWyS5GAZ8653qYWUGgaPbLFhHxzoYNuxgzJpYJE1azZ8/Rk5ZXr16aO++sT9++tShbtpgHFaYXyBF9AyDOObcJwMymAV0Bf9A753YAO8ysU9oNzawk0Ay4NXW9BODk32lERPK5hIRkZs1az5gxsSxa9NtJy6OiIrj++ssYPDiaa6+9MEeHR2ZXIEFfAUj7O0k8cHWA+78I2Am8Y2a1gVhgqHPuUMYVzWwQMAigcuXKAe5eRCR3xcfv5/XXv2X8+NXs2HFSdFGlSinuvLM+AwbUoVy53Bk1k12BBH1m/y25M9h/PeAfzrkVZjYKGAY8ftIOnRsLjAWIjo4OdP8iIrni++//5OWXlzF16g8kJaWkWxYRYXTpUp3Bg+vTps3FeXph9WwEEvTxQKU0rysCfwS4/3gg3jm3IvX1THxBLyKS7zjn+PzzTbz00lI+/3zTScsrVCjBHXfUY+DAelSsWNKDCs9OIEG/EqhmZlWBrUAvoE8gO3fO/WlmW8ysunNuA9CKNOf2RUTyg4SEZKZN+4GXXlrK2rU7TlrerNmF3HdfQ6677lIKFAi+uSCzDHrnXJKZ3QPMxze8crxzbp2ZDU5dPtrMzgdigJJAipndC9R0zu0H/gFMSR1xswkYkEt9ERE5I/v2HWXs2FhGjVrB1q0H0i2LiDB69KjJAw80okGDCh5VmDMCGkfvnJsLzM3QNjrN93/iO6WT2bargehs1CgikqM2b97HqFHLGTfuu5NubipaNIqBA+ty770Nueiicz2qMGfpzlgRCRsbNuxixIjFTJ26luTk9GM+ypUrxpAhVzN4cDTnnVfEowpzh4JeREJeXNxunnrqayZPXkNKSvqAv+yyMjz4YGP69LmSwoVDMxJDs1ciIsCmTXsYMeJrJk78/qQj+ObNq/DPfzamfftL8v3wyOxS0ItIyPn99708/fRi3nln9Ulj4Nu1u5gnnmhOw4aZXlYMSQp6EQkZW7bs45lnFvP226tITEwf8K1aVeXJJ5vTpEn43XmvoBeRoPfHHwd49tnFjB37HQkJ6R/Nd+21F/Lvf7egWbMLParOewp6EQla27cf5NlnlzB6dMxJz15t2rQy//53c1q0qOpNcfmIgl5Egs6RI4mMHLmcZ59dwsGD6cfBN2pUkX//uwWtWlXNVzNIeklBLyJBwznH1Kk/MGzYArZs2Z9uWYMGFXjyyea0a3exAj4DBb2IBIVly7Zw333zWbFia7r2yy8vy3PPtaZTp2oK+FNQ0ItIvvbbb3sZNmwB06evS9detmxRnnqqBQMH1gvKicbykoJeRPKl/fuP8eyzixk5cnm6C60FC0Zy330NefTRayhZspCHFQYPBb2I5CvJySm8/fYqHn980UlPdLrppst57rlWVK0aGpON5RUFvYjkG59//gsPPPDZSXPCN2hQgZEj29G4caVTbCmno6AXEc9t3bqfIUM+5YMP1qdrr1ixJM8914reva8M+flocpOCXkQ8k5ycwptvxvDoowvTzQtfrFgUw4Y15f77G1G0aJSHFYYGBb2IeGL16j+5886P+fbb9MMlb721Ds8805Ly5Ut4VFnoUdCLSJ46dCiBJ574kpEjl6ebOrh69dKMGXMd115bxbviQpSCXkTyzNy5G/n73z/h99/3+dsKFoxk+PBrePjhJhQqpEjKDfqpikiu27btAEOHfsqMGT+ma2/evAqjR3eievUyHlUWHhT0IpJrUlIcY8bEMGzYQvbvP+ZvL126CC+/3Jb+/Wtr2oI8oKAXkVyxdu12Bg36mOXL49O133JLbV56qS1lyhT1qLLwE9AEEWbW3sw2mFmcmQ3LZHkNM1tmZsfM7MFMlkea2Soz+zgnihaR/CspKYWnn/6aevXGpgv5atXO44sv+jNhQjeFfB7L8ojezCKB14E2QDyw0szmOOfSnmzbDQwBup1iN0OB9UDJ7JUrIvlZXNxu+vWblS7go6IieOSRpjzyyDUULqyTCF4I5Ii+ARDnnNvknEsApgFd067gnNvhnFsJJGbc2MwqAp2At3KgXhHJh5zznYuvXXt0upBv2LAi338/mCefbKGQ91AgP/kKwJY0r+OBq8/gPf4LPATo7geRELRt2wFuv/0j5s7d6G8rUCCCJ59szkMPNdEUwvlAIEGf2SVxl0nbyRuaXQfscM7FmlnzLNYdBAwCqFw5/J7SLhKM3n//R+6882P++uuIv61mzbJMmtSdevXKe1iZpBXIf7XxQNop4yoCfwS4/yZAFzP7Dd8pn5ZmNjmzFZ1zY51z0c656LJlywa4exHxwr59R7nlltn06DEjXcjfd19DYmLuUMjnM4Ec0a8EqplZVWAr0AvoE8jOnXOPAI8ApB7RP+ic63t2pYpIfvDll79xyy2z2bz5xN2tlSqVZMKEbrRsWdXDyuRUsgx651ySmd0DzAcigfHOuXVmNjh1+WgzOx+IwTeqJsXM7gVqOuf2n3LHIhJUjh5NYvjwhYwcuRyX5uRt3761ePXVDpQqVdi74uS0zLmATrfnqejoaBcTE+N1GSKSavXqP+nb9wPWrdvpbzvvvCKMGXMdPXrU9LAyOc7MYp1z0Zkt03gnETkl5xxjx8YydOin6Z7b2r79JYwf30VTCQcJBb2IZOrgwQQGD/6YKVPW+tuKFo3i5Zfbcued9TVHTRBR0IvISdav38kNN7zH+vW7/G21apVjxowbufTS0h5WJmdDdzKISDrvvruWq64aly7kBw6sy/LlAxXyQUpH9CIC+EbV3Hffp4weHetvK1KkAG++2YlbbqnjYWWSXQp6EeHXX/fQo8cMvvtum7/t0ktLM3PmjVx5ZTkPK5OcoKAXCXNz5mzglltms3fvUX/bTTddzrhxnSlZspCHlUlOUdCLhKnExGSGD/+CF19c6m+LiorgP/9px913X6VRNSFEQS8ShrZu3U+vXu+zZMlmf1vlyucwY8aNNGhQwcPKJDco6EXCzJdf/kbPnjPZseOQv61jx2pMnNiN0qX15KdQpOGVImFkzJgY2rSZ5A/5iAjjmWda8tFHvRXyIUxH9CJhICkphQcemM8rr3zrbytXrhjTpvWgefMq3hUmeUJBLxLi9u49Sq9eM5k//xd/W7165fnww15UrKjHOIcDBb1ICIuL203nzlP56acTd7n26FGTCRO6UqxYQQ8rk7ykc/QiIWrRol9p0GBcupB//PFmTJ/eQyEfZnRELxKCxoyJ4Z575pGUlAJA4cIFeOedrvTqdYXHlYkXFPQiISQpKYX775/Pq6+euOhavnxxZs/upfHxYUxBLxIi9u49Ss+eM/nsM110lfQU9CIhIC5uN9dd9y4bNvzlb9NFVzlOF2NFgtwXX/guuqYNeV10lbR0RC8SxCZO/J6BA+fooqucloJeJAg553jhhW8YNmyhv00XXeVUFPQiQSYlxXH//fMZNWqFv61WrXJ88kkfXXSVTAV0jt7M2pvZBjOLM7NhmSyvYWbLzOyYmT2Ypr2SmS0ys/Vmts7MhuZk8SLh5tixJG6++YN0IX/ttRfy9de3KuTllLI8ojezSOB1oA0QD6w0sznOuR/TrLYbGAJ0y7B5EvCAc+47MysBxJrZ5xm2FZEA7N9/jOuvn87Chb/622644TImT76ewoX1y7mcWiBH9A2AOOfcJudcAjAN6Jp2BefcDufcSiAxQ/s259x3qd8fANYDOoEocoa2bz9I8+YT0oX83/8ezfTpPRTykqVAgr4CsCXN63jOIqzNrApQF1hxiuWDzCzGzGJ27tx5prsXCVlxcbtp3Hg8q1b96W8bMaIFr73WkchIjZCWrAXytySzB0e6M3kTMysOvA/c65zbn9k6zrmxzrlo51x02bJlz2T3IiErNvYPmjQZz6ZNewDfg0Leeqszw4c30zNdJWCB/M4XD1RK87oi8Eegb2BmUfhCfopz7oMzK08kfH3++S9cf/17HDyYAPjGyE+f3oMuXap7XJkEm0CO6FcC1cysqpkVBHoBcwLZufkOOd4G1jvn/nP2ZYqEl6lT19Kp07v+kD/33MIsWNBPIS9nJcsjeudckpndA8wHIoHxzrl1ZjY4dfloMzsfiAFKAilmdi9QE6gF9APWmtnq1F0+6pybmwt9EQkJI0cu4/77P/O/rlixJPPn96VmTZ3SlLMT0OX61GCem6FtdJrv/8R3SiejJWR+jl9EMnDOMWzYAl54Yam/7fLLy/Lpp301Rl6yReOyRPKBlBTHP/4xlzfeiPG3NW1amTlzenHuuUU8rExCgYJexGPJySkMGvQR48ev9rd16VKdadNuoEiRKA8rk1ChoBfxUFJSCrfcMpt3313rb+vd+womTuxOgQIaIy85Q0Ev4pGEhGT69Hmf999f728bMKAO48Z11o1QkqMU9CIeOHo0iRtvnMHHH//sb7vrrmhee60jEREavyA5S0EvkscOH06kW7dpfP75Jn/bffc15OWX2+puV8kV+v1QJA8dOHCMDh2mpAv5Rx9tqpCXXKUjepE8snfvUTp0mMLy5fH+tqeeasFjjzXzsCoJBwp6kTzw11+Hadt2Mt99t83f9tJLbXjggcYeViXhQkEvksu2bz9ImzaTWLt2h7/ttdc6cPfdDTysSsKJgl4kF23dup/WrSfx00+7ADCDceM6M3BgPY8rk3CioBfJJZs376Nly//xyy8n5pL/3/+60bdvLY8rk3CjoBfJBZs376N58wn8+uteAAoUiGDq1Bvo0aOmx5VJOFLQi+SwrVv307Ll//whX7BgJDNn3kjnzppLXryhoBfJQdu2HaBFixOna6KiIpg1qycdO1bzuDIJZ7phSiSHbN9+kJYtJ7Jx427Ad7rm/fdvUsiL5xT0Ijlg585DtGw50T+6JjLSeO+9HjpdI/mCgl4km3btOkyrVhP58cedgC/kp069ge7dL/O4MhEfBb1INuzefSTdzVAREcakSd258cbLPa5M5AQFvchZ2rv3KG3bTmL16j8B381QEyZ0pXfvKz2uTCQ9Bb3IWdi37yjt2k0mNvbE3DVvv92Ffv1qe1iVSOYU9CJn6PhUw99+u9XfNnbsdQwYUNfDqkROLaCgN7P2ZrbBzOLMbFgmy2uY2TIzO2ZmD57JtiLB5ODBBDp2fJdly05MNfzGGx254476HlYlcnpZBr2ZRQKvAx2AmkBvM8t4H/duYAjw0llsKxIUDh9OpHPnqSxZstnf9sor7bnrrqs8rEoka4Ec0TcA4pxzm5xzCcA0oGvaFZxzO5xzK4HEM91WJBgcOZJIly5T+fLL3/xtL7/cln/842rvihIJUCBBXwHYkuZ1fGpbIALe1swGmVmMmcXs3LkzwN2L5L6EhGR69JjBwoW/+tuef74199/fyMOqRAIXSNBn9iBLF+D+A97WOTfWORftnIsuW7ZsgLsXyV1JSSn07fsBc+du9LeNGNGChx5q4mFVImcmkKCPByqleV0R+CPA/WdnWxFPpaQ47rjjI2bM+NHfNnz4NQwfrme8SnAJJOhXAtXMrKqZFQR6AXMC3H92thXxjHOOoUPnMWHCan/bkCENeOqpFh5WJXJ2spym2DmXZGb3APOBSGC8c26dmQ1OXT7azM4HYoCSQIqZ3QvUdM7tz2zb3OqMSE4ZPvwLXnttpf/1bbfVYeTI9phldjZSJH8z5wI93Z53oqOjXUxMjNdlSJh69tnFPProF/7XPXtezpQp1xMZqfsLJf8ys1jnXHRmy/Q3VySNV19dkS7kr7vuUiZN6q6Ql6Cmv70iqd55ZxVDhnzqf92yZVVmzLiRqKhID6sSyT4FvQgwY8Y6br/9I//rRo0q8uGHvShcWE/blOCnoJew98knP9OnzwekpPiuV9Wpcz5z595M8eIFPa5MJGco6CWsLVr0Kzfc8B5JSSkA1KhRhs8+60upUoU9rkwk5yjoJWwtXx5P585TOXYsGYCqVUuxYEE/ypYt5nFlIjlLQS9h6fvv/6RDhykcOuSbh++CC0qwcGF/KlQo6XFlIjlPQS9hZ8OGXbRtO5m9e48CUKZMURYs6EfVqud6XJlI7lDQS1jZvHkfbdpMYseOQwCcc04hPvusL5ddpon0JHQp6CVsbN9+kNatJ7Jly34AihaNYu7cm6lbt7zHlYnkLgW9hIU9e47Qtu1kNm7cDUDBgpF8+GEvGjeulMWWIsFPQS8h7/hzXtes2Q5AZKQxbdoNtG59kceVieQNBb2EtKNHk+jWbRrLl594mPf48V3p3v0yD6sSyVsKeglZSUkp9O79frpHAL76agf696/tYVUieU9BLyEpJcVx220fMnv2T/62ESNacM89DTysSsQbCnoJOc45hgyZx6RJa/xt//xnYx599BoPqxLxjoJeQs5jj33B66+feDrUoEH1eP751no6lIQtBb2ElBde+IZnnlnif92r1xW88UYnhbyENQW9hIwxY2J4+OEF/tedOlVj4sRuejqUhD39C5CQMHXqWu666xP/6+bNq+jpUCKpFPQS9D7++Gf69ZvF8efcX3XVBcyZ04siRaK8LUwknwiJ56Q555g3L44lSzbz7bdb+eSTPhQqFBJdkywsWLCJHj3eIznZl/KXX16WefNupkSJQh5XJpJ/hEQamhlDh35KXJxvHpPY2G2awyQMLFmyma5dp/kfHHLRRefy2Wf9KF26qMeVieQvAZ26MbP2ZrbBzOLMbFgmy83MXkldvsbM6qVZdp+ZrTOzH8xsqpnlyjPamjQ5EexLlmzOjbeQfGTlyq107DiFw4d9Dw6pWLEkCxf254ILSnhcmUj+k2XQm1kk8DrQAagJ9DazmhlW6wBUS/0aBLyZum0FYAgQ7Zy7AogEeuVY9WmkDfpvvtmSG28h+cSaNdtp124yBw4kAFCuXDEWLuxPlSqlPK5MJH8K5Ii+ARDnnNvknEsApgFdM6zTFZjofJYDpczs+CTfBYAiZlYAKAr8kUO1p9OkSWX/90uXbsEdvzInIeWnn3bRuvVE9uzxPR2qdOkiLFjQn0svLe1xZSL5VyBBXwFIe4gcn9qW5TrOua3AS8BmYBuwzzn3WWZvYmaDzCzGzGJ27twZaP1+NWqU4dxzfWeFdu06zM8//3XG+5D87ZdfdtOq1UR27jwMHH86VD+uuOJvHlcmkr8FEvSZ3VKY8XA503XM7Fx8R/tVgQuAYmbWN7M3cc6Ndc5FO+eiy5Y988e6RURYuguwOn0TWrZs2UerVhP5448DABQrFsW8eTdTr56eDiWSlUCCPh5IO4SlIieffjnVOq2BX51zO51zicAHQOOzL/f0mjY9cfrmm290QTZU/PnnQVq1msjvv+8DoHDhAnz0UW8aNdLIKpFABBL0K4FqZlbVzAriu5g6J8M6c4D+qaNvGuI7RbMN3ymbhmZW1HyTjbQC1udg/enogmzo2bXrMK1bT/Q/AjAqKoJZs3rSokVVjysTCR5ZjqN3ziWZ2T3AfHyjZsY759aZ2eDU5aOBuUBHIA44DAxIXbbCzGYC3wFJwCpgbG50BCA6+gKioiJITExhw4a/2LnzEGXLFsutt5NctnfvUdq1m8y6db5rNpGRxvTpPWjf/hKPKxMJLpYfR6dER0e7mJiYs9q2UaO3/Y+Nmz27J1271sjJ0iSPHDyYQNu2k1i2zPdZmsGUKdfTu/eVHlcmkj+ZWaxzLjqzZSE3141O3wS/I0cS6dJlqj/kAd56q4tCXuQsKeglXzl2LIkbbniPRYt+87e9+moHbrutrndFiQS5kAv6tEMsY2L+4OjRJA+rkTNxPOTnzYvztz33XCs951Ukm0Iu6MuVK84ll5wHQEJCMrGxuXIjruSwo0eT6N59Op98stHf9n//14yHH27qYVUioSHkgh4yjqfX6Zv87siRRLp1m5buSP7RR5vyxBPNvStKJISEZNDrPH3wOHw4kS5dpjF//i/+tv/7v2aMGNFSz3kVySEhMR99RumDfjPOOYVGPnToUAJdukzjiy9+9bc98cS1/Otfzb0rSiQEheQRffXqZTjvvCIA/PXXETZs0ARn+c3Bgwl06vRuupB/6qkWCnmRXBCSQX/yBGea9yY/OXDgGB07TuGrr373tz3zTEsee6yZh1WJhK6QDHrQefr8av/+Y3ToMIXFi0/85/vCC6155JFrPKxKJLSF5Dl6UNDnR/v2HaV9+yn+KSoAXn65Lfff38jDqkRCX8ge0R+f4Azg5599E5yJd/buPUrbtpPThfx//9tOIS+SB0I26IsUiSI6+gL/66VLdVTvlT17jtCmzSS+/Xarv+3VVzswdGhDD6sSCR8hG/Sg0zf5we7dR2jdehIxMSfuUH7jjY6a1kAkD4V40OsOWS/9+edBWrb8H999t83fNmbMddx111UeViUSfkL2YixkPsFZ4cIh3eV8Y+PGv2jXbjK//roX8M0nP25cZwYOrOdxZSLhJ6SP6P/2t2JUq3ZigrO0pw8k93z77VYaNx7vD/mICGP8+K4KeRGPhHTQQ8bTN7pxKrfNnbuRFi3+x65dhwEoUqQAs2f35NZb63hcmUj4CoOg1wXZvDJhwmq6dJnK4cOJAJx3XhEWLuxP587VPa5MJLyFVdAvXbqF/PiM3GDnnOPpp79mwIAPSU72/XwvvPAcvvnmNho1qpTF1iKS20I+6GvUKEPp0prgLLckJ6dwzz1zeeyxRf622rXLsXTpQGrUKONhZSJyXMgHvZkmOMstR48mcdNNM3njjRh/W4sWVfjqq1u54IIS3hUmIukEFPRm1t7MNphZnJkNy2S5mdkrqcvXmFm9NMtKmdlMM/vJzNabWZ7f867z9Dlvz54jtG07iQ8+WO9v69nzcubNu5lzzinsYWUiklGWQW9mkcDrQAegJtDbzGpmWK0DUC31axDwZpplo4BPnXM1gNrAevJY2pE3S5boiD67tmzZxzXXvJNuBsp7772ad9+9gUKFdJ+CSH4TyBF9AyDOObfJOZcATAO6ZlinKzDR+SwHSplZeTMrCTQD3gZwziU45/bmYP0BiY6+gIIFIwHYuHE3O3ZogrOztW7dDho1ept163b62158sQ3/+U87IiL0FC+R/CiQoK8ApD3fEZ/aFsg6FwE7gXfMbJWZvWVmxTJ7EzMbZGYxZhazc+fOzFY5a4ULF6B+/fL+15rg7Ox89dVvNG36Dlu3HgCgQIEIJk/uzoMPNtajGkXysUCCPrN/wRnHKJ5qnQJAPeBN51xd4BBw0jl+AOfcWOdctHMuumzZsgGUdWYyPkdWAuecY9So5bRuPYm9e48CULx4QebO7cPNN9fyuDoRyUogQR8PpB0MXRHIOJfAqdaJB+KdcytS22fiC/48pwnOzs7Bgwn06fMB9947n6SkFMA3tcRXX91KmzYXe1ydiAQikKBfCVQzs6pmVhDoBczJsM4coH/q6JuGwD7n3Dbn3J/AFjM7fmtkK+DHnCr+TKQ9oo+N3cbRo0lelBFUfv75Lxo2fItp037wt0VHX8DKlXdQr17502wpIvlJlkHvnEsC7gHm4xsx855zbp2ZDTazwamrzQU2AXHAOODvaXbxD2CKma0B6gDP5GD9AStbthiXXloa0ARngZg1az3R0WPTXXQdNKgeixcPoHLlczysTETOVEBj4Zxzc/GFedq20Wm+d8Ddp9h2NRCdjRpzTJMmlfj5Z9+dsUuWbKZp08pZbBF+kpJSePzxL3juuW/8bYUKRfLmm50YMKCuh5WJyNkK+Ttj09KNU6e3c+ch2refnC7kq1QpxdKlAxXyIkEsrO5uSXtBdunSLaSkOI39TvXtt1u54Yb3iI/f729r3/4SJk/uTunSRT2sTESyK6yO6KtXL+2f4Gz37iNs2LDL44q855xjzJgYrrnmnXQh/69/XcvHH/dWyIuEgLAK+pMnOAvv0zdHjiRy221zGDz4ExISkgEoVaowH3/cmyeeaE5kZFj99RAJWWH3L1nn6X02bNhFkybjmTBhtb+tdu1yxMTcQadOl3pYmYjktLAL+rQjbcLxDtnk5BRefPEb6tQZw6pVf/rb+/evzdKlA7n44vM8rE5EckNYXYwFqF/fN8FZQkKyf4Kzv/0t0+l3Qs769TsZMOBDVqzY6m+Liopg1Kj2DB4crflqREJU2B3RFy5cgOjoC/yvw+GoPikpheeeW0LdumPShXzduuezcuUd3HXXVQp5kRAWdkEP4XWe/tVXHQwAAAmFSURBVIcffNMKP/LIQo4d811wjYqK4KmnWrBixe3Urn2+xxWKSG5T0Ido0CcmJvP0019Tr96YdNM91K9fntjYQTz2WDOioiI9rFBE8krYnaMH0g2xjI39gyNHEilSJMrDinLWmjXbufXW2ekuthYsGMmTTzbnwQcbU6BAWP7/LhK2wvJffNoJzhITU0JmgrOEhGSefPJL6tcfmy7kr766AqtW3cmwYU0V8iJhKGz/1Yfa6ZtVq7bRoME4nnjiK/+88YUKRfLii2345pvbqFkz5x/mIiLBIWyDPv14+uAN+q1b9zNw4IdER4/j+++3+9sbNarI6tWDefDBxrrDVSTMhW0CpD2i/+yzX5gzZ4OH1Zy5AweO8fjjX1Ct2quMH7+alBTf0x2LFCnAf/7TlsWLB1CjRhmPqxSR/CBsg/7SS0tTp45vaGFCQjLXXz893ZOU8qvExGTeeGMlF1/8CiNGLObIkRNPymrf/hK+/34w993XSEfxIuIXtmlgZsya1ZOLLz4XgORkR58+7/P22995XFnmnHPMnv0TV1zxJnffPZedOw/7l9Wpcz4LFvRj3rybqVattIdVikh+FLZBD76HaixePMB/odI5uP32jxg1arnHlaW3fHk8zZpNoHv36f4nZAFUqlSSiRO7ERs7iFatLvKwQhHJz8I66AHKly/BV1/dmu5h1/feO5+nn/4a3xMSvfPLL7u56aYZNGr0NkuWnJiq4ZxzCvH8863ZsOEe+vWrrYeniMhphX3QA5QpU5Qvvuif7gLtY48t4pFHFnoS9tu3H+Teez/lssteZ8aMH/3tUVERDB16NXFxQ3jooSYhdZOXiOSesLwzNjPnnFOY+fP70q3bdBYs2ATA889/w8GDCbzySodcP2pOSXEsWvQrY8bEMmvWT/6x8MfdeGNNnn22laYRFpEzpqBPo1ixgnz0UW969pzpH275+usrOXQokXHjOufKXaW7dh1mwoTVjB0by8aNu09a3qRJJV56qS0NG1bM8fcWkfCgoM+gcOECzJx5I/37z/YPt5wwYTUHDyYwZcr1FCyY/YnAnHMsXryZMWNimTnzR/9j/NJq3LgS//xnY7p2ra4phEUkWwIKejNrD4wCIoG3nHPPZVhuqcs7AoeBW51z36VZHgnEAFudc9flUO25JioqksmTu1O8eBRvvbUKgJkzf+Tw4URmzrzxrM+N79lzhIkTv2fMmFjWrz/5weQlSxaif/9aDBpUnyuvLJetPoiIHJdl0KeG9OtAGyAeWGlmc5xzP6ZZrQNQLfXrauDN1D+PGwqsB0rmUN25LjIygrFjO1OsWEFGjVoBwNy5G+nU6V0+/LAXJUoUOuW2iYnJHDiQwIEDx9i//xjbtx9i8uQ1TJ++jqNHk05av0GDCtx5Z3169rycYsUK5lqfRCQ8BXJE3wCIc85tAjCzaUBXIG3QdwUmOt8QleVmVsrMyjvntplZRaAT8DRwf86Wn7vMjJEj21GiREFGjFgMwKJFv9G06TtcdlkZDhxIYP/+Yxw4cCxdsB9/wMfpFC9ekJtvvpI776xP3brls1xfRORsBRL0FYC0s37Fk/5o/VTrVAC2Af8FHgJKnO5NzGwQMAigcuXKp1s1T5kZTz3VkhIlCvHwwwsA33zva9Zsz2LLzNWtez6DB0fTu/cVp/2tQEQkpwQS9JldCcw4uDzTdczsOmCHcy7WzJqf7k2cc2OBsQDR0dHe3qmUiYceakLx4gW5++65Wa4bEWGUKFGQEiUK+f+sU6cct99ej+joC3RxVUTyVCBBHw9USvO6IpDxSR2nWqcH0MXMOgKFgZJmNtk51/fsS/bO3/9+FfXrl2fNmu0UL54+yEuUKEjJkoUoUaIQRYoUUJiLSL5hWd35aWYFgJ+BVsBWYCXQxzm3Ls06nYB78I26uRp4xTnXIMN+mgMPBjLqJjo62sXExJxZT0REwpiZxTrnojNbluURvXMuyczuAebjG1453jm3zswGpy4fDczFF/Jx+IZXDsip4kVEJHuyPKL3go7oRUTOzOmO6DWpmYhIiFPQi4iEOAW9iEiIU9CLiIQ4Bb2ISIjLl6NuzGwn8HuapjLAydM9BrdQ61Oo9QdCr0+h1h8IvT5lpz8XOufKZrYgXwZ9RmYWc6phQ8Eq1PoUav2B0OtTqPUHQq9PudUfnboREQlxCnoRkRAXLEE/1usCckGo9SnU+gOh16dQ6w+EXp9ypT9BcY5eRETOXrAc0YuIyFlS0IuIhDjPg97M2pvZBjOLM7NhmSw3M3sldfkaM6sX6LZeyGZ/fjOztWa22szyzfSdAfSphpktM7NjZvbgmWzrhWz2J1g/o5tT/76tMbOlZlY70G29kM3+BOtn1DW1P6vNLMbMmga6bZacc5594Zvf/hfgIqAg8D1QM8M6HYF5+B5X2BBYEei2wdSf1GW/AWW87MNZ9ulvwFX4HgD/4JlsG0z9CfLPqDFwbur3HULg31Gm/Qnyz6g4J66b1gJ+yqnPyOsj+gZAnHNuk3MuAZgGdM2wTldgovNZDpQys/IBbpvXstOf/CrLPjnndjjnVgKJZ7qtB7LTn/wqkD4tdc7tSX25HN/jPgPa1gPZ6U9+FUifDrrUZAeKceLZ3Nn+jLwO+grAljSv41PbAlknkG3zWnb6A74P9jMzizWzQblW5ZnJzs85WD+j0wmFz2ggvt8qz2bbvJCd/kAQf0Zm1t3MfgI+AW47k21PJ5CHg+emzJ6gnXG856nWCWTbvJad/gA0cc79YWZ/Az43s5+cc1/naIVnLjs/52D9jE4nqD8jM2uBLxiPn/8N6s8ok/5AEH9GzrlZwCwzawY8BbQOdNvT8fqIPh6olOZ1ReCPANcJZNu8lp3+4Jw7/ucOYBa+X9m8lp2fc7B+RqcUzJ+RmdUC3gK6Ouf+OpNt81h2+hPUn9Fxqf8xXWxmZc5021Pt0MsLFAWATUBVTlxkuDzDOp1If/Hy20C3DbL+FANKpPl+KdDey/6c6c8ZeIL0F2OD8jM6TX+C9jMCKgNxQOOz/XkESX+C+TO6hBMXY+sBW1NzItufkaedT+1QR+BnfFeVh6e2DQYGp35vwOupy9cC0afb1uuvs+0Pvivq36d+rcsv/QmwT+fjO+rYD+xN/b5kEH9GmfYnyD+jt4A9wOrUr5jTbev119n2J8g/o4dTa14NLAOa5tRnpCkQRERCnNfn6EVEJJcp6EVEQpyCXkQkxCnoRURCnIJeRCTEKehFREKcgl5EJMT9PxifalHytdrNAAAAAElFTkSuQmCC\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.plot(h, Int_bias_h + Int_var_h, color = col_navy, lw=3)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Cross-validation\n",
    "\n",
    "Estimate \n",
    "\\begin{equation}\n",
    "\\hat{f}_{h,-i\\,}(x_i) = \\sum_{j\\ne i} \\frac{K_h(x_i-x_j)}{{\\sum_{k\\ne i} K_h(x_i-x_k)}}\\,y_j\n",
    "\\end{equation}"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "metadata": {},
   "outputs": [],
   "source": [
    "r.seed(201053)      # fix seed to have bias and variance calculated on the same random set\n",
    "\n",
    "for j in range(J_cv):\n",
    "    y = gendata(f_x, sigma, n)\n",
    "    for i in range( len(h) ):\n",
    "        f_l = 0\n",
    "        for l in range( len(x) ):\n",
    "            x_l = np.delete(x, l,  axis = 0)\n",
    "            y_l = np.delete(y, l,  axis = 0)\n",
    "            u = fh([x[l]], x_l, y_l, h[i]) \n",
    "            f_l = f_l + (y[l] - u)**2 \n",
    "          \n",
    "            # end l loop over obs\n",
    "        CV[i] = f_l\n",
    "        \n",
    "    # end h loop\n",
    "    \n",
    "        CV_all[j, i] = CV[i]/n\n",
    "# end J loop\n",
    "CV_mean = np.mean( CV_all, axis = 0 )"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[<matplotlib.lines.Line2D at 0x7fc04e04acd0>]"
      ]
     },
     "execution_count": 19,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAXoAAAD4CAYAAADiry33AAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADh0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uMy4yLjIsIGh0dHA6Ly9tYXRwbG90bGliLm9yZy+WH4yJAAAgAElEQVR4nO3da3RV5b3v8e8/CYQQSAKEexISuYoKFQNYtVp31Xprqbt2V+3l7F4223br7hlndI92jPPmjHaccXbHOGecXo6t5Vhb666b01axtl7Qumtta8WEuyBIIAFC0EAkQEhCbs958ay4VpIFWUnWylxrrt9nDIZrPnNOeCZLfkye+cz/Y845REQkvHKC7oCIiKSWgl5EJOQU9CIiIaegFxEJOQW9iEjI5QXdgXhKS0tdZWVl0N0QEckYW7duPemcmxlvX1oGfWVlJbW1tUF3Q0QkY5jZ4Qvt09CNiEjIKehFREJOQS8iEnIKehGRkFPQi4iEnIJeRCTkFPQiIumgpR1eqPP/TbK0nEcvIpIVznXBtuPwxjE4eMq3OQe3LU7qL6OgFxEZT929sLvZh/ueZugdtCbIG8fg1kVglrRfUkEvIpJqfQ4OtPgQ3/4OdPYMPSbHYFkprJkPDkhezivoRURSwjk4dtaHe20TtHbGP66yBFbPg6vmQVF+SrqioBcRSaa2LtjSCH9thKaz8Y+ZORlWz/cBP3tKyrukoBcRGSvn4NAp+NMR/3C1p2/oMVMmwlVz/dBMZUlSx+CHo6AXERmtjm4/NPOnI/Hv3ifkwMo5/s59+UzIDWZGu4JeRGSkjpyGPx32Y+/ne4furyiGD1X4cfdJwcds8D0QEckE53tg63Ef8IdPD90/MdffuV9XAQtKxr9/F6GgFxG5mOZz8EqDf8DaEWda5Lyp/u59zXwomDDu3UuEgl5EJJ6GVnjpIOx4x89rj5WXA6vm+oC/ZNq4PlgdDQW9iEi/PuffVv39ITjw3tD9swr90MzVZX4WTYZQ0IuI9PRBzTEf8Mfbhu5fPhM+UgVLS/0brBlGQS8i2aujG/58BP7QMPTN1RyD6nlw0yVQVhRI95JFQS8i2ae1E/5Q7+e/D647k5/rh2durILpBcH0L8kU9CKSPU62w/MH/EtOg6tGFuXDjZXwoQUwOT1nz4yWgl5Ewq+10wf8a0eHBvzsQj88s2Y+TMgNpn8ppqAXkfA6ex42H4RXDw+tP7NwGty8EC6flZEPWEdCQS8i4dPe7efAv9IwtETBounwsSWweEYgXQuCgl5EwqOzxz9k/f2hoW+xLiiGjy2FS0vT/gWnZFPQi0jm6+qFPzbAiwfhXPfAffOm+jv4FbOzLuD7KehFJHN19/oHrC/UwenzA/fNKoQ7l/hSBSEfgx+Ogl5EMo9zsOtdeOotONE+cN+MArh9sZ9FE1D993STUNCb2a3A94Bc4BHn3L8O2v9h4DdAfaTpKefctxI5V0RkRJrOwq/3wr6TA9uL8+G2xXBNuS86Ju8bNujNLBd4CLgZaARqzOwZ59zeQYf+yTl35yjPFRG5uHNd8OwBP1WyL2Yu/OQJcOsiuH6BrwkvQyRyR78GqHPOHQIws43AOiCRsB7LuSIi0Nvn69H87u2BD1oN/xbrnUsyqpJkEBIJ+vnA0ZjtRmBtnOM+aGY7gSbg6865PSM4FzNbD6wHqKioSKBbIhJ6+0/Cr/YOXY91yQz41HKYn9nFxsZLIkEf73H14DL824AFzrk2M7sdeBpYnOC5vtG5DcAGgOrq6rjHiEiWONnuH7TueGdg+4wC+NtL4QNzsnaq5GgkEvSNQHnMdhn+rv19zrkzMZ+fM7MfmllpIueKiLyvswc218HL9QNLFkzM9ePwH6kKbT2aVEok6GuAxWZWBRwD7gHuiz3AzOYA7zrnnJmtAXKAFqB1uHNFRHAOtr8Dv9ozdD78mvnwiWVQMimYvoXAsEHvnOsxsweAzfgpko865/aY2f2R/Q8DdwNfMbMeoAO4xznngLjnpuhaRCQTtXbCxjf9vPhYC4rhU5f5NVllTMzncXqprq52tbW1QXdDRFKpz/m3Wp96a+DiH0X5sG4prC3L+jdaR8LMtjrnquPt05uxIjL+ms/BL3YNXYD72nK469LQLfwRNAW9iIyf3j7/oPXZt6E75mHrzMnwmRV+2qQknYJeRMbH0dPwb7vg6JloW4751Z1uX6y3WlNIQS8iqdXVC88d8DXiY0sXlBf5u/iK4uD6liUU9CKSOm+3wBO7/Zh8vwk5cMcSPyde1SXHhYJeRJKvs8fPpvnzkYHti6f7u/hZhcH0K0sp6EUkuRpa4afbB9aJn5TnSxdcU64pkwFQ0ItIcvQ5X77g2QMDx+JXzIZ7LtebrQFS0IvI2LW0w2M7oS5mXvykPPj0Zb6EgQqQBUpBLyJjU3PMlzDoiHm79ZJp8PcfgNLJwfVL3qegF5HR6eiG/7cH3jgWbcsxPyf+ows1oyaNKOhFZOQOvgc/2wEtHdG2GQXwhStVhCwNKehFJHG9ffB8HTx/YOASQleX+RWfClSjJh0p6EUkMSfb/bTJ+tZoW0Ee3HsFVM8Lrl8yLAW9iAxvS6N/4Hq+N9q2aLp/4Dq9ILh+SUIU9CJyYd298Ou98KeYN1xzDO5cArcs1MtPGUJBLyLxvdcBj2zzb7r2m1Xo7+IrS4Lrl4yYgl5Ehtp3En6yDc51R9tWzYXPrvAvQklG0TcmIlF9Dl48CL/dH51Vk2O+Ts2NlXrDNUMp6EXEa++Gx3bA7uZoW1E+fHmVf/AqGUtBLyLQeAb+79aBFScXTYcvXQnFKkaW6RT0Itnu9Ub4990D13D9SBV8YpnKGISEgl4kW8WbOpmfC59b6R+8Smgo6EWyUbypk3OmwD+sgrlTg+uXpISCXiTb7D8JP9kObV3RNk2dDDV9qyLZ5M9HfCmD/hWgNHUyKyjoRbJBn4NNb8HL9dE2TZ3MGgp6kbDr7IFHt8ObMfPjy4vg/mqYpoJk2UBBLxJmLe3wo1poOhttWznb16vJ1x//bKFvWiSs6k/Bw7VwNuah6y0L4eNLVXUyyyjoRcKo5hg8vgt6Ii9B5RrcdwV8sDzYfkkgFPQiYeIcPHsAnjsQbSucAP9YrYeuWUxBLxIWXb3w+E7YejzaNmcKfKUaZhYG1y8JXEKFLMzsVjPbb2Z1ZvbNixy32sx6zezumLYGM9ttZjvMrDYZnRaRQU53wndfHxjyl5bCv1yjkJfh7+jNLBd4CLgZaARqzOwZ59zeOMd9B9gc56e50Tl3Mgn9FZHBGs/Aj2rgVGe07YYFcPdyFSUTILGhmzVAnXPuEICZbQTWAXsHHfcg8CSwOqk9FJEL29Psa9b0L9ptwKcugw9XBtkrSTOJ/HU/Hzgas90YaXufmc0H7gIejnO+A140s61mtn60HRWRQV5v9HPk+0N+Uh58dbVCXoZI5I4+3oRbN2j7u8A3nHO9NrRexrXOuSYzmwW8ZGb7nHOvDvlF/F8C6wEqKioS6JZIlnIOXjoET++Ltk0v8CE/T5UnZahEgr4RiJ18WwY0DTqmGtgYCflS4HYz63HOPe2cawJwzjWb2Sb8UNCQoHfObQA2AFRXVw/+i0REwNes+fVeeKUh2jZ/KvzTGijRSlASXyJBXwMsNrMq4BhwD3Bf7AHOuar+z2b2M+B3zrmnzawQyHHOnY18vgX4VrI6L5JVunvhsZ2wLWZmzeLpvmZNwYTg+iVpb9igd871mNkD+Nk0ucCjzrk9ZnZ/ZH+8cfl+s4FNkTv9POAJ59wLY++2SJbp6PblDA68F21bNRf+00qYkBtcvyQjmHPpN0pSXV3tams15V4EgNZOeOgNOBZTmOzDlX76pGrWSISZbXXOVcfbpzdjRdLZu23wgzf80n/9PrEMbr5EC4VIwhT0Iumq/hT8sAbOdfvtHPPL/V1dFmy/JOMo6EXS0e53/YtQ3ZHqkxNz/cLdl80Ktl+SkRT0IunmtaPwxO7ouq5TJvo58pUlwfZLMpaCXiRdOAcv1MFv3462zSiAB9fCLBUmk9FT0IukA+fgqUGLd5cX+Tv5Yr0IJWOjoBcJWp+Df98Nf4kpKbV0Bqy/Si9CSVIo6EWC1NsHP9sxsI78ytnwxSv1IpQkjYJeJChdvX5mzZvN0ba18/0UStWRlyRS0IsEobPHLxYSW9LghgW+lrzedpUkU9CLjLdzXfBQDTS0Rts+uhA+vlRvu0pKKOhFxtPpTl/SoCmmbs0nlsEtC4Prk4Segl5kvLS0w/e3wIl2v23Apy+H6xcE2i0JPwW9yHh4tw2+t8VXogQ/Dv+5FbBWdWsk9RT0Iql29LQfrmnr8tt5OfClK2HlnGD7JVlDQS+SSodO+VryHT1+e2KuXxFqWWmw/ZKsoqAXSZV9J+HHtXC+128X5PmSBgunB9svyToKepFUeLMZNmyFnkiZ4SkT4cE1UF4cbL8kKynoRZJt5zv+jdfeSJnhkknwtbUwe0qw/ZKspaAXSaatTfDTHdFa8jMK4GtXQ+nkYPslWU1BL5IsWxrh5zshkvHMnOxDfnpBoN0SUdCLJMNrR+EXu6IhP7vQh3yJaslL8BT0ImP16mHY+GZ0e95UPyY/NT+4PonEUNCLjMXLh+DJt6Lb5UV+6b8pE4Prk8ggCnqR0dpcB7/ZH92uLIEH1sBkrQol6UVBLzJSzsFzB+DZA9G2hdP8y1Ba+k/SkIJeZCScg2f2w+aD0bYlM+Ar1ZCvP06SnvR/pkiinPPj8f9RH21bPtMv4j1R67tK+lLQiySiz8Ev9/gZNv2umAVfXqVFvCXtKehFhtPn4Indfq58vyvnwBeu9CWHRdKcgl7kYvocPL4TthyLtq2eB59fCbkKeckMCnqRC+ntg8d2Qm1TtO3qMvjsCr9ClEiGUNCLxNPTBz/dDtvfibZdWw73XqGQl4yjoBcZrLsXfrIddr0bbbthAXzqMoW8ZKSEBhnN7FYz229mdWb2zYsct9rMes3s7pGeK5IWunr9giGxIf83VfB3CnnJXMMGvZnlAg8BtwHLgXvNbPkFjvsOsHmk54qkha5eeLgW9pyItt2yED55KZhCXjJXInf0a4A659wh51wXsBFYF+e4B4EngeZRnCsSrM4ev4j3vpPRttsXw7qlCnnJeIkE/XwgZgIxjZG295nZfOAu4OGRnhvzc6w3s1ozqz1x4kS8Q0RSo6Pbh/yB96JtH1sCdy5RyEsoJBL08f5Pd4O2vwt8wznXO4pzfaNzG5xz1c656pkzZybQLZEkaO+GH7wBB09F2+5aBrctDq5PIkmWyKybRqA8ZrsMaBp0TDWw0fzdTylwu5n1JHiuSDDauuAHW+DomWjb3cv9w1eREEkk6GuAxWZWBRwD7gHuiz3AOff+nwwz+xnwO+fc02aWN9y5IoE4ex6+vwWOnY223XM5XL8guD6JpMiwQe+c6zGzB/CzaXKBR51ze8zs/sj+wePyw56bnK6LjNLpTh/yx9v8tgH3XQHXVgTaLZFUMefiDpkHqrq62tXW1gbdDQmjlnYf8ifa/bbh69asLQu0WyJjZWZbnXPV8fbpzVjJHu+2wfe2QGun384x+PsPQPW8YPslkmIKeskOjWfg/7wBZ8777bwcX0t+xexg+yUyDhT0En4NrT7k27v99sRcuL8alpUG2y+RcaKgl3A70AI/rIHzkVc8JuXBP62GhdOD7ZfIOFLQS3jtafYFyrr7/HbhBHhwLVQUB9svkXGmoJdw2n4cHt0OvZFZZcX5PuTnTQ22XyIBUNBL+GxphMd3+WUAAaYXwD+vhVmFwfZLJCAKegmXVw/Dxjej27MKfchPLwiuTyIBU9BLeLx0EDbti27Pm+pDvig/uD6JpAEFvWQ+5+DZA/DcgWhbZYmfXVM4Mbh+iaQJBb1ktj4HT+6FPzRE2xZNh6+u9lMpRURBLxmsuxce2wnbjkfbls+E9Vf5l6JEBFDQS6Zq74Yf1w5cFerKOb52zQSFvEgsBb1knlMdvqRBf5lhgA9X+kVDcrT0n8hgCnrJLE1nfcj3V6AEv/TfTZdofVeRC1DQS+Z4u8UP13T0+O1cg8+thDVx15sXkQgFvWSGbcfhZzugJ1K3ZlKef+iqCpQiw1LQS/r7j3o/hbJ/MbSifD9HvlzFyUQSoaCX9NXn4Ol98PtD0bbZhfDAGpgxObh+iWQYBb2kp+5eX5istinadsk0v2DIFL3tKjISCnpJPx3dvo78/pZo28rZ8IUr9SKUyCiEK+i7euHoaX/np6l2mem9DvhRDRw7G237UAV8+nLNkRcZpfAE/f/+Kxw85cd1v32jxnAz0dst8Mg2aOuKtn18KXx0of7iFhmD8AR9jkUXmqhvVdBnEud8UbKn3op+hzkGn10BV5cF2jWRMMgJugNJU1kS/dzQGlw/ZGS6ev38+F/vjYb81Inwn69WyIskSXju6KumRT/XnwquH5K4k+3+oWvjmWhbZYl/EapkUnD9EgmZEAV9zB390TP+Dcq88PyDJXT2noCfbodz3dG2a8vh7y5T9UmRJAtP0E/NhxkF0NLhQ77xzMDhHEkPzsGLB+GZ/dE3XfNyfMBfVxFo10TCKjxBD374pqXDf64/paBPN5098PhO2P5OtK043w/VxA69iUhShSvoK0uib1LqgWx6aT7nK0/G1pBfOA2+vAqKNR4vkkrhCvoqzbxJS7vfhZ/u8Hf0/W5YAJ9crucoIuMgXEFfVuSDo6cPTrT7F29UFyU4vX3wfB08dyDaNiEH7r1CUydFxlG4gn5Crg/7/rv5+lNwxexg+5Stms7Cz3fCkdPRtukFfjy+QuWFRcZTQv9uNrNbzWy/mdWZ2Tfj7F9nZrvMbIeZ1ZrZdTH7Gsxsd/++ZHY+Lg3fBKu3D16og3/988CQXzoDvnGtQl4kAMPe0ZtZLvAQcDPQCNSY2TPOub0xh70MPOOcc2a2AvglsCxm/43OuZNJ7PeFxc60qVfQj6vjkbv4wzEBn5cDdy6Bj1RBrsbjRYKQyNDNGqDOOXcIwMw2AuuA94PeORczlYJCojOkx1/sNL2GVv9avaoeplaf84uD/O7t6FJ/AAuK4fMrYe7U4PomIgkF/XzgaMx2I7B28EFmdhfwP4BZwB0xuxzwopk54MfOuQ2j724CZhT4B7BtXX6Wx7ttCppUeqfN38XHDpPlGtyxBG6+RHfxImkgkaCPdzs85I7dObcJ2GRm1wPfBm6K7LrWOddkZrOAl8xsn3Pu1SG/iNl6YD1ARcUY3pA08+P0u5v9dkOrgj4V+hy8fAh+O+guviJyFz9Pv+ci6SKR261GoDxmuwxousCxREJ8oZmVRrabIv9tBjbhh4LinbfBOVftnKueOXNmgt2/gAEFzjROn3TvtsH/eg027YuGfK7Bx5bAv1yjkBdJM4nc0dcAi82sCjgG3APcF3uAmS0CDkYexq4CJgItZlYI5DjnzkY+3wJ8K6lXEI9KFqdGn4M/1Ps6Nd0xd/HlRfC5lX5qq4iknWGD3jnXY2YPAJuBXOBR59weM7s/sv9h4JPA582sG+gAPh0J/dn44Zz+X+sJ59wLKbqWqAXFfsDJAcfO+LH6SeF6ZWDcvXUCnt7nK4P2yzG4fbFfAUpj8SJpK6H0c849Bzw3qO3hmM/fAb4T57xDwMox9nHkCibAnCm+rorDz+deMmPcuxEKR0/7IZp9g2bHlhX5sXjdxYukvfDe5lZNixbQamhV0I9US7t/0PrGsYHtE3Lgo4t0Fy+SQUIc9CXwWmRWqFacSlxbF2yugz8eHjibxoBryv20Sa3+JJJRwhv0g9+Qdc5PvZT4unrhlQYf8h09A/etmA3rlmqaqkiGCm/Qz50K+blwvhfOnIdTnb6olgzU52BLox+mae0cuK+qBO66FBZND6ZvIpIU4Q36HIMFJfB2i99uaFXQx3IO3myG3+z3lSZjzSr0d/AfmKN/BYmEQHiDHvwdaX/Q15+CVXOD7U866HOw4x2/bmtsdUmAonw/XfLacj1oFQmRkAe93pB9X0+fn0Hz4kG/rF+s/Fy4eSH8TZXeNxAJoXD/qY59IHv0tA+7bFu6rrMH/nIEXq4fOgY/IQeurYBbF/m7eREJpXAHfVG+r2bZ0uFf2T92xo/bZ4O2Lvhjg59Jc6574L6CPLh+gb+Dn6qAFwm7cAc9+Lv6lg7/uaE1/EHf2ulrw//liJ9xFKso34f7hyr828MikhXCH/RV02Drcf+5vhVuCLY7KdN8zo+/b2mE3kFVpGcU+DH4D5b5dXVFJKtkQdDHvjgVwjdkT3XAcwfgr41+Rk2seVN9qYJVczWLRiSLhT/oy4p8rfReByfa/dj1lIlB92rsznXB5oN+HD62ZDDAJdN8wF8+S/PgRSQLgn5CLpQXR+vSN7T6AMxU53vgDw3w0sGhpQqWzvC1aPQmq4jECH/Qg38gm+lB39PnH7A+X+dLOsSqKIZPLINlpcH0TUTSWnYEfVUJvBL5nGnj9H0Otjb5WjQn2wfum1UIH18KV6pUgYhcWJYEfcwbsg2tPjxz0jwYnYO9J3wtmsYzA/eVTII7FsPVZXrIKiLDyo6gn1HgH8C2dflx7eZzfgWqdHX0NPxqL9S9N7B98gT/kPWGSpioaZIikpjsCHozP07/ZrPfrj+VnkHvnH+T9am3Bs6Fn5gLN1b6ufCT9aKTiIxMdgQ9+HH6/qBvaIUPlgfbn8E6uuHfdsH2d6JtOQbXVcBti6BYqzqJyOhkUdCncSXLI6fhkW0DH7aWF8EXr4TZafgvDxHJKNkT9AuK/bqnDl/c7HwP5Ad8+c75tVmfemvg+qw3LIC/vVTlCkQkKbIn6Asm+HH5420+7I+chsUzgutPvKGaSXnw2RVaIEVEkip7gh78A9njbf5zfWtwQX+hoZovrfJz40VEkii7gr5qmi/+BdE3ZceTc/DqYXhy0FDN9QvgkxqqEZHUyK6grxxUydK58XujtKMbfrEbth2Ptk3Kg89cAVfNG58+iEhWyq6gnzfVr496vhdOn4dTnTC9IPW/7pHT8JNtvnpmv7Ii+LKGakQk9bIr6HPMFwA7EHnjtKE19UG/pdHfyccO1XyoAu5erqEaERkX2VcoZcB8+hQXOHvpIDy2Mxryk/L83Ph7r1DIi8i4ya47ehg4Tp+qB7J9Dja9BS/XR9vmToF/rNZQjYiMu+wL+tilBY+cht6+5FaA7OmDx3dCTVO0bdF0uL9adWpEJBDZN3RTPCk6Lt/dB8fOJu/n7uyBH9UMDPmVs+GBNQp5EQlM9gU9DJ1mmQxnz8N3X4e3TkbbrquAf7hKJYVFJFDZGfSxwzfJKHB2sh3+52t+KKjfHYvh3svTf4ETEQm9hILezG41s/1mVmdm34yzf52Z7TKzHWZWa2bXJXpuIJL5QPboaR/y/XPkDbjncr9It5b3E5E0MOzDWDPLBR4CbgYagRoze8Y5tzfmsJeBZ5xzzsxWAL8EliV47vgrL4Zc84t7NJ/zK09NmTjyn2f/SfjxVj82D5CX46dPfmBOcvsrIjIGidzRrwHqnHOHnHNdwEZgXewBzrk251z/kkiF+PqQCZ0biIm5ML8oun14FHf1247DQzXRkC/IgwfXKORFJO0kEvTzgaMx242RtgHM7C4z2wc8C3xxJOdGzl8fGfapPXHiRCJ9H5uxjNO/0uBLGvS/CFWcD//lg8GWPRYRuYBE5tHHG2h2Qxqc2wRsMrPrgW8DNyV6buT8DcAGgOrq6rjHJFXVNL/oB8DrjXDiXHRM3fCfY/9L5PO5roE15GcX+umTMyanvMsiIqORSNA3ArELrJYBTRc4Fufcq2a20MxKR3ruuIp9IPteh/8xmp/jq6tHN74vIjJOEhm6qQEWm1mVmU0E7gGeiT3AzBaZ+dthM1sFTARaEjk3MDMnDwz7kbpsJnxtrUJeRNLesHf0zrkeM3sA2AzkAo865/aY2f2R/Q8DnwQ+b2bdQAfw6cjD2bjnpuhaRsYM/nkt7DsJXb2+Nn3/gFH/5wu1lU6GS2dqjryIZASLTpZJH9XV1a62tjboboiIZAwz2+qcq463LzvfjBURySIKehGRkFPQi4iEnIJeRCTkFPQiIiGnoBcRCTkFvYhIyKXlPHozOwEcjmkqBU5e4PBMFbZrCtv1QPiuKWzXA+G7prFczwLn3Mx4O9Iy6Aczs9oLvQiQqcJ2TWG7HgjfNYXteiB815Sq69HQjYhIyCnoRURCLlOCfkPQHUiBsF1T2K4HwndNYbseCN81peR6MmKMXkRERi9T7uhFRGSUFPQiIiEXeNCb2a1mtt/M6szsm3H2m5l9P7J/V2QFq4TODcIYr6fBzHab2Q4zS5uC/Alc0zIz+6uZnTezr4/k3CCM8Xoy9Tv6TOT/t11m9pqZrUz03CCM8Xoy9TtaF7meHWZWa2bXJXrusJxzgf3Arzp1ELgEv/zgTmD5oGNuB57HL9F9NbAl0XMz6Xoi+xqA0iCvYZTXNAtYDfx34OsjOTeTrifDv6NrgGmRz7eF4M9R3OvJ8O9oCtHnpiuAfcn6joK+o18D1DnnDjnnuoCNwLpBx6wDfu6814ESM5ub4LnjbSzXk66GvSbnXLNzrgboHum5ARjL9aSrRK7pNefcqcjm60BZoucGYCzXk64SuaY2F0l2oJDoQqZj/o6CDvr5wNGY7cZIWyLHJHLueBvL9YD/Yl80s61mtj5lvRyZsfw+Z+p3dDFh+I6+hP9X5WjOHQ9juR7I4O/IzO4ys33As8AXR3LuxQy7OHiKxVtde/B8zwsdk8i5420s1wNwrXOuycxmAS+Z2T7n3KtJ7eHIjeX3OVO/o4vJ6O/IzG7EB2P/+G9Gf0dxrgcy+Dtyzm0CNpnZ9cC3gZsSPfdigr6jbwTKY7bLgKYEj0nk3PE2luvBOdf/32ZgE/6fbEEby+9zpn5HF5TJ35GZrQAeAdY551pGcu44G8v1ZPR31C/yF9NCMysd6bkX+gmDfECRBxwCqog+ZLhs0DF3MPDh5RuJnpth11MITI35/Bpwa5DXM9LfZ+C/MfBhbEZ+Rxe5noz9joAKoA64ZrS/HxlyPZn8HS0i+jB2FXAskhNj/o4CvfjIBd0OvI1/qvxfI233A/dHPhvwUGT/bhfTA0IAAAB8SURBVKD6YucG/WO014N/or4z8mNPulxPgtc0B3/XcQZojXwuyuDvKO71ZPh39AhwCtgR+VF7sXOD/jHa68nw7+gbkT7vAP4KXJes70glEEREQi7oMXoREUkxBb2ISMgp6EVEQk5BLyIScgp6EZGQU9CLiIScgl5EJOT+P5GWdo1t7nTmAAAAAElFTkSuQmCC\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.plot(h, CV_mean, lw = 3,  color = col_pink)\n",
    "#plt.plot(h,L_emp)\n",
    "#plt.plot(h,[sigma**2]*len(h))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "metadata": {},
   "outputs": [],
   "source": [
    "r.seed(201053)\n",
    "\n",
    "for j in range(J_cv):\n",
    "    y = gendata(f_x, sigma, n)\n",
    "    for i in range(  len(h) ):\n",
    "        f_h[:,i] = fh(x, x, y, h[i])\n",
    "        L_emp[j,i] = np.mean( (y - f_h[:,i])**2 , axis = 0)\n",
    "    # end h loop\n",
    "\n",
    "# end J loop  \n",
    "\n",
    "L_emp = L_emp/J_cv\n",
    "L_emp_mean = np.mean( L_emp, axis = 0 )"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Generalized Cross-validation"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "\\begin{equation}\n",
    "L^\\mathrm{GCV}(\\hat{f}_h) = GCV(h)=\n",
    "\\frac{1}{\\left[1-\\mathrm{tr}\\{S(h)\\}/n\n",
    "\\right]^2}\\,L_\\mathrm{emp}(\\hat{f}_h) \\rightarrow \\min_{h\\in \\mathbb{R}} \n",
    "\\end{equation}\n",
    "\n",
    "for the hat matrix $$\\mathrm{tr}\\{S(h)\\} = \\sum_{i=1}^n s_{ii}= n^{-1} K_h(0) \\sum_{i=1}^n\\{n^{-1} {\\sum_{k=1}^n K_h(x_i-x_k)}\\}^{-1}$$\n",
    "We know that \n",
    "$\\mathrm{tr}\\{S(h)\\} \\approx h^{-1} K(0)$\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "metadata": {},
   "outputs": [],
   "source": [
    "inner_t = (norm.pdf(0)/h)/n      # tr( S(h)/n )\n",
    "\n",
    "L_GCV = L_emp/( (1-inner_t)**2 )\n",
    "\n",
    "L_GCV_mean = np.mean( L_GCV, axis = 0 )"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[<matplotlib.lines.Line2D at 0x7fc04e187700>]"
      ]
     },
     "execution_count": 22,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.plot(h, sig2 + np.zeros( len(h) ), color=col_pink, lw=3)\n",
    "plt.plot(h, CV_mean,linestyle='dashdot', color=col_blue, lw=3)\n",
    "plt.plot(h, L_emp_mean, color= col_red, lw=3)\n",
    "plt.plot(h, L_GCV_mean, linestyle='dotted', color = 'green', lw=3)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.4"
  },
  "toc": {
   "base_numbering": 1,
   "nav_menu": {},
   "number_sections": true,
   "sideBar": true,
   "skip_h1_title": false,
   "title_cell": "Table of Contents",
   "title_sidebar": "Contents",
   "toc_cell": false,
   "toc_position": {},
   "toc_section_display": true,
   "toc_window_display": false
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}

```

automatically created on 2021-05-18