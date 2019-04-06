#####################################################################
# Author: Alex Marsh                                                #
# Title: Module for Bayesian Linear Regression                      #
# Date: 03/24/2019                                                  #
# Description: A module for functions to implement bayesian linear  # 
# regression (BLR).                                                 #
# Notes:                                                            #
# -As of right now, only supports conjugate Gaussian priors.        #
# -Assumes the variance of the errors is fixed i.e. non-random.     #
# -The optimal hyperparameter tau^2 (beta's prior variance) and     #
#  sigma^2 (the fixed variance of the errors) are estimated with    #
#  expectation maximizaton.                                         #
# -I hope to include other common priors used for BLR as well as    #
#  allow for the variance of the errors to also be a random         #
#  variable in the future.                                          #
#####################################################################

#===================================================================#
#                      Required Libraries                           #
#===================================================================#
import numpy as np

from numpy.linalg import inv
from numpy.linalg import eig
from numpy.linalg import norm

from math import sqrt

#===================================================================#
#                     Functions in Module                           #
#===================================================================#

def t(X):
    """
    Simple function to transpose matrix is the style of R's 
    transpose function.
    """
    return X.transpose()


def est_beta_bayes(y,X,lambd):
    """
    Estimate coefficient vector using Bayesian linear regression 
    assuming sigma^2 (the variance of the errors) is fixed. 
    Assumes Gaussian conjugate prior. Requires lambda parameter, 
    which is the the ratio of the (fixed) variance of the errors 
    (sigma^2) to the hyperparameter tau^2, which is beta's prior
    variance. So lambda=sigma^2/tau^2. Alternatively, if we 
    define alpha to be the prior precision of beta and gamma to be 
    the precision of the errors, then lambda = alpha/gamma.
    """

    k = len(X[0])
    XtX = matmul(t(X),X)
    beta_hat = matmul(inv(XtX+lambd*np.identity(k)),matmul(t(X),y))
    return beta_hat


def find_sig2_0(y,X):
    """
    Estimates sigma^2_hat (sample variance of linear regression 
    residuals) using vanilla OLS. This is used to get starting values 
    for the EM algorithm for finding the optimal hyperparameter tau^2 a
    nd optimal variance of the errors, sigma^2
    """
    k = len(X[0])
    XtX = matmul(t(X),X)
    beta_OLS = matmul(inv(XtX),matmul(t(X),y))
    es = y - matmul(X,beta_OLS)
    sig2_0 = np.var(es)
    return sig2_0

def find_tau2_sig2(y,X):
    """
    Finds the optimal hyperparameter tau^2 and optimal variance of 
    the errors sigma^2 for Bayesian linear regression. Estimation 
    is through the Expectation Maximization algorithm.
    
    Remember, alpha=1/tau^2 and gamma=1/sigma^2. The starting value 
    for both alpha and gamma is alpha_0=gamma_0=1/sqrt(sigma^2_OLS)
    
    The algorithm is described in detail in Pattern Recognition &
    Machine Learning by Christopher Bishop (2006), section 3.5.2
    "Maximizing the evidence function" starting on page 168.
    """
    k = len(X[0])
    n = len(X)
    sig2_0 = find_sig2_0(y,X)
    
    gamma_opt = sqrt(1/sig2_0) #set starting 
    alpha_opt = gamma_opt
    XtX = matmul(t(X),X)
    mean_post = est_beta_bayes(y,X,alpha_opt/gamma_opt)
    e_vals0 = eig(XtX)[0]
    theta = np.array([1/alpha_opt,1/gamma_opt],dtype=float)
    theta_prev = np.array([0,0],dtype=float)
    theta_diff = norm(theta-theta_prev)
    tol = 0.000000000001
    N = 1.0
    while theta_diff>tol and N<1000:
        theta_prev = theta
        e_vals = gamma_opt*e_vals0
        gamma = np.sum(e_vals/(np.add(e_vals,[alpha_opt]*k)))
        alpha_opt = gamma/np.dot(mean_post,mean_post)
        gamma_opt = np.sum(np.square(y-matmul(X,t(mean_post))))
        gamma_opt = (n-gamma)/gamma_opt
        mean_post = est_beta_bayes(alpha_opt/gamma_opt,y,X)
        theta = np.array([1/alpha_opt,1/gamma_opt],dtype=float)
        theta_diff = norm(theta-theta_prev)
        N=N+1
    return theta
