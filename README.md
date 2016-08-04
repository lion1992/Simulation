# Simulation
##This repository contains self written R code for the 6 simulation methods for different purposes.
##1. HMC is a function for sampling the coefficients of the predicting variables, which is the β vector, with logistic likelihood and Normal prior, using Hamiltonian Monte Carlo algorithm with leapfrogs method.
##2. MH is a function for sampling the coefficients of the predicting variables, which is the β vector, with logistic likelihood and improper prior 1, using Metropolis-Hastings algorithm.
##3. slice is a function for sampling the coefficients of the predicting variables, which is the β vector, with logistic likelihood and improper prior 1, using Slice sampler. 
##4. rejection is a function for sampling a normal distribution truncated from left at 0 using a Gamma distribution as proposal distribution. 
##5. importance is a function for finding a user-specified order of moment of a normal distribution truncated from left at 0 using a Gamma distribution as proposal distribution. 
##6. gibbs is a function for sampling from a posterior distribution, which is a normal distribution, with priors for mean and variance as well as the likelihood both being normal too.

