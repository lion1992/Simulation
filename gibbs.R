gibbs <- function(m, y, burnin){
  
  # m is the desired size of the sample.
  # y is the data (vector) of the variable whose posterior distribution that the user wish to sample.
  # burnin is the number of samples obtained in the beginning that the user wants to exclude.
  
  n = length(y); # get the length of the vector.
  # Define the mean and variance of the normal prior for the mean u.
  u0 = mean(y);     
  tao0 = 1;
  # Define the degree of freedom and the scale parameter of the scaled inversed chisqaured prior for the variance sigma.
  v0 = 1;
  sigma0 = var(y);
  # Initiate a value of the mean of the posterior distribution using the prior normal distribution.
  MU=rnorm(1, u0, tao0);
  # Initiate a value of the variance of the posterior distribution using prior gamma distribution (parameters transformed from those of scaled inversed chi-sqaured).
  rho = rgamma(1, v0/2, 2/(v0*sigma0));
  SIGMA = 1/rho;
  
  for(i in 1:m){
    # update the mean u of the posterior normal distribution for the mean MU, conditional on variance SIGMA.
    u = (u0/tao0 + n * mean(y)/SIGMA[i]) / (1/tao0 + n/SIGMA[i]); 
    # update the variance tao of the posterior normal distribution for the mean MU, conditional on variance SIGMA.
    tao = 1 / (1/tao0 + n/SIGMA[i]); 
    # generate a new sample of MU from the posterior normal distribution and combine it with the samples from previous iterations.
    MU = cbind(MU, rnorm(1, u, sqrt(tao)));
    
    # update the shape parameter of the posterior gamma distribution for the variance SIGMA, conditional on mean MU.
    alpha = (v0+n)/2;
    # update the scale parameter of the posterior gamma distribution for the variance SIGMA, conditional on mean MU.
    beta = 2/(v0*sigma0 + sum((y-MU[i+1])^2));
    # generate a new sample of SIGMA from the posterior gamma distribution and combine it with the samples from previous iterations.
    rho = rgamma(1, alpha, 1/beta);
    SIGMA = cbind(SIGMA, 1/rho);
  }
  
  par(mfrow = c(2,1)) # set up the space for two plots.
  # Plot the trace plot of samples of posterior mean MU.
  plot(c(burnin:m), MU[burnin:m+1], type="l", ylab=expression(mu),xlab="Iteration");
  # Plot the trace plot of samples of posterior variance SIGMA.
  plot(c(burnin:m), SIGMA[burnin:m+1], type="l", ylab=expression(sigma^2), xlab="Iteration");
  return(list(MU, SIGMA));
}