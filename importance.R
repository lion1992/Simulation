importance <- function(n, mu, sigma, alpha, beta, moment){
  # n is the desired size of the sample.
  # mu is the mean of the normal distribution that is truncated;
  # sigma is the variance of the normal distribution that is truncated;
  # alpha is the shape parameter of the proposal gamma distribution;
  # alpha is the rate parameter of the proposal gamma distribution;
  # moment is the moment that the user want to get.
 
  w = matrix(0, nrow = n, ncol = 1);  # Set up space for the weights.
  hx = matrix(0, nrow = n, ncol = 1); # Set up space for x^(moment).
  for(i in 1:n){
    x = rgamma(1, alpha, beta);  # Generate data point x randomly from the proposal gamma distribution.
    gx = dgamma(x, alpha, beta); # get the gamma density of the data point.
    fx = dnorm(x, mu, sigma); # get the targeted normal density of the data point.
    w[i] = fx/gx; # compute weight.
    hx[i] =  x^moment; # compute x^(moment)
  }
  EXsquared = (t(w) %*% hx)/sum(w); # compute the weighted sum of the x^(moment) obtained in the loop.
  return (EXsquared);
}