rejection <- function (n, mu,sigma, alpha, beta, c){
 
  # n is the desired size of the sample.
  # mu is the mean of the normal distribution that is truncated;
  # sigma is the variance of the normal distribution that is truncated;
  # alpha is the shape parameter of the proposal gamma distribution;
  # alpha is the rate parameter of the proposal gamma distribution;

  sample = matrix(0, nrow = n, ncol = 1); # Set up space for the sample.
  i = 1;
  count = 0; # set up counter to record the number of iterations.
  while(i <= n){
    x = rgamma(1, alpha, beta); # randomly generate a data point from the proposal gamma distribution.
    gx = dgamma(x, alpha, beta); # get the density of the data point generated.
    u = runif(1, 0, c*gx); # randomly select a value u from interval (0, cg(x))
    # Check if u is within the "envelope" of the density value of data point x in the normal distribution to be sampled.
    # if yes, then x is accepted as a new sample;
    if(u <= dnorm(x, mu, sigma)){
      sample[i] = x;
      i = i + 1;
    }
    count = count + 1;  # Update counter.
  }
  return (list(sample, count))  # Return both the sample and the counter.
}
