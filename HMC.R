HMC <- function (m, X, y, mu, sigma, L,e){
  
  # m is the desired size (number of time points) of the sample.
  # X is the data set of entries of the predicting variables
  # y is the target.
  # mu is the vector of means of the prior Normal distribution of the variables.
  # sigma is the vector of standard deviations of the prior Normal distribution of the variables.
  # L is the number of leapfrog steps.
  # e is the step size.
  
  n = length(y) ;
  d = length(X[1,]);
  beta0 = rep(0, d); # Initiate a beta vector.
  BETA = beta0; 
  t=1;

  while(t<=m){ 
    # Initialize p
    p0 = rep(0, d);
    for(np in 1:d)
    {
      p0[np] = rnorm(1, 0, 1);
    }
   
    betaprop=beta0;
    pprop = p0;
    # Leapfrogs method
    for(l in 1:L)
    {
      # propose the predicting variables (betas) and the momentum variables (p's) in turn in one leapfrog step;
      for(j in 1:d)
      {
        # compute the gradient of potential energy (negative loglikelihood of the posterior distribution) wrt betaj;
        gradL = -t(X[,j]) %*% (y - 1/(1+exp(-X %*% betaprop))) + (betaprop[j] - mu[j])/(sigma[j]^2);
        # update the jth element in p for half a leapfrog step;
        pprop[j] = pprop[j] - e/2 * gradL;
        # propose a new BETA vector by proposing a new betaj; 
        betaprop[j] = betaprop[j] + e*pprop[j];
        # compute the gradient of the potential energy wrt the newly proposed betaj;
        gradL = -t(X[,j]) %*% (y - 1/(1+exp(-X %*% betaprop))) + (betaprop[j] - mu[j])/(sigma[j]^2);
        # update the jth element in p for the full leapfrog step;
        pprop[j] = pprop[j] - e/2 * gradL;
      }
    }
    
    # compute the negative of the Hamiltonian functions with original beta and with newly proposed beta respectively;
    nh0 = (t(y) %*% X %*% beta0-sum(log(1+exp(X %*% beta0)))-1/2*t(p0) %*% p0);
    nhprop = (t(y) %*% X %*% betaprop-sum(log(1+exp(X %*% betaprop)))-1/2*t(pprop) %*% pprop);
    for (j in 1:d)
    {
      nh0 = nh0 + log(dnorm(beta0[j], mu[j], sigma[j]));
      nhprop = nhprop + log(dnorm(betaprop[j], mu[j], sigma[j]));
    }
    
    # Compute the acceptance probability.
    a = min(1, exp(nhprop - nh0));
    # randomly generate u from Uniform(0, 1);
    u = runif(1, 0, 1);
    # check if the acceptance probability is bigger than u
    if(u < a){
      beta0 = betaprop; # if yes, accept the proposal as the new sample;
    }
    else{
      betaprop = beta0; # if no, take the current state as the new sample;
    }
    BETA = rbind(BETA, beta0); # combine the new sample with those from previous iterations.
    t=t+1;
  }
  # Plot the traceplots of samples for each beta.
  par(mfrow=c(d,1));
  for (k in 1:d){
    s=sigma[k];
    plot(c(0:m), BETA[,k], type="l", main=substitute(paste(sigma, " = ", s), list(s=sigma[k])),ylab=substitute(paste(beta,p),list(p=k) ), xlab="Interation");
  }
  return (BETA);
}