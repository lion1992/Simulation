MH <- function (m, X, y, sigma){
  
  # m is the desired size of the sample.
  # X is the data set of entries of the predicting variables
  # y is the target.
  # sigma is the vector of standard deviations of the normal distribution for the random walk process (one for each predicting variable).
  
 n = length(y) ;
 d = length(X[1,]);
 beta0 = rep(0, d); # Initiate a beta vector.
 BETA = beta0; 
 i=1;
 while(i<=m){
   j=1;
   betaprop=beta0;
   # Sample the variables (betas) in turn;
   while (j <= d){
     # Propose a new beta vector by randomly generating a new value for the jth variable in the beta vector from the normal distribution with current beta[j] as mean. 
     betaprop[j] = rnorm(1, beta0[j], sigma[j]); 
     # Compute the loglikelihood of the current beta vector.
     lf0 = (t(y) %*% X %*% beta0-sum(log(1+exp(X %*% beta0))));
     # Compute the loglikelihood of the newly proposed beta vector. 
     lfprop = (t(y) %*% X %*% betaprop-sum(log(1+exp(X %*% betaprop))));
     # Compute the acceptance probability.
     a = min(1, exp((lfprop + log(dnorm(beta0[j], betaprop[j], sigma[j])))- (lf0 + log(dnorm(betaprop[j], beta0[j], sigma[j])))));
     # randomly generate u from Uniform(0, 1);
     u = runif(1, 0, 1);
     # check if the acceptance probability is bigger than u
     if(u < a){
       beta0[j] = betaprop[j]; # if yes, accept the proposal as the new sample;
     }
     else{
       betaprop[j] = beta0[j]; # if no, take the current state as the new sample;
     }
     j = j+1;
   }
   BETA = rbind(BETA, beta0); # combine the new sample with those from previous iterations.
   i=i+1;
 }
 # Plot the traceplots of samples for each beta.
 par(mfrow=c(d,1));
 for (k in 1:d){
   s=sigma[k];
   plot(c(0:m), BETA[,k], type="l", main=substitute(paste(sigma, " = ", s), list(s=sigma[k])),ylab=substitute(paste(beta,p),list(p=k) ), xlab="Interation");
 }
 return (BETA);
 }