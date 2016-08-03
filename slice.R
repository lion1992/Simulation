slice <- function(n, X, y, w, m){
  # m is the desired size of the sample.
  # X is the data set of entries of the predicting variables
  # y is the target.
  # w is an estimate of the typical size of a slice;
  # m is an integer limiting the size of a slice to mw; 
  d = length(X[1,]);
  beta0 = rep(0, d);
  BETA = beta0;
  i=1;
  while(i<=n){
    j=1;
    # Sample the variables (betas) in turn
    while (j <= d){
      lf = (t(y) %*% X %*% beta0-sum(log(1+exp(X %*% beta0)))); # compute the loglikelihood of the current status. 
      ys = runif(1, 0, exp(lf)); # uniformly sample a new point from the interval from 0 to the likelihood of current status.
      z = log(ys); # compute the auxiliary variable and define the slice as  S = {x : z < g(x)}.
      
      # Stepping out
      U = runif(1, 0, 1);
      L = beta0[j] - w*U;
      R = L+ w;
      V = runif(1, 0 ,1);
      J = floor(m*V);
      K = (m-1) - J;
      betaLR = beta0;
      betaLR[j] = L;
      
      ## stop the stepping out loop for the lower bound (left bound) until J hits 0 
      ### or z gets bigger than the loglikelihood of beta at the lower bound;
      while(J>0 && z < (t(y) %*% X %*% betaLR-sum(log(1+exp(X %*% betaLR))))){
        L = L-w;
        J = J-1;
        betaLR[j] = L;
      }

      betaLR[j] = R;
      ## stop the stepping out loop for the upper bound (right bound) until K hits 0 
      ### or z gets bigger than the loglikelihood of beta at the upper bound;
      while(K>0 && z < (t(y) %*% X %*% betaLR-sum(log(1+exp(X %*% betaLR))))){
        R = R+w;
        K = K-1;
        betaLR[j] = R;
      }
      
      # Shrinkage
      Ls = L;
      Rs = R;
      done= FALSE;
      while (done == FALSE){
        # uniformly sample a new beta[j], x1, from the lower to the upper bound;
        U = runif(1, 0, 1);
        x1 = Ls + U * (Rs-Ls); 
        betaLR[j] = x1;
        # check if z is smaller than the loglikelihood of the newly sampled beta vector (which contains x1);
        if(z < (t(y) %*% X %*% betaLR-sum(log(1+exp(X %*% betaLR))))){
          # if yes, accept the newly sampled beta as the new sample and break the loop.
          beta0[j] = x1;
          done = TRUE;
          break;
        }
        # if not, check if x1 is smaller than the current beta[j];
        else if(x1 < beta0[j]){
          # if yes, make x1 as the new lower bound;
          Ls = x1;
        }
        else {
          # if not, make x1 as the new upper bound;
          Rs = x1;
        }
      }
      
      j=j+1;
    }
    BETA = rbind(BETA, beta0); # combine the new sample with those from previous iterations.
    i=i+1;
  }
  # Plot the traceplots of samples for each beta.
  par(mfrow=c(d,1));
  for (k in 1:d){
    plot(c(0:n), BETA[,k], type="l", ylab=substitute(paste(beta,p),list(p=k)), xlab="Interation");
  }
  return (BETA);
}