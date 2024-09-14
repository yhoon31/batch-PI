library(MASS)

##sampling distributions for simulation

f_X_unif <- function(n,p)
{
  return(matrix(runif(n*p),nrow=n))
}

f_X_normal <- function(n,mu_x,Sigma_x)
{
  return(mvrnorm(n,mu_x,Sigma_x))
}

f_Y_normal <- function(X,beta1,beta2,beta3)
{
  n <- dim(X)[1]
  return(X%*%beta1 + (X%*%beta2)^2+rnorm(n, sd=abs(X%*%beta3)))
}

f_Y_normal_truncated <- function(X,beta1,beta2,beta3,x_min,x_max)
{
  n <- dim(X)[1]
  Y <- rtruncnorm(n,a=x_min,b=x_max,mean = X%*%beta1 + (X%*%beta2)^2,sd=abs(X%*%beta3))
  return(Y)
}

f_Y_exp <- function(X,beta0, sigma)
{
  n <- dim(X)[1]
  return(log(1+exp(X%*%beta0 + sigma*rnorm(n))))
}

#function for weighted quantile
weighted_quantile <- function(probs,values,alpha)
{
  ##inputs
  #probs : vector of nonnegative number that sum to 1 / probability masses of discrete distribution
  #values : vector of real numbers / support of discrete distribution
  #alpha : value in (0,1), target level of inference
  
  ##output
  #(1-alpha)-th quantile of the discrete distribution
  
  values_sorted <- sort(values)
  q <- min(values_sorted[sapply(values_sorted,function(x) sum(probs[values <= x]))>=1-alpha])
  return(q)
}

#batch_PI for one quantile (one-sided)
q_quantile <- function(n,m,zeta,alpha)
{
  ##inputs
  #n : calibration size
  #m : test size
  #zeta : value in {1,...,m} : target rank
  #alpha : value in (0,1) : target level of inference
  
  #output : rank z s.t. S_{(z)} bounds the target quantile 
  
  size_total <- choose(n+m,m)
  probs <- sapply(1:(n+1), function(k) choose(k+zeta-2,zeta-1)*choose(n+m-k-zeta+1,m-zeta)/size_total)
  q <- weighted_quantile(probs, 1:(n+1), alpha)
  return(q)
}

#size of level set for two ranks
two_quantile_level <- function(n,m,r_1,r_2,zeta_1,zeta_2)
{
  return(choose(r_1+zeta_1-2,zeta_1-1)*choose(r_2-r_1+zeta_2-zeta_1-1,zeta_2-zeta_1-1)*choose(n+m-r_2-zeta_2+1,m-zeta_2))
}

#probability of p <= R_1 <= R_2 <= q
two_quantile_lower_upper_prob <- function(n,m,p,q,zeta_1,zeta_2)
{
  total <- choose(n+m,m)
  prob <- 0
  for (r_1 in max(1,p+1):min(n+1, q))
  {
    for (r_2 in r_1:min(n+1, q))
    {
      prob <- prob+two_quantile_level(n,m,r_1,r_2,zeta_1,zeta_2)/total
    }
  }
  return(prob)
}

#batch PI bounds for inference on two quartiles
two_quantiles_lower_upper_bound <- function(n,m,zeta_1,zeta_2,alpha)
{
  t <- 1
  for(l in 1:n)
  {
    p <- two_quantile_lower_upper_prob(n,m,max(round(0.25*n)-t+1,1)-1,min(round(0.75*n)+t,n+1),zeta_1,zeta_2)
    p#rint(p)
    if(p > 1-alpha) break
    t <- t+1
  }
  return(c(max(round(0.25*n)-t+1,1)-1+1,min(round(0.75*n)+t,n+1)+1))
}


#dynamic programming procedure (for inference on the mean)
maximize_sum <- function(n, m, q, s) {
  
  dp <- matrix(-Inf, nrow = m + 1, ncol = q + 1)
  dp[1, 1] <- 0
  
  for (i in 1:m) {
    for (j in 1:q) {
      for (r in 1:n) {
        if (j >= r) {
          dp[i + 1, j + 1] <- max(dp[i + 1, j + 1], dp[i, j + 1 - r] + s[r])
        }
      }
    }
  }
  
  return(max(dp[m + 1, ]))
}


