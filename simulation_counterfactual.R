setwd("C:/Wharton/batch predictive inference/code")
source("functions.R")


#Inference on counterfactual median
p <- 20      #dimension
n <- 200    #calibration size
m <- 40     #test size
zeta <- 20  #target quantile

alpha_set <- seq(0.05,0.2,0.025)

set.seed(1234)
mu_x <- runif(p,-2,2)
Sigma_x <- 5*diag(p)
f_X <- f_X_unif

#treatment assignment
beta_A <- runif(p,-0.2,0)

#parameter for counterfactual distribution
beta_Y <- runif(p,0,1)/p

tn <- 500 #trial number
n_super <- 1000

coverage_median <- matrix(0,nrow=tn,length(alpha_set))
pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
  #generate calibration & test data of size n & m
  X_super <- f_X_unif(n_super,p)
  p_A_super <- exp(X_super%*%beta_A)/(1+exp(X_super%*%beta_A))
  A_super <- rbinom(n_super,1,p_A_super)
  Y_0_super <- rbeta(n_super,1+X_super%*%beta_Y,1-X_super%*%beta_Y)
  Y_1_super <- rbeta(n_super,1-X_super%*%beta_Y,1+X_super%*%beta_Y)
  Y_super <- (1-A_super)*Y_0_super + A_super*Y_1_super
  ind_cal <- which(A_super==0)[1:n]
  ind_test <- which(A_super==1)[1:m]
  
  X_cal <- X_super[ind_cal,]
  Y_cal <- Y_super[ind_cal]
  Y_test_counterfactual <- Y_0_super[ind_test]
  
  #rejection sampling
  p_A_cal <- exp(X_cal%*%beta_A)/(1+exp(X_cal%*%beta_A))
  B <- rbinom(n,1,p_A_cal/(1-p_A_cal))
  ind_rej <- which(B==1)
  Y_rej <- Y_cal[ind_rej]
  
  n_rej <- length(ind_rej)
  score_cal <- c(0,sort(Y_rej),1)
  score_test <- sort(Y_test_counterfactual)
  
  #batch PI
  q_U_set <- sapply(alpha_set, function(alpha) q_quantile(n_rej,m,zeta,alpha/2))+1
  q_L_set <- sapply(alpha_set, function(alpha) q_quantile(n_rej,m,zeta,1-alpha/2))
  
  target <- score_test[zeta]
  bound_L <- score_cal[q_L_set]
  bound_U <- score_cal[q_U_set]
  
  coverage_median[i,] <- (bound_L <= target)*(target <= bound_U)
    
  setTxtProgressBar(pb, i)
  
}

mean_median <- colMeans(coverage_median)
se_median <- apply(coverage_median,2,sd)/sqrt(tn)

#Inference on counterfactual quartiles
alpha_set <- seq(0.05,0.2,0.025)
p <- 20              
n <- 200
m <- 40
zeta_1 <- 10
zeta_2 <- 30


tn <- 500 #trial number
coverage_quartiles <- matrix(0,nrow=tn,ncol=length(alpha_set))
n_super <- 1000

pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
  #generate calibration & test data of size n & m
  X_super <- f_X_unif(n_super,p)
  p_A_super <- exp(X_super%*%beta_A)/(1+exp(X_super%*%beta_A))
  A_super <- rbinom(n_super,1,p_A_super)
  Y_0_super <- rbeta(n_super,1+X_super%*%beta_Y,1-X_super%*%beta_Y)
  Y_1_super <- rbeta(n_super,1-X_super%*%beta_Y,1+X_super%*%beta_Y)
  Y_super <- (1-A_super)*Y_0_super + A_super*Y_1_super
  ind_cal <- which(A_super==0)[1:n]
  ind_test <- which(A_super==1)[1:m]
  
  X_cal <- X_super[ind_cal,]
  Y_cal <- Y_super[ind_cal]
  Y_test_counterfactual <- Y_0_super[ind_test]
  
  #rejection sampling
  p_A_cal <- exp(X_cal%*%beta_A)/(1+exp(X_cal%*%beta_A))
  B <- rbinom(n,1,p_A_cal/(1-p_A_cal))
  ind_rej <- which(B==1)
  Y_rej <- Y_cal[ind_rej]

  n_rej <- length(ind_rej)
  score_cal <- c(0,sort(Y_rej),1)
  score_test <- sort(Y_test_counterfactual)
  
  #batch PI
  index_quartiles <- sapply(alpha_set, function(alpha) two_quantiles_lower_upper_bound(n_rej,m,zeta_1,zeta_2,alpha))

  bound_set_quartiles <- apply(index_quartiles,2,function(v) score_cal[v])

  coverage_quartiles[i,] <- apply(bound_set_quartiles,2,function(v) (score_test[zeta_1] >= v[1])*(score_test[zeta_2] <= v[2]))
 
  setTxtProgressBar(pb, i)
  
}


mean_quartiles <- colMeans(coverage_quartiles)
se_quartiles <- apply(coverage_quartiles,2,sd)/sqrt(tn) 


par(mfrow=c(1,2),cex.lab = 1.6, cex.main=1.6, cex.axis=1.6,mar=c(2,3,2,3), xpd=NA, oma = c(3, 3, 1, 1))
plot(1-alpha_set,mean_median,pch=16,ylim=c(0.76,1),lwd=2,xlab=expression(1-alpha),ylab="Coverage rate",main="Median")
arrows(1-alpha_set, mean_median-se_median, 1-alpha_set, mean_median+se_median, length=0.02, angle=90, code=3,lwd=2)
segments(0.8,0.8,0.95,0.95,lty=2,lwd=2)

plot(1-alpha_set,mean_quartiles,pch=16,ylim=c(0.76,1),lwd=2,xlab=expression(1-alpha),ylab="",main="Quartiles")
arrows(1-alpha_set, mean_quartiles-se_quartiles, 1-alpha_set, mean_quartiles+se_quartiles, length=0.02, angle=90, code=3,lwd=2)
segments(0.8,0.8,0.95,0.95,lty=2,lwd=2)


#Inference on counterfactual mean
alpha_set <- seq(0.05,0.2,0.025)
p <- 20              
n <- 100
m <- 5


tn <- 500
coverage_mean <- matrix(0,nrow=tn,ncol=length(alpha_set))

set.seed(1234)
n_super <- 1000
n_sample <- 20000
pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
  #generate calibration & test data of size n & m
  X_super <- f_X_unif(n_super,p)
  p_A_super <- exp(X_super%*%beta_A)/(1+exp(X_super%*%beta_A))
  A_super <- rbinom(n_super,1,p_A_super)
  Y_0_super <- rbeta(n_super,1+X_super%*%beta_Y,1-X_super%*%beta_Y)
  Y_1_super <- rbeta(n_super,1-X_super%*%beta_Y,1+X_super%*%beta_Y)
  Y_super <- (1-A_super)*Y_0_super + A_super*Y_1_super
  ind_cal <- which(A_super==0)[1:n]
  ind_test <- which(A_super==1)[1:m]
  
  X_cal <- X_super[ind_cal,]
  Y_cal <- Y_super[ind_cal]
  Y_test_counterfactual <- Y_0_super[ind_test]
  
  #rejection sampling
  p_A_cal <- exp(X_cal%*%beta_A)/(1+exp(X_cal%*%beta_A))
  B <- rbinom(n,1,p_A_cal/(1-p_A_cal))
  ind_rej <- which(B==1)
  Y_rej <- Y_cal[ind_rej]
  
  n_rej <- length(ind_rej)
  score_cal <- c(sort(Y_rej),1)
  score_test <- sort(Y_test_counterfactual)
  
  #batch PI
  ind_sample <- t(sapply(1:n_sample,function(i) sort(sample(1:(n_rej+1),m,replace=T))))
  q_set <- sapply(alpha_set,function(alpha) quantile(rowSums(ind_sample),1-alpha))
 
  bound_set <- sapply(q_set,function(q) maximize_sum(n_rej+1, m, q, score_cal)/m)

  coverage_mean[i,] <- sapply(bound_set,function(b) mean(score_test) <= b)
  


  setTxtProgressBar(pb, i)
  
}

mean_mean <- colMeans(coverage_mean)
se_mean <- apply(coverage_mean,2,sd)/sqrt(tn) 


par(mfrow=c(1,1),cex.lab = 1.6, cex.main=1.6, cex.axis=1.6,mar=c(2,3,2,3), xpd=NA, oma = c(3, 3, 1, 1))
plot(1-alpha_set,mean_mean,pch=16,ylim=c(0.76,1),lwd=2,xlab=expression(1-alpha),ylab="Coverage rate",main="Mean (n=100, m=10)")
arrows(1-alpha_set, mean_mean-se_mean, 1-alpha_set, mean_mean+se_mean, length=0.02, angle=90, code=3,lwd=2)
segments(0.8,0.8,0.95,0.95,lty=2,lwd=2)


##------------------------------------------------------------------------------------------------------------

#batch PI for the mean, under different score distributions

          
n <- 40
m <- 10
alpha_set <- seq(0.05,0.2,0.025)

tn <- 500
coverage_1 <- matrix(0,nrow=tn,ncol=length(alpha_set))
coverage_2 <- matrix(0,nrow=tn,ncol=length(alpha_set))
coverage_3 <- matrix(0,nrow=tn,ncol=length(alpha_set))
coverage_4 <- matrix(0,nrow=tn,ncol=length(alpha_set))
coverage_5 <- matrix(0,nrow=tn,ncol=length(alpha_set))

set.seed(1234)
n_sample <- 20000

ind_sample <- t(sapply(1:n_sample,function(i) sort(sample(1:(n+1),m,replace=T))))
q_set <- sapply(alpha_set,function(alpha) quantile(rowSums(ind_sample),1-alpha))

pb <- txtProgressBar(min = 0, max = tn, style = 3)
for(i in 1:tn)
{
 
  score_cal_1 <- c(sort(runif(n)),1)
  score_test_1 <- sort(runif(m))
  bound_set_1 <- sapply(q_set,function(q) maximize_sum(n+1, m, q, score_cal_1)/m)
  coverage_1[i,] <- sapply(bound_set,function(b) mean(score_test_1) <= b)
  score_cal_2 <- c(sort(rbeta(n,0.8,0.8)),1)
  score_test_2 <- sort(rbeta(m,0.8,0.8))
  bound_set_2 <- sapply(q_set,function(q) maximize_sum(n+1, m, q, score_cal_2)/m)
  coverage_2[i,] <- sapply(bound_set,function(b) mean(score_test_2) <= b)
  score_cal_3 <- c(sort(rbeta(n,0.5,0.5)),1)
  score_test_3 <- sort(rbeta(m,0.5,0.5))
  bound_set_3 <- sapply(q_set,function(q) maximize_sum(n+1, m, q, score_cal_3)/m)
  coverage_3[i,] <- sapply(bound_set,function(b) mean(score_test_3) <= b)
  score_cal_4 <- c(sort(rbeta(n,1.2,1.2)),1)
  score_test_4 <- sort(rbeta(m,1.2,1.2))
  bound_set_4 <- sapply(q_set,function(q) maximize_sum(n+1, m, q, score_cal_4)/m)
  coverage_4[i,] <- sapply(bound_set,function(b) mean(score_test_4) <= b)
  score_cal_5 <- c(sort(rbeta(n,1.5,1.5)),1)
  score_test_5 <- sort(rbeta(m,1.5,1.5))
  bound_set_5 <- sapply(q_set,function(q) maximize_sum(n+1, m, q, score_cal_5)/m)
  coverage_5[i,] <- sapply(bound_set,function(b) mean(score_test_5) <= b)
  
  setTxtProgressBar(pb, i)
  
}


cov_1 <- colMeans(coverage_1)
cov_2 <- colMeans(coverage_2)
cov_3 <- colMeans(coverage_3)
cov_4 <- colMeans(coverage_4)
cov_5 <- colMeans(coverage_5)
se_1 <- apply(coverage_1,2,sd)/sqrt(tn)
se_2 <- apply(coverage_2,2,sd)/sqrt(tn)
se_3 <- apply(coverage_3,2,sd)/sqrt(tn)
se_4 <- apply(coverage_4,2,sd)/sqrt(tn)
se_5 <- apply(coverage_5,2,sd)/sqrt(tn)

col1 <- rgb(70/255, 130/255, 180/255, alpha = 1) 
col2 <- rgb(0/255, 0/255, 139/255, alpha = 1)   
col3 <- rgb(0/255, 128/255, 128/255, alpha = 1)   
col4 <- rgb(178/255, 34/255, 34/255, alpha = 1)   
col5 <- rgb(240/255, 128/255, 128/255, alpha = 1) 
col1_1 <- rgb(70/255, 130/255, 180/255, alpha = 0.6) 
col2_1 <- rgb(0/255, 0/255, 139/255, alpha = 0.6)   
col3_1 <- rgb(0/255, 128/255, 128/255, alpha = 0.6)   
col4_1 <- rgb(178/255, 34/255, 34/255, alpha = 0.6)   
col5_1 <- rgb(240/255, 128/255, 128/255, alpha = 0.6) 

ttt <- seq(0.01,0.99,0.01)


par(mfrow=c(1,2),cex.lab = 1.6, cex.main=1.6, cex.axis=1.6,mar=c(2,3,2,3), xpd=NA, oma = c(3, 3, 1, 9))
matplot(ttt,cbind(dunif(ttt),dbeta(ttt,0.8,0.8),dbeta(ttt,0.5,0.5),dbeta(ttt,1.2,1.2),dbeta(ttt,1.5,1.5))
,type="l",ylim=c(0,3.5),ylab="Density",xlab="",lwd=2,lty=1,col=c(col1,col2,col3,col4,col5))


matplot(1-alpha_set,cbind(cov_1,cov_2,cov_3,cov_4,cov_5),type="l",ylim=c(0.79,1),ylab="Coverage rate",xlab=expression(1-alpha),lwd=2,lty=1,col=c(col1,col2,col3,col4,col5))
arrows(1-alpha_set, cov_1-se_1, 1-alpha_set, cov_1+se_1, length=0.02, angle=90, code=3,lwd=1,col=col1_1)
arrows(1-alpha_set, cov_2-se_2, 1-alpha_set, cov_2+se_2, length=0.02, angle=90, code=3,lwd=1,col=col2_1)
arrows(1-alpha_set, cov_3-se_3, 1-alpha_set, cov_3+se_3, length=0.02, angle=90, code=3,lwd=1,col=col3_1)
arrows(1-alpha_set, cov_4-se_4, 1-alpha_set, cov_4+se_4, length=0.02, angle=90, code=3,lwd=1,col=col4_1)
arrows(1-alpha_set, cov_5-se_5, 1-alpha_set, cov_5+se_5, length=0.02, angle=90, code=3,lwd=1,col=col5_1)
segments(0.8,0.8,0.95,0.95,lty=2,lwd=2)

legend1 = "Unif([0,1])"
legend2 = "Beta(0.8,0.8)"
legend3 = "Beta(0.5,0.5)"
legend4 = "Beta(1.2,1.2)"
legend5 = "Beta(1.5,1.5)"
legend(0.97,0.95, legend = c(legend1, legend2, legend3, legend4, legend5), col = c(col1,col2,col3,col4,col5),lty=1,lwd=2,box.lty=0,cex=1.2,bty="n",)
