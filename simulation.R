source("functions.R")

library(MASS)
library(randomForest)

###Simulataneous predictive inference with PAC guarantee

p <- 20               #dimension
n_train <- 200        #training size
n <- 200              #calibration size
m <- 100             #test size

alpha_conformal <- 0.1 #level of split conformal
delta_batch <- 0.1     #levels of batch PI
alpha_batch_0 <- 0.1
alpha_batch_1 <- 0.05
alpha_batch_2 <- 0.01

set.seed(1234)
mu_x <- runif(p,-2,2)
Sigma_x <- 5*diag(p)

beta1 <- runif(p,-1,1)
beta2 <- runif(p,-0.2,0.2)
beta3 <- runif(p,-0.6,0.6)

f_X <- f_X_normal
f_Y <- f_Y_normal

#generate training data & fit regression
X_train <- f_X(n_train,mu_x,Sigma_x)
Y_train <- f_Y(X_train,beta1,beta2,beta3)

rf <- randomForest(X_train,Y_train)

#repeat simulations with multiple generations of calibration/test data

K=500 #trial number
coverage_conformal <- rep(1,K)
width_conformal <- rep(1,K)
coverage_batch_0 <- rep(1,K)
width_batch_0 <- rep(1,K)
coverage_batch_1 <- rep(1,K)
width_batch_1 <- rep(1,K)
coverage_batch_2 <- rep(1,K)
width_batch_2 <- rep(1,K)


for(k in 1:K)
{
  X <- f_X(n,mu_x,Sigma_x)
  Y <- f_Y(X,beta1,beta2,beta3)
  X_test <- f_X(m,mu_x,Sigma_x)
  Y_test <- f_Y(X_test,beta1,beta2,beta3)
  muhat_cal <- predict(rf,X)
  muhat_test <- predict(rf,X_test)
  names(muhat_cal) <- NULL
  names(muhat_test) <- NULL
  
  score_cal <- c(sort(abs(Y-muhat_cal)),100000)
  score_test <- abs(Y_test - muhat_test)
  
  rank_conformal <- ceiling((n+1)*(1-alpha_conformal))
  rank_batch_0 <- q_quantile(n,m,ceiling((1-delta_batch)*m),alpha_batch_0)
  rank_batch_1 <- q_quantile(n,m,ceiling((1-delta_batch)*m),alpha_batch_1)
  rank_batch_2 <- q_quantile(n,m,ceiling((1-delta_batch)*m),alpha_batch_2)
  
  bound_conformal <- score_cal[rank_conformal]
  bound_batch_0 <- score_cal[rank_batch_0]
  bound_batch_1 <- score_cal[rank_batch_1]
  bound_batch_2 <- score_cal[rank_batch_2]
  
  coverage_conformal[k] <- mean(score_test <= bound_conformal)
  coverage_batch_0[k] <- mean(score_test <= bound_batch_0)
  coverage_batch_1[k] <- mean(score_test <= bound_batch_1)
  coverage_batch_2[k] <- mean(score_test <= bound_batch_2)
  width_conformal[k] <- 2*bound_conformal
  width_batch_0[k] <- 2*bound_batch_0
  width_batch_1[k] <- 2*bound_batch_1
  width_batch_2[k] <- 2*bound_batch_2
}

#histograms
col1 <- rgb(1, 192/255, 203/255, alpha = 0.4)
colb_0 <- rgb(255/255, 187/255, 120/255, alpha=0.4) #######
colb_1 <- rgb(135/255, 206/255, 250/255, alpha=0.4)
colb_2 <- rgb(200/255, 150/255, 250/255, alpha=0.3)

legend1 <- "split conformal"
legendb_0 <- expression(paste("batch-PI (",alpha,"=0.1)"))
legendb_1 <- expression(paste("batch-PI (",alpha,"=0.05)"))
legendb_2 <- expression(paste("batch-PI (",alpha,"=0.01)"))

cov_min <- min(c(coverage_conformal,coverage_batch_0, coverage_batch_1,coverage_batch_2))-0.01
cov_max <- 1
width_min <- min(c(width_conformal,width_batch_0,width_batch_1,width_batch_2))-1
width_max <- max(c(width_conformal,width_batch_0,width_batch_1,width_batch_2))+1

breaks_coverage <- seq(cov_min,cov_max,0.01)
breaks_width <- seq(width_min,width_max,1)

par(mfrow=c(1,2), cex.lab = 1.5, cex.axis=1.5,mar=c(5,5,2,0),oma = c(0, 0, 0, 14),xpd=NA)
hist(coverage_conformal, breaks = breaks_coverage, col = col1, main = "", xlab = "Test coverage rate", ylab = "", xlim=c(cov_min,cov_max),freq=F,ylim=c(0,25))
hist(coverage_batch_0, breaks = breaks_coverage, col = colb_0, freq=F, add = TRUE)
hist(coverage_batch_1, breaks = breaks_coverage, col = colb_1, freq=F, add = TRUE)
hist(coverage_batch_2, breaks = breaks_coverage, col = colb_2, freq=F, add = TRUE)

lines(x = c(0.9, 0.9), y = c(0, 25), col = "black", lty = 2, lwd = 2) 

hist(width_conformal, breaks = breaks_width, col = col1, main = "", xlab = "Prediction interval width", ylab = "", xlim=c(width_min,width_max),freq=F,ylim=c(0,0.3))
hist(width_batch_0, breaks = breaks_width, col = colb_0, freq=F, add = TRUE)
hist(width_batch_1, breaks = breaks_width, col = colb_1, freq=F, add = TRUE)
hist(width_batch_2, breaks = breaks_width, col = colb_2, freq=F, add = TRUE)

legend(width_max+0.2,0.13, legend = c(legend1, legendb_0,legendb_1, legendb_2), fill = c(col1,colb_0,colb_1,colb_2),box.lty=0,cex=1.5,bg="transparent", ncol = 1)

coverage <- rbind(coverage_conformal,coverage_batch_0,coverage_batch_1,coverage_batch_2)
width <- rbind(width_conformal,width_batch_0,width_batch_1,width_batch_2)

#table

table_mean <- cbind(rowMeans(coverage), apply(coverage,1,function(v) mean(v >= 0.9)), rowMeans(width))
table_sd <- cbind(apply(coverage,1,function(v) sd(v)), apply(coverage,1,function(v) sd(v >= 0.9)), apply(width,1,function(v) sd(v)))/sqrt(K)

table_mean
table_sd


###Selection with false discovery control
p <- 20               #dimension
n_train <- 500        #training size
n <- 1000              #calibration size
m <- 100               #test size
alpha <- 0.1          #target level

set.seed(1234)
mu_x <- runif(p,-2,2)
Sigma_x <- 5*diag(p)

beta <- runif(p,-1,1)
sigma <- 5

X_train <- f_X_normal(n_train,mu_x,Sigma_x)
Y_train <- f_Y_exp(X_train,beta,sigma)

rf <- randomForest(X_train,Y_train)

eta_set <- seq(0,10,2)  #target bound for false selections


#repeat simulations with multiple generations of calibration/test data
K <- 500 #trial number
c <- 5 #threshold

false_pos <- matrix(0,nrow=length(eta_set), ncol=K)
power <- matrix(0,nrow=length(eta_set), ncol=K)

pb <- txtProgressBar(min = 0, max = length(eta_set), style = 3)

t <- 1
for(eta in eta_set)
{
  num_rej <- rep(1,K)
  false_positive_trial <- rep(1,K)
  power_trial <- rep(1,K)
  for(k in 1:K)
  {
    X <- f_X_normal(n,mu_x,Sigma_x)
    Y <- f_Y_exp(X,beta,sigma)
    X_test <- f_X_normal(m,mu_x,Sigma_x)
    Y_test <- f_Y_exp(X_test,beta,sigma)
    muhat_cal <- predict(rf,X)
    muhat_test <- predict(rf,X_test)
    names(muhat_cal) <- NULL
    names(muhat_test) <- NULL
    
    score_cal <- c(sort(muhat_cal*(Y <= c)),100000)
    
    rank_batch <- q_quantile(n,m,m-eta,alpha)
    rejection_threshold <- score_cal[rank_batch]
    
    num_rej[k] <- sum(muhat_test > rejection_threshold)
    false_positive_trial[k] <- sum((muhat_test > rejection_threshold)*(Y_test <= c))
    power_trial[k] <- sum((muhat_test > rejection_threshold)*(Y_test > c))
  }
  false_pos[t,] <- false_positive_trial
  power[t,] <- power_trial
  setTxtProgressBar(pb, t)
  t <- t+1
}

false_pos_prob <- sapply(1:length(eta_set), function(t) mean(false_pos[t,] > eta_set[t]))
false_pos_prob_sd <- sapply(1:length(eta_set), function(t) sd(false_pos[t,] > eta_set[t])/sqrt(K))

power_mean <- rowMeans(power)
power_sd <- apply(power,1,sd)/sqrt(K)

#plots
par(mfrow=c(1,3), cex.lab = 2, cex.axis=2,mar=c(5,6,2,0), xpd=NA, oma = c(0, 0, 0, 4), mgp = c(3, 1, 0))
boxplot(t(false_pos),names=eta_set,xlab=expression(eta),ylab="Number of false discoveries")
boxplot(t(power),names=eta_set,xlab=expression(eta),ylab="Number of true discoveries")
plot(eta_set,false_pos_prob,ylab=expression(paste("P(false discovery > ",eta,")")),xlab=expression(eta), ylim=c(0,0.15),pch=16)
arrows(eta_set, false_pos_prob-false_pos_prob_sd, eta_set, false_pos_prob+false_pos_prob_sd, length=0.05, angle=90, code=3,lwd=2)
segments(0,alpha,10,alpha,lty=2)

#table
fd_mean <- rowMeans(false_pos)
fd_sd <- apply(false_pos,1,sd)/sqrt(K)
table_selection <- rbind(fd_mean,fd_sd,false_pos_prob, false_pos_prob_sd, power_mean, power_sd)
print(round(table_selection,4))



