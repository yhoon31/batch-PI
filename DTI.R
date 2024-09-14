source("functions.R")
y_cal <- read.table("y_cal.txt")$V1
y_test <- read.table("y_test.txt")$V1
muhat_cal <- read.table("muhat_cal.txt")$V1
muhat_test <- read.table("muhat_test.txt")$V1
set.seed(1)
perm <- sample(1:length(y_test)) #shuffle test points before grouping
y_test <- y_test[perm]
muhat_test <- muhat_test[perm]


##Predictive inference with PAC guarantee
residual_cal <- y_cal-muhat_cal
residual_test <- y_test-muhat_test

n <- 500
m <- 100
num_test_set <- 160
score_cal <- abs(residual_cal[sort(sample(1:length(residual_cal),n,replace=F))]) #construct calibration set of size n
score_cal_sorted <- sort(score_cal)

delta <- 0.1
alpha_set <- seq(0.05,0.3,by=0.05)

rank_bpi <- sapply(alpha_set,function(alpha) q_quantile(n,m,ceiling((1-delta)*m),alpha))
bpi_bound <- score_cal_sorted[rank_bpi]
conformal_bound <- score_cal_sorted[ceiling((1-delta)*(n+1))]
coverage_bpi <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))
coverage_conformal <- rep(0,num_test_set)

for(t in 1:num_test_set)
{
  score_test <- abs(residual_test[(m*(t-1)+1):(m*t)])
  coverage_bpi[t,] <- sapply(bpi_bound, function(b) mean(score_test <= b))
  coverage_conformal[t] <- mean(score_test <= conformal_bound)
}

coverage_bpi_pac = (coverage_bpi >= 1-delta)
cov_bpi_pac_mean <- colMeans(coverage_bpi_pac)
cov_bpi_pac_se <- apply(coverage_bpi_pac,2,sd)/sqrt(num_test_set)
cov_bpi_mean <- colMeans(coverage_bpi)
cov_bpi_se <- apply(coverage_bpi,2,sd)/sqrt(num_test_set)
cov_conformal_pac <- mean(coverage_conformal >= 1-delta)
cov_conformal_pac_se <- sd(coverage_conformal >= 1-delta)/sqrt(num_test_set)
cov_conformal_mean <- mean(coverage_conformal)
cov_conformal_se <- sd(coverage_conformal)/sqrt(num_test_set)

col1 <- rgb(30/255, 144/255, 255/255, alpha = 1)
col2 <- rgb(255/255, 20/255, 147/255, alpha = 1)

par(mfrow=c(1,2), cex.lab = 1.8, cex.main=1.8, cex.axis=1.8,mar=c(3,3,3,3), xpd=NA, oma = c(3, 3, 1, 12), mgp = c(3.5, 1.25, 0))
plot(1-alpha_set,cov_bpi_pac_mean,xlim=c(0.7,1),ylim=c(0.5,1),pch=16,xlab=expression(1-alpha),ylab="P(coverage > 0.9)",cex=0.8,col=col1)
segments(0.7,0.7,1,1,lty=2,lwd=3)
segments(0.7,cov_conformal_pac,1,cov_conformal_pac,lty=1,lwd=3,col=col2)
arrows(1-alpha_set, cov_bpi_pac_mean-cov_bpi_pac_se, 1-alpha_set, cov_bpi_pac_mean+cov_bpi_pac_se, length=0.05, angle=90, code=3,lwd=3,col=col1)

plot(1-alpha_set,cov_bpi_mean,xlim=c(0.7,1),ylim=c(0.87,0.96),pch=16,xlab=expression(1-alpha),ylab="Mean coverage rate",cex=0.8,col=col1)
segments(0.7,1-delta,1,1-delta,lty=2,lwd=3)
segments(0.7,cov_conformal_mean,1,cov_conformal_mean,lty=1,lwd=3,col=col2)
arrows(1-alpha_set, cov_bpi_mean-cov_bpi_se, 1-alpha_set, cov_bpi_mean+cov_bpi_se, length=0.05, angle=90, code=3,lwd=3,col=col1)

legend1 <- "batch PI"
legend2 <- "split conformal"
legend(1.05,0.93, legend = c(legend1, legend2),pch=c(16,NA),lty=c(NA,1), col=c(col1,col2), lwd=3, box.lty=0,cex=1.5,bg="transparent")

print(round(rbind(cov_bpi_pac_mean,cov_bpi_pac_se,cov_bpi_mean,cov_bpi_se),4))
round(c(cov_conformal_pac,cov_conformal_pac_se,cov_conformal_mean,cov_conformal_se),4)


##---------------------------------------------------------------------------------------------------------------------------------------------------------------------

##Selection with FWER control

##Predictive inference with PAC guarantee
set.seed(1)

n <- 2000
m <- 100
num_test_set <- 160
c <- 7
cal_ind <- sort(sample(1:length(y_cal),n,replace=F))
score_cal <- muhat_cal[cal_ind]*(y_cal[cal_ind ] <= c)
score_cal_sorted <- sort(score_cal)

alpha_set <- seq(0.05,0.3,by=0.05)

eta_set <- c(0,3,5)
rank_bpi_1 <- sapply(alpha_set,function(alpha) q_quantile(n,m,m-eta_set[1],alpha))
rank_bpi_2 <- sapply(alpha_set,function(alpha) q_quantile(n,m,m-eta_set[2],alpha))
rank_bpi_3 <- sapply(alpha_set,function(alpha) q_quantile(n,m,m-eta_set[3],alpha))
bpi_bound_1 <- score_cal_sorted[rank_bpi_1]
bpi_bound_2 <- score_cal_sorted[rank_bpi_2]
bpi_bound_3 <- score_cal_sorted[rank_bpi_3]
selection_bpi_1 <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))
selection_bpi_2 <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))
selection_bpi_3 <- matrix(0,nrow=num_test_set,ncol=length(alpha_set))

for(t in 1:num_test_set)
{
  score_test <- muhat_test[(m*(t-1)+1):(m*t)]
  y_test_trial <- y_test[(m*(t-1)+1):(m*t)]
  selection_bpi_1[t,] <- sapply(bpi_bound_1, function(b) sum((score_test > b)*(y_test_trial <= c)))
  selection_bpi_2[t,] <- sapply(bpi_bound_2, function(b) sum((score_test > b)*(y_test_trial <= c)))
  selection_bpi_3[t,] <- sapply(bpi_bound_3, function(b) sum((score_test > b)*(y_test_trial <= c)))
}

fd_1 <- colMeans(selection_bpi_1>eta_set[1])
fd_se_1 <- apply(selection_bpi_1>eta_set[1],2,sd)/sqrt(num_test_set)
fd_2 <- colMeans(selection_bpi_2>eta_set[2])
fd_se_2 <- apply(selection_bpi_2>eta_set[2],2,sd)/sqrt(num_test_set)
fd_3 <- colMeans(selection_bpi_3>eta_set[3])
fd_se_3 <- apply(selection_bpi_3>eta_set[3],2,sd)/sqrt(num_test_set)

par(mfrow=c(1,3), cex.lab = 1.8, cex.main=1.8, cex.axis=1.8,mar=c(3,3,3,3), xpd=NA, oma = c(3, 3, 1, 1), mgp = c(3.5, 1.25, 0))
plot(alpha_set,fd_1,ylab="P(# false discovery > 0)",ylim=c(0,0.4),pch=16,xlab=expression(alpha))
arrows(alpha_set, fd_1-fd_se_1, alpha_set, fd_1+fd_se_1, length=0.05, angle=90, code=3,lwd=3,col=col1)
segments(0.05,0.05,0.3,0.3,lty=2,lwd=3)
plot(alpha_set,fd_2,ylab="P(# false discovery > 3)",ylim=c(0,0.4),pch=16,xlab=expression(alpha))
arrows(alpha_set, fd_2-fd_se_2, alpha_set, fd_2+fd_se_2, length=0.05, angle=90, code=3,lwd=3,col=col1)
segments(0.05,0.05,0.3,0.3,lty=2,lwd=3)
plot(alpha_set,fd_3,ylab="P(# false discovery > 5)",ylim=c(0,0.4),pch=16,xlab=expression(alpha))
arrows(alpha_set, fd_3-fd_se_3, alpha_set, fd_3+fd_se_3, length=0.05, angle=90, code=3,lwd=3,col=col1)
segments(0.05,0.05,0.3,0.3,lty=2,lwd=3)


round(rbind(fd_1,fd_se_1,fd_2,fd_se_2,fd_3,fd_se_3),4)
