#
# Replication script of the simulation study
# for "High-level Causation and Causal Inference"
# by Lorenzo Casini and Alessio Moneta
#
# The script replicates the study illustrated in Figure 6 of the paper
#
#
# These packages must be installed and uploaded before running the script:
#
library(MASS)
#library(ppcor)
#library(dHSIC)
#library(OOmisc)
library(sm)
#
#  Function to replicate data for Model I
# i.e.  model with ill-defined variable
#
nid_ill_2ex<-function(ssize, mruns, sseed=46){ # two experiments (two pairs of coefficients)
  cc<-1:mruns
  alphas<-seq(from=1, to=3, length=2)
  for(i in 1:mruns){
    k<-4
    set.seed(i)
    E<-mvrnorm(ssize,rep(0,k), Sigma= diag(k))
    I<-E[,1]
    set.seed(i)
    alpha<-sample(alphas, size=1, replace=TRUE)
    X1<-alpha*I + E[,2]
    X2<-(1-alpha)*I +E[,3]
    X= X1 + X2 
    Y<- X1 - X2 +E[,4]
    rm<-lm(Y~X)
    c0<-rm$coefficients[2]
    cc[i]<-c0
  }
  cc
}
#
#  Function to replicate data for Model II
#  i.e. model where X is well-defined and 
# X1 and X2 are common causes of both X and Y
#
wd_mod2_2ex<-function(ssize, mruns, sseed=46){ 
  cc<-1:mruns
  #alphas<-seq(from=0, to= 0.5, length=2)
  alphas<-seq(from=1/3, to= 5/3, length=2)
  for(i in 1:mruns){
    k<-5
    set.seed(i)
    E<-mvrnorm(ssize,rep(0,k), Sigma= diag(k))
    I<-E[,1]
    set.seed(i)
    alpha<-sample(alphas, size=1, replace=TRUE)
    X1<-alpha*I + E[,2]
    X2<-alpha*I +E[,3]
    X=  I +E[,5]
    Y<- X1 + X2 +E[,4]
    rm<-lm(Y~X)
    c0<-rm$coefficients[2]
    cc[i]<-c0
  }
  cc
}

#### PLOT FIGURE 6 LEFT #######

cc<-nid_ill_2ex(ssize=100, mruns=1000, sseed = 46)
sm100<-sm.density(cc)
cc<-nid_ill_2ex(ssize=500, mruns=1000, sseed = 46)
sm500<-sm.density(cc, display="none")
cc<-nid_ill_2ex(ssize=1000, mruns=1000, sseed = 46)
sm1000<-sm.density(cc, display="none")
cc<-nid_ill_2ex(ssize=200, mruns=1000, sseed = 46)
sm200<-sm.density(cc, display="none")
pdf("plot1B.pdf", width=4, height=4)
matplot(sm100$eval.points, sm100$estimate, t="l", 
        xlim=c(-0.5,3), ylim=c(0,0.8),
        xlab=expression(hat(beta)),
        #  ylab="probability density function", 
        ylab="",
        main="X ill-defined wrt Y (Model I)", cex.main=1, lty=3)
matplot(sm200$eval.points, sm200$estimate, t="l", col="red",add=TRUE, lty=2)
matplot(sm1000$eval.points, sm1000$estimate, t="l",add=TRUE, col="blue", lty=1)
legend(-0.5,0.8,c("n=1000", "n=200", "n=100"),lty=c(1,2,3),bty="n", cex=1, text.col=c("blue", "red", 1), col=c("blue", "red", 1))
dev.off()

#### PLOT FIGURE 6 RIGHT #######

cc<-wd_mod2_2ex(ssize=100, mruns=1000, sseed = 46)
sm100<-sm.density(cc)
cc<-wd_mod2_2ex(ssize=200, mruns=1000, sseed = 46)
sm200<-sm.density(cc)
cc<-wd_mod2_2ex(ssize=500, mruns=1000, sseed = 46)
sm500<-sm.density(cc)
cc<-wd_mod2_2ex(ssize=1000, mruns=1000, sseed = 46)
sm1000<-sm.density(cc)

pdf("plot1C.pdf", width=4, height=4)
matplot(sm100$eval.points, sm100$estimate, t="l", 
        #  xlim=c(-0.5,1), ylim=c(0,2),
        xlim=c(-0.5,3), ylim=c(0,0.8), 
        xlab=expression(hat(beta)),
        #  ylab="probability density function", 
        ylab="",
        main="X well-defined wrt Y (Model II) ", cex.main=1, lty=3)
matplot(sm200$eval.points, sm200$estimate, t="l", col="red",add=TRUE, lty=2)
matplot(sm1000$eval.points, sm1000$estimate, t="l",add=TRUE, col="blue", lty=1)
legend(-0.5,0.8,c("n=1000", "n=200", "n=100"),lty=c(1,2,3),bty="n", cex=1, text.col=c("blue", "red", 1), col=c("blue", "red", 1))
#legend(-0.5,2,c("n=1000", "n=200", "n=100"),lty=c(1,2,3),bty="n", cex=1, text.col=c("blue", "red", 1), col=c("blue", "red", 1))
dev.off()


