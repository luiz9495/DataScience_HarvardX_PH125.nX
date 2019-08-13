############################################
# T-STAT, P-VALUE, Type-I and II error,POWER
############################################

library(downloader)
library(dplyr)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)

N=5
set.seed(1)
dat.ns <- sample(bwt.nonsmoke , N)
dat.s <- sample(bwt.smoke , N)

X.ns <- mean(dat.ns)
sd.ns <- sd(dat.ns)

X.s <- mean(dat.s)
sd.s <- sd(dat.s)

sd.diff <- sqrt(sd.ns^2/N+sd.s^2/N)
tval <- (X.s - X.ns)/sd.diff
tval

pval <- 1-pnorm(abs(tval))+pnorm(-abs(tval))
pval

qnorm(0.995) * sd.diff
qt(0.995,2*N-2) * sd.diff

t.test(dat.ns,dat.s)

N=5
set.seed(1)
foo <- function(N, alpha=0.05) {
      dat.ns <- sample(bwt.nonsmoke , N)
      dat.s <- sample(bwt.smoke , N)
      pval <- t.test(dat.ns,dat.s)$p.value
      pval < alpha
}

res <- replicate(10000, foo(N))
mean(res)

for (N in c(30,60,90,120)) {
      res <- replicate(10000, foo(N, 0.01))
      print(mean(res))
}

Ns=c(30,60,90,120)
res <- sapply(Ns, function(N){
      set.seed(1)
      rejects <- replicate(10000,{
            dat.ns <- sample(bwt.nonsmoke , N)
            dat.s <- sample(bwt.smoke , N)
            t.test(dat.s, dat.ns)$p.value < 0.01
      })
      mean(rejects)
})
Ns[ which.min( abs( res - .8) ) ] 


################################
# MC Simulations
################################

set.seed(1)
N <- 5
B <- 1000
tstats <- replicate(B, {
      dat <- rnorm(N)
      sqrt(N) * mean(dat) / sd(dat)
      })
mean(tstats > 2)

1-pt(2,df=4)


library(rafalib)
mypar(3,2)

set.seed(1)
Ns<-seq(5,30,5)
B <- 1000
mypar(3,2)
LIM <- c(-4.5,4.5)
for(N in Ns){
      ts <- replicate(B, {
            X <- rnorm(N)
            sqrt(N)*mean(X)/sd(X)
      })
      ps <- seq(1/(B+1),1-1/(B+1),len=B)
      qqplot(qt(ps,df=N-1),ts,main=N,
             xlab="Theoretical",ylab="Observed",
             xlim=LIM, ylim=LIM)
      abline(0,1)
} 


set.seed(1)
Ns<-seq(5,30,5)
B <- 1000
mypar(3,2)
LIM <- c(-4.5,4.5)
for(N in Ns){
      ts <- replicate(B,{
            x <- rnorm(N)
            y <- rnorm(N)
            t.test(x,y, var.equal = TRUE)$stat
      })
      ps <- seq(1/(B+1),1-1/(B+1),len=B)
      qqplot(qt(ps,df=2*N-2),ts,main=N,
             xlab="Theoretical",ylab="Observed",
             xlim=LIM, ylim=LIM)
      abline(0,1)
}

set.seed(1)
N <- 15
B <- 10000
tstats <- replicate(B,{
      X <- sample(c(-1,1), N, replace=TRUE)
      sqrt(N)*mean(X)/sd(X)
})
ps=seq(1/(B+1), 1-1/(B+1), len=B) 
qqplot(qt(ps,N-1), tstats, xlim=range(tstats))
abline(0,1)
#The population data is not normal thus the theory does not apply.
#We check with a Monte Carlo simulation. The qqplot shows a large tail. 
#Note that there is a small but positive chance that all the X are the same.
##In this case the denominator is 0 and the t-statistics is not defined


set.seed(1)
N <- 1000
B <- 10000
tstats <- replicate(B,{
      X <- sample(c(-1,1), N, replace=TRUE)
      sqrt(N)*mean(X)/sd(X)
})
ps=seq(1/(B+1), 1-1/(B+1), len=B) 
qqplot(qt(ps,N-1), tstats, xlim=range(tstats))
abline(0,1)
#With N=1000, CLT kicks in and the t-statistic is approximated with normal 0,1
##Furthermore, t-distribution with df=999 and normal are practically the same.


set.seed(1)
Ns <- seq(5,45,5)
library(rafalib)
mypar(3,3)
for(N in Ns){
      medians <- replicate(10000, median ( rnorm(N) ) )
      title <- paste("N=",N,", avg=",round( mean(medians), 2) , ", sd*sqrt(N)=", round( sd(medians)*sqrt(N),2) )
      qqnorm(medians, main = title )
      qqline(medians)
}
##there is an asymptotic result that says SD is sqrt(N*4*dnorm(0)^2)

#############################################
# PERMUTATIONS
#############################################

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist


N=10
set.seed(1)
nonsmokers <- sample(bwt.nonsmoke , N)
smokers <- sample(bwt.smoke , N)
obs <- mean(smokers) - mean(nonsmokers)

B <- 1000
set.seed(1)
null <- replicate (B, {
      dat <- c(smokers,nonsmokers)
      shuffle <- sample( dat )
      smokersstar <- shuffle[1:N]
      nonsmokersstar <- shuffle[(N+1):(2*N)]
      mean(smokersstar)-mean(nonsmokersstar)
      })
( sum( abs(null) >= abs(obs)) +1 ) / ( length(null)+1 ) 
##we add the 1s to avoid p-values=0 but we also accept:


#######################################################
#Association Tests
#######################################################

d = read.csv("assoctest.csv")

tab = table(d$allele, d$case)

chisq.test(tab) 
fisher.test(tab)
