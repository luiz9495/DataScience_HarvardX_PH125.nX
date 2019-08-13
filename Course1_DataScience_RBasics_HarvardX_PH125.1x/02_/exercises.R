n <- 5
p <- 1/2
var <- p*(1-p)/n
set.seed(1)
foo <- function(x = sample(1:1/p, n, replace=TRUE), z = (mean(x==1/p) - p) / sqrt(var)) z
res <- replicate(10000, foo())
mean(abs(res)>2)
hist(res)
qqnorm(res)
abline(0,1)

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
if(!file.exists("femaleMiceWeights.csv")) download(url,destfile=filename)
dat <- read.csv(filename)

X <- filter(dat, Diet=="chow") %>% select(Bodyweight) %>% unlist
Y <- filter(dat, Diet=="hf") %>% select(Bodyweight) %>% unlist

mean(X)