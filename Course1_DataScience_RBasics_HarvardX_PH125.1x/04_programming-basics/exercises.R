#library(downloader)
#library(dplyr)

############################################
# QQ-plot Exercises
############################################

load("skew.RData")

par(mfrow = c(3,3))

for (i in 1:9) {
      qqnorm(dat[,i])
      qqline(dat[,i])
}

par(mfrow=c(1,1))


############################################
# Boxplot Exercises
############################################

head(InsectSprays)
summary(InsectSprays[InsectSprays$spray=='A',]$count)

# Try out two equivalent ways of drawing boxplots in R
# 1) using a formula
boxplot(count ~ spray, data = InsectSprays,
        xlab = "Type of spray", ylab = "Insect count",
        main = "InsectSprays data", varwidth = TRUE, col = "lightgray")

# 2) using split
boxplot(split(InsectSprays$count, InsectSprays$spray))

library(UsingR)
library(dplyr)
data(nym.2002, package="UsingR")
str(nym.2002)
time_male <- nym.2002[nym.2002$gender=="Male",]
time_fem <- nym.2002[nym.2002$gender=="Female",]
boxplot(time ~ gender, data = nym.2002,
        xlab = "gender", ylab = "finish time",
        main = "NY 2002 Marathon", varwidth = TRUE, col = "lightgray")
hist(time_male$time)
hist(time_fem$time)
summary(time_male$time)
summary(time_fem$time)

mypar(1,3)
boxplot(time_fem, time_male)
hist(time_fem,xlim=c(range( nym.2002$time)))
hist(time_male,xlim=c(range( nym.2002$time)))


############################################
# Scatterplot
############################################

# For males, what is the Pearson correlation between age and time to finish?
cor.test(time_male$age, time_male$time)
cor(time_fem$age, time_fem$time)

# If we interpret these correlations without visualizing the data, we would 
# conclude that the older we get, the slower we run marathons, regardless of 
# gender. Look at scatterplots and boxplots of times stratified by age groups 
# (20-25, 25-30, etc..). After examining the data, what is a more reasonable 
# conclusion
groups <- split(time_male$time,round(time_male$age)) 
boxplot(groups)


############################################
# Log Ratios
############################################

time = sort(nym.2002$time)

min(time) / median(time)
max(time) / median(time)

plot(time/median(time), ylim=c(1/4,4))
abline(h=c(1/2,1,2))

plot(log2(time/median(time)),ylim=c(-2,2))
abline(h=-1:1)


############################################
# Median, MAD, and Spearman Correlation
############################################

data(ChickWeight)
head(ChickWeight)
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)
chick = reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time",
                direction="wide")
chick = na.omit(chick)
head(chick)

# add an outlier 3000g: check impact on mean
mean(c(chick$weight.4,3000)) / mean(chick$weight.4)

# same exercise, using median
median(c(chick$weight.4,3000)) / median(chick$weight.4)

# same exercise, with stddev instead of mean
sd(c(chick$weight.4,3000)) / sd(chick$weight.4)

# same exercise, using median absolute deviation
mad(c(chick$weight.4,3000)) / mad(chick$weight.4)

plot( chick$Time, chick$weight, row=chick$Diet)

# same exercise, with Pearson correlation
cor(c(chick$weight.4,3000),c(chick$weight.21,3000)) / cor(chick$weight.4,chick$weight.21)


############################################
# Mann-Whitney-Wilcoxon
############################################

# Observation vs ranks: correct p-values via Z-score
x <- chick[chick$Diet==1, "weight.4"]
y <- chick[chick$Diet==4, "weight.4"]

# t-test
t.test(x,y)

# Wilcoxon test
wilcox.test(x,y)

# t-test with outlier
t.test(c(x,200),y)$p.value

# Wilcoxon test with outlier
wilcox.test(c(x,200),y)$p.value

# shifting series mean
library(rafalib)
mypar(1,3)
boxplot(x,y)
boxplot(x,y+10)
boxplot(x,y+100)

# difference of t-test
t.test(x,y+10)$statistic - t.test(x,y+100)$statistic

# difference of Wilcoxon test
wilcox.test(x,y+10)
wilcox.test(x,y+100)
# => Because the Wilcoxon works on ranks, once the two groups show complete 
# separation, that is all points from group 'y' are above all points from group
#'x', the statistic will not change, regardless of how large the difference 
# grows. Likewise, the p-value has a minimum value, regardless of how far apart
# the groups are. This means that the Wilcoxon test can be considered less 
# powerful than the t-test in certain contexts. In fact, for small sample sizes,
# the p-value can't be very small, even when the difference is very large.
# Example:
wilcox.test(c(1,2,3),c(4,5,6))$p.value
wilcox.test(c(1,2,3),c(400,500,600))$p.value
