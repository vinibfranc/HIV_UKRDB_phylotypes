#install.packages("betareg")
library(betareg)
library(dplyr)
library(ggplot2)

packageVersion("betareg")
# [1] ‘3.2.1’

# 11 seqs close relative to PT40, only one close relative to PT69 (the VOI with signig VL and CD4), and 3 close relative to PT133
bvoi = c(2.6, 3.4, 3.8, 4, 17.1, 20.9, 20.9, 50.5, 86.9, 89.1, 93.8,
									50.3, 
									40.7, 63.4, 67
)/100
median(bvoi)#;IQR(bvoi)
# 0.407
calc_iqr(bvoi) 
# 0.11, 0.65
# proportion of sequences X-tropic
length(bvoi[bvoi<0.2]) / length(bvoi)
# 0.3333333

bother =c(0.1, 0.7, 0.7, 1.7, 1.7, 1.7, 1.7, 1.7, 2.5, 2.6, 2.6, 2.9, 4.7, 4.8, 5, 5, 5.3, 
										5.3, 6.9, 6.9, 6.9, 6.9, 9.6, 10.1, 11.4, 11.4, 13.2, 13.2, 15, 17.1, 17.1, 17.1, 
										18, 19.5, 20.8, 24.7, 25.2, 30.1, 34.6, 35.3, 35.6, 36.2, 36.9, 36.9, 37.7, 38, 38.8,
										39.3, 39.6, 39.7, 39.8, 41.4, 42.2, 42.8, 43, 43, 44.2, 44.2, 46.8, 47.8, 49, 50.2, 56,
										57.7, 58.3, 58.6, 59.2, 60.5, 67.3, 68.9, 73.7, 74, 74.6, 75.9, 79.5, 82, 86.2, 86.5, 86.5, 
										89, 90.4, 90.9, 93.6, 95.2, 95.7, 95.7, 98.9, 99, 99.4)/100 
median(bother);#IQR(bother)
# 0.377
calc_iqr(bother)
# 0.1, 0.59

# number of sequences X-tropic
length(bother[bother<0.2]) / length(bother)
# 0.3820225

wilcox.test( bvoi, bother )
# 	Wilcoxon rank sum test with continuity correction
# 
# data:  bvoi and bother
# W = 699.5, p-value = 0.7707
# alternative hypothesis: true location shift is not equal to 0

d = data.frame( score = c( bvoi, bother ), group = c( rep('VOI', length(bvoi)), rep('Non-VOI', length(bother))))
head(d)
f = betareg( score ~ group, data = d )
summary(f) 
# Call:
# betareg(formula = score ~ group, data = d)
# 
# Quantile residuals:
# 	Min      1Q  Median      3Q     Max 
# -2.2167 -0.8062  0.0480  0.5379  2.4941 
# 
# Coefficients (mean model with logit link):
# 	Estimate Std. Error z value Pr(>|z|)   
# (Intercept) -0.38279    0.12704  -3.013  0.00259 **
# 	groupVOI     0.03539    0.32848   0.108  0.91420  
# 
# Phi coefficients (precision model with identity link):
#       Estimate Std. Error z value Pr(>|z|)    
# (phi)    1.515      0.174   8.706   <2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
# 
# Type of estimator: ML (maximum likelihood)
# Log-likelihood:  11.36 on 3 Df
# Pseudo R-squared: 0.0001225
# Number of iterations: 13 (BFGS) + 1 (Fisher scoring) 

hist(d$score[d$group=="VOI"], breaks=100)
hist(d$score[d$group=="Non-VOI"], breaks=100)

ggplot(d, aes(x = score, fill = group)) +
	geom_histogram(binwidth = 0.05, position = "dodge", color = "black", alpha = 0.7) +
	scale_fill_manual(values = c("VOI" = "blue", "Non-VOI" = "red")) +
	scale_x_continuous(breaks=seq(0,1,by=0.1)) +
	labs(x = "False Positivity Rate (FPR)", y = "Count", fill = "Group") +
	theme_classic() +
	geom_vline(aes(xintercept = 0.2), color = "red", linetype = "dashed")