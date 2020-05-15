#---- Shrinkage Regularization ------------
#---- rCCA with shrinkage regularization-----
# This evaluation requires R- package ‘‘mixOmics’’
# and installation of XQuartz app in Mac OS
# or X11 in Windows
require(mixOmics)
require(ggplot)


start.time = Sys.time()
# Output of Regularized CCA
rcca_out = rcc(v1, v2, method = "shrinkage")
Obtained canonical correlations
rcca_cor = rcca_out$cor
# Obtained first degree projection vectors (loadings)
x_load = rcca_out$loadings$X[,1]
y_load = rcca_out$loadings$Y[,1]
# recall target genes and influenced traits

par(mfrow=c(1,2))
# Plot of view 1 loading with target genes
# marked in red
x = seq(1,1000,1)
plot(x, x_load, pch = ifelse(x%in%c(1:10), 17, 1),col =
ifelse(x%in%c(1:10), "red", "black"),
cex = 0.8, xlab = "View 1 variables",
ylab = "Canonical Loadings",
main = "Canonical loadings of View 1 Variables")
abline(h = 0, col = "green")
# Plot of absolute value of view 1 loading, target genes
# marked in blue
y = seq(1,10,1)
plot(y, y_load,pch = ifelse(x%in%c(1:5), 17, 1) ,col =
ifelse(y%in%c(1:5), "blue", "black"),
cex = 0.8, xlab = "View 2 variables",
ylab = "Canonical Loadings",
main = "Canonical loadings of View 2 Variables")
abline(h = 0, col = "green")
end.time = Sys.time()

# -----Assessment of performance ------------
# view 1
# mean level of effect size of noise variable in view 1
mean_noise_1 = mean(abs(x_load[-s1]))
# return Distinctiveness of interactive variables
abs(x_load[s1])>mean_noise_1
# mean level of effect size of target variables in view 1
mean_int_1 = mean(abs(x_load[s1]))
# Degree of separation in View 1
dos_1 = mean_noise_1/mean_int_1
dos_1
any(x_load == 0)
# mean and standard deviation of noise variable in view 1
mean(x_load[s1])
sd(x_load[s1])
# view 2
# mean level of effect size of noise variable in view 2
63
mean_noise_2 = mean(abs(y_load[-s2]))
# return Distinctiveness of interactive variables
abs(y_load[s2])>mean_noise_2
# mean level of effect size of target variables in view 2
mean_int_2 = mean(abs(y_load[s2]))
# degree of separation in view 2
dos_2 = mean_noise_2/mean_int_2
any(y_load == 0)
# mean and standard deviation of noise variables in view 2
mean(y_load[s2])
sd(y_load[s2])
# running time
run.time = end.time-start.time
64
#---- Evaluation of Regularized CCA - via CrossValidation Regularization
#---- Regularized Canonical Correlation Analysis--------
#---- via Cross-Validation Regularization -----
require(CCA)
start.time = Sys.time()
# Cross-Validation Regularization via
# estim.regul() required for high dimensional data
# format - estim.regul(X, Y, grid1 = NULL,
# grid2 = NULL, plt = TRUE)
# grid: if NULL - grid1, grid2 vector use
# seq(0.001, 1, length = 5) as default otherwise specify
# grid values ie. c(0.01,0.5)
# plt: logic, whether the CV heatmap should be plotted
# Regularization parameters
# reg_par = estim.regul(v1, v2)
lam1 = reg_par$lambda1
lam2 = reg_par$lambda2
# Implement rCCA with previously obtained
# regularization parameters
rcca_out = rcc(v1, v2, lam1, lam2)
65
rcca_out
# Obtained canonical correlations
rcca_cor = rcca_out$cor
rcca_cor
# Obtained first degree projection vectors (loadings)
x_load = rcca_out$xcoef[,1]
y_load = rcca_out$ycoef[,1]
# recall target genes and influenced traits
s1
s2
# plot layout specification for visual inspection
par(mfrow=c(1,2))
# Plot of absolute value of view 1 loading, target genes
# marked in red
x = seq(1,1000,1)
plot(x, x_load, pch = ifelse(x%in%c(1:10), 17, 1),col =
ifelse(x%in%c(1:10), "red", "black"),
cex = 0.8, xlab = "View 1 variables",
ylab = "Canonical Loadings",
main = "Canonical loadings of View 1 Variables")
abline(h = 0, col = "green")
66
# Plot of absolute value of view 1 loading, target genes
# marked in blue
y = seq(1,10,1)
plot(y, y_load,pch = ifelse(x%in%c(1:5), 17, 1) ,col =
ifelse(y%in%c(1:5), "blue", "black"),
cex = 0.8, xlab = "View 2 variables",
ylab = "Canonical Loadings",
main = "Canonical loadings of View 2 Variables")
abline(h = 0, col = "green")
end.time = Sys.time()
# -----Assessment of performance ------------
# view 1
# mean level of effect size of noise variable in view 1
mean_noise_1 = mean(abs(x_load[-s1]))
# return Distinctiveness of interactive variables
abs(x_load[s1])>mean_noise_1
# mean level of effect size of target variables in view 1
mean_int_1 = mean(abs(x_load[s1]))
# Degree of separation in View 1
dos_1 = mean_noise_1/mean_int_1
67
dos_1
any(x_load == 0)
# mean and standard deviation of noise variable in view 1
mean(x_load[-s1])
sd(x_load[-s1])
# view 2
# mean level of effect size of noise variable in view 2
mean_noise_2 = mean(abs(y_load[-s2]))
# return Distinctiveness of interactive variables
abs(y_load[s2])>mean_noise_2
# mean level of effect size of target variables in view 2
mean_int_2 = mean(abs(y_load[s2]))
# degree of separation in view 2
dos_2 = mean_noise_2/mean_int_2
any(y_load == 0)
# mean and standard deviation of noise variables in view 2
mean(y_load[-s2])
68
sd(y_load[-s2])
# running time
run.time = end.time-start.time
69
B.5 Evaluation of Sparse CCA
require(PMA)
# timing starts
start.time = Sys.time()
# Sparse CCA output
sparse_out = CCA(v1,v2,"standard","ordered",
standardize = TRUE)
#Obtained canonical correaltions
sparse_cor = sparse_out$cors
sparse_cor
# Obtained first degree view 1 and view 2 loadings
x_load = sparse_out$u
y_load = sparse_out$v
# graph layout specification
par(mfrow=c(1,2))
# plot of absolute value of view 1 loading, target
# genes marked in red
x = seq(1,1000,1)
70
plot(x, x_load, col = ifelse(x %in% c(1:10), "red", "black"),
cex = 0.8, xlab = "Genes",
ylab = "Canonical loadings",
main = "Canonical loadings of View 1")
abline(h = 0, col = "green")
# plot of absolute value of view 2 loading, target
# genes marked in blue
y = seq(1,10,1)
plot(y, y_load, col = ifelse(y %in% c(1:5), "blue", "black"),
cex = 0.8, xlab = "Traits",
ylab = "Canonical loadings",
main = "Canonical loadings of view 2")
# horizontal reference
abline(h = 0, col = "green")
#timing ends
end.time = Sys.time()
# -----Assessment of performance ------------
# view 1
# mean level of effect size of noise variable in view 1
mean_noise_1 = mean(abs(x_load[-s1]))
71
# return Distinctiveness of interactive variables
abs(x_load[s1])>mean_noise_1
# mean level of effect size of target variables in view 1
mean_int_1 = mean(abs(x_load[s1]))
# Degree of separation in View 1
dos_1 = mean_noise_1/mean_int_1
dos_1
any(x_load == 0)
# mean and standard deviation of noise variable in view 1
mean(x_load[-s1])
sd(x_load[-s1])
# view 2
# mean level of effect size of noise variable in view 2
mean_noise_2 = mean(abs(y_load[-s2]))
# return Distinctiveness of interactive variables
abs(y_load[s2])>mean_noise_2
# mean level of effect size of target variables in view 2
mean_int_2 = mean(abs(y_load[s2]))
72
# degree of separation in view 2
dos_2 = mean_noise_2/mean_int_2
any(y_load == 0)
# mean and standard deviation of noise variables in view 2
mean(y_load[-s2])
sd(y_load[-s2])
# running time
run.time = end.time-start.time