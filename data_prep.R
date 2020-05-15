#---- Data Preparation for Evaluation of Correlation Clustering
#----------Simulated data generation using SimGen() ---------------
# Simulated data generation engine for Cluster 1
# A Multi-view data matrix generator
# using MVN(mu, sig)
# MASS package
require(MASS)
# Specifications of input parameters
# n <- number of observations
# p <- number of features for view 1
# q<- number of features for view 2
# n_gene <- number of genes that influences
# the traits been studied
# n_trait number of traits actually influenced
# the selected genes
simgen = function(n, p, q) {
# renaming variables
ng1 = 10
74
nt1 = 5
# view 1 feature means
mu_1 = rep(3, p) # fix mu
# view 1 feature std dev
vars_1 = rep(0.1,p)
sig_1 = diag(vars_1)
# view 1 matrix
v1 = mvrnorm(n, mu_1, sig_1, tol = 1e-6,
empirical = FALSE,
EISPACK = FALSE)
# indexing the genes that influences the traits
s1 = c(1:ng1) # first ng1 variables as target
# indexing the traits actually under influence
s2 = c(1:nt1) # first nt1 variables as influenced
# artificial effect size matrix,
# each row represents the coefficient vector for one of
# the influenced traits
eff = matrix(c(0,-1,1,-1,1,
75
0,-1,-1,-1,1,
1,0,1,-1,1,
1,0,-1,1,-1,
1,1,0,1,-1,
1,1,0,1,-1,
-1,-1,1,0,1,
-1,-1,-1,0,-1,
-1,1,1,-1,0,
-1,1,-1,1,0
), nrow = nt1, ncol = ng1)
# background noise of traits domiain
# view 2 feature means
mu_2 = rep(5, q)
# view 1 feature std dev
vars_2 = rep(0.1, q)
sig_2 = diag(vars_2)
# view 1 matrix
v2 = mvrnorm(n, mu_2, sig_2, tol = 1e-6,
empirical = FALSE,
EISPACK = FALSE)
# replace influenced traits under mapping
for (j in s2){
v2[,j] = v1[,s1] %*% eff[match(j,s2),]
76
+ rnorm(n, 0, 0.1)
}
# capturing output
view1 = v1
view2 = v2
total_view = as.data.frame(cbind(view1, view2))
# adding ID column to each row
total_view$ID = seq.int(nrow(total_view))
# renaming variables
s1 -> index_v1
s2 -> index_v2
eff -> betas
# Function Output
return(out = list(total_view, n, p, q,
index_v1, index_v2, view1, view2, betas))
}
# Data generation engine for Cluster 2
simgen_2 = function(n, p, q) {
# renaming variables
ng1 = 10
77
nt1 = 5
# view 1 feature means
mu_1 = rep(3, p) # fix mu
# view 1 feature std dev
vars_1 = rep(0.1,p)
sig_1 = diag(vars_1)
# view 1 matrix
v1 = mvrnorm(n, mu_1, sig_1, tol = 1e-6,
empirical = FALSE,
EISPACK = FALSE)
# indexing the genes that influences the traits
s1 = c((p-ng1+1):p) # last ng1 variables as target
# indexing the traits actually under influence
s2 = c((q-nt1+1):q) # last nt1 variables as influenced
# artificial effect size matrix,
# each row represents the coefficient vector for one of
# the influenced traits
78
eff = matrix(c(1,-1,1,-1,0,
-1,1,1,-1,0,
1,-1,1,0,1,
-1,1,1,0,-1,
1,-1,0,-1,1,
-1,1,0,-1,-1,
1,0,-1,1,1,
-1,0,-1,1,-1,
0,-1,-1,1,1,
0,1,-1,1,-1
), nrow = nt1, ncol = ng1)
# background noise of traits domiain
# view 2 feature means
mu_2 = rep(5, q)
# view 1 feature std dev
vars_2 = rep(0.1, q)
sig_2 = diag(vars_2)
# view 1 matrix
v2 = mvrnorm(n, mu_2, sig_2, tol = 1e-6,
empirical = FALSE,
EISPACK = FALSE)
# replace influenced traits under mapping
for (j in s2){
79
v2[,j] = v1[,s1] %*% eff[match(j,s2),]
+ rnorm(n, 0, 0.1)
}
# capturing output
view1 = v1
view2 = v2
total_view = as.data.frame(cbind(view1, view2))
# adding ID column to each row
total_view$ID = seq.int(nrow(total_view))
# renaming variables
s1 -> index_v1
s2 -> index_v2
eff -> betas
# Function Output
return(out = list(total_view, n, p, q
index_v1, index_v2, view1, view2, betas))
}
# Data generation for a two cluster multi-view dataset
set.seed(2018)
80
# Specifying multi-view dataset dimensions
n = 100
p = 1000
q = 10
ng = 10
nt = 5
# intrinsic cluster 1
sim_1 = simgen(n, p , q)
# capturing output view_1, view_2 matrices for cluster 1
v1_1 = as.data.frame(sim_1[7])
v2_1 = as.data.frame(sim_1[8])
mv_1 = cbind(v1_1, v2_1)
# capture the indices of target genes and influenced traits
s1_1 = c(1:ng)
s2_1 = c(1:nt)
# intrinsic cluster 2
sim_2 = simgen_2(n, p , q)
# capture output view_1, view_2 matrices for cluster 2
v1_2 = as.data.frame(sim_2[7])
v2_2 = as.data.frame(sim_2[8])
81
mv_2 = cbind(v1_2, v2_2)
# capture the indices of target genes and influenced traits
s1_2 = c((p-ng+1):q) # last ng1 variables as target
s2_2 = c((q-nt+1):q) # last nt1 variables as influenced
# combining two clusters to obtain simulated multi-view data
mvdata = as.data.frame(rbind(mv_1, mv_2))
# add tag for the intrinsic cluster of each row
mvdata[, "ID"] = c(1:nrow(mvdata))
# dimension check
dim(mvdata)