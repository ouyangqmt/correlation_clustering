#---- Genetic expression simulator---------------
#---- as function 'sim_gen()'-------------------

#---- required packages -------

require(MASS)
require(ggplot)

#---- set seed ----
set.seed(2018)

#-----------------------------
# Parameters specification

# n: number of observations
# p: number of features for view 1
# q: number of features for view 2
# ng: number of genes that influences the interested traits
# nt: number of traits influenced the truly associated genes
SimGen = function(n, p, q) {
# view 1 feature means
mu_1 = rep(3, p)
# view 1 feature std dev
vars_1 = rep(0.1,p)
sig_1 = diag(vars_1)
# view 1 data matrix
v1 = mvrnorm(n, mu_1, sig_1, tol = 1e-6,
empirical = FALSE,
EISPACK = FALSE)
# indexing the genes that influences the traits
55
s1 = c(1:ng) # first ng variables as target
# indexing the traits actually under influence
s2 = c(1:nt) # first nt variables as influenced
# artificial effect size matrix, each row represents the
# coefficient vector for one of the influenced traits
eff = matrix(c(0,-1,1,-1,1,
0,-1,-1,-1,1,
1,0,1,-1,1,
1,0,-1,1,-1,
1,1,0,1,-1,
1,1,0,1,-1,
-1,-1,1,0,1,
-1,-1,-1,0,-1,
-1,1,1,-1,0,
-1,1,-1,1,0
), nrow = nt, ncol = ng)
# background noise of traits domain
# view 2 feature means
mu_2 = rep(5, q)
# view 1 feature std dev
vars_2 = rep(0.1, q)

sig_2 = diag(vars_2)
# view 2 data matrix
v2 = mvrnorm(n, mu_2, sig_2, tol = 1e-6,
empirical = FALSE,
EISPACK = FALSE)
# replace influenced traits under mapping
for (j in s2){
v2[,j] = v1[,s1] %*% eff[match(j,s2),]
+ rnorm(n, 0, 0.1)
}
# combining view 1 and view 2 for the multi-view dataset
total_view = as.data.frame(cbind(view1, view2))
# adding row ID
total_view$ID = seq.int(nrow(total_view))
# Renaming variables
s1 -> index_v1
s2 -> index_v2
eff -> betas
v1 -> view_1
v2 -> view_2
# Function Output

return(out = lsit(total_view, n, p, q
index_v1, index_v2, v1, v2, betas))
}

