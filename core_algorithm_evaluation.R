
#------ rCCA correlation clustering ---------
#------ core algorithm and evaluation -------
require(CCA)
require(ggplot)
# regularization via cross validation
# Extracting View 1 and View 2 matrices
v1 = mvdata[, 1:1000]
v2 = mvdata[, 1001:1010]
# Cross validation function
reg_par = estim.regul(v1, v2)
# Regularization parameters
lam1 = reg_par$lambda1
lam2 = reg_par$lambda2
# rCCA correlation clustering engine
reg_clust = function(mvdata, p, q, k, iter){
# The true clustering scheme
truth = c(rep(1,100), rep(2,100))
# creating array for error_rate
error_rate = rep(0, iter)
83
mvdata$tg = sample(1:k, nrow(mvdata), replace = TRUE )
group = vector("list",length = k)
cca_out = vector("list",length = k)
U = vector(list, length = k)
V = vector(list, length = k)
slr_out = vector("list", length = k)
# the slope of V~U fit
a = vector("list", length = k)
# the intercept of V~U fit
b = vector("list", length = k)
for (i in c(1:iter)){
for (j in c(1:k)){
group[[j]] = mvdata[which(mvdata$tag ==j), ]
cca_out[[j]] = rcc(group[j][,1:p], group[[j]][
(p+1):(p+q)], lam1, lam2)
U[[j]] = as.matrix(group[[j]][, 1:p]) %*%
cca_out[[j]]$xcoef[,1]
V[[j]] = as.matrix(group[[j]][,(p+1):(p+q)]) %*%
cca_out[[j]]$ycoef[,1]
slr_out[[j]] = lm(V[[j]] ~ U[[j]])
84
b[[j]] = as.numeric(cef(slr_out[[j]])[1]) # intercept
a[[j]] = as.numeric(cef(slr_out[[j]])[2]) # slope
}
for (id in c(1:nrow(mvdata))){
item_v1 = mvdata[id,1:p]
item_v2 = mvdata[id, (p+1):(p+q)]
item_U = vector("list", length = k)
item_V = vector("list", length = k)
V_hat = vector("list", length = k)
dist = vector("list", length = k)
for (cl in c(1:k)){
item_U[[cl]] = as.matrix(item_v1) %*%
as.matrix( cca_out[[cl]]$xcoef[,1] )
item_V[[cl]] = as.matrix(item_v2) %*%
as.matrix( cca_out[[cl]]$ycef[,1] )
V_hat[[cl]] = a[[cl]]%*%item_U[[cl]] + b[[cl]]
dist[[cl]] = (V_hat[cl]- item_V[[cl]] )^2
85
}
new_tag_id = match( min(unlist(dist)), dist)
mvdata[id, ncol(mvdata)] = new_tag_id
} # reassign end
label_out = table(mvdata$tag == truth)
error_rate[i] = as.vector(label_out)[1]/nrow(mvdata)
} # iteration end
return(list(mvdata, cca_out, mvdata$tag, error_rate))
}
# Specifications of function output
# mvdata - the original multi-view dataset
# cca_out - the local cca models on each cluster
# mvdata$tag - the resulted clustering scheme
# error_rate - the error rate of clustering scheme
86
B.8 Evaluation of sCCA correlation clustering
# Spase CCA clustering with optimized penalties
require(PMA)
sparse_clust_calipena = function(mvdata, p, q, k, iter){
# the Truth clustering scheme
truth = c(rep(1,100), rep(2,100))
# creating array for error_rate
error_rate = rep(0, iter)
#1 randomly assign instances into k clusters
mvdata$tag = sample(1:k, nrow(mvdata), replace = TRUE)
group = vector("list",length = k)
cca_out = vector("list",length = k)
U = vector("list", length = k) # canonical covariates
V = vector("list", length = k)
slr_out = vector("list", length = k)
a = vector(list, length = k) # the slope of V~U fit
b = vector(list, length = k) # the intercept of V~U fit
87
for (i in c(1:iter)){
for(j in c(1:k)){
group[[j]] = mvdata[which(mvdata$tag ==j), ]
par = CCA.permute( group[[j]][,1:p], group[[j]][ ,(p+1):(p+q)] )
cca_out[[j]] = CCA(group[[j]][,1:p],group[[j]][,(p+1):(p+q)]
penaltyx = par$bestpenaltyx,
penaltyz = par$bestpenaltyz,
standardize = TRUE )
U[[j]] = as.matrix(group[[j]][ , 1:p]) %*% cca_out[[j]]$u
V[[j]] = as.matrix(group[[j][[ (p+1):(p+q)]) %*% cca_out[[j]]$v
slr_out[[j]] = lm(V[[j]] ~ U[[j]])
b[[j]] = as.numeric(coef(slr_out[[j]])[1]) # intercept
a[[j]] = as.numeric(coef(slr_out[[j]])[2]) # slope
}
# reassignment
for (id in c(1:nrow(mvdata))){
item_v1 = mvdata[id, 1:p]
item_v2 = mvdata[id, (p+1):(p+q)]
88
V_hat = vector("list", length = k)
item_U = vector("list", length = k)
item_V = vector("list", length = k)
dist = vector("list", length = k)
for (cl in c(1:k)){
item_U[[cl]] = as.matrix(item_v1) %*% cca_out[[cl]]$u
item_V[[cl]] = as.matrix(item_v2) %*% cca_out[[cl]]$v
V_hat[[cl]] = a[[cl]]%*%item_U[[cl]] + b[[cl]]
dist[[cl]] = (V_hat[[cl]]- item_V[[cl]] )^2
}
new_tag_id = match( min(unlist(dist)), dist)
mvdata[id, ncol(mvdata)] = new_tag_id
} # reassignment end
label_out = table(mvdata$tag = truth)
error_rate[i] = as.vector(label_out)[1]/nrow(mvdata)
} # iteration end
89
return(list(mvdata, cca_out, mvdata$tag, error_rate))
} # function end
# Specifications of function output
# mvdata - the original multi-view dataset
# cca_out - the local cca models on each cluster
# mvdata$tag - the resulted clustering scheme
# error_rate - the error rate of clustering scheme