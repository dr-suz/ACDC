######################################
## Example setting: Cauchy regression 
######################################
set.seed(99)
n <- 20 # also try 50, 100
N <- 100 # monte carlo sample size 
  
x1<- rep(NA, n) 
for(i in 1:n){
  if(i==1){
    x1[i] <- 0.5 + runif(1,20,30)}
  else{x1[i] <- x1[i-1] + runif(1,20,30)}
}

x2 <- rep(NA,n)
for(i in 1:n){
  if(i==1){
    x2[i] <- 0.2 * i + runif(1,40,60)}
  else{x2[i] <- 0.2*i + 0.3*x2[i-1] + runif(1,40,60)}
}

b0 <- 10; b1 <- 2; b2 <- 1 
response <- b0 + b1*x1 + b2*x2 + (rnorm(n,0,1)/rnorm(n,0,1)) # errors are Cauchy(0,1)

response; x1; x2
design <- cbind(rep(1,n),x1,x2)

beta_hat_obs <- solve(t(design) %*% design) %*% t(design) %*% response 

## Estimation via method of scoring/iterative max likelihood 
## (Alternatively, can also try estimation via minimum sum absolute errors, maybe for smaller n)

######################################
## divide into subsets and compute subset_est
# to replace divide() function from file :Divide_Recombine_Parallelise.R
######################################
divide_matrix<-function(obj,group_no,ini_index=0){
  obj_ind<-1:length(obj[,1])
  obj_ind_split<-split(obj_ind,ceiling(obj_ind/(length(obj_ind)/group_no))+ini_index)
  return(lapply(obj_ind_split,function(x)obj[x,]))
}

obs_data <- cbind(x1,x2,response)
k <- floor(sqrt(n))
subgroup <- divide_matrix(obs_data,k)

LSE <- function(data_matrix, n){
  # Input: data_matrix an n x 3 dimensional matrix with column 1 = x1, column 2 = x2, column 3 = y 
  design_tmp <- cbind(rep(1,n), data_matrix[,1], data_matrix[,2])
  response_tmp <- data_matrix[,3]
  return(solve(t(design_tmp) %*% design_tmp) %*% t(design_tmp) %*% response_tmp)
}

subgroup_est<-matrix(unlist(lapply(subgroup,LSE)), ncol=3, byrow=TRUE)
