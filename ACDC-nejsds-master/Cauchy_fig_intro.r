library(matrixStats) # For function 'rowMads'
library(snowfall) # For parallelisation
library(parallel) # For parallelisation
library(pbapply) # For function 'pblapply'
library(mcmc)
library(tidyverse)
library(gridExtra)


remove(list=ls())

#setwd('/Users/sthornt1/Downloads/core-functions-main')
#setwd('/Users/m/Downloads/W. Li Code/core-functions-main')
source('Auxiliary_Functions.R')
source('Divide_Recombine_Parallelise.R') 
source('Functions_of_RejABC.R')
source('Functions_of_RegABC.R')
source('Functions_of_ABC_PMC.R')

#setwd('/Users/m/Downloads')
#setwd("/Users/sthornt1/Downloads/cauchy")
source('Cauchy_example_functionsv2.R')



parameter_true<-c(location=10,scale=0.55)
platform = 'Mac'
tolerance_percentile = 0.005 

n = 100
N = 5e4 
N_MCMC<-1e5 
prior_choice = 'Cauchy'

# generate observed data 
#set.seed(404)  #ISSBC_params[[1]], post_samp[[1]]
set.seed(400)
obs_dat = rcauchy(n, parameter_true[1], parameter_true[2])

group_no<-sqrt(n)
subgroup<-divide(obs_dat,group_no)

data_i = 1 
prior_df = 1 

test_size = 1 

ISABC_params = ABC_params = post_samps = MCMC_accept = list()

logposterior_Cauchy<-function(params,obs){ 
	loglik_val<-sum(dcauchy(obs,location=params[1],scale=params[2],log=TRUE))
	if(is.na(loglik_val)) loglik_val<--Inf
	prior_val<-1/params[2]
	if(prior_val>0) logprior_val<-log(prior_val)
	if(prior_val<=0) logprior_val<--Inf
	return(loglik_val+logprior_val)
}


for(ii in 1:2){ #parameter setting
parameter_setting <- ii

tested_parameters<-matrix(rep(parameter_true,test_size),ncol=2,byrow=T) # a test_size*2 matrix

if(parameter_setting==1) prior_location<-parameter_true[1]
if(parameter_setting==2) prior_location<-parameter_true[1]

prior_scale<-rep(10,length(prior_location))

# Simulate the parameter values and calculate their prior and proposal densities, conditional on the given parameter setting and proposal distribution
if(parameter_setting==1){
				param_est<-1
				subgroup_medians<-unlist(lapply(subgroup,median)) # a N*1 subgroup_medians
				simulated_parameters<-cbind(simulate_from_r2(N,subgroup_medians),tested_parameters[data_i,2])
				proposal_densities<-density_r2(simulated_parameters[,1],subgroup_medians)
				if(prior_choice=='jeffery') prior_densities<-rep(1,N)
				if(prior_choice=='Cauchy') prior_densities<-dt((simulated_parameters[,1]-prior_location)/prior_scale,df=prior_df)
				tested_summary <- median(obs_dat)
}
if(parameter_setting==2){
				param_est<-1
				subgroup_means<-unlist(lapply(subgroup,mean)) # a N*1 vector
				simulated_parameters<-cbind(simulate_from_r2(N,subgroup_means),tested_parameters[data_i,2])
				proposal_densities<-density_r2(simulated_parameters[,1],subgroup_means)
				if(prior_choice=='jeffery') prior_densities<-rep(1,N)
				if(prior_choice=='Cauchy') prior_densities<-dt((simulated_parameters[,1]-prior_location)/prior_scale,df=prior_df)
				tested_summary <- mean(obs_dat)
}

simulated_observations<-simulate_pseudo_datasets(simulated_parameters,n) 
simulated_summaries<-cal_summaries(simulated_observations,setting=parameter_setting,platform=platform)


reference_tables<-list(parameters=simulated_parameters,summaries=simulated_summaries)
weights_of_sample<-prior_densities/proposal_densities

############rej ABC (rejABC)#################
method_i<-1 # 
results_rejABC<-rejABC(reference_tables,t(tested_summary),
					tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform)
ABC_params[[ii]]  = results_rejABC$set_1$parameters %>% data.frame
names(ABC_params[[ii]]) = c("theta", "tau")


###########importance sampling regression ABC (ISregABC)###################
method_i<-2 # 
results_ISregABC<-regABC(reference_tables,t(tested_summary),
					tolerance_percentile=tolerance_percentile,paral=FALSE,platform=platform,
					weights_of_sample=weights_of_sample)
ISABC_params_tmp  = results_ISregABC$set_1$parameters %>% data.frame
names(ISABC_params_tmp) = c("theta", "tau")
ISABC_params[[ii]] = ISABC_params_tmp %>% mutate(weights = results_ISregABC$set_1$weights/sum(results_ISregABC$set_1$weights))    ### maybe I need to weight the parameters? 


############target posterior #####################
MCMC_results<-metrop(logposterior_Cauchy,initial=parameter_true,scale=diag(c(0.03,0.02)),nbatch=N_MCMC,obs=obs_dat)
MCMC_accept[[ii]]<-MCMC_results$accept
post_samps[[ii]]<-MCMC_results$batch[-(1:1e4),] %>% data.frame
names(post_samps[[ii]]) = c("theta", "tau")
}


##save.image('ABC-median.Rdata')
#load('ABC-median.Rdata')
#save.image('ABC-mean.Rdata')
#load('ABC-mean.Rdata')
#abc <- ABC_params[[2]]
#post <- post_samps[[2]]
#remove(list=setdiff(ls(),c('abc', 'post')))
###################################################
#dev.off()
#p1 <- ggplot() + geom_density(data=ABC_params[[1]],aes(x=theta)) +
#	  			 geom_density(data=post_samps[[1]],aes(x=theta), col="gray") +
#				 xlim(9,11) + labs(x=expression(theta), y='Density') + 
#				 theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
#p2 <- ggplot() +  geom_density(data=ABC_params[[2]],aes(x=theta)) +
#				  geom_density(data=post_samps[[2]],aes(x=theta), col="gray") +
#				  xlim(5,12) + labs(x=expression(theta), y=' ') + 
#				  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

#pdf(paste0("Cauchy_ABC_fig_n",n, ".pdf"))
#grid.arrange(p1,p2,ncol=2)
#dev.off()

## try using stat_density and different kernel than gaussian 
###################################################
#dev.off()
p1 <- ggplot() + 
         geom_density(data=ABC_params[[1]],aes(x=theta), adjust=5, kernel='triangular') +
	  		 geom_density(data=post_samps[[1]],aes(x=theta), adjust=5, col="gray", kernel='triangular') +
				 xlim(8,12) + labs(x=expression(theta), y='Density') + theme_classic() +
				 theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) 
p2 <- ggplot() +  
          #geom_density(data=abc,aes(x=theta), adjust = 3/2, kernel='triangular') +
          geom_density(data=ABC_params[[2]],aes(x=theta), adjust = 3/2, kernel='triangular') +
				  #geom_density(data=post,aes(x=theta), col="gray", kernel='triangular') +
          geom_density(data=post_samps[[2]],aes(x=theta), adjust=3/2, col="gray", kernel='triangular') +        
          xlim(8,12) + labs(x=expression(theta), y=' ') + theme_classic() +
				  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) 

#png(paste0("Cauchy_ABC_fig_n",n, ".png"))
grid.arrange(p1,p2,ncol=2)
#dev.off()







