median((raw_results_all$setting_1$wid_vol_all$rejACC)[1,])/median((raw_results_all$setting_1$wid_vol_all$ISABC)[1,])

#------
# 5/16/23

## Create for examples in Thorton, Li and Xie (2018) 

library(matrixStats) # For function 'rowMads'
library(snowfall) # For parallelisation
library(parallel) # For parallelisation
library(pbapply) # For function 'pblapply'
library(mcmc)

setwd('/Users/wentao/bitbucket/')
# setwd('D:\\bitbucket\\core-functions')
source('Auxiliary_Functions.R')
source('Divide_Recombine_Parallelise.R') 
source('Functions_of_RejABC.R')
source('Functions_of_RegABC.R')
source('Functions_of_ABC_PMC.R')

source('Cauchy_example_functionsv2.r')


##########################
# Experiment Setting
##########################
### Setup for parallelisation ###
platform<-'Mac'

### Experiment variables ###
test_size<-500; N<-5*10^5; 
# test_size<-300; N<-1e5
n_all<-400
# n_all<-50

p_all<-signif(c(0.005,0.05,0.1,0.2,0.4), digits=2)  ## five acceptance proportions 

parameter_true<-c(location_cubicrt=10^(1/3),scale_sqr=0.55^2); CI_alpha<-0.95 # Data model is Cauchy(location_cubicrt^3,scale_sqr^1/2)

 prior_scale_all<-c(0.5,1.5,1.5); 
 prior_location_all<-parameter_true[2]+c(0,0,prior_scale_all[3]); prior_location_ind<-c('True','True','False')
 number_prior<-length(prior_scale_all)+1

set.seed(100)
tested_observations_all<-simulate_pseudo_datasets(matrix(rep(parameter_true,test_size),ncol=2,byrow=T),n_all[1]) # a test_size*n matrix

save(tested_observations_all, file=paste0("tested_obs_all_n", n_all, ".rdata"))

##########################

res <- Example1_main(test_size,parameter_true,tested_observations_all,
              parameter_setting=1,
              posterior_sample=NULL,N,platform,n_all,p_all,
              CI_alpha=0.05,divide=FALSE,divide_paral=FALSE,
              prior_choice='Cauchy',  #'jeffrey'
              prior_location=parameter_true[1],
              prior_scale=0.5)

##########################
# Experiment 
##########################
# Parameter settings include: 1 -- location unknown and scale known with median as summary, 2 -- location unknown and scale known with mean as summary, 3 -- location known and scale unknown with median absolute deviation (MAD) as summary, 4 -- both unknown.
# Parameter proposals include: 1 -- KDE of certain point estimates, e.g. median and MAD, 2 -- Cauchy distribution with given location and scale.
# Choices of prior include: Student t distribution with degree of freedom 4 and scale parameter, 'prior_scale' 1:4; Jeffrey prior for location and scale family, i.e. uniform over all location and uniform over all log-scale.

###########################################################################
# Cauchy example: N1e5, test_size 500 


raw_results_all<-list(list(0,0),list(0),list(0),list(0),list(0))
dev.new()
for(i in 1:4){
#  if(i!=2) prior_set<-number_prior
#  if(i==2) prior_set<-1:number_prior
  if(i==1) {parameter_setting<-1; parameter_proposal<-1; prior_location<-parameter_true[1]}
  if(i==2) {parameter_setting<-2; parameter_proposal<-2; prior_location<-parameter_true[1]}
  if(i==3) {parameter_setting<-3; parameter_proposal<-1; prior_location<-parameter_true[2]}
  if(i==4) {parameter_setting<-4; parameter_proposal<-1; prior_location<-parameter_true}
  if(i==5) {parameter_setting<-5; parameter_proposal<-2; prior_location<-parameter_true}
  
  
  names(raw_results_all)[parameter_setting]<-paste0('setting_',parameter_setting)
  names(raw_results_all[[parameter_setting]])[parameter_proposal]<-paste0('proposal_',parameter_proposal)
  raw_results_all[[parameter_setting]][[parameter_proposal]]<-list(0)
  
  for(prior_i in prior_set){
    ##########################
    # Set Random Seed 
    ##########################
    set.seed(100)
    
    if(prior_i<number_prior){
      prior_scale<-prior_scale_all[prior_i]
      if(i==4) prior_scale<-c(prior_scale_all[prior_i],prior_scale_all[prior_i])
      prior_location<-prior_location_all[prior_i]
      tmp<-Example1_main(test_size,parameter_true,tested_observations_all=tested_observations_all,parameter_setting,parameter_proposal,N=N,platform=platform,n_all=n_all,p_all=p_all,CI_alpha=CI_alpha,divide=TRUE,PMC_refine=FALSE,prior_choice='t4',prior_location=prior_location,prior_scale=prior_scale)
      raw_results_all[[parameter_setting]][[parameter_proposal]][[prior_i]]<-tmp
      names(raw_results_all[[parameter_setting]][[parameter_proposal]])[prior_i]<-paste0('prior_scale',prior_scale_all[prior_i],'_location',prior_location_ind[prior_i])
    }
    if(prior_i==number_prior){ 
      tmp<-Example1_main(test_size,parameter_true,tested_observations_all=tested_observations_all,parameter_setting,parameter_proposal,N=N,platform=platform,n_all=n_all,p_all=p_all,CI_alpha=CI_alpha,divide=TRUE,PMC_refine=FALSE)
      raw_results_all[[parameter_setting]][[parameter_proposal]][[prior_i]]<-tmp
      names(raw_results_all[[parameter_setting]][[parameter_proposal]])[prior_i]<-paste0('uniform')
    }
  }
}

rm(sample_posterior)
save.image(file=paste0('Cauchy_results_',test_size,'rep_N',N,'_n',n_all,'_alpha',CI_alpha*100,'.RData'))

#############################################
# Cauchy example: N1e5, test_size 300 
# parameter_setting 1:4,6
# prior t4; prior_scale_used<-3
#############################################

##########################
# Results
##########################
# Under each setting, the raw result is a list with four elements: the posterior mean estimates, posterior variance estimates, coverage rate of the true parameter value in the test_size runs, and the tolerance values for all runs. Each element is a list with two levels. Level one is for the data sizes n_all, level two is for the acceptance proportions p_all, and the basic element is a test_size*2 matrix, for moments estimate, or a length test_size vector.

number_methods<-6 
results_all<-list(list(0,0),list(0),list(0),list(0),list(0)) # Corresponds to five experiment settings
tmp_matrix<-matrix(0,number_methods,length(p_all))
colnames(tmp_matrix)<-paste0('p=',p_all)
rownames(tmp_matrix)<-c('ISABC','ISregABC','ISpi_ACC','ISpi_regACC','rejACC','regACC')
for(i in c(1:4,6)){
  twoD_coverage<-FALSE
  if(i!=3) prior_set<-number_prior
  if(i==3) prior_set<-1:number_prior
  if(i==1) {parameter_setting<-1; parameter_proposal<-1; param_est<-1; number_coverage_type<-1; type_name<-'location'}
  if(i==2) {parameter_setting<-1; parameter_proposal<-2; param_est<-1; number_coverage_type<-1; type_name<-'location'}
  if(i==3) {parameter_setting<-2; parameter_proposal<-1; param_est<-1; number_coverage_type<-1; type_name<-'location'}
  if(i==4) {parameter_setting<-3; parameter_proposal<-1; param_est<-2; number_coverage_type<-1; type_name<-'scale'}
  if(i==5) {parameter_setting<-4; parameter_proposal<-1; param_est<-2; number_coverage_type<-1; type_name<-'scale'}
  if(i==6) {parameter_setting<-5; parameter_proposal<-1; param_est<-1:2; number_coverage_type<-3; type_name<-c('location','scale','both')}
  names(results_all)[parameter_setting]<-paste0('setting_',parameter_setting)
  names(results_all[[parameter_setting]])[parameter_proposal]<-paste0('proposal_',parameter_proposal)
  
  tmp<-list(0)
  for(n_i in 1:length(n_all)){
    tmp[[n_i]]<-list(0); param_count<-1
    names(tmp)[n_i]<-paste0('n=',n_all[n_i])
    for(param_i in param_est){
      tmp[[n_i]][[param_count]]<-list(tmp_matrix,tmp_matrix,tmp_matrix,tmp_matrix)
      if(param_i==1) names(tmp[[n_i]])[param_count]<-c('location')
      if(param_i==2) names(tmp[[n_i]])[param_count]<-c('scale')
      names(tmp[[n_i]][[param_count]])<-c('mean_avg','mean_var','var_avg','var_var')
      param_count<-param_count+1
    }
    param_count<-param_count-1
    for(coverge_i in 1:number_coverage_type){
      tmp[[n_i]][[param_count+coverge_i]]<-tmp_matrix
      names(tmp[[n_i]])[param_count+coverge_i]<-paste0(type_name[coverge_i],'_coverage')
    }
    for(coverge_i in 1:number_coverage_type){
      tmp[[n_i]][[param_count+number_coverage_type+coverge_i]]<-tmp_matrix
      names(tmp[[n_i]])[param_count+number_coverage_type+coverge_i]<-paste0(type_name[coverge_i],'_CIwidth')
    }
  }
  results_all[[parameter_setting]][[parameter_proposal]]<-list(0)
  for(prior_i in prior_set){
    tmp_raw<-raw_results_all[[parameter_setting]][[parameter_proposal]][[prior_i]]
    output_all<-list(n=n_all,p=p_all,post_means=tmp_raw$post_means_all,post_variances=tmp_raw$post_variances_all,post_coverage=(tmp_raw$post_coverage_all)$twoside,post_coverage_twoD=(tmp_raw$post_coverage_all)$twoD,CIbounds=tmp_raw$CIbounds_all,results=tmp)
    results_all[[parameter_setting]][[parameter_proposal]][[prior_i]]<-Process_all_output(1:test_size,output_all,param_est=param_est,number_coverage_type=number_coverage_type)
    names(results_all[[parameter_setting]][[parameter_proposal]])[prior_i]<-names(raw_results_all[[parameter_setting]][[parameter_proposal]])[prior_i]
  }
}


#------- is any of the below needed? 
#? load('Cauchy_true_posterior_sample_300rep_n400.RData')
#? p_i_table<-c(1,5,17); method_adj<-c(2,4,6); method_credible<-2
#? tmp_matrix<-matrix(0,nrow=6,ncol=8)
#? rownames(tmp_matrix)<-paste0('p=',rep(p_all[p_i_table],rep(2,3)),' ',rep(c('confidence','credible'),3))
#? colnames(tmp_matrix)<-c('rej-ABC','width','IS-ABC','width','piACC','width','rACC','width')
tmp_matrix[c(2,4,6),5:8]<-NaN; 
results_table<-list(0)
for(i in c(1:4,6)){
  if(i!=3) prior_set<-number_prior
  if(i==3) prior_set<-1:number_prior
  if(i==1) {parameter_setting<-1; parameter_proposal<-1; param_est<-1; number_coverage_type<-1;setting_name<-'location/KDE/median'; sub_name<-'(i)'; type_name<-'location'}
  if(i==2) {parameter_setting<-1; parameter_proposal<-2; param_est<-1; number_coverage_type<-1;setting_name<-'location/Cauchy/median'; sub_name<-'(ii)'; type_name<-'location'}
  if(i==3) {parameter_setting<-2; parameter_proposal<-1; param_est<-1; number_coverage_type<-1;setting_name<-'location/KDE/mean'; sub_name<-'(iii)'; type_name<-'location'}
  if(i==4) {parameter_setting<-3; parameter_proposal<-1; param_est<-2; number_coverage_type<-1;setting_name<-'scale/KDE/MAD'; sub_name<-'(iv)'; type_name<-'scale'}
  if(i==5) {parameter_setting<-4; parameter_proposal<-1; param_est<-2; number_coverage_type<-1;setting_name<-'scale/KDE/HL'; sub_name<-'(v)'; type_name<-'scale'}
  if(i==6) {parameter_setting<-5; parameter_proposal<-1; param_est<-1:2; number_coverage_type<-3;setting_name<-'both/KDE'; sub_name<-'(vi)'; type_name<-c('location','scale','both')}
  
  prior_count<-1; results_table[[i]]<-list(0)
  names(results_table)[i]<-paste0(setting_name,' ',sub_name)
  for(coverge_i in 1:number_coverage_type){
    results_table[[i]][[coverge_i]]<-list(0)
    names(results_table[[i]])[coverge_i]<-paste0(type_name[coverge_i],'_coverage')
    for(prior_i in prior_set){
      results_tmp<-results_all[[parameter_setting]][[parameter_proposal]][[prior_i]][[1]]
      results_table[[i]][[coverge_i]][[prior_count]]<-tmp_matrix
      results_table[[i]][[coverge_i]][[prior_count]][c(1,3,5),c(3,5,7)]<-signif(t(results_tmp[[length(param_est)+coverge_i]][method_adj,p_i_table]),3)
      results_table[[i]][[coverge_i]][[prior_count]][c(1,3,5),c(4,6,8)]<-signif(t(results_tmp[[length(param_est)+number_coverage_type+coverge_i]][method_adj,p_i_table]),6)
      if(length(param_est)==1){
        tmp_credible<-matrix(0,length(p_i_table),length(method_credible))
        CIbounds_all<-raw_results_all[[parameter_setting]][[parameter_proposal]][[prior_i]]$CIbounds_all
        CIupper_all<-CIbounds_all$twosideupper
        CIlower_all<-CIbounds_all$twosidelower
        for(credible_i in 1:length(method_credible)){
          method_i<-method_credible[credible_i]
          for(table_i in 1:length(p_i_table)){ 
            p_i<-p_i_table[table_i]
            tmp<-pblapply(as.list(1:test_size),function(ind_i)sum((sample_posterior[[ind_i]][,param_est]<=CIupper_all[[param_est]][[method_i]][p_i,ind_i])&(sample_posterior[[ind_i]][,param_est]>=CIlower_all[[param_est]][[method_i]][p_i,ind_i]))/nrow(sample_posterior[[ind_i]]))
            tmp_credible[table_i,credible_i]<-median(unlist(tmp))
          }
        }
      }
      results_table[[i]][[coverge_i]][[prior_count]][c(2,4,6),3]<-signif(tmp_credible,3)
      names(results_table[[i]][[coverge_i]])[prior_count]<-names(raw_results_all[[parameter_setting]][[parameter_proposal]])[prior_i]
      prior_count<-prior_count+1
    }		
  }
}
results_table

