####################################################################################
### Code for: "Global relationships in tree functional traits" by Maynard et al. 
### Contact Dan Maynard (dan.s.maynard@gmail) if you have any questions.
### 
### This script fits the random forest models with 
### leave-out-out buffered cross validation
###
###
### ****NOTE**** Because this file used the raw TRY trait data which is not open-access
### it is not executable using the data in the Github repo. Please contact Dan or
### the TRY database for access to these data. 
###
####################################################################################

rm(list=ls())

library(feather)
library(doParallel)
library(ranger)
library(tidyverse)

# source the functions
source("code/R_functions.R")

# register parallel cores
registerDoParallel(min(detectCores()-2, 72))

# read in trait names and variables
use_traits <- read_csv("data/raw_data/Trait_names.csv") %>% arrange(TraitID)

# read in the plot-level trait data. Note that values are already logged and aggregated to the pixel/species level
dt <- read_feather("data/raw_data/TRY_trait_data.feather")

# specify the base directory
my_dir_base <- paste0(getwd(), "/data/model_fits/")

# read in angio/gymno
ag <- read_csv("data/raw_data/ANGIO_GYMNO_lookup.csv") %>% 
	mutate(angio = ifelse(group == "Angiosperms", 1, ifelse(group == "Gymnosperms",0, NA))) %>% 
	filter(!is.na(genus))

# make sure the directories exist
dir.create(my_dir_base, showWarnings = FALSE)
dir.create(paste0(my_dir_base,"/phy"), showWarnings = FALSE)
dir.create(paste0(my_dir_base,"/both"), showWarnings = FALSE)
dir.create(paste0(my_dir_base,"/env"), showWarnings = FALSE)

# number of threads for parallel computation
nthreads <- 72

for(ijk in 1:2){
	
	# are we doing double imputation?	
	if(ijk==2){
		
		# create lists for previous fits
		fit_prev_both  <- fit_prev_phy <- fit_prev_env <- vector(mode = "list", length = nrow(use_traits))
		
		# load the previous fits
		for(i in 1:nrow(use_traits)){
			mymod_b <- readRDS(paste0(paste0(my_dir_base,"both/"),use_traits$trait_label[i],file_suffix,".RDS"))
			fit_prev_both[[i]] <- mymod_b[[1]]$model
			names(fit_prev_both)[i] <- names(mymod_b)
			
			mymod_p <- readRDS(paste0(paste0(my_dir_base,"phy/"),use_traits$trait_label[i],file_suffix,".RDS"))
			fit_prev_phy[[i]] <- mymod_p[[1]]$model
			names(fit_prev_phy)[i] <- names(mymod_p)
			
			mymod_e <- readRDS(paste0(paste0(my_dir_base,"env/"),use_traits$trait_label[i],file_suffix,".RDS"))
			fit_prev_env[[i]] <- mymod_e[[1]]$model
			names(fit_prev_env)[i] <- names(mymod_e)
		}
		# edit the suffic for double imputation
		file_suffix <- "_DOUBLE"
	}else{file_suffix <- ""}
	
	# cycle through the traits and fit the random forest models
	for(j in 1:nrow(use_traits)){
		
		# print the current trait
		cat("----------------------------------------------\n")
		cat(j,"of",nrow(use_traits),"--",use_traits$trait_label[j],"--",use_traits$trait[j],"\n")
		
		# are we doing quantile or mean random forest
		summ_func  <-  ifelse(use_traits$use_max[j] == 1, max, mean)
		
		# subset the data to the current trait
		my_X <- dt %>% filter(TraitID == use_traits$TraitID[j]) %>% select(-TraitID)

		# get the final data
		use_X <- my_X %>% 
			select(LAT, LON, accepted_bin, y, tax_bin, env_clust, contains("Phy_"), contains("Env_")) %>% 
			group_by(LAT, LON, accepted_bin) %>% summarize_all(.funs = summ_func, na.rm=T) %>% ungroup
		
		# the subset for phylogeny+environmental coviariates
		df_both <- use_X %>% left_join(ag, by = "accepted_bin") %>% filter(complete.cases(.)) 
		
		# for non-georeferenced phylogeny
		df_phy <- use_X %>% select(accepted_bin, y, tax_bin, contains("Phy_")) %>% 
			group_by(accepted_bin) %>% summarize_all(.funs = summ_func, na.rm=T) %>% ungroup %>%
			left_join(ag, by = "accepted_bin") %>% filter(complete.cases(.))
		
		# georeferenced phylogeny only. not used
		df_env <- use_X %>% left_join(ag, by = "accepted_bin") %>% filter(complete.cases(.)) %>% select(-contains("Phy_"))
		
		# just init to full model if no extrra non-georeferenced species
		if(nrow(df_phy) == 0){
			df_phy <- df_env
		}
		
		if(ijk==2){
			# if this is the second time around, impute the traits to use as predictors
			for(i in 1:nrow(use_traits)){
				if(names(fit_prev_both)[i]!=use_traits$trait_label[j]){
					df_both <- df_both %>% bind_cols(tibble(var = predict(fit_prev_both[[i]], data = df_both)$predictions) %>% setNames(paste0("BOTH_",names(fit_prev_both)[i])))
				}
				if(names(fit_prev_env)[i]!=use_traits$trait_label[j]){
					df_env <- df_env %>% bind_cols(tibble(var = predict(fit_prev_env[[i]], data = df_env)$predictions) %>% setNames(paste0("ENV_",names(fit_prev_env)[i])))
				}
				if(names(fit_prev_phy)[i]!=use_traits$trait_label[j]){
					df_phy <- df_phy %>% bind_cols(tibble(var = predict(fit_prev_phy[[i]], data = df_phy)$predictions) %>% setNames(paste0("PHY_",names(fit_prev_phy)[i])))
				}
				
			}
		}
		
		# are we fitting kfold or just full model?
		do_kfold_each <- I(ijk==2 & use_traits$focal_trait[j]==1)
		
		# ### fit the random forest spatial k-fold. call the wrapper function to make sure changes are reloaded
		res_all <- ranger_kfold_clust(df_both = df_both, df_phy = df_phy, df_env = df_env, dm = dm, do_kfold = do_kfold_each, nfolds = 1000,
									  summ_func = summ_func, seed = 10, crown = I(use_traits$TraitID[j]%in%c(324, 773, 281)),  
									  num_threads = nthreads, quantreg = I(use_traits$use_max[j] == 1))
		
		# save the separate models
		res_env <- res_all$env
		res_phy <- res_all$phy
		res_both <- res_all$both
		
		# save the model results on the 2nd pass through
		if(do_kfold_each & use_traits$focal_trait[j] == 1){
			
			r2_df <- tibble(
				TraitID = use_traits$TraitID[j],
				trait = use_traits$trait[j],
				use_max = use_traits$use_max[j],
				r2_both = round(res_both$r2,2),
				r2_phy =  round(res_phy$r2,2),
				rel_exp_both = round(res_both$obs_pred %>% mutate(rel_err = abs((exp(y) - exp(value))/exp(y))) %>% summarize(rel_err = median(rel_err)) %>% select(rel_err) %>% unlist(),2) %>% unlist(),
				rel_exp_phy = round(res_phy$obs_pred %>% mutate(rel_err = abs((exp(y) - exp(value))/exp(y))) %>% summarize(rel_err = median(rel_err)) %>% select(rel_err) %>% unlist(),2) %>% unlist())
			
			# append to file if it exists
			write_csv(r2_df, paste0(my_dir_base,"RF_summary_results.csv"), append = I(j>1))
			
		}
		# save the models
		saveRDS(list(res_both) %>% setNames(use_traits$trait_label[j]), paste0(paste0(my_dir_base,"both/"),use_traits$trait_label[j],file_suffix,".RDS"))
		saveRDS(list(res_phy) %>% setNames(use_traits$trait_label[j]), paste0(paste0(my_dir_base,"phy/"),use_traits$trait_label[j],file_suffix,".RDS"))
		saveRDS(list(res_env) %>% setNames(use_traits$trait_label[j]), paste0(paste0(my_dir_base,"env/"),use_traits$trait_label[j],file_suffix,".RDS"))
	}
}	

