####################################################################################
### Code for: "Global relationships in tree functional traits" by Maynard et al. 
### Contact Dan Maynard (dan.s.maynard@gmail) if you have any questions.
###
### This script takes the random forest models and imputes the traits
### for all of the missing observations. The output from this code is 
### for all subsequent analysis. 
###
### ****NOTE**** Because this file used the raw TRY trait data which is not open-access
### it is not executable using the data in the Github repo. Please contact Dan or
### the TRY database for access to these data. 
###
##############################################################

rm(list = ls())

library(feather)
library(stringr)
library(doParallel)
library(ranger)
library(geosphere)
library(nlme)
library(ape)
library(Metrics)
library(tidyverse)

# source the functions
source("code/R_functions.R")

# read in trait names
use_traits <- read_csv("data/raw_data/Trait_names.csv") %>% 
	mutate(low = ifelse(is.na(low), -Inf, log(low)), high = ifelse(is.na(high), Inf, log(high)))

# extract quantile trait names
quant_traits <- use_traits %>% filter(use_max == 1) %>% select(trait_label) %>% unlist()

# read in the model names
files_both <- str_sort(paste0(getwd(), "/data/model_fits/both/", list.files("data/model_fits/both/")), numeric = TRUE)
files_phy <- str_sort(paste0(getwd(), "/data/model_fits/phy/", list.files("data/model_fits/phy/")), numeric = TRUE)

# specify single fit models
files_both_single <- files_both[!grepl("DOUBLE", files_both)]
files_phy_single <- files_phy[!grepl("DOUBLE", files_phy)]

# specify double fit models
files_both_double <- files_both[grepl("DOUBLE", files_both)]
files_phy_double <- files_phy[grepl("DOUBLE", files_phy)]

# get the variables used in the model
mod1 <- readRDS(files_both[1])
used_vars <- names(mod1[[1]]$model$variable.importance)

# read in trait data. Note, values are logged
tr_obso <- read_feather("data/raw_data/TRY_trait_data.feather")

# unlog and extract variables
tr_obs <- tr_obso %>% 
	select(y, LAT, LON, accepted_bin, TraitID) %>% distinct() %>% 
	mutate(trait = paste0("trait",TraitID)) %>% rename(value = y) %>% select(-TraitID)

# read in scaled phylogeny
phy <- read_csv("data/raw_data/Phylogenetic_eigenvectors.csv",  col_types = cols()) %>% 
	select(accepted_bin, names(.)[names(.)%in%used_vars]) %>% distinct() %>% filter(complete.cases(.))

# read in the composite
env <- read_csv("data/raw_data/Environmental_covariates.csv") %>% 
	select(LAT, LON, names(.)[names(.)%in%used_vars]) %>% distinct() %>% filter(complete.cases(.))


# get the observed trait values
tr_both_obs <- tr_obs %>% filter(!is.na(LAT)) 

tr_phy_obs <- tr_obs %>% group_by(accepted_bin, trait) %>%
	summarize(mean_value = mean(value), max_value = max(value), .groups = "drop_last") %>% ungroup %>%
	rowwise() %>% mutate(value = ifelse(trait%in%quant_traits, max_value, mean_value)) %>% ungroup %>% select(accepted_bin, trait, value)

# add in phylo and env data to traits
tr <- tr_obs %>% filter(!is.na(LAT)) %>% select(LAT, LON, accepted_bin) %>% distinct() %>% 
	left_join(phy, by = "accepted_bin") %>% left_join(env, by = c("LAT", "LON")) %>% filter(complete.cases(.))

# register cores
registerDoParallel(min(detectCores()-2, 72))

# single imputation, phy + env
fit_both_single <- foreach(i = 1:length(files_both_single), .combine = bind_rows, .inorder = FALSE, .multicombine = TRUE, .packages = c("tidyverse", "ranger"))%dopar%{
	outtr <- tr %>% bind_cols(get_prediction(my_ranger_fit = files_both_single[i], my_data = (.), type = "single")) 
	return(outtr)
	
}

# single imputation. phy only
fit_phy_single <- foreach(i = 1:length(files_phy_single), .combine = bind_rows, .inorder = FALSE, .multicombine = TRUE, .packages = c("tidyverse", "ranger"))%dopar%{
	outtr <- phy %>% bind_cols(get_prediction(my_ranger_fit = files_phy_single[i], my_data = (.), type = "single")) 
	return(outtr)
}

# spread the data
res_both_spread <- fit_both_single %>% spread(trait, value)
res_phy_spread <- fit_phy_single %>% spread(trait, value)


# fit the double imputation models. phy+env 
fit_both_double <- foreach(i = 1:length(files_both_double), .combine = bind_rows, .inorder = FALSE, .multicombine = TRUE, .packages = c("tidyverse", "ranger"))%dopar%{
	outtr <- res_both_spread %>% bind_cols(get_prediction(my_ranger_fit = files_both_double[i], my_data = (.), type = "double", qlevel = 0.9)) %>% 
		select(accepted_bin, LAT, LON, trait, value) 
	return(outtr)
}

# fit the double imputation models. phy only
fit_phy_double <- foreach(i = 1:length(files_phy_double), .combine = bind_rows, .inorder = FALSE, .multicombine = TRUE, .packages = c("tidyverse", "ranger"))%dopar%{
	outtr_phy <- res_phy_spread %>% bind_cols(get_prediction(my_ranger_fit = files_phy_double[i], my_data = (.), type = "double", qlevel = 0.9)) %>% 
		select(accepted_bin, trait, value) 
	return(outtr_phy)
}

# add in the lat/long/accepted bin, and the observed values
res_both_double <- fit_both_double %>% left_join(tr_both_obs %>% rename(obs_value = value), by = c("accepted_bin", "LAT", "LON", "trait")) 
res_phy_double <- fit_phy_double %>% left_join(tr_phy_obs %>% rename(obs_value = value), by = c("accepted_bin", "trait"))

# combine
res_comb <- res_both_double %>% mutate(fit = "both") %>% bind_rows(res_phy_double %>% mutate(fit = "phy")) %>% 
	mutate(TraitID = as.numeric(gsub("trait", "", trait))) 

# get the bounds
bound_traits <- use_traits %>% filter(focal_trait == 1) %>% select(high, low, TraitID) 

# replace imputed traits with observed traits, excluding outliers
tr_both <- res_comb %>%
	left_join(bound_traits) %>%
	filter(!is.na(high)) %>%
	rowwise() %>%
	mutate(obs_value = ifelse(!is.na(obs_value), ifelse(obs_value<low | obs_value>high, NA, obs_value), NA)) %>%
	mutate(obs_value = ifelse(trait%in%quant_traits, NA, obs_value)) %>%
	ungroup %>%
	rename(pred_value = value) %>% 
	select(accepted_bin, LAT, LON, fit, TraitID, pred_value, obs_value) %>% 
	arrange(TraitID, accepted_bin, LAT, LON)


# write file, for direct use in PCA code
write_feather(tr_both, "data/results/Estimated_trait_table.feather")





