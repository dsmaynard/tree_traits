####################################################################################
### Code for: "Global relationships in tree functional traits" by Maynard et al. 
### Contact Dan Maynard (dan.s.maynard@gmail) if you have any questions.
###
### This script calculates phylogenetic conservatism and 
### phylogenetic-independent contrasts. It relies on the 
### output from the PCA plotting code, and its output is
### required to run the dendrogram code. 
###
### ****NOTE**** Because this file used the raw TRY trait data which is not open-access,
### the last section (calculating phylogenetic signal) is not executable using the data
### in the Github repo. Please contact Dan or the TRY database for access to these data. 
###
##############################################################


rm(list=ls())

library(feather)
library(phylosignal)
library(phytools)
library(tidyverse)

# source the functions
source("code/R_functions.R")

# read in the PCA data
tr_num <- read_feather("data/results/PCA_df_for_corr_BOTH.feather")

# read in the trait names
tn <- read_csv("data/raw_data/Trait_names.csv") %>% filter(focal_trait == 1)

# read in the cleaned dichotomous species tree 
phm <- read.tree("data/raw_data/dichotomous_phylogenetic_tree.tree")
phm$tip.label <- gsub("_", " ", phm$tip.label)

# get species level averages
tr_mean <- tr_num %>% select(accepted_bin, PC1, PC2, PC3, PC4, all_of(tn$trait_short)) %>% group_by(accepted_bin) %>% summarize_all(.funs = mean) %>% ungroup

###########################
### Independent contrasts

picdf <- tibble()

# calculate phylogenetic contrast for each trait
for(i in 1:ncol(tr_mean)){
	if(names(tr_mean)[i]!="accepted_bin"){
		tri <- tr_mean[,i] %>% unlist
		names(tri) <- tr_mean$accepted_bin
		if(nrow(picdf)>0){
			picdf <- picdf %>% bind_cols(tibble(x = pic(tri, phy = phm)) %>% setNames(names(tr_mean)[i]))
		}else{
			picdf <- tibble(x = pic(tri, phy = phm)) %>% setNames(names(tr_mean)[i])
		}
	}
}

# calculate correlation between contrasts
cm <- matrix(0, ncol(picdf), ncol(picdf))
for(i in 1:ncol(picdf)){
	for(j in 1:ncol(picdf)){
		cm[i,j] <- cm[j,i] <- cor(picdf[,i] %>% unlist, picdf[,j] %>% unlist, method = "spearman")
	}
}

write_csv(cm %>% data.frame() %>% setNames(names(picdf)) %>% as_tibble(), paste0("data/results/PIC_contrasts.csv"))


###################################
### phylogenetic signal

# read in the raw trait data, get the species-level average, and exponentiate
trb <- read_feather("data/raw_data/TRY_trait_data.feather") %>% 
	select(accepted_bin, TraitID, y) %>% group_by(accepted_bin, TraitID) %>% summarize(y = mean(y)) %>% ungroup %>% 
	mutate(trait = paste0("trait", TraitID)) %>% select(-TraitID) %>% #mutate(y = exp(y)) %>% 
	filter(trait%in%tn$trait_label)

# unique traits
utr <- str_sort(unique(trb$trait), numeric = TRUE)

res <- tibble()
set.seed(10)

# calculate phylogenetic signal for each trait
for(i in 1:length(utr)){
	
	print(i)
	
	# get the focal trait name
	my_trait <- utr[i]
	
	# subset the traits
	my_tr <- trb %>% filter(trait == my_trait)
	
	# subset the tree
	my_ph <- drop.tip(phm, trim.internal = TRUE, tip = phm$tip.label[!phm$tip.label%in%my_tr$accepted_bin])
	
	# filter to only species in tree
	my_tr <- my_tr %>% filter(accepted_bin%in%my_ph$tip.label)
	
	# get a named vector
	my_vec <- my_tr$y %>% unlist()
	names(my_vec) <- my_tr$accepted_bin
	
	# calculate signal
	my_res <- phylosig(my_ph, my_vec, method = "lambda", test = TRUE, nsim = 500)
	
	# add to the output data
	res <- res %>% bind_rows(tibble(trait = my_trait, value = my_res$lambda, p = round(my_res$P, 8), method = "lambda"))
}
	
write_csv(res, "data/results/phylo_conservatism.csv")
