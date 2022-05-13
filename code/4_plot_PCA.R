####################################################################################
### Code for: "Global relationships in tree functional traits" by Maynard et al. 
### Contact Dan Maynard (dan.s.maynard@gmail) if you have any questions.
###
### This script creates the plots in Figure 2. 
### The output from this code is used for Figs. 3 and 4
###
###
### ***NOTE*** The data being read in this file (Imputed_trait_data.csv) does
### not match the output of 2_Estimate_traits.R (Estimated_trait_table.feather) as
### that file contains raw TRY data that is not open access. The file being loaded 
### in here has already been processed and cleaned to remove this data. Please 
### contact Dan or TRY for the raw data.
###
##############################################################

rm(list = ls())

library(Redmonder)
library(ggrepel)
library(feather)
library(aroma.light)
library(wCorr)
library(colorspace)
library(cowplot)
library(tidyverse)

# source the functions
source("code/R_functions.R")


# read in the list of angio vs gymno species
ag0 <- read_csv("data/raw_data/ANGIO_GYMNO_lookup.csv")

ag <- ag0 %>% select(accepted_bin, group) %>%
	mutate(angio = ifelse(group == "Angiosperms", 1, ifelse(group == "Gymnosperms",0, NA))) %>%
	select(accepted_bin, angio)

# read in the trait summaries
use_traits <- read_csv("data/raw_data/Trait_names.csv") %>% 
	filter(focal_trait == 1) %>% 
	mutate(low = log(low), high = log(high))


# Load in the imputed traits. This data has been cleaned to remove raw data.
tr_both <- read_csv("data/results/Imputed_trait_data.csv")


# add in the species weights
tr_both_wt <- tr_both %>% 
	gather(trait, value, -accepted_bin, -LAT, -LON) %>% 
	group_by(trait) %>% 
	mutate(value = (value - mean(value))/sd(value)) %>% 
	ungroup %>% 
	spread(trait, value) %>% 
	left_join((.) %>% group_by(accepted_bin) %>% tally() %>% ungroup %>% mutate(spp_n = n) %>% select(-n), by = "accepted_bin") %>%
	group_by(accepted_bin) %>% mutate(spp_wt = (1/spp_n)/sum(1/spp_n)) %>% ungroup %>% 
	mutate(wt = spp_wt) %>% 
	mutate(wt = wt/sum(wt)) %>% 
	left_join(ag0 %>% select(accepted_bin, order, family)) %>%
	select(accepted_bin, LAT, LON, all_of(use_traits$trait_short), wt) 


# weighted pca, using all data. See "R_functions.R" for ggpca_weighted function
# note that the random seed in R might flip the axes. adjust flip_y, flip_x to fix. 
gb <- gb0 <- ggpca_weighted(tr_use = tr_both_wt, 
							trait_vars = use_traits$trait_short,
							angio_list = ag, line_wd_big = 1,show_labs = TRUE,
							num_show = 3, alpha = 0, label_size = 3, flip_coord = FALSE,
							show_legend = TRUE, flip_x = 1, flip_z = -1, flip_zz = -1, 
							flip_y = 1, title = NULL)

# just angiosperms
gba <- ggpca_weighted(tr_use = tr_both_wt %>% left_join(ag) %>% filter(angio == 1) %>% select(-angio), 
					  trait_vars = use_traits$trait_short,
					  angio_list = ag, line_wd_big = 0.75, show_labs = TRUE, alpha_const = 1.5, arrow_head = 4,
					  num_show = 3, alpha = 0, label_size = 2, flip_coord = FALSE,
					  show_legend = FALSE, flip_x = 1, use_vars = gb$vars,
					  flip_y = -1, title = NULL)

# just gymnosperms
gbg <- ggpca_weighted(tr_use = tr_both_wt %>% left_join(ag) %>% filter(angio == 0) %>% select(-angio), 
					  trait_vars = use_traits$trait_short,
					  angio_list = ag,line_wd_big = 0.75, show_labs = TRUE, alpha_const = 1.25, arrow_head = 4, 
					  num_show = 3, alpha = 0, label_size = 2, flip_coord = FALSE, use_vars = gb$vars,
					  show_legend = FALSE, flip_x = -1, 
					  flip_y = -1, title = NULL)


# save the results for use in other figures
write_feather(tr_both_wt, paste0("data/results/PCA_data_for_dendro_BOTH.feather"))
write_feather(gb$tr, paste0("data/results/PCA_df_for_corr_BOTH.feather"))
write_csv(tibble(value = round(gb$axis_perc*100, 1)), paste0("data/results/PCA_axis_perc_BOTH.csv"))
write_feather(gb$loadings %>% gather(pc, value, -trait) %>% mutate(value = round(value,2)) %>% spread(pc, value), paste0("data/results/PCA_loadings_BOTH.feather"))


# plot
gcomb <- plot_grid(gba$plot, gbg$plot, nrow = 2)
plot_grid(gb$plot, gcomb, nrow = 1, rel_widths = c(1, 0.6), rel_heights = c(1, 0.8))


