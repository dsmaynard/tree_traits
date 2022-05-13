####################################################################################
### Code for: "Global relationships in tree functional traits" by Maynard et al. 
### Contact Dan Maynard (dan.s.maynard@gmail) if you have any questions.
###
### This script creates the plots in Figure 3. 
### It relies on the output from the PCA code. 
###
##############################################################

rm(list = ls())

library(colorspace)
library(feather)
library(doParallel)
library(foreach)
library(viridis)
library(latex2exp)
library(cowplot)
library(ranger)
library(fastshap)
library(ggbeeswarm)
library(tidyverse)

# source the functions
source("code/R_functions.R")

# read in angio vs gymno
ag <- read_csv("data/raw_data/ANGIO_GYMNO_lookup.csv") %>% select(accepted_bin, group) %>%
	mutate(angio = ifelse(group == "Angiosperms", 1, ifelse(group == "Gymnosperms",0, NA))) %>% 
	select(accepted_bin, angio)

# read in the environmental covariates
env <- read_csv("data/raw_data/Environmental_covariates.csv")

# read in the environmental names
enames <- read_csv("data/raw_data/env_names.csv") 

# specify the focal covariates 
env_clust <- c(
	"Env_CHELSA_BIO_Precipitation_Seasonality",
	"Env_CHELSA_BIO_Annual_Precipitation",
	"Env_CHELSA_BIO_Annual_Mean_Temperature",
	"Env_SG_SOC_Content_015cm",
	"Env_SG_Sand_Content_015cm",
	"Env_SG_Depth_to_bedrock",
	"Env_EarthEnvTopoMed_Northness",
	"Env_CHELSA_BIO_Mean_Diurnal_Range",
	"Env_EarthEnvTopoMed_Elevation",
	"Env_EsaCci_BurntAreasProbability"
)

# read in the PCA data and scale the PC variables
tr_num <- 
	read_feather("data/results/PCA_df_for_corr_BOTH.feather") %>%
	left_join(env) %>% left_join(ag) %>% 
	mutate(PC1 = scale11(PC1), PC2 = scale11(PC2), PC3 = scale11(PC3), PC4 = scale11(PC4)) 


# filter the environemntal names to those used for fitting
enames <- enames %>% filter(var3%in%names(tr_num)) %>% mutate(show_units = ifelse(is.na(show_units), "", show_units))

# scale variables to correct units and log the depths
tr_num$Env_CHELSA_BIO_Annual_Mean_Temperature <- tr_num$Env_CHELSA_BIO_Annual_Mean_Temperature/10
tr_num$Env_CHELSA_BIO_Mean_Diurnal_Range <- (tr_num$Env_CHELSA_BIO_Mean_Diurnal_Range) + 2731.5
tr_num$Env_SG_Depth_to_bedrock <- log(max(tr_num$Env_SG_Depth_to_bedrock, na.rm=T) - tr_num$Env_SG_Depth_to_bedrock + 1)
tr_num$Env_SG_Depth_to_bedrock <- max(tr_num$Env_SG_Depth_to_bedrock, na.rm=T) - tr_num$Env_SG_Depth_to_bedrock
tr_num$Env_EarthEnvTopoMed_Elevation <- log((tr_num$Env_EarthEnvTopoMed_Elevation-min(tr_num$Env_EarthEnvTopoMed_Elevation)))

# remove any NAs from the dataset. add in row id and get inverse-species weights
treg_use <- tr_num %>% select(LAT, LON, PC1, PC2, accepted_bin, wt, all_of(enames$var3)) %>% left_join(ag) %>%
	filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %>% 
	filter(complete.cases(.)) %>% 
	mutate(id = 1:nrow(.)) %>% group_by(accepted_bin) %>% mutate(wt = length(accepted_bin)) %>% ungroup %>% mutate(wt = 1/wt)  


############################################
#### Shapley values
#############################################


# specify the number of out-of-fit test points, and number of simulations for fastshap. the more sims the better
ntest <- 5000
nsim <- 100

# register parallel cores
registerDoParallel(min(detectCores()-2, 72))

# prediction function for fastshap
pfun <- function(object, newdata) {
	predict(object, data = newdata)$predictions
}

set.seed(15)

# create the test/train split
test_id <- treg_use %>% slice(sample(1:nrow(.), ntest, replace = FALSE)) %>% select(id) %>% unlist() %>% as.numeric()
train_data <- treg_use %>% filter(!id%in%test_id)
test_data <- treg_use %>% filter(id%in%test_id) 

# Fit the random forest models to the training data
r1 <- ranger(PC1~., data = train_data %>% select(PC1, contains("Env_")), case.weights = train_data$wt)
r2 <- ranger(PC2~., data = train_data %>% select(PC2, contains("Env_")), case.weights = train_data$wt)

#  Estimate the shapley values on the test data
fs1 <- fastshap::explain(r1, X = test_data %>% select(contains("Env_")) %>% as.matrix(), pred_wrapper = pfun, nsim = nsim, .parallel = TRUE, adjust = TRUE)
fs2 <- fastshap::explain(r2, X = test_data %>% select(contains("Env_")) %>% as.matrix(), pred_wrapper = pfun, nsim = nsim, .parallel = TRUE, adjust = TRUE)

# combine the models and merge in the environmental covariates
shap_df <- bind_rows(fs1 %>% as_tibble() %>% mutate(id = test_data$id) %>% gather(var, shap, -id) %>% mutate(type = "shap", axis = 1),
					 fs2 %>% as_tibble() %>% mutate(id = test_data$id) %>% gather(var, shap, -id) %>% mutate(type = "shap", axis = 2)) %>%
	left_join(test_data %>% select(LAT, LON, accepted_bin, angio, id, contains("Env_")) %>% gather(var, value, -id, -angio, -accepted_bin, -LAT, -LON)) 

############################################
#### Plotting
#############################################

# baseline alpha transparency values
al <- c(0.15, 0.5)

# ylims for beeswarm plot
ylim <- c(-0.25, 0.25)

# blank plots for top 3 predictors of each axis
gg <- vector(mode = "list", length = 2*3)

# blank plots for beeswarm plots
varimp <- vector(mod = "list", length = 2)

# colors and shapes for plotting
col1 <- "#1A85FF"
col2 <- "#D41159"
my_bins <- 15
my_shape <- 23
my_cols <- c("red", "darkblue", "orange", "darkslategray3")

# initialize plotting index
count <- 0

# plot for each PC axis
for(j in 1:2){
	
	# select the axis and clean up some of the variables
	my_df <- shap_df %>% filter(axis == j) %>% 
		filter(var != "angio") %>% 
		filter(var%in%env_clust) %>% 
		left_join(enames %>% select(var3, var_short) %>% rename(var = var3), by = "var") %>% 
		select(-var) %>% rename(var = var_short) %>% 
		mutate(value = ifelse(grepl("Elevation", var), exp(value), value)) %>% 
		mutate(angio = ifelse(angio == 1, "Angiosperm", "Gymnosperm"))
	
	# sort the variables by sum(abs(shap value))
	shap_ord <- my_df %>% group_by(var) %>% summarize(feature_imp = mean(abs(shap))) %>% ungroup %>% arrange(desc(feature_imp)) %>% select(var) %>% unlist() %>% as.character()

	# create the beeswarm plot
	varimp[[j]] <- my_df %>% group_by(var) %>% mutate(value = symmetric_scale(value)) %>% ungroup %>%
		mutate(var = factor(var, levels = rev(shap_ord))) %>% 
		ggplot(aes(x = var, y = shap, col = value))+
		geom_quasirandom(dodge.width = 0, bandwidth = 0.2, cex = 0.03, groupOnX = TRUE, alpha = 0.5)+coord_flip()+
		ylab(paste0("Influence on PC Axis ", j, " (Shapley value)"))+
		theme_bw()+
		theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10))+
		scale_color_viridis(option = "viridis", limits = c(-1,1), breaks = seq(-1, 1, length = 10), labels = c("Low", rep("", 10-2), "High"))+
		labs(color = "Feature value")+ylim(ylim)

	# get the three strongest variables for plotting the full relationship
	var_short <- shap_ord[1:3]
	
	# cycle through the top 3 and creat the subplot
	for(i in 1:3){
		# keep track of plotting index
		count <- count+1
		
		# subset the data to the current variable
		sub_df <- my_df %>% filter(var == var_short[i])

		# get the x range and create the x sequence
		x_range <- sub_df %>% select(value) %>% distinct() %>% summarize(max = max(value), min = min(value)) %>% unlist %>% as.numeric %>% sort()
		
		# add in the units for the x axis, format for latex
		xnm <- enames %>% rename(vs = var_short) %>% filter(all_of(var_short[i]) == vs) %>% mutate(var = ifelse(!is.na(show_units), paste(vs, show_units), vs)) %>% select(var) %>% unlist() %>% as.character()
		
		# get the jitter offset based on 50 bins
		offset <- (x_range[2] - x_range[1])/50
		if(grepl("Elevation", var_short[i], ignore.case = TRUE)){
			offset <- (log(x_range[2]) - log(x_range[1]))/50
		}
		
		# get symmetrical ylimits
		ylims <- sub_df %>% filter(value>x_range[1], value<x_range[2]) %>% summarize(max(abs(max(shap)), abs(min(shap)))) %>% unlist()
		
		# get the proportion of angio/gymno, for scaling colors
		my_props <- my_df %>% group_by(angio) %>% tally() %>% mutate(n = n/sum(n)) %>% arrange(angio) %>% select(n) %>% unlist()
		
		# scale the alpha transparencies by prop of angio/gymno
		aluse <- al
		aluse[1] <- al[2]/(my_props[1]/my_props[2])
		
		# create the plot
		ggtmp <- 
			ggplot(data = sub_df %>% filter(value>x_range[1], value<x_range[2]), aes(x = value, y = shap))+
			geom_jitter(width = offset, size = 0.75, aes(color = factor(angio), shape = factor(angio), fill = factor(angio), alpha = factor(angio))) +
			geom_hline(yintercept = 0, linetype = 2, color = "gray50")+
			theme_bw()+
			scale_y_continuous(expand = c(0, 0), limits = c(-ylims, ylims), name = paste0("Influence on PC Axis ",j," (Shapley value)")) + 
			scale_shape_manual(name="Legend", values = c(22, 24))+
			scale_color_manual(name="Legend", values = darken(c(col1, col2), 0.3))+
			scale_fill_manual(name="Legend", values = lighten(c(col1, col2), 0))+
			scale_alpha_manual(values = aluse)
		
		# if looking at elevation, use the log scale
		if(grepl("Elevation", var_short[i], ignore.case = TRUE)){
			ggtmp <- ggtmp + scale_x_log10(expand = c(0, 0), name = TeX(xnm)) 
		}else{
			ggtmp <- ggtmp + scale_x_continuous(expand = c(0, 0), name = TeX(xnm))  
		}
		
		# remove the legend
		gg[[count]] <- ggtmp + theme(legend.position = "none",
									 plot.margin = unit(c(3,20,3,3), "pt"))
		
	}
}


# creat the subplots
p1 <- plot_grid(plotlist = varimp, ncol = 1, byrow = FALSE) 
p2 <- plot_grid(plotlist = gg, nrow = 2, byrow = TRUE)

# plot Fig 3
plot_grid(p1, NULL, p2, ncol = 3, rel_widths = c(0.65, 0.05, 1))

