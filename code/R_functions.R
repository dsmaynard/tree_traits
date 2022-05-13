####################################################################################
### Code for: "Global relationships in tree functional traits" by Maynard et al. 
### For Review purposes only. 
### 
### This script contains additional functions called by the main scripts
###
###########################################################################

# calculate coefficient of determination relative to the 1:1 line, also called "VEcv", Li 2017
coef_det <- function(xtrue, xpred){
	return(1-sum((xtrue-xpred)^2)/sum((xtrue-mean(xtrue))^2))
}

# scale a variable from -1 to 1
scale11<-function(x){
	if(sd(x) == 0){
		return(x)
	}else{
		mym <- max(c(max(abs(x[x>0])), max(abs(x[x<0]))))
		x <- x/mym
		return(x)
	}
}

# scale function to return a vector, rather than built in scale()
my_scale <- function(x){
	return((x-mean(x, na.rm=T))/sd(x,na.rm=T))
}



# symmetrically scale the values to -1, 1, centered on the median, removing outliers for plotting
symmetric_scale <- function(x, trim = 0.025){
	x <- x  - median(x)
	x[x<0] <- x[x<0]/abs(quantile(x, trim))
	x[x>0] <- x[x>0]/abs(quantile(x, 1-trim))
	x[x>1] <- 1
	x[x< (-1)] <- -1
	return(x)
}


# main code for fitting the random forest trait prediction models and implementing
# spatially + phylogenetically buffered leave one out cross validation
# Note that this also fits georeferenced phylogeny-only predictions, though not used in this analysis
ranger_kfold_clust <- function(df_both, df_phy, df_env, nfolds, num_threads, quantreg, family_genus, crown, do_kfold = FALSE, ...){
	
	set.seed(10)
	
	
	# get all possible trait variables
	en_all <- df_both %>% select(contains("ENV_trait", ignore.case = FALSE), contains("Env_", ignore.case = FALSE)) %>% names()
	pn_all <- df_both %>% select(contains("Phy_", ignore.case = FALSE)) %>% names()
	bn_all <- df_both %>% select(contains("BOTH_trait", ignore.case = FALSE), contains("Env_", ignore.case = FALSE), contains("Phy_", ignore.case = FALSE)) %>% names()
	
	
	set.seed(10)
	# get the unique species and reduce folds if needed
	spp_fold <- sample(unique(df_both$accepted_bin), nfolds, replace = TRUE)
	spp_fold_phy <- sample(unique(df_phy$accepted_bin), nfolds, replace = TRUE)
	nfolds <- min(nfolds, length(spp_fold))
	
	# can run kfold or just fit the full model
	if(do_kfold){
		
		cat("Running k-fold cross validation\n")
		
		obs_pred <- r2e <- r2p <- r2b <- my_r2 <- tibble()
		
		for(ii in 1:nfolds){
			
			# pick a specific species cluster within a spatial cluster
			test_both <- df_both %>% filter(accepted_bin==spp_fold[ii])
			test_env <- df_env %>% filter(accepted_bin==spp_fold[ii])
			test_phy <- df_phy %>% filter(accepted_bin==spp_fold_phy[ii])
			
			if(crown){
				# get all the points either far away phylogenetically or spatially
				train_both <- df_both %>% filter(accepted_bin!=spp_fold[ii], !env_clust%in%test_both$env_clust | !family%in%test_both$family)
				train_env <- df_env %>% filter(accepted_bin!=spp_fold[ii], !env_clust%in%test_both$env_clust | !family%in%test_both$family) #df_env %>% filter(accepted_bin!=spp_fold[ii]) #, !family%in%test_env$family) 
				train_phy <- df_phy %>% filter(accepted_bin!=spp_fold_phy[ii]) 
			}else{
				train_both <- df_both %>% filter(accepted_bin!=spp_fold[ii], !env_clust%in%test_both$env_clust | !genus%in%test_both$genus)
				train_env <- df_env %>% filter(accepted_bin!=spp_fold[ii], !env_clust%in%test_both$env_clust | !genus%in%test_both$genus) # df_env %>% filter(accepted_bin!=spp_fold[ii]) #, !genus%in%test_env$genus) 
				train_phy <- df_phy %>% filter(accepted_bin!=spp_fold_phy[ii])  
			}
			
			if(nrow(train_both)>10){
				
				# fit the random forest models
				rb <- ranger(formula = y~., data=train_both %>% select(y, all_of(bn_all)), seed = 9, 
							 quantreg = quantreg, num.trees = 50, num.threads = num_threads)
				
				rp <- ranger(formula = y~., data=train_phy %>% select(y, all_of(pn_all)), seed = 9, 
							 quantreg = quantreg, num.trees = 50, num.threads = num_threads)
				
				re <- ranger(formula = y~., data=train_env %>% select(y, all_of(en_all)), seed = 9,
							 quantreg = quantreg, num.trees = 50, num.threads = num_threads)
				
				# get the out-of-fit prediction
				my_pred <- test_both %>% mutate(value = predict(rb, (.))$predictions) %>% select(y, value, accepted_bin) %>% mutate(fit = "pred_b", n = nrow(test_both)) %>% 
					bind_rows(test_phy %>% mutate(value = predict(rp, (.))$predictions) %>% select(y, value, accepted_bin) %>% mutate(fit = "pred_p", n = nrow(test_phy))) %>%
					bind_rows(test_env %>% mutate(value = predict(re, (.))$predictions) %>% select(y, value, accepted_bin) %>% mutate(fit = "pred_e", n = nrow(test_env))) %>% 
					group_by(fit, accepted_bin) %>% summarize(y = mean(y), value = mean(value), .groups = "drop_last") %>% ungroup
				
				
				
				# get the obs vs pred, calculate R2
				obs_pred <- obs_pred %>% bind_rows(my_pred) 
				my_r2 <- obs_pred %>% group_by(fit) %>% dplyr::summarize(r2 = coef_det(y, value), .groups = "drop_last") 
			}
			
			# print running tally
			if(nrow(my_r2)>0 & ii%%10==0){
				print(paste(ii, "of", nfolds, "folds; R2_both = ", round(my_r2$r2[my_r2$fit=="pred_b"], 2),
							";  R2_phy = ", round(my_r2$r2[my_r2$fit=="pred_p"], 2),
							";  R2_env = ", round(my_r2$r2[my_r2$fit=="pred_e"], 2)))
			}
			
			
		}
		
		# add the results and recalculate r2
		if(nrow(obs_pred)>0){
			ope <- obs_pred %>% filter(fit == "pred_e") %>% select(-fit)
			r2e <- ope %>% summarize(r2 = coef_det(y, value)) %>% unlist()
			
			opp <- obs_pred %>% filter(fit == "pred_p") %>% select(-fit)
			r2p <- opp %>% summarize(r2 = coef_det(y, value)) %>% unlist()
			
			opb <- obs_pred %>% filter(fit == "pred_b") %>% select(-fit)
			r2b <- opb %>% summarize(r2 = coef_det(y, value)) %>% unlist()
		}else{
			r2e <- r2p <- r2b <- r2r <- -Inf
			ope <- opp <- opb <- opr <- tibble(y = NA, value = NA)
		}
		cat("Fitting full model\n")
	}else{
		cat("Fitting full model ONLY\n")
		r2e <- r2p <- r2b <- r2r <- ope <- opp <- opb <- opr <- NULL
	}
	
	set.seed(10)

	# env + phylo full model
	full_mod_b <- ranger(formula = y~., data=df_both %>% select(y, all_of(bn_all)) %>% setNames(gsub("BOTH_trait", "trait", names(.))),
						 quantreg = quantreg, importance = "permutation", num.trees = 500, scale.permutation.importance = TRUE, num.threads = num_threads)
	
	# phylogeny only
	full_mod_p <- ranger(formula = y~., data=df_phy %>% select(y, all_of(pn_all)) %>% setNames(gsub("PHY_trait", "trait", names(.))),
						 quantreg = quantreg, importance = "permutation", num.trees = 500, scale.permutation.importance = TRUE, num.threads = num_threads)
	
	# georeferenced phylogeny only. not used.
	full_mod_e <- ranger(formula = y~., data=df_env %>% select(y, all_of(en_all)) %>% setNames(gsub("ENV_trait", "trait", names(.))),
						 quantreg = quantreg, importance = "permutation", num.trees = 500, scale.permutation.importance = TRUE, num.threads = num_threads)
	
	
	# return the model and the fit
	return(list(
		env = list(model = full_mod_e, r2 =  r2e, obs_pred = ope, quantreg = quantreg, fit_type = quantreg, obs = df_env$y, pred = full_mod_e$predictions),
		phy = list(model = full_mod_p, r2 =  r2p, obs_pred = opp, quantreg = quantreg, fit_type = quantreg, obs = df_phy$y, pred = full_mod_p$predictions),
		both = list(model = full_mod_b, r2 =  r2b, obs_pred = opb, quantreg = quantreg, fit_type = quantreg, obs = df_both$y, pred = full_mod_b$predictions))) 
	
}


# impute traits based on a random forest fit
get_prediction <- function(my_ranger_fit, my_data, type, qlevel = 0.90){
	
	# read in the file
	my_mod <- readRDS(my_ranger_fit)
	
	# get the model names
	my_name <- names(my_mod)
	my_name_se <- paste0(my_name, "_se")
	
	# get the model
	my_mod <- my_mod[[1]]
	
	# is it quantile regression or not
	if(my_mod$quantreg & type == "double"){
		# get the prediction, se,  and wrap into df
		pred <- predict(my_mod$model, data = my_data, type = "quantiles", quantiles = c(qlevel))
		pred_df <- tibble(value = pred$predictions[,1], trait = my_name)
	}else{
		# get the prediction, se,  and wrap into df				
		pred <- predict(my_mod$model, data = my_data)
		pred_df <- tibble(value = pred$predictions, trait = my_name)
	}
	
	# output
	return(pred_df)
}


# plot the species-weighted principal components
ggpca_weighted <- function(tr_use, trait_vars, num_show, axes = c(1,2), alpha = 0.1, flip_coord = FALSE, arrow_head = 6, 
						   flip_x = 1, flip_y = 1, use_vars = NULL, label_size = 5, show_labs = TRUE, alpha_const = 1,
						   alpha_all = 0, line_wd_big = 1.5, all_alpha = 0.2, angio_list = NULL, flip_z = 1, flip_zz = 1,
						   show_legend = TRUE, title = ""){
	
	
	# browser()
	# set seed to ensure coordinates are consistent
	set.seed(10)
	
	# flip the coordinates if needed
	if(flip_coord){
		tmp <- flip_x
		flip_x <- flip_y
		flip_y <- tmp
	}
	
	# keep only the unique observations, and add in some noise to allow for merge later on
	tr_unique <- tr_use %>% select(wt, all_of(trait_vars)) %>%
		mutate(rowid = 1:nrow(.)) %>% 
		gather(trait, value, -wt, -rowid) %>% 
		rowwise() %>% 
		mutate(value = value+runif(1, -1e-12,1e-12)) %>%
		ungroup %>% 
		spread(trait, value) %>% arrange(rowid) %>% select(-rowid)
	
	# keep only the traits
	tr_use <- tr_use %>% select(-all_of(trait_vars), -wt) %>% bind_cols(tr_unique)
	
	# get the trait matrix and weights
	x_pca <- tr_unique %>% select(-wt) %>% as.matrix()
	w_pca <- tr_unique$wt
	
	# run the weighted pca
	pca_use <- wpca(x = x_pca, w = w_pca, center=TRUE, scale=TRUE)
	
	# get the mean value to scale the eigenvectors to mean unit variance
	sdmean <- sqrt(mean(pca_use$d^2))
	
	# extract the objects, scale the eigenvalues and x matrix
	sdev <- pca_use$d/sdmean
	rotmat <- t(pca_use$vt)
	xmat <- pca_use$pc/sdmean
	
	# flip the rotations to align the plots
	if(flip_x == -1){
		xmat[,1] <- -xmat[,1]
		rotmat[,1] <- -rotmat[,1]
	}
	if(flip_y == -1){
		xmat[,2] <- -xmat[,2]
		rotmat[,2] <- -rotmat[,2]		
	}
	if(flip_z == -1){
		xmat[,3] <- -xmat[,3]
		rotmat[,3] <- -rotmat[,3]		
	}
	if(flip_zz == -1){
		xmat[,4] <- -xmat[,4]
		rotmat[,4] <- -rotmat[,4]		
	}
	
	# add in names
	rownames(rotmat) <- colnames(x_pca)
	colnames(rotmat) <- paste0("PC",1:ncol(rotmat))
	colnames(xmat) <- paste0("PC",1:ncol(xmat))
	
	if(num_show<0){
		num_show <- nrow(rotmat)
	}
	
	# make sure wearen't displaying more traits than we have
	num_show <- min(num_show, nrow(rotmat))
	
	# get the percentage explained, for plotting on the axes
	percvar <- round((sdev^2)/sum(sdev^2), 3)*100
	
	# get the dataframe for plotting points
	pto <- tr_use %>% left_join(bind_cols(tr_unique %>% select(-wt), xmat %>% data.frame() %>% as_tibble() %>% select(contains("PC", ignore.case = FALSE))), by = trait_vars) %>% mutate(genus = word(accepted_bin, 1)) 
	
	# scale the axes to allow direct overlay of the vectors
	pt_scaled <- pto %>% select(-wt) %>% left_join(angio_list) %>% 
		gather(pc, value, -accepted_bin, -angio, -genus, -LAT, -LON) %>% 
		group_by(pc) %>% 
		mutate(value = ifelse(grepl("PC", pc, ignore.case = FALSE), scale11(value), value)) %>% ungroup %>% 
		spread(pc, value) %>% select(-LAT, -LON)
	
	# the final dataset for plotting, with angio vs. gymno added
	pt_use <-  pt_scaled %>% 
		rename(Genus = genus) %>% mutate(angio = ifelse(angio == 1, "Angiosperm", "Gymnosperm"))
	
	# get the loadings
	loads <- loads_orig <- rotmat%*%diag(sdev) %>% data.frame() %>%
		setNames(paste0("PC",1:ncol(.))) %>% mutate(trait = rownames(.)) %>% 
		as_tibble() %>% select(trait, names(.)) 
	
	# get the data frame for the arrows
	arrow_df <- loads %>% select(trait, paste0("PC",axes)) %>% setNames(c("trait","x","y")) %>% mutate(Genus = NA) 
	
	# if we haven't specified which variables to plot (Fig 2a), calculate the ones with the highest loadings
	if(is.null(use_vars)){
		# variable closest to the x or y axis
		best <- arrow_df %>% rowwise() %>%
			mutate(y_rat = min(dist(rbind(c(x,y), c(0,1))), dist(rbind(c(x,y), c(0,-1))))) %>% 
			mutate(x_rat = min(dist(rbind(c(x,y), c(1,0))), dist(rbind(c(x,y), c(-1,0))))) %>% 
			arrange(y_rat) %>% 
			add_column(axis1 = c(rep(1, num_show), rep(0, nrow(.) - num_show))) %>% 
			arrange(x_rat) %>% 
			add_column(axis2 = c(rep(2, num_show), rep(0, nrow(.) - num_show))) %>% 
			rowwise() %>% 
			mutate(axis = ifelse((axis1 + axis2) == 0, 0, ifelse((axis1 + axis2)%in%c(1,2), max(axis1, axis2), min(axis1, axis2)))) %>% 
			ungroup %>% 
			select(axis, trait)
	}else{
		best <- tibble(trait = use_vars, axis = axes[1])
	}
	
	# get the colors
	my_color <- my_fill <- c(redmonder.pal(8,"qPBI")[c(1,3, 5,8)],redmonder.pal(8,"qMSOMed")[c(4)])
	my_shp <- c(15,16,17,18,19)
	
	# create a copy for plotting
	plot_dt <- pt_use 
	
	# scale the transparency
	alpha_scale <- plot_dt  %>% group_by(Genus) %>% tally() %>% ungroup %>% mutate(n = n%/%1000) %>% mutate(al = 0.1) %>% select(Genus, al) %>% deframe
	
	# create the base plot
	g1 <- ggplot(data = plot_dt %>% 
				 	select(angio, PC1, PC2) %>% filter(complete.cases(.), abs(PC1)< 1, abs(PC2)<1) %>% rename(Genus = angio), aes(x = PC1, y = PC2))+
		geom_point(aes(color = Genus, fill = Genus, shape = Genus, alpha = Genus), size = 0.35)+
		theme_bw()
	
	# create the transparancy value based on number of observations (note that Genus is set to angio vs. gymno in the main figure)
	alv <- plot_dt %>% select(angio, PC1, PC2) %>% filter(complete.cases(.), abs(PC1)< 1, abs(PC2)<1) %>% rename(Genus = angio) %>% 
		group_by(Genus) %>% tally() %>% ungroup %>%
		rowwise() %>% mutate(n = max(n, 15001)) %>% ungroup %>% 
		mutate(al = 1/(n%/%1000)/alpha_const) %>% arrange(Genus) %>% select(al) %>% unlist() %>% as.numeric()
	
	# plotting colors/scales
	my_color <- darken(c("#1A85FF", "#D41159"), 0.3)
	my_fill <- lighten(c("#1A85FF", "#D41159"), 0)
	my_shp <- c(22, 24) #my_shp[c(1,2)]
	alpha_scale <- alv
	names(my_color) <- names(my_fill) <- names(my_shp) <- c("Angiosperm", "Gymnosperm")
	
	# set the transparency depending on if we're plotting both angio and gymno, or just one
	if(length(alpha_scale)==1){
		names(alpha_scale) <- unique(plot_dt$angio)
	}else{
		names(alpha_scale) <- c("Angiosperm", "Gymnosperm")
	}
	
	if(num_show < nrow(rotmat)){
		#plot
		g1 <- g1 + 
			geom_hline(yintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
			geom_vline(xintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
			xlab(paste0("PC Axis ",axes[1]," (",percvar[axes[1]],"%)"))+
			ylab(paste0("PC Axis ",axes[2]," (",percvar[axes[2]],"%)"))+
			geom_segment(data = arrow_df %>% left_join(best, by = "trait") %>% filter(axis%in%axes), aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(arrow_head, "pt")), size = line_wd_big, color = "black") + 
			geom_text_repel(data = arrow_df %>% left_join(best, by = "trait") %>% filter(axis%in%axes),  aes(x = x, y = y, label = trait), fontface = "bold", box.padding = 0.75, segment.size = 0, size = label_size*1)+
			geom_segment(data = arrow_df %>% left_join(best, by = "trait") %>% filter(!axis%in%axes), aes(x = 0, y = 0, xend = x, yend = y), size = 0.7, linetype = 1, color  = "gray40", alpha = 0.4)+
			scale_fill_manual(values = my_fill) +
			scale_color_manual(values = my_color) +
			scale_shape_manual(values = my_shp) +
			scale_alpha_manual(values = alpha_scale)
	}else{
		g1 <- g1 + 
			geom_hline(yintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
			geom_vline(xintercept = 0, linetype = 2, size = 0.3, color = "gray50")+
			xlab(paste0("PC Axis ",axes[1]," (",percvar[axes[1]],"%)"))+
			ylab(paste0("PC Axis ",axes[2]," (",percvar[axes[2]],"%)"))+
			geom_segment(data = arrow_df %>% left_join(best, by = "trait"), aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(arrow_head, "pt")), size = line_wd_big, color = "black") + 
			geom_text_repel(data = arrow_df %>% left_join(best, by = "trait"),  aes(x = x, y = y, label = trait), fontface = "bold", box.padding = 0.75, segment.size = 0, size = label_size*1)+
			geom_segment(data = arrow_df %>% left_join(best, by = "trait"), aes(x = 0, y = 0, xend = x, yend = y), size = 0.7, linetype = 1, color  = "gray40", alpha = 0.4)+
			scale_fill_manual(values = my_fill) +
			scale_color_manual(values = my_color) +
			scale_shape_manual(values = my_shp) +
			scale_alpha_manual(values = alpha_scale)
	}
	# add in a title if desired
	if(!is.null(title)){
		g1 <- g1 + ggtitle(title)
	}
	
	# add the legend
	if(!show_legend){
		g1 <- g1 + theme(legend.position = "none", plot.title = element_text(face = "bold"))
	}else{
		g1 <- g1 + theme(legend.position = c(0.8,0.1), 
						 legend.title = element_blank(), 
						 legend.text = element_text(size = 11),
						 legend.background = element_rect(colour = 'black', fill = NA, linetype='solid', size = 0.2),
						 legend.spacing.y = unit(0, "mm"))+
			guides(fill = guide_legend(override.aes = list(size = 2)), alpha = FALSE)
		
	}
	
	# flip the coordinates to align the figures
	if(flip_coord){
		g1 <- g1 + coord_flip()
	}
	
	# plot
	show(g1)
	
	return(list(plot = g1, vars = best %>% filter(axis%in%axes) %>% select(trait) %>% unlist(), pca = pca_use, tr = pto, loadings = loads_orig, axis_perc = sdev^2/sum(sdev^2)))
}



