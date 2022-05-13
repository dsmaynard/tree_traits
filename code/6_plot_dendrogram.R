####################################################################################
### Code for: "Global relationships in tree functional traits" by Maynard et al. 
### Contact Dan Maynard (dan.s.maynard@gmail) if you have any questions.
###
### This script creates Figure 4. It uses the output from the PCA plotting
### script, as well as the phylogenetic contrasts/conservatism script. 
###
#############################################################################

rm(list = ls())

library(feather)
library(dendextend)
library(RColorBrewer)
library(factoextra)
library(wCorr)
library(corrplot)
library(tidyverse)

# source the functions
source("code/R_functions.R")

# read in the trait names
tn <- read_csv("data/raw_data/Trait_names.csv") %>% filter(focal_trait == 1)

# read in angio/gymno
ag <- read_csv("data/raw_data/ANGIO_GYMNO_lookup.csv") %>% select(accepted_bin, group) %>%
	mutate(angio = ifelse(group == "Angiosperms", 1, ifelse(group == "Gymnosperms",0, NA))) %>% 
	select(accepted_bin, angio)

# read in the weighted trait data, outputted from the PCA code
tr_wt <- read_feather("data/results/PCA_data_for_dendro_BOTH.feather")

# get the trait names
my_traits <- names(tr_wt)[names(tr_wt)%in%tn$trait_short]


#########################################
## Calculate species-weighted correlation

cm <- p_value <- matrix(0, length(my_traits), length(my_traits))
for(i in 1:length(my_traits)){
	for(j in i:length(my_traits)){
		if(j>i){
			cm[j,i] <- cm[i,j] <- weightedCorr(x = tr_wt %>% select(all_of(my_traits[i])) %>% unlist(), tr_wt %>% select(all_of(my_traits[j])) %>% unlist(), 
											   method = "Spearman", ML = TRUE, fast = FALSE, weights = tr_wt$wt)
			p_value[j,i] <- p_value[i,j] <- cor.test(x = tr_wt %>% select(all_of(my_traits[i])) %>% unlist(), tr_wt %>% select(all_of(my_traits[j])) %>% unlist(), 
								method = "spearman", exact = FALSE)$p.value
		}
	}
}
rownames(cm) <- colnames(cm) <- my_traits

# get non signficiant values, if interested
cbind(my_traits[which(p_value>0.05/(18*(17)/2), arr.ind = TRUE)[,1]], my_traits[which(p_value>0.05/(18*(17)/2), arr.ind = TRUE)[,2]], round(cm[which(p_value>0.001, arr.ind = TRUE)], 2)) %>% 
	data.frame() %>% 
	filter(nchar(X1)>=nchar(X2))

# save a duplicate for the correlation matrix plotting, below
cmat0 <- cm

# convert the correlation to a distance matrix, ignoring sign of the correlation
cm <- as.dist(1-abs(cm))

##############################################
### Figure 4, dendrogram and cluster analysis

method <- "average"

# cluster analysis of the correlation matrix, using average absolute value correlation
hc <- hclust(cm, method = method)

# run the hierarchical clustering to get optimal number of clusters
(clust_num <- find_k(hc)$nc)

# convert to dendrogram for plotting
dend <- as.dendrogram(hc)

# specify colors for plotting
cols <- c(brewer.pal(8, "Dark2")[1:6], c(brewer.pal(3, "Set1")[c(2)]), c(brewer.pal(3, "Set2")[3]))
col_ord <- 1:clust_num

# great the dendrogram
tc <- fviz_dend(dend, main = method,
				k_colors = cols[col_ord],
				type = "rectangle",
				rect = TRUE, 
				horiz = TRUE)

# get the names of the tips
name_ord <- rev(my_traits[order.dendrogram(dend)])

# redorder the correlation matrix
(cmat <- cmat0[name_ord, name_ord])

# plot the clustered dendrogram
col_ord <- c(7, 2, 1, 8, 4,6,3,5) 
par(mar=c(3,1,1,10))


par(mar=c(3,1,1,10))
dend %>%
	color_branches(k = clust_num, col = cols[col_ord]) %>%
	color_labels(k = clust_num, col = cols[col_ord]) %>%
	set("branches_lwd",3) %>%
	set("nodes_col", 1) %>%
	set("nodes_pch", c(19)) %>%
	set("labels_cex", 1.4) %>%
	plot(horiz = TRUE, type = "rectangle")


##############################################
## Figure 4, correlation matrix

cmat2 <- cmat


# replace the lower triangle with phylogenetic independent contrasts
pic <- read_csv("data/results/PIC_contrasts.csv") %>% as.matrix()
rownames(pic) <- colnames(pic)
pic <- pic[name_ord, name_ord]
cmat2[lower.tri(cmat2)] <- pic[lower.tri(pic)]

# ensure that the upper and lower show the same correlations
cmat2[lower.tri(cmat2)][is.na(t(cmat2)[lower.tri(cmat2)])] <- NA

# remove the idagonal
diag(cmat2) <- NA

corrplot(corr = cmat2, cl.pos = "n", method = "circle", type = "full", diag = TRUE, tl.pos = "n", tl.col = "black", na.label = "  ", tl.cex = 0.7, yaxt = 'n', outline= TRUE, col = colorRampPalette(c("white","black"))(2))


#########################################################
### Figure 4, plot of pc loadings

# read in the PCA loadings
tr_comb <- read_feather("data/results/PCA_loadings_BOTH.feather")


# sort the columns to match the dendrogram
tr_both <- tr_comb %>% arrange(match(trait, colnames(cmat))) %>% 
	arrange(desc(abs(PC1)))

# add the row/col names
pc <- tr_both %>% select(paste0("PC", 1:2)) %>% as.matrix()
rownames(pc) <- tr_both$trait

# reorder to match the dendrogram
pc <- pc[name_ord,]

# plot the loading
corrplot(pc, method = "circle", type = "full", cl.pos="n", is.corr = FALSE, diag = TRUE,  tl.col = "black", na.label = "  ", tl.cex = 0.7, yaxt = 'n', outline= TRUE, col = colorRampPalette(c("white","black"))(2))


#####################################3
### Figure 4, plot of Pagel's Lambda

# read in the traits
phy <- read_csv("data/results/phylo_conservatism.csv") %>%  
	left_join(tn %>% select(trait_label, trait_short) %>% rename(trait = trait_label)) %>% 
	filter(trait%in%tn$trait_label) %>% filter(method == "lambda")

# turn it into a correlation matrix and reorder to match the dendrogram
cm <- phy %>% select(value) %>% data.frame()
rownames(cm) <- phy$trait_short
cm <- cm[name_ord, ]
cm <- matrix(cm, ncol = 1)
rownames(cm) <- name_ord

# plot Pagel's lambda
corrplot(cm, method = "circle", type = "full", diag = TRUE, tl.col = "black", na.label = "  ", tl.cex = 0.7, outline= TRUE, col = colorRampPalette(c("white","black"))(2))#rev(colorspace::diverge_hcl(100)),"darkblue")) #colorRampPalette(c(cols[1],"white",cols[2]))(100))



