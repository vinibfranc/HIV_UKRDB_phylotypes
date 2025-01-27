libs_load <- c("glue","dplyr","forcats", "ape", "lubridate", "ggplot2",  "phytools", "data.table", "motmot", "geiger", "phylolm")
invisible( lapply(libs_load, library, character.only=TRUE) )

### 1. Internode lengths, used as estimates of maximum transmission intervals, based on: ###
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000590 and
# https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050050

RES_PATH <- "results"
RDS_PATH <- "rds"
# use merge_para_timetr from 04 that resolve paraphyletic PTs by tracing the MRCA of all sequences in the phylotype until a monophyletic group was formed

merge_para_timetr <- readRDS(glue("{RDS_PATH}/merge_para_timetr.rds"))
options(scipen=999)

# for mcs=30 and subtype B only
subtype_b_vois_ids_only <- c(40,69,133)
timetr_pt <- timetr_pt_only_internal <- list() #node_heights <- internode_intervals <- internode_intervals_df <- list()
internode_lengths <- internode_intervals <- internode_intervals_df <- list()
for(i in 1:length(merge_para_timetr[[1,4]])) {
	print(i)
	if(!is.null(merge_para_timetr[[1,4]][[i]])) {
		timetr_pt[[i]] <- merge_para_timetr[[1,4]][[i]]
		# get only node to node edges (ignore node to tip edges)
		timetr_pt_only_internal[[i]] <- timetr_pt[[i]]$edge[ timetr_pt[[i]]$edge[,2] >= (Ntip(timetr_pt[[i]])+1), ]
		internode_lengths[[i]] <- timetr_pt[[i]]$edge.length[ timetr_pt[[i]]$edge[,2] >= (Ntip(timetr_pt[[i]])+1) ]
		
		internode_intervals_df[[i]] <- data.frame(phylotype=i, internode_intervals=internode_lengths[[i]], status=ifelse(i %in% subtype_b_vois_ids_only, yes="VOI",no="Non-VOI"))
	}
}

# tree with tips
plot(merge_para_timetr[[1,4]][[69]]); axisPhylo(root.time=max( node.depth.edgelength(merge_para_timetr[[1,4]][[69]]) ), backward = T)

ii_all_df <- rbindlist(internode_intervals_df)

# extract mean, median, and max internode intervals for each phylotype & VOI status
ii_summary_pt <- ii_all_df %>% group_by(phylotype) %>% summarise(mean=mean(internode_intervals), median=median(internode_intervals), max=max(internode_intervals), status=unique(status))
# convert to months instead of years
ii_all_df$internode_intervals_months <- ii_all_df$internode_intervals * 12
ii_summary_pt2 <- ii_all_df %>% group_by(phylotype) %>% summarise(mean=mean(internode_intervals_months), median=median(internode_intervals_months), max=max(internode_intervals_months), status=unique(status))

res_ii_path <- glue("{RESULTS_PATH}/08_internode_intervals")
system(glue("mkdir -p {res_ii_path}"))

# plot distribution (of means, medians, and max) of the above for VOI vs non-VOI phylotypes
plot_ii_summary <- function(df, var, unit_ii) {
	pl <- ggplot(df, aes(x = !!sym(var), fill = status)) +
		geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
		labs(x = glue("Internode intervals ({unit_ii})"), y = "Count") +
		theme_classic()
	ggsave(plot=pl, file=glue("{res_ii_path}/{var}_distr_{unit_ii}.jpg"), width=8, height=6.5, dpi=300)
}
# years
plot_ii_summary(ii_summary_pt, "mean", "years")
plot_ii_summary(ii_summary_pt, "median", "years")
plot_ii_summary(ii_summary_pt, "max", "years")
# months
plot_ii_summary(ii_summary_pt2, "mean", "months")
plot_ii_summary(ii_summary_pt2, "median", "months")
plot_ii_summary(ii_summary_pt2, "max", "months")

# summary based on grouping by status only across all observations inside each phylotype (months)
ii_summary_voi_st <- ii_all_df %>% group_by(status) %>% summarise(mean=mean(internode_intervals_months), median=median(internode_intervals_months), max=max(internode_intervals_months))
print(ii_summary_voi_st)
# status   mean median   max
# Non-VOI  17.8   6.90  321.
# VOI      18.7   8.28  132.
calc_iqr2 <- function(data) {
	lower_bound <- quantile(data, 0.25, na.rm=T)  # 25th percentile
	upper_bound <- quantile(data, 0.75, na.rm=T)  # 75th percentile
	print(paste0("Lower Bound (25th Percentile):", round(lower_bound,3)))
	print(paste0("Upper Bound (75th Percentile):", round(upper_bound,3)))
}

# to use alongside medians
calc_iqr2(ii_all_df$internode_intervals_months[ii_all_df$status == "Non-VOI"]) #0.002-23.778
calc_iqr2(ii_all_df$internode_intervals_months[ii_all_df$status == "VOI"]) #0.002-26.72

# to use with means
sd(ii_all_df$internode_intervals_months[ii_all_df$status == "Non-VOI"]) #29.0
sd(ii_all_df$internode_intervals_months[ii_all_df$status == "VOI"]) #26.6

# in previous work median transmission intervals were 14 mo for MSM and 27 mo for het
# "2% of transmissions in this population were estimated to occurr within 6 months or less"
# Here: means 18 (non-VOIs) and 19 (VOIs), medians 7 (non-VOIs) and 8 (VOIs), max 321 (~27 y, non-VOIs) and 132 (~11y for VOIs)
# Comparison needs to be against MSM, our means are similar to their medians, but our medians show much lower values
# Anyway values for VOIs seem a bit lower, suggesting faster transmission (but maybe should compare against non-VOIs with matched demog data instead)
# PERHAPS could plot internode intervals (y) over time coloured by VOI status OR ii (y) by phylotype

### 2. Pagel's Lambda / phylo-informed glm testing association between genetic distance and viral load ###

# Pagel’s lambda is a measure of phylogenetic ‘signal’ in which the degree to which shared history of taxa has driven trait distributions at tips. 
# In this model, internal branch lengths are transformed by the lambda parameter value. When the parameter lambda equals 1, 
# branches are transformed by multiplying by 1 and so the model is equal to Brownian motion (high phylogenetic signal). 
# Values of lambda under 1 suggest there has been less influence of shared history on trait values at the tips. 
# Finally, a lambda value of 0 indicates no phylogenetic influence on trait distributions, 
# and is equivalent to a ‘star phylogeny’ with no shared branch lengths.

# estimate lambda and plot out with CIs for VOIs/non-VOIs coloured
res_pagels_lambda_mm <- glue("{RES_PATH}/08_lambda_mm")
res_pagels_lambda_geig_fc <- glue("{RES_PATH}/08_lambda_geig_fc")
res_pagels_lambda_phylolm <- glue("{RES_PATH}/08_lambda_phylolm")
system(glue("mkdir -p {res_pagels_lambda_mm}"))
system(glue("mkdir -p {res_pagels_lambda_geig_fc}"))
system(glue("mkdir -p {res_pagels_lambda_phylolm}"))
NCPU <- 1

b_ml_tree <- readRDS(glue("{RDS_PATH}/tree_b_adj2.rds"))# ML tree
vl_b <- readRDS(glue("{RDS_PATH}/vl_subtypes_comb2.rds"))[[1,4]] # mean VLs for treestructure phylotypes
vl_b_trait_df <- vl_b %>% dplyr::select(taxon, mean_log_vl_pat, cluster)
vl_b_trait_df2 <- as.data.frame(vl_b_trait_df)

# match tips with VL info
vl_b_trait_df_pt <- vl_b_trait_df[vl_b_trait_df2$taxon %in% b_ml_tree$tip.label,] #8810
tree_with_vl <- keep.tip(b_ml_tree, tip=vl_b_trait_df_pt$taxon)
vl_b_trait_mx <- as.matrix(vl_b_trait_df_pt)
rownames(vl_b_trait_mx) <- vl_b_trait_mx[,1]
vl_b_trait_mx <- vl_b_trait_mx[,-1]
rownames_saved <- rownames(vl_b_trait_mx)
vl_b_trait_mx <- apply(vl_b_trait_mx, 2, as.numeric)
rownames(vl_b_trait_mx) <- rownames_saved
trait_data <- sortTraitData(phy=tree_with_vl, y=vl_b_trait_mx, data.name="mean_log_vl_pat", log.trait = F)
phy <- trait_data$phy
vl <- trait_data$trait

# fit model on continuous trait (VL) - bayes option geiger::fitContinuousMCMC (only BM and other models, but not lambda)
geig_fc <- geiger::fitContinuous(phy=phy, dat = vl, model = "lambda", ncores=NCPU) # hessian should return CIs (but doesn't work): niter=1000, hessian=TRUE, CI=0.95 
geig_fc_df <- data.frame(lambda_mean=geig_fc$opt$lambda)
print(geig_fc_df) #lambda 0.13

### using phylolm for directly estimating parameters like lambda and integrating it into linear models

vl_b_trait_df3 <- vl_b_trait_df2
rownames(vl_b_trait_df3) <- vl_b_trait_df3$taxon

phy2 <- phy
phy2$edge.length[phy2$edge.length == 0] <- 1e-8
			
# lambda from viral load: crude estimate of heritability in the tree
phylolm_lambda <- tryCatch({
	phylolm(mean_log_vl_pat ~ 1, data = vl_b_trait_df3, phy = phy2, model = "lambda", boot=100) #sexid + ethnicityid + exposureid + age_group
}, error = function(e) {
	cat("phylolm error ", e$message, "\n")
	return(list(NA))
})

# phylolm r2
phylolm_lambda_cis <- data.frame(adj.r.squared = phylolm_lambda$adj.r.squared, 
																																				lambda_bootmean=unname(phylolm_lambda$bootmean["lambda"]),
																																				lambda_lower=unname(phylolm_lambda$bootconfint95[1,"lambda"]),
																																				lambda_upper=unname(phylolm_lambda$bootconfint95[2,"lambda"]),
																																				n = phylolm_lambda$n )
phylolm_lambda_cis
# 0.1332465    0.1054192    0.1571757