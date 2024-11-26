libs_load <- c("ggplot2","dplyr","ggpubr","data.table","glue","lme4","lubridate","viridis")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"

extracted_clusters <- readRDS(glue("{RDS_PATH}/extracted_clusters.rds"))
subtype_choices <- c("A (A1)","CRF 02_AG","C", "B")

### STEP 9: PREPARE VL AND CD4 DATA ###
# CD4 decline data preparation
demog_md <- readRDS(glue("{RDS_PATH}/demog_md.rds"))
subtype_md <- readRDS(glue("{RDS_PATH}/subtype_md.rds"))
demog_md_subtype_match <- readRDS(glue("{RDS_PATH}/demog_md_subtype_match.rds")) #100591
demog_md_subtype_match <- demog_md_subtype_match[ (demog_md_subtype_match$status=="Naïve") & 
																																																			(demog_md_subtype_match$rega3subtype %in% c("A (A1)","CRF 02_AG","C","B")) & 
																																																			(demog_md_subtype_match$exposureid != "Not known"), ] #53382

### 9.1. Viral loads ###
vl_md <- read.csv(glue("{DATA_PATH}/vloads.csv", header=T)) #995260 rows,34757 unique patients
vl_md <- vl_md[vl_md$vl > 0,]; print(nrow(vl_md)); print(length(unique(vl_md$patientindex))) #976418 rows, 34740 patients
vl_md$log_vl <- log10(vl_md$vl) 
vl_md$vl_md_ymd <- as.Date(gsub("\\/", "15", vl_md$vl_date_my), "%m%d%Y")
saveRDS(vl_md, glue("{RDS_PATH}/vl_md.rds"))
print(nrow(vl_md[vl_md$onartflag==1,])) #781806
print(nrow(vl_md[vl_md$onartflag==0,])) #194612
# Get only viral load before treatment
vl_md_preart <- vl_md[vl_md$onartflag == 0,]
vl_md_preart$vl_decimal_date <- decimal_date(as.Date(vl_md_preart$vl_md_ymd))
print(nrow(vl_md_preart)); print(length(unique(vl_md_preart$patientindex))) #194612, 31091

# below 200 -> 2.301 log10 (viral suppression) -> 18778 rows, 5949 patients
vl_md_preart_less_200 <- vl_md_preart[vl_md_preart$vl < 200,]; print(nrow(vl_md_preart_less_200)); print(length(unique(vl_md_preart_less_200$patientindex)))
# >= 10,000,000 -> 7 log10 (above level of viral quant assay precision) -> 148 rows, 130 patients
vl_md_preart_more_10m <- vl_md_preart[vl_md_preart$vl > 10000000,]; print(nrow(vl_md_preart_more_10m)); print(length(unique(vl_md_preart_more_10m$patientindex)))

# filter out these probably wrong measurements (175686 rows, 30477 patients)
vl_md_preart <- vl_md_preart[vl_md_preart$vl >= 200 & vl_md_preart$vl <= 10000000,]
print(nrow(vl_md_preart)); print(length(unique(vl_md_preart$patientindex)))

saveRDS(vl_md_preart, glue("{RDS_PATH}/vl_md_preart.rds"))

# Get viral loads between date of diagnosis and 2 years after for diff subtypes
get_vl <- function(demog_md_choice, subtype_choice) {
	
	demog_md_choice <- demog_md_choice[demog_md_choice$rega3subtype == subtype_choice,]
	
	demog_md_choice$hivpos_ymd <- as.Date(gsub("\\/", "15", demog_md_choice$hivpos_my), "%m%d%Y")
	demog_md_choice$hivpos_decimal_date <- decimal_date(as.Date(demog_md_choice$hivpos_ymd))
	print("nrows md")
	print(nrow(demog_md_choice))
	print("unique patients")
	print(length(unique(demog_md_choice$patientindex)))
	
	# Old "spvl" extraction
	#spvl_md <- demog_md_choice %>% left_join(vl_md_preart, by="patientindex", relationship="many-to-many") %>% filter( (vl_decimal_date - hivpos_decimal_date >= 0.5) & (vl_decimal_date - hivpos_decimal_date) <= 2)
	spvl_md <- demog_md_choice %>% left_join(vl_md_preart, by="patientindex", relationship="many-to-many") %>% filter( (vl_decimal_date - hivpos_decimal_date >= 0) & (vl_decimal_date - hivpos_decimal_date) <= 2)
	#print(nrow(spvl_md))
	print("unique patient vl pre-art")
	print(length(unique(spvl_md$patientindex)))
	
	lookup_pat_tests <- unique( spvl_md[ , c('patientindex','testindex') ] )
	lookup_pat_tests$testindex <- paste0("t.",lookup_pat_tests$testindex)
	
	spvl_md_cl <- spvl_md %>% distinct(patientindex, vl_decimal_date, .keep_all = TRUE)
	
	spvl_md_cl <- spvl_md_cl[ , !names(spvl_md_cl) %in%  c("testindex","dbsample_date")]
	
	list(spvl_md_cl, lookup_pat_tests)
}

min_cl_size_choices <- c(30, 50, 100) #250,500
tree_names <- c("A_A1","CRF_02_AG","C","B")

vl_subtypes <- list()
for(i in 1:length(tree_names)) {
	print(tree_names[[i]])
	vl_subtypes[[i]] <- get_vl(demog_md_subtype_match, subtype_choices[i])
}
#saveRDS(vl_subtypes, glue("{RDS_PATH}/vl_subtypes.rds"))
vl_subtypes <- readRDS("rds/vl_subtypes.rds")

# Mean VL measurements before treatment (without excluding phylotypes with less than e.g. 10 measurements)
vl_subtypes_all <- list(vl_subtypes[[1]][[1]], vl_subtypes[[2]][[1]], vl_subtypes[[3]][[1]], vl_subtypes[[4]][[1]])
vl_subtypes_all <- bind_rows(vl_subtypes_all); print(nrow(vl_subtypes_all)); print(length(unique(vl_subtypes_all$patientindex))) #50667 rows, 15702 unique patients 
spvl_measur <- vl_subtypes_all %>% add_count(patientindex, name="n_spvl") %>% group_by(patientindex) %>% filter(row_number() >= (n())) #filter(n() == 1)
mean(spvl_measur$n_spvl) # 3.226
sd(spvl_measur$n_spvl) # 2.533
median(spvl_measur$n_spvl) # 2
IQR(spvl_measur$n_spvl); table(spvl_measur$n_spvl); hist(spvl_measur$n_spvl) # 4

treestruct_min_cl_size_res_yes_sup <- readRDS(glue("{RDS_PATH}/treestruct_min_cl_size_res_yes_sup.rds")) 

# create matrix for each mcs and subtype + summarise log_vl to get only one per patient (treestructure phylotypes)
vl_subtypes_comb1 <- vl_subtypes_comb2 <- vl_subtypes_mean_p_pat <- vl_subtypes_comb3 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		print("treestructure matching")
		print(nrow(treestruct_min_cl_size_res_yes_sup[[i,j]][[2]]))
		print(nrow(vl_subtypes[[j]][[2]]))
		vl_subtypes_comb1[[i,j]] <- inner_join(treestruct_min_cl_size_res_yes_sup[[i,j]][[2]], vl_subtypes[[j]][[2]], by=c("taxon"="testindex"))
		print(nrow(vl_subtypes_comb1[[i,j]]))
		vl_subtypes_mean_p_pat[[i,j]] <- vl_subtypes[[j]][[1]] %>% group_by(patientindex) %>% summarise(mean_log_vl_pat=mean(log_vl))
		print(nrow(vl_subtypes_mean_p_pat[[i,j]]))
		vl_subtypes_comb2[[i,j]] <- inner_join(vl_subtypes_comb1[[i,j]], vl_subtypes_mean_p_pat[[i,j]], by="patientindex")
		print(nrow(vl_subtypes_comb2[[i,j]]))
		print(nrow(vl_subtypes_comb2[[i,j]][ is.na(vl_subtypes_comb2[[i,j]]$mean_log_vl_pat), ]))
		print("===")
	}
}
saveRDS(vl_subtypes_comb2, glue("{RDS_PATH}/vl_subtypes_comb2.rds"))
# A1: 757 patients, CRF_02_AG: 935, C: 2371, B: 8810

# Estimating group-level differences using a multi-level Bayesian model (thanks Chris Wymant for the advice: https://htmlpreview.github.io/?https%3A//github.com/ChrisHIV/teaching/blob/main/other_topics/Stan_example_hierarchical_model_false_positives.html)
libs_load2 <- c("tidyverse","rstan","ggforce", "Rcpp", "bayesplot")
invisible( lapply(libs_load2, library, character.only=TRUE) )

options(mc.cores = parallel::detectCores()) # parallelise
rstan_options(auto_write = TRUE)            # avoid re-compiling stan code
theme_set(theme_classic())
set.seed(123946)

ITER <- 4000 #TODO change to 10000
CHAINS <- 2 #TODO change to 2
model_compiled <- stan_model(file="scripts/stan/vl_model.stan")

vl_bayes_model <- function(vl_subtypes_c, mcs_subtype_choice, partition_method_lbl, model_compiled) {
	# Parameters that we condition on
	num_groups <- length(unique(vl_subtypes_c$cluster))
	print("num groups")
	print(num_groups)
	
	# Parameters that we estimate and that are 'top level' (they do not have
	# distributions controlled by other parameters to be estimated)
	y_mean_pop <- mean(vl_subtypes_c$mean_log_vl_pat)
	print("Mean pop")
	print(y_mean_pop)
	y_sd <- sd(vl_subtypes_c$mean_log_vl_pat)
	print("SD pop")
	print(y_sd)
	# y_sd_group <- 0.5???
	
	# TODO uncomment: exploratory plots
	pl1 <- ggplot(vl_subtypes_c, aes(x=as.factor(cluster), y=mean_log_vl_pat)) +
		geom_violin() + geom_sina() + labs(x = "phylotype") + theme(axis.text=element_text(size=6),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	system(glue("mkdir -p {RESULTS_PATH}/05_vl_{partition_method_lbl}/"))
	ggsave(plot=pl1, filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_expl.jpg"), width=10, height=8, dpi=300)

	# The following Stan code (stored as a string in R) estimates those parameters of the data-generating process that we don’t
	# condition on. The likelihood it uses matches the one we used to generate the data.
	
	# parameters {
	# 	// Top-level parameters
	# 	real<lower = 2, upper = 7> y_mean_pop;
	# 	real<lower = 0,   upper = 2> y_sd;
	# 	real<lower = 0,   upper = 5> y_sd_group;
	# 	// Lower-level parameters (with a non-centred parameterisation for numerical
	# 																												// efficiency, see e.g.
	# 																												// https://mc-stan.org/docs/2_18/stan-users-guide/reparameterization-section.html
	# 																												vector[num_groups] group_effects_unscaled;
	# }

	vl_subtypes_c$cluster <- as.integer(vl_subtypes_c$cluster)
	
	gc()

	fit <- sampling(model_compiled,
																	data = list(
																		num_groups = max(vl_subtypes_c$cluster,na.rm=T),
																		num_y = nrow(vl_subtypes_c),
																		group = vl_subtypes_c$cluster,
																		y = vl_subtypes_c$mean_log_vl_pat,
																		sample_prior_only=0),
																	iter = ITER,
																	chains = CHAINS,
																	seed = 123)

	gc()
	
	fit_summary <- summary(fit)
	print("head oh Rhat and n_eff")
	print(head(fit_summary$summary[, "Rhat" ])) #R-hat < 1.01 -> Good convergence
	print(head(fit_summary$summary[, "n_eff" ]))
	print("Any Rhat > 1.01?")
	print(any(fit_summary$summary[, "Rhat"] > 1.01))
	print("Any n_eff < 400?")
	print(any(fit_summary$summary[, "n_eff"] < 400))

	# TODO printing empty PDF
	# pdf(glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_traceplot_model.pdf"), width = 10, height = 7)
	# traceplot(fit, pars=names(fit)[1:10], inc_warmup=TRUE) # there are a lot more, just plotting 10 to explore
	# dev.off()
	
	print("Getting df_fit")

	df_fit <- fit %>%
		as.data.frame() %>%
		as_tibble() %>%
		mutate(sample = row_number()) %>%
		pivot_longer(-sample, names_to = "param")
	
	# posterior predictive checks (PPCs)
	# Extract posterior predictive samples
	print("Doing and plotting posterior predictive checks (PPCs)")
	posterior_samples <- extract(fit)
	y_rep <- posterior_samples$y_rep
	
	rm(fit)
	gc()
	
	# Overlay observed and simulated densities
	# obs data
	observed_density <- ggplot() +
		geom_density(aes(x = vl_subtypes_c$mean_log_vl_pat), fill = "blue", alpha = 0.5) +
		labs(title = "Posterior Predictive Check", x = "y", y = "Density")
	
	# sim data (50 first iters)
	y_rep_long <- data.frame(value = as.vector(y_rep[1:50, ]))
	pl2 <- observed_density +
		geom_density(data = y_rep_long, aes(x = value), color = "red", alpha = 0.4)
	
	ggsave(plot=pl2, filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_ppc_sim_dens_.jpg"), width=8, height=6, dpi=300)
	
	# PPC with bayesplot
	ppc_dens_overlay(vl_subtypes_c$mean_log_vl_pat, y_rep[1:50, ])
	ggsave(filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_ppc_dens_overlay_bp_.jpg"), width=8, height=6, dpi=300)
	ppc_stat(vl_subtypes_c$mean_log_vl_pat, y_rep, stat = "mean")
	ggsave(filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_ppc_stat_mean_bp_.jpg"), width=8, height=6, dpi=300)
	ppc_stat(vl_subtypes_c$mean_log_vl_pat, y_rep, stat = "sd")
	ggsave(filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_ppc_stat_sd_bp_.jpg"), width=8, height=6, dpi=300)
	
	data_sample_prior <- list(
		num_groups = max(vl_subtypes_c$cluster,na.rm=T),
		num_y = nrow(vl_subtypes_c),
		group = vl_subtypes_c$cluster,
		y = vl_subtypes_c$mean_log_vl_pat,
		sample_prior_only=1)
	
	gc()
	
	print("Sampling from the prior only")
	
	# Sample from the prior
	fit_prior <- sampling(model_compiled, data = data_sample_prior, iter = ITER, chains = CHAINS)
	
	gc()
	
	print("Getting df_prior")
	df_prior <- fit_prior %>%
		as.data.frame() %>%
		as_tibble() %>%
		mutate(sample = row_number()) %>%
		pivot_longer(-sample, names_to = "param")
	rm(fit_prior)
	gc()
	
	# Combine prior and posterior into one data frame
	df_prior$source <- "Prior"
	#View(df_prior)
	df_fit$source <- "Posterior"
	#View(df_fit)
	
	df_combined <- bind_rows(df_prior, df_fit)
	#View(df_combined)
	print(length(unique(df_combined$param)))
	df_combined2 <- subset(
		df_combined, 
		param != "lp__" & 
			!grepl("^group_effects_unscaled", param) & 
			!grepl("^y_mean\\[", param) & 
			!grepl("^y_rep", param)
	)
	print(length(unique(df_combined2$param)))
	
	# For each parameter, plot its marginal posterior overlayed with prior
	# plots of posterior distributions such as those above should always show the prior distribution as well,
	# to show how much information came from the data and how much from the prior: we always want to know this
	print("For each parameter, plotting its marginal posterior overlayed with prior")
	pl3 <- ggplot(df_combined2) +
		geom_density(aes(value, fill = source), alpha = 0.5) +
		facet_wrap(~param, scales = "free") +
		labs(x = "Parameter value", y = "Density", fill = "Source")
	
	ggsave(plot=pl3, filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_posterior.jpg"), width=18, height=12, dpi=300)

	# For each group, calculate the two-tailed Bayesian p value that it is different from the population-level mean, i.e.
	# twice the posterior mass for the group effect being less than zero (for positive effects) or greater than zero (for negative effects)
	print("Calculating the two-tailed Bayesian p value that it is different from the population-level mean")
	df_p <- df_fit %>%
		filter(startsWith(param, "group_effects_unscaled")) %>%
		summarise(.by = param,
												p = 2 * min(sum(value > 0), sum(value < 0)) / n()) %>%
		mutate(group = str_match(param, "group_effects_unscaled\\[([0-9]+)\\]")[,2])
	#View(df_p)

	# Calculate the estimated group effects - median and 95% credible intervals
	quantiles <- c(0.025, 0.5, 0.975)
	df_estimate <- df_fit %>%
		filter(startsWith(param, "y_by_group")) %>%
		group_by(param) %>%
		reframe(quantile = quantiles,
										value = quantile(value, quantiles)) %>%
		pivot_wider(values_from = value, names_from = quantile, names_prefix = "quantile_") %>%
		mutate(group = str_match(param, "y_by_group\\[([0-9]+)\\]")[,2])
	#View(df_estimate)
	# Now plot the estimated group effects in red, and data as black circles, ordered by the Bayesian p value for the group.

	print("Plotting estimates for each phylotype with 95% CIs")
	pl4 <- ggplot(data = df_estimate %>%
																	inner_join(df_p %>% select(group, p), by = "group") %>%
																	mutate(group = fct_reorder(group, p), significance = if_else(p < 0.05, "*", ""))) +
		geom_sina(data = vl_subtypes_c %>%
													mutate(group = as.character(cluster)) %>%
													left_join(df_p %>% select(group, p), by = "group") %>%
													mutate(group = fct_reorder(group, p)),
												aes(group, mean_log_vl_pat, ), alpha=0.15) + #y
		geom_errorbar(aes(x = group, ymin = quantile_0.025, ymax = quantile_0.975), col = "blue") +
		geom_point(aes(x = group, y = quantile_0.5), col = "blue") +
		geom_text(aes(x=group, label = significance, y = quantile_0.975 + 0.2), size = 10, color = "blue", fontface = "bold") +
		#geom_point(aes(group, y_true), shape = 4, size = 4) +
		labs(x = "phylotype") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	ggsave(plot=pl4, filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_effects_cis.jpg"), width=15, height=8, dpi=300)

	vl_subtypes_c$group <- vl_subtypes_c$cluster
	vl_subtypes_c$y <- vl_subtypes_c$mean_log_vl_pat

	fit_pairwise_lm <- function(group_) {
		lm_ <- lm(data = vl_subtypes_c, mean_log_vl_pat ~ group == group_)
		summary_ <- summary(lm_)
		confint_ <- confint(lm_)
		list(freq_estimate = summary_$coefficients[1,1] + summary_$coefficients[2,1],
							freq_p = summary_$coefficients[2,4],
							freq_lower = summary_$coefficients[1,1] + confint_[2,1],
							freq_upper = summary_$coefficients[1,1] + confint_[2,2])
	}

	print("Plotting comparison against simpler lm model")
	df_compare <- tibble(group = unique(vl_subtypes_c$cluster) %>% as.character) %>% #1:max(vl_subtypes_c$group,na.rm=T)
		mutate(result = map(group, fit_pairwise_lm)) %>%
		unnest_wider(result) %>%
		inner_join(df_estimate, by = "group") %>%
		inner_join(df_p %>%
														select(group, p) %>%
														rename(bayes_p = p), by = "group")

	pl5 <- pl4 + geom_errorbar(data = df_compare,
																aes(x = group, ymin = freq_lower, ymax = freq_upper),
																col = "red") +
		geom_point(data = df_compare, aes(x = group, y = freq_estimate), col = "red")

	ggsave(plot=pl5, filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_effectsBayes_VS_lm_.jpg"), width=15, height=8, dpi=300)

	print("Extracting significant phylotypes")
	df_table <- df_estimate %>% inner_join(df_p %>% select(group, p), by = "group")
	#print(head(df_table))
	df_table$y_mean_pop <- y_mean_pop
	# from Bayesian model, extract significant phylotypes, both HIGHER VL and LOWER VL
	sig_higher <- df_table %>% filter(p < 0.05 & quantile_0.5 > y_mean_pop)
	print("signif higher VL")
	print(sig_higher)
	write.csv(sig_higher, file=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_sig_higher_vl.csv"), quote=F, row.names = F)
	print("signif lower VL")
	sig_lower <- df_table %>% filter(p < 0.05 & quantile_0.5 < y_mean_pop)
	print(sig_lower)
	write.csv(sig_lower, file=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_sig_lower_vl.csv"), quote=F, row.names = F)

	print("Saving fits to rds if needed")
	saveRDS(df_fit, glue("{RDS_PATH}/df_fit_{partition_method_lbl}_{mcs_subtype_choice}.rds"))
 rm(df_fit)
 gc()
	saveRDS(df_prior, glue("{RDS_PATH}/df_prior_{partition_method_lbl}_{mcs_subtype_choice}.rds"))
	rm(df_prior)
	gc()
 #list(df_fit=df_fit, df_prior=df_prior, df_table=df_table, sig_higher=sig_higher, sig_lower=sig_lower, pl2=pl2, pl3=pl3, pl4=pl4, pl5=pl5)
	#list(df_fit=df_fit, df_table=df_table, sig_higher=sig_higher, sig_lower=sig_lower)
	list(df_table=df_table, sig_higher=sig_higher)
}

res_vl_model_ts <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) { #length(min_cl_size_choices)
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		res_vl_model_ts[[i,j]] <- vl_bayes_model(vl_subtypes_comb2[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}"), "treestructure", model_compiled)
	}
}

saveRDS(res_vl_model_ts, glue("{RDS_PATH}/res_vl_model.rds"))

# rm(res_vl_model)
# gc()

# top10 median estim now: 40, 101, 69, 125, 79, 20, 138, 122, 142, 137 (4 match previous SPVL-VOIs)
# previous SPVL-VOIs (not in order): 5,15,20(ok),24,26,40(ok),44,53,69(ok),81,115,121,137(ok)
# lm model seems to overestimate median and overall effect

# no signif diff for any mcs on the other 3 subtypes
# B30: 4 (40[ok],101[no],69[ok],20[ok]) higher, 5 lower
# B50: 3 higher (14[ok],27[ok],52[ok]), 2 lower
# B100: 1 higher (11[no]), 2 lower

# same model for fastbaps partitions
vl_subtypes_comb_fb <- vl_subtypes_comb_fb2 <- vl_subtypes_mean_p_pat_fb <- list()
for(j in 1:length(tree_names)) {
	print(glue("{tree_names[j]}"))
	print("fastbaps matching")
	print(nrow(fastbaps_df_list[[j]]))
	print(nrow(vl_subtypes[[j]][[2]]))
	vl_subtypes_comb_fb[[j]] <- inner_join(fastbaps_df_list[[j]], vl_subtypes[[j]][[2]], by=c("id"="testindex"))
	print(nrow(vl_subtypes_comb_fb[[j]]))
	vl_subtypes_mean_p_pat_fb[[j]] <- vl_subtypes[[j]][[1]] %>% group_by(patientindex) %>% summarise(mean_log_vl_pat=mean(log_vl))
	print(nrow(vl_subtypes_mean_p_pat_fb[[j]]))
	vl_subtypes_comb_fb2[[j]] <- inner_join(vl_subtypes_comb_fb[[j]], vl_subtypes_mean_p_pat_fb[[j]], by="patientindex")
	print(nrow(vl_subtypes_comb_fb2[[j]]))
	print(nrow(vl_subtypes_comb_fb2[[j]][ is.na(vl_subtypes_comb_fb2[[j]]$mean_log_vl_pat), ]))
	print("===")
}

res_vl_model_fb <- list()#length(min_cl_size_choices)
for(j in 1:length(tree_names)) {
	print(glue("{tree_names[j]}"))
	if(length(unique(vl_subtypes_comb_fb2[[j]]$cluster)) != 1) {
		res_vl_model_fb[[j]] <- vl_bayes_model(vl_subtypes_comb_fb2[[j]], glue("{tree_names[j]}"), "fastbaps", model_compiled)
	}
}

# no partition with signif higher or lower VL as in treestructure for A1
# for B, 3 partition have higher VL and 5 have lower VL
saveRDS(res_vl_model_fb, glue("{RDS_PATH}/res_vl_model_fb.rds"))

# Extract signif higher partitions in trestruct and fastbaps and check if same clusters or not
get_tips_cluster_id <- function(part_file, seq_id, vl_model_signif) { #cluster_id
	signif_cls <- as.integer(str_extract(vl_model_signif$param, "(?<=\\[)[0-9]+(?=\\])"))
	print(signif_cls)
	part_file$cluster <- as.integer(part_file$cluster)
	#print( part_file[[seq_id]])
	tips_clusters <- list()
	for(i in 1:length(signif_cls)) {
		print(signif_cls[i])
		#tips_clusters[[i]] <- part_file[[seq_id]][ part_file$cluster == signif_cls[i] ]
		tips_clusters[[i]] <- data.frame(tips=part_file[[seq_id]][ part_file$cluster == signif_cls[i] ], cluster=signif_cls[i])
	}
	#names(tips_clusters) <- signif_cls
	tips_clusters_all <- rbindlist(tips_clusters)
	print(tips_clusters_all)
}

res_vl_model_ts <- readRDS("rds/res_vl_model.rds")
higher_vl_B_ts <- get_tips_cluster_id(vl_subtypes_comb2[[1,4]], "taxon", res_vl_model_ts[[1,4]]$sig_higher)
# PT101 has 20 meas; PT20 has 69 meas; PT40 has 38 meas; PT69 has 26 meas
print(nrow(higher_vl_B_ts)) # 153

res_vl_model_fb <- readRDS("rds/res_vl_model_fb.rds")
higher_vl_B_fb <- get_tips_cluster_id(vl_subtypes_comb_fb2[[4]], "id", res_vl_model_fb[[4]]$sig_higher)
# PT118 has 42 meas; PT21 has 32 meas; PT3 has 40 meas
print(nrow(higher_vl_B_fb)) # 114
inner_join_vl_vois <- inner_join(higher_vl_B_ts, higher_vl_B_fb, by="tips") # 84 matches
View(inner_join_vl_vois)
# n=20 PT101 matches PT21 (+12 in fb), n=38 PT40 matches PT3 (all 38 match, +2 in fb), n=26 PT69 matches PT118 (+16 in fb), 
# PT20 not matching signif from fb
# trestruct VOIs estimates (PT101, PT20, PT40, PT69): 4.90, 4.81 (no correspondence), 4.91, 4.86
# fastbaps VOIs estimates (PT21, x, PT3, PT118):      4.79, no correspondence,        4.85, 4.79
# the PT that does not match (PT20 from trestruct is paraphyletic!)