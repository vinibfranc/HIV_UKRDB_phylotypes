libs_load <- c("ggplot2","dplyr","ggpubr","ggforce","data.table","glue","lme4","lubridate","viridis", "stringr")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"

extracted_clusters <- readRDS(glue("{RDS_PATH}/extracted_clusters.rds"))
subtype_choices <- c("A (A1)","CRF 02_AG","C", "B")
tree_names <- c("A_A1","CRF_02_AG","C","B")
min_cl_size_choices <- c(30, 50, 100)

### PREPARE VL AND CD4 DATA ###
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
	demog_md_choice$artstart_ymd <- as.Date(gsub("\\/", "15", demog_md_choice$artstart_my), "%m%d%Y")
	demog_md_choice$artstart_decimal_date <- decimal_date(as.Date(demog_md_choice$artstart_ymd))
	print("nrows md")
	print(nrow(demog_md_choice))
	print("unique patients")
	print(length(unique(demog_md_choice$patientindex)))
	
	# Old "spvl" extraction
	#spvl_md <- demog_md_choice %>% left_join(vl_md_preart, by="patientindex", relationship="many-to-many") %>% filter( (vl_decimal_date - hivpos_decimal_date >= 0.5) & (vl_decimal_date - hivpos_decimal_date <= 2))
	spvl_md <- demog_md_choice %>% left_join(vl_md_preart, by="patientindex", relationship="many-to-many") %>% 
		filter( (vl_decimal_date - hivpos_decimal_date >= 0) & (vl_decimal_date - hivpos_decimal_date <= 2) & (vl_decimal_date <= artstart_decimal_date)) 
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
saveRDS(vl_subtypes, glue("{RDS_PATH}/vl_subtypes.rds"))
vl_subtypes <- readRDS("rds/vl_subtypes.rds")

treestruct_min_cl_size_res_yes_sup <- readRDS(glue("{RDS_PATH}/treestruct_min_cl_size_res_yes_sup.rds")) 

# create matrix for each mcs and subtype + summarise log_vl (mean) to get only one per patient (treestructure phylotypes)
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
		vl_subtypes_comb3[[i,j]] <- inner_join(vl_subtypes_comb1[[i,j]], vl_subtypes[[j]][[1]], by="patientindex")
		print("===")
	}
}
saveRDS(vl_subtypes_comb1, glue("{RDS_PATH}/vl_subtypes_comb1.rds"))
saveRDS(vl_subtypes_comb2, glue("{RDS_PATH}/vl_subtypes_comb2.rds"))

# Mean VL measurements before treatment (without excluding phylotypes with less than e.g. 10 measurements)
vl_subtypes_all <- list(vl_subtypes_comb3[[1,1]], vl_subtypes_comb3[[1,2]], vl_subtypes_comb3[[1,3]], vl_subtypes_comb3[[1,4]])
vl_subtypes_all <- bind_rows(vl_subtypes_all); print(nrow(vl_subtypes_all)); print(length(unique(vl_subtypes_all$patientindex))) #41185 rows, 12833 unique patients 
vl_measur <- vl_subtypes_all %>% add_count(patientindex, name="n_vl") %>% group_by(patientindex) %>% filter(row_number() >= (n())) #filter(n() == 1)
mean(vl_measur$n_vl) # 3.21
sd(vl_measur$n_vl) # 2.53
median(vl_measur$n_vl) # 2
#IQR(vl_measur$n_vl);
calc_iqr(vl_measur$n_vl) #1-5
table(vl_measur$n_vl); hist(vl_measur$n_vl)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
# 4699 2188 1383 1127  924  779  706  517  262  128   57   31   17    7    4    3    1 

# For table 1
tab1_vl_ss <- readRDS(glue("{RDS_PATH}/vl_subtypes_comb2.rds"))
tab1_vl_ss[[1,1]]$subtype <- tree_names[1]; tab1_vl_ss[[1,2]]$subtype <- tree_names[2] 
tab1_vl_ss[[1,3]]$subtype <- tree_names[3]; tab1_vl_ss[[1,4]]$subtype <- tree_names[4] 
tab1_vl_ss_all <- rbind( tab1_vl_ss[[1,1]], tab1_vl_ss[[1,2]], tab1_vl_ss[[1,3]], tab1_vl_ss[[1,4]] )
print(nrow(tab1_vl_ss_all))
tab1_vl_ss_all %>% group_by(subtype) %>% summarise(n=n())
# A1: 755 patients, CRF_02_AG: 935, C: 2356, B: 8789

### Model 1: Estimating group-level differences using a multi-level Bayesian model
# Huge thanks to Chris Wymant for the advice and code: 
# https://htmlpreview.github.io/?https%3A//github.com/ChrisHIV/teaching/blob/main/other_topics/Stan_example_hierarchical_model_false_positives.html
libs_load2 <- c("tidyverse","rstan","ggforce", "Rcpp", "bayesplot")
invisible( lapply(libs_load2, library, character.only=TRUE) )

options(mc.cores = parallel::detectCores()) # parallelise
rstan_options(auto_write = TRUE)            # avoid re-compiling stan code
theme_set(theme_classic())
set.seed(123946)

ITER <- 4000
CHAINS <- 2
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
	
	pl1 <- ggplot(vl_subtypes_c, aes(x=as.factor(cluster), y=mean_log_vl_pat)) +
		geom_violin() + geom_sina() + labs(x = "phylotype") + theme(axis.text=element_text(size=6),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	system(glue("mkdir -p {RESULTS_PATH}/05_vl_{partition_method_lbl}/"))
	ggsave(plot=pl1, filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_expl.jpg"), width=10, height=8, dpi=300)

	# The following Stan code (stored as a string in R) estimates those parameters of the data-generating process that we don’t
	# condition on. The likelihood it uses matches the one we used to generate the data

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
	posterior_samples <- rstan::extract(fit)
	y_rep <- posterior_samples$y_rep
	#print("by group")
	print(names(posterior_samples))
	# print("mean")
	# print(head(posterior_samples$y_mean))
	
	# compute the variance explained in viral load (by this rstan model), we need to calculate the proportion of the total variance in
	# y that is explained by the group effects and the overall population-level effects
	# 1. Extract posterior samples for y_sd, y_sd_group, and group_effects_unscaled
	y_sd <- posterior_samples$y_sd  # Residual SD
	y_sd_group <- posterior_samples$y_sd_group  # Group-level SD
	# 2. Compute variances
	var_group <- y_sd_group^2  # Variance explained by group effects
	var_residual <- y_sd^2     # Residual variance
	var_total <- var_group + var_residual  # Total variance
	# 3. calculate R-squared for each posterior sample
	r_squared <- var_group / var_total
	# 4. summarize R-squared (mean, median, CIs)
	mean_r2 <- mean(r_squared)
	median_r2 <- median(r_squared)
	ci_r2 <- quantile(r_squared, probs = c(0.025, 0.975))
	print(paste("Mean R-squared:", mean_r2))
	print(paste("Median R-squared:", median_r2))
	print("95% Credible Interval for R-squared:")
	df_r2 <- data.frame(median_r2=median_r2, lower_ci=unname(ci_r2[1]), upper_ci=unname(ci_r2[2]))
	print(df_r2)
	
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
	ggsave(filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_ppc_dens_overlay_bp_.jpg"), width=8, height=6, dpi=300, bg="white")
	ppc_stat(vl_subtypes_c$mean_log_vl_pat, y_rep, stat = "mean")
	ggsave(filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_ppc_stat_mean_bp_.jpg"), width=8, height=6, dpi=300, bg="white")
	ppc_stat(vl_subtypes_c$mean_log_vl_pat, y_rep, stat = "sd")
	ggsave(filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_ppc_stat_sd_bp_.jpg"), width=8, height=6, dpi=300, bg="white")
	
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
							freq_upper = summary_$coefficients[1,1] + confint_[2,2], 
							adj_r_squared = summary_$adj.r.squared)
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
	#saveRDS(df_fit, glue("{RDS_PATH}/df_fit_{partition_method_lbl}_{mcs_subtype_choice}.rds"))
 rm(df_fit)
 gc()
	#saveRDS(df_prior, glue("{RDS_PATH}/df_prior_{partition_method_lbl}_{mcs_subtype_choice}.rds"))
	rm(df_prior)
	gc()

	list(df_table=df_table, sig_higher=sig_higher, df_compare=df_compare, df_r2=df_r2)
}

res_vl_model_ts <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) { #length(min_cl_size_choices)
	for(j in 1:length(tree_names)) { #length(tree_names)
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		res_vl_model_ts[[i,j]] <- vl_bayes_model(vl_subtypes_comb2[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}"), "treestructure", model_compiled)
	}
}

saveRDS(res_vl_model_ts, glue("{RDS_PATH}/res_vl_model.rds"))

# Table S7: bayes vl model estimates
vl_tables <- readRDS( glue("{RDS_PATH}/res_vl_model.rds") )
vl_tables[[1,1]]$df_table$subtype <- tree_names[1]; vl_tables[[1,2]]$df_table$subtype <- tree_names[2] 
vl_tables[[1,3]]$df_table$subtype <- tree_names[3]; vl_tables[[1,4]]$df_table$subtype <- tree_names[4] 
vl_table <- rbind( vl_tables[[1,1]]$df_table, vl_tables[[1,2]]$df_table, vl_tables[[1,3]]$df_table, vl_tables[[1,4]]$df_table )
vl_table <- vl_table %>% dplyr::select( subtype, group, quantile_0.025, quantile_0.5, quantile_0.975, y_mean_pop, p  )
vl_table <- vl_table %>% arrange( p, desc(quantile_0.5) )
vl_table$quantile_0.025 <- round(vl_table$quantile_0.025,3)
vl_table$quantile_0.5 <- round(vl_table$quantile_0.5,3)
vl_table$quantile_0.975 <- round(vl_table$quantile_0.975,3)
vl_table$y_mean_pop <- round(vl_table$y_mean_pop,3)
vl_table$p <- signif(vl_table$p, digits = 2)
options(scipen=999)
write.csv( vl_table, file=glue("{RESULTS_PATH}/tables/tableS7.csv"), quote=F, row.names = F )

# Fig. 2A: estimates with CIs for significantly lower and higher VLs
subtype_b_vois_ids_only <- c(40,69,133)
vl_table_B <- vl_tables[[1,4]]$df_table
mean_vl_pop <- unique(vl_table_B$y_mean_pop)
vl_table_B <- vl_table_B %>% dplyr::select( subtype, group, quantile_0.025, quantile_0.5, quantile_0.975, p )
backbone_cl_control <- readRDS(glue("{RDS_PATH}/backbone_cl_control.rds"))
vl_table_B_sel <- vl_table_B %>% filter( group == backbone_cl_control[[1,4]] | p <= 0.05 | group %in% subtype_b_vois_ids_only ) # last condition to include PT133
vl_table_B_sel$phylotype <- paste0("PT.B.",vl_table_B_sel$group,".UK")
vl_table_B_sel$phylotype[vl_table_B_sel$phylotype == "PT.B.153.UK"] <- "Backbone"
custom_order <- c("Backbone", "PT.B.8.UK", "PT.B.12.UK", "PT.B.20.UK", 
																		"PT.B.27.UK", "PT.B.40.UK", "PT.B.69.UK", "PT.B.83.UK", "PT.B.101.UK", "PT.B.133.UK")
ordered_vl_sig <- vl_table_B_sel$phylotype[match(custom_order, vl_table_B_sel$phylotype)]
vl_table_B_sel$phylotype <- factor(vl_table_B_sel$phylotype, levels = ordered_vl_sig)

vl_pal_vois <- c("PT.B.40.UK"="#56B4E9", "PT.B.69.UK"="#D55E00", "PT.B.133.UK"="#009E73", "Backbone"="#999999", "Other PTs"="grey25")

vl_table_B_sel$phylotype_display <- ifelse( vl_table_B_sel$phylotype %in% names(vl_pal_vois), as.character(vl_table_B_sel$phylotype), "Other PTs")

f2a_vl <- ggplot(data=vl_table_B_sel, aes(color=phylotype_display)) +
	#geom_errorbar(aes(x = phylotype, ymin = quantile_0.025, ymax = quantile_0.975), col = "blue") +
	geom_hline(yintercept = mean_vl_pop, linetype = "dashed") +
	geom_pointrange(aes(x=phylotype, y=quantile_0.5, ymin = quantile_0.025, ymax = quantile_0.975), size = 0.25) + #color = "grey25"
	scale_color_manual(values=vl_pal_vois, name="Phylotype") + #+ scale_fill_manual(values=ne_pt_pal, name="Phylotype") +
	#scale_y_continuous(limits=c(4.2,5.2)) +
	labs(x = "Phylotype", y=bquote(bold("Viral load estimate and 95% CIs ("~log[10]~"copies / mL)"))) + theme_classic() + 
	theme(axis.text.x = element_text(angle = 45,hjust = 1, size=8), axis.text=element_text(size=8)) + leg #legend.position = c(0.201,0.75)
saveRDS(f2a_vl, glue("{RDS_PATH}/f2a_vl.rds"))

### same model for fastbaps partitions (only doing for subtype B as none other PT signif in VL analysis in other subtypes) ###
fastbaps_df_list <- readRDS(glue("{RDS_PATH}/fastbaps_df_list.rds"))
vl_subtypes_comb_fb <- vl_subtypes_comb_fb2 <- vl_subtypes_mean_p_pat_fb <- list()
res_vl_model_fb <- list()
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

### same model for treecluster partitions gd=0.075, method=avg_clade (only doing for subtype B as none other PT signif in VL analysis in other subtypes) ###
res_treecluster <- glue("{RESULTS_PATH}/03_treecluster/B/")
vl_b_tc <- read.table(glue("{res_treecluster}/B_gd_0.075_avg_clade.log"), header=T) #24095
vl_b_tc <- vl_b_tc[vl_b_tc$ClusterNumber != -1,] # singletons removed
print(nrow(vl_subtypes[[4]][[2]])) #13349
vl_b_comb_tc <- inner_join(vl_b_tc, vl_subtypes[[4]][[2]], by=c("SequenceName"="testindex"))
print(nrow(vl_b_comb_tc)) # 8809
vl_b_mean_p_pat_tc <- vl_subtypes[[4]][[1]] %>% group_by(patientindex) %>% summarise(mean_log_vl_pat=mean(log_vl))
print(nrow(vl_b_mean_p_pat_tc)) #10903
vl_b_comb_tc2 <- inner_join(vl_b_comb_tc, vl_b_mean_p_pat_tc, by="patientindex")
print(nrow(vl_b_comb_tc2)) #8809
colnames(vl_b_comb_tc2) <- c("testindex","cluster","patientindex", "mean_log_vl_pat")

# Run VL model (subtype B only)
res_vl_model_tc <- vl_bayes_model(vl_b_comb_tc2, "B", "treecluster", model_compiled)
saveRDS(res_vl_model_tc, glue("{RDS_PATH}/res_vl_model_tc.rds"))
res_vl_model_tc <- readRDS(glue("{RDS_PATH}/res_vl_model_tc.rds"))
higher_vl_B_tc <- get_tips_cluster_id(vl_b_comb_tc2, "testindex", res_vl_model_tc$sig_higher) # 5 clusters, 568 sequences
inner_join_vl_vois_ts_tc <- inner_join(higher_vl_B_ts, higher_vl_B_tc, by="tips") # 153 matches, all from ts
# n=20 PT101 matches PT127 (+86 in tc), n=38 PT40 matches PT113 (+131 in tc), n=26 PT69 matches PT8 (+16 in tc), 
# n=69 PT20 matches PT20 (+48 in tc); one phylotype with >VL from tc not from ts: PT57 (142 seqs)
# trestruct VOIs estimates   (PT101, PT20, PT40, PT69, x):    4.90, 4.81, 4.91, 4.86, x
# treecluster VOIs estimates (PT127, PT20, PT113, PT8, PT57): 4.72, 4.73, 4.83, 4.77, 4.70
res_vl_model_tc <- readRDS(glue("{RDS_PATH}/res_vl_model_tc.rds"))

# plot correlation of bayes and lm model mean estimates for partitioning methods separately
plot_corr_vl_models_bayes_lm <- function(estims_df, mcs_subtype_choice, partition_method_lbl) {
	corr <- cor(estims_df$freq_estimate, estims_df$quantile_0.5)
	pl <- ggplot(estims_df, aes(x = freq_estimate, y = quantile_0.5)) +
		#geom_segment(aes(x = y_true, y = y_true, xend = quantile_0.5, yend = freq_estimate)) +
		geom_point(aes(x = quantile_0.5, y = freq_estimate)) +
		geom_abline(linetype = "dashed") +
		geom_smooth(method = "lm", color = "blue", se = FALSE) +
		annotate("text", x = max(estims_df$freq_estimate) - 1, y = max(estims_df$quantile_0.5), 
											label = paste("r =", round(corr, 2)), 
											size = 5, color = "black") +  # Annotate with correlation coefficient
		theme_classic()
	ggsave(plot=pl, filename=glue("{RESULTS_PATH}/05_vl_{partition_method_lbl}/{mcs_subtype_choice}_corr_bayes_lm.jpg"), width=10, height=7, dpi=300)
}

# treestructure mcs=30, subtype B
plot_corr_vl_models_bayes_lm(res_vl_model_ts[[1,4]]$df_compare, "30-B", "treestructure")
# fastbaps optimise.symmetric, subtype B
plot_corr_vl_models_bayes_lm(res_vl_model_fb[[4]]$df_compare, "B", "fastbaps")
# treecluster gd=0.075, method=avg_clade, subtype B
plot_corr_vl_models_bayes_lm(res_vl_model_tc$df_compare, "B", "treecluster")

# variance explained in different partition methods
res_vl_model_ts <- readRDS("rds/res_vl_model.rds")
res_vl_model_fb <- readRDS("rds/res_vl_model_fb.rds")
res_vl_model_tc <- readRDS("rds/res_vl_model_tc.rds")
res_vl_model_ts[[1,4]]$df_r2$method <- "treestructure (min. clade size = 30)"
res_vl_model_ts[[2,4]]$df_r2$method <- "treestructure (min. clade size = 50)"
res_vl_model_ts[[3,4]]$df_r2$method <- "treestructure (min. clade size = 100)"
res_vl_model_ts_all <- rbind(res_vl_model_ts[[1,4]]$df_r2, res_vl_model_ts[[2,4]]$df_r2, res_vl_model_ts[[3,4]]$df_r2)
res_vl_model_fb[[4]]$df_r2$method <- "fastbaps"
res_vl_model_tc$df_r2$method <- "treecluster"
comb_var_expl <- rbind( res_vl_model_ts_all,  res_vl_model_fb[[4]]$df_r2, res_vl_model_tc$df_r2)
View(comb_var_expl)
ggplot(comb_var_expl) + geom_errorbar(aes(x = method, ymin = lower_ci, ymax = upper_ci)) + geom_point(aes(x = method, y = median_r2)) + 
	labs(x="Method", y = expression(R^2)) + theme_classic() + theme(axis.text=element_text(family=helv, color="black"), axis.title=element_text(family=helv, color="black"))

# variance explained by lm models (adj_r_squared)
res_vl_model_ts[[1,4]]$df_compare$method <- "treestructure (min. clade size = 30)"
res_vl_model_ts[[2,4]]$df_compare$method <- "treestructure (min. clade size = 50)"
res_vl_model_ts[[3,4]]$df_compare$method <- "treestructure (min. clade size = 100)"
res_vl_model_ts_all2 <- rbind(res_vl_model_ts[[1,4]]$df_compare, res_vl_model_ts[[2,4]]$df_compare, res_vl_model_ts[[3,4]]$df_compare)
res_vl_model_fb[[4]]$df_compare$method <- "fastbaps"
res_vl_model_tc$df_compare$method <- "treecluster"

comb_var_expl_lm <- rbind( res_vl_model_ts_all2,  res_vl_model_fb[[4]]$df_compare, res_vl_model_tc$df_compare)
#View(comb_var_expl_lm)
ggplot(comb_var_expl_lm) + geom_point(aes(x = method, y = adj_r_squared)) + theme_classic() #geom_errorbar(aes(x = method, ymin = lower_ci, ymax = upper_ci))

# overlay them
comb_var_expl$model <- "Bayesian joint model"
comb_var_expl_lm$model <- "Linear model"
combined_data <- bind_rows(
	comb_var_expl %>% rename(y = median_r2),
	comb_var_expl_lm %>% rename(y = adj_r_squared)
)
combined_data$method <- factor(combined_data$method, levels = c("treestructure (min. clade size = 30)", "treestructure (min. clade size = 50)", "treestructure (min. clade size = 100)", "fastbaps", "treecluster"))
pal_bayes_lm_var_expl <- c("Bayesian joint model"="#999999", "Linear model"="black")
plot_var_expl_vl <- ggplot(combined_data) +
	geom_errorbar(data = subset(combined_data, model == "Bayesian joint model"), aes(x = method, ymin = lower_ci, ymax = upper_ci, color = model),width = 0.2) +
	geom_point(aes(x = method, y = y, color = model), size = 1) + scale_color_manual(values=pal_bayes_lm_var_expl, name="Viral load model") +
	labs(x = "Partitioning method", y = bquote("Variance explained (" ~ R^2 ~ ")") ) + #color = "Viral load model"
	theme_classic() +
	theme(axis.text = element_text(family = helv, color = "black"),axis.title = element_text(family = helv, color = "black"),
							axis.text.x = element_text(angle = 60, hjust = 1))
	
# Figure S4
ggsave(plot = plot_var_expl_vl, filename = glue("{RESULTS_PATH}/figs/figS4.jpg"), dpi = 600, width = 8, height = 6, bg = "white")
ggsave(plot = plot_var_expl_vl, filename = glue("{RESULTS_PATH}/figs/figS4.eps"), dpi = 600, width = 8, height = 6, bg = "white")

### Model 2: Welch's T-test with confints based on lm and adjusted for multiple testing

# vl_subtypes_comp_bef_ttest_clusts have clusters and vl_subtypes_comp_bef_ttest_ref have backbone mean_vls

vl_subtypes_comp_bef_ttest_clusts <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
vl_subtypes_comp_bef_ttest_ref <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) { #length(min_cl_size_choices)
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		vl_subtypes_comp_bef_ttest_clusts[[i,j]] <- split(vl_subtypes_comb2[[i,j]], vl_subtypes_comb2[[i,j]]$cluster)
		for(k in 1:length(vl_subtypes_comp_bef_ttest_clusts[[i,j]])) {
			if( (k != as.integer(backbone_cl_control[[i,j]])) & (!is.null(extracted_clusters[[i,j]][[1]][[k]]$tip.label)) ) {
				vl_subtypes_comp_bef_ttest_clusts[[i,j]][[k]] <- vl_subtypes_comb2[[i,j]][vl_subtypes_comb2[[i,j]]$cluster == k,]
			} else {
				vl_subtypes_comp_bef_ttest_clusts[[i,j]][[k]] <- NULL
			}
		}
		vl_subtypes_comp_bef_ttest_ref[[i,j]][[1]] <- vl_subtypes_comb2[[i,j]][vl_subtypes_comb2[[i,j]]$cluster == as.integer( backbone_cl_control[[i,j]] ),]
	}
}


t_test_diff_means_vl <- function(vl_subtypes_c, vl_subtypes_ref, extr_clusters, backbone_control, subtype_choice) {
	ttest <- ttest_res <- lm_vl <- confints <-  list() #sds
	for(k in 1:length(extr_clusters[[1]])) {
		print(k)
		if(length(unique(vl_subtypes_c[[k]]$patientindex)) < 10) {
			ttest[[k]] <- NULL
			print("Less than 10 patients, not computing!")
		} else {
			ttest[[k]] <- t.test(vl_subtypes_c[[k]]$mean_log_vl_pat, vl_subtypes_ref[[1]]$mean_log_vl_pat, alternative="greater")
			combined_clust_backb <- rbind(vl_subtypes_c[[k]], vl_subtypes_ref[[1]])
			ttest_res[[k]] <- data.frame(phylotype=k, rows_x=nrow(vl_subtypes_c[[k]]), estimate_x=round(ttest[[k]]$estimate[[1]],3), 
																																#lower_ci=(summ$coefficients[1,1] + confints[[k]][2,1]), upper_ci=(summ$coefficients[1,1] + confints[[k]][2,2]),
																																estimate_y=round(ttest[[k]]$estimate[[2]],3), p_value=ttest[[k]]$p.value, sd_x=round(sd(vl_subtypes_c[[k]]$mean_log_vl_pat),3), 
																																sd_y=round(sd(vl_subtypes_ref[[1]]$mean_log_vl_pat),3))
			ttest_res[[k]]$p_value <- signif(ttest_res[[k]]$p_value, digits = 2)
		}
		ttest_res_list <- rbindlist(ttest_res)
		
		# results un-adjusted for multiple testing
		system(glue("mkdir -p {RESULTS_PATH}/05_vl_ttest_treestructure/all_unadj/"))
		write.csv(ttest_res_list, file=glue("{RESULTS_PATH}/05_vl_ttest_treestructure/all_unadj/{subtype_choice}_vl.csv"), quote=F, row.names=F) # join all of those for mcs=30 for table S8
		
		ttest_res_list_sig <- ttest_res_list[ttest_res_list$p_value <= 0.05, ]
		system(glue("mkdir -p {RESULTS_PATH}/05_vl_ttest_treestructure/sig_unadj/"))
		write.csv(ttest_res_list_sig, file=glue("{RESULTS_PATH}/05_vl_ttest_treestructure/sig_unadj/{subtype_choice}_vl.csv"), quote=F, row.names=F)
		
		# results adjusted for multiple testing
		ttest_res_list_p <- ttest_res_list$p_value
		
		.adjust_methods <- function(ps, method) {
			res_adj <- p.adjust(ps, method = method)
			res_adj
		}
		
		bonf <- .adjust_methods(ttest_res_list_p, "bonferroni")
		fdrr <- .adjust_methods(ttest_res_list_p, "fdr")
		
		ttest_res_list_adj <- ttest_res_list
		ttest_res_list_adj$p_adj_bonf <- bonf
		ttest_res_list_adj$p_adj_fdr <- fdrr
		
		# all adjusted
		system(glue("mkdir -p {RESULTS_PATH}/05_vl_ttest_treestructure/all_adj/"))
		write.csv(ttest_res_list_adj, file=glue("{RESULTS_PATH}/05_vl_ttest_treestructure/all_adj/{subtype_choice}_vl.csv"), quote=F, row.names=F)
		
		# bonferroni corrected
		system(glue("mkdir -p {RESULTS_PATH}/05_vl_ttest_treestructure/sig_adj_bonf/"))
		ttest_res_list_adj_bonf <- ttest_res_list_adj[ttest_res_list_adj$p_adj_bonf <= 0.05,]
		write.csv(ttest_res_list_adj_bonf, file=glue("{RESULTS_PATH}/05_vl_ttest_treestructure/sig_adj_bonf/{subtype_choice}_vl.csv"), quote=F, row.names=F)
		
		# fdr corrected
		system(glue("mkdir -p {RESULTS_PATH}/05_vl_ttest_treestructure/sig_adj_fdr/"))
		ttest_res_list_adj_fdr <- ttest_res_list_adj[ttest_res_list_adj$p_adj_fdr <= 0.05,]
		write.csv(ttest_res_list_adj_fdr, file=glue("{RESULTS_PATH}/05_vl_ttest_treestructure/sig_adj_fdr/{subtype_choice}_vl.csv"), quote=F, row.names=F)
		
		#View(ttest_res_list)
	}
	list(all_adj=ttest_res_list_adj, sig_unadj=ttest_res_list_sig, sig_adj_bonf=ttest_res_list_adj_bonf, sig_adj_fdr=ttest_res_list_adj_fdr, ttest_res_list=ttest_res_list)
}

rm_paraphyletic_pt_alns <- readRDS("rds/rm_paraphyletic_pt_alns.rds")
extracted_clusters <- readRDS("rds/extracted_clusters.rds")
vl_st_test <- vl_vois_incl_parap <- vl_vois_excl_parap <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(subtype_choices))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		vl_st_test[[i,j]] <- t_test_diff_means_vl(vl_subtypes_comp_bef_ttest_clusts[[i,j]], vl_subtypes_comp_bef_ttest_ref[[i,j]], extracted_clusters[[i,j]], backbone_cl_control[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		# matrix with vl VOIs (incl parap)
		#vl_vois_incl_parap[[i,j]] <- as.integer(unique(vl_st_test[[i,j]]$phylotype))
		# matrix with vl VOIs (excl parap)
		#vl_vois_excl_parap[[i,j]] <- setdiff(vl_vois_incl_parap[[i,j]], as.integer(rm_paraphyletic_pt_alns[[i,j]]) )
	}
}
sum(vl_st_test[[1,4]]$ttest_res_list$rows_x) # number of patients in t-test analysis (after removing PTs with <10 patients)
vl_st_test[[1,4]]$sig_unadj$phylotype #|> length() #29
vl_st_test[[1,4]]$sig_adj_bonf$phylotype #PTs 20,40 (overlap bayes ones, but 2 missing [69 and 101])
vl_st_test[[1,4]]$sig_adj_fdr$phylotype #|> length() #7 PTs: 4,20,40,69,79,125,137 (overlap with 20,40 from bayes; 69 from before and bayes; 137 from before; 4, 79 and 125 new)
vl_st_test[[1,4]]$all_adj$phylotype #|> length() #115
saveRDS(vl_st_test, glue("{RDS_PATH}/vl_st_test.rds"))

# Table S6: T-test on viral loads corrected by multiple testing
vl_tables_ttest <- readRDS(glue("{RDS_PATH}/vl_st_test.rds"))
vl_tables_ttest[[1,1]]$all_adj$subtype <- tree_names[1]; vl_tables_ttest[[1,2]]$all_adj$subtype <- tree_names[2] 
vl_tables_ttest[[1,3]]$all_adj$subtype <- tree_names[3]; vl_tables_ttest[[1,4]]$all_adj$subtype <- tree_names[4] 
vl_table_ttest <- rbind( vl_tables_ttest[[1,1]]$all_adj, vl_tables_ttest[[1,2]]$all_adj, vl_tables_ttest[[1,3]]$all_adj, vl_tables_ttest[[1,4]]$all_adj )
vl_table_ttest <- vl_table_ttest %>% dplyr::select( subtype, phylotype, rows_x, estimate_x, sd_x, estimate_y, sd_y, p_value, p_adj_bonf, p_adj_fdr )
vl_table_ttest <- vl_table_ttest %>% arrange( p_adj_fdr, desc(estimate_x) )
vl_table_ttest$p_value <- signif(vl_table_ttest$p_value, digits = 2)
vl_table_ttest$p_adj_fdr <- signif(vl_table_ttest$p_adj_fdr, digits = 2)
vl_table_ttest$p_adj_bonf <- signif(vl_table_ttest$p_adj_bonf, digits = 2)
options(scipen=999)
write.csv( vl_table_ttest, file=glue("{RESULTS_PATH}/tables/tableS6.csv"), quote=F, row.names = F )