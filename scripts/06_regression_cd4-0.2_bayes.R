libs_load <- c("ggplot2","dplyr","ggpubr","data.table","glue","lubridate","sjPlot","viridis", "rstan", "brms","bayesplot",  "ape")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"
system(glue("mkdir -p {RESULTS_PATH}/06_cd4_bayes_model"))

# Bayesian CD4 model
residuals_removal_mx <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))

fit_bayes_cd4_model <- function(outcome, cd4_df, mcs_subtype_choice_cd4_transf) {
	
	formula <- as.formula(paste(outcome, "~ years_since_1cd4 + age_group + 
    years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
    (1 + years_since_1cd4 | patientindex) + (1 + years_since_1cd4 | phylotype)"))
	
	# TODO later extract and modify priors
	# gp <- get_prior(formula, data=cd4_df, warmup = WARMUP, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
	# print("priors:")
	# print(gp)
	#gp$prior[ grepl(gp$coef, pattern = 'years_since_1cd4:phylotype') ] <- glue("normal(0, {se_cd4_prior[[i,j]]})")
	#print(gp)
	
	# TODO add prior=gp after data below when changing priors
	bayes_fit <- brm(formula, data=cd4_df, warmup = WARMUP, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
	#saveRDS(bayes_fit, glue("{RDS_PATH}/bayes_fit_{mcs_subtype_choice_cd4_transf}.rds"))
	
	# Extract diagnostics from the Stan fit object
	stan_summary <- rstan::summary(bayes_fit$fit)$summary
	
	# R-hat
	rhat_values <- stan_summary[, "Rhat"]
	# Effective sample sizes
	n_eff_values <- stan_summary[, "n_eff"]
	
	# Combine
	diagnostics <- data.frame(
		Parameter = rownames(stan_summary),
		Rhat = rhat_values,
		n_eff = n_eff_values
	)
	View(diagnostics)
	
	print("Any Rhat > 1.01?")
	print(any(diagnostics$Rhat > 1.01))
	print("Any n_eff < 400?")
	print(any(diagnostics$n_eff < 400))
	
	# First approach: random effects only
	
	print("Extracting phylotype coefs")
	pattern_pts <- "^r_phylotype\\[\\d+,years_since_1cd4\\]$"
	phylotype_coefs <- stan_summary[grepl(pattern_pts, rownames(stan_summary)), ]
	phylotype_coefs <- as.data.frame(phylotype_coefs)
	phylotype_coefs$phylotype <- rownames(phylotype_coefs)
	phylotype_coefs$phylotype <- as.factor(gsub(".*\\[(\\d+),.*", "\\1", phylotype_coefs$phylotype))
	
	phylotype_coefs$significance <- ifelse(phylotype_coefs$`2.5%` < 0 & phylotype_coefs$`97.5%` < 0, yes="*", no="")
	
	print("Extracting significant phylotypes")
	cd4_signifs <- phylotype_coefs[ phylotype_coefs$significance == "*", ]
	cd4_signifs <- cd4_signifs[ order(cd4_signifs$mean), ] #`50%`
	write.table(cd4_signifs, glue("results/06_cd4_bayes_model/{mcs_subtype_choice_cd4_transf}_signif_phylotype_slopes.tsv"), sep="\t", quote=F, row.names=F)
	print(cd4_signifs)
	
	print("Plotting phylotype slopes with 95% CIs")
	p1 <- ggplot(data=phylotype_coefs) +
		geom_errorbar(aes(x = phylotype, ymin = `2.5%`, ymax = `97.5%`), col = "blue") +
		geom_point(aes(x = phylotype, y = mean), col = "blue") + #`50%`
		geom_text(aes(x=phylotype, label = significance, y = `97.5%` + 0.2), size = 10, color = "blue", fontface = "bold") +
		labs(x = "phylotype", y="estimate") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	ggsave(plot=p1, filename=glue("{RESULTS_PATH}/06_cd4_bayes_model/{mcs_subtype_choice_cd4_transf}_phylotype_slopes.jpg"), width=15, height=8, dpi=300)
	
	print("Extract 'b_' and 'r_phylotype' parameters to plot posteriors")
	phylotype_coefs_poster <- stan_summary[grepl("^b_", rownames(stan_summary)) | grepl(pattern_pts, rownames(stan_summary)), ]
	phylotype_coefs_poster <- as.data.frame(phylotype_coefs_poster)
	phylotype_coefs_poster$parameter <- rownames(phylotype_coefs_poster)
	posterior_array <- as.array(bayes_fit)
	p2 <- mcmc_areas(posterior_array,pars = rownames(phylotype_coefs_poster), prob = 0.95 ) + theme_minimal() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA))
	print("Plotting posterior of parameters")
	ggsave(plot=p2, filename=glue("{RESULTS_PATH}/06_cd4_bayes_model/{mcs_subtype_choice_cd4_transf}_posterior_params.jpg"), width=8, height=15, dpi=300)
	
	phylotype_coefs_poster$significance <- ifelse(phylotype_coefs_poster$`2.5%` < 0 & phylotype_coefs_poster$`97.5%` < 0, yes="*", no="")
	print("Plotting fixed effects")
	p3 <- ggplot(data=phylotype_coefs_poster) +
		geom_errorbar(aes(x = parameter, ymin = `2.5%`, ymax = `97.5%`), col = "blue") +
		geom_point(aes(x = parameter, y = mean), col = "blue") + #`50%`
		geom_text(aes(x=parameter, label = significance, y = `97.5%` + 0.2), size = 10, color = "blue", fontface = "bold") +
		labs(x = "parameter", y="estimate") + theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	ggsave(plot=p3, filename=glue("{RESULTS_PATH}/06_cd4_bayes_model/{mcs_subtype_choice_cd4_transf}_fixed_effects.jpg"), width=15, height=8, dpi=300)
	
	print("Extracting significant parameters among all fixed and random effects")
	phylotype_coefs_poster2 <- phylotype_coefs_poster[ phylotype_coefs_poster$significance == "*", ]
	phylotype_coefs_poster2 <- phylotype_coefs_poster2[ order(phylotype_coefs_poster2$mean), ] # `50%`
	write.table(phylotype_coefs_poster2, glue("results/06_cd4_bayes_model/{mcs_subtype_choice_cd4_transf}_signif_among_all_fixed_rand_effs.tsv"), sep="\t", quote=F, row.names=F)
	print(phylotype_coefs_poster2)
	# print to screen only the signif (phylotype_coefs_poster2) but return to function all (phylotype_coefs_poster)
	
	list(model_summary=stan_summary, phylotype_coefs=phylotype_coefs, cd4_signifs=cd4_signifs, phylotype_coefs_poster=phylotype_coefs_poster) #random_effects=random_effects
}
# for tests 10, 100, 1, 1
# for first try: 1000, 10000, 10, 4
WARMUP <- 1000 #2500
ITER <- 10000 #25000
THIN <- 10 #25
CORES_CHAINS <- 4

# call bayes CD4 model
#bmm <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
#bmm_sqrt <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
#bmm4 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
bmm3 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:1) { #length(min_cl_size_choices)
	for(j in 4:length(tree_names)) { #1:length(tree_names)
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		if(nrow(residuals_removal_mx[[i,j]]$cd4_df_filt) > 0) {
			#start_time <- Sys.time()
			# CD4 only basic filter
			# bmm[[i,j]] <- fit_bayes_cd4_model("cd4", residuals_removal_mx[[i,j]]$cd4_df_filt, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4"))
			# end_time <- Sys.time(); time_taken <- round(end_time - start_time,2)
			# print("Time to fit model and return results")
			# print(time_taken)
			
			# CD4 only basic filter
			# bmm[[i,j]] <- fit_bayes_cd4_model("cd4", residuals_removal_mx[[i,j]]$cd4_df_filt, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4"))
			# end_time <- Sys.time(); time_taken <- round(end_time - start_time,2)
			# print("Time to fit model and return results")
			# print(time_taken)
			
			# TODO try only with removal of short follow-up patients later?
			
			# CD4 only advanced filter (short follow-ups and >90% quantile slopes removed)
			start_time <- Sys.time()
			#bmm4[[i,j]] <- fit_bayes_cd4_model("cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt2, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4-filt_short_high_slopes"))
			bmm3[[i,j]] <- fit_bayes_cd4_model("cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt1, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4-filt_short"))
			end_time <- Sys.time(); time_taken <- round(end_time - start_time,2)
			print("Time to fit model and return results")
			print(time_taken)
			
			# sqrt CD4 only basic filter
			# start_time <- Sys.time()
			# bmm_sqrt[[i,j]] <- fit_bayes_cd4_model("sqrt_cd4", residuals_removal_mx[[i,j]]$cd4_df_filt, glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4"))
			# print("Time to fit model and return results")
			# end_time <- Sys.time(); time_taken <- round(end_time - start_time,2)
			# print(time_taken)
			
			# TODO get tables and plot for combined random phylotype + fixed effects & random phylotype + fixed years_s1cd4
			# When to report each:
			# 1) Only random effects: Use to highlight relative differences between phylotypes compared to the population average. When the main interest is in variability across phylotypes, rather than absolute effects.
			# 2) Combined random + fixed: These represent the absolute effect for each phylotype, accounting for both global and specific contributions. Use when phylotype-specific absolute slopes or intercepts are of primary interest.
			# 3) Comparison of Phylotype Slopes: If your primary focus is on how slopes differ between phylotypes (e.g., how years_since_1cd4 impacts CD4 for each phylotype)
			# Use the slope term r_phylotype[, years_since_1cd4] combined with the fixed effect (b_years_since_1cd4).
			# Highlight phylotypes with particularly positive or negative slopes.
		}
	}
}
#saveRDS(bmm_sqrt, glue("{RDS_PATH}/bmm_sqrt.rds"))
#saveRDS(bmm, glue("{RDS_PATH}/bmm.rds"))
#saveRDS(bmm4, glue("{RDS_PATH}/bmm4.rds"))
saveRDS(bmm3, glue("{RDS_PATH}/bmm3.rds"))

# TODO print summary of demog variables for PTs with both higher SPVL and CD4 counts?
# TODO: if verifying if paraphyletic VOIs have monophyletic VOIs within it: is any internal node of A the MRCA of B? then copy code from bkp of bayes_model
# TODO to plot average trajectory of PTs, get previous code from bkp of bayes_model