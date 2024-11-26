libs_load <- c("ggplot2","dplyr","ggpubr","data.table","glue","lme4","lubridate","sjPlot","viridis", "boot")
invisible( lapply(libs_load, library, character.only=TRUE) )

RDS_PATH="rds"
RESULTS_PATH="results"
DATA_PATH="data"
system(glue("mkdir -p {RESULTS_PATH}/06_cd4_ml_model"))

min_cl_size_choices <- c(30, 50, 100)
tree_names <- c("A_A1","CRF_02_AG","C","B")

treestruct_min_cl_size_res_yes_sup <- readRDS(glue("{RDS_PATH}/treestruct_min_cl_size_res_yes_sup.rds")) 
cd4_before_art_subtypes <- readRDS(glue("{RDS_PATH}/cd4_before_art_subtypes.rds"))

# include phylotypes in dfs
cd4_subtypes_comb1 <- cd4_subtypes_comb2 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names)) #cd4_subtypes_mean_p_pat <-
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		print("treestructure matching")
		print(nrow(treestruct_min_cl_size_res_yes_sup[[i,j]][[2]]))
		print(nrow(cd4_before_art_subtypes[[j]][[1]]))
		# first match seq ids
		cd4_subtypes_comb1[[i,j]] <- inner_join(treestruct_min_cl_size_res_yes_sup[[i,j]][[2]], cd4_before_art_subtypes[[j]][[2]], by=c("taxon"="testindex"))
		print("nrow and unique patients cd4_subtypes_comb1")
		print(nrow(cd4_subtypes_comb1[[i,j]]))
		print(length(unique(cd4_subtypes_comb1[[i,j]]$patientindex)))
		# then match back to get all CD4s (multiple per patient)
		cd4_subtypes_comb2[[i,j]] <- inner_join(cd4_subtypes_comb1[[i,j]], cd4_before_art_subtypes[[j]][[1]], by="patientindex")
		# get measurement ids ordering by cd4_date per patient
		cd4_subtypes_comb2[[i,j]] <- cd4_subtypes_comb2[[i,j]] %>% group_by(patientindex) %>% 
			arrange(cd4_decimal_date, .by_group = TRUE) %>% 
			mutate(measurement_id = row_number()) %>% ungroup()
		# add years since first CD4
		cd4_subtypes_comb2[[i,j]] <- cd4_subtypes_comb2[[i,j]] %>% group_by(patientindex) %>% mutate(years_since_1cd4 = round(cd4_decimal_date - cd4_decimal_date[measurement_id==1], 3) )
		# add age_group
		cd4_subtypes_comb2[[i,j]]$hiv_diag_decimal_date <- decimal_date(as.Date(cd4_subtypes_comb2[[i,j]]$hivpos_ymd))
		cd4_subtypes_comb2[[i,j]] <- cd4_subtypes_comb2[[i,j]] %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y))
		cd4_subtypes_comb2[[i,j]] <- cd4_subtypes_comb2[[i,j]] %>% mutate(
			age_group = dplyr::case_when(age_diag<=29 ~ "<29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
			age_group = factor(age_group,level = c("<29","30-39","40-49","50-59","60+")))
		# phylotype
		cd4_subtypes_comb2[[i,j]]$phylotype <- cd4_subtypes_comb2[[i,j]]$cluster
		print("nrow and unique patients cd4_subtypes_comb2")
		print(nrow(cd4_subtypes_comb2[[i,j]]))
		print(length(unique(cd4_subtypes_comb2[[i,j]]$patientindex)))
		print("===")
	}
}

hist(cd4_subtypes_comb2[[1,4]]$years_since_1cd4, breaks=50)

backbone_cl_control <- readRDS(glue("{RDS_PATH}/backbone_cl_control.rds"))

# CD4 model for cd4 decline
# cd4 ~ years_since_1cd4 + age_group + years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
# 	(years_since_1cd4 | patientindex)  + (years_since_1cd4 | phylotype),
fit_cd4_model <- function(outcome, cd4_df, mcs_subtype_choice_cd4_transf, do_boot) {
	# distribution for the slope and intercept random effects as a two-dimensional joint normal, 
	# rather than two independent one-dimensional normals (did not change crude estims of variance explained)
	
	# previously: (years_since_1cd4 | patientindex) + (years_since_1cd4 | phylotype)
	
	formula <- as.formula(paste(outcome, "~ years_since_1cd4 + age_group + 
    years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
    (1 + years_since_1cd4 | patientindex) + (1 + years_since_1cd4 | phylotype)")) #(years_since_1cd4 | patientindex) + (years_since_1cd4 | phylotype)
	fit <- lmer(formula, data=cd4_df, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) # adding optim removed 'fail to converge' warning
	#print(summary( fit ))
	
	# extract random effects of phylotype
	randeff_pt <- ranef( fit )$phylotype
	colnames(randeff_pt) <- c("intercept_randeff", "years_since_1cd4_randeff")
	randeff_pt$phylotype <- rownames(randeff_pt)
	#View(randeff_pt)
	# extract combined fixed effects and random effects for phylotype
	comb_rand_fixed_pt1 <- coef( fit )$phylotype
	comb_rand_fixed_pt <- comb_rand_fixed_pt1 %>% dplyr::select("(Intercept)", "years_since_1cd4")
	# TODO bring back fixed effects later
	colnames(comb_rand_fixed_pt) <- c("intercept_fixedeff", "slope_phylotype")
	comb_rand_fixed_pt$phylotype <- rownames(comb_rand_fixed_pt)
	combined_coeffs <- inner_join(randeff_pt, comb_rand_fixed_pt, by="phylotype")
	#View(combined_coeffs)
	#View(summary(fit)$coefficients[,1])
	
	if(do_boot) {
		# function to bootstrapping random effects
		.random_effect_se <- function(data, indices) {
			data <- data %>% dplyr::select(patientindex, phylotype, cd4, sqrt_cd4, years_since_1cd4, age_group, sexid, exposureid)
			boot_data <- data[indices, ]
			#print(boot_data)
			formula <- as.formula("cd4 ~ years_since_1cd4 + age_group + 
	    years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
	    (1 + years_since_1cd4 | patientindex) + (1 + years_since_1cd4 | phylotype)")
			#boot_fit <- lmer(formula, data=boot_data, REML=FALSE, control = lmerControl(optimizer = "bobyqa"))
			
			# Fit model and handle failures
			boot_fit <- tryCatch(
				lmer(formula, data = boot_data, REML = FALSE, control = lmerControl(optimizer = "bobyqa")),
				error = function(e) return(rep(NA, length(unique(data$phylotype))))  # Return consistent NA vector if model fails
			)
			
			# Extract random effects for all phylotypes
			if (is.numeric(boot_fit)) {
				return(boot_fit)  # Return NA vector if model fitting failed
			}
			
			ranef_vals <- coef(boot_fit)$phylotype %>% dplyr::select(years_since_1cd4)
			
			# Ensure all phylotype levels are present
			all_phylotypes <- unique(data$phylotype)
			ranef_vals_vec <- sapply(all_phylotypes, function(ptype) {
				if (ptype %in% rownames(ranef_vals)) {
					ranef_vals[ptype, "years_since_1cd4"]
				} else {
					NA  # Fill missing levels with NA
				}
			})
			
			#print(length(ranef_vals_vec))
			print("rep")
			ranef_vals_vec
		}
		# call bootstrapping
		
		NREP <- 100
		NCPU <- 4
		
		print("Running boot")
		boot_res <- boot::boot(data=cd4_df, statistic=.random_effect_se, R = NREP, ncpus=NCPU)
		print(boot_res)
		# Calculate standard errors for random effects
		random_effect_ses <- apply(boot_res$t, 2, sd)
		
		# Calculate t-values for random effects
		comb_rand_fixed_pt2 <- comb_rand_fixed_pt1 %>% dplyr::select(years_since_1cd4)
		colnames(comb_rand_fixed_pt2) <- c("phylotype_coef")
		comb_rand_fixed_pt2$phylotype <- rownames(comb_rand_fixed_pt1)
		comb_rand_fixed_pt2 <- comb_rand_fixed_pt2 %>% dplyr::select(phylotype, phylotype_coef)
		comb_rand_fixed_pt2$t_value <- comb_rand_fixed_pt2[, "phylotype_coef"] / random_effect_ses
		comb_rand_fixed_pt2$se <- random_effect_ses
		comb_rand_fixed_pt2$p_value <- 2 * (1 - pnorm(abs(comb_rand_fixed_pt2$t_value)))
		View(comb_rand_fixed_pt2)
		#comb_rand_fixed_pt2 <- comb_rand_fixed_pt2 %>% mutate(significance = if_else(p_value < 0.05, "*", ""))
		comb_rand_fixed_pt2 <- comb_rand_fixed_pt2[order(comb_rand_fixed_pt2$phylotype_coef), ]
		
		# Plot estimates with confidence intervals
		estim_ml_boot <- ggplot(comb_rand_fixed_pt2, aes(x = phylotype, y = phylotype_coef)) +
			geom_point() +
			geom_errorbar(aes(ymin = phylotype_coef - 1.96 * se, ymax = phylotype_coef + 1.96 * se)) +
			#geom_text(aes(x=phylotype, label = significance, y = (phylotype_coef + 1.96 * se) + 0.2), size = 10, color = "black", fontface = "bold") +
			theme(axis.text=element_text(size=5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
			labs(title = "Phylotype Effects on CD4 Decline", y = "Effect Estimate")
		ggsave(plot=estim_ml_boot, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice_cd4_transf}_boot_estim.jpg"), width=15, height=8, dpi=300)
		return(	list(fit=fit, combined_coeffs=combined_coeffs, table_pt_estims_p=comb_rand_fixed_pt2))
	} else {
		return(list(fit=fit, combined_coeffs=combined_coeffs))
	}
}

lmm1 <- lmm_sqrt1 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names)) 
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
			if(nrow(cd4_subtypes_comb2[[i,j]]) > 0) {
				print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
				# define as ref category for fixed effects: backbone phylotype, age 30-39, sexid Male, exposureid Homo-bisexual
				# untransformed cd4 model
				cd4_subtypes_comb2[[i,j]] <- cd4_subtypes_comb2[[i,j]] %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[i,j]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
				lmm1[[i,j]] <- fit_cd4_model("cd4", cd4_subtypes_comb2[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4"), do_boot=F)
				# sqrt model
				lmm_sqrt1[[i,j]] <- fit_cd4_model("sqrt_cd4", cd4_subtypes_comb2[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4"), do_boot=F)
			} else {
				next
			}
		}
}

saveRDS(cd4_subtypes_comb2, glue("{RDS_PATH}/cd4_subtypes_comb2.rds"))

# ranef() Returns only the random effects for phylotype (deviations from the fixed effects)
# coef() Combines the fixed effects with the random effects for phylotype. Provides the actual coefficients for each phylotype, including the specific intercept and slope for years_since_1cd4.
# For the random effect (years_since_1cd4 | phylotype), no specific reference group is needed because the random effects capture deviations for each phylotype from the overall population-level effects
# For fixed effect specified backbone, people in thirties, and male as reference
# Random Effects Covariance Matrix: 
# PATIENTINDEX -> variance: 45678.98, sd: 213.726 -> variation in intercepts across patients high
# slope variance: ~900, sd: ~30 also high, the rate of change in CD4 levels varies by 30 units per year across patients
# correlation random intercept vs slope: -0.54; patients with higher intercepts tend to have slower declines in CD4 levels over time

# PHYLOTYPE -> variance 592.41, sd 24.340 -> variation in intercepts lower than the one of paatientindex
# slope: variance 54.24, sd ~7.4 -> variation in slopes much smaller
# correlation: +0.49: higher intercepts within PT associated with faster decline in CD4 levels

# residual: variance and sd (21k and 145): still substantial variance not explained by fixed and random effects
# IMPORTANT: demonstrates heritability in a crude way

# Three "simple/crude" measures of the fraction of variance explained by the lmm.
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#how-do-i-compute-a-coefficient-of-determination-r2-or-an-analogue-for-glmms
r2.corr.mer <- function(m) {
	print("corr 1: ")
	lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
	print(summary(lmfit)$r.squared)
	print("corr 2")
	print( 1 - var(residuals(m)) / var(model.response(model.frame(m))) )
	print("corr 3")
	print( cor(model.response(model.frame(m)), predict(m, type = "response"))^2 )
}

r2_corrs <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names)) 
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		print("untrans CD4")
		r2.corr.mer(lmm1[[i,j]]$fit)
		print("sqrt CD4")
		r2.corr.mer(lmm_sqrt1[[i,j]]$fit)
	}
}

# for B30: # same when using random effects as indep normals or joint normals: cd4 ~0.78
# sqrt cd4: ~0.87

# Visualise and remove outlying individuals in a distribution of maximum-likehood slopes and intercepts + outlying measurements
vis_residuals <- function(fit, cd4_df, mcs_subtype_choice) {
	# Extract random effects
	random_effects <- ranef(fit)
	# For patientindex
	individual_effects_df <- as.data.frame(random_effects$patientindex)
	individual_effects_df$patientindex <- as.integer( rownames(individual_effects_df) )
	#print(individual_effects_df)
	
	# Scatterplot
	sc <- ggplot(individual_effects_df, aes(x = `(Intercept)`, y = `years_since_1cd4`)) +
		geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") +
		labs(x = "Random Intercept", y = "Random Slope", title = "Random Effects: Intercept vs. Slope")
	ggsave(plot=sc, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_scatter_residuals.jpg"), width=10, height=8, dpi=300)
	# Distr plots (intercepts and slopes)
	distr_interc <- ggplot(individual_effects_df, aes(x = `(Intercept)`)) +
		geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
		labs(x = "Random Intercept", title = "Distribution of Random Intercepts")
	ggsave(plot=distr_interc, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_distr_interc.jpg"), width=10, height=8, dpi=300)
	distr_slopes <- ggplot(individual_effects_df, aes(x = `years_since_1cd4`)) +
		geom_histogram(bins = 30, fill = "green", alpha = 0.7) +
		labs(x = "Random Slope", title = "Distribution of Random Slopes")
	ggsave(plot=distr_slopes, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_distr_slopes.jpg"), width=10, height=8, dpi=300)
	
	# Calculate z-scores for intercepts and slopes to identify outliers
	# `scale` standardises a vector or a matrix by subtracting the mean and dividing by the standard deviation
	zs <- individual_effects_df %>% mutate(intercept_z = scale(`(Intercept)`), slope_z = scale(`years_since_1cd4`))
	# if >3 the data point is quite different from the other data points
	outl_indivs <- zs %>% filter(abs(intercept_z) > 3 | abs(slope_z) > 3)
	print("Outliers detected (individuals):")
	print(nrow(outl_indivs))
	#View(outl_indivs)
	
	# Scatterplot again to see if looks better
	individual_effects_rm_outl_df <- individual_effects_df[!individual_effects_df$patientindex %in% outl_indivs$patientindex,]
	sc_outl_rm <- ggplot(individual_effects_rm_outl_df, aes(x = `(Intercept)`, y = `years_since_1cd4`)) +
		geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") +
		labs(x = "Random Intercept", y = "Random Slope", title = "Random Effects: Intercept vs. Slope")
	ggsave(plot=sc_outl_rm, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_scatter_residuals_outl_rm.jpg"), width=10, height=8, dpi=300)
	
	# usual CD4 counts up to 1800 cells/mm^3 (sqrt = 42.4); allowing up to sqrt = 45
	sc_meas <- ggplot(cd4_df, aes(x = `cd4_decimal_date`, y = `sqrt_cd4`)) +
		geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") +
		labs(x = "Date of CD4", y = "CD4 value")
	ggsave(plot=sc_meas, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_scatter_meas.jpg"), width=10, height=8, dpi=300)
	
	# remove outlier measurements
	outl_meas <- cd4_df[cd4_df$sqrt_cd4 > 45,]
	print("Outliers detected (measurements):")
	print(nrow(outl_meas))
	
	# remove outlier individuals and measurements to later fit model without them
	cd4_df_filt <- cd4_df[!cd4_df$patientindex %in% outl_indivs,]
	print("nrow after removing individuals: ")
	print(nrow(cd4_df_filt))
	cd4_df_filt <- cd4_df[cd4_df$sqrt_cd4 <= 45,]
	print("nrow after removing measurements: ")
	print(nrow(cd4_df_filt))
	cd4_df_filt <- cd4_df_filt[cd4_df_filt$patientindex %in% unique(cd4_df_filt$patientindex[duplicated(cd4_df_filt$patientindex)]),]
	print("nrow after double-checking if by removing a measurement a patient ended up with n<2: ")
	print(nrow(cd4_df_filt))
	
	list(outl_indivs=outl_indivs, outl_meas=outl_meas, cd4_df_filt=cd4_df_filt)
}

# Run removal of outlying individuals and measurements
vl_subtypes_comb2 <- readRDS(glue("{RDS_PATH}/vl_subtypes_comb2.rds"))
residuals_removal_mx <-  ij_outliers_pts <- ij_outliers_meas <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		residuals_removal_mx[[i,j]] <- vis_residuals(lmm1[[i,j]]$fit, cd4_subtypes_comb2[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}")) # 130 outlier indivs, 11 outlier measurements
		# inner_join of outlier individuals with cluster membership to see if removing VL VOI patients
		ij_outliers_pts[[i,j]] <- inner_join(vl_subtypes_comb2[[1,4]], residuals_removal_mx[[i,j]]$outl_indivs, by="patientindex")
		print("Patients with outlying slopes or intercepts")
		print(length(unique(ij_outliers_pts[[i,j]]$patientindex))) # 98 patients with outlying slopes or intercepts
		print("Distribution of outliers (slope and intercept) in phylotypes")
		print( table(ij_outliers_pts[[i,j]]$cluster) ) # none of subtype B VL VOIs (PTs 20, 40, 69, and 101) have a patient removed
		# same as above but for individual measurements
		ij_outliers_meas[[i,j]] <- inner_join(vl_subtypes_comb2[[1,4]], residuals_removal_mx[[i,j]]$outl_meas, by="patientindex")
		print("Patients with outlying measurements")
		print(length(unique(ij_outliers_meas[[i,j]]$patientindex))) # 10 patients
		print("Distribution of outliers (measurements) in phylotypes")
		print( table(ij_outliers_meas[[i,j]]$cluster.x) ) # no outlier measurements from potential VL VOIs as well
	}
}

#saveRDS(residuals_removal_mx, glue("{RDS_PATH}/residuals_removal_mx.rds"))
residuals_removal_mx <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))

# plot distribution of follow-up times and positive slopes (potentially leading to biased estimates) + filter original CD4 df to exclude short follow-ups & >90% quantile slopes
plot_distr_cd4_follow_ups <- function(cd4_df, n_measurements_to_consider, mcs_subtype_choice) {
	# calculate slopes
	cd4_dt <- as.data.table(cd4_df)
	cd4_dt_slopes <- cd4_dt[,list(phylotype=phylotype, intercept=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[2]), by=patientindex]
	#View(cd4_dt_slopes)
	pl <- ggplot(cd4_dt_slopes, aes(slope)) + #follow_up_cat
		geom_histogram()
	ggsave(plot=pl, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_inspect_cd4_slopes.jpg"), width=10, height=8, dpi=300)
	
	print("max slope")
	print(max(cd4_dt_slopes$slope))
	print("min slope")
	print(min(cd4_dt_slopes$slope))
	print("median slope")
	quantile_50_slopes <- quantile(cd4_dt_slopes$slope, probs = 0.5)
	print(quantile_50_slopes) #-48
	print("95% quantile slopess")
	quantile_95_slopes <- quantile(cd4_dt_slopes$slope, probs = 0.95)
	print(quantile_95_slopes) #348
	quantile_90_slopes <- quantile(cd4_dt_slopes$slope, probs = 0.90)
	print(quantile_90_slopes) #99

	# follow-up time distributions
	cd4_df_slopes <- as.data.frame(cd4_dt_slopes)
	if(n_measurements_to_consider == "all") {
		cd4_df_follup <- cd4_df %>% group_by(patientindex) %>% summarise(follow_up_time = years_since_1cd4[measurement_id==max(measurement_id)], phylotype=phylotype) %>% ungroup()
	}
	else if(n_measurements_to_consider == 2) {
		cd4_df_follup <- cd4_df %>% group_by(patientindex) %>% filter(max(measurement_id) <= 2) %>% summarise(follow_up_time = years_since_1cd4[measurement_id==max(measurement_id)], phylotype=phylotype) %>% ungroup()
	}
	
	# custom categories for follow-up time
	cd4_df_follup <- cd4_df_follup %>% mutate(
		follow_up_cat = dplyr::case_when(follow_up_time>0 & follow_up_time <= 0.085 ~ "[0-1 month]", follow_up_time>0.085 & follow_up_time<=0.252 ~ "(1-3 months]",
																																			follow_up_time>0.252 & follow_up_time<=0.504 ~ "(3-6 months]", follow_up_time>0.504 & follow_up_time<=1.002 ~ "(6 months - 1 year]",
																																			follow_up_time>1.002 & follow_up_time<=3.002 ~ "(1-3 years]", follow_up_time>3.002 & follow_up_time<=5.002 ~ "(3-5 years]",
																																			follow_up_time>5.002 & follow_up_time<=10.002 ~ "(5-10 years]", follow_up_time>10.002 ~ ">10 years"),
		follow_up_cat = factor(follow_up_cat,level = c("[0-1 month]","(1-3 months]","(3-6 months]","(6 months - 1 year]","(1-3 years]","(3-5 years]","(5-10 years]",">10 years")))
	
	# summarise follow-ups for all data
	cd4_df_follup_summ <- cd4_df_follup %>% group_by(follow_up_cat) %>% summarise(n=n())
	print(cd4_df_follup_summ)
	pl3 <- ggplot(cd4_df_follup_summ, aes(x = follow_up_cat, y = n)) + #follow_up_cat
		geom_bar(stat = "identity", fill = "skyblue", color = "black") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
	ggsave(plot=pl3, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_inspect_cd4_follow_ups_overall_{n_measurements_to_consider}_meas.jpg"), width=10, height=8, dpi=300)

	# summarise for each phylotype
	cd4_df_follup_summ_pt <- cd4_df_follup %>% group_by(follow_up_cat, phylotype) %>% summarise(n=n())
	pl4 <- ggplot(cd4_df_follup_summ_pt, aes(x = follow_up_cat, y = n)) +
		geom_bar(stat = "identity", fill = "skyblue", color = "black") +
		facet_wrap(~phylotype, scales = "free") + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=4))
	ggsave(plot=pl4, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_inspect_cd4_follow_ups_phylotypes_{n_measurements_to_consider}_meas.jpg"), width=20, height=15, dpi=300)

	# summarise follow-ups for slopes >90% quantile (e.g. for B30 = 99) across phylotypes
	cd4_df_follup_slopes <- inner_join(cd4_df_follup, cd4_df_slopes, by="patientindex")
	cd4_df_follup_slopes_higher <- cd4_df_follup_slopes[cd4_df_follup_slopes$slope >= as.numeric(unname(quantile_90_slopes)),]
	#View(cd4_df_follup_slopes_higher)
	cd4_df_follup_slope_pt <- cd4_df_follup_slopes_higher %>% group_by(follow_up_cat, phylotype.x) %>% summarise(n=n())
	pl5 <- ggplot(cd4_df_follup_slope_pt, aes(x = follow_up_cat, y = n)) +
		geom_bar(stat = "identity", fill = "skyblue", color = "black") +
		facet_wrap(~phylotype.x, scales = "free") + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=4))
	ggsave(plot=pl5, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice}_inspect_cd4_follow_ups_phylotypes_higher_slopes_{n_measurements_to_consider}_meas.jpg"), width=20, height=15, dpi=300)
	
	# List of patients to remove: 1) only 2 measurements, (2) one only one month apart from the other, (3) and in >90% quantile of slopes (high positive value) - this last one regardless of n meas or time of follow-up
	
	if(n_measurements_to_consider == 2) {
		cd4_df_follup_short <- unique( cd4_df_follup[cd4_df_follup$follow_up_cat == "[0-1 month]",]$patientindex )
		print("patients removed 1")
		print(length(cd4_df_follup_short))
		
		cd4_df_slopes <- as.data.frame(cd4_dt_slopes)
		cd4_df_unlik_high_slopes <- unique( cd4_df_slopes[cd4_df_slopes$slope >= as.numeric(unname(quantile_90_slopes)),]$patientindex )
		print("patients removed 2")
		print(length(cd4_df_unlik_high_slopes))
		
		print("patients removed both")
		union_rm <- union(cd4_df_follup_short, cd4_df_unlik_high_slopes)
		print(length(union_rm))
		
		cd4_df_filt1 <- cd4_df[!(cd4_df$patientindex %in% cd4_df_follup_short),]
		cd4_df_filt2 <- cd4_df[!(cd4_df$patientindex %in% union_rm),]
		return(list(cd4_df_filt1=cd4_df_filt1, cd4_df_filt2=cd4_df_filt2))
	}
	
}

cd4_meas_filt_short_fups <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		# all CD4 measurements
		plot_distr_cd4_follow_ups(residuals_removal_mx[[i,j]]$cd4_df_filt, n_measurements_to_consider = "all", glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		# Only patients that have 2 measurements
		cd4_meas_filt_short_fups[[i,j]] <- plot_distr_cd4_follow_ups(residuals_removal_mx[[i,j]]$cd4_df_filt, n_measurements_to_consider = 2, glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
	}
}

saveRDS(cd4_meas_filt_short_fups, glue("{RDS_PATH}/cd4_meas_filt_short_fups.rds"))

# check model again after removing outlying indivs and measurements
# IMPORTANT: not doing boot for A1, CRF, and C because of singularFits error
#lmm2 <- lmm_sqrt2 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
#lmm3 <- lmm_sqrt3 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
lmm4 <- lmm_sqrt4 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		if(j <= 3) {
			print("model CD4")
			#lmm2[[i,j]] <- fit_cd4_model("cd4", residuals_removal_mx[[i,j]]$cd4_df_filt, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4"), do_boot=FALSE)
			#lmm3[[i,j]] <- fit_cd4_model("cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt1, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4-filt_short"), do_boot=FALSE)
			lmm4[[i,j]] <- fit_cd4_model("cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt2, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4-filt_short_high_slopes"), do_boot=FALSE)
			print("model sqrt CD4")
			#lmm_sqrt2[[i,j]] <- fit_cd4_model("sqrt_cd4", residuals_removal_mx[[i,j]]$cd4_df_filt, glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4"), do_boot=FALSE)
			#lmm_sqrt3[[i,j]] <- fit_cd4_model("sqrt_cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt1, glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4-filt_short"), do_boot=FALSE)
			lmm_sqrt4[[i,j]] <- fit_cd4_model("sqrt_cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt2, glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4-filt_short_high_slopes"), do_boot=FALSE)
		} else {
			print("model CD4")
			#lmm2[[i,j]] <- fit_cd4_model("cd4", residuals_removal_mx[[i,j]]$cd4_df_filt, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4"), do_boot=TRUE)
			#lmm3[[i,j]] <- fit_cd4_model("cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt1, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4-filt_short"), do_boot=FALSE)
			lmm4[[i,j]] <- fit_cd4_model("cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt2, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4-filt_short_high_slopes"), do_boot=FALSE)
			print("model sqrt CD4")
			#lmm_sqrt2[[i,j]] <- fit_cd4_model("sqrt_cd4", residuals_removal_mx[[i,j]]$cd4_df_filt, glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4"), do_boot=TRUE)
			#lmm_sqrt3[[i,j]] <- fit_cd4_model("sqrt_cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt1, glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4-filt_short"), do_boot=FALSE)
			lmm_sqrt4[[i,j]] <- fit_cd4_model("sqrt_cd4", cd4_meas_filt_short_fups[[i,j]]$cd4_df_filt2, glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4-filt_short_high_slopes"), do_boot=FALSE)
		}
	}
}

# saveRDS(lmm2, glue("{RDS_PATH}/lmm2.rds"))
# saveRDS(lmm_sqrt2, glue("{RDS_PATH}/lmm_sqrt2.rds"))
# saveRDS(lmm3, glue("{RDS_PATH}/lmm3.rds"))
# saveRDS(lmm_sqrt3, glue("{RDS_PATH}/lmm_sqrt3.rds"))
#saveRDS(lmm4, glue("{RDS_PATH}/lmm4.rds"))
#saveRDS(lmm_sqrt4, glue("{RDS_PATH}/lmm_sqrt4.rds"))

r2_corrs2 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names)) 
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		print("untrans CD4 - outlying filter")
		r2.corr.mer(lmm2[[1,4]]$fit)
		print("untrans CD4 - short follow-up filter")
		r2.corr.mer(lmm3[[1,4]]$fit)
		print("untrans CD4 - short follow-up + high slope filter")
		r2.corr.mer(lmm4[[1,4]]$fit)
		print("sqrt CD4 - outlying filter")
		r2.corr.mer(lmm_sqrt2[[1,4]]$fit)
		print("sqrt CD4 - short follow-up filter")
		r2.corr.mer(lmm_sqrt3[[1,4]]$fit)
		print("sqrt CD4 - short follow-up + high slope filter")
		r2.corr.mer(lmm_sqrt4[[1,4]]$fit)
	}
}

# For B30
# ~0.85 (previously 0.78)
# # ~0.88 (previosly 0.87)

### lmm1 random effects
# Groups       Name             Variance Std.Dev. Corr 
# patientindex (Intercept)      45678.98 213.726       
# years_since_1cd4   900.22  30.004  -0.54
# phylotype    (Intercept)        592.41  24.340       
# years_since_1cd4    54.24   7.365  0.49 
# Residual                      21061.44 145.126       
# Number of obs: 31562, groups:  patientindex, 8948; phylotype, 151

### lmm2 random effects
# Groups       Name             Variance Std.Dev. Corr 
# patientindex (Intercept)      46720.16 216.148       
# years_since_1cd4  1348.78  36.726  -0.48
# phylotype    (Intercept)        657.62  25.644       
# years_since_1cd4    65.05   8.066  0.72 
# Residual                      13379.87 115.671       
# Number of obs: 31551, groups:  patientindex, 8948; phylotype, 151

# variance of patient and phylotype random effects increase, correlations similar, but overall variance decreased
