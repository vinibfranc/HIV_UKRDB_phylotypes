libs_load <- c("ggplot2","dplyr","ggpubr","data.table","glue","lme4","lubridate","viridis", "boot", "parallel") #,"sjPlot"
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

NREP <- 10
NCPU <- 6

fit_cd4_model <- function(outcome, fostr, cd4_df, mcs_subtype_choice_cd4_transf, do_boot) {
	
	formula <- as.formula(paste(outcome, fostr)) 
	fit <- lmer(formula, data=cd4_df, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
	
	randeff_pt <- ranef(fit)$phylotype
	colnames(randeff_pt) <- c("intercept_randeff", "years_since_1cd4_randeff")
	randeff_pt$phylotype <- rownames(randeff_pt)
	
	comb_rand_fixed_pt1 <- coef(fit)$phylotype
	comb_rand_fixed_pt <- comb_rand_fixed_pt1 %>% dplyr::select("(Intercept)", "years_since_1cd4")
	colnames(comb_rand_fixed_pt) <- c("intercept_fixedeff", "slope_phylotype")
	comb_rand_fixed_pt$phylotype <- rownames(comb_rand_fixed_pt)
	combined_coeffs <- inner_join(randeff_pt, comb_rand_fixed_pt, by="phylotype")
	
	resid_sd <- sigma(fit)
	
	# variance-covariance matrix
	rand_eff_sd_df <- as.data.frame(lme4::VarCorr(fit))
	rand_eff_sd_df <- rand_eff_sd_df[, c("grp", "var1", "var2", "vcov", "sdcor")]
	colnames(rand_eff_sd_df) <- c("grouping_factor", "effect1", "effect2", "std_dev", "correlation")
	
	# bootstrap data to get std errors of random phylotype slopes, etc
	if (do_boot) {
		.random_effect_se <- function(data, fo, indices) {
  data <- data %>% dplyr::select(patientindex, phylotype, cd4, sqrt_cd4, years_since_1cd4, age_group, sexid, exposureid)
  boot_data <- data[indices, ]
  
  boot_fit <- tryCatch(
    lmer(fo, data = boot_data, REML = FALSE, control = lmerControl(optimizer = "bobyqa")),
    error = function(e) return(rep(NA, length(unique(data$phylotype)) + 8))  # updated length to match added parameters
  )
  
  if (is.numeric(boot_fit)) return(boot_fit)

  # Extract phylotype random slope effects
  ranef_vals <- ranef(boot_fit)$phylotype %>% dplyr::select(years_since_1cd4)
  all_phylotypes <- unique(data$phylotype)
  ranef_vals_vec <- sapply(all_phylotypes, function(ptype) {
    if (ptype %in% rownames(ranef_vals)) {
      ranef_vals[ptype, "years_since_1cd4"]
    } else {
      NA
    }
  })

  # Extract SDs and correlation from VarCorr
  vc_all <- VarCorr(boot_fit)
  
  # From phylotype
  vc_pt <- vc_all$phylotype
  sd_intercept_pt <- attr(vc_pt, "stddev")[1]
  sd_slope_pt <- attr(vc_pt, "stddev")[2]
  correlation_pt <- attr(vc_pt, "correlation")[1, 2]
  
  # From patientindex
  vc_patient <- vc_all$patientindex
  sd_intercept_patient <- attr(vc_patient, "stddev")[1]
  sd_slope_patient <- attr(vc_patient, "stddev")[2]
  correlation_patient <- attr(vc_patient, "correlation")[1, 2]
  
  # Residual SD
  resid_sd <- sigma(boot_fit)
  
  # Combine all in one vector
  c(ranef_vals_vec,
    sd_intercept_pt, sd_slope_pt, correlation_pt,
    sd_intercept_patient, sd_slope_patient, correlation_patient,
    resid_sd)
	}
		
		print("Running boot")
		boot_res <- boot::boot(data=cd4_df, fo=formula, statistic=.random_effect_se, R=NREP, parallel="multicore", ncpus=NCPU)
		print(boot_res)
		
		# Extract random slope SEs
		n_phylotypes <- length(unique(cd4_df$phylotype))
		random_effect_ses <- apply(boot_res$t[, 1:n_phylotypes, drop=FALSE], 2, sd)
		
		comb_rand_fixed_pt2 <- randeff_pt %>% dplyr::select(years_since_1cd4_randeff)
		colnames(comb_rand_fixed_pt2) <- c("phylotype_coef")
		comb_rand_fixed_pt2$phylotype <- rownames(comb_rand_fixed_pt2)
		comb_rand_fixed_pt2 <- comb_rand_fixed_pt2 %>% dplyr::select(phylotype, phylotype_coef)
		comb_rand_fixed_pt2$t_value <- comb_rand_fixed_pt2$phylotype_coef / random_effect_ses
		comb_rand_fixed_pt2$se <- random_effect_ses
		comb_rand_fixed_pt2$p_value <- 2 * (1 - pnorm(abs(comb_rand_fixed_pt2$t_value)))
		comb_rand_fixed_pt2 <- comb_rand_fixed_pt2[order(comb_rand_fixed_pt2$phylotype_coef), ]
		
		estim_ml_boot <- ggplot(comb_rand_fixed_pt2, aes(x = phylotype, y = phylotype_coef)) +
			geom_point() +
			geom_errorbar(aes(ymin = phylotype_coef - 1.96 * se, ymax = phylotype_coef + 1.96 * se)) +
			theme(axis.text=element_text(size=5), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
			labs(title="Phylotype Effects on CD4 Decline", y="Effect Estimate")
		
		ggsave(plot=estim_ml_boot, filename=glue("{RESULTS_PATH}/06_cd4_ml_model/{mcs_subtype_choice_cd4_transf}_boot_estim.jpg"), width=15, height=8, dpi=300)
		
		# Collect additional bootstrapped SDs and correlation
		boot_sd_info <- boot_res$t[, (n_phylotypes + 1):(n_phylotypes + 7)]
		colnames(boot_sd_info) <- c("sd_intercept_pt", "sd_slope_pt", "correlation_pt", "sd_intercept_patient", "sd_slope_patient", "correlation_patient","resid_sd")
		boot_sd_summary <- apply(boot_sd_info, 2, function(x) c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE)))
		
		return(list(
			fit = fit,
			combined_coeffs = combined_coeffs,
			table_pt_estims_p = comb_rand_fixed_pt2,
			residual_sd = resid_sd,
			rand_eff_summary = as.data.frame(t(boot_sd_summary))
		))
	} else {
		return(list(
			fit = fit,
			combined_coeffs = combined_coeffs,
			residual_sd = resid_sd,
			rand_eff_summary = rand_eff_sd_df
		))
	}
}

# IMPORTANT: lme4 includes random intercept by default (no need to add e.g. "(1 + years_since_1cd4 | patientindex)", just "(years_since_1cd4 | patientindex)" is fine)
# Default behaviour is also to assume joint-normal distribution for intercept and slope and that they are correlated
cd4_model_form <- paste(" ~ years_since_1cd4 + age_group +  years_since_1cd4 * sexid + years_since_1cd4:exposureid + (years_since_1cd4 | patientindex) + (years_since_1cd4 | phylotype)")

lmm1 <- lmm_sqrt1 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
			if(nrow(cd4_subtypes_comb2[[i,j]]) > 0) {
				print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
				# define as ref category for fixed effects: backbone phylotype, age 30-39, sexid Male, exposureid Homo-bisexual
				# untransformed cd4 model: without bootstrap (raw first estim)
				cd4_subtypes_comb2[[i,j]] <- cd4_subtypes_comb2[[i,j]] %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[i,j]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
				lmm1[[i,j]] <- fit_cd4_model("cd4", cd4_model_form, cd4_subtypes_comb2[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4"), do_boot=F)
				# sqrt model: without bootstrap (raw first estim)
				lmm_sqrt1[[i,j]] <- fit_cd4_model("sqrt_cd4", cd4_model_form, cd4_subtypes_comb2[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4"), do_boot=F)
			} else {
				next
			}
		}
}
saveRDS(lmm1, glue("{RDS_PATH}/lmm1.rds"))
saveRDS(lmm_sqrt1, glue("{RDS_PATH}/lmm_sqrt1.rds"))

# inspect B30 models closer (indices 1,4 of matrix)
imcs <- 1 
jtree <- 4 
da <- cd4_subtypes_comb2[[imcs,jtree]]
da <- da %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[imcs,jtree]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
'
- Is intercept correlated with slope?? 
	* yes, but cor is reversed and much smaller with fixed eff...
- Recreate distribution of fixed effect slopes and rand eff slopes; compare variance 
	* FE has large outliers 
- Is fixed effect slopes centered zero & approx normal? 
	* centered, but large outliers
- excl backbone & repeat 
'

# random effect of phylotype (as two-dimensional joint normal)
form1 <- as.formula(paste("cd4 ~ years_since_1cd4 + age_group + 
    years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
    (1 + years_since_1cd4 | patientindex) + (1 + years_since_1cd4 | phylotype)")) 
f1 <- lmer(form1, data=da, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
f1reff <-  coef(f1)$phylotype 
f1slope  <- f1reff$years_since_1cd4 
plot( f1reff[,1], f1reff$years_since_1cd4 )
cor( cbind(f1reff[,1], f1reff$years_since_1cd4 ))

f1reff2 <- ranef(f1)$phylotype 
f1slope2  <- f1reff2$years_since_1cd4
plot( f1reff2[,1], f1reff2$years_since_1cd4 )

# fixed effect of phylotype 
form2 <- as.formula(paste("cd4 ~ years_since_1cd4*phylotype + age_group + 
    years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
    (1 + years_since_1cd4 | patientindex) ")) 
f2 <- lmer(form2, data=da, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
f2coef <- summary(f2)$coefficients 
f2intercept <- f2coef[ grepl(rownames(f2coef), patt='^phylotype') , 1]
f2slope <- f2coef[ grepl(rownames(f2coef), patt=':phylotype') , 1]
plot( f2intercept, f2slope )
cor( f2intercept, f2slope )
summary( lm( f2intercept ~ f2slope ) )
# [1] -0.3532388
plot( f1slope[-1], f2slope )
plot( f1slope[-1], f2slope, ylim=c(-50,10) )
cor( f1slope[-1],  f2slope )
# [1] 0.249739

# gives same result as f1 (random eff not as two-dim joint normal), so using this simpler version
form3 <- as.formula(paste("cd4 ~ years_since_1cd4 + age_group + 
    years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
    (years_since_1cd4 | patientindex) + (years_since_1cd4 | phylotype)")) 
f3 <- lmer(form3, data=da, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
f3reff <-  coef(f3)$phylotype 
f3slope  <- f3reff$years_since_1cd4 
f3intercept <- f3reff[,1] 
plot( f3intercept, f3slope )

# drop backbone 
da$backbone <- da$phylotype == 153
d1 <- da[!da$backbone,]
form4 <- as.formula(paste("cd4 ~ years_since_1cd4 + age_group + 
    years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
    (years_since_1cd4 | patientindex) + (years_since_1cd4 | phylotype)")) 
f4 <- lmer(form4, data=d1, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
f4reff <-  coef(f4)$phylotype 
f4slope  <- f4reff$years_since_1cd4 
f4intercept <- f4reff[,1] 
plot( f4intercept, f4slope )
## very tight cor! 
cor( f4intercept, f4slope )
# [1] 0.9781044

# fixed eff backbone, doesn't work 
form5 <- as.formula(paste("cd4 ~ years_since_1cd4 + age_group + 
    years_since_1cd4 * sexid + years_since_1cd4:exposureid + 
    	years_since_1cd4:backbone + 
    (years_since_1cd4 | patientindex) + (years_since_1cd4 | phylotype)")) 
f5 <- lmer(form5, data=da, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) # singular fit...
summary( f5)
f5reff <-  coef(f5)$phylotype 
head(f5reff) 
f5slope  <- f5reff$years_since_1cd4 
f5intercept <- f5reff[,1] 
plot( f5intercept, f5slope )

AIC(f1)
# [1] 424831.4
AIC(f2)
# [1] 425018.8
AIC(f3)
# [1] 424831.4
AIC(f4)
# [1] 169701.5; not incl backbone 
AIC(f5)
# [1] 424829.3

# Visualise and remove outlying individuals in a distribution of maximum-likehood slopes and intercepts + outlying measurements (thanks Chris Wymant for advice)
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
	#View(zs)
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
	print("nrow BEFORE removing individuals: ")
	print(nrow(cd4_df))
	cd4_df_filt <- cd4_df[!cd4_df$patientindex %in% outl_indivs$patientindex,]
	print("nrow after removing individuals: ")
	print(nrow(cd4_df_filt))
	cd4_df_filt <- cd4_df_filt[cd4_df_filt$sqrt_cd4 <= 45,] # same as outl_meas
	print("nrow after removing measurements: ")
	print(nrow(cd4_df_filt))
	cd4_df_filt <- cd4_df_filt[cd4_df_filt$patientindex %in% unique(cd4_df_filt$patientindex[duplicated(cd4_df_filt$patientindex)]),]
	print("nrow after double-checking if by removing a measurement a patient ended up with n<2: ")
	print(nrow(cd4_df_filt))
	
	list(outl_indivs=outl_indivs, outl_meas=outl_meas, cd4_df_filt=cd4_df_filt, zs=zs)
}

# Run removal of outlying individuals and measurements
cd4_subtypes_comb2 <- readRDS(glue("{RDS_PATH}/cd4_subtypes_comb2.rds"))
lmm1 <- readRDS(glue("{RDS_PATH}/lmm1.rds"))
residuals_removal_mx <-  ij_outliers_pts <- ij_outliers_meas <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		residuals_removal_mx[[i,j]] <- vis_residuals(lmm1[[i,j]]$fit, cd4_subtypes_comb2[[i,j]], glue("{min_cl_size_choices[i]}-{tree_names[j]}")) # 130 outlier indivs, 11 outlier measurements
		# inner_join of outlier individuals with cluster membership to see if removing cd4 VOI patients
		ij_outliers_pts[[i,j]] <- inner_join(cd4_subtypes_comb2[[i,j]], residuals_removal_mx[[i,j]]$outl_indivs, by="patientindex")
		print("Distribution of outliers (slope and intercept) in phylotypes")
		print( table(ij_outliers_pts[[i,j]]$cluster) ) # none of subtype B cd4 VOIs (PTs 20, 40, 69, and 101) have a patient removed
	}
}

saveRDS(residuals_removal_mx, glue("{RDS_PATH}/residuals_removal_mx.rds"))

# run boostrap on random effs for both cd4 and sqrt_cd4 for all subtypes
residuals_removal_mx <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))

NREP <- 1000
NCPU <- 6

lmm2 <- lmm_sqrt2 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		if(nrow(cd4_subtypes_comb2[[i,j]]) > 0) {
			print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
			# define as ref category for fixed effects: backbone phylotype, age 30-39, sexid Male, exposureid Homo-bisexual
			# untransformed cd4 model: without bootstrap (raw first estim)
			db <- residuals_removal_mx[[i,j]]$cd4_df_filt
			db <- db %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[i,j]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
			lmm2[[i,j]] <- fit_cd4_model("cd4", cd4_model_form, db, glue("{min_cl_size_choices[i]}-{tree_names[j]}-cd4"), do_boot=TRUE)
			gc()
			# sqrt model: without bootstrap (raw first estim)
			#lmm_sqrt2[[i,j]] <- fit_cd4_model("sqrt_cd4", cd4_model_form, db, glue("{min_cl_size_choices[i]}-{tree_names[j]}-sqrt_cd4"), do_boot=FALSE) #TRUE
			gc()
		} else {
			next
		}
	}
}

saveRDS(lmm2, glue("{RDS_PATH}/lmm2.rds"))
#saveRDS(lmm_sqrt2, glue("{RDS_PATH}/lmm_sqrt2.rds"))

# join signif together
lmm2 <- readRDS(glue("{RDS_PATH}/lmm2.rds"))
#lmm_sqrt2 <- readRDS(glue("{RDS_PATH}/lmm_sqrt2.rds"))
cd4_pot_cd4_vois <- sqrt_cd4_pot_cd4_vois <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		cd4_pot_cd4_vois[[i,j]] <- lmm2[[i,j]]$table_pt_estims_p[lmm2[[i,j]]$table_pt_estims_p$t_value < 0 & lmm2[[i,j]]$table_pt_estims_p$p_value <= 0.05,]
		if(nrow(cd4_pot_cd4_vois[[i,j]]) > 0) cd4_pot_cd4_vois[[i,j]]$mcs_subtype <- glue("{min_cl_size_choices[i]}-{tree_names[j]}")
		#sqrt_cd4_pot_cd4_vois[[i,j]] <- lmm_sqrt2[[i,j]]$table_pt_estims_p[lmm_sqrt2[[i,j]]$table_pt_estims_p$t_value < 0 & lmm_sqrt2[[i,j]]$table_pt_estims_p$p_value <= 0.05,]
		#if(nrow(sqrt_cd4_pot_cd4_vois[[i,j]]) > 0) sqrt_cd4_pot_cd4_vois[[i,j]]$mcs_subtype <- glue("{min_cl_size_choices[i]}-{tree_names[j]}")
	}
}
cd4_pot_cd4_vois_df <- do.call(rbind, cd4_pot_cd4_vois)
cd4_pot_cd4_vois_df <- cd4_pot_cd4_vois_df[!is.na(cd4_pot_cd4_vois_df$phylotype_coef),]
saveRDS(cd4_pot_cd4_vois_df, glue("{RDS_PATH}/cd4_pot_cd4_vois_df.rds")) #24
cd4_pot_cd4_vois_df_summ <- cd4_pot_cd4_vois_df %>% group_by(mcs_subtype) %>% summarise(n=n())
print(cd4_pot_cd4_vois_df_summ)
#sqrt_cd4_pot_cd4_vois_df <- do.call(rbind, sqrt_cd4_pot_cd4_vois)
#sqrt_cd4_pot_cd4_vois_df <- sqrt_cd4_pot_cd4_vois_df[!is.na(sqrt_cd4_pot_cd4_vois_df$phylotype_coef),]
#saveRDS(sqrt_cd4_pot_cd4_vois_df, glue("{RDS_PATH}/sqrt_cd4_pot_cd4_vois_df.rds")) # 29
#sqrt_cd4_pot_cd4_vois_df_summ <- sqrt_cd4_pot_cd4_vois_df %>% group_by(mcs_subtype) %>% summarise(n=n())
#print(sqrt_cd4_pot_cd4_vois_df_summ)

# agreements and disagreements cd4 and sqrt_cd4
nrow(cd4_pot_cd4_vois_df) #24
#nrow(sqrt_cd4_pot_cd4_vois_df) #29
#agr_cd4_sqrt <- inner_join(cd4_pot_cd4_vois_df, sqrt_cd4_pot_cd4_vois_df, by=c("mcs_subtype","phylotype")) # 20 matches
#cd4_but_not_sqrt <- anti_join(cd4_pot_cd4_vois_df, sqrt_cd4_pot_cd4_vois_df, by=c("mcs_subtype","phylotype")) # 3
#sqrt_but_not_cd4 <- anti_join(sqrt_cd4_pot_cd4_vois_df, cd4_pot_cd4_vois_df, by=c("mcs_subtype","phylotype")) # 9

# Table S8: ML CD4 model with random effect on phylotype
cd4_ml_rand_eff_tables <- readRDS(glue("{RDS_PATH}/lmm2.rds"))
cd4_ml_rand_eff_tables[[1,1]]$table_pt_estims_p$subtype <- tree_names[1]; cd4_ml_rand_eff_tables[[1,2]]$table_pt_estims_p$subtype <- tree_names[2] 
cd4_ml_rand_eff_tables[[1,3]]$table_pt_estims_p$subtype <- tree_names[3]; cd4_ml_rand_eff_tables[[1,4]]$table_pt_estims_p$subtype <- tree_names[4] 
cd4_ml_rand_eff_table <- rbind( cd4_ml_rand_eff_tables[[1,1]]$table_pt_estims_p, cd4_ml_rand_eff_tables[[1,2]]$table_pt_estims_p, cd4_ml_rand_eff_tables[[1,3]]$table_pt_estims_p, cd4_ml_rand_eff_tables[[1,4]]$table_pt_estims_p )
cd4_ml_rand_eff_table <- cd4_ml_rand_eff_table %>% dplyr::select( subtype, phylotype, phylotype_coef, se, t_value, p_value )
cd4_ml_rand_eff_table <- cd4_ml_rand_eff_table[!is.na(cd4_ml_rand_eff_table$se),]
cd4_ml_rand_eff_table <- cd4_ml_rand_eff_table %>% arrange( p_value, phylotype_coef )
cd4_ml_rand_eff_table$phylotype_coef <- round(cd4_ml_rand_eff_table$phylotype_coef, 3)
cd4_ml_rand_eff_table$se <- round(cd4_ml_rand_eff_table$se, 3)
cd4_ml_rand_eff_table$t_value <- round(cd4_ml_rand_eff_table$t_value, 3)
cd4_ml_rand_eff_table$p_value <- signif(cd4_ml_rand_eff_table$p_value, digits = 2)
options(scipen=999)
cd4_ml_rand_eff_table <- cd4_ml_rand_eff_table[order(cd4_ml_rand_eff_table$p_value, cd4_ml_rand_eff_table$phylotype_coef),]
write.csv( cd4_ml_rand_eff_table, file=glue("{RESULTS_PATH}/tables/tableS8.csv"), quote=F, row.names = F )
saveRDS(cd4_ml_rand_eff_table, glue("{RDS_PATH}/cd4_ml_rand_eff_table.rds"))

# For table S8 as well, extract residual sd, sds interc and slopes, and corr interc and slope
cd4_ml_rand_eff_tables_sds <- readRDS(glue("{RDS_PATH}/lmm2.rds"))
cd4_ml_rand_eff_tables_sds <- matrix(
	data = lapply(seq_along(cd4_ml_rand_eff_tables_sds), function(i) {
		cd4_ml_rand_eff_tables_sds[[i]]$rand_eff_summary}),
	nrow = nrow(cd4_ml_rand_eff_tables_sds),ncol = ncol(cd4_ml_rand_eff_tables_sds),byrow = FALSE)
cd4_ml_rand_eff_tables_sds[[1,1]]$subtype <- tree_names[1]; cd4_ml_rand_eff_tables_sds[[1,2]]$subtype <- tree_names[2] 
cd4_ml_rand_eff_tables_sds[[1,3]]$subtype <- tree_names[3]; cd4_ml_rand_eff_tables_sds[[1,4]]$subtype <- tree_names[4] 
cd4_ml_rand_eff_tables_sds_mcs30 <- do.call(rbind, lapply(1:4, function(j) {
	cd4_ml_rand_eff_tables_sds[[1, j]]
}))
cd4_ml_rand_eff_tables_sds_mcs30$mean <- round(cd4_ml_rand_eff_tables_sds_mcs30$mean, 3)
cd4_ml_rand_eff_tables_sds_mcs30$sd <- round(cd4_ml_rand_eff_tables_sds_mcs30$sd, 3)
write.csv( cd4_ml_rand_eff_tables_sds_mcs30, file=glue("{RESULTS_PATH}/tables/tableS8_2.csv"), quote=F, row.names = T )

### Fig. 1 VL heatmap: After identyfing outliers using random effects model, fit individual fixed eff model for each phylotype
# backbone = most likely to reflect contain ancestral state of the virus
fit_cd4_model_fixed_eff_pts_vs_backbone <- function(outcome, fostr, cd4_df) {
	
	# assumes df already contain VOI, backbone and non-voi phylotypes
	formula <- as.formula(paste(outcome, fostr)) 
	fit <- lmer(formula, data=cd4_df, REML=FALSE, control = lmerControl(optimizer = "bobyqa"))
	
	summ_df <- data.frame(phylotype=names(summary(fit)$coefficients[,1]), estimate=round( unname(summary(fit)$coefficients[,1]), 3), 
																							stderr=round( unname(summary(fit)$coefficients[,2]), 3), t_val=round( unname(summary(fit)$coefficients[,3]), 3))
	summ_df <- summ_df[order(summ_df$estimate), ]
	summ_df$norm <- summ_df$t_val |> pnorm() # note left tail only (negative t values)
	summ_df$p_value <- summ_df$norm / 2 # one tail test
	summ_df$p_value <- signif(summ_df$p_value, digits = 2)
	summ_df_pts <- summ_df[ grepl(summ_df$phylotype, pattern = 'years_since_1cd4:phylotype'), ] #years_since_1cd4:phylotype
	summ_df_pts$phylotype <- gsub("years_since_1cd4:phylotype", "", as.character(summ_df_pts$phylotype))
	
	summ_df_pts_interc <- summ_df[ grepl(summ_df$phylotype, pattern = '^phylotype'), ]
	summ_df_pts_interc$phylotype <- gsub("^phylotype", "", as.character(summ_df_pts_interc$phylotype))
	
	return(list(all_covar=summ_df, pt_effs=summ_df_pts, pt_interc=summ_df_pts_interc))
}

# For fig. 1 (tree+vl & cd4 estimates): run fixed effects model on each phylotype and non-VOIs against backbone (mcs=30, subtype B)
cd4_model_form_fixed_eff <- cd4_model_form_with_intercept_pt <- " ~ years_since_1cd4 + age_group + years_since_1cd4 * sexid + years_since_1cd4:exposureid + years_since_1cd4*phylotype + (years_since_1cd4 | patientindex)"
d_fixedeff_all_b <- d_fixedeff_all_b_ref_bb <- res_fixeff_all_b_vs_bb <- list()
for(i in 1:154) {
	print(i)
	d_fixedeff2 <- residuals_removal_mx[[ 1,4 ]]$cd4_df_filt
	d_fixedeff_all_b[[i]] <- d_fixedeff2
	d_fixedeff_all_b[[i]]$phylotype <- ifelse( d_fixedeff2$cluster ==  as.integer(backbone_cl_control[[ 1,4 ]]),
																																												yes=as.character(backbone_cl_control[[ 1,4 ]]),
																																												no=ifelse( d_fixedeff2$cluster == i,
																																																							yes= i,
																																																							no=as.character("non-VOI") ))
	d_fixedeff_all_b[[i]] <- d_fixedeff_all_b[[i]][ d_fixedeff_all_b[[i]]$phylotype != "non-VOI", ]
	#print(length(unique(d_fixedeff_all_b[[i]]$cluster)))
	if(i == as.integer(backbone_cl_control[[ 1,4 ]]) | length(unique(d_fixedeff_all_b[[i]]$cluster)) < 2) {
		next
	} else {
		d_fixedeff_all_b[[i]]$phylotype <- as.factor(d_fixedeff_all_b[[i]]$phylotype)
		print(table(d_fixedeff_all_b[[i]]$phylotype))
		d_fixedeff_all_b_ref_bb[[i]] <- d_fixedeff_all_b[[i]] %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[ 1,4 ]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
		res_fixeff_all_b_vs_bb[[i]] <- fit_cd4_model_fixed_eff_pts_vs_backbone("cd4",cd4_model_form_fixed_eff, d_fixedeff_all_b_ref_bb[[i]])
	}
}

# join all phylotype slopes together
res_fixeff_all_b_vs_bb_combined <- do.call(rbind, lapply(res_fixeff_all_b_vs_bb, function(x) x$pt_effs))
#res_fixeff_all_b_vs_bb_combined <- res_fixeff_all_b_vs_bb_combined[res_fixeff_all_b_vs_bb_combined$phylotype != "non-VOI",]
saveRDS(res_fixeff_all_b_vs_bb_combined, glue("{RDS_PATH}/res_fixeff_all_b_vs_bb_combined.rds") )
