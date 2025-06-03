libs_load <- c("ggplot2","dplyr","ggpubr","data.table","glue","lubridate","sjPlot","viridis", "rstan", "brms","bayesplot",  "ape", "forcats","pbmcapply", "posterior")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"
system(glue("mkdir -p {RESULTS_PATH}/06_cd4_bayes_model"))

min_cl_size_choices <- c(30, 50, 100)
tree_names <- c("A_A1","CRF_02_AG","C","B")

residuals_removal_mx <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))
backbone_cl_control <- readRDS(glue("{RDS_PATH}/backbone_cl_control.rds"))

rstan_options(auto_write = TRUE)

# Sensitivity analysis: Bayesian version of model with random effects on phylotype
fit_bayes_cd4_model_randeff_pts <- function(outcome, fostr, warmup_, cd4_df, change_prior=TRUE) {
	
	formula <- as.formula(paste(outcome, fostr))
	
	gp <- get_prior(formula, data=cd4_df, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
	print(gp)
	
	if(change_prior) {
		gp$prior[gp$coef == "years_since_1cd4" & gp$class == "sd" & gp$group == "phylotype"] <- "normal(0, 20)"
		gp$prior[gp$coef == "Intercept" & gp$class == "sd" & gp$group == "phylotype"] <- "normal(0, 50)"
		gp$prior[gp$class == "cor" & gp$group == "phylotype"] <- "lkj(2)"
		gp$prior[gp$class == "cor" & gp$group == "patientindex"] <- "lkj(2)"
		print(gp)
		print("changing prior")
		bayes_fit <- brm(formula, data=cd4_df, prior=gp, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
	} else {
		bayes_fit <- brm(formula, data=cd4_df, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
	}
	loo_model <- loo(bayes_fit) #moment_match = TRUE
	print("========")
	print("LOO:")
	print(loo_model)
	print("========")
	
	summary_fit <- summary(bayes_fit)
	print(str(summary_fit))
	
	components <- list(
		fixed = summary_fit$fixed,
		spec_pars = summary_fit$spec_pars,
		random_patient = summary_fit$random$patientindex,
		random_phylotype = summary_fit$random$phylotype
	)
	
	# Combine into one data.frame
	diagnostics_df <- do.call(rbind, lapply(names(components), function(name) {
		df <- components[[name]]
		if (nrow(df) > 0) {
			df$source <- name
			return(df[, c("Rhat", "Bulk_ESS", "source")])
		}
	}))
	
	# Drop rows with all NA Rhat/Bulk_ESS (optional)
	diagnostics_df <- diagnostics_df[!(is.na(diagnostics_df$Rhat) & is.na(diagnostics_df$Bulk_ESS)), ]
	
	# Check conditions
	all_rhat_ok <- all(diagnostics_df$Rhat <= 1.01, na.rm = TRUE)
	print("All Rhat values are <= 1.01?")
	print(all_rhat_ok)
	all_ess_ok <- all(diagnostics_df$Bulk_ESS >= 400, na.rm = TRUE)
	print("All effective sample sizes are >= 400")
	print(all_ess_ok)
	
	gc()
	
	.get_p_val <- function(df) {
		df$p_value_freq <- 2 * (1 - pnorm(abs(df$Estimate) / df$Est.Error))
		df
	}
	
	ranef_pt <- ranef(bayes_fit)$phylotype
	p_pt_only <- as.data.frame(ranef_pt[, , "years_since_1cd4"])
	p_pt_only$phylotype <- rownames(ranef_pt)
	rownames(p_pt_only) <- NULL
	p_pt_only <- .get_p_val(p_pt_only)
	
	# get posterior draws for random slopes
	posterior_draws <- as_draws_df(bayes_fit)
	slope_cols <- grep("^r_phylotype\\[.*?,years_since_1cd4\\]", names(posterior_draws), value = TRUE)
	
	# calc Bayes p-values and classify based on posterior probability
	.bayes_p_vals <- lapply(slope_cols, function(col) {
		values <- posterior_draws[[col]]
		
		# posterior probabilities
		p_neg <- mean(values < 0)
		p_pos <- mean(values > 0)
		
		# two-sided Bayesian p-value (symmetric tail test)
		p_val <- 2 * min(p_neg, p_pos)
		
		# classification based on posterior probability
		evidence_label <- case_when(
			p_neg > 0.90 ~ "Strong negative slope",
			p_neg > 0.80  ~ "Moderate negative slope",
			p_neg < 0.10 ~ "Strong positive slope",
			p_neg < 0.20  ~ "Moderate positive slope",
			TRUE          ~ "No clear evidence"
		)
		
		# extract phylotype name from column name
		phylotype_name <- gsub("r_phylotype\\[(.*?),years_since_1cd4\\]", "\\1", col)
		
		data.frame(phylotype = phylotype_name, p_value_bayes = p_val, prob_negative = p_neg, evidence_label = evidence_label)
	})
	bayes_p_df <- do.call(rbind, .bayes_p_vals)
	# merge with your existing phylotype-level data
	p_pt_only <- merge(p_pt_only, bayes_p_df, by = "phylotype", all.x = TRUE)
	
	p_pt_intercepts <- as.data.frame(ranef_pt[, , "Intercept"])
	p_pt_intercepts$phylotype <- rownames(ranef_pt)
	rownames(p_pt_intercepts) <- NULL
	
	p_all_rows <- fixef(bayes_fit)
	p_all_rows <- as.data.frame(p_all_rows)
	p_all_rows$param <- rownames(p_all_rows)
	rownames(p_all_rows) <- NULL 
	p_all_rows <- p_all_rows %>% dplyr::select(param, Estimate, Est.Error, Q2.5, Q97.5)
	#p_all_rows <- .get_p_val(p_all_rows)
	
	print(summary_fit$spec_pars)
	print(summary_fit$random$patientindex)
	print(summary_fit$random$phylotype)
	sigma_sds <- data.frame(
		Param = c("sigma", "sd(Intercept) - patientindex", "sd(years_since_1cd4) - patientindex",
												"sd(Intercept) - phylotype", "sd(years_since_1cd4) - phylotype"),
		Estimate = c(summary_fit$spec_pars$Estimate, summary_fit$random$patientindex$Estimate[1:2],summary_fit$random$phylotype$Estimate[1:2]),
		Q2.5 = c(summary_fit$spec_pars$`l-95%`,summary_fit$random$patientindex$`l-95%`[1:2],summary_fit$random$phylotype$`l-95%`[1:2]),
		Q97.5 = c(summary_fit$spec_pars$`u-95%`,summary_fit$random$patientindex$`u-95%`[1:2],summary_fit$random$phylotype$`u-95%`[1:2])
	)
	
	# Helper function to extract correlation with CIs from a group
	extract_cor_with_ci <- function(group_summary, group_name) {
		cor_rows <- grep("^cor\\(", rownames(group_summary))
		if (length(cor_rows) == 0) return(NULL)
		
		cor_df <- as.data.frame(group_summary[cor_rows, ])
		print("cor_df")
		print(cor_df)
		cor_df$group <- group_name
		cor_df$param <- rownames(group_summary)[cor_rows]
		rownames(cor_df) <- NULL
		cor_df <- cor_df[, c("group", "param", "Estimate", "l-95% CI", "u-95% CI")]
		return(cor_df)
	}
	
	# Extract from each random effect group
	cor_patient <- extract_cor_with_ci(summary_fit$random$patientindex, "patientindex")
	cor_phylotype <- extract_cor_with_ci(summary_fit$random$phylotype, "phylotype")
	
	# Combine
	cor_all <- rbind(cor_patient, cor_phylotype)
	
	print(p_pt_only[ order(p_pt_only$Estimate), ])
	
	return(list(all_covar=p_all_rows, pt_effs=p_pt_only, p_pt_intercepts=p_pt_intercepts, p_cors=cor_all, sigma_sds=sigma_sds, loo_model=loo_model, rhat_values=diagnostics_df$Rhat, ess_values=diagnostics_df$Bulk_ESS)) # sigma=sigma, random_effects_sd=random_effects_sd,
}

ITER <- 10000
WARMUP10 <- ITER/10
WARMUP50 <- ITER/2
THIN <- 10
CORES_CHAINS <- 2

cd4_model_form_with_intercept_pt <- " ~ years_since_1cd4 + age_group + years_since_1cd4 * sexid + years_since_1cd4:exposureid + (years_since_1cd4 | phylotype) + (years_since_1cd4 | patientindex)"

bmm_res_ranef_pts <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		db <- residuals_removal_mx[[i,j]]$cd4_df_filt
		db <- db %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[i,j]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
		bmm_res_ranef_pts[[i,j]] <- fit_bayes_cd4_model_randeff_pts("cd4",cd4_model_form_with_intercept_pt,WARMUP50, db, change_prior=TRUE) #priors_slope_interc = priors_bayes
	}
}
saveRDS(bmm_res_ranef_pts, glue("{RDS_PATH}/bmm_res_ranef_pts.rds"))

# join all pt effects together for mcs=30
cd4_ranef_pts <- do.call(rbind, lapply(1:4, function(j) { 
	ss <- bmm_res_ranef_pts[[1, j]]$pt_effs 
	ss$subtype <- tree_names[j]
	ss
}))

# order by estimate
cd4_ranef_pts <- cd4_ranef_pts[order(cd4_ranef_pts$Estimate), ]
# order by posterior probability as none signif (conservative model)
cd4_ranef_pts2 <- cd4_ranef_pts[order(cd4_ranef_pts$prob_negative, decreasing = T), ]

cd4_ranef_pts2_table <- cd4_ranef_pts2
cd4_ranef_pts2_table <- cd4_ranef_pts2_table %>% dplyr::select(subtype, phylotype, Q2.5, Estimate, Q97.5, Est.Error, p_value_bayes, prob_negative) #evidence_label
cd4_ranef_pts2_table <- cd4_ranef_pts2_table %>%	mutate(across(c(Q2.5, Estimate, Q97.5, Est.Error, prob_negative), ~ round(.x, 3)))
cd4_ranef_pts2_table <- cd4_ranef_pts2_table %>%	mutate(across(c(p_value_bayes), ~ signif(.x, 2)))
cd4_ranef_pts2_table <- cd4_ranef_pts2_table[order(-cd4_ranef_pts2_table$prob_negative, cd4_ranef_pts2_table$Estimate), ]
write.csv( cd4_ranef_pts2_table, file=glue("{RESULTS_PATH}/tables/tableS9.csv"), quote=F, row.names = F )
saveRDS(cd4_ranef_pts2_table, glue("{RDS_PATH}/cd4_ranef_pts2_table.rds"))

# For table S9 as well, extract residual sd, sds interc and slopes, and corr interc and slope
cd4_ranef_pts2_tables_sds <- readRDS(glue("{RDS_PATH}/bmm_res_ranef_pts.rds"))
cd4_ranef_pts2_tables_sds <- do.call(rbind, lapply(1:4, function(j) { 
	ss <- cd4_ranef_pts2_tables_sds[[1, j]]$sigma_sds 
	ss$subtype <- tree_names[j]
	ss$group2 <- NA
	ss2 <- cd4_ranef_pts2_tables_sds[[1, j]]$p_cors
	ss2$subtype <- tree_names[j]
	ss2$group2 <- ss2$group
	ss2$group <- NULL
	colnames(ss2) <- c("Param","Estimate","Q2.5","Q97.5","subtype","group2" )
	ss3 <- bind_rows(ss, ss2)
	ss3
}))
cd4_ranef_pts2_tables_sds <- cd4_ranef_pts2_tables_sds %>% dplyr::select(subtype, Param, Q2.5, Estimate, Q97.5, group2)
cd4_ranef_pts2_tables_sds <- cd4_ranef_pts2_tables_sds %>%	mutate(across(c(Q2.5, Estimate, Q97.5), ~ round(.x, 3)))
write.csv( cd4_ranef_pts2_tables_sds, file=glue("{RESULTS_PATH}/tables/tableS9_2.csv"), quote=F, row.names = T )

# extract the ones with post prob of being negative slope > 0.8
cd4_ranef_pts2_suspVOIs <- cd4_ranef_pts2[cd4_ranef_pts2$prob_negative > 0.8,] #20

# load ML random eff 0.1.1 model to compare
cd4_ml_rand_eff_table <- readRDS(glue("{RDS_PATH}/cd4_ml_rand_eff_table.rds"))
cd4_ml_rand_eff_table_p <- cd4_ml_rand_eff_table[cd4_ml_rand_eff_table$p_value <= 0.05 & cd4_ml_rand_eff_table$phylotype_coef < 0,] #22
# inner join 
overlap_ml_bay_cd4_randeff <- inner_join(cd4_ranef_pts2_suspVOIs, cd4_ml_rand_eff_table_p, by=c("subtype","phylotype")) #10
saveRDS(overlap_ml_bay_cd4_randeff, glue("{RDS_PATH}/overlap_ml_bay_cd4_randeff.rds"))
overlap_ml_bay_cd4_randeff[,c("subtype","phylotype")]
# subtype phylotype
# 1     A_A1         3
# 2        B       133 *
# 3        B        90
# 4        B       118
# 5        B        24
# 6        B       137
# 7     A_A1         8
# 8        B        69 *
# 9        B        84
# 10       B        62

# in ML but not in Bayes
# 710 CRF_02_AG         7
# 77          B        77