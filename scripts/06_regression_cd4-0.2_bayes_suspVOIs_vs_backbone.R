libs_load <- c("ggplot2","dplyr","ggpubr","data.table","glue","lubridate","sjPlot","viridis", "lme4", "performance", "brms","bayesplot",  "ape", "forcats","pbmcapply") # "rstan"
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"
system(glue("mkdir -p {RESULTS_PATH}/06_cd4_bayes_model"))

min_cl_size_choices <- c(30, 50, 100)
tree_names <- c("A_A1","CRF_02_AG","C","B")

# Bayesian CD4 model
residuals_removal_mx <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))
backbone_cl_control <- readRDS(glue("{RDS_PATH}/backbone_cl_control.rds"))

# After identyfing outliers using random effects model, fit individual fixed eff model for each VOI (consider the n=10 from last ML CD4 analysis), non-VOI & backbone to get effect size
fit_bayes_cd4_model_fixed_eff_vois_vs_backbone <- function(outcome, fostr, warmup_, cd4_df, prior_phylotype_eff="flat", prior_val=NULL) { #priors_slope_interc,pt_id
	
	formula <- as.formula(paste(outcome, fostr))
	
	gp <- get_prior(formula, data=cd4_df, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
	#print(gp)
	if(prior_phylotype_eff=="flat") {
		bayes_fit <- brm(formula, data=cd4_df, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
	} else if(prior_phylotype_eff=="normal_conservative") {
		gp$prior[ grepl(gp$coef, pattern = 'years_since_1cd4:phylotype') ] <- glue("normal(0,20)")
		gp$prior[ grepl(gp$coef, pattern = '^phylotype') ] <- glue("normal(0,20)")
		print(gp)
		bayes_fit <- brm(formula, data=cd4_df, prior=gp, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
	} else if(prior_phylotype_eff=="r2d2") {
		print("R2D2")
		bayes_fit <- brm(formula, data=cd4_df, prior=prior_val, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
		#gp <- get_prior(formula, data=cd4_df, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
 } else {
		stop("Not a valid prior option")
	}
	print(gp)
	
	loo_model <- loo(bayes_fit)
	print("========")
	print("LOO:")
	print(loo_model)
	print("========")
	
	summary_fit <- summary(bayes_fit)
	rhat_values <- summary_fit$rhats
	print("head rhats")
	print(head(rhat_values))
	ess_values <- summary_fit$ess_bulk
	print("head ess")
	print(ess_values)
	
	if (length(which(rhat_values > 1.01)) == 0) {
		print("All Rhat values are <= 1.01")
	}
	if (length(which(ess_values < 400)) == 0) {
		print("All effective sample sizes are >= 400")
	}
	
	gc()
	
	p_all_rows <- fixef(bayes_fit)
	p_all_rows <- as.data.frame(p_all_rows)
	p_all_rows$phylotype <- rownames(p_all_rows)
	rownames(p_all_rows) <- NULL 
	p_all_rows <- p_all_rows %>% dplyr::select(phylotype, Estimate, Est.Error, Q2.5, Q97.5)
	p_pt_only <- p_all_rows
	p_pt_only <- p_pt_only[ grepl(p_pt_only$phylotype, pattern = 'years_since_1cd4:phylotype'), ]
	p_pt_only$phylotype <- gsub("years_since_1cd4:phylotype", "", as.character(p_pt_only$phylotype))
	
	p_pt_intercepts <- p_all_rows[grepl(p_all_rows$phylotype, pattern = '^phylotype'),]
	
	sigma_df <- as.data.frame(summary_fit$spec_pars)
	sigma_df <- sigma_df["sigma", c("Estimate", "l-95% CI", "u-95% CI")]
	sigma_df$param <- "sigma"; sigma_df$group <- "residual"
	
	patient_sd_cor <- as.data.frame(summary_fit$random$patientindex)
	
	# Get posterior draws
	posterior_draws <- as_draws_df(bayes_fit)
	
	# Find fixed effect interaction terms between phylotype and years_since_1cd4
	slope_cols <- grep("b_years_since_1cd4:phylotype", names(posterior_draws), value = TRUE)
	
	# Calculate Bayesian p-values and classify based on posterior probability
	.bayes_p_vals <- lapply(slope_cols, function(col) {
		values <- posterior_draws[[col]]
		
		# Posterior probabilities
		p_neg <- mean(values < 0)
		p_pos <- mean(values > 0)
		
		# Two-sided Bayesian p-value
		p_val <- 2 * min(p_neg, p_pos)
		
		# Extract phylotype name from the column name
		phylotype_name <- gsub("b_years_since_1cd4:phylotype", "\\1", col)
		
		data.frame(phylotype = phylotype_name, p_value_bayes = p_val, prob_negative = p_neg)
	})
	
	bayes_p_df <- do.call(rbind, .bayes_p_vals)
	
	p_pt_only <- merge(p_pt_only, bayes_p_df, by = "phylotype", all.x = TRUE)
	
	return(list(all_covar=p_all_rows, pt_effs=p_pt_only, p_pt_intercepts=p_pt_intercepts, sigma_df=sigma_df, patient_sd_cor=patient_sd_cor, loo_model=loo_model)) # df_p=df_p
}

empirical_bayes_r2d2_prior <- function(model) {
	# Calculate R2 values
	r2p <- performance::r2_nakagawa(model)
	print(r2p)
	
	# Check if R2_marginal exists and is not NA
	if (!("R2_marginal" %in% names(r2p)) || is.na(r2p$R2_marginal)) {
		warning("R2_marginal not found or is NA. Using default value of 0.3")
		r2_marg <- 0.63  # Default fallback value
	} else {
		r2_marg <- unname(r2p$R2_marginal)
	}
	
	print(paste("Using R2 marginal value:", r2_marg))
	
	# Create the prior specification string - THIS IS THE KEY CHANGE
	prior_string <- paste0("R2D2(mean_R2 = ", r2_marg, ", prec_R2 = 2, cons_D2 = 0.5, autoscale = TRUE)")
	
	# Create and return the prior
	bm1prior <- brms::set_prior(prior_string, class = "b")
	print(bm1prior)
	return(bm1prior)
}

ITER <- 5000 #10000
WARMUP40 <- ITER/2.5
WARMUP10 <- ITER/10
WARMUP50 <- ITER/2
THIN <- 5 #10
CORES_CHAINS <- 2

cd4_df_filt <- do.call(rbind, lapply(1:4, function(j) { 
	ss <- residuals_removal_mx[[1, j]]$cd4_df_filt
	ss$subtype <- tree_names[j]
	ss
}))

backbone_cl_control2 <- do.call(rbind, lapply(1:4, function(j) { 
	ss <- data.frame(subtype=tree_names[j], backbone=backbone_cl_control[[1, j]])
	ss
}))

# suspected VOIs from ML random effect model]
cd4_ml_rand_eff_table <- readRDS(glue("{RDS_PATH}/cd4_ml_rand_eff_table.rds"))
cd4_ml_rand_eff_table_p <- cd4_ml_rand_eff_table[cd4_ml_rand_eff_table$p_value <= 0.05 & cd4_ml_rand_eff_table$phylotype_coef < 0,] #12
cd4_ml_rand_eff_table_p_min <- cd4_ml_rand_eff_table_p %>% dplyr::select(subtype, phylotype)

cd4_model_form_with_intercept_pt <- " ~ years_since_1cd4 + age_group + years_since_1cd4 * sexid + years_since_1cd4:exposureid + years_since_1cd4*phylotype + (years_since_1cd4 | patientindex)"
d_bmm_fixedeff_vois <- d_bmm_fixedeff_vois_ref_bb <- bmm_res_fixeff_vois_vs_bb <- list()
for(i in 1:nrow(cd4_ml_rand_eff_table_p_min)) { #intersect_cd4_susp_VOIs
	print(glue("n={i}"))
	
	target_subtype <- cd4_ml_rand_eff_table_p_min$subtype[i]
	print(target_subtype)
	target_phylotype <- cd4_ml_rand_eff_table_p_min$phylotype[i]
	print(target_phylotype)
	backbone_phylotype <- backbone_cl_control2$backbone[backbone_cl_control2$subtype==target_subtype]
	print(backbone_phylotype)
	
	# Subset from cd4_df where both subtype and phylotype match
	match_pair <- cd4_df_filt[cd4_df_filt$subtype == target_subtype & cd4_df_filt$phylotype == target_phylotype, ]
	print("nrow match pair")
	print(nrow(match_pair))
	
	# Subset from cd4_df where only phylotype matches backbone
	match_backbone <- cd4_df_filt[cd4_df_filt$phylotype == backbone_phylotype, ]
	print("nrow match_backbone")
	print(nrow(match_backbone))
	
	# Combine both subsets
	d_bmm_fixedeff_vois[[i]] <- rbind(match_pair, match_backbone)
	d_bmm_fixedeff_vois_ref_bb[[i]] <- d_bmm_fixedeff_vois[[i]] %>% mutate(phylotype = relevel(phylotype, ref=backbone_phylotype), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )

	bmm_res_fixeff_vois_vs_bb[[i]] <- fit_bayes_cd4_model_fixed_eff_vois_vs_backbone("cd4",cd4_model_form_with_intercept_pt,WARMUP40, d_bmm_fixedeff_vois_ref_bb[[i]], prior_phylotype_eff="normal_conservative", prior_val = NULL ) #flat
	print(bmm_res_fixeff_vois_vs_bb[[i]]$pt_effs)
}

saveRDS(bmm_res_fixeff_vois_vs_bb, glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb_NORMAL_P_CONS.rds"))

# summarise phylotype slope effects
cd4_fixeff_suspVOIs_vs_bb <- do.call(rbind, lapply(bmm_res_fixeff_vois_vs_bb, function(x) x$pt_effs))
cd4_fixeff_suspVOIs_vs_bb_table <- cd4_fixeff_suspVOIs_vs_bb[order(cd4_fixeff_suspVOIs_vs_bb$p_value_bayes), ]
# add subtype back
cd4_fixeff_suspVOIs_vs_bb_table$subtype <- cd4_ml_rand_eff_table_p_min$subtype[
	match(cd4_fixeff_suspVOIs_vs_bb_table$phylotype, cd4_ml_rand_eff_table_p_min$phylotype)
]
cd4_fixeff_suspVOIs_vs_bb_table <- cd4_fixeff_suspVOIs_vs_bb_table %>% dplyr::select(subtype, phylotype, Q2.5, Estimate, Q97.5, Est.Error, p_value_bayes, prob_negative) #evidence_label
cd4_fixeff_suspVOIs_vs_bb_table <- cd4_fixeff_suspVOIs_vs_bb_table %>%	mutate(across(c(Q2.5, Estimate, Q97.5, Est.Error, prob_negative), ~ round(.x, 3)))
cd4_fixeff_suspVOIs_vs_bb_table <- cd4_fixeff_suspVOIs_vs_bb_table %>%	mutate(across(c(p_value_bayes), ~ signif(.x, 2)))
cd4_fixeff_suspVOIs_vs_bb_table <- cd4_fixeff_suspVOIs_vs_bb_table[order(cd4_fixeff_suspVOIs_vs_bb_table$p_value_bayes, cd4_fixeff_suspVOIs_vs_bb_table$Estimate),]
write.csv( cd4_fixeff_suspVOIs_vs_bb_table, file=glue("{RESULTS_PATH}/tables/tableS11.csv"), quote=F, row.names = F )
saveRDS(cd4_fixeff_suspVOIs_vs_bb_table, glue("{RDS_PATH}/cd4_fixeff_suspVOIs_vs_bb_table.rds"))

# For table S11 as extract residual sd, sds interc and slopes, and corr interc and slope
cd4_fixeff_suspVOIs_vs_bb_sds <- bmm_res_fixeff_vois_vs_bb #bmm_res_fixeff_vois_vs_bb
cd4_fixeff_suspVOIs_vs_bb_sds <- do.call(rbind, lapply(seq_along(cd4_fixeff_suspVOIs_vs_bb_sds), function(i) {
	x <- cd4_fixeff_suspVOIs_vs_bb_sds[[i]]
	ss <- x$sigma_df
	ss$param <- rownames(ss); rownames(ss) <- NULL; ss$Est.Error <- NULL; ss$param <- NULL; ss$group <- NULL
	ss$phylotype <- x$pt_effs$phylotype
	ss2 <- x$patient_sd_cor
	ss2$param <- rownames(ss2); rownames(ss2) <- NULL; ss2$Rhat <- NULL; ss2$Bulk_ESS <- NULL; ss2$Tail_ESS <- NULL
	ss2$phylotype <- x$pt_effs$phylotype
	ss3 <- dplyr::bind_rows(ss, ss2)
	ss3
}))
cd4_fixeff_suspVOIs_vs_bb_sds <- cd4_fixeff_suspVOIs_vs_bb_sds %>% dplyr::select(phylotype, param, `l-95% CI`, Estimate, `u-95% CI`)
cd4_fixeff_suspVOIs_vs_bb_sds <- cd4_fixeff_suspVOIs_vs_bb_sds %>%	mutate(across(c(`l-95% CI`, Estimate, `u-95% CI`), ~ round(.x, 3)))
cd4_fixeff_suspVOIs_vs_bb_sds$subtype <- cd4_ml_rand_eff_table_p_min$subtype[
	match(cd4_fixeff_suspVOIs_vs_bb_sds$phylotype, cd4_ml_rand_eff_table_p_min$phylotype)
]
# mean and sd across models (n=10) separated by subtype
cd4_fixeff_suspVOIs_vs_bb_sds2 <- cd4_fixeff_suspVOIs_vs_bb_sds %>%
	group_by(subtype, param) %>% 
	summarise(mean = round(mean(Estimate, na.rm = TRUE), 3),sd   = round(sd(Estimate, na.rm = TRUE), 3),.groups = "drop")
write.table( cd4_fixeff_suspVOIs_vs_bb_sds2, file=glue("{RESULTS_PATH}/tables/tableS11_2.tsv"), sep = "\t", quote=F, row.names = F )

# Table 2
# vl_t-test_adj_fdr (one-tailed test)
vl_t_test_adj_fdrs <- readRDS(glue("{RDS_PATH}/vl_st_test.rds"))
vl_t_test_adj_fdrs[[1,1]]$all_adj$subtype <- tree_names[1]; vl_t_test_adj_fdrs[[1,2]]$all_adj$subtype <- tree_names[2]
vl_t_test_adj_fdrs[[1,3]]$all_adj$subtype <- tree_names[3]; vl_t_test_adj_fdrs[[1,4]]$all_adj$subtype <- tree_names[4]
vl_t_test_adj_fdr <- rbind(vl_t_test_adj_fdrs[[1,1]]$all_adj, vl_t_test_adj_fdrs[[1,2]]$all_adj, vl_t_test_adj_fdrs[[1,3]]$all_adj, vl_t_test_adj_fdrs[[1,4]]$all_adj)
vl_t_test_adj_fdr$vl_ttest_adj_fdr_estim <- vl_t_test_adj_fdr$estimate_x
vl_t_test_adj_fdr$vl_ttest_adj_fdr_p <- vl_t_test_adj_fdr$p_adj_fdr
vl_t_test_adj_fdr$p_value1 <- signif(vl_t_test_adj_fdr$vl_ttest_adj_fdr_p, digits = 2)
vl_t_test_adj_fdr <- vl_t_test_adj_fdr %>% dplyr::select(subtype, phylotype, vl_ttest_adj_fdr_estim, p_value1, estimate_y) %>% arrange(p_value1) #156 rows (<10 patients not computed)
vl_t_test_adj_fdr$phylotype <- as.character(vl_t_test_adj_fdr$phylotype)
vl_t_test_adj_fdr_sig <- vl_t_test_adj_fdr[vl_t_test_adj_fdr$p_value1 <= 0.05 & vl_t_test_adj_fdr$vl_ttest_adj_fdr_estim >= vl_t_test_adj_fdr$estimate_y,] #9
vl_t_test_adj_fdr_sig <- vl_t_test_adj_fdr_sig %>% dplyr::select(subtype, phylotype, vl_ttest_adj_fdr_estim, p_value1)
colnames(vl_t_test_adj_fdr_sig) <- c("subtype", "phylotype_id","vl_ttest_adj_fdr_estim","p_value1") 
vl_t_test_adj_fdr <- vl_t_test_adj_fdr  %>% dplyr::select(subtype, phylotype, vl_ttest_adj_fdr_estim, p_value1)
colnames(vl_t_test_adj_fdr) <- c("subtype", "phylotype_id","vl_ttest_adj_fdr_estim","p_value1") 

# vl_bayes_joint (two-tailed test; should make one-tailed to be consistent with the ones below?)
vl_bayes_joints <- readRDS(glue("{RDS_PATH}/res_vl_model.rds"))
vl_bayes_joints[[1,1]]$df_table$subtype <- tree_names[1]; vl_bayes_joints[[1,2]]$df_table$subtype <- tree_names[2]
vl_bayes_joints[[1,3]]$df_table$subtype <- tree_names[3]; vl_bayes_joints[[1,4]]$df_table$subtype <- tree_names[4]
vl_bayes_joint <- rbind(vl_bayes_joints[[1,1]]$df_table, vl_bayes_joints[[1,2]]$df_table, vl_bayes_joints[[1,3]]$df_table, vl_bayes_joints[[1,4]]$df_table)
rownames(vl_bayes_joint) <- NULL
vl_bayes_joint$vl_bayes_joint_estim <- vl_bayes_joint$quantile_0.5
vl_bayes_joint$vl_bayes_joint_p <- vl_bayes_joint$p
vl_bayes_joint$phylotype_id <- vl_bayes_joint$group
vl_bayes_joint$p_value2 <- signif(vl_bayes_joint$vl_bayes_joint_p, digits = 2)
vl_bayes_joint <- vl_bayes_joint %>% dplyr::select(subtype, phylotype_id, vl_bayes_joint_estim, p_value2, y_mean_pop) %>% arrange(p_value2) # 207 rows
# filter only signif p and higher vl
vl_bayes_joint_sig <- vl_bayes_joint[vl_bayes_joint$p_value2 <= 0.05 & vl_bayes_joint$vl_bayes_joint_estim >= vl_bayes_joint$y_mean_pop,] 
vl_bayes_joint_sig <- vl_bayes_joint_sig %>% dplyr::select(subtype, phylotype_id, vl_bayes_joint_estim, p_value2)
vl_bayes_joint <- vl_bayes_joint %>% dplyr::select(subtype, phylotype_id, vl_bayes_joint_estim, p_value2)
colnames(vl_bayes_joint) <- c("subtype", "phylotype_id","vl_bayes_joint_estim","p_value2") 

# vl (sample sizes): for bayes model (for t-test if pt has <10 patients they are removed and not computed)
vl_ss <- readRDS(glue("{RDS_PATH}/vl_subtypes_comb2.rds"))[1,]
names(vl_ss) <- tree_names
vl_ss <- rbindlist(vl_ss, idcol = "subtype")
vl_ss <- vl_ss %>% group_by(subtype, cluster) %>% summarise(vl_npatients=n())
colnames(vl_ss) <- c("subtype", "phylotype_id","vl_npatients") #205

# cd4_ml_random_eff_ml (one-tailed test)
lmm2 <- readRDS(glue("{RDS_PATH}/lmm2.rds"))[1,] # mcs=30
lmm2 <- lapply(lmm2, `[[`, 3)
names(lmm2) <- tree_names
cd4_ml_rand_eff <- rbindlist(lmm2, idcol = "subtype")
rownames(cd4_ml_rand_eff) <- NULL
cd4_ml_rand_eff$cd4_ml_random_eff_estim <- cd4_ml_rand_eff$phylotype_coef
cd4_ml_rand_eff$cd4_ml_random_eff_estim_p <- cd4_ml_rand_eff$p_value
cd4_ml_rand_eff$p_value3 <- signif(cd4_ml_rand_eff$cd4_ml_random_eff_estim_p, digits = 2)
cd4_ml_rand_eff <- cd4_ml_rand_eff %>% dplyr::select(subtype, phylotype, cd4_ml_random_eff_estim, p_value3) %>% arrange(p_value3) #cd4_ml_random_eff_estim
colnames(cd4_ml_rand_eff) <- c("subtype", "phylotype_id","cd4_ml_random_eff_estim","p_value3")
cd4_ml_rand_eff_sig <- cd4_ml_rand_eff[cd4_ml_rand_eff$p_value3 <= 0.05 & cd4_ml_rand_eff$cd4_ml_random_eff_estim < 0,] #12
colnames(cd4_ml_rand_eff_sig) <- c("subtype", "phylotype_id","cd4_ml_random_eff_estim","p_value3") # 204 rows

bmm_res_fixeff_vois_vs_bb <- readRDS(glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb.rds"))

# cd4_bayes_fixed_eff_vois_vs_backone 
cd4_bayes_fixed_eff <- readRDS(glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb_NORMAL_P_CONS.rds"))
cd4_bayes_fixed_eff <- do.call(rbind, lapply(cd4_bayes_fixed_eff, function(x) x$pt_effs))
cd4_bayes_fixed_eff$subtype <- cd4_ml_rand_eff_table_p_min$subtype[
	match(cd4_bayes_fixed_eff$phylotype, cd4_ml_rand_eff_table_p_min$phylotype)
]
cd4_bayes_fixed_eff$cd4_bayes_fixed_eff_estim <- cd4_bayes_fixed_eff$Estimate
cd4_bayes_fixed_eff$p_value4 <- signif(cd4_bayes_fixed_eff$p_value, digits = 2)
cd4_bayes_fixed_eff <- cd4_bayes_fixed_eff %>% dplyr::select(subtype, phylotype, cd4_bayes_fixed_eff_estim, p_value4) %>% arrange(p_value4) #cd4_bayes_fixed_eff
colnames(cd4_bayes_fixed_eff) <- c("subtype", "phylotype_id","cd4_bayes_fixed_eff_estim","p_value4") #12 rows
cd4_bayes_fixed_eff_sig <- cd4_bayes_fixed_eff[cd4_bayes_fixed_eff$p_value <= 0.05,] # 6 rows
colnames(cd4_bayes_fixed_eff_sig) <- c("subtype", "phylotype_id","cd4_bayes_fixed_eff_estim","p_value4")

# cd4 (sample size)
cd4_ss <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))[1,] # mcs=30
cd4_ss <- lapply(cd4_ss, `[[`, 3)
names(cd4_ss) <- tree_names
cd4_ss <- rbindlist(cd4_ss, idcol = "subtype")
cd4_ss <- cd4_ss %>% group_by(subtype, cluster) %>% summarise(cd4_npatients=length(unique(patientindex)))
colnames(cd4_ss) <- c("subtype", "phylotype_id","cd4_npatients") # 204 rows

# merge all
dfs_to_merge <- list(vl_t_test_adj_fdr, vl_bayes_joint, vl_ss, cd4_ml_rand_eff, cd4_bayes_fixed_eff, cd4_ss)
combined_vl_cd4_estims <- Reduce(function(x, y) merge(x, y, by = c("subtype","phylotype_id"), all = TRUE), dfs_to_merge)
combined_vl_cd4_estims$vl_ttest_adj_fdr_estim <- round(combined_vl_cd4_estims$vl_ttest_adj_fdr_estim, 3)
combined_vl_cd4_estims$vl_bayes_joint_estim <- round(combined_vl_cd4_estims$vl_bayes_joint_estim, 3)
combined_vl_cd4_estims$cd4_ml_random_eff_estim <- round(combined_vl_cd4_estims$cd4_ml_random_eff_estim, 3)
combined_vl_cd4_estims$cd4_bayes_fixed_eff_estim <- round(combined_vl_cd4_estims$cd4_bayes_fixed_eff_estim, 3)
saveRDS(combined_vl_cd4_estims, glue("{RDS_PATH}/combined_vl_cd4_estims.rds"))
# order vl bayes
combined_vl_cd4_estims1 <- combined_vl_cd4_estims %>% arrange( p_value2 )
# order cd4 ml random eff
combined_vl_cd4_estims3 <- combined_vl_cd4_estims %>% arrange( p_value3 )
# order cd4 bayes fixed eff
combined_vl_cd4_estims4 <- combined_vl_cd4_estims %>% arrange( p_value4 )

# order by subtype (B, C, A1, CRF)
order_subtypes_show <- c("B","C","A_A1", "CRF_02_AG")
combined_vl_cd4_estims$subtype <- factor(combined_vl_cd4_estims$subtype, levels=order_subtypes_show, labels=order_subtypes_show)
combined_vl_cd4_estims$phylotype_id <- fct_relevel(combined_vl_cd4_estims$phylotype_id,function(x){as.character(sort(as.integer(x)))})
combined_vl_cd4_estims_f <- combined_vl_cd4_estims %>% arrange(subtype,phylotype_id)
write.csv(combined_vl_cd4_estims_f, glue("{RESULTS_PATH}/tables/table2_all_subtypeB_pts.csv"), quote = F, row.names = F, na = "") # 207 rows
# only show the ones that are signif in at least one analysis
combined_vl_cd4_estims_f2 <- combined_vl_cd4_estims_f %>% filter( (phylotype_id %in% vl_bayes_joint_sig$phylotype_id & subtype %in% vl_bayes_joint_sig$subtype) |
																																																																		(phylotype_id %in% vl_t_test_adj_fdr_sig$phylotype_id & subtype %in% vl_t_test_adj_fdr_sig$subtype) |
																																																																		(phylotype_id %in% cd4_ml_rand_eff_sig$phylotype_id & subtype %in% cd4_ml_rand_eff_sig$subtype) |
																																																																		(phylotype_id %in% cd4_bayes_fixed_eff_sig$phylotype_id & subtype %in% cd4_bayes_fixed_eff_sig$subtype))
write.csv(combined_vl_cd4_estims_f2, glue("{RESULTS_PATH}/tables/table2.csv"), quote = F, row.names = F, na = "") # 30 rows

# Shortlist
combined_vl_cd4_estims <- as.data.frame(combined_vl_cd4_estims)
i0 = with( combined_vl_cd4_estims, (is.na(p_value1) | p_value1 <.05) | (is.na(p_value2) | p_value2 <.05) | (is.na(p_value3) | p_value3 <.05) | (is.na(p_value4) | p_value4 <.05))
pvnames <- c('p_value1', 'p_value2', 'p_value3', 'p_value4')
i0 = apply( combined_vl_cd4_estims[,pvnames]|>as.matrix(), MAR=1, FUN=function(x) {
	if( all( is.na(x))) return(FALSE)
	( min( na.omit(x)) < .05)
})
d0 <- combined_vl_cd4_estims[ i0, ] #nrow(d0): 34
d0

i1 <- combined_vl_cd4_estims$vl_bayes_joint_estim > median(combined_vl_cd4_estims$vl_bayes_joint_estim )
d1 <- combined_vl_cd4_estims[ i0 & i1, ]  #nrow(d1): 19

cor( d1[ , pvnames], use='na.or.complete')

with( d1, plot( p_value1, p_value3, log='xy' ))
j0 <- with( d1, p_value1<.05 & p_value3<.05 )
h0 <- d1[j0, ] #nrow(h0): 2 (69, 137)

j1 <- with( d1, p_value2<.05 & p_value3<.05 )
h1 <- d1[j1, ] #nrow(h1): 1 (69)

j2 <- with( d1, p_value2<.05 & p_value4<.05 )
h2 <- d1[j2, ] #nrow(h1): 1 (69)

# VOI 1 = PT69 because significant in all tests

# order vl_ttest by increasing p
vl_ttest_sl <- combined_vl_cd4_estims %>% filter( p_value1<.05 ) %>% arrange(p_value1)
vl_ttest_sl$phylotype_id # VOI 2: PT40 (smallest p and higher vl estimate) 20  79 125 137	69  14  18	4

# order vl_bayes by increasing p
vl_bayes_sl <- combined_vl_cd4_estims %>% filter( vl_bayes_joint_estim > median(vl_bayes_joint_estim ) & p_value2<.05 ) %>% arrange(p_value2)
vl_bayes_sl$phylotype_id # VOI 2: 40 (smallest p and higher vl estimate), 101, 20, 69

# order cd4_bayes_fixed_eff by increasing p
cd4_bayes_fe_sl <- combined_vl_cd4_estims %>% filter( p_value4<.05 ) %>% arrange(p_value4)
cd4_bayes_fe_sl$phylotype_id # 133 (smallest p and fastest estim/decline among B) 3	8	90 118	69
cd4_bayes_fe_sl %>% select(subtype, phylotype_id, cd4_bayes_fixed_eff_estim)
# subtype phylotype_id cd4_bayes_fixed_eff_estim
# subtype phylotype_id cd4_bayes_fixed_eff_estim
# 1       B          133                   -31.861
# 2    A_A1            3                   -18.771
# 3    A_A1            8                   -22.834
# 4       B           90                   -25.673
# 5       B          118                   -21.484
# 6       B           69                   -22.320

# order cd4_bayes_rand_eff by increasing p
cd4_bayes_re_sl <- combined_vl_cd4_estims %>% filter( cd4_ml_random_eff_estim < 0 & p_value3<.05 ) %>% arrange(p_value3)
cd4_bayes_re_sl$phylotype_id #133 (smallest p and fastest estim/decline among B)  69  90 118 7	3	8	137  84  24  77
cd4_bayes_re_sl %>% select(subtype, phylotype_id, cd4_ml_random_eff_estim, p_value3)
# subtype phylotype_id cd4_ml_random_eff_estim p_value3
# 1          B          133                 -10.182  0.00004
# 2          B           69                  -7.978  0.00240
# 3          B           90                  -7.987  0.00330
# 4          B          118                  -8.773  0.00390
# 5  CRF_02_AG            7                  -9.919  0.00600
# 6       A_A1            3                 -15.564  0.01900
# 7       A_A1            8                  -7.336  0.01900
# 8          B          137                  -8.610  0.02000
# 9          B           84                  -8.173  0.03000
# 10         B           24                  -7.638  0.03300
# 11         B           77                  -3.687  0.03400

# Table S13: coeffs of covariates for VOIs in fixed effects model
cd4_bayes_fixed_eff_tables <- readRDS(glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb_NORMAL_P_CONS.rds"))

#cd4_bayes_fixed_eff_table2 <- do.call(rbind, lapply(cd4_bayes_fixed_eff_tables, function(x) data.frame(x$all_covar, mcs=x$pt_effs$mcs, subtype=x$pt_effs$subtype))) #322
cd4_bayes_fixed_eff_table2 <- do.call(rbind, lapply(seq_along(cd4_bayes_fixed_eff_tables), function(i) {
	x <- cd4_bayes_fixed_eff_tables[[i]]
	data.frame(
		x$all_covar,
		subtype = cd4_ml_rand_eff_table_p_min[i,"subtype"],
		phylotype_id = cd4_ml_rand_eff_table_p_min[i,"phylotype"]
	)
}))

cd4_bayes_fixed_eff_table2.1 <- cd4_bayes_fixed_eff_table2[!(cd4_bayes_fixed_eff_table2$phylotype %in%  c("Intercept","years_since_1cd4","years_since_1cd4:exposureidOther")),] #132, s"years_since_1cd4:phylotypenonMVOI"
# remove PTs with non-sig slope effect from these tables
ns_pts <- cd4_bayes_fixed_eff[cd4_bayes_fixed_eff$p_value > 0.05,]
cd4_bayes_fixed_eff_table2.2 <- cd4_bayes_fixed_eff_table2.1[!(cd4_bayes_fixed_eff_table2.1$phylotype_id %in% ns_pts$phylotype_id),]
cd4_bayes_fixed_eff_table2.2 <- cd4_bayes_fixed_eff_table2.2 %>% arrange(Estimate)
cd4_bayes_fixed_eff_table2.2 <- cd4_bayes_fixed_eff_table2.2 %>% dplyr::select( subtype, phylotype_id, phylotype, Q2.5, Estimate, Q97.5, Est.Error )
colnames(cd4_bayes_fixed_eff_table2.2) <- c("subtype", "phylotype", "parameter", "Q2.5", "Estimate", "Q97.5", "Est.Error")
cd4_bayes_fixed_eff_table2.2 <- cd4_bayes_fixed_eff_table2.2 %>% mutate(across(c(Q2.5, Estimate, Q97.5, Est.Error), ~ round(.x, 3)))
write.csv( cd4_bayes_fixed_eff_table2.2, file=glue("{RESULTS_PATH}/tables/tableS13.csv"), quote=F, row.names = F )

# median effects by age_group60P (intercept)
cd4_bayes_fixed_eff_table2.2_60_older <- cd4_bayes_fixed_eff_table2.2 %>% filter(parameter == "age_group60P") %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))
# median effects by age_group50M59 (intercept)
cd4_bayes_fixed_eff_table2.2_50s <- cd4_bayes_fixed_eff_table2.2 %>% filter(parameter == "age_group50M59") %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))
# median effects by females (intercept)
cd4_bayes_fixed_eff_table2.2_females <- cd4_bayes_fixed_eff_table2.2 %>% filter(parameter == "sexidFemale") %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))
# median effects by age_group40M49 (intercept)
cd4_bayes_fixed_eff_table2.2_40s <- cd4_bayes_fixed_eff_table2.2 %>% filter(parameter == "age_group40M49") %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))
# median effects of phylotype (intercept)
cd4_bayes_fixed_eff_table2.2_pt_interc <- cd4_bayes_fixed_eff_table2.2 %>% filter(grepl("^phylotype", parameter)) %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))

# Fig. 2B: Plot average CD4 trajectory of VOIs
# IMPORTANT: even though PT40 is not a "CD4" VOI, run Bayes fixed effects model on it to include in Fig. and provide comparison against the "2 CD4 VOIs"
residuals_removal_mx <- readRDS( glue("{RDS_PATH}/residuals_removal_mx.rds") )
d_pt40_run <- residuals_removal_mx[[ 1,4 ]]$cd4_df_filt
d_pt40_run$phylotype <- ifelse( d_pt40_run$cluster ==  backbone_cl_control[[ 1,4 ]], yes=backbone_cl_control[[ 1,4 ]],no=
																																	ifelse( d_pt40_run$cluster == "40", yes= "40", no=as.character("non-VOI") ))
d_pt40_run <- d_pt40_run[d_pt40_run$phylotype != "non-VOI",]
d_pt40_run$phylotype <- as.factor(d_pt40_run$phylotype)
print(table(d_pt40_run$phylotype))
d_pt40_run_ref_bb <- d_pt40_run %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[ 1,4 ]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
# run bayes
bmm_res_pt40_run_ref_bb <- fit_bayes_cd4_model_fixed_eff_vois_vs_backbone("cd4", cd4_model_form_with_intercept_pt, WARMUP40, d_pt40_run_ref_bb, prior_phylotype_eff="normal_conservative", prior_val = NULL)
saveRDS(bmm_res_pt40_run_ref_bb, glue("{RDS_PATH}/bmm_res_pt40_run_ref_bb.rds"))

bmm_res_pt40_run_ref_bb$pt_effs$subtype <- "B"
#bmm_res_pt40_run_ref_bb$pt_effs <- bmm_res_pt40_run_ref_bb$pt_effs[ bmm_res_pt40_run_ref_bb$pt_effs$phylotype_coef != "nonMVOI", ]
pt40_other_coeffs <- bmm_res_pt40_run_ref_bb$all_covar
pt40_other_coeffs$subtype <- "B"; pt40_other_coeffs$phylotype_id <- "40"
colnames(pt40_other_coeffs) <- c("param","Estimate","Est.Error","Q2.5","Q97.5","subtype","phylotype_id")

bmm_vois_combined <- readRDS(glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb_NORMAL_P_CONS.rds"))
bmm_vois_combined2 <- do.call(rbind, lapply(seq_along(bmm_vois_combined), function(i) {
	bmm_vois_combined[[i]]$pt_effs
	#data.frame(x$pt_effs, phylotype_id = cd4_ml_rand_eff_table_p_min[i,"phylotype"])
}))

bmm_vois_combined2$subtype <- cd4_ml_rand_eff_table_p_min$subtype[
	match(bmm_vois_combined2$phylotype, cd4_ml_rand_eff_table_p_min$phylotype)
]

subtype_b_vois_ids_only <- c(40,69,133)
# bmm_vois_combined2 <- bmm_vois_combined2[bmm_vois_combined2$phylotype_coef != "nonMVOI",]
bmm_vois_combined2 <- bmm_vois_combined2[(bmm_vois_combined2$phylotype %in% subtype_b_vois_ids_only) & bmm_vois_combined2$subtype=="B", ]
bmm_vois_combined3 <- rbind(bmm_vois_combined2, bmm_res_pt40_run_ref_bb$pt_effs)

# Standardise initial (baseline) CD4 at 500 and use Bayesian regression coeffs to represent decline of CD4 counts over time
bmm_vois_combined3$phylotype_id <- bmm_vois_combined3$phylotype
bmm_vois_combined3$phylotype <- paste0("PT.B.",bmm_vois_combined3$phylotype_id,".UK")

colnames(cd4_bayes_fixed_eff_table2) <- c("param","Estimate","Est.Error","Q2.5","Q97.5","subtype","phylotype_id")
cd4_bayes_fixed_eff_table2.a <- rbind(cd4_bayes_fixed_eff_table2, pt40_other_coeffs)
cd4_bayes_fixed_eff_table2.b <- cd4_bayes_fixed_eff_table2.a[ (cd4_bayes_fixed_eff_table2.a$phylotype_id %in% subtype_b_vois_ids_only) & cd4_bayes_fixed_eff_table2.a$subtype == "B", ]
cd4_bayes_fixed_eff_table2.b$phylotype <- paste0("PT.B.",cd4_bayes_fixed_eff_table2.b$phylotype_id,".UK")
# usual rate of decline
rate_decl <- cd4_bayes_fixed_eff_table2.b[ grepl(cd4_bayes_fixed_eff_table2.b$param, pattern = '^years_since_1cd4$'), ]#$Estimate
rate_decl <- rate_decl %>% summarise(phylotype="Backbone", Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5=median(Q97.5))

# add usual decline per year together with phylotype effects
bmm_vois_combined4 <- bmm_vois_combined3
bmm_vois_combined4$Estimate <- rate_decl$Estimate + bmm_vois_combined4$Estimate
bmm_vois_combined4$Q2.5 <- rate_decl$Q2.5 + bmm_vois_combined4$Q2.5
bmm_vois_combined4$Q97.5 <- rate_decl$Estimate + bmm_vois_combined4$Q97.5
bmm_vois_combined4_bb <- bind_rows(bmm_vois_combined4, rate_decl)

# usual intercepts
intercs <- cd4_bayes_fixed_eff_table2.b[ grepl(cd4_bayes_fixed_eff_table2.b$param, pattern = 'Intercept'), ]
intercs <- intercs %>% summarise(phylotype="Backbone", Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5=median(Q97.5))
intercs$Estimate_interc <- intercs$Estimate
intercs$Q2.5_interc <- intercs$Q2.5
intercs$Q97.5_interc <- intercs$Q97.5

# VOI intercepts
bmm_vois_combined_interc <- do.call(rbind, lapply(seq_along(bmm_vois_combined), function(i) {
	x <- bmm_vois_combined[[i]]$p_pt_intercepts
	data.frame(x, phylotype_id = cd4_ml_rand_eff_table_p_min[i,"phylotype"])
}))
bmm_vois_combined_interc <- bmm_vois_combined_interc[(bmm_vois_combined_interc$phylotype_id %in% subtype_b_vois_ids_only), ] #& bmm_vois_combined_interc$mcs == "30" & bmm_vois_combined_interc$subtype=="B"
bmm_res_pt40_run_ref_bb$p_pt_intercepts$phylotype_id <- "40"
bmm_vois_combined_interc <- rbind(bmm_vois_combined_interc, bmm_res_pt40_run_ref_bb$p_pt_intercepts)

# add baseline CD4 at VOIs
bmm_vois_combined_interc$phylotype <- paste0("PT.B.",bmm_vois_combined_interc$phylotype_id,".UK")
bmm_vois_interc <- bmm_vois_combined_interc
bmm_vois_interc$Estimate_interc <- intercs$Estimate + bmm_vois_interc$Estimate
bmm_vois_interc$Q2.5_interc  <- intercs$Q2.5 + bmm_vois_interc$Q2.5
bmm_vois_interc$Q97.5_interc <- intercs$Estimate + bmm_vois_interc$Q97.5
bmm_vois_interc_bb <- bind_rows(bmm_vois_interc, intercs)

# create df dropping year 0.1 by 0.1 units
subtype_b_vois <- c("PT.B.40.UK","PT.B.69.UK","PT.B.133.UK")
new_data <- data.frame(
	years_since_1cd4 = rep(seq(from=0,to=6,by=0.1), times=4), #to=6
	phylotype = rep(c('Backbone', c(subtype_b_vois)), each = 61) #each=61
)

# merge new df with regression coeffs
bmm_vois_sl_interc <- inner_join(bmm_vois_combined4_bb, bmm_vois_interc_bb, by="phylotype")
new_data_coeffs <- left_join(new_data, bmm_vois_sl_interc, by="phylotype")

# when initial CD4 = 500, overlapping CI ribbons
# get CD4 values based on estimated declines for each 0.1 time period (initial CD4 = 500)
new_data_coeffs_cd4s <- new_data_coeffs %>% mutate(cd4_count = Estimate_interc + (Estimate.x) * years_since_1cd4, cd4_lower = Q2.5_interc + (Q2.5.x) * years_since_1cd4, cd4_upper = Q97.5_interc + (Q97.5.x) * years_since_1cd4)
new_data_coeffs_cd4s$phylotype <- factor(new_data_coeffs_cd4s$phylotype, levels = c("Backbone", subtype_b_vois)) #"Backbone",c
saveRDS(new_data_coeffs_cd4s, glue("{RDS_PATH}/new_data_coeffs_cd4s.rds"))

cd4_regr_pt_pal <- c("PT.B.40.UK"="#56B4E9", "PT.B.69.UK"="#D55E00", "PT.B.133.UK"="#009E73", "Backbone"="#999999")

leg <- theme(text=element_text(family="Helvetica"), axis.text=element_text(size=10, color="black"), axis.title=element_text(size=10, color="black", face="bold"))

new_data_coeffs_cd4s_list <- list()
for(i in 1:length(subtype_b_vois)) {
	new_data_coeffs_cd4s_list[[i]] <- new_data_coeffs_cd4s %>% filter(phylotype == "Backbone" | phylotype == subtype_b_vois[i] ) %>% mutate(comp_group=paste0("Backbone vs ", subtype_b_vois[i]))
}

new_data_coeffs_cd4s_list_all <- rbindlist(new_data_coeffs_cd4s_list)
new_data_coeffs_cd4s_list_all$comp_group <- factor(new_data_coeffs_cd4s_list_all$comp_group, levels = c("Backbone vs PT.B.40.UK","Backbone vs PT.B.69.UK","Backbone vs PT.B.133.UK"))

# Plot
f2b <- ggplot(new_data_coeffs_cd4s_list_all, aes(x = years_since_1cd4, y = cd4_count)) +
	geom_line(aes(color = phylotype), size = 0.75) +
	geom_ribbon(aes(ymin = cd4_lower, ymax = cd4_upper, fill = phylotype), alpha = 0.3) +
	scale_color_manual(values = cd4_regr_pt_pal, name = "Phylotype") +
	scale_fill_manual(values = cd4_regr_pt_pal, name = "Phylotype") +
	ylab(bquote(bold("Pre-treatment CD4 count (cells /"~mm^3~")"))) +
	xlab("Years since first CD4") +
	theme_classic() +
	coord_cartesian(xlim = c(0, 6), ylim = c(0, 500), expand = FALSE) +
	scale_x_continuous(breaks = seq(0, 6, 1)) +
	#scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
	theme(
		plot.margin = margin(0.5, 0.8, 0.5, 0.5, "cm"),
		legend.key.size = unit(1, 'cm'),
		legend.key.height = unit(0.3, 'cm'),
		legend.key.width = unit(0.3, 'cm'),
		legend.text = element_text(size = 8),
		legend.title = element_text(size=10),
		legend.position = "top", #c(0.75, 0.20),
		strip.background = element_blank(),
		strip.text = element_blank(),
		panel.spacing = unit(1, "lines")
	) + guides(color = guide_legend(ncol = 2,title.position = "top"), fill = guide_legend(ncol = 2,title.position = "top")) +
	geom_hline(yintercept = 350, linetype = "longdash", color = "black", linewidth = 0.3) +
	facet_wrap(~comp_group, ncol = 1, scales = "fixed") +
	leg

saveRDS(f2b, glue("{RDS_PATH}/f2b.rds"))

# time to reach 350 cells/mm3
new_data_coeffs_cd4s_time_to_350 <- new_data_coeffs_cd4s %>%
	filter(cd4_count <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350[,c("phylotype", "years_since_1cd4")]

new_data_coeffs_cd4s_time_to_350_lower <- new_data_coeffs_cd4s %>%
	filter(cd4_lower <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_lower[,c("phylotype", "years_since_1cd4")]

new_data_coeffs_cd4s_time_to_350_upper <- new_data_coeffs_cd4s %>%
	filter(cd4_upper <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_upper[,c("phylotype", "years_since_1cd4")]

# given that intercept estimates are wide and overlap zero, consider all phylotypes starting at baseline 469 cells/mm3 for
# Fig. S10
baseline_ref_group <- round(intercs$Estimate)
new_data_coeffs_cd4s_common_intercept <- new_data_coeffs
new_data_coeffs_cd4s_common_intercept$Estimate_interc <- baseline_ref_group
#new_data_coeffs_cd4s_common_intercept <- new_data_coeffs %>% mutate(cd4_count = baseline_ref_group + (Estimate.x) * years_since_1cd4, cd4_lower = baseline_ref_group + (Q2.5.x) * years_since_1cd4, cd4_upper = baseline_ref_group + (Q97.5.x) * years_since_1cd4)
new_data_coeffs_cd4s_common_intercept <- new_data_coeffs_cd4s_common_intercept %>% mutate(cd4_count = Estimate_interc + (Estimate.x) * years_since_1cd4, cd4_lower = Estimate_interc + (Q2.5.x) * years_since_1cd4, cd4_upper = Estimate_interc + (Q97.5.x) * years_since_1cd4)
s_ci <- ggplot(new_data_coeffs_cd4s_common_intercept, aes(x = years_since_1cd4, y = cd4_count)) + #fill = phylotype, color = phylotype
	geom_ribbon(aes(ymin = cd4_lower, ymax = cd4_upper, fill = phylotype), alpha = 0.3) + #color = phylotype
	geom_line(aes(color = phylotype), size = 0.75) +
	scale_color_manual(values=cd4_regr_pt_pal, name="Phylotype") + scale_fill_manual(values=cd4_regr_pt_pal, name="Phylotype") +
	ylab(bquote(bold("Pre-treatment CD4 count (cells /"~mm^3~")"))) + xlab("Years since first CD4") + theme_classic() + 
	coord_cartesian(xlim=c(0,6), ylim = c(0,500), expand=FALSE) +
	scale_x_continuous(breaks=seq(from=0, to=6, by=1)) + 
	#scale_y_continuous(limits=c(0,500), breaks=seq(from=0,to=500,by=50)) + coord_cartesian(expand = c(0, 0)) +
	theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'),
							legend.key.width = unit(0.5, 'cm'), legend.text=element_text(size=10), legend.position = c(0.25,0.3)) + leg +
	geom_hline(yintercept=350, linetype="longdash", color="black", linewidth=0.3)
ggsave(s_ci, file=glue("{RESULTS_PATH}/figs/figS10.eps"), device=cairo_ps, dpi=600, width=8, height=6, bg="white")
ggsave(s_ci, file=glue("{RESULTS_PATH}/figs/figS10.jpg"), dpi=600, width=8, height=6, bg="white")

new_data_coeffs_cd4s_time_to_350_comm <- new_data_coeffs_cd4s_common_intercept %>%
	filter(cd4_count <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_comm[,c("phylotype","years_since_1cd4")]

new_data_coeffs_cd4s_time_to_350_lower_comm <- new_data_coeffs_cd4s_common_intercept %>%
	filter(cd4_lower <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_lower_comm[,c("phylotype","years_since_1cd4")]

new_data_coeffs_cd4s_time_to_350_upper_comm <- new_data_coeffs_cd4s_common_intercept %>%
	filter(cd4_upper <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_upper_comm[,c("phylotype","years_since_1cd4")]

# difference in estimates when comparing estimated baseline CD4 and common baseline CD4
new_data_coeffs_cd4s_time_to_350$years_since_1cd4 - new_data_coeffs_cd4s_time_to_350_comm$years_since_1cd4
new_data_coeffs_cd4s_time_to_350_lower$years_since_1cd4 - new_data_coeffs_cd4s_time_to_350_lower_comm$years_since_1cd4
new_data_coeffs_cd4s_time_to_350_upper$years_since_1cd4 - new_data_coeffs_cd4s_time_to_350_upper_comm$years_since_1cd4