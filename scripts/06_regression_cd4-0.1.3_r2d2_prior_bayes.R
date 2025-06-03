libs_load <- c("ggplot2","dplyr","ggpubr","data.table","glue","lme4","lubridate", "performance", "brms")
invisible( lapply(libs_load, library, character.only=TRUE) )

# Sensitivity analysis: Bayesian model with fixed effects on phylotype slopes and R2D2 prior 
# R2D2 shrinkage prior 
# There are some advantages to this approach: 
# 1) No assumption of normality of the effect sizes 
# 2) Reasonable prior for shrinking PT effects towards zero w/o the drawback of assuming normality as in 0.1.2 model 
# 3) The resulting effect sizes are larger and closer to what you got with the individual fixed effect models 
# 4) It still makes use of the lmer model, because the prior is calibrated against the R^2 from that fit. 

min_cl_size_choices <- c(30, 50, 100)
tree_names <- c("A_A1","CRF_02_AG","C","B")

RDS_PATH="rds"
RESULTS_PATH="results"
residuals_removal_mx <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))
backbone_cl_control <- readRDS(glue("{RDS_PATH}/backbone_cl_control.rds"))

fit_bayes_cd4_model_r2d2_prior <- function(cd4_df, r2d2_prior) {
	print(r2d2_prior)
	
	bayes_fit <- brm(cd4 ~ years_since_1cd4 + age_group + years_since_1cd4:sexid
																		+ years_since_1cd4:exposureid + phylotype + years_since_1cd4:phylotype
																		+ (years_since_1cd4 | patientindex), data=cd4_df, prior = r2d2_prior, 
																		iter = ITER_, warmup=WARMUP50_, chains  = 2, thin=THIN_, cores = CORES_CHAINS_, seed=1111)
	
	loo_model <- loo(bayes_fit) #moment_match = TRUE
	print("========")
	print("LOO:")
	print(loo_model)
	print("========")
	
	summary_fit <- summary(bayes_fit)
	rhat_values <- summary_fit$rhats
	ess_values <- summary_fit$ess_bulk
	
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
		
		# Classification
		evidence_label <- dplyr::case_when(
			p_neg > 0.90 ~ "Strong negative slope",
			p_neg > 0.80 ~ "Moderate negative slope",
			p_neg < 0.10 ~ "Strong positive slope",
			p_neg < 0.20 ~ "Moderate positive slope",
			TRUE         ~ "No clear evidence"
		)
		
		# Extract phylotype name from the column name
		phylotype_name <- gsub("b_years_since_1cd4:phylotype", "\\1", col)
		
		data.frame(phylotype = phylotype_name, p_value_bayes = p_val, prob_negative = p_neg, evidence_label = evidence_label)
	})
	
	bayes_p_df <- do.call(rbind, .bayes_p_vals)
	
	p_pt_only <- merge(p_pt_only, bayes_p_df, by = "phylotype", all.x = TRUE)
	
	return(list(all_covar=p_all_rows, pt_effs=p_pt_only, p_pt_intercepts=p_pt_intercepts, sigma_df=sigma_df, patient_sd_cor=patient_sd_cor, loo_model=loo_model))
}

# ML LMM: R2D2 prior is calibrated against the R^2 from this fit that has random eff on phylotype slopes
lmmfo  <- "cd4 ~ years_since_1cd4 + age_group +  years_since_1cd4 * sexid + years_since_1cd4:exposureid + (years_since_1cd4 | patientindex) + (years_since_1cd4 | phylotype)"

## Empirical bayes for tuning R2D2 prior 
## Empirical Bayes prior: basing predicted R2 on lmer model 
empirical_bayes_r2d2_prior <- function(model) {
	r2p <- r2_nakagawa(model)
	# marginal R^2 quantifies the proportion of variance explained by fixed effects alone, 
	# while conditional R^2 captures the total variance explained by both fixed and random effects
	print(r2p)
	r2_marg <- r2p$R2_marginal
	bm1prior <- prior(R2D2(mean_R2 = r2_marg, prec_R2 = 2,cons_D2 = 0.5,autoscale = TRUE), class = b)
	bm1prior
}

ITER_ <- 5000 #3e3
WARMUP50_ <- ITER_/2.5
THIN_ <- 5
CORES_CHAINS_ <- 2
bmm_res_feff_r2d2_pts <- m0 <- eb_r2d2_prior <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		db <- residuals_removal_mx[[i,j]]$cd4_df_filt
		db <- db %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[i,j]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
		m0[[i,j]] <- lmer( lmmfo, data = db, REML=FALSE, control = lmerControl(optimizer = "bobyqa") )
		eb_r2d2_prior[[i,j]] <- empirical_bayes_r2d2_prior(m0[[i,j]])
		bmm_res_feff_r2d2_pts[[i,j]] <- fit_bayes_cd4_model_r2d2_prior(db, eb_r2d2_prior[[i,j]]) #priors_slope_interc = priors_bayes
	}
}

saveRDS(bmm_res_feff_r2d2_pts, glue("{RDS_PATH}/bmm_res_feff_r2d2_pts.rds"))
#write.csv( bmm_res_feff_r2d2_pts[[1,4]]$pt_effs[ order(bmm_res_feff_r2d2_pts[[1,4]]$pt_effs$Estimate), ], file=glue("{RESULTS_PATH}/tables/r2d2_B30.csv"), quote=F, row.names = F )

# join all pt effects together for mcs=30
cd4_fixeff_pts <- do.call(rbind, lapply(1:4, function(j) { 
	ss <- bmm_res_feff_r2d2_pts[[1, j]]$pt_effs ; ss$subtype <- tree_names[j]
	ss
}))

# order by estimate
cd4_fixeff_pts <- cd4_fixeff_pts[order(cd4_fixeff_pts$Estimate), ]
# order by posterior probability
cd4_fixeff_pts2 <- cd4_fixeff_pts[order(cd4_fixeff_pts$prob_negative, decreasing = T), ]

cd4_fixeff_pts2_table <- cd4_fixeff_pts2
cd4_fixeff_pts2_table <- cd4_fixeff_pts2_table %>% dplyr::select(subtype, phylotype, Q2.5, Estimate, Q97.5, Est.Error, p_value_bayes, prob_negative) #evidence_label
cd4_fixeff_pts2_table <- cd4_fixeff_pts2_table %>%	mutate(across(c(Q2.5, Estimate, Q97.5, Est.Error, prob_negative), ~ round(.x, 3)))
cd4_fixeff_pts2_table <- cd4_fixeff_pts2_table %>%	mutate(across(c(p_value_bayes), ~ signif(.x, 2)))
cd4_fixeff_pts2_table <- cd4_fixeff_pts2_table[order(-cd4_fixeff_pts2_table$prob_negative, cd4_fixeff_pts2_table$Estimate), ]
write.csv( cd4_fixeff_pts2_table, file=glue("{RESULTS_PATH}/tables/tableS10.csv"), quote=F, row.names = F )
saveRDS(cd4_fixeff_pts2_table, glue("{RDS_PATH}/cd4_fixeff_pts2_r2d2_table.rds"))

# For table S10 as well, extract residual sd, sds interc and slopes, and corr interc and slope
cd4_fixeff_pts2_table_sds <- readRDS(glue("{RDS_PATH}/bmm_res_feff_r2d2_pts.rds"))
cd4_fixeff_pts2_table_sds <- do.call(rbind, lapply(1:4, function(j) { 
	ss <- cd4_fixeff_pts2_table_sds[[1, j]]$sigma_df
	ss$param <- rownames(ss); rownames(ss) <- NULL; ss$subtype <- tree_names[j]; ss$`Est.Error` <- NULL; ss$param <- NULL; ss$group <- NULL
	ss2 <- cd4_fixeff_pts2_table_sds[[1, j]]$patient_sd_cor
	ss2$param <- rownames(ss2); rownames(ss2) <- NULL; ss2$subtype <- tree_names[j];	ss2$Rhat <- NULL; ss2$Bulk_ESS <- NULL; ss2$Tail_ESS <- NULL
	ss3 <- bind_rows(ss, ss2)
	ss3
}))
cd4_fixeff_pts2_tables_sds <- cd4_fixeff_pts2_table_sds %>% dplyr::select(subtype, param, `l-95% CI`, Estimate, `u-95% CI`)
cd4_fixeff_pts2_tables_sds <- cd4_fixeff_pts2_tables_sds %>%	mutate(across(c(`l-95% CI`, Estimate, `u-95% CI`), ~ round(.x, 3)))
write.table( cd4_fixeff_pts2_tables_sds, file=glue("{RESULTS_PATH}/tables/tableS10_2.tsv"), sep = "\t", quote=F, row.names = F )

# Read ML (0.1.1) & Bayes (0.1.2) results to see overlap in suspected VOIs
# ML
cd4_ml_rand_eff_table <- readRDS(glue("{RDS_PATH}/cd4_ml_rand_eff_table.rds"))
cd4_ml_rand_eff_table_p <- cd4_ml_rand_eff_table[cd4_ml_rand_eff_table$p_value <= 0.05 & cd4_ml_rand_eff_table$phylotype_coef < 0,] #12
# Bayes
cd4_ranef_pts2_table <- readRDS(glue("{RDS_PATH}/cd4_ranef_pts2_table.rds"))
cd4_ranef_pts2_table_p <- cd4_ranef_pts2_table[cd4_ranef_pts2_table$prob_negative > 0.8,] #20
# This model (R2D2)
cd4_fixeff_pts2_r2d2_table <- readRDS(glue("{RDS_PATH}/cd4_fixeff_pts2_r2d2_table.rds"))
cd4_fixeff_pts2_r2d2_table_p <- cd4_fixeff_pts2_r2d2_table[cd4_fixeff_pts2_r2d2_table$prob_negative > 0.8,] #30

# Overlap between ML & Bayes was n=10 as seen in 0.1.2 script
# Overlap ML & R2D2
overlap_ml_r2d2_cd4 <- inner_join(cd4_ml_rand_eff_table_p, cd4_fixeff_pts2_r2d2_table_p, by=c("subtype","phylotype")) #10
overlap_ml_r2d2_cd4[,c("subtype","phylotype")]
# subtype phylotype
# 1        B       133 *
# 2        B        69 *
# 3        B        90
# 4        B       118
# 5     A_A1         8
# 6     A_A1         3
# 7        B       137
# 8        B        84
# 9        B        24
# 10       B        62

# Overlap between Bayes & R2D2
overlap_bay_r2d2_cd4 <- inner_join(cd4_ranef_pts2_table_p, cd4_fixeff_pts2_r2d2_table_p, by=c("subtype","phylotype")) #20
overlap_bay_r2d2_cd4[,c("subtype","phylotype")]
# subtype phylotype
# 1     A_A1         3
# 2        B       133 *
# 3        B        83
# 4        B        90
# 5        B       118
# 6        B        28
# 7        B        24
# 8        B        54
# 9        B       137
# 10       B       151
# 11    A_A1         8
# 12       B        69 *
# 13       B       132
# 14       B        84
# 15       B       150
# 16       B        62
# 17       B         5
# 18       B        37
# 19       B        40
# 20       B       154

# Get union
ml_sub <- cd4_ml_rand_eff_table_p[, c("subtype", "phylotype")]
bay_sub <- cd4_ranef_pts2_table_p[, c("subtype", "phylotype")]
r2d2_sub <- cd4_fixeff_pts2_r2d2_table_p[, c("subtype", "phylotype")]
union_all <- unique(rbind(ml_sub, bay_sub, r2d2_sub)) #32
rownames(union_all) <- NULL
union_all[,c("subtype","phylotype")]
# subtype phylotype
# 1          B       133
# 2          B        69
# 3          B        90
# 4          B       118
# 5  CRF_02_AG         7
# 6       A_A1         8
# 7       A_A1         3
# 8          B       137
# 9          B        84
# 10         B        24
# 11         B        77
# 12         B        62
# 13         B        83
# 14         B        28
# 15         B        54
# 16         B       151
# 17         B       132
# 18         B       150
# 19         B         5
# 20         B        37
# 21         B        40
# 22         B       154
# 23         C        22
# 24         C        21
# 25         B        43
# 26         B        21
# 27      A_A1         4
# 28         B        65
# 29         B       124
# 30         B        95
# 31         B        55
# 32         B        67

# Intersection
inter12 <- inner_join(ml_sub, bay_sub, by = c("subtype", "phylotype")) #10
inter123 <- inner_join(inter12, r2d2_sub, by = c("subtype", "phylotype")) #10
inter123[,c("subtype","phylotype")]
intersect_cd4_susp_VOIs <- inter123
saveRDS(intersect_cd4_susp_VOIs, glue("{RDS_PATH}/intersect_cd4_susp_VOIs.rds"))
# subtype phylotype
# 1        B       133 *
# 2        B        69 *
# 3        B        90
# 4        B       118
# 5     A_A1         8
# 6     A_A1         3
# 7        B       137
# 8        B        84
# 9        B        24
# 10       B        62

# Extract mean and SD of negative slope coeffs from models
# ML
cd4_ml_rand_eff_table_p_neg <- cd4_ml_rand_eff_table_p[cd4_ml_rand_eff_table_p$phylotype_coef < 0,]
mean(cd4_ml_rand_eff_table_p_neg$phylotype_coef); sd(cd4_ml_rand_eff_table_p_neg$phylotype_coef) #-8.5, 2.8
hist(cd4_ml_rand_eff_table_p_neg$phylotype_coef)
# Bayes
mean(cd4_ranef_pts2_table_p$Estimate); sd(cd4_ranef_pts2_table_p$Estimate) # -7.8, 1.7
hist(cd4_ranef_pts2_table_p$Estimate)
# R2D2
mean(cd4_fixeff_pts2_r2d2_table_p$Estimate); sd(cd4_fixeff_pts2_r2d2_table_p$Estimate) # -12.7, 5.3
hist(cd4_fixeff_pts2_r2d2_table_p$Estimate)

# intercepts
cd4_fixeff_interc <- do.call(rbind, lapply(1:4, function(j) { 
	ss <- bmm_res_feff_r2d2_pts[[1, j]]$p_pt_intercepts ; ss$subtype <- tree_names[j]
	ss$phylotype <- gsub("^[^0-9]+", "", ss$phylotype)
	ss
}))
cd4_fixeff_interc_suspVOIs <- inner_join(cd4_fixeff_interc, intersect_cd4_susp_VOIs, by=c("subtype", "phylotype"))
median(cd4_fixeff_interc_suspVOIs$Estimate); mean(cd4_fixeff_interc_suspVOIs$Est.Error) # -2.8, 17.4
hist(cd4_fixeff_interc_suspVOIs$Estimate)

# Bayesian LMM (0.1.2) vs ML LMM (0.1.1): only mcs=30 and subtype B
db <- residuals_removal_mx[[1,4]]$cd4_df_filt
db <- db %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[i,j]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
if ( FALSE )
{
	library( brms )
	bm0 <- brm( lmmfo, data = db, iter = 1e4, chains = 4, cores = NCPU ) # 218 min 
	summary( bm0)
	
	## Comparison
	# The brms‐fitted hyperparameters for both groupings almost exactly match the lmer estimates:
	#
	#     Patient‐level
	#
	#         sd(Intercept): 216.23 (Brms) vs 216.15 (Lmer)
	#
	#         sd(years): 36.79 vs 36.73
	#
	#         cor: –0.48 vs –0.48
	#
	#     Phylotype‐level
	#
	#         sd(Intercept): 26.68 vs 25.65
	#
	#         sd(years): 8.46 vs 8.07
	#
	#         cor: 0.57 vs 0.72
	#
	# The fixed‐effect slopes and intercepts (e.g. overall intercept ~494, years_since_1cd4 ~–41) are also essentially identical. 
	
	## check ESS's. They look ok
	bm0pteffs <- ranef(bm0, prob = c(.025,0.975) )$phylotype
	bm0ptci <- bm0pteffs[, "years_since_1cd4", c("Q2.5", "Q97.5")]
	bm0pteffs1  <- bm0pteffs[, ,2 ]
}
