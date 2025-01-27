libs_load <- c("ggplot2","dplyr","ggpubr","data.table","glue","lubridate","sjPlot","viridis", "rstan", "brms","bayesplot",  "ape", "forcats","pbmcapply")
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

#priors_bayes <- c(mean_slope, sd_slope, mean_interc, sd_interc)
priors_bayes <- c(-20, 20, 0, 50)

rstan_options(auto_write = TRUE)
# After identyfing outliers using random effects model, fit individual fixed eff model for each VOI (consider the n=10 from last ML CD4 analysis), non-VOI & backbone to get effect size
fit_bayes_cd4_model_fixed_eff_vois_vs_backbone_vs_nv <- function(outcome, fostr, warmup_, cd4_df, change_prior=TRUE, priors_slope_interc) { #pt_id
	
	formula <- as.formula(paste(outcome, fostr))
	
	if(change_prior) {
		gp <- get_prior(formula, data=cd4_df, warmup = warmup_, iter = ITER, thin=THIN, cores = CORES_CHAINS, chains = CORES_CHAINS, seed = 1111)
		gp$prior[ grepl(gp$coef, pattern = 'years_since_1cd4:phylotype') ] <- glue("normal({priors_slope_interc[1]}, {priors_slope_interc[2]})")
		gp$prior[ grepl(gp$coef, pattern = '^phylotype') ] <- glue("normal({priors_slope_interc[3]}, {priors_slope_interc[4]})")
		#print(gp)
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
	p_all_rows$phylotype_coef <- rownames(p_all_rows)
	rownames(p_all_rows) <- NULL 
	p_all_rows <- p_all_rows %>% dplyr::select(phylotype_coef, Estimate, Est.Error, Q2.5, Q97.5)
	p_pt_only <- p_all_rows
	p_pt_only <- p_pt_only[ grepl(p_pt_only$phylotype_coef, pattern = 'years_since_1cd4:phylotype'), ]
	p_pt_only$phylotype_coef <- gsub("years_since_1cd4:phylotype", "", as.character(p_pt_only$phylotype_coef))
	
	p_pt_intercepts <- p_all_rows[grepl(p_all_rows$phylotype_coef, pattern = '^phylotype'),]
	
	return(list(all_covar=p_all_rows, pt_effs=p_pt_only, p_pt_intercepts=p_pt_intercepts, loo_model=loo_model)) # df_p=df_p
}

ITER <- 10000
WARMUP10 <- ITER/10
WARMUP50 <- ITER/2
THIN <- 10
CORES_CHAINS <- 2

# vs backbone
residuals_removal_mx2 <- residuals_removal_mx 
rownames(residuals_removal_mx2) <- c("30","50","100")
colnames(residuals_removal_mx2) <- c("A_A1","CRF_02_AG","C", "B")
backbone_cl_control2 <- backbone_cl_control
rownames(backbone_cl_control2) <- c("30","50","100")
colnames(backbone_cl_control2) <- c("A_A1","CRF_02_AG","C", "B")

cd4_model_form_no_intercept_pt <- " ~ years_since_1cd4 + age_group + years_since_1cd4 * sexid + years_since_1cd4:exposureid + years_since_1cd4:phylotype + (years_since_1cd4 | patientindex)"
cd4_model_form_with_intercept_pt <- " ~ years_since_1cd4 + age_group + years_since_1cd4 * sexid + years_since_1cd4:exposureid + years_since_1cd4:phylotype + phylotype + (years_since_1cd4 | patientindex)"
d_bmm_fixedeff_vois <- d_bmm_fixedeff_vois_ref_bb <- bmm_res_fixeff_vois_vs_bb <- list()
for(i in 1:nrow(cd4_ml_fixeff_pot_vois)) {
	print(glue("n={i}"))
	print(paste0( cd4_ml_fixeff_pot_vois[i,"mcs"],"-", cd4_ml_fixeff_pot_vois[i,"subtype"], "-", cd4_ml_fixeff_pot_vois[i,"phylotype"] ))
	d_bmm_fixedeff <- residuals_removal_mx2[[ cd4_ml_fixeff_pot_vois[i, "mcs"], cd4_ml_fixeff_pot_vois[i, "subtype"] ]]$cd4_df_filt
	
	print(nrow(d_bmm_fixedeff))
	d_bmm_fixedeff_vois[[i]] <- d_bmm_fixedeff
	d_bmm_fixedeff_vois[[i]]$phylotype <- ifelse( d_bmm_fixedeff$cluster ==  as.integer(backbone_cl_control2[[ cd4_ml_fixeff_pot_vois[i, "mcs"], cd4_ml_fixeff_pot_vois[i, "subtype"] ]]), 
																																											yes=as.character(backbone_cl_control2[[ cd4_ml_fixeff_pot_vois[i, "mcs"], cd4_ml_fixeff_pot_vois[i, "subtype"] ]]),
																																											no=ifelse( d_bmm_fixedeff$cluster == as.integer(cd4_ml_fixeff_pot_vois$phylotype[ i ]), 
																																																						yes= as.integer(cd4_ml_fixeff_pot_vois$phylotype[ i ]), 
																																																						no=as.character("non-VOI") ))
	d_bmm_fixedeff_vois[[i]] <- d_bmm_fixedeff_vois[[i]][ d_bmm_fixedeff_vois[[i]]$phylotype != "non-VOI", ]
	d_bmm_fixedeff_vois[[i]]$phylotype <- as.factor(d_bmm_fixedeff_vois[[i]]$phylotype)
	print(table(d_bmm_fixedeff_vois[[i]]$phylotype))
	d_bmm_fixedeff_vois_ref_bb[[i]] <- d_bmm_fixedeff_vois[[i]] %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control2[[ cd4_ml_fixeff_pot_vois[i, "mcs"], cd4_ml_fixeff_pot_vois[i, "subtype"] ]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
	#View(d_bmm_fixedeff_vois_ref_bb[[i]])
	# IMPORTANT: NOW changing prior
	bmm_res_fixeff_vois_vs_bb[[i]] <- fit_bayes_cd4_model_fixed_eff_vois_vs_backbone_vs_nv("cd4",cd4_model_form_with_intercept_pt,WARMUP50, d_bmm_fixedeff_vois_ref_bb[[i]], change_prior=TRUE, priors_slope_interc = priors_bayes) #as.integer(sigs[i]), glue("{min_cl_size_choices[imcs]}-{tree_names[jtree]}-cd4-bayes-fixed_eff_voi_{sigs[i]}_vs_bb")
	# bmm_res_fixeff_vois_vs_bb[[i]] <- fit_bayes_cd4_model_fixed_eff_vois_vs_backbone_vs_nv("cd4",
	# 																																																																																							cd4_model_form_no_intercept_pt,
	# 																																																																																							WARMUP50, d_bmm_fixedeff_vois_ref_bb[[i]], change_prior=TRUE, priors_slope_interc = cd4_ml_fixeff_pot_vois$prior_bayes_se[i])
	bmm_res_fixeff_vois_vs_bb[[i]]$pt_effs$p_value <- 2 * (1 - pnorm(abs(bmm_res_fixeff_vois_vs_bb[[i]]$pt_effs$Estimate) / bmm_res_fixeff_vois_vs_bb[[i]]$pt_effs$Est.Error))
	print(bmm_res_fixeff_vois_vs_bb[[i]]$pt_effs)
	bmm_res_fixeff_vois_vs_bb[[i]]$pt_effs$mcs <- cd4_ml_fixeff_pot_vois$mcs[i]
	bmm_res_fixeff_vois_vs_bb[[i]]$pt_effs$subtype <- cd4_ml_fixeff_pot_vois$subtype[i]
	gc()
}

saveRDS(bmm_res_fixeff_vois_vs_bb, glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb.rds"))

# join all pt effects together
cd4_bayes_pt_effs_vs_bb <- do.call(rbind, lapply(bmm_res_fixeff_vois_vs_bb, function(x) x$pt_effs))
# join intercepts
cd4_bayes_interc_effs_vs_bb <- do.call(rbind, lapply(bmm_res_fixeff_vois_vs_bb, function(x) x$p_pt_intercepts))
# extract signif
cd4_bayes_pt_effs_vs_bb_signif <- cd4_bayes_pt_effs_vs_bb[cd4_bayes_pt_effs_vs_bb$Q2.5 < 0 & cd4_bayes_pt_effs_vs_bb$Q97.5 < 0 & cd4_bayes_pt_effs_vs_bb$phylotype_coef != "nonMVOI",]
# 11/12 signif

# Table S9: Bayesian CD4 model with fixed effect on phylotype (vs backbone)
cd4_bayes_fixed_eff_tables <- readRDS(glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb.rds"))
cd4_bayes_fixed_eff_table <- do.call(rbind, lapply(cd4_bayes_fixed_eff_tables, function(x) x$pt_effs))
#cd4_bayes_fixed_eff_table <- cd4_bayes_fixed_eff_table[cd4_bayes_fixed_eff_table$phylotype_coef != "nonMVOI",] #23
cd4_bayes_fixed_eff_table <- cd4_bayes_fixed_eff_table[cd4_bayes_fixed_eff_table$mcs == 30,] #12
table(cd4_bayes_fixed_eff_table$subtype)
# A_A1      B   CRF_02_AG 
# 2         9         1
cd4_bayes_fixed_eff_table <- cd4_bayes_fixed_eff_table %>% dplyr::select( subtype, phylotype_coef, Q2.5, Estimate, Q97.5, Est.Error, p_value )
cd4_bayes_fixed_eff_table <- cd4_bayes_fixed_eff_table %>% arrange( p_value, Estimate )
cd4_bayes_fixed_eff_table$Q2.5 <- round(cd4_bayes_fixed_eff_table$Q2.5, 3)
cd4_bayes_fixed_eff_table$Estimate <- round(cd4_bayes_fixed_eff_table$Estimate, 3)
cd4_bayes_fixed_eff_table$Q97.5 <- round(cd4_bayes_fixed_eff_table$Q97.5, 3)
cd4_bayes_fixed_eff_table$Est.Error <- round(cd4_bayes_fixed_eff_table$Est.Error, 3)
cd4_bayes_fixed_eff_table$p_value <- signif(cd4_bayes_fixed_eff_table$p_value, digits = 2)
options(scipen=999)
write.csv( cd4_bayes_fixed_eff_table, file=glue("{RESULTS_PATH}/tables/tableS9.csv"), quote=F, row.names = F )

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
#colnames(vl_t_test_adj_fdr) <- c("subtype", "phylotype_id","vl_ttest_adj_fdr_estim","p_value1") 
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
vl_bayes_joint$phylotype_id <- vl_bayes_joint$group #as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", vl_bayes_joint$param))
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

# join all pt effects together
cd4_bayes_fixed_eff <- readRDS(glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb.rds"))
cd4_bayes_fixed_eff <- do.call(rbind, lapply(cd4_bayes_fixed_eff, function(x) x$pt_effs))
cd4_bayes_fixed_eff <- cd4_bayes_fixed_eff[cd4_bayes_fixed_eff$mcs == 30,] # mcs=30
cd4_bayes_fixed_eff <- cd4_bayes_fixed_eff[cd4_bayes_fixed_eff$phylotype_coef !=  "nonMVOI",] #12
rownames(cd4_bayes_fixed_eff) <- NULL
cd4_bayes_fixed_eff$cd4_bayes_fixed_eff_estim <- cd4_bayes_fixed_eff$Estimate
cd4_bayes_fixed_eff$p_value4 <- signif(cd4_bayes_fixed_eff$p_value, digits = 2)
cd4_bayes_fixed_eff <- cd4_bayes_fixed_eff %>% dplyr::select(subtype, phylotype_coef, cd4_bayes_fixed_eff_estim, p_value4) %>% arrange(p_value4) #cd4_bayes_fixed_eff
colnames(cd4_bayes_fixed_eff) <- c("subtype", "phylotype_id","cd4_bayes_fixed_eff_estim","p_value4") #12 rows
cd4_bayes_fixed_eff_sig <- cd4_bayes_fixed_eff[cd4_bayes_fixed_eff$p_value <= 0.05,] # 9 rows
colnames(cd4_bayes_fixed_eff_sig) <- c("subtype", "phylotype_id","cd4_bayes_fixed_eff_estim","p_value4") #12 rows

# cd4 (sample size): careful if doing for other subtypes as well, because this one is only for B30
cd4_ss <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))[1,] # mcs=30
cd4_ss <- lapply(cd4_ss, `[[`, 3) #$cd4_df_filt
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
																																																																		(phylotype_id %in% cd4_bayes_fixed_eff_sig$phylotype_id & subtype %in% cd4_bayes_fixed_eff_sig$subtype)) #23 rows
write.csv(combined_vl_cd4_estims_f2, glue("{RESULTS_PATH}/tables/table2.csv"), quote = F, row.names = F, na = "") # 30 rows
#combined_vl_cd4_estims <- merge(, by = "phylotype_id", all = TRUE)

# Shortlist
combined_vl_cd4_estims <- as.data.frame(combined_vl_cd4_estims)
i0 = with( combined_vl_cd4_estims, (is.na(p_value1) | p_value1 <.05) | (is.na(p_value2) | p_value2 <.05) | (is.na(p_value3) | p_value3 <.05) | (is.na(p_value4) | p_value4 <.05))
pvnames <- c('p_value1', 'p_value2', 'p_value3', 'p_value4')
i0 = apply( combined_vl_cd4_estims[,pvnames]|>as.matrix(), MAR=1, FUN=function(x) {
	if( all( is.na(x))) return(FALSE)
	( min( na.omit(x)) < .05)
})
d0 <- combined_vl_cd4_estims[ i0, ] #nrow(d0): 35
d0

i1 <- combined_vl_cd4_estims$vl_bayes_joint_estim > median(combined_vl_cd4_estims$vl_bayes_joint_estim )
d1 <- combined_vl_cd4_estims[ i0 & i1, ]  #nrow(d1): 20

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
vl_ttest_sl$phylotype_id # VOI 2: PT40 (smallest p and higher vl estimate)  40  20  79  69 125 137   5  14  18

# order vl_bayes by increasing p
vl_bayes_sl <- combined_vl_cd4_estims %>% filter( vl_bayes_joint_estim > median(vl_bayes_joint_estim ) & p_value2<.05 ) %>% arrange(p_value2)
vl_bayes_sl$phylotype_id # VOI 2: 40 (smallest p and higher vl estimate), 101, 20, 69

# order cd4_bayes_fixed_eff by increasing p
cd4_bayes_fe_sl <- combined_vl_cd4_estims %>% filter( p_value4<.05 ) %>% arrange(p_value4)
cd4_bayes_fe_sl$phylotype_id # 133 (smallest p and fastest estim/decline among B) 90 69 118
cd4_bayes_fe_sl %>% select(subtype, phylotype_id, cd4_bayes_fixed_eff_estim)
# subtype phylotype_id cd4_bayes_fixed_eff_estim
# 1       A_A1            3                   -30.802
# 2          B          133                   -39.180 (VOI 3: PT133)
# 3       A_A1            8                   -30.278
# 4          B           90                   -32.215
# 5          B          118                   -25.543
# 6          B           69                   -28.386
# 7          B          137                   -25.031
# 8          B           84                   -31.113
# 9          B           24                   -24.093
# 10         B           62                   -22.940
# 11 CRF_02_AG            7                   -29.878

# order cd4_bayes_rand_eff by increasing p
cd4_bayes_re_sl <- combined_vl_cd4_estims %>% filter( cd4_ml_random_eff_estim < 0 & p_value3<.05 ) %>% arrange(p_value3)
cd4_bayes_re_sl$phylotype_id #133 (smallest p and fastest estim/decline among B)  69  90 118 137  84  24  77  24 62
cd4_bayes_re_sl %>% select(subtype, phylotype_id, cd4_ml_random_eff_estim, p_value3)
# subtype phylotype_id cd4_ml_random_eff_estim p_value3
# 1          B          133                 -10.182 0.000059 (VOI 3: PT133)
# 2  CRF_02_AG            7                  -9.919 0.000860
# 3          B           69                  -7.978 0.002600
# 4          B           90                  -7.987 0.003100
# 5          B          118                  -8.773 0.003700
# 6          B          137                  -8.610 0.013000
# 7       A_A1            8                  -7.336 0.018000
# 8       A_A1            3                 -15.564 0.022000
# 9          B           84                  -8.173 0.026000
# 10         B           77                  -3.687 0.027000
# 11         B           24                  -7.638 0.033000
# 12         B           62                  -5.907 0.042000

# Table S11: coeffs of covariates for VOIs in fixed effects model
cd4_bayes_fixed_eff_tables <- readRDS(glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb.rds"))
cd4_ml_fixeff_pot_vois <- readRDS(glue("{RDS_PATH}/cd4_ml_fixeff_pot_vois.rds"))
#cd4_bayes_fixed_eff_table2 <- do.call(rbind, lapply(cd4_bayes_fixed_eff_tables, function(x) data.frame(x$all_covar, mcs=x$pt_effs$mcs, subtype=x$pt_effs$subtype))) #322
cd4_bayes_fixed_eff_table2 <- do.call(rbind, lapply(seq_along(cd4_bayes_fixed_eff_tables), function(i) {
	x <- cd4_bayes_fixed_eff_tables[[i]]
	data.frame(
		x$all_covar, 
		mcs = x$pt_effs$mcs, 
		subtype = x$pt_effs$subtype, 
		phylotype_id = cd4_ml_fixeff_pot_vois[i,"phylotype"]
	)
}))

cd4_bayes_fixed_eff_table2 <- cd4_bayes_fixed_eff_table2[ cd4_bayes_fixed_eff_table2$mcs == 30 & cd4_bayes_fixed_eff_table2$subtype == "B", ] #126
cd4_bayes_fixed_eff_table2.1 <- cd4_bayes_fixed_eff_table2[!(cd4_bayes_fixed_eff_table2$phylotype_coef %in%  c("Intercept","years_since_1cd4","years_since_1cd4:exposureidOther")),] #99, s"years_since_1cd4:phylotypenonMVOI"
# remove PT77 because p-values not signif
cd4_bayes_fixed_eff_table2.2 <- cd4_bayes_fixed_eff_table2.1[!(cd4_bayes_fixed_eff_table2.1$phylotype_id %in% c("77")),]
cd4_bayes_fixed_eff_table2.2 <- cd4_bayes_fixed_eff_table2.2 %>% arrange(Estimate)
cd4_bayes_fixed_eff_table2.2 <- cd4_bayes_fixed_eff_table2.2 %>% dplyr::select( subtype, phylotype_id, phylotype_coef, Q2.5, Estimate, Q97.5, Est.Error )
cd4_bayes_fixed_eff_table2.2 <- cd4_bayes_fixed_eff_table2.2 %>% mutate(across(c(Q2.5, Estimate, Q97.5, Est.Error), ~ round(.x, 3)))
write.csv( cd4_bayes_fixed_eff_table2.2, file=glue("{RESULTS_PATH}/tables/tableS11.csv"), quote=F, row.names = F )
#View(cd4_bayes_fixed_eff_table2.2)
# median effects by age_group60P (intercept)
cd4_bayes_fixed_eff_table2.2_60_older <- cd4_bayes_fixed_eff_table2.2 %>% filter(phylotype_coef == "age_group60P") %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))
# median effects by age_group50M59 (intercept)
cd4_bayes_fixed_eff_table2.2_50s <- cd4_bayes_fixed_eff_table2.2 %>% filter(phylotype_coef == "age_group50M59") %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))
# median effects by females (intercept)
cd4_bayes_fixed_eff_table2.2_females <- cd4_bayes_fixed_eff_table2.2 %>% filter(phylotype_coef == "sexidFemale") %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))
# median effects by age_group40M49 (intercept)
cd4_bayes_fixed_eff_table2.2_40s <- cd4_bayes_fixed_eff_table2.2 %>% filter(phylotype_coef == "age_group40M49") %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))
# median effects of phylotype (intercept)
cd4_bayes_fixed_eff_table2.2_pt_interc <- cd4_bayes_fixed_eff_table2.2 %>%filter(grepl("^phylotype", phylotype_coef)) %>% 
	summarise(Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5 = median(Q97.5))

# Fig. 2B: Plot average trajectory of VOIs
# IMPORTANT: even though PT40 is not a "CD4" VOI, run Bayes fixed effects model on it to include in Fig. and provide comparison against the "2 CD4 VOIs"
residuals_removal_mx <- readRDS( glue("{RDS_PATH}/residuals_removal_mx.rds") )
d_pt40_run <- residuals_removal_mx[[ 1,4 ]]$cd4_df_filt
d_pt40_run$phylotype <- ifelse( d_pt40_run$cluster ==  backbone_cl_control[[ 1,4 ]], yes=backbone_cl_control[[ 1,4 ]],no=
																																	ifelse( d_pt40_run$cluster == "40", yes= "40", no=as.character("non-VOI") ))
d_pt40_run <- d_pt40_run[d_pt40_run$phylotype != "non-VOI",]
d_pt40_run$phylotype <- as.factor(d_pt40_run$phylotype)
print(table(d_pt40_run$phylotype))
d_pt40_run_ref_bb <- d_pt40_run %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[ 1,4 ]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
# first fit ML model and get SE to use as prior in bayes (NOT CHANCING PRIOR ANYMORE)
ml_pt40 <- fit_cd4_model_fixed_eff_vois_vs_backbone_vs_nv("cd4", cd4_model_form_with_intercept_pt, d_pt40_run_ref_bb)
#ml_pt40_prior_bayes_se <- 2 * ml_pt40$pt_effs$stderr[ml_pt40$pt_effs$phylotype == "40"]
# run bayes
bmm_res_pt40_run_ref_bb <- fit_bayes_cd4_model_fixed_eff_vois_vs_backbone_vs_nv("cd4", cd4_model_form_with_intercept_pt, WARMUP50, d_pt40_run_ref_bb, change_prior=TRUE, priors_slope_interc = priors_bayes)
saveRDS(bmm_res_pt40_run_ref_bb, glue("{RDS_PATH}/bmm_res_pt40_run_ref_bb.rds"))
bmm_res_pt40_run_ref_bb$pt_effs$p_value <- 2 * (1 - pnorm(abs(bmm_res_pt40_run_ref_bb$pt_effs$Estimate) / bmm_res_pt40_run_ref_bb$pt_effs$Est.Error))
bmm_res_pt40_run_ref_bb$pt_effs$mcs <- "30"; bmm_res_pt40_run_ref_bb$pt_effs$subtype <- "B"; bmm_res_pt40_run_ref_bb$pt_effs$phylotype_id <- "40"
#bmm_res_pt40_run_ref_bb$pt_effs <- bmm_res_pt40_run_ref_bb$pt_effs[ bmm_res_pt40_run_ref_bb$pt_effs$phylotype_coef != "nonMVOI", ]
pt40_other_coeffs <- bmm_res_pt40_run_ref_bb$all_covar; pt40_other_coeffs$mcs <- "30"; pt40_other_coeffs$subtype <- "B"; pt40_other_coeffs$phylotype_id <- "40"

bmm_vois_combined <- readRDS(glue("{RDS_PATH}/bmm_res_fixeff_vois_vs_bb.rds"))
bmm_vois_combined2 <- do.call(rbind, lapply(seq_along(bmm_vois_combined), function(i) {
	x <- bmm_vois_combined[[i]]
	data.frame(x$pt_effs, phylotype_id = cd4_ml_fixeff_pot_vois[i,"phylotype"])
}))
bmm_vois_combined2 <- bmm_vois_combined2[bmm_vois_combined2$phylotype_coef != "nonMVOI",]
bmm_vois_combined2 <- bmm_vois_combined2[(bmm_vois_combined2$phylotype_id %in% subtype_b_vois_ids_only) & bmm_vois_combined2$mcs == "30" & bmm_vois_combined2$subtype=="B", ]
bmm_vois_combined3 <- rbind(bmm_vois_combined2, bmm_res_pt40_run_ref_bb$pt_effs)

# Standardise initial (baseline) CD4 at 500 and use Bayesian regression coeffs to represent decline of CD4 counts over time
bmm_vois_combined3$phylotype <- paste0("PT.B.",bmm_vois_combined3$phylotype_id,".UK")

cd4_bayes_fixed_eff_table2.a <- rbind(cd4_bayes_fixed_eff_table2, pt40_other_coeffs)
cd4_bayes_fixed_eff_table2.b <- cd4_bayes_fixed_eff_table2.a[ (cd4_bayes_fixed_eff_table2.a$phylotype_id %in% subtype_b_vois_ids_only) & cd4_bayes_fixed_eff_table2.a$mcs == "30" & cd4_bayes_fixed_eff_table2.a$subtype == "B", ]
cd4_bayes_fixed_eff_table2.b$phylotype <- paste0("PT.B.",cd4_bayes_fixed_eff_table2.b$phylotype_id,".UK")
# usual rate of decline
rate_decl <- cd4_bayes_fixed_eff_table2.b[ grepl(cd4_bayes_fixed_eff_table2.b$phylotype_coef, pattern = '^years_since_1cd4$'), ]#$Estimate
rate_decl <- rate_decl %>% summarise(phylotype="Backbone", Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5=median(Q97.5))

# add usual decline per year together with phylotype effects
bmm_vois_combined4 <- bmm_vois_combined3
bmm_vois_combined4$Estimate <- rate_decl$Estimate + bmm_vois_combined4$Estimate
bmm_vois_combined4$Q2.5 <- rate_decl$Q2.5 + bmm_vois_combined4$Q2.5
bmm_vois_combined4$Q97.5 <- rate_decl$Estimate + bmm_vois_combined4$Q97.5
bmm_vois_combined4_bb <- bind_rows(bmm_vois_combined4, rate_decl)

# usual intercepts
intercs <- cd4_bayes_fixed_eff_table2.b[ grepl(cd4_bayes_fixed_eff_table2.b$phylotype_coef, pattern = '(Intercept)'), ]
intercs <- intercs %>% summarise(phylotype="Backbone", Estimate=median(Estimate), Q2.5=median(Q2.5), Q97.5=median(Q97.5))
intercs$Estimate_interc <- intercs$Estimate
intercs$Q2.5_interc <- intercs$Q2.5
intercs$Q97.5_interc <- intercs$Q97.5

# VOI intercepts
bmm_vois_combined_interc <- do.call(rbind, lapply(seq_along(bmm_vois_combined), function(i) {
	x <- bmm_vois_combined[[i]]
	data.frame(x$p_pt_intercepts, phylotype_id = cd4_ml_fixeff_pot_vois[i,"phylotype"])
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

#bmm_pred_summary <- inner_join(comb_cd4_vois, comb_cd4_vois_non, by=c("phylotype"="phylotype.x"))
f2b <- ggplot(new_data_coeffs_cd4s, aes(x = years_since_1cd4, y = cd4_count)) + #fill = phylotype, color = phylotype
	#geom_ribbon(aes(ymin = cd4_lower, ymax = cd4_upper, fill = phylotype, color = phylotype), alpha = 0.15) + #color=NA
	geom_line(aes(color = phylotype), size = 0.75) +
	scale_color_manual(values=cd4_regr_pt_pal, name="Phylotype") + scale_fill_manual(values=cd4_regr_pt_pal, name="Phylotype") +
	ylab(bquote(bold("Pre-treatment CD4 count (cells /"~mm^3~")"))) + xlab("Years since first CD4") + theme_classic() + coord_cartesian(xlim=c(0,6), ylim = c(0,500)) +
	scale_x_continuous(breaks=seq(from=0, to=6, by=1)) + scale_y_continuous(limits=c(0,500), breaks=seq(from=0,to=500,by=50)) + coord_cartesian(expand = c(0, 0)) +
	theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'),
							legend.key.width = unit(0.5, 'cm'), legend.text=element_text(size=10), legend.position = c(0.25,0.3)) + leg +
	geom_hline(yintercept=350, linetype="longdash", color="black", linewidth=0.3)
saveRDS(f2b, glue("{RDS_PATH}/f2b.rds"))

# time to reach 350 cells/mm3
new_data_coeffs_cd4s_time_to_350 <- new_data_coeffs_cd4s %>%
	filter(cd4_count <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350$years_since_1cd4

new_data_coeffs_cd4s_time_to_350_lower <- new_data_coeffs_cd4s %>%
	filter(cd4_lower <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_lower$years_since_1cd4

new_data_coeffs_cd4s_time_to_350_upper <- new_data_coeffs_cd4s %>%
	filter(cd4_upper <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_upper$years_since_1cd4

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
	ylab(bquote(bold("Pre-treatment CD4 count (cells /"~mm^3~")"))) + xlab("Years since first CD4") + theme_classic() + coord_cartesian(xlim=c(0,6), ylim = c(0,500)) +
	scale_x_continuous(breaks=seq(from=0, to=6, by=1)) + scale_y_continuous(limits=c(0,500), breaks=seq(from=0,to=500,by=50)) + coord_cartesian(expand = c(0, 0)) +
	theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'),
							legend.key.width = unit(0.5, 'cm'), legend.text=element_text(size=10), legend.position = c(0.25,0.3)) + leg +
	geom_hline(yintercept=350, linetype="longdash", color="black", linewidth=0.3)
ggsave(s_ci, file=glue("{RESULTS_PATH}/figs/figS10.eps"), device=cairo_ps, dpi=600, width=8, height=6, bg="white")
ggsave(s_ci, file=glue("{RESULTS_PATH}/figs/figS10.jpg"), dpi=600, width=8, height=6, bg="white")

new_data_coeffs_cd4s_time_to_350_comm <- new_data_coeffs_cd4s_common_intercept %>%
	filter(cd4_count <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_comm$years_since_1cd4

new_data_coeffs_cd4s_time_to_350_lower_comm <- new_data_coeffs_cd4s_common_intercept %>%
	filter(cd4_lower <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_lower_comm$years_since_1cd4

new_data_coeffs_cd4s_time_to_350_upper_comm <- new_data_coeffs_cd4s_common_intercept %>%
	filter(cd4_upper <= 350) %>% group_by(phylotype) %>% slice_min(order_by = years_since_1cd4, n = 1, with_ties = FALSE) %>% ungroup()
new_data_coeffs_cd4s_time_to_350_upper_comm$years_since_1cd4

# difference in estimates when comparing estimated baseline CD4 and common baseline CD4
new_data_coeffs_cd4s_time_to_350$years_since_1cd4 - new_data_coeffs_cd4s_time_to_350_comm$years_since_1cd4
new_data_coeffs_cd4s_time_to_350_lower$years_since_1cd4 - new_data_coeffs_cd4s_time_to_350_lower_comm$years_since_1cd4
new_data_coeffs_cd4s_time_to_350_upper$years_since_1cd4 - new_data_coeffs_cd4s_time_to_350_upper_comm$years_since_1cd4