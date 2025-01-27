libs_load <- c( "ggplot2", "dplyr", "glue", "ggtree", "reshape2", "lubridate", "data.table", "ggpubr", "scales")
invisible( lapply(libs_load, library, character.only=TRUE) )

# similar CD4 decline model as before
RDS_PATH <- "rds"

cd4_md <- readRDS(glue("{RDS_PATH}/cd4_md.rds"))

# patients with 2 meas regardless of subtype and treatment status
cd4_md_mult_meas <- cd4_md[cd4_md$patientindex %in% unique(cd4_md$patientindex[duplicated(cd4_md$patientindex)]),]

keep_cd4_after_art <- function(demog_md_choice, subtype_choice) {
	demog_md_choice <- demog_md_choice[demog_md_choice$rega3subtype == subtype_choice,]
	print("rows subtype md:")
	print(nrow(demog_md_choice))
	demog_md_choice$artstart_decimal_date <- decimal_date( as.Date(gsub("\\/", "15", demog_md_choice$artstart_my), "%m%d%Y") )
	demog_md_cd4_merged <- demog_md_choice %>% left_join(cd4_md_mult_meas, by="patientindex", relationship="many-to-many")
	suppressWarnings( cd4_after_art <- demog_md_cd4_merged %>% filter(cd4_decimal_date > artstart_decimal_date) )
	
	# for the df above I have for each sample all measurements of CD4
	lookup_pat_tests <- unique( cd4_after_art[ , c('patientindex','testindex') ] )
	lookup_pat_tests$testindex <- paste0("t.",lookup_pat_tests$testindex)
	
	# Here removes duplicates of patient and cd4 date
	cd4_after_art_cl <- cd4_after_art %>% distinct(patientindex, cd4_decimal_date, .keep_all = TRUE)
	cd4_after_art_cl <- cd4_after_art_cl[ , !names(cd4_after_art_cl) %in%  c("testindex","dbsample_date")]
	cd4_after_art_cl <- cd4_after_art_cl[cd4_after_art_cl$patientindex %in% unique(cd4_after_art_cl$patientindex[duplicated(cd4_after_art_cl$patientindex)]),]
	#View(cd4_after_art_cl)
	print("nrows cd4 after art")
	print(nrow(cd4_after_art_cl))
	print("unique patients after art")
	print(length(unique(cd4_after_art_cl$patientindex)))
	
	print("unique patients lookup_pat_tests")
	lookup_pat_tests <- lookup_pat_tests[lookup_pat_tests$patientindex %in% cd4_after_art_cl$patientindex,]
	print(length(unique(lookup_pat_tests$patientindex)))
	
	list(cd4_after_art_cl, lookup_pat_tests)
}

demog_md_subtype_match <- readRDS(glue("{RDS_PATH}/demog_md_subtype_match.rds")) #100591
demog_md_subtype_match <- demog_md_subtype_match[ (demog_md_subtype_match$status=="Naïve") & 
																																																			(demog_md_subtype_match$rega3subtype %in% c("A (A1)","CRF 02_AG","C","B")) & 
																																																			(demog_md_subtype_match$exposureid != "Not known"), ] #53382

subtype_choices <- c("A (A1)","CRF 02_AG","C", "B")
cd4_after_art_subtypes <- list()
for(i in 1:length(subtype_choices)) {
	print(subtype_choices[i])
	cd4_after_art_subtypes[[i]] <- keep_cd4_after_art(demog_md_subtype_match, subtype_choices[i])
}

# subtype B: 2544 rows, 2342 patients
treestruct_min_cl_size_res_yes_sup <- readRDS(glue("{RDS_PATH}/treestruct_min_cl_size_res_yes_sup.rds")) 
backbone_cl_control <- readRDS(glue("{RDS_PATH}/backbone_cl_control.rds")) 
cd4_after_art_subtypes_B <- cd4_after_art_subtypes[[4]]
cd4_after_art_subtypes_B_pts1 <- inner_join(treestruct_min_cl_size_res_yes_sup[[1,4]][[2]], cd4_after_art_subtypes_B[[2]], by=c("taxon"="testindex"))
cd4_after_art_subtypes_B_pts2 <- inner_join(cd4_after_art_subtypes_B_pts1, cd4_after_art_subtypes_B[[1]], by="patientindex") #3802 rows, 1738 patients
table(cd4_after_art_subtypes_B_pts2$cluster) # 32 rows PT40, 2 rows PT69, 6 rows PT133, (9 rows PT137), 2408 backbone
length(unique( cd4_after_art_subtypes_B_pts2[cd4_after_art_subtypes_B_pts2$cluster == 40,]$patientindex )) #12
length(unique( cd4_after_art_subtypes_B_pts2[cd4_after_art_subtypes_B_pts2$cluster == 69,]$patientindex )) #1
length(unique( cd4_after_art_subtypes_B_pts2[cd4_after_art_subtypes_B_pts2$cluster == 133,]$patientindex )) #3
length(unique( cd4_after_art_subtypes_B_pts2[cd4_after_art_subtypes_B_pts2$cluster == 153,]$patientindex )) #1105
cd4_after_art_subtypes_B_pts2$phylotype <- cd4_after_art_subtypes_B_pts2$cluster
# adjust variables
cd4_after_art_subtypes_B_pts2 <- cd4_after_art_subtypes_B_pts2 %>% group_by(patientindex) %>% arrange(cd4_decimal_date, .by_group = TRUE) %>% mutate(measurement_id = row_number()) %>% ungroup()
# IMPORTANT: years_since_1cd4 is actually years since first CD4 post-treatment
cd4_after_art_subtypes_B_pts2 <- cd4_after_art_subtypes_B_pts2 %>% group_by(patientindex) %>% mutate(years_since_1cd4 = round(cd4_decimal_date - cd4_decimal_date[measurement_id==1], 3) )
cd4_after_art_subtypes_B_pts2$hiv_diag_decimal_date <- decimal_date(as.Date(cd4_after_art_subtypes_B_pts2$hivpos_ymd))
cd4_after_art_subtypes_B_pts2 <- cd4_after_art_subtypes_B_pts2 %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y))
cd4_after_art_subtypes_B_pts2 <- cd4_after_art_subtypes_B_pts2 %>% mutate(
	age_group = dplyr::case_when(age_diag<=29 ~ "0-29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
	age_group = factor(age_group,level = c("0-29","30-39","40-49","50-59","60+")))
table(cd4_after_art_subtypes_B_pts2$age_group)
# 0-29 30-39 40-49 50-59   60+ 
# 838  1382   931   413   177 

# model 1: ML random effect model of post-treatment data to find potential outlier phylotypes (smaller increase OR even decrease)
cd4_after_art_subtypes_B_pts2 <- cd4_after_art_subtypes_B_pts2 %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[1,4]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
cd4_model_form <- paste(" ~ years_since_1cd4 + age_group +  years_since_1cd4 * sexid + years_since_1cd4:exposureid + (years_since_1cd4 | patientindex) + (years_since_1cd4 | phylotype)")
# IMPORTANT: load function "fit_cd4_model" from 06_regression_cd4-0.1_ml.R
NREP <- 1000; NCPU <- 6
lmm_reff_after_art_B <- fit_cd4_model("cd4", cd4_model_form, cd4_after_art_subtypes_B_pts2, glue("30-B-cd4-POST-ART"), do_boot=T)
# PT40: 17 (se: 67), PT69: -16 (se: 78), PT133: 1 (se: NA)
hist(lmm_reff_after_art_B$table_pt_estims_p$phylotype_coef)
lmm_reff_after_art_B$table_pt_estims_p$phylotype[lmm_reff_after_art_B$table_pt_estims_p$phylotype_coef > 10] # only PT40 from VOIs here
lmm_reff_after_art_B$table_pt_estims_p$phylotype[lmm_reff_after_art_B$table_pt_estims_p$phylotype_coef > 20] #"84"  "148" "39"  "117" "27"
sort(lmm_reff_after_art_B$table_pt_estims_p$p_value) # all p-values >0.05

# model 2: Bayes fixed effect model of combined VOIs (PT40+PT69+PT133) against backbone
subtype_b_vois_ids_only <- c(40,69,133)
cd4_after_art_subtypes_B_pts2_feff <- cd4_after_art_subtypes_B_pts2
cd4_after_art_subtypes_B_pts2_feff$phylotype <- ifelse(cd4_after_art_subtypes_B_pts2_feff$phylotype == as.integer(backbone_cl_control[[1, 4]]), yes = as.character(backbone_cl_control[[1, 4]]),
																																																							no = ifelse(cd4_after_art_subtypes_B_pts2_feff$phylotype %in% subtype_b_vois_ids_only, yes = "VOI", no = "non-VOI"))
cd4_after_art_subtypes_B_pts2_feff2 <- cd4_after_art_subtypes_B_pts2_feff
cd4_after_art_subtypes_B_pts2_feff <- cd4_after_art_subtypes_B_pts2_feff[cd4_after_art_subtypes_B_pts2_feff$phylotype != "non-VOI",]
cd4_after_art_subtypes_B_pts2_feff$phylotype <- as.factor(cd4_after_art_subtypes_B_pts2_feff$phylotype)
cd4_after_art_subtypes_B_pts2_feff <- cd4_after_art_subtypes_B_pts2_feff %>% mutate(phylotype = relevel(phylotype, ref=backbone_cl_control[[1,4]]), age_group = relevel(age_group, ref="30-39"), sexid = relevel(sexid, ref="Male"), exposureid = relevel(exposureid, ref="Homo/bisexual") )
cd4_model_form_with_intercept_pt <- " ~ years_since_1cd4 + age_group + years_since_1cd4 * sexid + years_since_1cd4:exposureid + years_since_1cd4:phylotype + phylotype + (years_since_1cd4 | patientindex)"
# IMPORTANT: load function "fit_cd4_model" from 06_regression_cd4-0.2_bayes.R
ITER <- 10000; WARMUP50 <- ITER/2; THIN <- 10; CORES_CHAINS <- 2
priors_slope_interc_post_treat <- c(0, 20, 0, 100)
bmm_feff_after_art_B <- fit_bayes_cd4_model_fixed_eff_vois_vs_backbone_vs_nv("cd4", cd4_model_form_with_intercept_pt,
																																																					WARMUP50, cd4_after_art_subtypes_B_pts2_feff, change_prior=TRUE, priors_slope_interc = priors_slope_interc_post_treat)
saveRDS(bmm_feff_after_art_B, glue("{RDS_PATH}/bmm_feff_after_art_B.rds"))

tab_post_art <- bmm_feff_after_art_B$all_covar
tab_post_art <- tab_post_art %>% dplyr::select( phylotype_coef, Q2.5, Estimate, Q97.5, Est.Error )
tab_post_art <- tab_post_art %>% mutate(across(c(Q2.5, Estimate, Q97.5, Est.Error), ~ round(.x, 3)))
tab_post_art <- tab_post_art %>% arrange(desc(Estimate))
options(scipen=999)
write.csv( tab_post_art, file=glue("{RESULTS_PATH}/tables/tableS15.csv"), quote=F, row.names = F )

# CIs too big and sample sizes too small, so just plot slopes labelled by VOI, non-VOI and backbone
cd4_after_art_subtypes_B_pts2_feff2$phylotype <- as.character(cd4_after_art_subtypes_B_pts2_feff2$phylotype)
cd4_after_art_subtypes_B_pts2_feff2$phylotype[cd4_after_art_subtypes_B_pts2_feff2$phylotype=="153"] <- "Backbone"
cd4_after_art_subtypes_B_pts2_feff2$phylotype <- as.factor(cd4_after_art_subtypes_B_pts2_feff2$phylotype)
cd4_after_art_subtypes_B_pts2_feff2$age_group <- factor(cd4_after_art_subtypes_B_pts2_feff2$age_group, levels = c("0-29", "30-39", "40-49", "50-59", "60+"))
cd4_after_art_subtypes_B_pts2_feff2 <- cd4_after_art_subtypes_B_pts2_feff2[!is.na(cd4_after_art_subtypes_B_pts2_feff2$age_group),]

transpar <- 0.5
age_group_palette <- c("0-29" = alpha("#E69F00", transpar),"30-39" = alpha("#CC79A7", transpar),"40-49" = alpha("#009E73", transpar),	"50-59" = alpha("#F0E442", transpar),"60+" = alpha("#0072B2", transpar))
f_postart_1 <- ggplot(cd4_after_art_subtypes_B_pts2_feff2, aes(x = cd4_decimal_date, y = cd4, group=patientindex, color=age_group)) + #group=patientindex, color=phylotype
	geom_smooth(method = "lm", fill = NA) + facet_wrap(~phylotype) + scale_color_manual(values = age_group_palette, name = "Age group") + theme_classic() +
	ylim(c(0,2500)) + ylab("CD4 count")+ xlab("Year") + 
	theme(axis.title = element_text(size = 10, family = helv, face = "bold"),axis.text = element_text(color="black"), 
							strip.text = element_text(size = 10, family = helv, face = "bold"),legend.title = element_text(family=helv, size=12), legend.text=element_text(family=helv,size=10),
							strip.background = element_blank(), strip.placement = "outside", strip.text.x = element_text(family=helv,size=12))

# distribution of slopes for backbone, non-voi, and VOI
# calculate slopes
cd4_after_art_subtypes_B_pts2_feff2_dt <- as.data.table(cd4_after_art_subtypes_B_pts2_feff2)
slopes_cd4_bb_nv_vois <- cd4_after_art_subtypes_B_pts2_feff2_dt[,list(intercept=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[2], phylotype=phylotype), by=patientindex] |> as.data.frame()
slopes_postart_pal <- c( "VOI"="#E69F00", "non-VOI"="black", "Backbone"="#999999" )
f_postart_2 <- ggplot(slopes_cd4_bb_nv_vois, aes(x = slope, color = phylotype, fill = phylotype)) + geom_density(alpha = 0.25) + scale_x_continuous(limits=c(-5000,5000), breaks=seq(from=-5000, to=5000, by=1000)) +
	theme_classic() + labs(x = "CD4 slopes",y = "Density",color = "Phylotype",fill = "Phylotype") + scale_color_manual(values=slopes_postart_pal, name="VOI status") + scale_fill_manual(values=slopes_postart_pal, name="VOI status") +
	theme(axis.title = element_text(size = 10, family = helv, face = "bold"),axis.text = element_text(color="black"), axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_text(family=helv, size=12), legend.text=element_text(family=helv,size=10))
f_postart_2
pl_postart <- ggarrange(f_postart_1, f_postart_2, nrow=2, ncol=1, labels=c("A","B"), font.label=list(family=helv, color="black",size=10), heights=c(2,1))
ggsave(plot=pl_postart, file=glue("{RESULTS_PATH}/figs/figS16.eps"), device=cairo_ps, dpi=600, width=8, height=9, bg="white")
ggsave(plot=pl_postart, file=glue("{RESULTS_PATH}/figs/figS16.jpg"), dpi=600, width=8, height=9, bg="white")
