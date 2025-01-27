libs_load <- c( "ggplot2", "dplyr","ggforce", "glue", "car", "tibble", "reshape2", "ggpubr", "gtools") #"multcomp"
invisible( lapply(libs_load, library, character.only=TRUE) )

RDS_PATH="rds"
RESULTS_PATH="results"
#dir.create( RESULTS_PATH )

# Figure 2D: age comparisons
subtype_b_vois_ids_only <- c(40,69,133)
subtype_b_vois <- c("PT.B.40.UK","PT.B.69.UK","PT.B.133.UK")
backbone_cl_control <- readRDS("rds/backbone_cl_control.rds")
# below from 07 script
d_demog_age <- d_demog[d_demog$age_diag > 0,]
d_demog_age_vois <- d_demog_age[d_demog_age$phylotype %in% subtype_b_vois_ids_only,]
d_demog_age_vois$phylotype <- paste0("PT.B.",d_demog_age_vois$phylotype,".UK")
d_demog_age_backbone <- d_demog_age[d_demog_age$phylotype == backbone_cl_control[[1,4]],]
d_demog_age_backbone$phylotype <- as.character(d_demog_age_backbone$phylotype)
d_demog_age_backbone$phylotype <- "Backbone"
d_demog_age_voi_non <- rbind(d_demog_age_vois, d_demog_age_backbone)
d_demog_age_voi_non$phylotype <- factor(d_demog_age_voi_non$phylotype, levels = c("Backbone", c(subtype_b_vois))) #"PT.B.5.UK","PT.B.137.UK"
table(d_demog_age_voi_non$phylotype)
# Backbone  PT.B.40.UK  PT.B.69.UK PT.B.133.UK 
# 13376          86          70          39 

# BEGIN COUNT PROPORTION OF SAMPLES IN VOIs SINCE 2015
count_vois_since_2015 <- d_demog
count_vois_since_2015$after_2015 <- ifelse(d_demog$dbsample_date >= 2015, yes="After 2015", no="Before 2015")
count_vois_since_2015 %>% group_by(after_2015) %>% summarise(n=n()) # 2999 after, 20776 before

count_vois_since_2015_vois <- count_vois_since_2015[count_vois_since_2015$phylotype %in% subtype_b_vois,]
tb <- table(count_vois_since_2015_vois$phylotype, count_vois_since_2015_vois$after_2015)
# format and include as table S11

# Table S12: Welch's t-test on age within each VOI phylotype vs backbone
pt_voi_lbls <- subtype_b_vois
#df_age_non_vois <- d100_ages_voi_non[!(d100_ages_voi_non$phylotype %in% pt_voi_lbls),]

df_age_vois <- ttest_ages <- ttest_ages_res <- list()
z <- 1
for(i in pt_voi_lbls) {
	df_age_vois[[z]] <- d_demog_age_voi_non[d_demog_age_voi_non$phylotype == i,] #d_demog_age_voi_non[d_demog_age_voi_non$phylotype == i,]
	ttest_ages[[z]] <- t.test(df_age_vois[[z]]$age_diag, d_demog_age_backbone$age_diag, alternative="two.sided")
	ttest_ages_res[[z]] <- data.frame(phylotype=i, rows_x=nrow(df_age_vois[[z]]), estimate_x=round(ttest_ages[[z]]$estimate[[1]],3), estimate_y=round(ttest_ages[[z]]$estimate[[2]],3), p_value=ttest_ages[[z]]$p.value, sd_x=round(sd(df_age_vois[[z]]$age_diag),3), sd_y=round(sd(d_demog_age_backbone$age_diag,na.rm=T),3))
	z <- z+1
}
ttest_ages_res_list <- rbindlist(ttest_ages_res)
# adjust for multiple testing
ttest_ages_res_list_p <- ttest_ages_res_list$p_value
.adjust_methods <- function(ps, method) {
	res_adj <- p.adjust(ps, method = method)
	res_adj
}
ttest_age_bonf_adj <- .adjust_methods(ttest_ages_res_list_p, "bonferroni")
ttest_age_fdr_adj <- .adjust_methods(ttest_ages_res_list_p, "fdr")
ttest_ages_res_adj <- ttest_ages_res_list
ttest_ages_res_adj$p_adj_bonf <- ttest_age_bonf_adj
ttest_ages_res_adj$p_adj_fdr <- ttest_age_fdr_adj
ttest_ages_res_adj <- ttest_ages_res_adj %>% mutate(across(c(p_value, p_adj_bonf, p_adj_fdr), ~ signif(.x, 2)))
ttest_ages_res_adj <- ttest_ages_res_adj %>% dplyr::select( phylotype, rows_x, estimate_x, sd_x, estimate_y, sd_y, p_adj_fdr )
# only PT133 characterised by younger age of patients (30.333 vs 36.241 backbone), p_adj_fdr=0.0003811728

# Table S12 csv
write.csv(ttest_ages_res_adj, file=glue("{RESULTS_PATH}/tables/tableS12.csv"), quote=F, row.names=F)

#age_diag_pt_pal <- c("PT.B.5.UK"="#0099B4","PT.B.40.UK"="#088c06","PT.B.69.UK"="#ED0000","PT.B.133.UK"="#00468B","PT.B.137.UK"="#580618", "Backbone"="#ADB6B6")
age_diag_pt_pal <- c("PT.B.40.UK"="#56B4E9", "PT.B.69.UK"="#D55E00", "PT.B.133.UK"="#009E73")

ttest_ages_res_adj$group1 <- ttest_ages_res_adj$phylotype
ttest_ages_res_adj$group2 <- "Backbone"
ttest_ages_res_adj$.y. <- "age_diag"
ttest_ages_res_adj$p <- signif(ttest_ages_res_adj$p_adj_fdr, digits=2)
ttest_ages_res_adj$p.signif <- stars.pval(as.numeric(ttest_ages_res_adj$p))

print(mean(ttest_ages_res_adj$estimate_x)) # mean age for VOIs combined
print(unique(ttest_ages_res_adj$estimate_y)) # mean age for backbone phylotype

helv <- "Helvetica"

f2_age_bxp <- ggplot(d_demog_age_voi_non, aes(x=phylotype, y=age_diag, fill=phylotype)) + #fill=hlt, alpha=hlt
	geom_violin(width=1, alpha=0.9) +
	geom_boxplot(width=0.1, color="grey10", fill="white", alpha=1, outlier.shape=NA) + stat_pvalue_manual(ttest_ages_res_adj, label="p = {p}", size=2.5, y.position = 80, step.increase = 0.04, hide.ns = TRUE) +
	
	scale_fill_manual(values=age_diag_pt_pal, name="Phylotype") + 
	labs(x="Phylotype", y="Age at diagnosis (years)") + scale_y_continuous(expand = c(0,0), limits=c(0,90), breaks=seq(0,80,by=10)) + scale_x_discrete(expand=c(0,0)) +
	theme_classic() + theme(axis.text.x = element_text(angle = 30,hjust = 1), plot.title = element_text(hjust=0.5),
																									legend.text=element_text(size=8,family=helv), legend.title=element_text(size=8,family=helv), legend.key.height=unit(.5, "cm"),legend.key.width=unit(.5, "cm")) + leg + theme(axis.text=element_text(size=8)) + guides(color="none", fill=guide_legend(ncol=1))
f2_age_bxp
saveRDS(f2_age_bxp, glue("{RDS_PATH}/f2_age_bxp.rds"))

# Figs S11 and S14
d_demog_eth_reg <- d_demog

d_demog_eth_reg$phylotype <- as.integer(d_demog_eth_reg$phylotype) 

eth_pal <- c("White"="#377EB8","Black-Caribbean"="#E69F00","Black-African"="#D55E00","Black-other/unspecified"="#CC79A7","Indian/Pakistani/Bangladeshi"="#009E73", "Other Asian/Oriental"="#56B4E9","Other/Mixed"="#999999","Other"="#F0E442","Not known"="#000000")
reg_pal <- c("South of England"="#56B4E9","Not Known"="#000000","London"="#377EB8","North of England"="#D55E00","Northern Ireland, Scotland, and Wales"="#009E73", "Midlands and East of England"="#E69F00")
expos_pal <- c("Homo/bisexual"="#377EB8","Heterosexual"="#E69F00","IDU"="#CC79A7","Blood products"="#009E73","Other"="#999999")
sex_pal <- c("Male"="#377EB8", "Female"="#D55E00", "Not known"="#000000")

# Boxplot or bar plot for age (numeric age, not age cat), initial cd4 value (related to time to diagnosis), maybe race, geographic region
demog_vars_across_clusters <- function(df2, threshold) { #df, 
	
	# Previously had boxplots of age and initial CD4, but info repeating with tables
	
	### BEGIN FIGURE S10 ###
	# Bar plot of ethnicity
	df_eth <- df2 %>% group_by(phylotype) %>% count(ethnicityid)
	pl3 <- ggplot(df_eth, aes(x=as.factor(phylotype), y=n, fill=ethnicityid)) +
		geom_bar(position="fill", stat="identity") + scale_fill_manual(values=eth_pal, name="Ethnicity") +
		labs(x="Phylotype", y="Proportion in phylotype") + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
		theme_classic() + leg + theme(axis.text.x=element_text(size=5.5,angle = 90, vjust = 0.5, hjust=1), legend.position = "top")
	### SOME PANELS OF FIGURE S10 ARE GENERATED LATER ###
	
	# Ethnicity excluding White & unknown
	df_exc_white <- df2[df2$ethnicityid != "White" & df2$ethnicityid != "Not known",]
	df_eth_exc_white <- df_exc_white %>% group_by(phylotype) %>% count(ethnicityid)
	pl4 <- ggplot(df_eth_exc_white, aes(x=as.factor(phylotype), y=n, fill=ethnicityid)) +
		geom_bar(position="fill", stat="identity") +
		labs(x="Phylotype", y="Ethnicity") + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
		theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7))
	system(glue("mkdir -p {RESULTS_PATH}/demog_comparisons/{threshold}/"))
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/demog_comparisons/{threshold}/ethnicity_excl_white.png"), plot=pl4, dpi=600, width=12, height=8, limitsize=FALSE, bg="white") )
	
	### BEGIN FIGURE S11 ###
	# Bar plot of region diagnosed
	df_reg <- df2 %>% group_by(phylotype) %>% count(PHE_regiondiagnosed)
	pl5 <- ggplot(df_reg, aes(x=as.factor(phylotype), y=n, fill=PHE_regiondiagnosed)) +
		geom_bar(position="fill", stat="identity") + scale_fill_manual(values=reg_pal, name="Region of diagnosis") +
		labs(x="Phylotype", y="Proportion in phylotype") + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
		theme_classic() + leg + theme(axis.text.x=element_text(size=5.5,angle = 90, vjust = 0.5, hjust=1), legend.position = "top")
	### SOME PANELS OF FIGURE S11 ARE GENERATED LATER ###
	
	system(glue("mkdir -p {RESULTS_PATH}/demog_comparisons/{threshold}/"))
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/demog_comparisons/{threshold}/PHE_region_diag.png"), plot=pl5, dpi=600, width=12, height=8, limitsize=FALSE, bg="white") )
	# Region excluding unknown
	df_exc_unkn <- df2[df2$PHE_regiondiagnosed != "Not Known",]
	df_reg_exc_unkn <- df_exc_unkn %>% group_by(phylotype) %>% count(PHE_regiondiagnosed)
	pl6 <- ggplot(df_reg_exc_unkn, aes(x=as.factor(phylotype), y=n, fill=PHE_regiondiagnosed)) +
		geom_bar(position="fill", stat="identity") +
		labs(x="Phylotype", y="Region diagnosed") + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
		theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7))
	system(glue("mkdir -p {RESULTS_PATH}/demog_comparisons/{threshold}/"))
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/demog_comparisons/{threshold}/PHE_region_diag_excl_unknown.png"), plot=pl6, dpi=600, width=12, height=8, limitsize=FALSE, bg="white") )
	
	### BEGIN FIGURE S12 ###
	# Bar plot of exposureid
	df_expos <- df2 %>% group_by(phylotype) %>% count(exposureid)
	pl6 <- ggplot(df_expos, aes(x=as.factor(phylotype), y=n, fill=exposureid)) +
		geom_bar(position="fill", stat="identity") + scale_fill_manual(values=expos_pal, name="Risk group") +
		labs(x="Phylotype", y="Proportion in phylotype") + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
		theme_classic() + leg + theme(axis.text.x=element_text(size=5.5,angle = 90, vjust = 0.5, hjust=1), legend.position = "top")
	### SOME PANELS OF FIGURE S12 ARE GENERATED LATER ###
	
	### BEGIN FIGURE S13 ###
	# Bar plot of sexid
	df_sex <- df2 %>% group_by(phylotype) %>% count(sexid)
	pl7 <- ggplot(df_sex, aes(x=as.factor(phylotype), y=n, fill=sexid)) +
		geom_bar(position="fill", stat="identity") + scale_fill_manual(values=sex_pal, name="Sex") +
		labs(x="Phylotype", y="Proportion in phylotype") + scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
		theme_classic() + leg + theme(axis.text.x=element_text(size=5.5,angle = 90, vjust = 0.5, hjust=1), legend.position = "top")
	### SOME PANELS OF FIGURE S12 ARE GENERATED LATER ###
	
	list(pl3, pl5, pl6, pl7)
}

res_age_eth_reg <- demog_vars_across_clusters(d_demog_eth_reg, "30-B")
saveRDS(res_age_eth_reg, glue("{RDS_PATH}/res_age_eth_reg.rds"))

#obs_exp_pal <- c("Expected"="#00468B7F", "Observed"="#ED00007F")
obs_exp_pal <- c("Expected"="#999999", "Observed"="#000000")

options(scipen=999)

stat_tests <- function(df_eth_reg, threshold) {

	### CONTINUE FIGURE S10 ####
	### ETHNICITY
	df_eth_reg$phylotype <- as.character(df_eth_reg$phylotype)
	df_eth_reg <- df_eth_reg[df_eth_reg$ethnicityid != "Not known",]
	chisq_ethnicity <- chisq.test(df_eth_reg$ethnicityid, df_eth_reg$phylotype)
	#print(chisq_ethnicity)
	chisq_ethnicity_expected <- as.data.frame(t(chisq_ethnicity$expected))
	chisq_ethnicity_expected <- tibble::rownames_to_column(chisq_ethnicity_expected, "df_eth_reg.phylotype")
	chisq_ethnicity_expected <- reshape2::melt(chisq_ethnicity_expected)
	#print(class(chisq_ethnicity_expected))
	chisq_ethnicity_expected <- chisq_ethnicity_expected  %>% dplyr::select(variable, df_eth_reg.phylotype, value)
	colnames(chisq_ethnicity_expected) <- c("df_eth_reg.ethnicityid", "df_eth_reg.phylotype", "Expected")

	chisq_ethnicity_df_exp <- reshape2::melt(chisq_ethnicity_expected, id.vars=c("df_eth_reg.ethnicityid","df_eth_reg.phylotype"))
	chisq_ethnicity$observed <- as.data.frame(chisq_ethnicity$observed)
	colnames(chisq_ethnicity$observed) <- c("df_eth_reg.ethnicityid", "df_eth_reg.phylotype", "Observed")
	chisq_ethnicity_df_obs <- reshape2::melt(chisq_ethnicity$observed, id.vars=c("df_eth_reg.ethnicityid","df_eth_reg.phylotype"))
	chisq_ethnicity_df_exp_obs <- rbind(chisq_ethnicity_df_exp, chisq_ethnicity_df_obs)

	
	chisq_eth_pl <- list()
	c <- 1
	
	for(j in subtype_b_vois_ids_only) {
		print(j)
		chisq_ethnicity_df_exp_obs2 <- chisq_ethnicity_df_exp_obs[chisq_ethnicity_df_exp_obs$df_eth_reg.phylotype == j,]
		chisq_eth_pl[[c]] <- ggplot(chisq_ethnicity_df_exp_obs2, aes(x = df_eth_reg.ethnicityid, y = value, fill = variable)) +
			geom_col(position = position_dodge()) + #labs(x="Ethnicity", y="Number of patients in phylotype") + 
			scale_fill_manual(values=obs_exp_pal, name="Comparison") + scale_y_continuous(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
			theme_classic() + leg + theme(axis.text.x = element_text(angle = 30,hjust = 1, size=8), axis.title.x=element_blank(), axis.title.y=element_blank(), plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"))
		c <- c+1
	}
	### END FIGURE S10 ####

	### CONTINUE FIGURE S11 ####
	### PHE region diagnosed
	df_eth_reg <- df_eth_reg[df_eth_reg$PHE_regiondiagnosed != "Not Known",]
	chisq_region <- chisq.test(df_eth_reg$PHE_regiondiagnosed, df_eth_reg$phylotype)

	chisq_region_expected <- as.data.frame(t(chisq_region$expected))
	chisq_region_expected <- tibble::rownames_to_column(chisq_region_expected, "df_eth_reg.phylotype")
	chisq_region_expected <- reshape2::melt(chisq_region_expected)
	chisq_region_expected <- chisq_region_expected  %>% dplyr::select(variable, df_eth_reg.phylotype, value)
	colnames(chisq_region_expected) <- c("df_eth_reg.PHE_regiondiagnosed", "df_eth_reg.phylotype", "Expected")

	chisq_region_df_exp <- reshape2::melt(chisq_region_expected, id.vars=c("df_eth_reg.PHE_regiondiagnosed","df_eth_reg.phylotype"))
	chisq_region$observed <- as.data.frame(chisq_region$observed)
	colnames(chisq_region$observed) <- c("df_eth_reg.PHE_regiondiagnosed", "df_eth_reg.phylotype", "Observed")
	chisq_region_df_obs <- reshape2::melt(chisq_region$observed, id.vars=c("df_eth_reg.PHE_regiondiagnosed","df_eth_reg.phylotype"))
	chisq_region_df_exp_obs <- rbind(chisq_region_df_exp, chisq_region_df_obs)

	chisq_reg_pl <- list()
	c <- 1
	for(j in subtype_b_vois) {
		chisq_region_df_exp_obs2 <- chisq_region_df_exp_obs[chisq_region_df_exp_obs$df_eth_reg.phylotype == j,]
		chisq_reg_pl[[c]] <- ggplot(chisq_region_df_exp_obs2, aes(x = df_eth_reg.PHE_regiondiagnosed, y = value, fill = variable)) +
			geom_col(position = position_dodge()) + #labs(x="Region of diagnosis", y="Number of patients in phylotype") +
			scale_fill_manual(values=obs_exp_pal, name="Comparison") + scale_y_continuous(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
			theme_classic() + leg + theme(axis.text.x = element_text(angle = 50,hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank(), plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"))
		c <- c+1
	}
	### END FIGURE S11 ####
	
	### CONTINUE FIGURE S12 ####
	### exposureid
	chisq_expos <- chisq.test(df_eth_reg$exposureid, df_eth_reg$phylotype)
	
	chisq_expos_expected <- as.data.frame(t(chisq_expos$expected))
	chisq_expos_expected <- tibble::rownames_to_column(chisq_expos_expected, "df_eth_reg.phylotype")
	chisq_expos_expected <- reshape2::melt(chisq_expos_expected)
	chisq_expos_expected <- chisq_expos_expected  %>% dplyr::select(variable, df_eth_reg.phylotype, value)
	colnames(chisq_expos_expected) <- c("df_eth_reg.exposureid", "df_eth_reg.phylotype", "Expected")
	
	chisq_expos_df_exp <- reshape2::melt(chisq_expos_expected, id.vars=c("df_eth_reg.exposureid","df_eth_reg.phylotype"))
	chisq_expos$observed <- as.data.frame(chisq_expos$observed)
	colnames(chisq_expos$observed) <- c("df_eth_reg.exposureid", "df_eth_reg.phylotype", "Observed")
	chisq_expos_df_obs <- reshape2::melt(chisq_expos$observed, id.vars=c("df_eth_reg.exposureid","df_eth_reg.phylotype"))
	chisq_expos_df_exp_obs <- rbind(chisq_expos_df_exp, chisq_expos_df_obs)
	
	chisq_expos_pl <- list()
	c <- 1
	for(j in subtype_b_vois) {
		chisq_expos_df_exp_obs2 <- chisq_expos_df_exp_obs[chisq_expos_df_exp_obs$df_eth_reg.phylotype == j,]
		chisq_expos_pl[[c]] <- ggplot(chisq_expos_df_exp_obs2, aes(x = df_eth_reg.exposureid, y = value, fill = variable)) +
			geom_col(position = position_dodge()) +
			scale_fill_manual(values=obs_exp_pal, name="Comparison") + scale_y_continuous(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
			theme_classic() + leg + theme(axis.text.x = element_text(angle = 40,hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank(), plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"))
		c <- c+1
	}
	### END FIGURE S12 ####
	
	### CONTINUE FIGURE S13 ####
	### exposureid
	chisq_sex <- chisq.test(df_eth_reg$sexid, df_eth_reg$phylotype)
	
	chisq_sex_expected <- as.data.frame(t(chisq_sex$expected))
	chisq_sex_expected <- tibble::rownames_to_column(chisq_sex_expected, "df_eth_reg.phylotype")
	chisq_sex_expected <- reshape2::melt(chisq_sex_expected)
	chisq_sex_expected <- chisq_sex_expected  %>% dplyr::select(variable, df_eth_reg.phylotype, value)
	colnames(chisq_sex_expected) <- c("df_eth_reg.sexid", "df_eth_reg.phylotype", "Expected")
	
	chisq_sex_df_exp <- reshape2::melt(chisq_sex_expected, id.vars=c("df_eth_reg.sexid","df_eth_reg.phylotype"))
	chisq_sex$observed <- as.data.frame(chisq_sex$observed)
	colnames(chisq_sex$observed) <- c("df_eth_reg.sexid", "df_eth_reg.phylotype", "Observed")
	chisq_sex_df_obs <- reshape2::melt(chisq_sex$observed, id.vars=c("df_eth_reg.sexid","df_eth_reg.phylotype"))
	chisq_sex_df_exp_obs <- rbind(chisq_sex_df_exp, chisq_sex_df_obs)
	
	chisq_sex_pl <- list()
	c <- 1
	for(j in subtype_b_vois) {
		chisq_sex_df_exp_obs2 <- chisq_sex_df_exp_obs[chisq_sex_df_exp_obs$df_eth_reg.phylotype == j,]
		chisq_sex_pl[[c]] <- ggplot(chisq_sex_df_exp_obs2, aes(x = df_eth_reg.sexid, y = value, fill = variable)) +
			geom_col(position = position_dodge()) +
			scale_fill_manual(values=obs_exp_pal, name="Comparison") + scale_y_continuous(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) +
			theme_classic() + leg + theme(axis.text.x = element_text(angle = 40,hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank(), plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"))
		c <- c+1
	}
	### END FIGURE S13 ####
	
	list(chisq_eth_pl, chisq_reg_pl, chisq_expos_pl, chisq_sex_pl)
}

d30_stats <- stat_tests(d_demog_eth_reg, "30-B")
saveRDS(d30_stats, glue("{RDS_PATH}/d30_stats.rds"))

s10_bd <- d30_stats[[1]]
s11_bd <- d30_stats[[2]]
s12_bd <- d30_stats[[3]]
s13_bd <- d30_stats[[4]]

system(glue("mkdir -p {RESULTS_PATH}/figs"))
generate_figs_demog_vars <- function(res_df, idx_bd, bd, fig_label) {
	s10_1 <- ggarrange(res_df[[idx_bd]], nrow=1, ncol=1, labels = c("A"), font.label=list(family="Helvetica", color="black",size=10))
	s10_2 <- ggarrange(bd[[1]], bd[[2]], bd[[3]], nrow=1, ncol=3, labels = c("B", "C", "D"), font.label=list(family="Helvetica", color="black",size=10), common.legend = T)
	s10_2_annot <- annotate_figure(s10_2, left = text_grob("Patients in phylotype", rot = 90, vjust = 0.5, hjust=0.3, size=10, family = helv, face="bold"), bottom = text_grob("Ethnicity", size=10, family = helv, face="bold"))
	ggarrange(s10_1, s10_2_annot, ncol=1, nrow=2)
	# Fig S10
	ggsave(file=glue("{RESULTS_PATH}/figs/{fig_label}.eps"), dpi=600, width=12, height=12, bg="white")
	ggsave(file=glue("{RESULTS_PATH}/figs/{fig_label}.jpg"), dpi=600, width=12, height=12, bg="white")
}

generate_figs_demog_vars(res_age_eth_reg, 1, s10_bd, "figS11") # ethnicity
generate_figs_demog_vars(res_age_eth_reg, 2, s11_bd, "figS12") # region
generate_figs_demog_vars(res_age_eth_reg, 3, s12_bd, "figS13") # exposureid
generate_figs_demog_vars(res_age_eth_reg, 4, s13_bd, "figS14") # sex