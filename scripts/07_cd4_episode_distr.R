libs_load <- c("ggplot2","dplyr", "lubridate", "ggsci", "ggpubr", "data.table","cowplot","grid", "glue") #"scales"
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"

min_cl_size_choices <- c(30, 50, 100) #250,500
tree_names <- c("A_A1","CRF_02_AG","C","B")

demog_md_subtype_match <- readRDS(glue("{RDS_PATH}/demog_md_subtype_match.rds"))
demog_md_subtype_match <- demog_md_subtype_match[ (demog_md_subtype_match$status=="Naïve") & (demog_md_subtype_match$rega3subtype %in% c("A (A1)","CRF 02_AG","C","B")) & (demog_md_subtype_match$exposureid != "Not known"), ]

demog_md_subtype_match$hiv_diag_decimal_date <- decimal_date(as.Date(demog_md_subtype_match$hivpos_ymd))
demog_md_subtype_match <- demog_md_subtype_match %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y))
demog_md_subtype_match <- demog_md_subtype_match %>% mutate(
	age_group = dplyr::case_when(age_diag<=29 ~ "0-29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
	age_group = factor(age_group,level = c("0-29","30-39","40-49","50-59","60+")))
demog_md_subtype_match$artstart_decimal_date <- decimal_date( as.Date(gsub("\\/", "15", demog_md_subtype_match$artstart_my), "%m%d%Y") )

extracted_clusters <- readRDS("rds/extracted_clusters.rds")
cd4_regr_subtypes <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))# [[1,4]]$cd4_df_filt # previously readRDS("rds/cd4_regr_subtypes.rds")

extr_clust_match_demog <- extr_clust_match_demog_df <- extr_clust_match_cd4 <- extr_clust_match_cd4_df <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		for(k in 1:length(extracted_clusters[[i,j]][[1]])) {
			print(glue("{min_cl_size_choices[i]}-{tree_names[j]}-PT{k}"))
			# Extract metadata for all patients within each phylotype (backbone included)
			if(!is.null( extracted_clusters[[i,j]][[1]][[k]]$tip.label)) {
				extr_clust_match_demog[[i,j]][[k]] <- demog_md_subtype_match[paste0("t.",demog_md_subtype_match$testindex) %in%  extracted_clusters[[i,j]][[1]][[k]]$tip.label, ]
				extr_clust_match_demog[[i,j]][[k]]$phylotype <- k
				# Extract metadata for patients with pre-ART CD4 (backbone included)
				extr_clust_match_cd4[[i,j]][[k]] <- cd4_regr_subtypes[[i,j]]$cd4_df_filt[ cd4_regr_subtypes[[i,j]]$cd4_df_filt$phylotype == k, ]
			}
		}
		extr_clust_match_demog_df[[i,j]] <- rbindlist(extr_clust_match_demog[[i,j]])
		extr_clust_match_cd4_df[[i,j]] <- rbindlist(extr_clust_match_cd4[[i,j]])
	}
}

d_demog <- extr_clust_match_demog_df[[1,4]] #mcs=30, subtype B
d_demog_ages <- d_demog[!is.na(d_demog$age_group),]
d_cd4 <- extr_clust_match_cd4_df[[1,4]]
d_cd4_ages <- d_cd4[!is.na(d_cd4$age_group),]
# adjust age group
levels(d_cd4_ages$age_group) <- c(levels(d_cd4_ages$age_group), "0-29")
d_cd4_ages$age_group[d_cd4_ages$age_group == "<29"] <- "0-29"
d_cd4_ages$age_group <- factor(d_cd4_ages$age_group, levels = c("0-29","30-39","40-49","50-59","60+"))

# CD4 slopes for each phylotype coloured by age_group
for(i in 1:length(unique(d_cd4_ages$phylotype))) {
	d_cd4_ages %>% filter(phylotype==i) %>% ggplot(aes(x = cd4_decimal_date, y = cd4, group=patientindex, color=age_group)) + #color=patientindex
		geom_smooth(method = "lm", fill = NA) + theme_classic() + coord_cartesian(ylim=c(0,1000)) + #+ theme(legend.position="none")
		ylab("CD4 counts")+ xlab("Year") + theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm")) +
		geom_hline(yintercept=350, linetype="longdash", color="black", size=0.3)
	system(glue("mkdir -p {RESULTS_PATH}/07_cd4regr_points/30-B/"))
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/07_cd4regr_points/30-B/phylotype_{i}.png"), dpi=600, width=5, height=4) )
}

###
# fig. 2B, average decline of VOIs with uncertainty
###

cd4_measur_palette <- c("1"="#56B4E9","2"="#377EB8","3"="#D55E00","4"="#009E73", "5"="#E69F00")

### BEGIN FIGURE S8 ###
# Stacked bar of nth measurements per year considering all phylotypes 
d_cd4$year <- lubridate::year(as.Date(date_decimal(d_cd4$cd4_decimal_date)))
d_cd4$year <- as.character(d_cd4$year)
print(length(unique(d_cd4$patientindex))) #8946
print(nrow(d_cd4)) #31503
d2 <- d_cd4 %>% group_by(measurement_id, year) %>% summarise(meas_year_n = n()) # measurement for each patient and group by it
d2$measurement_id <- as.factor(d2$measurement_id)
# Plot distribution of CD4 measurements order over time
options(scipen=999)
leg2 <- theme(text=element_text(family="Helvetica"), axis.text=element_text(size=7, color="black"), axis.title=element_text(size=10, color="black", face="bold"))
s8_1 <- ggplot(d2, aes(fill=measurement_id, y=year, x=meas_year_n)) + #scale_fill_viridis(discrete=T, name="CD4 measurement order") +
	scale_fill_manual(values=cd4_measur_palette, name="CD4 measurement episode") +
	geom_bar(position = position_stack(reverse = TRUE), stat="identity") + labs(x="CD4 measurements", y="Year") + theme_classic() + #position="stack"
	theme(axis.text=element_text(colour="black"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'),legend.title = element_text(size=10), legend.position="top") +
	leg2 + scale_x_continuous(expand = c(0,0)) + scale_y_discrete(expand = c(0,0))
s8_2 <- readRDS(glue("{RDS_PATH}/pl_c.rds"))
ggarrange(s8_1, s8_2 + leg2, nrow=2, ncol=1, labels=c("A","B"), font.label=list(family=helv, color="black",size=10), heights=c(3.5,5))
ggsave(file=glue("{RESULTS_PATH}/figs/figS8.eps"), dpi=600, width=10, height=12, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS8.jpg"), dpi=600, width=10, height=12, bg="white")
### END FIGURE S8 ###

# Same as above but 1 plot per phylotype
d3 <- d_cd4 %>% group_by(measurement_id, year, phylotype) %>% summarise(meas_year_n = n())
d3$measurement_id <- as.factor(d3$measurement_id)

for(i in 1:154) {
	if(i == 6 || i == 153) next
	d3 %>% filter(phylotype==i) %>% ggplot(aes(fill=measurement_id, y=year, x=meas_year_n)) + scale_fill_manual(values=cd4_measur_palette, name="CD4 measurement episode") +
		geom_bar(position = position_stack(reverse = TRUE), stat="identity") + labs(x="CD4 measurements", y="Year") + theme_bw() + 
		theme(axis.text=element_text(colour="black"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'),legend.title = element_text(size=10), legend.position="top")
	#system(glue("mkdir -p {RESULTS_PATH}/07_cd4regr_plots/30-B/stackbar_phylotypes/"))
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/07_cd4regr_plots/30-B/stackbar_phylotypes/phylotype_{i}_measur_years.png"), dpi=600, width=8, height=6) )
}

### BEGIN FIGURE S9A-C ###
s3_pl <- list()
c <- 1
subtype_b_vois_ids_only <- c(40,69,133)
for(i in subtype_b_vois_ids_only) {
	s3_pl[[c]] <- d3 %>% filter(phylotype==i) %>% ggplot(aes(fill=measurement_id, y=year, x=meas_year_n)) +
		geom_bar(position = position_stack(reverse = TRUE), stat="identity") + theme_classic() +
		scale_fill_manual(values=cd4_measur_palette, name="CD4 measurement\nepisode") + 
		#scale_x_continuous(expand = c(0,0), limits=c(0, 70), breaks=c(5,10,20,30,40,50,60,70)) + 
		scale_x_continuous(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
		theme(legend.position="top", legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'), 
								legend.title = element_text(family=helv, size=8), legend.text=element_text(size=7), plot.margin = margin(0.5, 0.25, 0.25, 0.25, "cm"), axis.title.x=element_blank(), axis.title.y=element_blank(), 
								axis.text.y=element_text(size=6)) + guides(fill = guide_legend(nrow = 2)) + leg #+
	c = c+1
}

s3_pl[[1]] <- s3_pl[[1]] + scale_x_continuous(expand = c(0,0), limits=c(0, 15), breaks=c(0,5,10,15))
s3_pl[[2]] <- s3_pl[[2]] + scale_x_continuous(expand = c(0,0), limits=c(0, 15), breaks=c(0,5,10,15))
s3_pl[[3]] <- s3_pl[[3]] + scale_x_continuous(expand = c(0,0), limits=c(0, 15), breaks=c(0,5,10,15))

### FIGURE S9D ###

age_group_palette <- c("0-29" = "#E69F00","30-39" = "#CC79A7", "40-49" = "#009E73", "50-59" = "#F0E442", "60+" = "#0072B2")

subtype_b_vois_ids_only <- c(40,69,133)
d_cd4_ages_filtered <- d_cd4_ages %>%
	filter(phylotype %in% subtype_b_vois_ids_only) %>%
	mutate(phylotype_label = factor(phylotype, 
																																	labels = c("PT.B.40.UK", "PT.B.69.UK", "PT.B.133.UK")))

# Create the plot using facet_wrap
f1c <- d_cd4_ages_filtered %>%
	ggplot(aes(x = cd4_decimal_date, y = cd4, group = patientindex, color = age_group)) +
	stat_smooth(geom = "line", method = "lm", alpha = 0.9, se = FALSE, linewidth = 0.75) +
	geom_hline(yintercept = 350, linetype = "longdash", color = "black", size = 0.3) +
	theme_classic() +
	coord_cartesian(ylim = c(0, 1000)) +
	scale_color_manual(values = age_group_palette, name = "Age group") +
	guides(color = guide_legend(ncol = 5)) +
	scale_x_continuous(limits = c(1995, 2015), breaks = c(1995, 2000, 2005, 2010, 2015)) +
	scale_y_continuous(limits=c(0,1000), breaks=seq(from=0,to=1000,by=100)) +
	facet_wrap(~ phylotype_label, nrow = 1) + # Use facet_wrap to create the subplots
	theme(plot.margin = margin(1.5, 0.2, 1.5, 0.2, "cm"),axis.title.x = element_text(size = 10, family = helv, face = "bold"),
							axis.title.y = element_text(size = 10, family = helv, face = "bold"),axis.text.x = element_text(angle = 60, hjust = 1),
							axis.text = element_text(color="black"), strip.text = element_text(size = 10, family = helv, face = "bold"), legend.position = c(0.5,0.9),
							legend.title = element_text(family=helv, size=8), legend.text=element_text(family=helv,size=6),
							strip.background = element_blank(), strip.placement = "outside", strip.text.x = element_text(family=helv,size=12)) +
	guides(color = guide_legend(nrow = 2)) +
	labs(x = "CD4 measurement year",y = bquote(bold("CD4 count (cells /"~mm^3~")")))

f1c
#saveRDS( f1c, glue("{RDS_PATH}/f1c.rds") )

s3 <- ggarrange(s3_pl[[1]], s3_pl[[2]], s3_pl[[3]], common.legend = T, nrow=3, ncol=1, labels = c("A","B","C"), font.label=list(family="Helvetica", color="black",size=10)) #lg2, labels=subtype_b_vois, 
s3_1 <- annotate_figure(s3, left = text_grob("Year", rot = 90, vjust = 1, size=10, family = helv, face="bold"), bottom = text_grob("CD4 measurements", size=10, family = helv, face="bold"))
s3_2 <- ggarrange(f1c, common.legend = F, nrow=1, ncol=1, labels = c("D"), font.label=list(family="Helvetica", color="black",size=10))
ggarrange(s3_1, s3_2, common.legend = F, nrow=1, ncol=2, widths=c(1,2))
ggsave(file=glue("{RESULTS_PATH}/figs/figS9.eps"), device=cairo_ps, dpi=600, width=7, height=6.5, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS9.jpg"), dpi=600, width=7, height=6.5, bg="white")
### END FIGURE S9 ###

# Begin part of table 1: CD4 summary stats (including indivs with 0 and 0 measurements, which are not in analysis)

# pick up cd4 df back to get indivs with 0 or 1 meas
cd4_counts_all_preart <- list()
for(i in 1:length(subtype_choices)) {
	print(subtype_choices[i])
	cd4_counts_all_preart[[i]] <- remove_cd4_after_art(demog_md_subtype_match, subtype_choices[i], cd4_md)
}
# all pre-art meas
cd4_counts_all_preart_df <- rbind(cd4_counts_all_preart[[1]][[3]], cd4_counts_all_preart[[2]][[3]], cd4_counts_all_preart[[3]][[3]], cd4_counts_all_preart[[4]][[3]])
# corresp patient and sample
cd4_counts_all_preart_corresp_df <- rbind(cd4_counts_all_preart[[1]][[2]], cd4_counts_all_preart[[2]][[2]], cd4_counts_all_preart[[3]][[2]], cd4_counts_all_preart[[4]][[2]])
# tips included
tips_incl <- unlist(lapply(timetrees, function(tree) tree$tip.label))
tips_incl_df <- data.frame(testindex=tips_incl)
# make sure get all testindex for incl tips
cd4_counts_all_preart_corresp_tips_incl_df <- left_join(tips_incl_df, cd4_counts_all_preart_corresp_df, by="testindex")
demog_md_subtype_match$testindex <- paste0("t.",demog_md_subtype_match$testindex)
# get all assoc patientindex
cd4_counts_all_preart_corresp_tips_incl_df <- cd4_counts_all_preart_corresp_tips_incl_df %>%
	left_join(demog_md_subtype_match %>% select(testindex, patientindex_match = patientindex), by = "testindex") %>%
	mutate(patientindex = ifelse(is.na(patientindex), patientindex_match, patientindex)) %>%
	select(-patientindex_match) 
cd4_counts_all_preart_corresp_tips_incl_all_meas_df <- left_join(cd4_counts_all_preart_df, cd4_counts_all_preart_corresp_tips_incl_df, by="patientindex", relationship = "many-to-many")
# 75283 rows + 2258 (0 meas) = 77541 meas including 
# count number of CD4 measurements per patient
cd4_counts_all_quantities <- cd4_counts_all_preart_corresp_tips_incl_all_meas_df %>% group_by(patientindex) %>% summarise(n_obs = n())
# add back ones with zero meas
cd4_counts_all_quantities2 <- data.frame(patientindex = seq(from = max(cd4_counts_all_quantities$patientindex),length.out = 40888-nrow(cd4_counts_all_quantities)),n_obs=0)
cd4_counts_all_quantities3 <- rbind(cd4_counts_all_quantities, cd4_counts_all_quantities2)

mean(cd4_counts_all_quantities3$n_obs) #1.84
sd(cd4_counts_all_quantities3$n_obs) #1.36
median(cd4_counts_all_quantities3$n_obs) # 1
calc_iqr(cd4_counts_all_quantities3$n_obs) #1-4
table(cd4_counts_all_quantities3$n_obs); hist(cd4_counts_all_quantities3$n_obs)
# 0     1     2     3     4     5     6     8    10 
# 2258 23471  3729  3003  6840  1571     3     9     4 
# End part of table 1

# Table S2: demog variable distributions after CD4 filters

# Gender
demog_md_subtype_match_sex_cd4 <- cd4_before_art_subtypes_all_df_ts2
demog_md_subtype_match_sex_cd4$subtype2 <- ifelse(demog_md_subtype_match_sex_cd4$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_sex_cd4$rega3subtype), "Others")
demog_md_subtype_match_sex_cd4 <- demog_md_subtype_match_sex_cd4 %>% group_by(sexid, subtype2) %>% summarise(n=n())
demog_md_subtype_match_sex_cd4 %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})"))
#View(demog_md_subtype_match_sex %>% group_by(sexid, subtype2) %>% summarise(n=n())) #%>% mutate(n_perc=glue("{n} ({n*100/})"))

# Ethnicity
demog_md_subtype_match_eth_cd4 <- cd4_before_art_subtypes_all_df_ts2
demog_md_subtype_match_eth_cd4$subtype2 <- ifelse(demog_md_subtype_match_eth_cd4$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_eth_cd4$rega3subtype), "Others")
demog_md_subtype_match_eth_cd4 <- demog_md_subtype_match_eth_cd4 %>% mutate(ethnicityid2 = case_when(
	(ethnicityid=="Black-Caribbean") | (ethnicityid=="Black-African") | (ethnicityid=="Black-other/unspecified") ~ "Black-Caribbean / African / other",
	(ethnicityid=="Indian/Pakistani/Bangladeshi") ~ "Indian/Pakistani/Bangladeshi", 
	(ethnicityid=="Other Asian/Oriental") ~ "Other Asian/Oriental", 
	(ethnicityid=="White") ~ "White",
	TRUE ~ "Other/mixed/NA"))
demog_md_subtype_match_eth_cd4 <- demog_md_subtype_match_eth_cd4 %>% group_by(ethnicityid2, subtype2) %>% summarise(n=n())
demog_md_subtype_match_eth_cd4 %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})"))

# Exposure (risk group)
demog_md_subtype_match_exposure_cd4 <- cd4_before_art_subtypes_all_df_ts2
demog_md_subtype_match_exposure_cd4$subtype2 <- ifelse(demog_md_subtype_match_exposure_cd4$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_exposure_cd4$rega3subtype), "Others")
demog_md_subtype_match_exposure_cd4 <- demog_md_subtype_match_exposure_cd4 %>% mutate(exposureid2 = case_when(
	(exposureid=="IDU") | (exposureid=="Blood products") ~ "IDU and blood products",
	(exposureid=="Homo/bisexual") ~ "Homo/bisexual", (exposureid=="Heterosexual") ~ "Heterosexual", TRUE ~ "Other/NA")) #(exposureid=="Not known") ~ "Not known"
demog_md_subtype_match_exposure_cd4 <- demog_md_subtype_match_exposure_cd4 %>% group_by(exposureid2, subtype2) %>% summarise(n=n())
demog_md_subtype_match_exposure_cd4 %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})"))

# Region
demog_md_subtype_match_region_cd4 <- cd4_before_art_subtypes_all_df_ts2
demog_md_subtype_match_region_cd4$subtype2 <- ifelse(demog_md_subtype_match_region_cd4$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_region_cd4$rega3subtype), "Others")
demog_md_subtype_match_region_cd4 <- demog_md_subtype_match_region_cd4 %>% group_by(PHE_regiondiagnosed, subtype2) %>% summarise(n=n())
View(demog_md_subtype_match_region_cd4 %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})")))

# Age at diagnosis
demog_md_subtype_match_age_cd4 <- cd4_before_art_subtypes_all_df_ts2
demog_md_subtype_match_age_cd4$subtype2 <- ifelse(demog_md_subtype_match_age_cd4$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_age_cd4$rega3subtype), "Others")
demog_md_subtype_match_age_cd4$hiv_diag_decimal_date <- decimal_date(as.Date(demog_md_subtype_match_age_cd4$hivpos_ymd))
demog_md_subtype_match_age_cd4 <- demog_md_subtype_match_age_cd4 %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y))
demog_md_subtype_match_age_cd4 <- demog_md_subtype_match_age_cd4 %>% mutate(
	age_group = dplyr::case_when(age_diag<=29 ~ "<29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
	age_group = factor(age_group,level = c("<29","30-39","40-49","50-59","60+")))
demog_md_subtype_match_age_cd4 <- demog_md_subtype_match_age_cd4 %>% group_by(age_group, subtype2) %>% summarise(n=n())
View(demog_md_subtype_match_age_cd4 %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})")))

# Table S3: VLs and CD4 measurements of individual PLWH pre-ART for defined VOIs
# VL
subtype_b_vois_ids_only <- c(40,69,133)
vls_indiv_vois <- readRDS(glue("{RDS_PATH}/vl_subtypes.rds"))[[4]][[1]]
vls_vois <- readRDS(glue("{RDS_PATH}/vl_subtypes_comb1.rds"))[[1,4]]
vls_vois <- vls_vois[vls_vois$cluster %in% subtype_b_vois_ids_only,]
vls_vois$cluster <- droplevels(vls_vois$cluster)
table(vls_vois$cluster) #81 patients
nrow(vls_vois); length(unique(vls_vois$taxon))
vls_vois_combined <- inner_join( vls_vois, vls_indiv_vois, by="patientindex" ) #320 meas
vls_vois_combined$log_vl <- round(vls_vois_combined$log_vl,3)
vls_vois_combined <- vls_vois_combined %>% dplyr::select(rega3subtype, cluster, taxon, hivpos_year, artstart_my, vl_date_my, log_vl) #artstart_ymd, vl_date_ymd
nrow(vls_vois_combined); length(unique(vls_vois_combined$taxon)) #307, 81

# CD4
cd4_indiv_vois <- readRDS(glue("{RDS_PATH}/residuals_removal_mx.rds"))[[1,4]]$cd4_df_filt
cd4_indiv_vois <- cd4_indiv_vois %>% ungroup() %>% dplyr::select(rega3subtype, cluster, taxon, hivpos_year, artstart_my, cd4_date_my, cd4) #artstart_ymd, cd4_date_ymd
print(nrow(cd4_indiv_vois)); print(length(unique(cd4_indiv_vois$taxon)))
cd4_indiv_vois <- cd4_indiv_vois[cd4_indiv_vois$cluster %in% subtype_b_vois_ids_only,] # 31503 to 258
cd4_indiv_vois$cluster <- droplevels(cd4_indiv_vois$cluster)
table(cd4_indiv_vois$cluster)
# measurements
# 40  69 133 
# 112  80  66
nrow(cd4_indiv_vois); length(unique(cd4_indiv_vois$taxon)) #258, 74

# merge both 307, 258
ij_vl_cd4_vois <- inner_join(vls_vois_combined, cd4_indiv_vois, by=c("taxon"="taxon", "vl_date_my"="cd4_date_my")) #155
fj_vl_cd4_vois <- full_join(vls_vois_combined, cd4_indiv_vois, by=c("taxon"="taxon", "vl_date_my"="cd4_date_my")) #410
fj_vl_cd4_vois$cluster.x <- as.integer(as.character(fj_vl_cd4_vois$cluster.x))
fj_vl_cd4_vois$cluster.y <- as.integer(as.character(fj_vl_cd4_vois$cluster.y))
fj_vl_cd4_vois$phylotype <- ifelse( is.na(fj_vl_cd4_vois$cluster.x), yes=fj_vl_cd4_vois$cluster.y, no=fj_vl_cd4_vois$cluster.x )
fj_vl_cd4_vois$subtype <- ifelse( is.na(fj_vl_cd4_vois$rega3subtype.x), yes=fj_vl_cd4_vois$rega3subtype.y, no=fj_vl_cd4_vois$rega3subtype.x )
fj_vl_cd4_vois$artstart_my <- ifelse( is.na(fj_vl_cd4_vois$artstart_my.x), yes=fj_vl_cd4_vois$artstart_my.y, no=fj_vl_cd4_vois$artstart_my.x )
fj_vl_cd4_vois <- fj_vl_cd4_vois %>% mutate(measurement_date = my(vl_date_my)) %>% arrange(phylotype, measurement_date, taxon)
# vl date is actually cd4 date when vls are empty because of full_join matching that drops that column
fj_vl_cd4_vois <- fj_vl_cd4_vois %>% dplyr::select( subtype, phylotype, taxon, artstart_my, vl_date_my, log_vl, cd4 )
colnames(fj_vl_cd4_vois) <- c("subtype", "phylotype_id", "sequence_id", "art_start_my", "measurement_date", "log10_vl", "cd4")
write.csv(fj_vl_cd4_vois, file=glue("{RESULTS_PATH}/tables/tableS3.csv"), quote=F, row.names=F, na="")