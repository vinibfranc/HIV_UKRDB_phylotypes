libs_load <- c("ggplot2","dplyr", "lubridate", "ggsci", "scales", "ggpubr", "data.table","cowplot","grid")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"
dir.create( RESULTS_PATH )

demog_md_subtype_match <- readRDS(glue("{RDS_PATH}/demog_md_subtype_match.rds"))
demog_md_subtype_match <- demog_md_subtype_match[ (demog_md_subtype_match$status=="Naïve") & (demog_md_subtype_match$rega3subtype %in% c("A (A1)","CRF 02_AG","C","B")) & (demog_md_subtype_match$exposureid != "Not known"), ]

demog_md_subtype_match$hiv_diag_decimal_date <- decimal_date(as.Date(demog_md_subtype_match$hivpos_ymd))
demog_md_subtype_match <- demog_md_subtype_match %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y))
demog_md_subtype_match <- demog_md_subtype_match %>% mutate(
	age_group = dplyr::case_when(age_diag<=29 ~ "<29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
	age_group = factor(age_group,level = c("<29","30-39","40-49","50-59","60+")))
demog_md_subtype_match$artstart_decimal_date <- decimal_date( as.Date(gsub("\\/", "15", demog_md_subtype_match$artstart_my), "%m%d%Y") )

saveRDS(demog_md_subtype_match, "rds/demog_md_subtype_match.rds")

extracted_clusters <- readRDS("rds/extracted_clusters.rds")
cd4_regr_subtypes <- readRDS("rds/cd4_regr_subtypes.rds")

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
				extr_clust_match_cd4[[i,j]][[k]] <- cd4_regr_subtypes[[i,j]][[1]][[k]][[1]]
			}
		}
		extr_clust_match_demog_df[[i,j]] <- rbindlist(extr_clust_match_demog[[i,j]])
		extr_clust_match_cd4_df[[i,j]] <- rbindlist(extr_clust_match_cd4[[i,j]])
	}
}

d_demog <- extr_clust_match_demog_df[[1,4]]
d_demog_ages <- d_demog[!is.na(d_demog$age_group),]
d_cd4 <- extr_clust_match_cd4_df[[1,4]]
d_cd4_ages <- d_cd4[!is.na(d_cd4$age_group),]

# CD4 slopes for each phylotype coloured by age_group
for(i in 1:length(unique(d_cd4_ages$phylotype))) {
	d_cd4_ages %>% filter(phylotype==i) %>% ggplot(aes(x = cd4_decimal_date, y = cd4, group=patientindex, color=age_group)) + #color=patientindex
		geom_smooth(method = "lm", fill = NA) + theme_classic() + coord_cartesian(ylim=c(0,1000)) + #+ theme(legend.position="none")
		ylab("CD4 counts")+ xlab("Year") + theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm")) +
		geom_hline(yintercept=350, linetype="longdash", color="black", size=0.3)
	system(glue("mkdir -p {RESULTS_PATH}/08_cd4regr_points/30-B/"))
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/08_cd4regr_points/30-B/phylotype_{i}.png"), dpi=600, width=5, height=4) )
}

### BEGIN FIGURE 2D ###
#show_col(pal_lancet("lanonc", alpha = 0.7)(9))
age_group_palette <- c("0-29"="#00468BB2", "30-39"="#ED0000B2", "40-49"="#42B540B2", "50-59"="#925E9FB2", "60+"="#AD002AB2")
#display.brewer.pal(8,"Accent")

# Get median slope patients within backbone phylotype
d_ages_backbone <- d_cd4_ages[d_cd4_ages$phylotype == 153,]
d_ages_backbone <- as.data.table(d_ages_backbone)
d_ages_backbone <- d_ages_backbone[,list(intercept=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[2], phylotype=153), by=patientindex]
d_ages_backbone <- as.data.frame(d_ages_backbone)
hist(d_ages_backbone$slope)
median_slp_b <- median(d_ages_backbone$slope)
print(median_slp_b)
# Find nearest slope to the median (with measurements within 10-year window and initial CD4 below 750)
d_ages_backbone_med_pat <- d_ages_backbone
# ID: 90742 has 13 years after 1 CD4 as last measurement, so not appropriate
#d_ages_backbone_med_pat <- d_ages_backbone[ which.min(abs(median_slp_b - d_ages_backbone$slope)) && (d_ages_backbone$years_since_1_cd4 <= 5), ]
# clostest to that slope value is patient 56402 (but has initial CD4 >1000)
# next closest is 81763
d_ages_backbone_med_pat <- d_ages_backbone[ d_ages_backbone$patientindex == 81763, ]
# Get all measurements from median slope patient back
d_ages_backbone_med_pat_ms <- inner_join(d_ages_backbone_med_pat, d_cd4_ages, by="patientindex")
d_ages_backbone_med_pat_ms

age_group_palette <- c("0-29"="#F7BC00", "30-39"="#F27F0C", "40-49"="#925E9F", "50-59"="#ff4654", "60+"="#6d1723")

leg <- theme(text=element_text(family="Times New Roman"), axis.text=element_text(size=10, color="black"), axis.title=element_text(size=10, color="black", face="bold"))
ax_rem_a <- theme(axis.text.x = element_blank())
ax_rem_b <- theme(axis.text.x = element_blank(),axis.text.y = element_blank())
ax_rem_d <- theme(axis.text.y = element_blank())

#phylotypes_interest <- c("PT.B.5.UK", "PT.B.40.UK", "PT.B.69.UK", "PT.B.133.UK", "PT.B.137.UK")
s1_pl <- list()
c <- 1
for(i in spvl_cd4_vois) { #d_ages_comb_vois_non
	s1_pl[[c]] <- d_cd4_ages %>% filter(phylotype==i) %>%
		ggplot() +
		stat_smooth(data = d_cd4_ages %>% filter(phylotype == i), aes(x = cd4_decimal_date, y = cd4, group=patientindex, color=age_group), geom="line", method="lm", alpha=0.9, se=FALSE, linewidth=0.45) +
		theme_classic() + coord_cartesian(ylim=c(0,1000)) +
		scale_color_manual(values=age_group_palette, name="Age group") + scale_x_continuous(limits=c(1995, 2015), breaks=c(1995, 2000, 2005, 2010, 2015)) +
		theme(plot.margin = margin(0.5, 0.2, 0.2, 0.2, "cm"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(angle = 60,hjust = 1)) + leg +
		geom_hline(yintercept=350, linetype="longdash", color="black", size=0.3)
	c = c+1
}

tnr <- "Times New Roman"
system(glue("mkdir -p {RESULTS_PATH}/figs/"))
lg <- as_ggplot(get_legend(s1_pl[[2]])) #cowplot::
for(i in 1:length(s1_pl)) { s1_pl[[i]] <- s1_pl[[i]] + theme(legend.position = "none") }
s1 <- ggarrange(s1_pl[[1]], s1_pl[[2]], s1_pl[[3]], s1_pl[[4]], s1_pl[[5]], lg, nrow=2, ncol=3, labels = c("PT.B.5.UK","PT.B.40.UK","PT.B.69.UK","PT.B.133.UK","PT.B.137.UK"), font.label=list(family="Times New Roman", color="black",size=10, face="plain")) #common.legend = T, legend="none"
f1c <- annotate_figure(s1, left = text_grob(bquote(bold("CD4 count (cells /"~mm^3~")")), rot = 90, vjust = 1, size=10, family = tnr, face="bold"), bottom = text_grob("CD4 measurement year", size=10, family = tnr, face="bold"))
f1c
### END FIGURE 2D ###

### BEGIN FIGURE 2B ###
# Get slopes for each patientd
d_cd4_ages_sl <- as.data.table(d_cd4_ages)
cd4_d_ages_list <- cd4_d_ages_list_comb <- list()

c <- 1
for(i in spvl_cd4_vois) { #seq_along()
	cd4_d_ages_list[[c]] <- d_cd4_ages_sl[d_cd4_ages_sl$phylotype == i,]
	#print(cd4_d_ages_list[[c]])
	cd4_d_ages_list[[c]] <- cd4_d_ages_list[[c]][,list(intercept=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[2], phylotype=i), by=patientindex]
	#print(head(cd4_d_ages_list[[c]]))
	median_slp <- median(cd4_d_ages_list[[c]]$slope)
	print(median_slp)
	cd4_d_ages_list_comb[[c]] <- as.data.frame(cd4_d_ages_list[[c]])
	# Find nearest slope to the median
	#cd4_d_ages_list_comb[[c]] <- cd4_d_ages_list_comb[[c]][cd4_d_ages_list_comb[[c]]$slope == median_slp,]
	cd4_d_ages_list_comb[[c]] <- cd4_d_ages_list_comb[[c]][which.min(abs(median_slp-cd4_d_ages_list_comb[[c]]$slope)),]
	#print(cd4_d_ages_list_comb[[c]])
	
	cd4_d_ages_list_comb[[c]] <- inner_join(cd4_d_ages_list_comb[[c]], d_cd4_ages_sl, by="patientindex")
	#print(cd4_d_ages_list_comb[[c]])
	
	f2_pl[[c]] <- ggplot(cd4_d_ages_list_comb[[c]], aes(x = cd4_decimal_date, y = cd4, group=age_group, color=age_group)) +
			geom_smooth(method = "glm", level=0.95, aes(fill=age_group), alpha=0.2) + theme_classic() + coord_cartesian(ylim=c(0,1000)) +
			#scale_x_continuous(limits=c(1995, 2015), breaks=c(1995, 2000, 2005, 2010, 2015)) + scale_x_continuous(limits=c(0,15), breaks=seq(0,15,3)) +
			coord_cartesian(ylim=c(0,750)) +
			scale_color_manual(values=age_group_palette, name="Age group") + scale_fill_manual(values=age_group_palette, name="Age group") +
			#ylab("CD4 counts")+ xlab("Years since first CD4") +
			theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"), axis.title.x=element_blank(), axis.title.y=element_blank()) + leg + #guides(fill="none") +
			geom_hline(yintercept=350, linetype="longdash", color="black", size=0.3)
	
	c <- c+1
}

# Selected patient with 'median slope' within each of the phylotypes of interest
f2_df <- rbindlist(cd4_d_ages_list_comb)
f2_df$phylotype.x <- paste0("PT.B.",f2_df$phylotype.x,".UK")

d_ages_backbone_med_pat_ms$phylotype.x <- as.character(d_ages_backbone_med_pat_ms$phylotype.x)
d_ages_backbone_med_pat_ms$phylotype.x <- "Backbone"

comb_cd4_vois_non <- rbind(f2_df, d_ages_backbone_med_pat_ms)
comb_cd4_vois_non$phylotype.x <- factor(comb_cd4_vois_non$phylotype.x, levels = c("Backbone","PT.B.5.UK","PT.B.40.UK","PT.B.69.UK","PT.B.133.UK","PT.B.137.UK"))

cd4_pt_pal <- c("PT.B.5.UK"="#0099B4","PT.B.40.UK"="#088c06","PT.B.69.UK"="#ED0000","PT.B.133.UK"="#00468B","PT.B.137.UK"="#580618", "Backbone"="#ADB6B6")

f2_pl <- ggplot(comb_cd4_vois_non, aes(x = years_since_1cd4, y = cd4, group=phylotype.x, color=as.factor(phylotype.x))) + #x=cd4_decimal_date
	geom_smooth(method = "glm", fill=NA) + theme_classic() + #coord_cartesian(ylim=c(0,1000)) +
	scale_x_continuous(limits=c(0,7), breaks=seq(0,6,1)) + 
	coord_cartesian(ylim=c(0,1000)) +
	scale_color_manual(values=cd4_pt_pal, name="Phylotype") + scale_fill_manual(values=cd4_pt_pal, name="Phylotype") +
	ylab(bquote(bold("CD4 count (cells /"~mm^3~")"))) + xlab("Years since first CD4") +
	theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'),
							legend.key.width = unit(0.5, 'cm'), legend.text=element_text(size=10)) + leg +
	geom_hline(yintercept=350, linetype="longdash", color="black", size=0.3) 
f2_pl
### END FIGURE 2B ###

# Standarise all lines starting at cd4=500 for discussion on decline
comb_cd4_vois_non_std <- comb_cd4_vois_non
comb_cd4_vois_non_std$new_cd4_std <- comb_cd4_vois_non_std$slope * comb_cd4_vois_non_std$years_since_1cd4 + 500
comb_cd4_vois_non_std$time_to_350 <- -150 / comb_cd4_vois_non_std$slope
mean(sort(unique(mean(comb_cd4_vois_non_std$time_to_350[comb_cd4_vois_non_std$phylotype.x == "Backbone"]))))
#sd(sort(unique(comb_cd4_vois_non_std$time_to_350[comb_cd4_vois_non_std$phylotype.x == "Backbone"])), na.rm=TRUE)

ggplot(comb_cd4_vois_non_std, aes(x = years_since_1cd4, y = new_cd4_std, group=phylotype.x, color=as.factor(phylotype.x))) + #x=cd4_decimal_date
	geom_smooth(method = "glm", fill=NA) + theme_classic() + coord_cartesian(ylim=c(0,1000)) +
	scale_x_continuous(limits=c(0,7), breaks=seq(0,7,1)) + 
	coord_cartesian(ylim=c(0,750)) +
	scale_color_manual(values=cd4_pt_pal, name="Phylotype") + scale_fill_manual(values=cd4_pt_pal, name="Phylotype") +
	ylab(bquote(bold("CD4 count (cells /"~mm^3~")"))) + xlab("Years since first CD4") +
	theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'),
							legend.key.width = unit(0.5, 'cm'), legend.text=element_text(size=10)) + leg +
	geom_hline(yintercept=350, linetype="longdash", color="black", size=0.3) 

sort(unique(mean(comb_cd4_vois_non_std$time_to_350[comb_cd4_vois_non_std$phylotype.x == "PT.B.137.UK"])))
sort(unique(mean(comb_cd4_vois_non_std$time_to_350[comb_cd4_vois_non_std$phylotype.x == "PT.B.133.UK"])))
sort(unique(mean(comb_cd4_vois_non_std$time_to_350[comb_cd4_vois_non_std$phylotype.x == "PT.B.69.UK"])))
sort(unique(mean(comb_cd4_vois_non_std$time_to_350[comb_cd4_vois_non_std$phylotype.x == "PT.B.40.UK"])))
sort(unique(mean(comb_cd4_vois_non_std$time_to_350[comb_cd4_vois_non_std$phylotype.x == "PT.B.5.UK"])))

# std all patients from backbone
# get all backbone patients for median time to AIDS calculation

#d_cd4_ages_sl_backbone <- d_cd4_ages_sl[d_cd4_ages_sl$phylotype == 153] #d_cd4_ages_sl
#d_cd4_ages_sl_backbone_slopes <- d_cd4_ages_sl_backbone[,list(intercept=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[2], phylotype="153"), by=patientindex]

d_cd4_ages_sl_backbone_slopes_ij <- inner_join(d_ages_backbone, d_cd4_ages, by="patientindex")
d_cd4_ages_sl_backbone_slopes_std <- d_cd4_ages_sl_backbone_slopes_ij
d_cd4_ages_sl_backbone_slopes_std$new_cd4_std <- d_cd4_ages_sl_backbone_slopes_std$slope * d_cd4_ages_sl_backbone_slopes_std$years_since_1cd4 + 500
d_cd4_ages_sl_backbone_slopes_std$time_to_350 <- -150 / d_cd4_ages_sl_backbone_slopes_std$slope
d_cd4_ages_sl_backbone_slopes_std_non_neg_up_20 <- d_cd4_ages_sl_backbone_slopes_std[d_cd4_ages_sl_backbone_slopes_std$time_to_350 >= -5 & d_cd4_ages_sl_backbone_slopes_std$time_to_350 <= 20,]
d_cd4_ages_sl_backbone_slopes_std_non_neg_up_20 <- d_cd4_ages_sl_backbone_slopes_std_non_neg_up_20 %>% distinct(patientindex, .keep_all = T)
hist(d_cd4_ages_sl_backbone_slopes_std_non_neg_up_20$time_to_350, breaks=1000)
median(d_cd4_ages_sl_backbone_slopes_std_non_neg_up_20$time_to_350) #1.4
mean(d_cd4_ages_sl_backbone_slopes_std_non_neg_up_20$time_to_350) #2.2
IQR(d_cd4_ages_sl_backbone_slopes_std_non_neg_up_20$time_to_350) #2.9
sd(d_cd4_ages_sl_backbone_slopes_std_non_neg_up_20$time_to_350) #3.36

cd4_measur_palette <- c("1"="#00468BB2", "2"="#ED0000B2", "3"="#42B540B2", "4"="#925E9FB2", "5"="#AD002AB2")

### BEGIN FIGURE S5 ###
# Stacked bar of nth measurements per year considering all phylotypes 
d_cd4$year <- lubridate::year(as.Date(date_decimal(d_cd4$cd4_decimal_date)))
d_cd4$year <- as.character(d_cd4$year)
d2 <- d_cd4 %>% group_by(cd4_episode, year) %>% summarise(episode_year_n = n()) # # measurement for each patient and group by it
d2$cd4_episode <- as.factor(d2$cd4_episode)
# Plot distribution of CD4 measurements order over time
ggplot(d2, aes(fill=cd4_episode, y=year, x=episode_year_n)) + #scale_fill_viridis(discrete=T, name="CD4 measurement order") +
	scale_fill_manual(values=cd4_measur_palette, name="CD4 measurement episode") +
	geom_bar(position = position_stack(reverse = TRUE), stat="identity") + labs(x="CD4 measurements", y="Year") + theme_classic() + #position="stack"
	theme(axis.text=element_text(colour="black"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'),legend.title = element_text(size=10), legend.position="top") +
	leg + scale_x_continuous(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) 
ggsave(file=glue("{RESULTS_PATH}/figs/figS5.svg"), dpi=600, width=8, height=6, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS5.jpg"), dpi=600, width=8, height=6, bg="white")
### END FIGURE S5 ###

# Same as above but 1 plot per phylotype
d3 <- d_cd4 %>% group_by(cd4_episode, year, phylotype) %>% summarise(episode_year_n = n())
d3$cd4_episode <- as.factor(d3$cd4_episode)

for(i in 1:154) {
	if(i == 6 || i == 153) next
	d3 %>% filter(phylotype==i) %>% ggplot(aes(fill=cd4_episode, y=year, x=episode_year_n)) + scale_fill_manual(values=cd4_measur_palette, name="CD4 measurement episode") +
		geom_bar(position = position_stack(reverse = TRUE), stat="identity") + labs(x="CD4 measurements", y="Year") + theme_bw() + 
		theme(axis.text=element_text(colour="black"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'),legend.title = element_text(size=10), legend.position="top")
	system(glue("mkdir -p {RESULTS_PATH}/08_cd4regr_plots/30-B/stackbar_phylotypes/"))
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/08_cd4regr_plots/30-B/stackbar_phylotypes/phylotype_{i}_measur_years.png"), dpi=600, width=8, height=6) )
}

### BEGIN FIGURE S6 ###
s3_pl <- list()
c <- 1
for(i in spvl_cd4_vois) {
	s3_pl[[c]] <- d3 %>% filter(phylotype==i) %>% ggplot(aes(fill=cd4_episode, y=year, x=episode_year_n)) +
		geom_bar(position = position_stack(reverse = TRUE), stat="identity") + theme_classic() +
		scale_fill_manual(values=cd4_measur_palette, name="CD4 measurement\nepisode") + 
		#scale_x_continuous(expand = c(0,0), limits=c(0, 70), breaks=c(5,10,20,30,40,50,60,70)) + 
		scale_x_continuous(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + 
		theme(legend.position="top", legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'), legend.text=element_text(size=10), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), axis.title.x=element_blank(), axis.title.y=element_blank()) + leg +
		guides(fill=guide_legend(ncol=1))
	#if(c==1 || c==2) s3_pl[[c]] <- s3_pl[[c]] + ax_rem_a
	c = c+1
}

s3_pl[[2]] <- s3_pl[[2]] + scale_x_continuous(expand = c(0,0), limits=c(0, 20), breaks=c(0,5,10,15,20))
s3_pl[[3]] <- s3_pl[[3]] + scale_x_continuous(expand = c(0,0), limits=c(0, 15), breaks=c(0,5,10,15))
s3_pl[[4]] <- s3_pl[[4]] + scale_x_continuous(expand = c(0,0), limits=c(0, 15), breaks=c(0,5,10,15))
s3_pl[[5]] <- s3_pl[[5]] + scale_x_continuous(expand = c(0,0), limits=c(0, 15), breaks=c(0,5,10,15))

lg2 <- as_ggplot(get_legend(s3_pl[[1]]))
for(i in 1:length(s3_pl)) { s3_pl[[i]] <- s3_pl[[i]] + theme(legend.position = "none") }
s3 <- ggarrange(s3_pl[[1]], s3_pl[[2]], s3_pl[[3]], s3_pl[[4]], s3_pl[[5]], lg2, nrow=2, ncol=3, labels = c("A","B","C","D","E"), font.label=list(family="Times New Roman", color="black",size=10)) #common.legend = T
annotate_figure(s3, left = text_grob("Year", rot = 90, vjust = 1, size=10, family = tnr, face="bold"), bottom = text_grob("CD4 measurements", size=10, family = tnr, face="bold"))
ggsave(file=glue("{RESULTS_PATH}/figs/figS6.svg"), dpi=600, width=8, height=6, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS6.jpg"), dpi=600, width=8, height=6, bg="white")
### END FIGURE S6 ###
