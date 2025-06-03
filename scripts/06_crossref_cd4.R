libs_load <- c("ggplot2","dplyr","ggpubr","glue","data.table","lubridate", "mgcv", "viridis")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"

cd4_md <- read.csv("data/cd4s.csv", header=T) #168,861 cd4s, 74,720 unique patients
print(nrow(cd4_md)); length(unique(cd4_md$patientindex))
cd4_md$cd4_date_ymd <- as.Date(gsub("\\/", "15", cd4_md$cd4_date_my), "%m%d%Y")

# Remove patients with 0 values
cd4_md <- cd4_md[cd4_md$cd4 > 0,]
print(nrow(cd4_md)) #168660
# sqrt transform cd4 values
cd4_md$sqrt_cd4 <- sqrt(cd4_md$cd4)
cd4_md$cd4_decimal_date <- decimal_date(as.Date(cd4_md$cd4_date_ymd))
cd4_md <- cd4_md[!is.na(cd4_md$cd4_date_my) & !is.na(cd4_md$cd4),]
print(nrow(cd4_md)) #168660
# check if CD4s before diagnosis
saveRDS(cd4_md, file=glue("{RDS_PATH}/cd4_md.rds"))

# these plots include all risk groups and the four considered subtypes
plot_cd4_measurement_distributions <- function(cd4_df, out_prefix) {
	# get only ids where more than one CD4 measurement was performed
	out_f <- "06_cd4_distr"
	cd4_md <- cd4_df[cd4_df$patientindex %in% unique(cd4_df$patientindex[duplicated(cd4_df$patientindex)]),]
	cd4_md_hist <- cd4_md %>% group_by(patientindex) %>% summarise(n=n())
	#View(cd4_md_hist)
	if(!dir.exists(glue("{RESULTS_PATH}/{out_f}/"))) dir.create(glue("{RESULTS_PATH}/{out_f}/"))
	ggplot(cd4_md_hist, aes(n)) + geom_bar(fill="#69b3a2") + labs(x="CD4 measurements", y="Count") #theme_classic()
	ggsave(file=glue("{RESULTS_PATH}/{out_f}/{out_prefix}_cd4_hist_measur.pdf"), dpi=600, width=6, height=5)
	
	cd4_md$year <- format(as.Date(cd4_md$cd4_date_ymd, format="%Y-%m-%d"),"%Y")
	cd4_md_id_occur_test <- cd4_md %>% mutate(measurement = rowid(patientindex))
	cd4_md_id_occur <- cd4_md %>% mutate(measurement = rowid(patientindex)) %>% group_by(measurement, year) %>% summarise(measur_year_n = n()) # # measurement for each patient and group by it
	cd4_md_2_more_cd4_year <- cd4_md %>% group_by(year) %>% summarise(year_n = n())
	cd4_md_measur_time <- merge(cd4_md_id_occur, cd4_md_2_more_cd4_year, by="year")
	cd4_md_measur_time <- cd4_md_measur_time %>% mutate(measurement_new = case_when(measurement==1 ~'1st', measurement==2 ~'2nd', measurement==3 ~'3rd', measurement==4 ~'4th', measurement==5 ~'5th', TRUE ~ 'Others'))
	cd4_md_measur_time$measurement_new <- as.factor(cd4_md_measur_time$measurement_new)
	# Plot distribution of CD4 measurements order over time
	ggplot(cd4_md_measur_time, aes(fill=measurement_new, y=year, x=measur_year_n)) + scale_fill_viridis(discrete=T, name="CD4 measurement order") +
		geom_bar(position="stack", stat="identity") + labs(x="CD4 measurements", y="Year") + #theme_classic() +
		theme(axis.text=element_text(colour="black"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'),legend.title = element_text(size=10), legend.position="top")
	ggsave(file=glue("{RESULTS_PATH}/{out_f}/{out_prefix}_cd4_measur_order_time.pdf"), dpi=600, width=6, height=5)
	
	# Plot raw number of CD4 measurements per year
	ggplot(cd4_md_2_more_cd4_year, aes(y=year, x=year_n)) + scale_fill_viridis(discrete=T, name="CD4 measurements per year") +
		geom_bar(stat="identity", fill="#69b3a2") + labs(x="CD4 measurements", y="Year") + #theme_classic() +
		theme(axis.text=element_text(colour="black"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.5, 'cm'),legend.title = element_text(size=10), legend.position="top")
	ggsave(file=glue("{RESULTS_PATH}/{out_f}/{out_prefix}_cd4_measur_raw_time.pdf"), dpi=600, width=6, height=5)
	
	write.csv(cd4_md_measur_time, file=glue("{RESULTS_PATH}/{out_f}/{out_prefix}_cd4_measur_order_time.csv"), quote=F, row.names=F) #col.names=T,
	write.csv(cd4_md_2_more_cd4_year, file=glue("{RESULTS_PATH}/{out_f}/{out_prefix}_cd4_measur_raw_time.csv"), quote=F, row.names=F)
	
	list(cd4_md_id_occur_test, cd4_md_measur_time, cd4_md_2_more_cd4_year)
}

# get only ids where more than one CD4 measurement was performed (118k to 89305 rows --> 38k patients to 16334)
cd4_more_2_measur <- plot_cd4_measurement_distributions(cd4_md, "all")

# from 74615 unique patients (initially) to 37037 with more than 1 cd4 measurement
cd4_md_mult <- cd4_md[cd4_md$patientindex %in% unique(cd4_md$patientindex[duplicated(cd4_md$patientindex)]),]
print(nrow(cd4_md_mult)); print(length(unique(cd4_md_mult$patientindex)))

remove_cd4_after_art <- function(demog_md_choice, subtype_choice, cd4_md_) {
	demog_md_choice <- demog_md_choice[demog_md_choice$rega3subtype == subtype_choice,]
	print("rows subtype md:")
	print(nrow(demog_md_choice))
	demog_md_choice$artstart_decimal_date <- decimal_date( as.Date(gsub("\\/", "15", demog_md_choice$artstart_my), "%m%d%Y") )
	#View(demog_md_choice)
	demog_md_cd4_merged <- demog_md_choice %>% left_join(cd4_md_, by="patientindex", relationship="many-to-many") #multiple="all"
	suppressWarnings( cd4_before_art <- demog_md_cd4_merged %>% filter(cd4_decimal_date <= artstart_decimal_date) )
	
	# for the df above I have for each sample all measurements of CD4
	lookup_pat_tests <- unique( cd4_before_art[ , c('patientindex','testindex') ] )
	lookup_pat_tests$testindex <- paste0("t.",lookup_pat_tests$testindex)
	
	# Here removes duplicates of patient and cd4 date
	cd4_before_art_cl <- cd4_before_art %>% distinct(patientindex, cd4_decimal_date, .keep_all = TRUE)
	cd4_before_art_cl <- cd4_before_art_cl[ , !names(cd4_before_art_cl) %in%  c("testindex","dbsample_date")]
	# keep all
	cd4_before_art_cl2 <- cd4_before_art_cl
	# only ones with >1 meas
	cd4_before_art_cl <- cd4_before_art_cl[cd4_before_art_cl$patientindex %in% unique(cd4_before_art_cl$patientindex[duplicated(cd4_before_art_cl$patientindex)]),]
	#View(cd4_before_art_cl)
	print("nrows cd4 before art")
	print(nrow(cd4_before_art_cl))
	print("unique patients before art")
	print(length(unique(cd4_before_art_cl$patientindex)))
	# Check if there are CD4 before diagnosis (perhaps before HIV infection)
	cd4_before_art_cl$hivpos_decimal_date <- decimal_date(as.Date(cd4_before_art_cl$hivpos_ymd))
	cd4_before_art_cl$diff_diag_cd4_dates <- cd4_before_art_cl$hivpos_decimal_date - cd4_before_art_cl$cd4_decimal_date
	
	# IMPORTANT: Measurements made in the days before diagnosis may still be valid for infection-related CD4 decrease, 
	# so considering up to one month ok (as considering mid-month precision, in the worst scenario would include CD4 from 1.5 month before diag)
	
	cd4_before_art_and_diag <- cd4_before_art_cl[cd4_before_art_cl$diff_diag_cd4_dates > 0.084,] # 1 month = 0.0833334 year, 3 months = 0.25 year, 6 months = 0.5 year
	print("CD4 more than 1 month diagnosis? Unique patients")
	print(length(unique(cd4_before_art_and_diag$patientindex)))
	
	# keeping the ones with NA diff_diag_cd4_dates because unknown hivpos_my (but known artstart_date)
	# but excluding the ones with CD4 measurements more than 1 month before diag
	
	cd4_before_art_cl <- cd4_before_art_cl[is.na(cd4_before_art_cl$diff_diag_cd4_dates) | cd4_before_art_cl$diff_diag_cd4_dates < 0.084,]
	print("unique patients before art removing the ones with KNOWN CD4 more than one month BEFORE diagnosis: ")
	print(length(unique(cd4_before_art_cl$patientindex)))
	print("same as above but now excluding the ones that now have less than 2 measurements (after removing the ones long before diag): ")
	cd4_before_art_cl <- cd4_before_art_cl[cd4_before_art_cl$patientindex %in% unique(cd4_before_art_cl$patientindex[duplicated(cd4_before_art_cl$patientindex)]),]
	print(length(unique(cd4_before_art_cl$patientindex)))
	#View(cd4_before_art_and_diag)
	
	print("unique patients lookup_pat_tests")
	lookup_pat_tests <- lookup_pat_tests[lookup_pat_tests$patientindex %in% cd4_before_art_cl$patientindex,]
	print(length(unique(lookup_pat_tests$patientindex)))
	
	list(cd4_before_art_cl, lookup_pat_tests, cd4_before_art_cl2)
}

demog_md_subtype_match <- readRDS(glue("{RDS_PATH}/demog_md_subtype_match.rds")) #100591
demog_md_subtype_match <- demog_md_subtype_match[ (demog_md_subtype_match$status=="Naïve") & 
																																																			(demog_md_subtype_match$rega3subtype %in% c("A (A1)","CRF 02_AG","C","B")) & 
																																																			(demog_md_subtype_match$exposureid != "Not known"), ] #53382

subtype_choices <- c("A (A1)","CRF 02_AG","C", "B")
cd4_before_art_subtypes <- list()
for(i in 1:length(subtype_choices)) {
	print(subtype_choices[i])
	cd4_before_art_subtypes[[i]] <- remove_cd4_after_art(demog_md_subtype_match, subtype_choices[i], cd4_md_mult)
}
# CD4 before diagnosis (perhaps before HIV infection) -> A1: 6 patients, CRF: 12, C: 35; B: 122
# if flagging only the ones more than 3 months before diagnosis -> A1: 4 patients, CRF: 9, C: 21; B: 90
# if flagging only the ones more than 3 months before diagnosis -> A1: 3 patients, CRF: 7, C: 14, B: 71
# if flagging only the ones more than 1 month before diagnosis -> A1: 5, CRF: 10, C: 32, B: 109
# patients before removal and after -> A1: 864 unchanged, CRF: 774 unchanged, C: 2520 to 2518 (2), B: 11011 to 11109 (2)
# same as above but RM the ones that now have less than 2 measurements: 864 to 863 (1), 774 to 768 (6), 2520 to 2505 (15), 11011 to 10971 (40)
saveRDS(cd4_before_art_subtypes, glue("{RDS_PATH}/cd4_before_art_subtypes.rds"))

give.n <- function(x){
	return(c(y = median(x)*1.05, label = length(x))) # experiment with the multiplier to find the perfect position
}

# initial exploratory analysis without formal modelling (comparying phylotype against backbone)
calc_regr_slopes <- function(cd4_before_art_subt, subtype_choice, extr_clusters, lookup_pat_tests, backbone_control) { #extr_clusters
	# Calculate individual regression for each test and get intercept, slope, r2, and median_time_measurement (FOR ALL PATIENTS WITH CD4 INFO WITHOUT MATCHING WITH SAMPLES)
	cd4_md_all <- as.data.table(cd4_before_art_subt)
	# Simple calculations not considering random effect of repeated measurements for each patient
	
	# Using sqrt CD4
	cd4_md_all_sqrtcd4 <- cd4_md_all[,list(intercept=coef(lm(sqrt_cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(sqrt_cd4~cd4_decimal_date, na.action=na.exclude))[2], med_time_cd4_measured=median(cd4_decimal_date)), by=patientindex]
	
	# Using raw CD4
	cd4_md_all_rawcd4 <- cd4_md_all[,list(intercept=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[2], med_time_cd4_measured=median(cd4_decimal_date)), by=patientindex]
	
	cd4_md_all <- as.data.frame(cd4_md_all)
	
	tree_clusters_tips_cd4_info <- md_clusters_tips_cd4_info <- regr_cd4_raw <- regr_cd4_sqrt <- list()
	for(k in 1:length(extr_clusters[[1]])) {
		print(paste("k=",k))
		tree_clusters_tips_cd4_info[[k]] <- data.frame(testindex=extr_clusters[[1]][[k]]$tip.label)
		print("ntips")
		print(nrow(tree_clusters_tips_cd4_info[[k]]))
		
		if(nrow(tree_clusters_tips_cd4_info[[k]]) != 0) {
			tree_clusters_tips_cd4_info[[k]] <- tree_clusters_tips_cd4_info[[k]] %>% inner_join(lookup_pat_tests, by="testindex")
			print("unique patients before cd4 matching")
			print(length(unique(tree_clusters_tips_cd4_info[[k]]$patientindex)))
			
			md_clusters_tips_cd4_info[[k]] <- tree_clusters_tips_cd4_info[[k]] %>% inner_join(cd4_before_art_subt, by="patientindex", multiple="all") #right_join
			# print("rows md"); print(nrow(md_clusters_tips_cd4_info[[k]]))
			print("unique patients after matching")
			print(length(unique(md_clusters_tips_cd4_info[[k]]$patientindex)))
			
			if(nrow(md_clusters_tips_cd4_info[[k]]) > 0) {
				md_clusters_tips_cd4_info[[k]] <- as.data.table(md_clusters_tips_cd4_info[[k]])
				# regression for each individual patient inside phylotype and get dt object with intercept and slope
				#View(md_clusters_tips_cd4_info[[k]])
				regr_cd4_sqrt[[k]] <- md_clusters_tips_cd4_info[[k]][,list(intercept=coef(lm(sqrt_cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(sqrt_cd4~cd4_decimal_date, na.action=na.exclude))[2], med_time_cd4_measured=median(cd4_decimal_date), phylotype=k), by=patientindex] #r2=summary(lm(sqrt_cd4~cd4_decimal_date, na.action=na.exclude))$adj.r.squared,
				regr_cd4_raw[[k]] <- md_clusters_tips_cd4_info[[k]][,list(intercept=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[1], slope=coef(lm(cd4~cd4_decimal_date, na.action=na.exclude))[2], med_time_cd4_measured=median(cd4_decimal_date), phylotype=k), by=patientindex]
				md_clusters_tips_cd4_info[[k]] <- as.data.frame(md_clusters_tips_cd4_info[[k]])
				md_clusters_tips_cd4_info[[k]] <- md_clusters_tips_cd4_info[[k]][ , !names(md_clusters_tips_cd4_info[[k]]) %in%  c("testindex")]
				md_clusters_tips_cd4_info[[k]]$phylotype <- k
				
				# Plot individual regressions inside each cluster
				p <- ggplot(md_clusters_tips_cd4_info[[k]], aes(x = cd4_decimal_date, y = sqrt_cd4, group=patientindex, color=patientindex)) + #color=patientindex
					geom_smooth(method = "lm", fill = NA) + theme_classic() + theme(legend.position="none") + ylim(c(0,50)) +
					ylab("Squared CD4 count")+ xlab("Year") + theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"))
				system(glue("mkdir -p {RESULTS_PATH}/06_regr/{subtype_choice}/"))
				suppressMessages( ggsave(file=glue("{RESULTS_PATH}/06_regr/{subtype_choice}/cluster_{k}_sqrt.png"), plot=p, dpi=300, width=3.85, height=2.8) )
				
				p2 <- ggplot(md_clusters_tips_cd4_info[[k]], aes(x = cd4_decimal_date, y = cd4, group=patientindex, color=patientindex)) + #color=patientindex
					geom_smooth(method = "lm", fill = NA) + theme_classic() + theme(legend.position="none") + ylim(c(0,2500)) +
					ylab("Raw CD4 count")+ xlab("Year") + theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"))
				system(glue("mkdir -p {RESULTS_PATH}/06_regr/{subtype_choice}/"))
				suppressMessages( ggsave(file=glue("{RESULTS_PATH}/06_regr/{subtype_choice}/cluster_{k}_raw.png"), plot=p2, dpi=300, width=3.85, height=2.8) )
			}
		}
		
	}
	
	regr_cd4_all <- rbindlist(regr_cd4_sqrt)
	
	# Add backbone/control cluster
	cd4_md_all_sqrtcd4$phylotype <- as.integer(backbone_control)
	cd4_md_all_sqrtcd4 <- cd4_md_all_sqrtcd4 %>% inner_join(lookup_pat_tests, by="patientindex", multiple="all")
	cd4_md_all_sqrtcd4 <- cd4_md_all_sqrtcd4[cd4_md_all_sqrtcd4$testindex %in% extr_clusters[[1]][[ as.integer(backbone_control) ]]$tip.label ,] #cd4_md_all_sqrtcd4[cd4_md_all_sqrtcd4$phylotype == backbone_control,]
	#print(colnames(cd4_md_all_sqrtcd4))
	#print(colnames(regr_cd4_all))
	regr_cd4_all <- rbind(cd4_md_all_sqrtcd4, regr_cd4_all, fill=TRUE)
	
	regr_cd4_all_count <- regr_cd4_all %>% group_by(phylotype) %>% summarise(count=n())
	
	regr_cd4_all <- regr_cd4_all[!is.na(regr_cd4_all$slope),]
	regr_cd4_all$phylotype <- as.factor(regr_cd4_all$phylotype)
	regr_cd4_all <- as.data.frame(regr_cd4_all)
	
	system(glue("mkdir -p {RESULTS_PATH}/06_boxplots_slopes/"))
	bxp <- ggboxplot(regr_cd4_all, x="phylotype", y="slope", fill="phylotype", outlier.shape=NA) + coord_cartesian(ylim=c(-20,20)) + labs(x="Phylotype", y="Slope of CD4+ cell counts") +
		stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
		theme_classic() + theme(legend.position="none")
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/06_boxplots_slopes/{subtype_choice}_phylotypes.png"), plot=bxp, dpi=300, width=14, height=7, limitsize=FALSE) )
	
	cd4s_across_clusters <- bind_rows(md_clusters_tips_cd4_info)
	uq_backbone <- unique(regr_cd4_all$patientindex[regr_cd4_all$phylotype == as.integer(backbone_control)])
	cd4_md_control <- cd4_md_all[cd4_md_all$patientindex %in% uq_backbone, ]
	cd4_md_control$phylotype <- "ref"
	
	cd4s_across_clusters$testindex <- NULL
	# print(colnames(cd4_md_control)); print(colnames(cd4s_across_clusters))
	cd4s_across_clusters_all <- rbind(cd4_md_control, cd4s_across_clusters)
	cd4s_across_clusters_all$phylotype <- as.factor(cd4s_across_clusters_all$phylotype)
	
	# Add years since diag and age categories
	cd4s_across_clusters_all$hiv_diag_decimal_date <- decimal_date(as.Date(cd4s_across_clusters_all$hivpos_ymd))
	# cd4 episodes
	cd4s_across_clusters_all <- cd4s_across_clusters_all %>% group_by(patientindex) %>% mutate(cd4_episode = seq_along(patientindex)) #mutate(cd4_episode = rowid(patientindex))
	# years from current cd4 to 1st cd4 episode (measurement)
	cd4s_across_clusters_all <- cd4s_across_clusters_all %>% group_by(patientindex) %>% mutate(years_since_1cd4 = round(cd4_decimal_date - cd4_decimal_date[cd4_episode==1], 3) )
	cd4s_across_clusters_all <- cd4s_across_clusters_all %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y)) #years_since_diagnosis=round(cd4_decimal_date-hiv_diag_decimal_date),
	# Age categories: https://www.science.org/doi/10.1126/science.abk1688
	cd4s_across_clusters_all <- cd4s_across_clusters_all %>% mutate(
		age_group = dplyr::case_when(age_diag<=29 ~ "0-29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
		age_group = factor(age_group,level = c("0-29","30-39","40-49","50-59","60+")
		)
	)
	
	cd4s_seqs_cluster <- cd4s_ref_ww_seqs_cl <- cd4s_both <- list()
	for(k in 1:length(extr_clusters[[1]])) {
		cd4s_seqs_cluster[[k]] <- cd4s_across_clusters_all[cd4s_across_clusters_all$phylotype == k,]
		# adjust indices of cd4 episodes
		cd4s_seqs_cluster[[k]] <- cd4s_seqs_cluster[[k]] %>% group_by(patientindex) %>% mutate(cd4_episode = seq_along(patientindex))
		cd4s_seqs_cluster[[k]] <- cd4s_seqs_cluster[[k]] %>% group_by(patientindex) %>% mutate(years_since_1cd4 = round(cd4_decimal_date - cd4_decimal_date[cd4_episode==1], 3) )
		
		cd4s_ref_ww_seqs_cl[[k]] <- cd4s_across_clusters_all[ (cd4s_across_clusters_all$phylotype == "ref") & !(cd4s_across_clusters_all$patientindex %in% cd4s_seqs_cluster[[k]]$patientindex) ,]
		
		cd4s_both[[k]] <- list(cd4s_seqs_cluster[[k]], cd4s_ref_ww_seqs_cl[[k]])
	}
	
	list(cd4s_both, regr_cd4_sqrt)
}
#View(tree_clusters_tips_cd4_info[[1]])

min_cl_size_choices <- c(30, 50, 100) #250,500
tree_names <- c("A_A1","CRF_02_AG","C","B")
extracted_clusters <- readRDS(glue("{RDS_PATH}/extracted_clusters.rds"))
backbone_cl_control <- readRDS(glue("{RDS_PATH}/backbone_cl_control.rds"))

cd4_regr_subtypes <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		print(glue("{min_cl_size_choices[i]}-{tree_names[j]}"))
		cd4_regr_subtypes[[i,j]] <- calc_regr_slopes(cd4_before_art_subtypes[[j]][[1]], glue("{min_cl_size_choices[i]}-{tree_names[j]}"), extracted_clusters[[i,j]], cd4_before_art_subtypes[[j]][[2]], backbone_cl_control[[i,j]])
	}
}
saveRDS(cd4_regr_subtypes, glue("{RDS_PATH}/cd4_regr_subtypes.rds"))