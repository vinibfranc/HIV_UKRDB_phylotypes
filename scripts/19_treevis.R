libs_load <- c("ggtree","ggplot2","dplyr","ggpubr","ape","ggnewscale","glue","lubridate", "treeio","forcats")
invisible( lapply(libs_load, library, character.only=TRUE) )

RDS_PATH <- "rds"

extracted_clusters_B30 <- readRDS(glue("{RDS_PATH}/extracted_clusters.rds"))[[1,4]]
cd4_regr_subtypes_B30 <- readRDS(glue("{RDS_PATH}/cd4_regr_subtypes.rds"))[[1,4]][[2]]
vl_annot_B30 <- readRDS(glue("{RDS_PATH}/vl_subtypes_comb2.rds"))[[1,4]] #readRDS(glue("{RDS_PATH}/vl_subtypes.rds"))[[4]][[1]]

legend_config <- theme(axis.text.x = element_text(color="black", size=5),
																							axis.title.x = element_text(color="black", size=8), #face="bold"
																							#legend.text = element_text(size=12),
																							legend.title = element_text(size=10),
																							plot.margin = unit(c(2,2,0,0), "lines"))

# [[mcs,subtype]][[1]] -> tree
# [[mcs,subtype]][[2]] -> sample_times

# variables to plot: 
# 1) age_group, sex, ethnicity, PHE_regiondiagnosed, years_since_diagnosis
# 2) mean vl per patient, mean CD4 slopes per patient
# 3) DRMs?

# Examples:
# https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html (4.15, 4.17, 4.18, 4.19)
# http://yulab-smu.top/treedata-book/chapter10.html#api-to-link-ggtree-and-ggplot2 (10.3, 10.5)

min_cl_size_choices <- c(30, 50, 100)
tree_names <- c("A_A1","CRF_02_AG","C","B")
subtype_choices <- c("A (A1)","CRF 02_AG","C", "B")

tree_annots <- readRDS("rds/demog_md_subtype_match.rds")
tree_annots <- tree_annots[ (tree_annots$status=="Naïve") & (tree_annots$rega3subtype=="B"), ]

cd4_slope_pal <- c("a. <-5"="#003f5c","b. -5 to -2.6"="#58508d","c. -2.5 to 0"="#bc5090","d. >0"="#ffa600", "Not known"="white")
vl_pal <- c("a. >5"="#003f5c","b. >4.5-5"="#58508d","c. >4-4.5"="#bc5090","d. 0-4"="#ffa600", "Not known"="white")
reg_pal <- c("South of England"="#56B4E9","Not Known"="#000000","London"="#377EB8","North of England"="#D55E00","Northern Ireland, Scotland, and Wales"="#009E73", "Midlands and East of England"="#E69F00")
age_group_palette <- c("0-29"="#0072B2", "30-39"="#E69F00", "40-49"="#56B4E9", "50-59"="#CC79A7", "60+"="#009E73")
helv <- "Helvetica"
leg <- theme(text=element_text(family="Helvetica"), axis.text=element_text(size=10, color="black"), axis.title=element_text(size=10, color="black", face="bold"))

### BEGIN FIGURE 4A and C ###
tree_annots_match_tips_selected <- function(demog_md_choice, subtype_choice, extr_clusters, phylotypes, vl_patients, cd4_patients) {
	demog_md_choice <- demog_md_choice[demog_md_choice$rega3subtype == subtype_choice,]

	demog_md_choice$testindex <- paste0("t.",demog_md_choice$testindex)

	demog_md_choice$hiv_diag_decimal_date <- decimal_date(as.Date(demog_md_choice$hivpos_ymd))
	demog_md_choice <- demog_md_choice %>% mutate(years_since_diagnosis=round(dbsample_date-hiv_diag_decimal_date), age_diag=round(hiv_diag_decimal_date-dob_y))
	demog_md_choice$years_since_diagnosis[demog_md_choice$years_since_diagnosis < 0] <- 0
	demog_md_choice <- demog_md_choice %>% mutate(
		age_group = dplyr::case_when(age_diag<=29 ~ "0-29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+", TRUE ~ "Not knwon"),
		age_group = factor(age_group,level = c("0-29","30-39","40-49","50-59","60+","Not known")
		)
	)
	demog_md_choice <- demog_md_choice %>% mutate(
		hivpos_group = dplyr::case_when(hivpos_year>=1975 & hivpos_year<1990 ~ "1975-1989", hivpos_year>=1990 & hivpos_year<2005 ~ "1990-2004", hivpos_year>=2005 & hivpos_year<2010 ~ "2005-2009", hivpos_year>=2010 & hivpos_year<2015 ~ "2010-2014", hivpos_year>=2015 & hivpos_year<=2020 ~ "2015-2020", TRUE ~ "Not Known"),
		hivpos_group = factor(hivpos_group,level = c("1975-1989","1990-2004","2005-2009","2010-2014","2015-2020","Not known")
		)
	)
	

	tree_annot_md <- vl_annot_md <- cd4_annot_md <- tree_annot_all <- tree_annot_all2 <- tree_annot_all2_dates <- p10 <- ggtr0 <- list() #tt <- ggtr0 
	ct <- 1
	for(k in phylotypes) {
		tree_annot_md[[ct]] <- demog_md_choice[demog_md_choice$testindex %in% extr_clusters[[1]][[k]]$tip.label,]
		tree_annot_md[[ct]] <- tree_annot_md[[ct]] %>% dplyr::select(testindex, dbsample_date, age_group, sexid, ethnicityid, PHE_regiondiagnosed, years_since_diagnosis, hivpos_group, patientindex)
		vl_pt <- vl_patients[vl_patients$cluster == k,]
		vl_annot_md[[ct]] <- tree_annot_md[[ct]] %>% inner_join(vl_pt, by="patientindex")
		cd4_annot_md[[ct]] <- tree_annot_md[[ct]] %>% inner_join(cd4_patients[[k]], by="patientindex")
		tree_annot_all[[ct]] <- left_join(tree_annot_md[[ct]], vl_annot_md[[ct]], by="patientindex")
		tree_annot_all[[ct]] <- left_join(tree_annot_all[[ct]], cd4_annot_md[[ct]], by="patientindex")
		
		tree_annot_all[[ct]] <- tree_annot_all[[ct]] %>% mutate(
			eth_cat = dplyr::case_when(ethnicityid.x == "White" ~ "White", TRUE ~ "Other ethnicities"),
			eth_cat = factor(eth_cat, level = c("White", "Other ethnicities")))

		tree_annot_all[[ct]] <- tree_annot_all[[ct]] %>% mutate(
			yrs_s_diag_cat = dplyr::case_when(years_since_diagnosis.x>=0 & years_since_diagnosis.x<=0.5 ~ "a. 0-0.5", years_since_diagnosis.x>0.5 & years_since_diagnosis.x<=2 ~ "b. 0.6-2", years_since_diagnosis.x>2 & years_since_diagnosis.x<=5 ~ "c. 2.1-5", years_since_diagnosis.x>5 ~ "d. >5", TRUE ~ NA), #"Not Known"
			yrs_s_diag_cat = factor(yrs_s_diag_cat, level = c("a. 0-0.5", "b. 0.6-2", "c. 2.1-5", "d. >5", "Not known")))

		tree_annot_all[[ct]] <- tree_annot_all[[ct]] %>% mutate(
			log_vl_cat = dplyr::case_when(mean_log_vl_pat>5 ~ "a. >5", mean_log_vl_pat>4.5 & mean_log_vl_pat<=5 ~ "b. >4.5-5", mean_log_vl_pat>4 & mean_log_vl_pat<=4.5 ~ "c. >4-4.5", mean_log_vl_pat>=0 & mean_log_vl_pat<=4 ~ "d. 0-4", TRUE ~ NA), #Not known
			log_vl_cat = factor(log_vl_cat, levels=c("a. >5", "b. >4.5-5", "c. >4-4.5", "d. 0-4", "Not known")))

		tree_annot_all[[ct]] <- tree_annot_all[[ct]] %>% mutate(
			cd4_slope_cat = dplyr::case_when(slope < -5 ~ "a. <-5", slope >= -5 & slope < -2.5 ~ "b. -5 to -2.6", slope >= -2.5 & slope<=0 ~ "c. -2.5 to 0", slope > 0 ~ "d. >0", TRUE ~ NA), #"Not known"
			cd4_slope_cat = factor(cd4_slope_cat, levels=c("a. <-5", "b. -5 to -2.6", "c. -2.5 to 0", "d. >0", "Not known")))
		
		tree_annot_all[[ct]] <- tree_annot_all[[ct]] %>% dplyr::select(testindex.x, patientindex, dbsample_date.x, age_group.x, sexid.x, ethnicityid.x, eth_cat, PHE_regiondiagnosed.x, yrs_s_diag_cat, log_vl_cat, cd4_slope_cat)
		tree_annot_all2[[ct]] <- tree_annot_all[[ct]]

		rownames(tree_annot_all[[ct]]) <- tree_annot_all[[ct]]$testindex
		tree_annot_all[[ct]] <- tree_annot_all[[ct]][ , !names(tree_annot_all[[ct]]) %in%  c("testindex.x","patientindex","sexid.x")] #age_group.x
		rn <- rownames(tree_annot_all[[ct]])

		d1 <- tree_annot_all[[ct]] %>% dplyr::select(eth_cat)
		d2 <- tree_annot_all[[ct]] %>% dplyr::select(PHE_regiondiagnosed.x)
		d3 <- tree_annot_all[[ct]] %>% dplyr::select(yrs_s_diag_cat)
		d4 <- tree_annot_all[[ct]] %>% dplyr::select(log_vl_cat)
		
		d5 <- tree_annot_all[[ct]] %>% dplyr::select(cd4_slope_cat)
		
		tree_annot_all2_dates[[ct]] <- tree_annot_all2[[ct]]
		tree_annot_all2_dates[[ct]]$date <- as.Date(date_decimal(tree_annot_all2_dates[[ct]]$dbsample_date.x))
		tree_annot_all2_dates[[ct]] <- tree_annot_all2_dates[[ct]] %>% dplyr::select(testindex.x, date)

		tr <- extr_clusters[[1]][[k]]
		
		tt <- ggtree(tr, mrsd=max(	tree_annot_all2_dates[[ct]]$date), as.Date=FALSE, color="grey20", size=0.2) +
			theme_tree2() + scale_x_ggtree() + xlab("Year")
		
		ggtr0[[ct]] <- tt %<+% tree_annot_all2[[ct]] +	geom_tippoint(aes(color=age_group.x)) + scale_color_manual(values=age_group_palette, name="Age group")
		
		p3 <- gheatmap(ggtr0[[ct]], d2, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("Region"), colnames_offset_y=1, family=helv, font.size=2.5, offset=0) + 
			scale_fill_manual(values=reg_pal, name = "Region\ndiagnosed", na.value="white", guide=guide_legend(order=2, ncol=2)) + coord_cartesian(clip = "off")
		p4 <- p3 + new_scale_fill()
		
		p7 <- gheatmap(p4, d4, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("vl"), colnames_offset_y=1,  family=helv,font.size=2.5, offset=1.75) +
			scale_fill_manual(values=vl_pal, name = "Viral load", na.value="white", guide=guide_legend(order=3, ncol=2)) + coord_cartesian(clip = "off")
		p8 <- p7 + new_scale_fill()
		p9 <- gheatmap(p8, d5, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("CD4"), colnames_offset_y=1, family=helv, font.size=2.5, offset=3.5) +
			scale_fill_manual(values=cd4_slope_pal, name = "CD4 slope", na.value="white", guide=guide_legend(order=4, ncol=2)) + coord_cartesian(clip = "off")
		p10[[ct]] <- p9 + theme(legend.text=element_text(size=7,family=helv), legend.title=element_text(size=9,family=helv), legend.key.height=unit(.5, "cm"),legend.key.width=unit(.5, "cm"),
																										legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0,0,0,0, unit="cm")) + leg + guides(color=guide_legend(order=1, ncol=2)) 
		
		#axis.text.x=element_text(size=4),axis.title.x = element_text(size=4))
		
		ct <- ct+1
	}

	list(ready_pl=p10, tr=ggtr0, md=tree_annot_all)
}

subtype_b_vois_ids_only <- c(40,69,133)
f4bd <- tree_annots_match_tips_selected(tree_annots, "B", extracted_clusters_B30, subtype_b_vois_ids_only, vl_annot_B30, cd4_regr_subtypes_B30)

f4_trees <- function(tr, md) {
	f4 <- tr
	md_dfs <- md
	md_df1 <- md_dfs %>% dplyr::select(age_group.x)
	md_df2 <- md_dfs %>% dplyr::select(PHE_regiondiagnosed.x)
	md_df4 <- md_dfs %>% dplyr::select(log_vl_cat)
	md_df5 <- md_dfs %>% dplyr::select(cd4_slope_cat)
	f4 <- f4  %<+% md_df1 + geom_tippoint(aes(color=age_group.x)) + scale_color_manual(values=age_group_palette, name="Age group", guide=guide_legend(order=1, ncol=1))
	p31 <- gheatmap(f4, md_df2, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("Region"), colnames_offset_y=0.5, family=helv, font.size=2.5, offset=0) + 
		scale_fill_manual(values=reg_pal, name = "Region diagnosed", na.value="white", guide=guide_legend(order=2, ncol=1)) + coord_cartesian(clip = "off") 
	p41 <- p31 + new_scale_fill()
	p71 <- gheatmap(p41, md_df4, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("VL"), colnames_offset_y=0.5, family=helv,font.size=2.5, offset=1.875) +
		scale_fill_manual(values=vl_pal, name = "Viral load", na.value="white", guide=guide_legend(order=3, ncol=1)) + coord_cartesian(clip = "off") 
	p81 <- p71 + new_scale_fill()
	p91 <- gheatmap(p81, md_df5, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("CD4"), colnames_offset_y=0.5, family=helv, font.size=2.5, offset=3.75) +
		scale_fill_manual(values=cd4_slope_pal, name = "CD4 slope", na.value="white", guide=guide_legend(order=4, ncol=1)) + coord_cartesian(clip = "off") 
	p101 <- p91 + new_scale_fill()
	p111 <- p101 + theme(legend.text=element_text(size=7.5,family=helv), legend.title=element_text(size=9,family=helv), legend.key.height=unit(.35, "cm"),legend.key.width=unit(.35, "cm"),
																						legend.position = "top", legend.direction = "vertical", legend.margin = margin(0,0,0,0, unit="cm")) + leg 
	
	return(p111)
}

# PT69 and PT133 are indices [[2]] and [[3]]
f4b_r <- f4_trees(f4bd$tr[[2]],f4bd$md[[2]])
f4b <- f4b_r + theme(legend.margin = margin(t=0.25, r=0, b=0.25, l=0, "cm"))
f4d_r <- f4_trees(f4bd$tr[[3]],f4bd$md[[3]])
f4d <- f4b_r + theme(legend.margin = margin(t=0.25, r=0, b=0.25, l=0, "cm"))

f4bd_pl <- ggarrange(f4b, f4d, nrow=2, ncol=1, labels=c("B","D"), common.legend=T, font.label=list(family=helv, color="black",size=10)) #common.legend=T
saveRDS(f4bd_pl, glue( "{RDS_PATH}/f4bd_pl.rds" ))
### END FIGURE 4B and 4D ###

fs17_trees <- function(tr, md, reduce_colnames=F) {
	fs17 <- tr
	md_dfs <- md
	md_df0 <- md_dfs %>% dplyr::select(age_group.x)
	md_df1 <- md_dfs %>% dplyr::select(eth_cat)
	md_df2 <- md_dfs %>% dplyr::select(PHE_regiondiagnosed.x)
	md_df3 <- md_dfs %>% dplyr::select(yrs_s_diag_cat)
	md_df4 <- md_dfs %>% dplyr::select(log_vl_cat)
	md_df5 <- md_dfs %>% dplyr::select(cd4_slope_cat)
	
	if(reduce_colnames) font_reduced=2
	else font_reduced=2.5
		
	fs17 <- fs17  %<+% md_df0 + geom_tippoint(aes(color=age_group.x)) + scale_color_manual(values=age_group_palette, name="Age group", guide=guide_legend(order=1, ncol=1))
	p11 <- gheatmap(fs17, md_df1, width = 0.1,colnames=T, colnames_position="top", custom_column_labels=c("Ethnicity"), colnames_offset_y=0.5, family=helv, font.size=font_reduced, offset=0) +
		scale_fill_viridis_d(option = "D", name = "Ethnicity", na.value="white", guide=guide_legend(order=2, ncol=1)) + coord_cartesian(clip = "off")
	p21 <- p11 + new_scale_fill()
	p31 <- gheatmap(p21, md_df2, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("Region"), colnames_offset_y=0.5, family=helv, font.size=font_reduced, offset=1.875) + 
		scale_fill_manual(values=reg_pal, name = "Region\ndiagnosed", na.value="white", guide=guide_legend(order=3, ncol=1)) + coord_cartesian(clip = "off")
	p41 <- p31 + new_scale_fill()
	p51 <- gheatmap(p41, md_df3, width = 0.1, colnames_position="top", custom_column_labels=c("Yrs diag."), colnames_offset_y=0.5, family=helv, font.size=font_reduced, offset=3.75) +
		scale_fill_viridis_d(option = "E", name = "Years\nsince\ndiagnosis", guide=guide_legend(order=4,ncol=1)) + coord_cartesian(clip = "off") 
	p61 <- p51 + new_scale_fill()
	p71 <- gheatmap(p61, md_df4, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("VL"), colnames_offset_y=0.5, family=helv, font.size=font_reduced, offset=5.625) +
		scale_fill_manual(values=vl_pal, name = "Viral\nload", na.value="white", guide=guide_legend(order=5, ncol=1)) + coord_cartesian(clip = "off") 
	p81 <- p71 + new_scale_fill()
	p91 <- gheatmap(p81, md_df5, width = 0.1, colnames=T, colnames_position="top", custom_column_labels=c("CD4"), colnames_offset_y=0.5, family=helv, font.size=font_reduced, offset=7.5) +
		scale_fill_manual(values=cd4_slope_pal, name = "CD4\nslope", na.value="white", guide=guide_legend(order=6, ncol=1)) + coord_cartesian(clip = "off")
	p101 <- p91 + new_scale_fill()
	p111 <- p101 + theme(legend.text=element_text(size=6,family=helv), legend.title=element_text(size=8,family=helv), legend.key.height=unit(.4, "cm"), legend.key.width=unit(.4, "cm"),
																						legend.position = "top", legend.direction = "horizontal", legend.margin = margin(-0.1,0,0,0, unit="cm")) + leg 
	
	return(p111)
}

# Figure S19
fs17_a <- fs17_trees(f4bd$tr[[2]], f4bd$md[[2]], reduce_colnames = T) #PT69
fs17_b <- fs17_trees(f4bd$tr[[3]], f4bd$md[[3]], reduce_colnames = T) #PT133
fs17 <- ggarrange(fs17_a + theme(legend.margin = margin(0.5, 0.1, 0.5, 0.1, "cm"), legend.text=element_text(size=7.5,family=helv)), 
																		fs17_b + theme(legend.margin = margin(1, 0.5, 0.5, 0.5, "cm"),legend.text=element_text(size=7.5,family=helv)), 
																		nrow=1, ncol=2, common.legend = T, labels=c("A","B"), font.label=list(family=helv, color="black",size=10))
ggsave(file=glue("{RESULTS_PATH}/figs/figS19_to_edit.svg"), plot=fs17, dpi=600, width=11.5, height=12.5, bg="white")

# PT40 is index 1
fs19 <- fs17_trees(f4bd$tr[[1]], f4bd$md[[1]], reduce_colnames = F)
# Figure S21
ggsave(file=glue("{RESULTS_PATH}/figs/figS21_to_edit.svg"), plot=fs19, dpi=600, width=8.5, height=6.5, bg="white")
### END FIGURES 4B, 4D, S21 ####