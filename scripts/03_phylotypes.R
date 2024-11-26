libs_load <- c("dplyr", "glue","ape","treedater","lubridate","treestructure","ggtree","ggpubr","ggplot2","viridis","data.table", "fastbaps")
invisible( lapply(libs_load, library, character.only=TRUE) )

### STEP 6: IDENTIFY CLUSTERS IN THE TIMETREE (TREESTRUCTURE) AND GET DF WITH METADATA FOR EACH ONE ###

NCPU <- 4
min_cl_size_choices <- c(30, 50, 100) #250,500
tree_names <- c("A_A1","CRF_02_AG","C","B")

timetree_a1_adj <- readRDS("rds/timetree_a1_adj.rds")
tree_a1_adj <- readRDS("rds/tree_a1_adj.rds")
all(timetree_a1_adj$tip.label == tree_a1_adj$tip.label)
class(timetree_a1_adj) <- "phylo"
timetree_a1_boot <- as.integer(tree_a1_adj$node.label)
timetree_a1_boot[is.na(timetree_a1_boot)] <- 95
hist(timetree_a1_boot)

timetree_crf02ag_adj <- readRDS("rds/timetree_crf02ag_adj.rds")
tree_crf02ag_adj <- readRDS("rds/tree_crf02ag_adj.rds")
all(timetree_crf02ag_adj$tip.label == tree_crf02ag_adj$tip.label)
class(timetree_crf02ag_adj) <- "phylo"
timetree_crf02ag_boot <- as.integer(tree_crf02ag_adj$node.label)
timetree_crf02ag_boot[is.na(timetree_crf02ag_boot)] <- 95
hist(timetree_crf02ag_boot)

timetree_c_adj <- readRDS("rds/timetree_c_adj.rds")
tree_c_adj <- readRDS("rds/tree_c_adj2.rds")
all(timetree_c_adj$tip.label == tree_c_adj$tip.label)
class(timetree_c_adj) <- "phylo"
timetree_c_boot <- as.integer(tree_c_adj$node.label)
timetree_c_boot[is.na(timetree_c_boot)] <- 95
hist(timetree_c_boot)

timetree_b_adj <- readRDS("rds/timetree_b_strict_adj.rds")
tree_b_adj <- readRDS("rds/tree_b_adj2.rds")
all(timetree_b_adj$tip.label == tree_b_adj$tip.label)
class(timetree_b_adj) <- "phylo"
timetree_b_boot <- as.integer(tree_b_adj$node.label)
timetree_b_boot[is.na(timetree_b_boot)] <- 95
hist(timetree_b_boot)

timetrees <- list(timetree_a1_adj, timetree_crf02ag_adj, timetree_c_adj, timetree_b_adj)
timetrees_boot <- list(timetree_a1_boot, timetree_crf02ag_boot, timetree_c_boot, timetree_b_boot)
ml_trees <- list(tree_a1_adj, tree_crf02ag_adj, tree_c_adj, tree_b_adj)

system("mkdir -p results/03_treestructure/struct/")
#p_signif <- seq(from=0.00125, 0.01, by=0.00125)
#n_sims <- seq(from=8000, to=1000, by=-1000)

# Treestructure with or without UB support
vary_treestructure_minCladeSize <- function(timetree, min_cl_size, consider_support=TRUE, supp_thr=90, boots) {
	if(consider_support)
		tree_struct_res <- trestruct(timetree, minCladeSize=min_cl_size, nodeSupportValues=boots, nodeSupportThreshold=supp_thr, level=0.01, ncpu=NCPU) # level=p, nsim=nsims
	else 
		tree_struct_res <- trestruct(timetree, minCladeSize=min_cl_size, level=0.01, ncpu=NCPU)
	structure_df <- as.data.frame(tree_struct_res)
	list(tree_struct_res, structure_df)
}

treestruct_min_cl_size_res_no_sup <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(timetrees))
treestruct_min_cl_size_res_yes_sup <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(timetrees))
bt_vals <- c(60,70,80,90)
for(j in 1:length(min_cl_size_choices)) {
	print(glue("min_cl_size_choices={min_cl_size_choices[j]}"))
	for(k in 1:length(timetrees)) {
		print(glue("tree_names={tree_names[k]}"))
		# for(i in 1:length(bt_vals)) {
		# 	print(glue("bt_vals={bt_vals[i]}"))
		# 	treestruct_min_cl_size_res_yes_sup[[j,k]] <- vary_treestructure_minCladeSize(timetrees[[k]], min_cl_size_choices[j], consider_support=TRUE, supp_thr=bt_vals[i], timetrees_boot[[k]]) #, p_signif[l], n_sims[l]
		# 	pl <- plot(treestruct_min_cl_size_res_yes_sup[[j,k]][[1]])  + ggtree::geom_tippoint()
		# 	system(glue("mkdir -p results/03_treestructure/struct/BT{bt_vals[i]}/"))
		# 	ggsave(glue('results/03_treestructure/struct/BT{bt_vals[i]}/mcs{min_cl_size_choices[j]}_subtype_{tree_names[k]}_SUPPORT{bt_vals[i]}.png'), plot=pl, width=8, height=10) #_{p_signif[l]}_{n_sims[l]}
		# }
		
		treestruct_min_cl_size_res_yes_sup[[j,k]] <- vary_treestructure_minCladeSize(timetrees[[k]], min_cl_size_choices[j], consider_support=TRUE, supp_thr=80, timetrees_boot[[k]])
		pl <- plot(treestruct_min_cl_size_res_yes_sup[[j,k]][[1]])  + ggtree::geom_tippoint()
		system(glue("mkdir -p results/03_treestructure/struct/BT80/"))
		ggsave(glue('results/03_treestructure/struct/BT80/mcs{min_cl_size_choices[j]}_subtype_{tree_names[k]}_SUPPORT80.png'), plot=pl, width=8, height=10)
		
		system(glue("mkdir -p results/03_treestructure/struct/noBT/"))
		treestruct_min_cl_size_res_no_sup[[j,k]] <- vary_treestructure_minCladeSize(timetrees[[k]], min_cl_size_choices[j], consider_support=FALSE, timetrees_boot[[k]])
		pl2 <- plot(treestruct_min_cl_size_res_no_sup[[j,k]][[1]])  + ggtree::geom_tippoint()
		ggsave(glue('results/03_treestructure/struct/noBT/mcs{min_cl_size_choices[j]}_subtype_{tree_names[k]}_NO_SUPPORT.png'), plot=pl2, width=8, height=10)
		# }
	}
}

# Fig S1
treestruct_min_cl_size_res_yes_sup <- readRDS("rds/treestruct_min_cl_size_res_yes_sup.rds")

plot_mcs30_b_trestruct <- treestruct_min_cl_size_res_yes_sup[[1,4]][[1]]
plot_mcs30_b_trestruct$clusterSets_bkp <- plot_mcs30_b_trestruct$clusterSets
plot_mcs30_b_trestruct$clusterSets <- list(NULL)

timetree_b_adj <- readRDS("rds/timetree_b_strict_adj.rds")
class(timetree_b_adj) <- "phylo"
plot_mcs30_b_trestruct$clusterSets[[1]] <- timetree_b_adj$tip.label[!(timetree_b_adj$tip.label %in% plot_mcs30_b_trestruct$clusterSets_bkp[[153]])]
plot_mcs30_b_trestruct$clusterSets[[2]] <- timetree_b_adj$tip.label[(timetree_b_adj$tip.label %in% plot_mcs30_b_trestruct$clusterSets_bkp[[153]])]
plot_mcs30_b_trestruct$data$cluster <- ifelse(plot_mcs30_b_trestruct$data$cluster == 153, yes="Backbone", no="Other phylotypes")
plot_mcs30_b_trestruct$data$partition <- ifelse(plot_mcs30_b_trestruct$data$partition == 21, yes="Backbone", no="Other phylotypes")

names(plot_mcs30_b_trestruct$clusterSets) <- c("Other phylotypes", "Backbone")
backbone_other_pal <- c("Other phylotypes"="#00468B7F", "Backbone"="#ED00007F")
tr_s1 <- plot(plot_mcs30_b_trestruct) + ggtree::geom_tippoint(ggplot2::aes_(color=~cluster, shape=NA), size=1.5) 
tr_s1 + ggtree::scale_color_manual(values=backbone_other_pal, name = "Phylotype") + guides(shape=FALSE) #aes(color=partition, shape=NA), size=1
ggsave(file=glue("{RESULTS_PATH}/figs/figS1.svg"), dpi=600, width=8, height=10, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS1.jpg"), dpi=600, width=8, height=10, bg="white")

# Table S2: demographics after outlier removal filters
demog_seq_filters_list <- list()
for(i in 1:length(tree_names)) {
	demog_seq_filters_list[[i]] <- demog_md_subtype_match[paste0("t.",demog_md_subtype_match$testindex) %in% timetrees[[i]]$tip.label,]
}
demog_seq_filters_df <- rbindlist(demog_seq_filters_list) # 40888 patients

# Gender
demog_md_subtype_match_sex_seqs <- demog_seq_filters_df
demog_md_subtype_match_sex_seqs$subtype2 <- ifelse(demog_md_subtype_match_sex_seqs$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_sex_seqs$rega3subtype), "Others")
demog_md_subtype_match_sex_seqs <- demog_md_subtype_match_sex_seqs %>% group_by(sexid, subtype2) %>% summarise(n=n())
View(demog_md_subtype_match_sex_seqs %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})")))

# Ethnicity
demog_md_subtype_match_eth_seqs <- demog_seq_filters_df
demog_md_subtype_match_eth_seqs$subtype2 <- ifelse(demog_md_subtype_match_eth_seqs$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_eth_seqs$rega3subtype), "Others")
demog_md_subtype_match_eth_seqs <- demog_md_subtype_match_eth_seqs %>% mutate(ethnicityid2 = case_when(
	(ethnicityid=="Black-Caribbean") | (ethnicityid=="Black-African") | (ethnicityid=="Black-other/unspecified") ~ "Black-Caribbean / African / other",
	(ethnicityid=="Indian/Pakistani/Bangladeshi") ~ "Indian/Pakistani/Bangladeshi", 
	(ethnicityid=="Other Asian/Oriental") ~ "Other Asian/Oriental", 
	(ethnicityid=="White") ~ "White",
	TRUE ~ "Other/mixed/NA"))
demog_md_subtype_match_eth_seqs <- demog_md_subtype_match_eth_seqs %>% group_by(ethnicityid2, subtype2) %>% summarise(n=n())
View(demog_md_subtype_match_eth_seqs %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})")))

demog_md_subtype_match_exposure_seqs <- demog_seq_filters_df
demog_md_subtype_match_exposure_seqs$subtype2 <- ifelse(demog_md_subtype_match_exposure_seqs$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_exposure_seqs$rega3subtype), "Others")
demog_md_subtype_match_exposure_seqs <- demog_md_subtype_match_exposure_seqs %>% mutate(exposureid2 = case_when(
	(exposureid=="IDU") | (exposureid=="Blood products") ~ "IDU and blood products",
	(exposureid=="Homo/bisexual") ~ "Homo/bisexual", (exposureid=="Heterosexual") ~ "Heterosexual", TRUE ~ "Other/NA")) #(exposureid=="Not known") ~ "Not known"
demog_md_subtype_match_exposure_seqs <- demog_md_subtype_match_exposure_seqs %>% group_by(exposureid2, subtype2) %>% summarise(n=n())
View(demog_md_subtype_match_exposure_seqs %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})")))

# Region
demog_md_subtype_match_region_seqs <- demog_seq_filters_df
demog_md_subtype_match_region_seqs$subtype2 <- ifelse(demog_md_subtype_match_region_seqs$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_region_seqs$rega3subtype), "Others")
demog_md_subtype_match_region_seqs <- demog_md_subtype_match_region_seqs %>% group_by(PHE_regiondiagnosed, subtype2) %>% summarise(n=n())
View(demog_md_subtype_match_region_seqs %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})")))

# Age at diagnosis
demog_md_subtype_match_age_seqs <- demog_seq_filters_df
demog_md_subtype_match_age_seqs$subtype2 <- ifelse(demog_md_subtype_match_age_seqs$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_age_seqs$rega3subtype), "Others")
demog_md_subtype_match_age_seqs$hiv_diag_decimal_date <- decimal_date(as.Date(demog_md_subtype_match_age_seqs$hivpos_ymd))
demog_md_subtype_match_age_seqs <- demog_md_subtype_match_age_seqs %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y))
demog_md_subtype_match_age_seqs <- demog_md_subtype_match_age_seqs %>% mutate(
	age_group = dplyr::case_when(age_diag<=29 ~ "<29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
	age_group = factor(age_group,level = c("<29","30-39","40-49","50-59","60+")))
demog_md_subtype_match_age_seqs <- demog_md_subtype_match_age_seqs %>% group_by(age_group, subtype2) %>% summarise(n=n())
View(demog_md_subtype_match_age_seqs %>% group_by(subtype2) %>% mutate(n_perc=glue("{n} ({round(n*100/sum(n),1)})")))

# get partitions using fastbaps (only one for each subtype)
seqs_folder_naive <- "data/subtype_seqs_naive"
fasta_paths <- c( glue("{seqs_folder_naive}/A_A1_curated_refB_aln_len_filter.fasta"), 
													glue("{seqs_folder_naive}/CRF_02_AG_curated_refB_aln_len_filter.fasta"), 
													glue("{seqs_folder_naive}/C_UK_final_aln_outliers_removed.fasta"), 
													glue("{seqs_folder_naive}/B_UK_final_aln_outliers_removed.fasta"))
fastas_match_tre <- fsparse <- fastbaps_parts <- fastbaps_df_list <- fastbaps_part_plot <- n_tips_each_cluster <- list()
for(i in 1:length(fasta_paths)) {
	fastas_match_tre[[i]] <- read.dna(fasta_paths[i], format="fasta")
	if(i == 1 | i == 2) rownames(fastas_match_tre[[i]]) <- paste0("t.",rownames(fastas_match_tre[[i]]))
	print(length(rownames(fastas_match_tre[[i]])))
	fastas_match_tre[[i]] <- fastas_match_tre[[i]][rownames(fastas_match_tre[[i]]) %in% ml_trees[[i]]$tip.label,]
	print(length(rownames(fastas_match_tre[[i]])))
	write.dna(fastas_match_tre[[i]], glue("{seqs_folder_naive}/{tree_names[i]}_match_tre.fasta"), format="fasta")
	
	# import sparse matrix
	f_aux <- read.FASTA(glue("{seqs_folder_naive}/{tree_names[i]}_match_tre.fasta"))
	fsparse[[i]] <- import_fasta_sparse_nt(f_aux, prior="baps")
	fsparse[[i]] <- optimise_prior(fsparse[[i]], type = "optimise.symmetric") #optimise.symmetric, optimise.baps, baps, symmetric
	# find partition
	fastbaps_parts[[i]] <- best_baps_partition(fsparse[[i]], ml_trees[[i]])
	
	fastbaps_df_list[[i]] <- data.frame(id = ml_trees[[i]]$tip.label, cluster = fastbaps_parts[[i]], stringsAsFactors = FALSE)
	gg <- ggtree(ml_trees[[i]])
	fastbaps_part_plot[[i]] <- facet_plot(gg, panel = "fastbaps", data = fastbaps_df_list[[i]], geom = geom_tile, aes(x = cluster), 
																		color = "blue")
	system("mkdir -p results/03_fastbaps/")
	ggsave(plot=fastbaps_part_plot[[i]], filename=glue("results/03_fastbaps/{tree_names[[i]]}_optim_symm.png"), width=10, height=8, dpi=300)
	
	n_tips_each_cluster[[i]] <- fastbaps_df_list[[i]] %>% group_by(cluster) %>% summarise(n=n())
	summary(n_tips_each_cluster[[i]])
	
	print(min(n_tips_each_cluster[[i]]$n))
	print(max(n_tips_each_cluster[[i]]$n))
	print(length(unique(fastbaps_df_list[[i]]$cluster)))
	print("=====")
	pdf(file=glue("results/03_fastbaps/{tree_names[[i]]}_hist_nclusters_optim_symm.pdf"))
	hist(n_tips_each_cluster[[i]]$n, breaks = 200)
	dev.off()
}

RDS_PATH <- "rds"
# TODO debug why fastbaps giving only one partition for CRF and C?
saveRDS(fastbaps_df_list, glue("{RDS_PATH}/fastbaps_df_list.rds")) # A1 26 partitions (17 with vl), CRF 1, C 1, B 126 (119 with VL)
