libs_load <- c("dplyr", "glue","ape","treedater","lubridate","treestructure","ggtree","ggpubr","ggplot2","viridis","data.table", "fastbaps", "forcats", "pbmcapply")
invisible( lapply(libs_load, library, character.only=TRUE) )

### IDENTIFY CLUSTERS IN THE TIMETREE (TREESTRUCTURE) AND GET DF WITH METADATA FOR EACH ONE ###

NCPU <- 4
min_cl_size_choices <- c(30, 50, 100) #250,500
tree_names <- c("A_A1","CRF_02_AG","C","B")

RESULTS_PATH="results"

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
plot_mcs30_b_trestruct$data$cluster <- ifelse(plot_mcs30_b_trestruct$data$cluster == 153, yes="Backbone", no="Different phylotypes")
plot_mcs30_b_trestruct$data$partition <- ifelse(plot_mcs30_b_trestruct$data$partition == 21, yes="Backbone", no="Different phylotypes")

names(plot_mcs30_b_trestruct$clusterSets) <- c("Different phylotypes", "Backbone")
backbone_other_pal <- c("Different phylotypes"="#56B4E9", "Backbone"="#D55E00")
tr_s1 <- plot(plot_mcs30_b_trestruct) + ggtree::geom_tippoint(ggplot2::aes_(color=~cluster, shape=NA), size=1.5) 
tr_s1 + ggtree::scale_color_manual(values=backbone_other_pal, name = "Phylotype status") + guides(shape=FALSE) #aes(color=partition, shape=NA), size=1
ggsave(file=glue("{RESULTS_PATH}/figs/figS1.eps"), device=cairo_ps, dpi=600, width=8, height=10, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS1.jpg"), dpi=600, width=8, height=10, bg="white")

# Table S2: demographics after sequence outlier removal filters
demog_seq_filters_list <- list()
for(i in 1:length(tree_names)) {
	demog_seq_filters_list[[i]] <- demog_md_subtype_match[paste0("t.",demog_md_subtype_match$testindex) %in% timetrees[[i]]$tip.label,]
}
demog_seq_filters_df <- rbindlist(demog_seq_filters_list) # 40888 patients
saveRDS(demog_seq_filters_df, file=glue("{RDS_PATH}/demog_seq_filters_df.rds"))

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

### FASTBAPS TREE PARTITIONING ###
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
# fastbaps giving only one partition for CRF and C
saveRDS(fastbaps_df_list, glue("{RDS_PATH}/fastbaps_df_list.rds")) # A1 26 partitions (17 with vl), CRF 1, C 1, B 126 (119 with VL)

### CLUSTER PICKER ###
# First of all, plot intra-phylotype genetic diversity (because clusterPicker relies on that and if too high might not be a good choice)
# IMPORTANT: this relies on phylotype fastas computed in step 04_aln_all_mono...

ALN_LEN <- 995

calculate_pi <- function(distance_matrix) {
	diag(distance_matrix) <- NA
	# Flatten the upper triangle of the matrix (excluding NAs)
	pairwise_distances <- distance_matrix[upper.tri(distance_matrix, diag = FALSE)]
	
	pairwise_distances <- na.omit(pairwise_distances)
	
	# Calculate the number of sequences (rows/columns in the matrix)
	n <- nrow(distance_matrix)
	
	# Calculate pi
	pi <- sum(pairwise_distances) / (n * (n - 1) / 2)
	
	return(pi)
}

plot_intraphylotype_genetic_diversity <- function(folder_fastas, color_limits = c(0, 100)) {
	# Get the list of files
	fold <- list.files(path = folder_fastas, full.names = TRUE, pattern = "*.fasta")
	fastas <- lapply(X = fold, FUN = read.FASTA)
	
	res_path <- glue("{RESULTS_PATH}/03_cluster_picker/intra_pt_distances_heatmaps/")
	system(glue("mkdir -p {res_path}"))
	
	# Process each file
	summary_list <- mapply(function(file_path, fasta_data) {
		#print(file_path)
		D_raw <- dist.dna(fasta_data, model = 'raw', as.matrix = TRUE, pairwise.deletion = TRUE) * ALN_LEN
		diag(D_raw) <- NA  # don't count zero distances on diagonal
		
		# Calculate statistics
		mean_value <- mean(D_raw, na.rm = TRUE)
		median_value <- median(D_raw, na.rm = TRUE)
		pi_value <- calculate_pi(D_raw)
		max_value <- max(D_raw, na.rm=TRUE)
		
		# Reshape matrix
		D_raw_mlt <- reshape2::melt(D_raw)
		D_raw_mlt <- D_raw_mlt[order(D_raw_mlt$Var1, D_raw_mlt$Var2), ]
		
		# Create the heatmap plot
		pl <- ggplot(data = D_raw_mlt, aes(x = as.factor(Var1), y = as.factor(Var2), fill = value)) +
			geom_tile() +
			labs(x = "", y = "") +
			theme_classic() +
			theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3.5),
									axis.text.y = element_text(angle = 90, vjust = 1, size = 3.5),
									legend.title = element_text(size = 10)) +
			scale_fill_distiller(palette = "Spectral", name = "Mutations",
																								limits = color_limits, oob = scales::squish)
		
		#print(pl)
		
		# Save the plot
		output_file <- glue("{res_path}/{basename(file_path)}.jpg")
		#ggsave(plot = pl, filename = output_file, dpi = 600, width = 10, height = 10, bg = "white")
		
		# Return a summary for this file
		list(
			file = file_path,
			mean = mean_value,
			median = median_value,
			pi_v = pi_value,
			max = max_value
		)
	}, fold, fastas, SIMPLIFY = FALSE)
	
	# Combine summaries into a single list
	summary_df <- do.call(rbind, lapply(summary_list, function(x) data.frame(file = x$file, mean = x$mean, median = x$median, pi_v = x$pi_v, max = x$max)))
	# summary_df$mean_prop <- summary_df$mean / ALN_LEN
	# summary_df$median_prop <- summary_df$median / ALN_LEN
	# summary_df$pi_prop <- summary_df$pi_v / ALN_LEN
	# View(summary_df)
	
	.plot_dist_summary_distr <- function(df, var, lbl) {
		gd_plot <- ggplot(df, aes(x = !!sym(var))) +
			geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
			theme_classic() + labs(title = glue("Distribution of {lbl} distances"), x = glue("{lbl} distance (nt mutations)"), y = "frequency")
		ggsave(plot = gd_plot, filename = glue("{res_path}/{lbl}_gen_dist_distribution.jpg"), dpi = 600, width = 8, height = 6, bg = "white")
	}
	
	.plot_dist_summary_distr(summary_df, "mean", "mean")
	.plot_dist_summary_distr(summary_df, "median", "median")
	.plot_dist_summary_distr(summary_df, "pi_v", "pi")
	.plot_dist_summary_distr(summary_df, "max", "max")
	
	# for each phylotype plot individual mean, median, and max values
	summary_df$phylotype <- as.factor(sub(".*_B_(\\d+)_.*", "\\1", summary_df$file))
	summary_df$phylotype <- fct_relevel(summary_df$phylotype,function(x){as.character(sort(as.integer(x)))})
	
	#summary_df$phylotype <- fct_reorder(summary_df$phylotype,as.integer(summary_df$phylotype))
	
	.plot_dist_summary_pts <- function(df, var, lbl) {
		pt_plot <- ggplot(df, aes(x = phylotype, y = !!sym(var))) +
			geom_bar(stat="identity") + theme_classic() + 
			theme(axis.text.x=element_text(size=5, angle=90, hjust = 1)) + labs(x = glue("phylotype"), y = glue("{lbl} genetic distance"))
		ggsave(plot = pt_plot, filename = glue("{res_path}/{lbl}_pt_vals.jpg"), dpi = 600, width = 15, height = 10, bg = "white")
	}
	
	.plot_dist_summary_pts(summary_df, "mean", "mean")
	.plot_dist_summary_pts(summary_df, "median", "median")
	.plot_dist_summary_pts(summary_df, "pi_v", "pi")
	.plot_dist_summary_pts(summary_df, "max", "max")
	
	return(summary_df)
}

# mean and pi are calculated differently, but in this case because diagonal is excluded and only unique pairs are considered (upper or lower triangle of the matrix) are considered anyway
pigd <- plot_intraphylotype_genetic_diversity("results/04_aln_all_mono_subtypes/mcs30/B")
# mean/pi around 20 muts (2%), upper values ~50 (5%)
# max around 50 muts (5%), upper values ~100 (10%)

### Run ClusterPicker ###
# exec java file: https://github.com/emmahodcroft/cluster-picker-and-cluster-matcher/raw/master/release/ClusterPicker_1.2.5.jar
HOME_DIR <- Sys.getenv("HOME")
CP_PATH <- glue("{HOME_DIR}/tools/clusterpicker/ClusterPicker_1.2.5.jar")
write_trees_gd <- function(ml_tree, out_fname) {
	write.tree(ml_tree, glue("{seqs_folder_naive}/{out_fname}"))
}
out_files <- c("A_A1_match_fasta.nwk", "CRF_02_AG_match_fasta.nwk", "C_match_fasta.nwk", "B_match_fasta.nwk")
for(i in 1:length(out_files)) write_trees_gd(ml_trees[i], out_files[i])

system(glue("java -jar {CP_PATH} {seqs_folder_naive}/B_match_tre.fasta {seqs_folder_naive}/B_match_fasta.nwk 0.8 0.8 0.1 0"))
# init_supp_thr main_supp_thr gen_dist_thr (max genetic dist allowed within clusters) large_clust_thr
# 0.8 0.8 <0.075> 30 -> 2110 clusters
# 0.8 0.8 <0.1> <0> -> 635 clusters
# IMPORTANT: not using because very large number of cluster (not comparable to treestructure)

### Run treecluster ###
# sudo pip install treecluster
gd_thrs <- seq(from=0.005, to=0.15, by=0.005) # 30 genetic distance thrs
clust_methods <- c("max_clade","avg_clade") # 2 methods
res_treecluster <- glue("{RESULTS_PATH}/03_treecluster")
system(glue("mkdir -p {res_treecluster}"))

run_treecluster <- function(params) {
	i <- params[1] # genetic dist threshold
	j <- params[2] # path to output with gd and clust_method
	k <- params[3] # subtype tree to consider
	print(glue("{i}-{j}-{k}"))
	
	system(glue("mkdir -p  {RESULTS_PATH}/03_treecluster/{tree_names[k]}"))
	command <- glue("TreeCluster.py -i {seqs_folder_naive}/{tree_names[k]}_match_fasta.nwk -o {RESULTS_PATH}/03_treecluster/{tree_names[k]}/{tree_names[k]}_gd_{gd_thrs[i]}_{clust_methods[j]}.log -t {gd_thrs[i]} -s 0.8 -m {clust_methods[j]}")
	
	system(command)
}

# combination of indices
param_combinations <- expand.grid(seq_along(gd_thrs), seq_along(clust_methods), seq_along(tree_names))

# Run treecluster in parallel
pbmclapply(seq_len(nrow(param_combinations)), function(idx) {
	print(idx)
	run_treecluster(as.numeric(param_combinations[idx, ]))
}, mc.cores = parallel::detectCores())

list_tc_outs <- list()
inspect_nclusts <- list()
for(i in 1:length(tree_names)) {
	list_tc_outs[[i]] <- list.files(path=glue("{RESULTS_PATH}/03_treecluster/{tree_names[i]}/"), full.names = T, pattern="*.log")
	inspect_nclusts[[i]] <- lapply(X=list_tc_outs[[i]], function(x) {
		#print("a")
		f <- read.csv(x, header=T, sep="\t", fileEncoding = "UTF-8")
		print(x)
		distr_sizes <- f %>% group_by(ClusterNumber) %>% summarise(n=n())
		distr_sizes <- distr_sizes[distr_sizes$ClusterNumber != -1,]
		#print(head(distr_sizes))
		
		# # plot distribution of cluster sizes
		# distr_sizes_plot <- ggplot(distr_sizes, aes(x = n)) +
		# 	geom_histogram(binwidth = 5, fill = "skyblue", color = "black", alpha = 0.7) +
		# 	theme_classic() + labs(x = glue("cluster size"), y = "frequency")
		# ggsave(plot = distr_sizes_plot, filename = glue("{res_treecluster}/{tree_names[i]}/{basename(x)}.jpg"), dpi = 600, width = 8, height = 6, bg = "white")
		
		print("n clusters")
		n_clusts <- max(f$ClusterNumber)
		print(n_clusts)
		print("max size among clusters")
		print(max(distr_sizes$n, na.rm=T))
		print("min size among clusters")
		print(min(distr_sizes$n, na.rm=T))
		
		gd_v <- as.numeric(sub(".*_gd_([0-9.]+)_.*", "\\1", basename(x)))
		method_v <- sub(".*_gd_[0-9.]+_(.*)\\.log", "\\1", basename(x))
		df_gd_method_nclust <- data.frame(gd=gd_v, method=method_v, n_clusts=n_clusts)
		print(df_gd_method_nclust)
		
		#return(df_gd_method_nclust=df_gd_method_nclust)
		#df_gd_method_nclust[[i]]
	})
}

# plot change of number of clusters based on genetic distance and method
inspect_nclusts_1 <- rbindlist(inspect_nclusts[[1]]); inspect_nclusts_2 <- rbindlist(inspect_nclusts[[2]])
inspect_nclusts_3 <- rbindlist(inspect_nclusts[[3]]); inspect_nclusts_4 <- rbindlist(inspect_nclusts[[4]])
inspect_nclusts_comb <- list(inspect_nclusts_1, inspect_nclusts_2, inspect_nclusts_3, inspect_nclusts_4)

pal_gd_thrs <- c("avg_clade"="#D55E00", "max_clade"="#56B4E9")
change_nclusts_gd_method <- function(gd_stats, range_y, treestructure_nclusts_mcs30, subtype_out) {
	print(subtype_out)
	plt <- ggplot(gd_stats, aes(x=gd, y=n_clusts, color=method)) +
		geom_line() + geom_point() + scale_color_manual(values = pal_gd_thrs, name="Genetic distance\nmethod") + scale_x_continuous(breaks=gd_thrs) + scale_y_continuous(breaks=seq(from=range_y[1],to=range_y[2],by=range_y[3])) + coord_cartesian(expand=FALSE) + geom_hline(yintercept = treestructure_nclusts_mcs30, linetype = "dashed") + 
		labs(x="Genetic distance threshold", y="Number of clusters") + #title="treecluster genetic distance vs method comparison",
		theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5,color="black"), axis.text=element_text(size=8,family=helv,color="black"), axis.title= element_text(family=helv),legend.title = element_text(family=helv), legend.text=element_text(family=helv,size=10), legend.position = c(.75,.5))
	ggsave(plot = plt, filename = glue("{res_treecluster}/{subtype_out}.jpg"), dpi = 600, width = 10, height = 8, bg = "white")
	return(list(plt=plt, gd_stats=gd_stats))
}


range_ys <- list( c(100, 1000, 100), c(100, 1000, 100), c(500, 3500, 500), c(500, 6000, 500) )
nclust_gds <- list()
for(i in 1:length(tree_names)) {
	print(i)
	current_range <- range_ys[[i]]
	#print(current_range)
	nclust_gds[[i]] <- change_nclusts_gd_method(inspect_nclusts_comb[[i]], current_range, length(unique(treestruct_min_cl_size_res_yes_sup[[1,i]][[2]]$cluster)), tree_names[i] )
}

# Figure S2
ggsave(plot = nclust_gds[[4]]$plt, filename = glue("{RESULTS_PATH}/figs/figS2.jpg"), dpi = 600, width = 8, height = 7, bg = "white")
ggsave(plot = nclust_gds[[4]]$plt, filename = glue("{RESULTS_PATH}/figs/figS2.eps"), dpi = 600, width = 8, height = 7, bg = "white")

# below focus on subtype B, mcs=30 (treestructure)

# ideal thresholds appear to be 0.075 (201 clusters) and 0.08 (120 clusters) for avg_clade
# for max_clade, only get ~200 clusters when increasing gd to 0.15, but distance too big
# clusterPicker is said to be comparable with treecluster max_clade, but for 0.05 threshold clusterpicker gives 635 clusters and treecluster gives 1361
# plot distribution of cluster sizes from treestructure 30B and see if comparable
distr_sizes_ts <- treestruct_min_cl_size_res_yes_sup[[1,4]][[2]] %>% group_by(cluster) %>% summarise(n=n())
# backbone (PT153) has 13567, removing that for vis purposes
distr_sizes_ts <- distr_sizes_ts[distr_sizes_ts$cluster != 153,]
distr_sizes_ts_plot <- ggplot(distr_sizes_ts, aes(x = n)) +
	geom_histogram(binwidth = 5, fill = "skyblue", color = "black", alpha = 0.7) +
	theme_classic() + labs(x = glue("cluster size"), y = "frequency")
ggsave(plot = distr_sizes_ts_plot, filename = glue("results/03_treestructure/30B_dist_sizes.jpg"), dpi = 600, width = 8, height = 6, bg = "white")

# distribution of cluster sizes for avg 0.075 looks more similar than 0.08

# compare visually treestructure mcs=30, fastbaps optim.symmetric, and treecluster 0.075_avg_clade
trestruct_b_comp <- treestruct_min_cl_size_res_yes_sup[[1,4]][[1]]
fastbaps_b_df <- readRDS(glue("{RDS_PATH}/fastbaps_df_list.rds"))[[4]]
treecluster_b_comp <- read.csv(glue("{res_treecluster}/B/B_gd_0.075_avg_clade.log"), sep="\t", fileEncoding = "UTF-8")
treecluster_b_comp <- treecluster_b_comp[treecluster_b_comp$ClusterNumber != -1,]

plot_comp12 <- inner_join(trestruct_b_comp$data, fastbaps_b_df, by=c("taxon"="id"))
plot_comp123 <- inner_join(plot_comp12, treecluster_b_comp, by=c("taxon"="SequenceName"))
colnames(plot_comp123) <- c("taxon","cluster_trestruct","partition_trestruct","cluster_fastbaps","cluster_treescluster")
plot_comp123$cluster_trestruct <- as.integer(plot_comp123$cluster_trestruct)
plot_comp123$cluster_fastbaps <- as.integer(plot_comp123$cluster_fastbaps)
plot_comp123$cluster_treescluster <- as.integer(plot_comp123$cluster_treescluster)

plot_comp123_cp <- plot_comp123

treeB <- readRDS("rds/tree_b_adj2.rds")
gg <- ggtree(treeB, color="grey20", size=0.1) + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + expand_limits(y = -25)

gg2 <- ggtree(treeB, color="grey20", size=0.1) %<+% plot_comp123_cp
gg3 <- gg2 + #%<+% plot_comp123_cp +
	geom_tippoint(aes(subset=(cluster_trestruct==40 & !is.na(cluster_trestruct))),fill="#56B4E9",size=2, stroke=0.15, shape=21) +
	geom_tippoint(aes(subset=(cluster_trestruct==69 & !is.na(cluster_trestruct))),fill="#D55E00",size=2, stroke=0.15, shape=21) +
	geom_tippoint(aes(subset=(cluster_trestruct==133 & !is.na(cluster_trestruct))),fill="#009E73",size=2, stroke=0.15, shape=21) +
	theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + expand_limits(y = -25)
ggsave(plot=gg3, file=glue("{RESULTS_PATH}/figs/figS3_aux.svg"), dpi=900, width=8, height=10, bg="white")

fcomp1 <- facet_plot(gg, panel = "treestructure", data = plot_comp123, geom = geom_tile, 
																					aes(x = cluster_trestruct), color = "black") #fill = highlight)

fcomp1 <- facet_plot(fcomp1, panel = "fastbaps", data = plot_comp123, geom = geom_tile, 
																					aes(x = cluster_fastbaps), color = "brown")
fcomp1 <- facet_plot(fcomp1, panel = "treecluster", data = plot_comp123, geom = geom_tile, 
																					aes(x = cluster_treescluster), color = "#999999")
fcomp1
ggsave(plot=fcomp1, file=glue("{RESULTS_PATH}/figs/figS3_to_edit.svg"), dpi=600, width=15, height=10, bg="white")
# Added red boxed manually in inkscape