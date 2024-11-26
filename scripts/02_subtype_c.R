# NOTE: after running timetrees and removing outliers for each individual subsample from folder `b_c_subsamples`, run here

libs_load <- c("dplyr","openxlsx", "glue","ape","treedater","lubridate","data.table","ggtree","ggpubr","ggplot2","viridis")
invisible( lapply(libs_load, library, character.only=TRUE) )

NCPU <- 7

seqs_folder_naive <- "data/subtype_seqs_naive"
output_iqtree_folder <- "results/iqtree"
fasta_c <- read.dna(glue("{seqs_folder_naive}/C_curated_refB_aln_len_filter.fasta"), format="fasta")

# Add "t" prefix to samples (tips in the tree)
trees_c_path <- c(); for(i in 1:5) { trees_c_path[i] <- glue("rds/tree_C{i}_adj.rds") }
timetrees_c_path <- c(); for(i in 1:5) { timetrees_c_path[i] <- glue("rds/timetree_C{i}_adj.rds") }
trees_c_adj <- timetrees_c_adj <- list()
for(f in 1:length(trees_c_path)) {
	trees_c_adj[[f]] <- readRDS(trees_c_path[[f]])
	timetrees_c_adj[[f]] <- readRDS(timetrees_c_path[[f]])
	#trees_c_adj[[f]]$tip.label <- paste0("t.",trees_c_adj[[f]]$tip.label)
}

# Gather together all tips that passed outlier removal
all_tips_filtered_c <- sapply(trees_c_adj, "[[",5) #drop levels
all_tips_filtered_c <- unlist(all_tips_filtered_c) # unlist to get only one vector (n=11700)
all_tips_filtered_c <- unique(all_tips_filtered_c) # remove duplicated names across diff replicates (reference seqs) (n=11696)
# 11711 - 11696 = 15 removed as outliers

#md_subtype_c_uk_only_dates_filtered <- md_subtype_c_uk_only_dates[md_subtype_c_uk_only_dates$sequence_name %in% all_tips_filtered_c,]
# Retrieve fasta aligned file with all sequences
rownames(fasta_c) <- paste0("t.",rownames(fasta_c))
fasta_c_outl_rm <- fasta_c[rownames(fasta_c) %in% all_tips_filtered_c,]
write.FASTA(fasta_c_outl_rm, glue("{seqs_folder_naive}/C_UK_final_aln_outliers_removed.fasta"))

##### NOTE: here will need to move C_UK_final_aln_outliers_removed.fasta and pbs file and run (probably try -bcor=0.80, 0.85, 0.90)

system(glue("{iqtree_bin} -s {seqs_folder_naive}/C_UK_final_aln_outliers_removed.fasta -m GTR+R -nt AUTO -ntmax 5 -B 1000 -nm 5000 -bcor 0.90"))
system(glue("mv {seqs_folder_naive}/*.bionj {seqs_folder_naive}/*.gz {seqs_folder_naive}/*.iqtree {seqs_folder_naive}/*.log {seqs_folder_naive}/*.mldist {seqs_folder_naive}/*.treefile {output_iqtree_folder}"))

tree_c_adj <- read.tree(glue("{output_iqtree_folder}/C_UK_final_aln_outliers_removed.fasta.treefile"))
demog_md_subtype_match_naive_uq_dates_list <- readRDS("rds/demog_md_subtype_match_naive_uq_dates_list.rds")
md_subtype_c_uk_only <- demog_md_subtype_match_naive_uq_dates_list$C
md_subtype_c_uk_only$testindex <- paste0("t.",md_subtype_c_uk_only$testindex)
md_subtype_c_uk_only_filtered <- md_subtype_c_uk_only[md_subtype_c_uk_only$testindex %in% tree_c_adj$tip.label,]

md_subtype_c_uk_only_filtered <- subset(md_subtype_c_uk_only_filtered, select=c("testindex", "dbsample_date")) #patientindex
colnames(md_subtype_c_uk_only_filtered) <- c("sequence_name","sample_date")
root_tip_b <- "t.Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
root_time_b <- 1983.5
md_subtype_c_uk_only_filtered[nrow(md_subtype_c_uk_only_filtered) + 1,] <- list(root_tip_b, root_time_b)
md_subtype_c_uk_only_filtered$sample_date <- date_decimal(md_subtype_c_uk_only_filtered$sample_date)
md_subtype_c_uk_only_filtered$sample_date <- as.Date(md_subtype_c_uk_only_filtered$sample_date)

subtype_c_sampleTimes <- decimal_date(md_subtype_c_uk_only_filtered$sample_date)
names(subtype_c_sampleTimes) <- md_subtype_c_uk_only_filtered$sequence_name

tree_c_adj <- ape::root(tree_c_adj, outgroup=glue("{root_tip_b}"), resolve.root=TRUE)

saveRDS(tree_c_adj, "rds/tree_c_adj.rds")
saveRDS(subtype_c_sampleTimes, "rds/subtype_c_sampleTimes.rds")

##### NOTE: here will need to move C_UK_final_aln_outliers_removed.fasta.treefile, subtype_c_sampleTimes.rds, raw_clock_rates_c.rds, R file, pbs file, s=995, ncpu=7 to remove PC (kingman)
raw_clock_c <- readRDS("rds/raw_clock_c.rds")

timetree_c <- dater( tree_c_adj, subtype_c_sampleTimes, s=995, clock='additive', omega0=unlist(raw_clock_c), ncpu=4 ) #omega0=median(unlist(raw_clock_rates_c))
saveRDS(timetree_c_adj, "rds/timetree_c.rds")

# Load resulting tree from remote
timetree_c <- readRDS("rds/timetree_c.rds")
print(timetree_c) # 1939, 12.3
rootToTipRegressionPlot(timetree_c) #0.054
outliers_c <- outlierTips(timetree_c, alpha=0.05)
# Remove all tips that don't have a high q-value (except subtype B reference sequence)
to_rm_c <- outliers_c[ (outliers_c$taxon != "t.Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455") & (outliers_c$q < 0.05) ,]
print("TO RM:")
print(nrow(to_rm_c)) #364
#subtype_c_tree <- read.tree(glue("{output_iqtree_folder}/C_UK_final_aln_outliers_removed.fasta.treefile"))
tree_c_adj <- readRDS("rds/tree_c_adj.rds")
tree_c_adj2 <- drop.tip(tree_c_adj, rownames(to_rm_c))
saveRDS(tree_c_adj2, "rds/tree_c_adj2.rds")

# Move tree_c_adj2.rds, pbs and R files to remote and run
timetree_c <- dater( tree_c_adj2, subtype_c_sampleTimes, s=995, clock='additive', omega0=unlist(raw_clock_c), ncpu=4 ) # done on HPC

# Load resulting tree without outliers from remote
timetree_c_adj <- readRDS("rds/timetree_c_adj.rds")
print(timetree_c_adj) # 1938.1, 8.31
rootToTipRegressionPlot(timetree_c_adj) #0.054
outliers_c_adj <- outlierTips(timetree_c_adj, alpha=0.05)
to_rm_c <- outliers_c_adj[ (outliers_c_adj$taxon != "t.Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455") & (outliers_c_adj$q < 0.05) ,]
print("TO RM:") 
print(nrow(to_rm_c)) #0