libs_load <- c("dplyr","openxlsx", "glue","ape","treedater","lubridate","data.table","ggtree","ggpubr","ggplot2","viridis","grid")
invisible( lapply(libs_load, library, character.only=TRUE) )

NCPU <- 7

seqs_folder_naive <- "data/subtype_seqs_naive"
output_iqtree_folder <- "results/iqtree"
fasta_b <- read.dna(glue("{seqs_folder_naive}/B_curated_refC_aln_len_filter.fasta"), format="fasta")

# Add "t" prefix to samples (tips in the tree)
trees_b_path <- c(); for(i in 1:10) { trees_b_path[i] <- glue("rds/tree_subsamples/tree_B{i}_adj.rds") }
timetrees_b_path <- c(); for(i in 1:10) { timetrees_b_path[i] <- glue("rds/tree_subsamples/timetree_B{i}_adj.rds") }
trees_b_adj <- timetrees_b_adj <- list()
for(f in 1:length(trees_b_path)) {
	trees_b_adj[[f]] <- readRDS(trees_b_path[[f]])
	timetrees_b_adj[[f]] <- readRDS(timetrees_b_path[[f]])
	#trees_b_adj[[f]]$tip.label <- paste0("t.",trees_b_adj[[f]]$tip.label)
}

# Gather together all tips that passed outlier removal
all_tips_filtered <- sapply(trees_b_adj, "[[",5) #drop levels (5th index is tip.label)
all_tips_filtered <- unlist(all_tips_filtered) # unlist to get only one vector (n=25210)
all_tips_filtered <- unique(all_tips_filtered) # remove duplicated names across diff replicates (reference seqs) (n=25201)
# 25282 - 25201 = 81 removed because outliers

#md_subtype_b_uk_only_dates_filtered <- md_subtype_b_uk_only_dates[md_subtype_b_uk_only_dates$sequence_name %in% all_tips_filtered,]
# Retrieve fasta aligned file with all sequences
rownames(fasta_b) <- paste0("t.",rownames(fasta_b))
fasta_b_outl_rm <- fasta_b[rownames(fasta_b) %in% all_tips_filtered,]
write.FASTA(fasta_b_outl_rm, glue("{seqs_folder_naive}/B_UK_final_aln_outliers_removed.fasta"))

# Estimate full tree using iqtree2

system(glue("{iqtree_bin} -s {seqs_folder_naive}/B_UK_final_aln_outliers_removed.fasta -m GTR+R -nt AUTO -ntmax 5 -B 1000 -nm 5000 -bcor 0.97"))
system(glue("mv {seqs_folder_naive}/*.bionj {seqs_folder_naive}/*.gz {seqs_folder_naive}/*.iqtree {seqs_folder_naive}/*.log {seqs_folder_naive}/*.mldist {seqs_folder_naive}/*.treefile {output_iqtree_folder}"))

# move treefile from server to local, load below
tree_b_adj <- read.tree(glue("{output_iqtree_folder}/B_UK_final_aln_outliers_removed.fasta.treefile")) #25201
demog_md_subtype_match_naive_uq_dates_list <- readRDS("rds/demog_md_subtype_match_naive_uq_dates_list.rds")
md_subtype_b_uk_only <- demog_md_subtype_match_naive_uq_dates_list$B
md_subtype_b_uk_only$testindex <- paste0("t.",md_subtype_b_uk_only$testindex)
md_subtype_b_uk_only_filtered <- md_subtype_b_uk_only[md_subtype_b_uk_only$testindex %in% tree_b_adj$tip.label,]

md_subtype_b_uk_only_filtered <- subset(md_subtype_b_uk_only_filtered, select=c("testindex", "dbsample_date")) #patientindex
colnames(md_subtype_b_uk_only_filtered) <- c("sequence_name","sample_date")
root_tip_c <- "t.Ref.C.ET.86.ETH2220.U46016"
root_time_c <- 1986.5
md_subtype_b_uk_only_filtered[nrow(md_subtype_b_uk_only_filtered) + 1,] <- list(root_tip_c, root_time_c)
md_subtype_b_uk_only_filtered$sample_date <- date_decimal(md_subtype_b_uk_only_filtered$sample_date)
md_subtype_b_uk_only_filtered$sample_date <- as.Date(md_subtype_b_uk_only_filtered$sample_date)

subtype_b_sampleTimes <- decimal_date(md_subtype_b_uk_only_filtered$sample_date) #25201
names(subtype_b_sampleTimes) <- md_subtype_b_uk_only_filtered$sequence_name

tree_b_adj <- ape::root(tree_b_adj, outgroup=glue("{root_tip_c}"), resolve.root=TRUE)

# move 3 following files to remove PC (kingman)
saveRDS(tree_b_adj, "rds/tree_b_adj.rds")
saveRDS(subtype_b_sampleTimes, "rds/subtype_b_sampleTimes.rds")

##### NOTE: had to use ape dev version >=5.7.1.2 because of issue with dist.nodes [tree too big], which is enhanced for larger trees in dev version
raw_clock_b <- readRDS("rds/raw_clock_b.rds")
# Unable to converge with additive model
#timetree_b <- dater( tree_b_adj, subtype_b_sampleTimes, s=n_sites_aln, clock='additive', omega0=unlist(raw_clock_b), ncpu=NCPU) #median(unlist(raw_clock_rates))
# doi.org/10.1093/ve/vex029 (rate 0.0014 - 0.0017 subset 1, 0.0012 - 0.0015 subset 2)
# range (raw_clock_b): 0.001072383 0.001240094
# median(raw_clock_b): 0.001152089
timetree_b <- dater( tree_b_adj, subtype_b_sampleTimes, s=n_sites_aln, clock='strict', quiet=FALSE, maxit=10, omega0=median(raw_clock_b), meanRateLimits=c(0.0010, 0.0013), parallel_foreach = TRUE, ncpu=NCPU)
saveRDS(timetree_b, "rds/timetree_b.rds")

# after first dating: Load resulting tree from HPC (B_full_timetree_strict.pbs)
timetree_b <- readRDS("rds/timetree_b_strict.rds")
print(timetree_b) #1953
#pdf("results/rtt/timetree_b_tips.pdf")
rootToTipRegressionPlot(timetree_b, show.tip.labels = F) #0.154
#dev.off()
tail_prob <- 0.0025
outliers_b <- outlierTips(timetree_b, alpha=tail_prob) #with alpha=0.05 gives 4055 outliers, 0.025 -> 3649, 0.01 -> 3265, 0.005 -> 3066, 0.0025 -> 1100
# Remove all tips that don't have a high q-value (except subtype C reference sequence)
to_rm_b <- outliers_b[ (outliers_b$taxon != "t.Ref.C.ET.86.ETH2220.U46016") & (outliers_b$q < tail_prob) ,]
print("TO RM:")
print(nrow(to_rm_b)) 

tree_b_adj <- readRDS("rds/tree_b_adj.rds")
tree_b_adj2 <- drop.tip(tree_b_adj, rownames(to_rm_b))
saveRDS(tree_b_adj2, "rds/tree_b_adj2.rds")
timetree_b_adj <- dater( tree_b_adj2, subtype_b_sampleTimes, s=n_sites_aln, clock='strict', quiet=FALSE, maxit=10, omega0=median(raw_clock_b), meanRateLimits=c(0.0010, 0.0013), parallel_foreach = TRUE, ncpu=NCPU) # done on HPC

# after second dating: Load resulting tree without outliers from HPC (B_full_timetree_strict_adj.pbs)
timetree_b_adj <- readRDS("rds/timetree_b_strict_adj.rds")
print(timetree_b_adj) 
rootToTipRegressionPlot(timetree_b_adj) #0.164
tail_prob <- 0.0025
outliers_b_adj <- outlierTips(timetree_b_adj, alpha=tail_prob)
to_rm_b <- outliers_b_adj[ (outliers_b_adj$taxon != "t.Ref.C.ET.86.ETH2220.U46016") & (outliers_b_adj$q < tail_prob) ,]
print("TO RM:") 
print(nrow(to_rm_b))
