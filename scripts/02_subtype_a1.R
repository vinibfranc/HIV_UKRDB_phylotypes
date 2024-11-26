libs_load <- c("dplyr","openxlsx", "glue","ape","treedater","lubridate","data.table","ggtree","ggpubr","ggplot2","viridis")
invisible( lapply(libs_load, library, character.only=TRUE) )

seqs_folder_naive <- "data/subtype_seqs_naive"
output_iqtree_folder <- "results/02_iqtree"

# Continue here for analysis of subtype A trees (after 01 script)
fasta_a1 <- read.dna(glue("{seqs_folder_naive}/A_A1_curated_refB_aln_len_filter.fasta"), format="fasta")
subtype_a1_tree <- read.tree(glue("{output_iqtree_folder}/A_A1_curated_refB_aln_len_filter.fasta.treefile"))
subtype_a1_tree$tip.label <- paste0("t.",subtype_a1_tree$tip.label)

root_tip_b <- "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
root_time_b <- 1983.5
demog_md_subtype_match_naive_uq_dates_list <- readRDS("rds/demog_md_subtype_match_naive_uq_dates_list.rds")
md_subtype_a1_uk_only <- demog_md_subtype_match_naive_uq_dates_list$`A (A1)`
md_subtype_a1_uk_only <- md_subtype_a1_uk_only[md_subtype_a1_uk_only$testindex %in% rownames(fasta_a1),]
md_subtype_a1_uk_only_dates <- subset(md_subtype_a1_uk_only, select=c("testindex", "dbsample_date")) #patientindex
colnames(md_subtype_a1_uk_only_dates) <- c("sequence_name","sample_date")
md_subtype_a1_uk_only_dates$sequence_name <- as.character(md_subtype_a1_uk_only_dates$sequence_name)
md_subtype_a1_uk_only_dates[nrow(md_subtype_a1_uk_only_dates) + 1,] <- list(root_tip_b, root_time_b)
md_subtype_a1_uk_only_dates$sequence_name <- paste0("t.",md_subtype_a1_uk_only_dates$sequence_name)
md_subtype_a1_uk_only_dates$sample_date <- date_decimal(md_subtype_a1_uk_only_dates$sample_date)
md_subtype_a1_uk_only_dates$sample_date <- as.Date(md_subtype_a1_uk_only_dates$sample_date)

out_trees <- "results/02_explore_trees"
md_subtype_a1_combined <- md_subtype_a1_uk_only_dates[md_subtype_a1_uk_only_dates$sequence_name %in% subtype_a1_tree$tip.label,]
subtype_a1_tree <- keep.tip(subtype_a1_tree, intersect(subtype_a1_tree$tip.label, md_subtype_a1_combined$sequence_name))
subtype_a1_sampleTimes <- decimal_date(md_subtype_a1_combined$sample_date)
names(subtype_a1_sampleTimes) <- md_subtype_a1_combined$sequence_name
ggtree(subtype_a1_tree) + theme_tree()
system(glue("mkdir -p {out_trees}"))
ggsave(glue("{out_trees}/A_A1.png"), width=10, height=15, dpi=300, bg="white", limitsize=FALSE)

# Root-to-tip regression to see if evolution could be well explained by molecular clock model
out_rtt_mltr <- "results/02_rtt"
subtype_a1_tree <- ape::root(subtype_a1_tree, outgroup=glue("t.{root_tip_b}"), resolve.root=TRUE) #TRUE
print(is.rooted(subtype_a1_tree))
d2root <- node.depth.edgelength( subtype_a1_tree )[ 1:length(subtype_a1_tree$tip.label) ]
if(!dir.exists(glue("{out_rtt_mltr}")))
	suppressWarnings( dir.create(glue("{out_rtt_mltr}")) )
png(glue('{out_rtt_mltr}/refB_root_A1.png'))
# a scatter plot with local regression line
scatter.smooth( x=subtype_a1_sampleTimes[subtype_a1_tree$tip.label], y=d2root, xlab="Sampling date", ylab="Root-to-tip distance" )
text(subtype_a1_sampleTimes[subtype_a1_tree$tip.label], d2root, labels=names(subtype_a1_sampleTimes[subtype_a1_tree$tip.label]), cex= 0.7, pos=3)
#print(summary( lm( d2root ~ subtype_a1_sampleTimes[subtype_a1_tree$tip.label] ) ))
raw_clock_a1 <- unname(coef( lm( d2root ~  subtype_a1_sampleTimes[subtype_a1_tree$tip.label] ) ))[2] #0.001406869
dev.off()
#saveRDS(subtype_a1_sampleTimes, "rds/subtype_a1_sampleTimes.rds")

ggtree(subtype_a1_tree) + theme_tree()
system(glue("mkdir -p {out_trees}"))
ggsave(glue("{out_trees}/ref_rerootB_A1.png"), width=10, height=15, dpi=300, bg="white", limitsize=FALSE)

# timetree and outlier removal
n_sites_aln_a1 <- ncol(fasta_a1)
tail_prob <- 0.05 #  tail probability used for classifying seqs as outliers
# doi.org/10.1093/ve/vex029 (rate 0.0014 - 0.0017)
timetree_a1 <- dater( subtype_a1_tree, subtype_a1_sampleTimes, s=n_sites_aln_a1, clock='additive', omega0=c(raw_clock_a1), ncpu=7 ) # 1959.13, 7.43
png(glue('{out_rtt_mltr}/timetree_a1_rtt.png'))
rootToTipRegressionPlot(timetree_a1) #0.051
dev.off()
outliers_a1 <- outlierTips(timetree_a1, alpha=tail_prob) #prev=0.1
# Remove all tips that don't have a high q-value (except subtype B reference sequence)
to_rm_a1 <- outliers_a1[ (outliers_a1$taxon != glue("t.{root_tip_b}")) & (outliers_a1$q < tail_prob) ,]
print("TO RM:")
print(nrow(to_rm_a1)) ##122
tree_a1_adj <- drop.tip(subtype_a1_tree, rownames(to_rm_a1))
timetree_a1_adj <- dater( tree_a1_adj, subtype_a1_sampleTimes, s=n_sites_aln_a1, clock='additive', omega0=c(raw_clock_a1), ncpu=7 ) #1959.65, 6.32
rootToTipRegressionPlot(timetree_a1_adj) #0.053
#class(timetree_a1_adj) <- "phylo"

saveRDS(tree_a1_adj, "rds/tree_a1_adj.rds")
saveRDS(subtype_a1_sampleTimes, "rds/subtype_a1_sampleTimes.rds")
saveRDS(timetree_a1_adj, "rds/timetree_a1_adj.rds")

library(treeio)
tr_a1 <- read.iqtree(glue("{output_iqtree_folder}/A_A1_curated_refB_aln_len_filter.fasta.treefile"))
ggtree(tr_a1) +  geom_text(aes(label = UFboot), hjust = 1, vjust = -0.4, size = 1, color="red") #geom_nodelab(aes(label = UFboot, color="red", vjust=-.5, size=2))
ggsave(glue("results/02_explore_trees/A1_UBoot.pdf"), width=12, height=20, dpi=300, bg="white")
