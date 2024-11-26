libs_load <- c("dplyr","openxlsx", "glue","ape","treedater","lubridate","data.table","ggtree","ggpubr","ggplot2","viridis")
invisible( lapply(libs_load, library, character.only=TRUE) )

seqs_folder_naive <- "data/subtype_seqs_naive"
output_iqtree_folder <- "results/iqtree"

# Continue here for analysis of subtype A trees (after ~ line 220 of main script)
fasta_crf02ag <- read.dna(glue("{seqs_folder_naive}/CRF_02_AG_curated_refB_aln_len_filter.fasta"), format="fasta")
subtype_crf02ag_tree <- read.tree(glue("{output_iqtree_folder}/CRF_02_AG_curated_refB_aln_len_filter.fasta.treefile"))
subtype_crf02ag_tree$tip.label <- paste0("t.",subtype_crf02ag_tree$tip.label)

root_tip_b <- "Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
root_time_b <- 1983.5
demog_md_subtype_match_naive_uq_dates_list <- readRDS("rds/demog_md_subtype_match_naive_uq_dates_list.rds")
md_subtype_crf02ag_uk_only <- demog_md_subtype_match_naive_uq_dates_list$`CRF 02_AG`
md_subtype_crf02ag_uk_only <- md_subtype_crf02ag_uk_only[md_subtype_crf02ag_uk_only$testindex %in% rownames(fasta_crf02ag),]
md_subtype_crf02ag_uk_only_dates <- subset(md_subtype_crf02ag_uk_only, select=c("testindex", "dbsample_date")) #patientindex
colnames(md_subtype_crf02ag_uk_only_dates) <- c("sequence_name","sample_date")
md_subtype_crf02ag_uk_only_dates$sequence_name <- as.character(md_subtype_crf02ag_uk_only_dates$sequence_name)
md_subtype_crf02ag_uk_only_dates[nrow(md_subtype_crf02ag_uk_only_dates) + 1,] <- list(root_tip_b, root_time_b)
md_subtype_crf02ag_uk_only_dates$sequence_name <- paste0("t.",md_subtype_crf02ag_uk_only_dates$sequence_name)
md_subtype_crf02ag_uk_only_dates$sample_date <- date_decimal(md_subtype_crf02ag_uk_only_dates$sample_date)
md_subtype_crf02ag_uk_only_dates$sample_date <- as.Date(md_subtype_crf02ag_uk_only_dates$sample_date)

out_trees <- "results/02_explore_trees"
md_subtype_crf02ag_combined <- md_subtype_crf02ag_uk_only_dates[md_subtype_crf02ag_uk_only_dates$sequence_name %in% subtype_crf02ag_tree$tip.label,]
subtype_crf02ag_tree <- keep.tip(subtype_crf02ag_tree, intersect(subtype_crf02ag_tree$tip.label, md_subtype_crf02ag_combined$sequence_name))
subtype_crf02ag_sampleTimes <- decimal_date(md_subtype_crf02ag_combined$sample_date)
names(subtype_crf02ag_sampleTimes) <- md_subtype_crf02ag_combined$sequence_name
ggtree(subtype_crf02ag_tree) + theme_tree()
system(glue("mkdir -p {out_trees}"))
ggsave(glue("{out_trees}/CRF_02_AG.png"), width=10, height=15, dpi=300, bg="white", limitsize=FALSE)

# Root-to-tip regression to see if evolution could be well explained by molecular clock model
out_rtt_mltr <- "results/02_rtt"
subtype_crf02ag_tree <- ape::root(subtype_crf02ag_tree, outgroup=glue("t.{root_tip_b}"), resolve.root=TRUE) #TRUE
print(is.rooted(subtype_crf02ag_tree))
d2root <- node.depth.edgelength( subtype_crf02ag_tree )[ 1:length(subtype_crf02ag_tree$tip.label) ]
if(!dir.exists(glue("{out_rtt_mltr}")))
	suppressWarnings( dir.create(glue("{out_rtt_mltr}")) )
png(glue('{out_rtt_mltr}/refB_root_CRF_02_AG.png'))
# a scatter plot with local regression line
scatter.smooth( x=subtype_crf02ag_sampleTimes[subtype_crf02ag_tree$tip.label], y=d2root, xlab="Sampling date", ylab="Root-to-tip distance" )
text(subtype_crf02ag_sampleTimes[subtype_crf02ag_tree$tip.label], d2root, labels=names(subtype_crf02ag_sampleTimes[subtype_crf02ag_tree$tip.label]), cex= 0.7, pos=3)
#print(summary( lm( d2root ~ subtype_crf02ag_sampleTimes[subtype_crf02ag_tree$tip.label] ) ))
raw_clock_crf02ag <- unname( coef( lm( d2root ~  subtype_crf02ag_sampleTimes[subtype_crf02ag_tree$tip.label] ) ))[2] #0.001050625
dev.off()
saveRDS(subtype_crf02ag_sampleTimes, "rds/subtype_crf02ag_sampleTimes.rds")

ggtree(subtype_crf02ag_tree) + theme_tree()
system(glue("mkdir -p {out_trees}"))
ggsave(glue("{out_trees}/ref_rerootB_CRF_02AG.png"), width=10, height=15, dpi=300, bg="white", limitsize=FALSE)

# timetree and outlier removal
n_sites_aln_crf02ag <- ncol(fasta_crf02ag)
tail_prob <- 0.05 #  tail probability used for classifying seqs as outliers
# doi.org/10.1093/ve/vex029 (rate 0.0006 - 0.0011)
timetree_crf02ag <- dater( subtype_crf02ag_tree, subtype_crf02ag_sampleTimes, s=n_sites_aln_crf02ag, clock='additive', omega0=c(raw_clock_crf02ag), ncpu=7 ) #1948.70, 7.03
png(glue('{out_rtt_mltr}/timetree_crf02ag_rtt.png'))
rootToTipRegressionPlot(timetree_crf02ag) #0.11
dev.off()
outliers_crf02ag <- outlierTips(timetree_crf02ag, alpha=tail_prob)
# Remove all tips that don't have a high q-value (except subtype B reference sequence)
to_rm_crf02ag <- outliers_crf02ag[ (outliers_crf02ag$taxon != glue("t.{root_tip_b}")) & (outliers_crf02ag$q < tail_prob) ,] 
print("TO RM:") #
print(nrow(to_rm_crf02ag)) # 101
tree_crf02ag_adj <- drop.tip(subtype_crf02ag_tree, rownames(to_rm_crf02ag))
timetree_crf02ag_adj <- dater( tree_crf02ag_adj, subtype_crf02ag_sampleTimes, s=n_sites_aln_crf02ag, clock='additive', omega0=c(raw_clock_crf02ag), ncpu=7 ) #1949.76, 7.22
rootToTipRegressionPlot(timetree_crf02ag_adj) #0.12
saveRDS(tree_crf02ag_adj, "rds/tree_crf02ag_adj.rds")
saveRDS(subtype_crf02ag_sampleTimes, "rds/subtype_crf02ag_sampleTimes.rds")
saveRDS(timetree_crf02ag_adj, "rds/timetree_crf02ag_adj.rds")

library(treeio)
tr_crf02ag <- read.iqtree(glue("{output_iqtree_folder}/CRF_02_AG_curated_refB_aln_len_filter.fasta.treefile"))
ggtree(tr_crf02ag) +  geom_text(aes(label = UFboot), hjust = 1, vjust = -0.4, size = 1, color="red") #geom_nodelab(aes(label = UFboot, color="red", vjust=-.5, size=2))
ggsave(glue("results/02_explore_trees/CRF02AG_UBoot.pdf"), width=12, height=20, dpi=300, bg="white")
