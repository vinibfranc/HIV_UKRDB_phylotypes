# Calculate distance between consensus sequences of phylotype
libs_load <- c("glue", "ape", "RColorBrewer", "reshape2", "ggplot2", "viridis", "readr", "ggtree", "magrittr", "ggnewscale", "viridis", "dplyr")
invisible( lapply(libs_load, library, character.only=TRUE) )

cons_folder <- "15_consensus"

# Deleted nucleotides after position 1302 and masked ends of aln using aliview (pt_consensus_adj.fasta)
# keep only subtype B consensus phylotype pol sequences
fasta_all_cons <- read.dna(glue("{RESULTS_PATH}/{cons_folder}/pt_consensus_adj.fasta"), format="fasta" ) 
fasta_b_cons <- fasta_all_cons[ grepl("^B@", rownames(fasta_all_cons)), ]
write.dna(fasta_b_cons, glue("{RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB.fasta"), format="fasta")

# manually get 1302 and 995 nt length sequences (pt_consensus_adj_subtypeB_1302.fasta and pt_consensus_adj_subtypeB_995.fasta)

RDS_PATH="rds"
RESULTS_PATH="results"
aln_len1 <- 1302
aln_len2 <- 995

# compute distance matrices
pt_cons1 <- read.FASTA(glue("{RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB_1302.fasta") ) 
D_raw1 <- dist.dna( pt_cons1, model = 'raw' , as.matrix = T, pairwise.deletion=TRUE) * aln_len1
diag(D_raw1) <- NA # don't count zero distances on diagonal
hist(D_raw1)
mean(D_raw1,na.rm=T) #66.07

D_raw1_mlt <- melt(D_raw1)

ggplot(data = D_raw1_mlt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
	labs(x="",y="") + theme_bw() + theme(axis.text.x = element_text(angle = 40,hjust = 1)) +
	scale_fill_viridis(direction = -1, name="nt muts")

pt_cons2 <- read.FASTA(glue("{RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB_995.fasta") ) 

D_raw2 <- dist.dna( pt_cons2, model = 'raw' , as.matrix = T, pairwise.deletion=TRUE) * aln_len2
diag(D_raw2) <- NA
hist(D_raw2)
mean(D_raw2,na.rm=T) #48.4
median(D_raw2,na.rm=T) #48

D_raw2_mlt <- melt(D_raw2)
D_raw2_mlt$Var1 <- gsub("B@", "", D_raw2_mlt$Var1); D_raw2_mlt$Var2 <- gsub("B@", "", D_raw2_mlt$Var2)
D_raw2_mlt$Var1 <- as.integer(D_raw2_mlt$Var1)
D_raw2_mlt$Var2 <- as.integer(D_raw2_mlt$Var2)
D_raw2_mlt <- D_raw2_mlt[order(D_raw2_mlt$Var1, D_raw2_mlt$Var2),]

### BEGIN FIGURE S7 ###
library(scales)
custom_palette <- gradient_n_pal(c("#00468B", "#56B4E9", "#009E73", "#F0E442", "#D55E00"))
ggplot(data = D_raw2_mlt, aes(x=as.factor(Var1), y=as.factor(Var2), fill=value)) + geom_tile() +
	labs(x="",y="") + theme_classic() +	theme(axis.text.x = element_text(angle = 90,hjust = 1, size=3.5, family=helv), axis.text.y=element_text(angle = 90,vjust = 1,size=3.5, family=helv), legend.title=element_text(family = helv, size=10)) +
	#scale_fill_distiller(palette="Spectral", name="Mutations" ) + scale_fill_viridis(option = "cividis") + scale_fill_distiller(palette = "RdBu") 
	scale_fill_gradientn(colours = custom_palette(seq(0, 1, length.out = 5)), name="# mutations") + leg
suppressMessages( ggsave(file=glue("{RESULTS_PATH}/figs/figS7.eps"), dpi=600, width=10, height=10, bg="white") )# limitsize=FALSE,
suppressMessages( ggsave(file=glue("{RESULTS_PATH}/figs/figS7.jpg"), dpi=600, width=10, height=10, bg="white") )#limitsize=FALSE, 
### END FIGURE S4 ###

# Build tree with consensus seqs of each phylotype + VB + refs B and C

# Generate consensus of VB clade
# First, align all 15 VB sequences
system(glue("mafft --thread 4 --auto {DATA_PATH}/refs/vb_clade.fasta > {DATA_PATH}/refs/vb_clade_aln.fasta"))

# read sequences
vb <- read.alignment(glue('{DATA_PATH}/refs/vb_clade_aln.fasta'), format = 'fasta')
vb_cons <- seqinr::consensus(vb,  method = 'majority', threshold = .8)
write.fasta( vb_cons, names=c("VB_consensus"), file = glue('{RESULTS_PATH}/{cons_folder}/vb_consensus.fasta'))

# align against ref B, VB clade sample, ref C
system(glue("mafft --6merpair --thread 4 --addfragments {RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB_995.fasta data/refs/HIV1_REF_2020_pol_DNA_SUBTYPE_B.fasta > {RESULTS_PATH}/14_consensus/pt_consensus_adj_subtypeB_995_refB.fasta"))
# trim ends manually before running below
system(glue("mafft --6merpair --thread 4 --addfragments {RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB_995_refB.fasta {RESULTS_PATH}/14_consensus/vb_consensus.fasta > {RESULTS_PATH}/14_consensus/pt_consensus_adj_subtypeB_995_refB_vb.fasta"))
# trim ends manually before running below
system(glue("mafft --6merpair --thread 4 --addfragments {RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB_995_refB_vb.fasta data/refs/HIV1_REF_2020_pol_DNA_SUBTYPE_C.fasta > {RESULTS_PATH}/14_consensus/pt_consensus_adj_subtypeB_995_refB_vb_refC.fasta"))
# trim ends manually before running below
# rename manually refs and VB to make names shorter (refC.1986.ETH2220, VB.consensus, refB.1983.HXB2)

all_cons <- read.FASTA(glue("{RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB_995_refB_vb_refC.fasta"))
D_raw3 <- dist.dna( all_cons, model = 'raw' , as.matrix = T, pairwise.deletion=TRUE) * aln_len2
diag(D_raw3) <- NA # don't count zero distances on diagonal
D_raw3["B@40","VB.consensus"] #dist VB = 56 (95% matches imply up to 50 muts)
mean(D_raw3["B@40",],na.rm=T) # mean distance to all: 53.75
D_raw3["B@69","VB.consensus"] # dist VB = 52 (93.5% matches imply up to 70 muts)
mean(D_raw3["B@69",],na.rm=T) # mean distance to all: 47.43
D_raw3["B@133","VB.consensus"] # dist VB = 49 (93.5% matches imply up to 70 muts)
mean(D_raw3["B@133",],na.rm=T) # mean distance to all: 45.8

# Build tree
system(glue("{iqtree_bin} -s {RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB_995_refB_vb_refC.fasta -m GTR+R -nt AUTO -B 1000 -nm 5000"))

### BEGIN FIGURE 1 ###
ml_cons_unr <- read.tree(glue("{RESULTS_PATH}/{cons_folder}/pt_consensus_adj_subtypeB_995_refB_vb_refC.fasta.treefile")) # edit treefile to have underscores instead of dots
ml_cons_r <- ape::root(ml_cons_unr, outgroup = "refC_1986_ETH2220", resolve.root=TRUE)

ggtree(ml_cons_r) + geom_treescale(x=0.12, y=20, fontsize=4, linesize=2, offset=2, width=0.01) + geom_tiplab() +
	geom_point2(aes(subset = as.numeric(sub("/.*", "", label))>75 & !isTip))

##Annotate branch support
#Extract node label
bootsv <-data.frame(ml_cons_r[["node.label"]])
ntip <- length(ml_cons_r[["tip.label"]])
nnode <- ml_cons_r[["Nnode"]]
#bootvals <- data.frame(lapply(bootsv, function(x) as.numeric(sub("/.*", "",x))))
bootsv$node <- c((ntip + 1) : (ntip + nnode))
colnames(bootsv) <- c("label","node")

bootsv$label <- cut(as.numeric(bootsv$label), breaks = c(0,50,75,90,100),right = F,include.lowest = TRUE)
rownames(bootsv) <- bootsv$node

bootsv$label <- bootsv$label[, drop=T]

#Put the results back into the tree.
ml_cons_r[["node.label"]] <- bootsv$label
boot_pal <- c("#E6A35E", "#3A8CC3", "#7B5EA4", "#6DAA75")

# get number of sequences for each phylotype
extracted_clusters <- readRDS("rds/extracted_clusters.rds")
n_seqs_phylot <- list()
for(i in 1:length(extracted_clusters[[1,4]][[1]])) {
	if(!is.null(extracted_clusters[[1,4]][[1]][i]))
	n_seqs_phylot[[i]] <- length(extracted_clusters[[1,4]][[1]][[i]]$tip.label)
}
n_seqs_phylot_df <- as.data.frame(unlist(n_seqs_phylot))
n_seqs_phylot_df$phylotype <- rownames(n_seqs_phylot_df)
colnames(n_seqs_phylot_df) <- c("n_seqs","phylotype")
n_seqs_phylot_df$tip <- paste0("B_", n_seqs_phylot_df$phylotype)
n_seqs_phylot_df <- n_seqs_phylot_df %>% add_row(tip="VB_consensus", phylotype="VB clade", n_seqs=1)
n_seqs_phylot_df <- n_seqs_phylot_df %>% add_row(tip="refB_1983_HXB2", phylotype="Reference B", n_seqs=1)

# get CD4 model regression coeffs (ML fixed effects model)
cd4_coeffs <- readRDS(glue("{RDS_PATH}/res_fixeff_all_b_vs_bb_combined.rds"))

contre_100_cd4 <- left_join(n_seqs_phylot_df, cd4_coeffs, by="phylotype")
contre_100_cd4 <- contre_100_cd4 %>% dplyr::select(tip, n_seqs, phylotype, estimate, stderr)
contre_100_cd4$estimate[contre_100_cd4$tip == "VB_consensus"] <- -49.3 # add VB CD4 from paper, Table S1
contre_100_cd4$stderr[contre_100_cd4$tip == "VB_consensus"] <- 15 # not sure, just inputting to make sure it is included in plot

vl_bayes_model <- readRDS(glue("{RDS_PATH}/res_vl_model.rds"))[[1,4]]$df_table
vl_bayes_model$tip <- paste0("B_",vl_bayes_model$group)
vl_bayes_model <- vl_bayes_model %>% select(tip, quantile_0.5) %>% add_row(tip="VB_consensus", quantile_0.5=5.33)
contre_100_all <- left_join(contre_100_cd4, vl_bayes_model, by=c("tip"="tip"))

contre_100_all_sizes <- contre_100_all %>% dplyr::select(tip, n_seqs)
rownames(contre_100_all_sizes) <- contre_100_all_sizes$tip
contre_100_all_sizes$tip <- NULL

f3_1 <- ggtree(ml_cons_r, color="grey20") +
	geom_point2(aes(subset = !isTip, color  = label), size=1,shape=19)+
	geom_tiplab(size=1.8, font=helv, offset=0.003) +# align=TRUE, linesize=0.25
	scale_size(limits=c(1,3500), breaks=c(1,100, 250, 500, 1000, 3500), name="Sequences") +
	scale_color_manual("Bootstrap support", values = boot_pal, na.translate = F,labels = levels(bootsv$label)) + #guide=guide_legend(order=1)
	ggplot2::xlim(0, 0.21) +
	theme_tree(legend.position=c(0.8,0.5)) +
	guides(color=guide_legend(ncol=1, order=1))
ggsave(file=glue("{RESULTS_PATH}/figs/fig1_aux.svg"), plot=f3_1, dpi=600, width=8, height=10, bg="white")

leg <- theme(text=element_text(family="Helvetica"), axis.text=element_text(size=10, color="black"), axis.title=element_text(size=10, color="black", face="bold"))

f3_2 <- f3_1 %<+% contre_100_all +
	geom_treescale(x=0.13, y=15, fontsize=4, linesize=1, offset=2, width=0.01)

contre_100_all_vl <- contre_100_all %>% dplyr::select(tip,quantile_0.5)
rownames(contre_100_all_vl) <- contre_100_all_vl$tip
contre_100_all_vl$tip <- NULL

range(contre_100_all_vl$quantile_0.5,na.rm=T)
min_vl <- 4.2
max_vl <- 5.4
breaks_vl <- seq(from=4.35, to=5.35, by=0.25)
f3_3 <- gheatmap(f3_2, contre_100_all_vl[, "quantile_0.5", drop=FALSE], width = 0.1, colnames=T, custom_column_labels=c("VL"), colnames_position = "top", colnames_offset_y=0.5, family=helv, offset = -0.09) +
	scale_fill_viridis_c(option="C", values = c(0.2, 0.5, 0.6,1), direction=-1, name = "Median viral load", na.value="white", limits=c(min_vl,max_vl), breaks = breaks_vl) + coord_cartesian(clip = "off") +
	guides(fill=guide_colorbar())

contre_100_all_regr <- contre_100_all %>% dplyr::select(tip, estimate, stderr)
rownames(contre_100_all_regr) <- contre_100_all_regr$tip
contre_100_all_regr$tip <- NULL #156 rows
contre_100_all_regr <- contre_100_all_regr[contre_100_all_regr$stderr <= 30,] # remove CD4 coefficients when stderr is above 40; 148 rows

range(contre_100_all_regr$estimate,na.rm=T)
max_cd4 <- 60
min_cd4 <- -60
breaks_cd4 <- seq(from=-60, to=60, by=30)
f3_4 <- f3_3 + new_scale_fill()
f3_5 <- gheatmap(f3_4, contre_100_all_regr[, "estimate", drop=FALSE], width = 0.1, colnames=T, custom_column_labels=c("CD4"), colnames_position = "top", colnames_offset_y=0.5, family=helv, offset= -0.0725) +
	scale_fill_viridis_c(option="C", values = c(0, 0.1, 0.25, 0.5, 1), direction=1, name = "CD4 regression\ncoefficient", na.value="white", limits=c(min_cd4,max_cd4), breaks = breaks_cd4) + coord_cartesian(clip = "off") +
	leg + theme(legend.text=element_text(size=8,family=helv), legend.title=element_text(size=10,family=helv), legend.key.height=unit(.9, "cm"),legend.key.width=unit(.7, "cm"),
													legend.position = c(0.675, 0.5)) + guides(fill=guide_colorbar())

f3_5 <- f3_5 +
	geom_text(aes(x=0.061, y=108.5), label = '*', check_overlap = TRUE, color = 'grey10', size = 6, family=helv) + #40
	geom_text(aes(x=0.049, y=29.65), label = '*', check_overlap = TRUE, color = 'grey10', size = 6, family=helv) + #69
	geom_text(aes(x=0.045, y=51.5), label = '*', check_overlap = TRUE, color = 'grey10', size = 6, family=helv) #+ #133
saveRDS(f3_5, glue("{RDS_PATH}/f3_5.rds"))
### END FIGURE 1 ###

# Quantify branch lengths (compare VOIs against non-VOIs because would not be accurate to compare against consensus of backbone)
bl_cons_tr <- node.depth.edgelength(ml_cons_r)
nt <- Ntip(ml_cons_r)
nn <- Nnode(ml_cons_r)

df_tips_bl <- data.frame(tip=ml_cons_r$tip.label, bl=bl_cons_tr[1:nt])
df_tips_bl <- df_tips_bl[!(df_tips_bl$tip %in% c("refC_1986_ETH2220", "refB_1983_HXB2","VB_consensus")),]

vois_bl <- c("B_40","B_69","B_133")
df_tips_bl$status <- ifelse(df_tips_bl$tip %in% vois_bl, "VOI", "Non-VOI") 
df_tips_bl$indiv_vois_vs_non <- ifelse(df_tips_bl$tip %in% vois_bl, df_tips_bl$tip, "Non-VOI")

# wilcox.test
gghistogram(df_tips_bl, x = "bl", y = "..density..", 
												fill = "steelblue",bins = 4, add_density = TRUE)

df_tips_bl %>% group_by(status) %>% get_summary_stats(bl, type = "median_iqr")
library(rstatix)
stat_test <- df_tips_bl %>% wilcox_test(bl ~ status) %>% add_significance() # ns
print(stat_test)

mean_bl_voi_non <- df_tips_bl %>% group_by(status) %>% summarise(median_bl=median(bl), iqr=IQR(bl))
mean_bl_voi_non

# status  median_bl     iqr
# Non-VOI    0.0363 0.0115 
# VOI        0.0378 0.00880

ggplot(data=df_tips_bl, aes(x=bl, y=tip, fill=status)) + geom_bar(stat="identity")
system(glue("mkdir -p {RESULTS_PATH}/16_brlen/"))
ggsave(file=glue("{RESULTS_PATH}/16_brlen/bl1_cons.png"), dpi=600, width=10, height=12, bg="white")
ggplot(data=df_tips_bl, aes(x=status, y=bl, fill=status)) + geom_boxplot()

df_tips_bl2 <- df_tips_bl
df_tips_bl2 <- df_tips_bl2 %>% group_by(indiv_vois_vs_non) %>% summarise(bl=median(bl))
df_tips_bl2$status <- ifelse(df_tips_bl2$indiv_vois_vs_non %in% vois_bl, "VOI", "Non-VOI")
df_tips_bl2

# indiv_vois_vs_non     bl status 
# B_133             0.0322 VOI    
# B_40              0.0498 VOI (higher)
# B_69              0.0378 VOI    
# Non-VOI           0.0363 Non-VOI

ggplot(data=df_tips_bl2, aes(x=indiv_vois_vs_non, y=bl, fill=status)) + geom_bar(stat="identity")
ggsave(file=glue("{RESULTS_PATH}/16_brlen/bl2_cons.png"), dpi=600, width=8, height=8, bg="white")