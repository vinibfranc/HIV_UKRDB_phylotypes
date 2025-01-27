libs_load <- c("glue", "ape", "ggtree", "ggplot2","magrittr", "ggnewscale", "dplyr", "tidyverse", "data.table", "lubridate", "ggpubr")
invisible( lapply(libs_load, library, character.only=TRUE) )

RDS_PATH <- "rds"

timetree_b_adj <- readRDS("rds/timetree_b_strict_adj.rds") 
class(timetree_b_adj) <- "phylo"

### BEGIN FIGURE S1 ###
# Include demographic and phylotype info
demog_md_subtype_match <- readRDS(glue("{RDS_PATH}/demog_md_subtype_match.rds"))
demog_md_subtype_match <- demog_md_subtype_match[ (demog_md_subtype_match$status=="Naïve") & 
																																																			(demog_md_subtype_match$rega3subtype %in% c("A (A1)","CRF 02_AG","C","B")) & 
																																																			(demog_md_subtype_match$exposureid != "Not known"), ] #53382
df_sb_tree <- demog_md_subtype_match
df_sb_tree$testindex <- paste0("t.",df_sb_tree$testindex)
df_sb_tree <- df_sb_tree[df_sb_tree$testindex %in% timetree_b_adj$tip.label,]

extracted_clusters <- readRDS(glue("{RDS_PATH}/extracted_clusters.rds"))
backbone_cl_control <- readRDS(glue("{RDS_PATH}/backbone_cl_control.rds"))

pt_ids <- list()
for(i in 1:length(extracted_clusters[[1,4]][[1]])) {
	if(!(i %in% as.integer(backbone_cl_control[[1,4]])) ) {
		if(!is.null(extracted_clusters[[1,4]][[1]][[i]]$tip.label)) {
			pt_ids[[i]] <- df_sb_tree[ df_sb_tree$testindex %in% extracted_clusters[[1,4]][[1]][[i]]$tip.label, ]
			pt_ids[[i]]$phylotype <- i
			pt_ids[[i]] <- pt_ids[[i]] %>% dplyr::select(testindex,phylotype)
		}
	}
}

pt_ids_all <- rbindlist(pt_ids)

df_sb_tree_pt <- left_join(df_sb_tree, pt_ids_all, by="testindex")

df_sb_tree_pt$date <- as.Date(date_decimal(df_sb_tree_pt$dbsample_date))
df_sb_tree_pt <- df_sb_tree_pt %>% dplyr::select(testindex, date, phylotype)

df_sb_tree_pt$phylotype <- paste0("PT.B.",df_sb_tree_pt$phylotype,".UK")
subtype_b_vois <- c("PT.B.40.UK","PT.B.69.UK","PT.B.133.UK")
df_sb_tree_pt$phylotype_annot <- ifelse(!(df_sb_tree_pt$phylotype %in% subtype_b_vois), NA_character_, df_sb_tree_pt$phylotype) #"Non-VOIs"

voi_pal <- c("PT.B.40.UK"="#56B4E9", "PT.B.69.UK"="#D55E00", "PT.B.133.UK"="#009E73")

# Subsample tree according to phylotype sizes (divide by 10) for visualisation purposes
df_sb_tree_pt_list <- df_sb_tree_pt_list_subsamp <- list()
dv <- 5
for(i in 1:154) {
	if(!(i %in% as.integer(backbone_cl_control[[1,4]])) ) {
		if(!is.null(extracted_clusters[[1,4]][[1]][[i]]$tip.label)) {
			df_sb_tree_pt_list[[i]] <- df_sb_tree_pt[df_sb_tree_pt$phylotype == paste0("PT.B.",i,".UK"),]
			nr <- nrow(df_sb_tree_pt_list[[i]])
			df_sb_tree_pt_list_subsamp[[i]] <- df_sb_tree_pt_list[[i]][ sample(nr,nr/dv, replace=FALSE), ]
		}
	}
}
df_sb_tree_pt_subsamp_comb <- rbindlist(df_sb_tree_pt_list_subsamp)
df_sb_tree_pt_subsamp_comb$phylotype_annot <- factor(df_sb_tree_pt_subsamp_comb$phylotype_annot, subtype_b_vois) #"PT.B.50.UK"
df_sb_tree_pt_subsamp_comb <- df_sb_tree_pt_subsamp_comb %>% add_row(testindex="t.Ref.C.ET.86.ETH2220.U46016",date=as.Date("1986-06-01"), phylotype=NA, phylotype_annot=NA)

timetree_b_subsamp <- keep.tip(timetree_b_adj, intersect(timetree_b_adj$tip.label, df_sb_tree_pt_subsamp_comb$testindex))

# ML tree
ml_tree_b <- readRDS("rds/tree_b_adj2.rds")
ml_tree_b_subsamp <- keep.tip(ml_tree_b, intersect( ml_tree_b$tip.label, df_sb_tree_pt_subsamp_comb$testindex))

df_sb_tree_pt_subsamp_comb <- as.data.frame(df_sb_tree_pt_subsamp_comb)

helv <- "Helvetica"

ml_tr_f1 <- ggtree(ml_tree_b_subsamp, color="grey20", size=0.1) + geom_treescale(x=0.12, y=80, fontsize=3, linesize=1, offset=10)
ml_tr_f2 <- ml_tr_f1 %<+% df_sb_tree_pt_subsamp_comb +
	geom_tippoint(aes(subset=(phylotype_annot=="PT.B.40.UK" & !is.na(phylotype_annot))),fill="#56B4E9",size=2.5, stroke=0.15, shape=21) +
	geom_tippoint(aes(subset=(phylotype_annot=="PT.B.69.UK" & !is.na(phylotype_annot))),fill="#D55E00",size=2.5, stroke=0.15, shape=21) +
	geom_tippoint(aes(subset=(phylotype_annot=="PT.B.133.UK" & !is.na(phylotype_annot))),fill="#009E73",size=2.5, stroke=0.15, shape=21) +
	theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + expand_limits(y = -25)

ml_tr_f3 <- ml_tr_f1 %<+% df_sb_tree_pt_subsamp_comb + geom_tippoint(aes(color = phylotype_annot, subset = !is.na(phylotype_annot)), size = 1.5, alpha = 1, na.rm = TRUE) +
	scale_color_manual(values=voi_pal, name="VOIs", guide=guide_legend(order=1, ncol=4)) + theme(legend.text=element_text(family = helv, size=8), legend.title=element_text(family = helv, size=10)) 
leg_f1 <- get_legend(ml_tr_f3)

# Figure S6: main tree subsampled according to phylotype size
ggarrange(ml_tr_f2, nrow=1, ncol=1, legend.grob = leg_f1, legend="top", font.label=list(family="Helvetica", color="black",size=10))
ggsave(file=glue("{RESULTS_PATH}/figs/figS6.jpg"), dpi=900, width=8, height=10, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS6.eps"), device = cairo_ps, family = "Helvetica", dpi=900, width=8, height=10, bg="white")
### END FIGURE S6 ###

# Timetree (does not add much info I think)
class(timetree_b_subsamp) <- "phylo"
tr_f1 <- ggtree(timetree_b_subsamp, mrsd=max(df_sb_tree_pt_subsamp_comb$date), as.Date=TRUE, color="grey40", size=0.1) + theme_tree2()
tr_f1 %<+% df_sb_tree_pt_subsamp_comb + geom_tippoint(aes(color = phylotype_annot, subset = !is.na(phylotype_annot)), size = 1.5, alpha = 1, na.rm = TRUE) +
	scale_color_manual(values=voi_pal, name="VOIs") + scale_x_date(limits=c(as.Date("1970-01-01"),as.Date("2020-01-01")), date_labels = "%Y",breaks="5 years") +
	xlab("Time") +
	theme(axis.text.x=element_text(family = helv, size=10, color="black"), legend.text=element_text(family = helv, size=8), legend.title=element_text(family = helv, size=10))

# test if branches leading to phylotype are bigger for VOI vs non-VOIs
bl_full_tr <- node.depth.edgelength(ml_tree_b)
nt2 <- ape::Ntip(ml_tree_b)
nn2 <- Nnode(ml_tree_b)

# get MRCA nodes leading to phylotypes
mrca_nodes <- c()
for(i in 1:154) {
	if(!(i %in% as.integer(backbone_cl_control[[1,4]])) ) {
		if(!is.null(extracted_clusters[[1,4]][[1]][[i]]$tip.label)) {
			mrca_nodes[i] <- getMRCA(ml_tree_b, extracted_clusters[[1,4]][[1]][[i]]$tip.label)
		}
	}
}
#mrca_nodes <- mrca_nodes[ !is.na(mrca_nodes) ]

bl_nodes_phylotypes <- c()
for(i in 1:154) { 
	bl_nodes_phylotypes[i] <- bl_full_tr[ mrca_nodes[i] ] 
} 

df_nodes_bl <- data.frame(node=mrca_nodes, bl=bl_nodes_phylotypes)
df_nodes_bl$phylotype <- 1:nrow(df_nodes_bl)

vois_bl2 <- c("40","69","133")
df_nodes_bl$status <- ifelse(df_nodes_bl$phylotype %in% vois_bl2, "VOI", "Non-VOI")
df_nodes_bl$indiv_vois_vs_non <- ifelse(df_nodes_bl$phylotype %in% vois_bl2, df_nodes_bl$phylotype, "Non-VOI")
df_nodes_bl <- df_nodes_bl[!is.na(df_nodes_bl$node),]

mean_bl_voi_non2 <- df_nodes_bl %>% group_by(status) %>% summarise(mean_bl=mean(bl), sd=sd(bl))
mean_bl_voi_non2

# status  mean_bl    sd
# 1 Non-VOI  0.0515 0.0147
# 2 VOI      0.0505 0.0111

median_bl_voi_non2 <- df_nodes_bl %>% group_by(status) %>% summarise(median_bl=median(bl), iqr=IQR(bl))
median_bl_voi_non2

# status  median_bl    iqr
# 1 Non-VOI    0.0508 0.0184
# 2 VOI        0.0445 0.00987

library(tidyverse)
library(rstatix)
library(ggpubr)

# wilcox.test
gghistogram(df_nodes_bl, x = "bl", y = "..density..", 
												fill = "steelblue",bins = 4, add_density = TRUE)

df_nodes_bl %>% group_by(status) %>% get_summary_stats(bl, type = "median_iqr")
stat_test2 <- df_nodes_bl %>% wilcox_test(bl ~ status) %>% add_significance() # ns
stat_test2 #ns

ggplot(data=df_nodes_bl, aes(x=bl, y=as.factor(phylotype), fill=status)) + geom_bar(stat="identity")
system(glue("mkdir -p {RESULTS_PATH}/14_brlen/"))
ggsave(file=glue("{RESULTS_PATH}/14_brlen/bl1_full.png"), dpi=600, width=10, height=12, bg="white")
ggplot(data=df_nodes_bl, aes(x=status, y=bl, fill=status)) + geom_boxplot()

df_nodes_bl2 <- df_nodes_bl
df_nodes_bl2 <- df_nodes_bl2 %>% group_by(indiv_vois_vs_non) %>% summarise(bl=median(bl))
df_nodes_bl2$status <- ifelse(df_nodes_bl2$indiv_vois_vs_non %in% vois_bl2, "VOI", "Non-VOI")
df_nodes_bl2

# indiv_vois_vs_non     bl status 
# 133               0.0436 VOI    
# 40                0.0633 VOI (slightly higher)   
# 69                0.0445 VOI    
# Non-VOI           0.0508 Non-VOI

ggplot(data=df_nodes_bl2, aes(x=indiv_vois_vs_non, y=bl, fill=status)) + geom_bar(stat="identity")
ggsave(file=glue("{RESULTS_PATH}/14_brlen/bl2_full.png"), dpi=600, width=8, height=8, bg="white")