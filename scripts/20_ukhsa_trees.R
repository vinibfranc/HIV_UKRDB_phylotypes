libs_load <- c("glue", "ape", "ggtree", "ggplot2", "ggpubr")
invisible( lapply(libs_load, library, character.only=TRUE) )

tree_files <- list.files(path="results/20_phylotype_ukhsa_trees", pattern = "\\.tree$", full.names = T)
voi_ids <- c("PT.B.133.UK", "PT.B.40.UK", "PT.B.69.UK") #"PT.B.137.UK","PT.B.5.UK"

read_and_plot_ukhsa_trees <- function(tf, voi_id) {
	tr <- read.tree(tf)
	df_tr <- data.frame(name=tr$tip.label)
	df_tr$dataset <- ifelse(grepl("t.", df_tr$name), yes=glue("UKRDB VOI {voi_id}"), no="INITIO and COMPARE studies")
	df_tr$dataset <- ifelse(grepl("B.", df_tr$name), yes="Subtype reference sequences", no=df_tr$dataset)
	
	pl <- ggtree(tr, layout = "rectangular", color="grey20") %<+% df_tr + geom_tippoint(aes(fill=dataset), size=1.5, stroke=0.5, color="black", shape=21, alpha=1) + 
		scale_fill_manual("Dataset", values=phylot_pal_wgs) + 
		geom_treescale(x=0.15, y=30, fontsize=4, linesize=1, offset=10, family=helv, width=0.05) + 
		theme(text=element_text(family=helv)) + expand_limits(y = -25)
	pl
}

phylot_pal_wgs <- c("UKRDB VOI PT.B.40.UK"="#56B4E9","UKRDB VOI PT.B.69.UK"="#D55E00","UKRDB VOI PT.B.133.UK"="#009E73", "INITIO and COMPARE studies"="#CC79A7", "Subtype reference sequences"="black") #92A0AD

ukhsa_trees <- list()
for(i in 1:length(tree_files)) {
	ukhsa_trees[[i]] <- read_and_plot_ukhsa_trees(tree_files[i], voi_ids[i])
}

library(gridExtra)

legend2 <- ggplotGrob(ukhsa_trees[[2]]); legend2_grob <- legend2$grobs[[which(legend2$layout$name == "guide-box")]]
legend3 <- ggplotGrob(ukhsa_trees[[3]]); legend3_grob <- legend3$grobs[[which(legend3$layout$name == "guide-box")]]
legend4 <- ggplotGrob(ukhsa_trees[[1]]); legend4_grob <- legend4$grobs[[which(legend4$layout$name == "guide-box")]]
combined_legend <- grid.arrange(legend2_grob, legend3_grob, legend4_grob, nrow = 1)

s5 <- grid.arrange(
	combined_legend,
	nrow = 2,
	heights = c(1, 7),
	arrangeGrob(ukhsa_trees[[2]] + theme(legend.position = "none"),
													ukhsa_trees[[3]] + theme(legend.position = "none"),
													ukhsa_trees[[1]] + theme(legend.position = "none"),
													nrow = 1)
)
ggsave(s5, file=glue("{RESULTS_PATH}/figs/figS5.eps"), device=cairo_ps, dpi=600, width=8, height=10, bg="white")
ggsave(s5, file=glue("{RESULTS_PATH}/figs/figS5.jpg"), dpi=600, width=8, height=10, bg="white")