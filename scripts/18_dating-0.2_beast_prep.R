libs_load <- c("ggtree","ggplot2","dplyr","ggpubr","mlesky","treedater","ape","glue","lubridate","RBeast","pbmcapply", "data.table")
invisible( lapply(libs_load, library, character.only=TRUE) )

extracted_clusters <- readRDS("rds/extracted_clusters.rds")
subtype_b_vois_ids_only <- c(40,69,133)
NCPU <- 6

size_label_adjust <- theme(axis.text=element_text(size=7), axis.title=element_text(size=8), plot.title=element_text(size=10))
seqs_folder_naive <- "data/subtype_seqs_naive"

trs <- sts <- dfs_seqs_times <- aln_files <- list()
system("mkdir -p results/18_beast/")
for(i in 1:length(subtype_b_vois_ids_only)) {
	print(i)
	
	trs[[i]] <- res_ml_dating[[i]]$td
	
	# Generate sequence_name, sample_date table for BEAST
	system("mkdir -p results/18_beast/")
	dfs_seqs_times[[i]] <- data.frame(sequence_name=trs[[i]]$tip.label)
	dfs_seqs_times[[i]]$sample_date <- trs[[i]]$sts[ match(dfs_seqs_times[[i]]$sequence_name, names(trs[[i]]$sts)) ]
	dfs_seqs_times[[i]]$sample_date <- as.Date(date_decimal(dfs_seqs_times[[i]]$sample_date))
	write.table(dfs_seqs_times[[i]], file=glue("results/18_beast/global_and_PT_B_{subtype_b_vois_ids_only[i]}_UK.tsv"), sep="\t", quote=F, row.names=F)
	
	# Extract sequences from alignment matching each cluster to run phylodynamic analysis in BEAST
	aln_files[[i]] <- read.dna( glue('{DB_PATH}/queries/aln_comb_phylotypes_global_refC_{subtype_b_vois_ids_only[i]}.fasta' ), format="fasta")
	aln_files[[i]] <- aln_files[[i]][ rownames(aln_files[[i]]) %in% dfs_seqs_times[[i]]$sequence_name, ]
	write.FASTA(aln_files[[i]], glue("results/18_beast/global_and_PT_B_{subtype_b_vois_ids_only[i]}_UK.fasta")) 
}

# Generate xmls in beast and run
# BEAST analysis: 
# GTR+G4 empirical (previously HKY+G4 estimated) -> URLC -> Bayesian skygrid (grid points and time at last point transition: max(trs[[i]]$sts) - trs[[i]]$timeOfMRCA
# Clock rate prior: previously 0.00115, now trs[[i]]$mean.rate (0.001, 0.001, 0.001, 0.001, )
# Skygrid parameters: PT40: 44, PT69: 44, PT133: 46
# Previously: PT40, 69 and 133 converged well with one chain of 200M states.
# Command: beast global_and_PT_B_{i}_UK_chain1.xml > chain1.log 2>&1 &

# Read posterior trees and extract tmrca of PTs
# 1) https://github.com/beast-dev/RBeast/blob/master/R/read_beast2_trees.R

#devtools::install_github("beast-dev/RBeast")
vois <- c("40_UK","69_UK","133_UK")
subtype_b_vois_ids_only <- c(40,69,133)

post_trees <- pbmcapply::pbmclapply(1:length(subtype_b_vois_ids_only), function(i) { read_beast2_trees(glue("results/18_beast/global_and_PT_B_{vois[i]}.trees")) }, mc.cores = NCPU)
#saveRDS(post_trees, "rds/post_trees.rds")
post_trees <- readRDS("rds/post_trees.rds")

extract_tmrcas_beast <- function(i, md, voi_id) {
	mrcas_beast <- rh <- mrst <- mrca_node <- c() # tmrca tmrcas_beast
	nde <- stimes <- shs <- nhs <- tmrca <- tmrcas_beast <- list()
	for(x in 1:length(post_trees[[i]])) {
		#print(x)
		nde[[x]] <- ape::node.depth.edgelength(post_trees[[i]][[x]])
		rh[x] <- max(nde[[x]][1:ape::Ntip(post_trees[[i]][[x]])]) # root heights
		stimes[[x]] <- nde[[x]][1:ape::Ntip(post_trees[[i]][[x]])]
		shs[[x]] <- rh[x] - stimes[[x]] # time to most recent sample
		nhs[[x]] <- rh[x] - nde[[x]] # node heights 	nhs[[x]] = nodeheights <- rh[x] - nde[[x]]
		md <- md[md$Accession %in% post_trees[[i]][[x]]$tip.label,]
		mrst[x] <- max(md$Sampling.Year) # most recent sample time
		tmrca[[x]] <- mrst[x] - nhs[[x]]
		#mrca_node[x] <- getMRCA(phy=post_trees[[i]][[x]], tip=intersect(md$Accession, post_trees[[i]][[x]]$tip.label) )
		mrcas_beast[x] <- getMRCA(phy=post_trees[[i]][[x]], tip=md$Accession[md$phylotype == glue("PT.B.{voi_id}.UK")] )
		#tmrcas_beast[x] <- ape::node.depth.edgelength(post_trees[[i]][[x]])[ mrcas_beast[x] ] + tmrca[ mrca_node[x] ]
		tmrcas_beast[[x]] <- tmrca[[x]][ mrcas_beast[x] ]
	}
	return(tmrcas_beast)
}

combined_seqs_md_adj <- readRDS("rds/combined_seqs_md_adj.rds")
pt69_beast <- read.dna("results/18_beast/global_and_PT_B_69_UK.fasta", format="fasta")
combined_seqs_md_adj[[2]] <- combined_seqs_md_adj[[2]][ combined_seqs_md_adj[[2]]$Accession %in% rownames(pt69_beast),  ]
subtype_b_vois_ids_only <- c(40, 69, 133)
subtype_b_vois <- c("PT.B.40.UK","PT.B.69.UK","PT.B.133.UK")

tmrca_post_beast <- tmrca_post_beast_res <- tmrca_post_beast_res2 <- list()
for(i in 1:length(subtype_b_vois_ids_only)) {
	print(glue("======{subtype_b_vois_ids_only[i]}======"))
	tmrca_post_beast[[i]] <- extract_tmrcas_beast(i, combined_seqs_md_adj[[i]], subtype_b_vois_ids_only[i]) # careful this (combined_seqs_md_adj) is in order from PTs 5 to 137 (1 to 5)
	tmrca_post_beast_res[[i]] <- unlist(tmrca_post_beast[[i]])
	tmrca_post_beast_res2[[i]] <- data.frame(repl=seq(1,length(tmrca_post_beast_res[[i]])), tmrca=tmrca_post_beast_res[[i]], phylotype=subtype_b_vois_ids_only[i], method="BEAST")
}

df_tmrcas_beast <- rbindlist(tmrca_post_beast_res2)
df_tmrcas_beast$phylotype <- paste0("PT.B.",df_tmrcas_beast$phylotype,".UK")
df_tmrcas_beast$phylotype <- factor(df_tmrcas_beast$phylotype, levels = subtype_b_vois)

df_tmrcas <- readRDS("rds/df_tmrcas.rds")

df_tmrcas_combined <- rbind(df_tmrcas, df_tmrcas_beast)
saveRDS(df_tmrcas_combined, glue("{RDS_PATH}/df_tmrcas_combined.rds"))

helv <- "Helvetica"
RESULTS_PATH="results"
method_pal <- c("BEAST"="#D55E00","Treedater parboot"="#377EB8")
ggplot(df_tmrcas_combined, aes(x=tmrca, y=phylotype, group = interaction(phylotype,method))) +
	geom_violin(width=1, position=position_dodge(1), alpha=0.9, bw=1.5, aes(fill=method)) +
	geom_boxplot(width=0.1, position=position_dodge(1), fill="white", alpha=1, outlier.shape=NA) + #fill="white"
	scale_fill_manual(values=method_pal, name="Estimation\nmethod") + 
	labs(x="Year", y="Phylotype") + scale_x_continuous(expand = c(0,0), limits=c(1975,2005), breaks=seq(1975,2005,by=5)) + #+ scale_x_discrete(expand=c(0,0)) +
	theme_classic() + theme(legend.text=element_text(size=10,family=helv), legend.title=element_text(size=11,family=helv)) #+ leg 
ggsave(file=glue("{RESULTS_PATH}/figs/figS17.eps"), device=cairo_ps, dpi=600, width=10, height=12, bg="white") 
ggsave(file=glue("{RESULTS_PATH}/figs/figS17.jpg"), dpi=600, width=10, height=12, bg="white") 

# Summary stats
df_tmrcas_combined_summ <- df_tmrcas_combined %>% group_by(method, phylotype) %>% summarise(min=min(tmrca), max=max(tmrca), median=date_decimal(median(tmrca)), lower_ci=date_decimal(quantile(tmrca, probs = 0.025 )), upper_ci=date_decimal(quantile(tmrca, probs = 0.975 ))) #iqr=stats::IQR(trmca)
View(df_tmrcas_combined_summ)
saveRDS(df_tmrcas_combined_summ, glue("{RDS_PATH}/df_tmrcas_combined_summ.rds"))