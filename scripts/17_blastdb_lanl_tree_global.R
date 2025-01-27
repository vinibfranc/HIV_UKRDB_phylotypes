#install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
libs_load <- c( "ape", "glue", "rBLAST", "pbmcapply", "dplyr","data.table", "ggtree", "ggplot2", "viridis", "ggnewscale", "patchwork", "caper")
invisible( lapply(libs_load, library, character.only=TRUE) )

# rBLAST version 0.99.2 used here, be aware that newest 0.99.4 used in other analysis changed columns (to qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) 

# Download all pol sequences from LANL (B, Pol CDS, min_length=950)
DB_PATH="data/blast_db_lanl"
REF_SEQS_PATH="data/refs/"
FASTA_FILE="all_global_b"  #"B_pol_lanl_OLD"
seqs_folder_naive <- "data/subtype_seqs_naive"
NCPU <- 6

# Remove gaps from fasta files and join
system(glue("mkdir -p {DB_PATH}/fasta_without_gaps/"))
system(glue("for f in {DB_PATH}/alns/*.fasta; do sed 's/-//g' $f > {DB_PATH}/fasta_without_gaps/$(basename $f); done"))
system(glue("cat {DB_PATH}/fasta_without_gaps/*.fasta > {DB_PATH}/all_global_b.fasta"))
# Delete manually 290 sequences with repeated IDs indicated by AliView

# Join metadata files 
system(glue("awk -F '\t' '(NR == 1) || (FNR > 1)' {DB_PATH}/mds/*.txt > {DB_PATH}/metadata_all_global_b.txt"))

# Build BLAST database
system(glue("makeblastdb -in {DB_PATH}/{FASTA_FILE}.fasta -parse_seqids -blastdb_version 5 -dbtype nucl"))

# Create fasta for phylotypes of interest
fasta_b <- read.dna(glue("{seqs_folder_naive}/B_UK_final_aln_outliers_removed_final.fasta"), format="fasta")
extracted_clusters <- readRDS("rds/extracted_clusters.rds")

create_fasta <- function(phylotype) {
	phylot <- fasta_b[rownames(fasta_b) %in% extracted_clusters[[1,4]][[1]][[phylotype]]$tip.label,]
	system(glue("mkdir -p {DB_PATH}/queries/"))
	write.FASTA(phylot, glue("{DB_PATH}/queries/phylotype_{phylotype}.fasta"))
}

subtype_b_vois <- c(40,69,133)
for(i in subtype_b_vois) {
	create_fasta(i)
}

# Replace ? by Ns and remove gaps
for(i in subtype_b_vois) {
	system(glue("sed 's/?/N/g' {DB_PATH}/queries/phylotype_{i}.fasta > {DB_PATH}/queries/phylotype_{i}_adj.fasta"))
	system(glue("sed 's/-//g' {DB_PATH}/queries/phylotype_{i}_adj.fasta > {DB_PATH}/queries/phylotype_{i}_adj2.fasta"))
}

# load db
bl <- blast(db=glue("{DB_PATH}/{FASTA_FILE}.fasta"))

query_seqs <- function(query_phylotype_id, percent_id_par=95) {
	seqs <- readDNAStringSet(glue("{DB_PATH}/queries/phylotype_{query_phylotype_id}_adj2.fasta"))
	print(glue("ntips: {length(seqs)}"))
	
	# Matches against each sequence
	# https://github.com/mhahsler/rBLAST
	cl = pbmcapply::pbmclapply( 1:length(seqs), function(i){
		predict(bl, seqs[i,], BLAST_args = glue("-num_threads {NCPU} -perc_identity {percent_id_par}"))
	}, mc.cores = NCPU)
	
	# Unique matches within the phylotype
	all_matches <- rbindlist(cl)
	
	# Order by percent identity
	all_matches$Perc.Ident <- as.numeric(all_matches$Perc.Ident)
	all_matches <- all_matches[order(all_matches$Perc.Ident, decreasing=T),]
	
	# Get all sequences (queries) matching a particular subject
	all_matches <- all_matches %>% group_by(SubjectID) %>% mutate(seqs_matching = paste0(na.omit(QueryID), collapse=",")) %>% ungroup()
	# Get unique close relatives
	all_matches_uniq <- all_matches[!duplicated(all_matches$SubjectID), ]
	#View(all_matches_uniq)
	print(glue("Number of unique closely relatives > {percent_id_par}% percent identity threshold"))
	print(nrow(all_matches_uniq))
	
	metadata_global_seqs <- read.table(glue("{DB_PATH}/metadata_all_global_b.txt"), header=T, sep="\t")
	metadata_global_seqs <- metadata_global_seqs[!duplicated(metadata_global_seqs$Accession), ] # remove duplicated accessions
	# Join metadata info with closely relatives
	all_matches_uniq_merged <- inner_join(all_matches_uniq, metadata_global_seqs, by=c("SubjectID"="Accession"))
	
	# Retrieve only non-UK sequences
	all_matches_uniq_nonuk <- all_matches_uniq_merged[all_matches_uniq_merged$Country != "UNITED KINGDOM",]
	print("Number of NON-UK unique closely relatives")
	print(nrow(all_matches_uniq_nonuk))
	
	return(list(all_matches_uniq_nonuk, all_matches_uniq_merged))
}

id_thr <- c(90, 95, 96, 97, 98, 99, 100)

blast_id_vars <- matrix(list(), nrow=length(subtype_b_vois), ncol=length(id_thr))
rownames(blast_id_vars) <- subtype_b_vois
colnames(blast_id_vars) <- id_thr
for(i in seq_along(subtype_b_vois)) {
	for(j in seq_along(id_thr)) {
		print(glue("PT={subtype_b_vois[i]}; ID={id_thr[j]}"))
		blast_id_vars[[i,j]] <- query_seqs(subtype_b_vois[i], id_thr[j])
	}
}
saveRDS(blast_id_vars, "rds/blast_id_vars.rds")

global_db <- read.dna(glue("{DB_PATH}/{FASTA_FILE}.fasta"), format="fasta")
md_global_db <- read.table(glue("{DB_PATH}/metadata_all_global_b.txt"), sep="\t", header=T)
md_global_db <- md_global_db[!duplicated(md_global_db$Accession), ]

output_iqtree_folder_blast <- glue("{DB_PATH}/queries/iqtree")

### BEGIN FIGURES 4B AND S13 ####
non_uk_matches_vois <- list() 
for(i in seq_along(subtype_b_vois))
	non_uk_matches_vois[[i]] <- blast_id_vars[[i,2]][[1]] # get non-UK matches from 95% perc id threshold

non_uk_matches_vois_uniq <- non_uk_matches_vois_sampled <- list()
for(i in 1:length(non_uk_matches_vois)) {
	# Exclude duplicates of Perc.Ident and S.start
	non_uk_matches_vois_uniq[[i]] <- non_uk_matches_vois[[i]] %>% distinct(Perc.Ident, S.start, .keep_all=TRUE)
	# Randomly sample 250 matches to add to the phylotype's tree
	non_uk_matches_vois_sampled[[i]] <- non_uk_matches_vois_uniq[[i]][sample(nrow(non_uk_matches_vois_uniq[[i]]), 250, replace=FALSE),]
}
hist(non_uk_matches_vois_sampled[[1]]$Perc.Ident) # inspect from 1 to 5

saveRDS(non_uk_matches_vois_sampled, "rds/non_uk_matches_vois_sampled.rds")

non_uk_matches_vois_sampled <- readRDS("rds/non_uk_matches_vois_sampled.rds")
# did not save seed generated

retrieve_cl_related_seqs_md <- function(phylotypes_list, qr_perc) {
	close_rel_fasta <- close_rel_md <- phylotype_md <- merged_dfs <- list() 
	for(i in 1:length(phylotypes_list)) {
		close_rel_fasta[[i]] <- global_db[names(global_db) %in% qr_perc$SubjectID] 
		print(length(names(close_rel_fasta[[i]])))
		close_rel_md[[i]] <- md_global_db[md_global_db$Accession %in% names(close_rel_fasta[[i]]),]
		close_rel_md[[i]] <- subset(close_rel_md[[i]], select=c("Accession", "Sampling.Year", "Country", "Georegion"))
		close_rel_md[[i]] <- close_rel_md[[i]][!is.na(close_rel_md[[i]]$Sampling.Year),]
		close_rel_md[[i]]$Sampling.Year <- paste0(close_rel_md[[i]]$Sampling.Year,".500") # Add .5 to define midpoint time for integer years
		close_rel_md[[i]]$Sampling.Year <- as.numeric(close_rel_md[[i]]$Sampling.Year)
		close_rel_md[[i]]$phylotype <- "Non-UK"
		close_rel_fasta[[i]] <- close_rel_fasta[[i]][names(close_rel_fasta[[i]]) %in% close_rel_md[[i]]$Accession]

		write.FASTA(close_rel_fasta[[i]], file=glue("{DB_PATH}/queries/cl_relatives_phylotype_{phylotypes_list[i]}.fasta"))
		
		phylotype_md[[i]] <- d_demog[d_demog$phylotype ==  phylotypes_list[i],]
		phylotype_md[[i]]$phylotype <- paste0("PT.B.",phylotype_md[[i]]$phylotype,".UK")
		phylotype_md[[i]]$Georegion <- "Europe"
		phylotype_md[[i]]$Country <- "UNITED KINGDOM"
		phylotype_md[[i]] <- subset(phylotype_md[[i]], select=c("testindex", "dbsample_date", "Country", "Georegion", "phylotype"))
		colnames(phylotype_md[[i]]) <- c("Accession", "Sampling.Year", "Country", "Georegion", "phylotype")
		merged_dfs[[i]] <- rbind(phylotype_md[[i]], close_rel_md[[i]])
		
		merged_dfs[[i]] <- as.data.frame(merged_dfs[[i]])
		# careful below!
		system(glue("cat {DB_PATH}/queries/phylotype_{phylotypes_list[i]}_adj2.fasta {DB_PATH}/queries/cl_relatives_phylotype_{phylotypes_list[i]}.fasta > {DB_PATH}/queries/comb_phylotypes_global_{phylotypes_list[i]}.fasta"))
		
	}
	merged_dfs
}

subtype_b_vois_ids_only <- c(40, 69, 133)

# Subsampled 95% identity close relatives with phylotype seqs
combined_seqs_md <- list()
for(i in 1:length(subtype_b_vois)) {
	combined_seqs_md[[i]] <- retrieve_cl_related_seqs_md(subtype_b_vois_ids_only[i], non_uk_matches_vois_sampled[[i]])
}

saveRDS(combined_seqs_md, glue("{RDS_PATH}/combined_seqs_md.rds"))

# align all VB clade samples against C
## system(glue("mafft --6merpair --thread 4 --addfragments {REF_SEQS_PATH}/vb_clade.fasta {REF_SEQS_PATH}/HIV1_REF_2020_pol_DNA_SUBTYPE_C.fasta > {DB_PATH}/queries/aln_all_VB_refC.fasta"))
# align VOIs + closely relatives against VB (not anymore) and refC
for(i in subtype_b_vois) {
	system(glue("mafft --6merpair --thread 4 --addfragments {DB_PATH}/queries/comb_phylotypes_global_{i}.fasta {REF_SEQS_PATH}/HIV1_REF_2020_pol_DNA_SUBTYPE_C.fasta > {DB_PATH}/queries/aln_comb_phylotypes_global_refC_{i}.fasta"))
}

# Trim ends of alns, add Ns at ends and remove gap-only columns
# PT40: removing gaps from only one seq and keeping if >1 seq has it

# Run IQTREE
for(i in subtype_b_vois) {
	system(glue("{iqtree_bin} -s {DB_PATH}/queries/aln_comb_phylotypes_global_refC_{i}.fasta -m GTR+R -nt AUTO -B 1000 -nm 5000")) #for B133 setting -nm 10000 -bcor 0.985 because BTnot converging quickly
}

system(glue("mkdir -p {output_iqtree_folder_blast}"))
system(glue("mv {DB_PATH}/queries/*.bionj {DB_PATH}/queries/*.gz {DB_PATH}/queries/*.iqtree {DB_PATH}/queries/*.log {DB_PATH}/queries/*.mldist {DB_PATH}/queries/*.treefile {DB_PATH}/queries/*.nex {DB_PATH}/queries/*.contree {output_iqtree_folder_blast}"))

# Palettes
phylot_pal <- c("PT.B.40.UK"="#56B4E9", "PT.B.69.UK"="#D55E00", "PT.B.133.UK"="#009E73", "Non-UK"="#999999")
continent_pal <- c("Africa" = "#E69F00","Asia" = "#377EB8","Europe" = "#73D055","North America" = "#F0E442","Oceania" = "#9467BD","South America" = "#CC79A7")

combined_seqs_md_adj <- unlist(combined_seqs_md, recursive=FALSE)

# Add "t." to UK seqs
for(i in 1:length(combined_seqs_md_adj)) {
	combined_seqs_md_adj[[i]]$Accession <- ifelse(combined_seqs_md_adj[[i]]$phylotype != "Non-UK", yes=paste0("t.",combined_seqs_md_adj[[i]]$Accession), no=as.character(combined_seqs_md_adj[[i]]$Accession))
}

saveRDS(combined_seqs_md_adj, "rds/combined_seqs_md_adj.rds")

combined_seqs_md_adj <- readRDS( glue("{RDS_PATH}/combined_seqs_md_adj.rds") )

sts_vois <- sts_vois_order <- trees_global_vois_r <- dfs_gheat <- dfs_gheat_stack <- p1 <- pl2 <- p3 <- pst <- list() #trees_global_vois
for(i in 1:length(subtype_b_vois_ids_only)) {
	rf <- read.tree(glue("{DB_PATH}/queries/iqtree/aln_comb_phylotypes_global_refC_{subtype_b_vois_ids_only[i]}.fasta.treefile"))
	combined_seqs_md_adj[[i]] <- combined_seqs_md_adj[[i]] %>% add_row(Accession="Ref.C.ET.86.ETH2220.U46016", Sampling.Year=1986.5, Country="Ethiopia", Georegion="Africa", phylotype="C reference")
	combined_seqs_md_adj[[i]] <- combined_seqs_md_adj[[i]][combined_seqs_md_adj[[i]]$Accession %in% rf$tip.label,] #names(rf)
	rf <- keep.tip(rf, intersect(combined_seqs_md_adj[[i]]$Accession, rf$tip.label) )
	sts_vois[[i]] <- combined_seqs_md_adj[[i]]$Sampling.Year
	names(sts_vois[[i]]) <- combined_seqs_md_adj[[i]]$Accession
	# Make sure sample times are in same order as tip in the tree
	sts_vois_order[[i]] <- sts_vois[[i]][order(match(names(sts_vois[[i]]),rf$tip.label))] #trees_global_vois[[i]]$tip.label
	# Root using outgroup
	trees_global_vois_r[[i]] <- root(rf, outgroup = "Ref.C.ET.86.ETH2220.U46016") #trees_global_vois[[i]]
	
	dfs_gheat[[i]] <- combined_seqs_md_adj[[i]][ !duplicated(combined_seqs_md_adj[[i]]$Accession), ]
	dfs_gheat[[i]] <- dfs_gheat[[i]] %>% mutate(
		Georegion = dplyr::case_when(Georegion == "Africa" ~ "Africa", Georegion == "Asia" ~ "Asia", Georegion == "Europe"  ~ "Europe", Georegion == "Caribbean"  ~ "North America", 
																															Georegion == "Central America"  ~ "North America", Georegion == "North America"  ~ "North America",Georegion == "Oceania"  ~ "Oceania", 
																															Georegion == "South America"  ~ "South America"),
		Georegion = factor(Georegion,level = c("Africa","Asia","Europe","North America","Oceania","South America")))
	# Considering Russia as Europe
	dfs_gheat[[i]]$Georegion <- ifelse(dfs_gheat[[i]]$Country == "RUSSIAN FEDERATION", yes="Europe", no=as.character(dfs_gheat[[i]]$Georegion))
	dfs_gheat_stack[[i]] <- dfs_gheat[[i]]
	dfs_gheat_stack[[i]]$year <- as.integer(dfs_gheat_stack[[i]]$Sampling.Year)
	dfs_gheat_stack[[i]] <- dfs_gheat_stack[[i]] %>% group_by(Georegion, year) %>% summarise(n=n())
	
	# Distribution of sequences over time across continents (global regions)
	pst[[i]] <- ggplot(dfs_gheat_stack[[i]], aes(fill=Georegion, y=n, x=year)) + geom_bar(position="stack", stat="identity") + labs(x="Year", y="Count") + 
		scale_fill_manual(values=continent_pal, name="Global region") + scale_x_continuous(expand = c(0,0),limits=c(1975,2025), breaks=seq(1980,2020,by=10)) + 
		scale_y_continuous(expand=c(0,0)) + theme_classic() + leg + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) #limits=c(0,50), breaks=seq(10,40,by=10)
	
	dfs_gheat[[i]] <- dfs_gheat[[i]][ , !names(dfs_gheat[[i]]) %in% c("Sampling.Year","Country","phylotype")]
	rownames(dfs_gheat[[i]]) <- dfs_gheat[[i]]$Accession
	dfs_gheat[[i]]$Accession <- NULL 
	
	print("ntips:")
	print(length(rf$tip.label))
	
	print("nrows dfs_gheat:")
	print(nrow(dfs_gheat[[i]]))
	
	pl <- ggtree(trees_global_vois_r[[i]], layout = "rectangular", color="grey20") %<+% combined_seqs_md_adj[[i]] + geom_tippoint(aes(color=phylotype), size=2, alpha=.8) + geom_treescale(x=0.1, y=20, fontsize=4, linesize=1, offset=1) 
	p1[[i]] <- gheatmap(pl, dfs_gheat[[i]], width = 0.2,colnames=FALSE) + #offset = 1#,colnames_offset_y = 1, font.size=1.35,custom_column_labels=c("Georegion")
		scale_fill_manual(values=continent_pal, name = "Global region", drop=FALSE) + scale_color_manual(values=phylot_pal, name="Label", drop=FALSE) + coord_cartesian(clip = "off") +
		theme(legend.text=element_text(size=6), legend.title=element_text(size=8), legend.key.height=unit(.5, "cm"),legend.key.width=unit(.5, "cm"),
								axis.text.x=element_blank(),axis.title.x = element_blank())
	
	pl2[[i]] <- ggtree(trees_global_vois_r[[i]]) %<+% combined_seqs_md_adj[[i]] + geom_tippoint(aes(color=Georegion), size=3, alpha=.75) + geom_treescale(x=0.12, y=40, fontsize=6, linesize=1, offset=1)
	
	pl3 <- ggtree(trees_global_vois_r[[i]]) %<+% combined_seqs_md_adj[[i]] + geom_tippoint(aes(color=Country), size=3, alpha=.75) + geom_treescale(x=0.12, y=40, fontsize=6, linesize=1, offset=1)
	p3[[i]] <- gheatmap(pl3, dfs_gheat[[i]], width = 0.2,colnames=FALSE) +
		scale_fill_manual(values=continent_pal, name = "Global region") + coord_cartesian(clip = "off") +
		theme(legend.text=element_text(size=6), legend.title=element_text(size=8), legend.key.height=unit(.25, "cm"),legend.key.width=unit(.25, "cm"),
								axis.text.x=element_blank(),axis.title.x = element_blank())
}

# Figure 4A and 4C (PT69 and 133 are indices 2 and 3)
f4a <- p1[[2]] + theme(legend.title=element_text(family = helv, size=10, color="black"),legend.text=element_text(family = helv, size=8, color="black")) +
	geom_text(aes(x=0.16, y=12), label = 'RefC.1986.\nETH2220', check_overlap = TRUE, color = 'grey10', size = 2.5, family=helv) + theme(legend.position = c(0.55, 0.75))

f4c <- p1[[3]] + theme(legend.title=element_text(family = helv, size=10, color="black"),legend.text=element_text(family = helv, size=8, color="black")) + 
	geom_text(aes(x=0.18, y=12), label = 'RefC.1986.\nETH2220', check_overlap = TRUE, color = 'grey10', size = 2.5, family=helv) + theme(legend.position = c(0.55, 0.75))

f4ac <- ggarrange(f4a, f4c, nrow=2, ncol=1, labels=c("A","C"), font.label=list(family=helv, color="black",size=10))
saveRDS( f4ac, glue("{RDS_PATH}/f4ac.rds") )

# Figure S18: distribution of sequences over time across global regions
blank_axis <- theme(axis.text.x=element_blank(),axis.title.x=element_blank(), axis.title.y=element_blank())
fs16 <- ggarrange(pst[[1]] + blank_axis, pst[[2]] + blank_axis, pst[[3]] + theme(axis.title.x=element_blank(),axis.title.y=element_blank()), 
																		nrow=3,  ncol=1, labels = c("A","B","C"), font.label=list(family=helv, color="black",size=10), common.legend = T, legend="top") #, heights = c(2,1)
annotate_figure(fs16, left = text_grob("Number of sequences", rot = 90, vjust = 1, size=10, family = helv, face="bold"), bottom = text_grob("Year", size=10, family = helv, face="bold"))
suppressMessages( ggsave(file=glue("{RESULTS_PATH}/figs/figS18.eps"), dpi=600, width=6, height=8, bg="white") )
suppressMessages( ggsave(file=glue("{RESULTS_PATH}/figs/figS18.jpg"), dpi=600, width=6, height=8, bg="white") )

# Figure S20
p1[[1]] + theme(legend.title=element_text(family = helv, size=10, color="black"),legend.text=element_text(family = helv, size=8, color="black")) +
	geom_text(aes(x=0.15, y=15), label = 'RefC.1986.\nETH2220', check_overlap = TRUE, color = 'grey10', size = 2.5, family=helv)
ggsave(file=glue("{RESULTS_PATH}/figs/figS20.jpg"), dpi=600, width=8, height=10, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS20.eps"), device = cairo_ps, dpi=600, width=8, height=10, bg="white")

### END FIGURES 4A, 4C, S16, S20 ####