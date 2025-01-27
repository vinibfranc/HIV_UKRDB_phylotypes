# Derive consensus sequence for each phylotype

libs_load <- c("glue", "ape", "seqinr", "lubridate")
invisible( lapply(libs_load, library, character.only=TRUE) )

RESULTS_PATH="results"

# Read from results/04_aln_all_mono_subtypes since it already resolves paraphyly by tracing the MRCA of all sequences in the phylotype until a monophyletic group was formed
d_cons <- d_cons_join <- list()
tree_names_adj <- c("A1","CRF02AG","C","B")
for(i in 1:length(tree_names_adj)) {
	for(j in 1:length(extracted_clusters[[1,i]][[1]])) {
		if(!(j %in% as.integer(backbone_cl_control[[1,i]])) ) {
			if(!is.null(extracted_clusters[[1,i]][[1]][[j]]$tip.label)) {
				print(glue("{i}-{j}"))
				r <- read.FASTA(glue("results/04_aln_all_mono_subtypes/mcs30/{tree_names_adj[i]}/PT_{tree_names_adj[i]}_{j}_UK.fasta", type="DNA"))
				d_cons[[j]] <- data.frame(testindex=names(r), phylotype=rep(glue("{tree_names_adj[i]}@{j}"), length(names(r))))
			}
		}
	}
	d_cons_join[[i]] <- rbindlist(d_cons)
}

d_cons_all <- rbindlist(d_cons_join)

# get seq id's for each phylotype
PTS <- sort( unique( d_cons_all$phylotype ))
sids <- lapply( PTS, function(pt){
	d_cons_all$testindex[ d_cons_all$testindex %in% d_cons_all$testindex[d_cons_all$phylotype==pt] ]
})

# read sequences 
s <- read.alignment(glue('data/seqs_20230518.fasta'), format = 'fasta')

cs <- lapply( sids, function(ptsids){
        s1 <- seqinr::as.alignment(nb=length(ptsids), nam=ptsids, seq=s$seq[match(ptsids, paste0("t.",s$nam)) ] )
        seqinr::consensus(s1,  method = 'majority', threshold = 0.8)
})

system(glue("mkdir -p {RESULTS_PATH}/15_consensus/"))
write.fasta( cs, names=PTS, file = glue('{RESULTS_PATH}/15_consensus/pt_consensus.fasta'))

# GenBank submission below
d_cons <- d_cons_join <- list()
tree_names_adj <- c("A1","CRF02AG","C","B")
for(i in 1:length(tree_names_adj)) {
	for(j in 1:length(extracted_clusters[[1,i]][[1]])) {
		if(!(j %in% as.integer(backbone_cl_control[[1,i]])) ) {
			if(!is.null(extracted_clusters[[1,i]][[1]][[j]]$tip.label)) {
				print(glue("{i}-{j}"))
				r <- read.FASTA(glue("results/04_aln_all_mono_subtypes/mcs30/{tree_names_adj[i]}/PT_{tree_names_adj[i]}_{j}_UK.fasta", type="DNA"))
				d_cons[[j]] <- data.frame(testindex=names(r), phylotype=rep(glue("PT.{tree_names_adj[i]}.{j}.UK"), length(names(r))))
			}
		}
	}
	d_cons_join[[i]] <- rbindlist(d_cons)
}

d_cons_all <- rbindlist(d_cons_join)

# get seq id's for each phylotype
PTS <- sort( unique( d_cons_all$phylotype ))
sids <- lapply( PTS, function(pt){
	d_cons_all$testindex[ d_cons_all$testindex %in% d_cons_all$testindex[d_cons_all$phylotype==pt] ]
})

# read sequences 
s <- read.alignment(glue('data/seqs_20230518.fasta'), format = 'fasta')

cs <- lapply( sids, function(ptsids){
	s1 <- seqinr::as.alignment(nb=length(ptsids), nam=ptsids, seq=s$seq[match(ptsids, paste0("t.",s$nam)) ] )
	seqinr::consensus(s1,  method = 'majority', threshold = 0.8)
})

system(glue("mkdir -p {RESULTS_PATH}/15_consensus/"))
write.fasta( cs, names=glue("{PTS}_phylotype_consensus [organism=Human immunodeficiency virus 1] pol gene, partial cds"), file = glue('{RESULTS_PATH}/15_consensus/pt_consensus_genbank.fasta'))

# trim end manually to keep 1302 nucleotides of length and delete all gaps in all sequences

# get rane of dates of seqs in phylotypes
seq_ranges_df <- readRDS(file=glue("{RDS_PATH}/demog_seq_filters_df.rds"))
seq_ranges_df$testindex <- paste0("t.",seq_ranges_df$testindex)
#seq_ranges_df <- seq_ranges_df[seq_ranges_df$testindex %in% d_cons_all$testindex[d_cons_all$phylotype==pt],]
dates <- lapply( PTS, function(pt){
	floor(range(seq_ranges_df$dbsample_date[ seq_ranges_df$testindex %in% d_cons_all$testindex[d_cons_all$phylotype==pt] ],na.rm=T))
})
dates <- sapply(dates, function(x) paste(x, collapse = "-"))

notes_consensus <- "These sequences were derived from a cluster/phylotype analysis performed on the UK Drug Resistace Database dataset (no original GenBank submission of complete dataset) for subtypes A1, B, C, and CRF02AG. Therefore, collection dates represent the range of dates of the individual sequences for which a phylotype was called. Phylotypes were identified using subtype-specific time-scaled trees estimated using treedater v0.5.3. Samples sizes of those trees were, respectively, 2714, 24100, 11331, and 2743. Phylotypes were estimated using treestructure v0.3.1 with a minimum clade size of 30 and a bootstrap support threshold of 80%. The consensus phylotype sequences were obtained using seqinr v4.2.8 with the site majority method (i.e. high-frequency character is returned as consensus) and a minimum relative frequency threshold of 80%. Initially, some phylotypes were paraphyletic (n=30 of the submitted consensuses, listed IDs in paper Table S5; preprint available at https://doi.org/10.2139/ssrn.4929798). We resolved this by tracing the most recent common ancestor (MRCA) of all sequences in the phylotype until a monophyletic group was formed. Then called the consensus as mentioned above. There is no consensus sequence uploaded for the backbone phylotype (paraphyletic) of each subtype, which represents the large (~50% of samples) and diverse ancestral type from which all of the other phylotypes were descended."
df_seqs <- data.frame(Sequence_ID=glue("{PTS}_phylotype_consensus"), Collection_date=dates, Country="United Kingdom", Host="Homo sapiens", Isolate=PTS, Note=notes_consensus)
write.table(df_seqs, file=glue("{RESULTS_PATH}/15_consensus/source_modifier_table_UK_phylotypes.txt"), sep = "\t", quote = F, row.names = F)
