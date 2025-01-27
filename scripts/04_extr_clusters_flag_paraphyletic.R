libs_load <- c("glue","ape", "rlist", "ggtree", "data.table")
invisible( lapply(libs_load, library, character.only=TRUE) )

treestruct_min_cl_size_res_yes_sup <- readRDS("rds/treestruct_min_cl_size_res_yes_sup.rds")

min_cl_size_choices <- c(30, 50, 100)
tree_names <- c("A_A1","CRF_02_AG","C","B")

treestruct_min_cl_size_res_yes_sup <- readRDS("rds/OLD_treestruct_min_cl_size_res_yes_sup.rds")
# IMPORTANT: After doing VL and CD4 analysis later in the workflow, the paraphyletic PT.B.5 (n=333, mcs=30) was found significantly different (VOI).
# We decided to join it with its descendants (matching exactly with PT.B.6, n=31) making it monophyletic
# Merging clusters/PTs 5 and 6 from subtype B (mcs=30, index=[[1,4]]) to subsequent analyses (so new cluster is PT.B.5.UK)
# treestruct_min_cl_size_res_yes_sup[[1,4]][[2]]$cluster[ treestruct_min_cl_size_res_yes_sup[[1,4]][[2]]$cluster == 6 ] <- 5
# PT6 will be empty now

# Flag paraphyletic phylotypes
treestruct_check_paraphyly <- treestruct_check_paraphyly_pt <- treestruct_check_paraphyly_pt2 <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(timetrees))
for(j in 1:length(min_cl_size_choices)) { #length(min_cl_size_choices)
	print(glue("min_cl_size_choices={min_cl_size_choices[j]}"))
	for(k in 1:length(timetrees)) {
		print(glue("tree_names={tree_names[k]}"))
		treestruct_check_paraphyly[[j,k]] <- base::split(treestruct_min_cl_size_res_yes_sup[[j,k]][[2]], treestruct_min_cl_size_res_yes_sup[[j,k]][[2]]$cluster)
		treestruct_check_paraphyly_pt[[j,k]] <- treestruct_check_paraphyly[[j,k]]
		print(names(treestruct_check_paraphyly_pt[[j,k]]))
		
		for(i in 1:length(treestruct_check_paraphyly_pt[[j,k]])) {
			print(glue("phylotype={i}"))
			treestruct_check_paraphyly_pt2[[j,k]][[i]] <- treestruct_check_paraphyly_pt[[j,k]][[i]]$taxon
			
			# if paraphyletic (non-monophyletic)
			if(!is.monophyletic(ml_trees[[k]], treestruct_check_paraphyly_pt2[[j,k]][[i]])) {
				print("Paraphyletic!")
				treestruct_check_paraphyly_pt[[j,k]][[i]]$paraphyletic <- ifelse(treestruct_check_paraphyly_pt[[j,k]][[i]]$cluster == i, yes="yes", no="no")
			} else {
				treestruct_check_paraphyly_pt[[j,k]][[i]]$paraphyletic <- ifelse(treestruct_check_paraphyly_pt[[j,k]][[i]]$cluster != i, yes="yes", no="no")
			}
		}
		treestruct_min_cl_size_res_yes_sup[[j,k]][[2]] <- rbindlist(treestruct_check_paraphyly_pt[[j,k]])
	}
}

#unique(treestruct_min_cl_size_res_yes_sup[[x,y]][[2]]$cluster [ treestruct_min_cl_size_res_yes_sup[[x,y]][[2]]$paraphyletic == "yes" ] )
saveRDS(treestruct_min_cl_size_res_yes_sup, "rds/treestruct_min_cl_size_res_yes_sup.rds")

# Split to get a data frame for each cluster
get_clusters_tree <- function(struct_data, struct_obj_tree, sampleTimes, out_folder) {
	#uniq_clusts <- unique(struct_data$cluster); uniq_clusts <- as.integer(uniq_clusts); uniq_clusts <- sort(uniq_clusts); print(uniq_clusts)
	all_clusts <- list()
	sample_times <- list()
	system(glue("mkdir -p results/03_treestructure/{out_folder}/"))
	for(u in 1: max(as.integer(struct_data$cluster)) ) { #uniq_clusts
		if(!identical(struct_data$taxon[ struct_data$cluster==u ], character(0))) { # if not empty (to make sure PT.6 does not raise an error)
			print(u)
			all_clusts[[u]] <- with(struct_data, ape::keep.tip(struct_obj_tree,  taxon[ cluster==u ])) #uniq_clusts[u]
			sample_times[[u]] <- sampleTimes[names(sampleTimes) %in% all_clusts[[u]]$tip.label]
			tr <- ggtree(all_clusts[[u]], color="grey30") + geom_tiplab(size=1)
			suppressMessages( ggsave(file=glue("results/03_treestructure/{out_folder}/phylotype_{u}.png"), plot=tr, dpi=600, limitsize=FALSE) ) #uniq_clusts[u]
		}
	}
	#list(all_clusts)
	list(all_clusts,sample_times)
}

sampleTimes_a1 <- readRDS("rds/subtype_a1_sampleTimes.rds")
sampleTimes_crf02ag <- readRDS("rds/subtype_crf02ag_sampleTimes.rds")
sampleTimes_c <- readRDS("rds/subtype_c_sampleTimes.rds")
sampleTimes_b <- readRDS("rds/subtype_b_sampleTimes.rds")
sampleTimes_all <- list(sampleTimes_a1, sampleTimes_crf02ag, sampleTimes_c, sampleTimes_b)

# Get all clusters: phylo object and no node.labels/bootstrap
extracted_clusters <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(timetrees))
for(j in 1:length(min_cl_size_choices)) { #length(min_cl_size_choices)
	print(glue("min_cl_size_choices={min_cl_size_choices[j]}"))
	for(k in 1:length(timetrees)) {
		print(glue("tree_names={tree_names[k]}"))
		# skip PT.6 for mcs=30 [[1,4]]
		extracted_clusters[[j,k]] <- get_clusters_tree(treestruct_min_cl_size_res_yes_sup[[j,k]][[2]], treestruct_min_cl_size_res_yes_sup[[j,k]][[1]]$tree, sampleTimes_all[[k]], glue('{tree_names[k]}/mcs{min_cl_size_choices[j]}'))
	}
}
saveRDS(extracted_clusters, "rds/extracted_clusters.rds")

seqs_folder_naive <- "data/subtype_seqs_naive"
fasta_c <- read.dna(glue("{seqs_folder_naive}/C_UK_final_aln_outliers_removed.fasta"), format="fasta")
fasta_c <- fasta_c[rownames(fasta_c) %in% timetrees[[3]]$tip.label,]
write.FASTA(fasta_c, glue("{seqs_folder_naive}/C_UK_final_aln_outliers_removed_final.fasta"))
fasta_b <- read.dna(glue("{seqs_folder_naive}/B_UK_final_aln_outliers_removed.fasta"), format="fasta")
fasta_b <- fasta_b[rownames(fasta_b) %in% timetrees[[4]]$tip.label,]
write.FASTA(fasta_b, glue("{seqs_folder_naive}/B_UK_final_aln_outliers_removed_final.fasta"))

fasta_files <- c(glue("{seqs_folder_naive}/A_A1_curated_refB_aln_len_filter.fasta"), 
																	glue("{seqs_folder_naive}/CRF_02_AG_curated_refB_aln_len_filter.fasta"),
																	glue("{seqs_folder_naive}/C_UK_final_aln_outliers_removed_final.fasta"),
																	glue("{seqs_folder_naive}/B_UK_final_aln_outliers_removed_final.fasta"))

# Remove paraphyletic clades from PT aln generation
treestruct_min_cl_size_res_yes_sup <- readRDS("rds/treestruct_min_cl_size_res_yes_sup.rds")
extracted_clusters <- readRDS("rds/extracted_clusters.rds")

# compute backbone phylotype
rm_paraphyletic_pt_alns <- backbone_cl_control <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(timetrees))
for(x in 1:length(min_cl_size_choices)) {
	aux_ntips_para <- c()
	print(glue("min_cl_size_choices={min_cl_size_choices[x]}"))
	for(y in 1:length(timetrees)) {
		print(glue("tree_names={tree_names[y]}"))
		rm_paraphyletic_pt_alns[[x,y]] <- as.character(unique(treestruct_min_cl_size_res_yes_sup[[x,y]][[2]]$cluster [ treestruct_min_cl_size_res_yes_sup[[x,y]][[2]]$paraphyletic == "yes" ] ))
		for(z in 1:length(rm_paraphyletic_pt_alns[[x,y]])) {
			print(glue("paraphyletic clade id={rm_paraphyletic_pt_alns[[x,y]][z]}"))
			#print(glue("ntips: {nrow(treestruct_min_cl_size_res_yes_sup[[x,y]][[2]][ treestruct_min_cl_size_res_yes_sup[[x,y]][[2]]$cluster == rm_paraphyletic_pt_alns[[x,y]][z], ])}"))
			#names(aux_ntips_para[[z]]) <- rm_paraphyletic_pt_alns[[x,y]][z]
			aux_ntips_para[z] <- nrow(treestruct_min_cl_size_res_yes_sup[[x,y]][[2]][ treestruct_min_cl_size_res_yes_sup[[x,y]][[2]]$cluster == rm_paraphyletic_pt_alns[[x,y]][z], ])
		}
		names(aux_ntips_para) <- rm_paraphyletic_pt_alns[[x,y]]
		print(aux_ntips_para)
		backbone_cl_control[[x,y]] <- names( aux_ntips_para [which(aux_ntips_para == max(aux_ntips_para))] ) #aux_ntips_para[which.max(aux_ntips_para))]
		print(glue("backbone phylotype: {backbone_cl_control[[x,y]]}"))
		aux_ntips_para <- c()
	}
}

saveRDS(backbone_cl_control, "rds/backbone_cl_control.rds")
saveRDS(rm_paraphyletic_pt_alns, "rds/rm_paraphyletic_pt_alns.rds")

# http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants<-function(tre,node,curr=NULL){
	if(is.null(curr)) curr<-vector()
	daughters<-tre$edge[which(tre$edge[,1]==node),2]
	curr<-c(curr,daughters)
	w<-which(daughters>=length(tre$tip))
	if(length(w)>0) for(i in 1:length(w))
		curr<-getDescendants(tre,daughters[w[i]],curr)
	return(curr)
}

# resolved paraphyletic phylotypes by tracing the most recent common ancestor (MRCA) of all sequences in the phylotype until a monophyletic group was formed
merge_paraphyletic_desc <- function(tre, extr_clust_para) {
	print(glue("ntips BEFORE: {length(extr_clust_para)}"))
	mrcaa <- ape::getMRCA(phy=tre, tip=extr_clust_para)
	descs <- getDescendants(tre, node=mrcaa)
	descs_indices <- descs[descs <= ape::Ntip(tre)] 
	desc_ids <- unlist( lapply( 1:(ape::Ntip(tre)+ape::Nnode(tre)) , function(x) { na.exclude(tre$tip.label[descs_indices[x]]) }) )
	tree_adj <- keep.tip(tre, desc_ids)
	print(glue("ntips AFTER: {length(tree_adj$tip.label)}"))
	tree_adj
}

merge_para <- merge_para_timetr <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		#for(k in as.integer( rm_paraphyletic_pt_alns[[i,j]]) ) {
		for(k in 1:length(extracted_clusters[[i,j]][[1]])) { 
			if(!(k %in% as.integer(backbone_cl_control[[i,j]])) ) {
				print(glue("{min_cl_size_choices[i]}-{tree_names[j]}-PT{k}"))
				#merge_para[[i,j]][[k]] <- merge_paraphyletic_desc(ml_trees[[j]], extracted_clusters[[i,j]][[1]][[k]]$tip.label)
				merge_para_timetr[[i,j]][[k]] <- merge_paraphyletic_desc(timetrees[[j]], extracted_clusters[[i,j]][[1]][[k]]$tip.label)
			}
		#}
		}
	}
}

#saveRDS(merge_para, "rds/merge_para.rds")
saveRDS(merge_para_timetr, "rds/merge_para_timetr.rds")