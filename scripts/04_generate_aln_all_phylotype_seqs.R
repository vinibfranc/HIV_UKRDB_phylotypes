# Generate alignments for each phylotype. If backbone does not generate, if paraphyletic then use merge_para, 
# which absorbs the monophyla nested inside them, to include all seqs that make clade monophyletic

libs_load <- c("ape","glue")
invisible( lapply(libs_load, library, character.only=TRUE) )

tree_names_aj <- c("A1", "CRF02AG", "C", "B")

folder_aln_pts <- "results/04_aln_all_mono_subtypes"
system(glue("mkdir -p {folder_aln_pts}"))
fastas_subtypes <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(timetrees))

for(i in 1:length(min_cl_size_choices)) {
	system(glue("mkdir -p {folder_aln_pts}/mcs{min_cl_size_choices[i]}"))
	print(min_cl_size_choices[i])
	for(j in 1:length(fasta_files)) {
		fastas_subtypes[[i,j]] <- read.dna(fasta_files[j], format="fasta")
		if(!(j %in% c(3,4))) { # subtypes C and B already have t.
			rownames(fastas_subtypes[[i,j]]) <- paste0("t.",rownames(fastas_subtypes[[i,j]]))
		}
		system(glue("mkdir -p {folder_aln_pts}/mcs{min_cl_size_choices[i]}/{tree_names_aj[j]}"))
		print(tree_names_aj[j])
		for(k in 1:length(extracted_clusters[[i,j]][[1]])) {
			if(!(k %in% as.integer(backbone_cl_control[[i,j]])) ) {
				if( k %in%  as.integer( rm_paraphyletic_pt_alns[[i,j]]) ) {
					fs <- fastas_subtypes[[i,j]][( rownames(fastas_subtypes[[i,j]]) %in% merge_para[[i,j]][[k]]$tip.label), ]
				} else {
					fs <- fastas_subtypes[[i,j]][( rownames(fastas_subtypes[[i,j]]) %in% extracted_clusters[[i,j]][[1]][[k]]$tip.label), ]
				}
				write.FASTA(fs, file=glue("{folder_aln_pts}/mcs{min_cl_size_choices[i]}/{tree_names_aj[j]}/PT_{tree_names_aj[j]}_{k}_UK.fasta"))
				#print(fs)
			}
		}
	}
}