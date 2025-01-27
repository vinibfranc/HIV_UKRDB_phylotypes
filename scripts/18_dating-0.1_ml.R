libs_load <- c("ggtree","ggplot2","ape","treedater","glue","caper","data.table")
invisible( lapply(libs_load, library, character.only=TRUE) )

NCPU <- 1
RDS_PATH <- "rds"

options(scipen=999)
ml_dating <- function(voi_id, md) {
	print("=========")
	print(glue("PT.B.{voi_id}.UK"))
	print("=========")
	tr <- read.tree( glue('{DB_PATH}/queries/iqtree/aln_comb_phylotypes_global_refC_{voi_id}.fasta.treefile' ))
	md <- md
	print(glue("nrow md: {nrow(md)}"))
	print(glue("ntips: {length(tr$tip.label)}"))
	tr <- keep.tip(tr, intersect(tr$tip.label, md$Accession) )
	sts <- md$Sampling.Year[ match( tr$tip.label, md$Accession ) ] |> setNames( tr$tip.label )
	
	# drop ref C
	tr <- drop.tip( tr, tr$tip.label[grepl(tr$tip.label, patt = '^Ref.C.*')] )
	# for PTs 40 and 133: clock = 'additive', omega0 = 0.001152089, meanRateLimits=c(0.0008, 0.0015)
	# for PT69: clock = 'uncorrelated' (NO omega0 and meanRateLimits)
	if(voi_id == 69) {
		td <- dater ( (tr), sts = sts, s = 995, clock = 'uncorrelated', ncpu = NCPU )
	} else {
		td <- dater ( (tr), sts = sts, s = 995, clock = 'additive', omega0 = 0.001152089, meanRateLimits=c(0.0008, 0.0015), ncpu = NCPU ) #omega0 = 0.001152089,meanRateLimits=c(0.0008, 0.0015), parallel_foreach = TRUE
	}
	print("treedater coef_var_rates before outlier removal:")
	print(td$coef_of_variation)
	plot( td, cex=0.7 )
	system("mkdir -p results/18_dating_ml/")
	png(glue('results/18_dating_ml/bef_outl_rm_rtt_PT_{voi_id}.png'))
	rootToTipRegressionPlot(td)
	dev.off()
	ot <- outlierTips( td )
	print("Mean q for seqs <.05 before removing outliers:")
	print(mean( ot$q <.05 ))
	tr.1 <- drop.tip( tr, ot$taxon [ ot$q < .05 ] )
	print(tr.1)
	if(voi_id == 69) {
		td.1 <- dater ( (tr.1), sts = sts, s = 995, clock = 'uncorrelated', ncpu = NCPU )
	} else {
		td.1 <- dater ( (tr.1), sts = sts, s = 995, clock = 'additive', omega0 = 0.001152089, meanRateLimits=c(0.0008, 0.0015), ncpu = NCPU ) #omega0 = 0.001152089,meanRateLimits=c(0.0008, 0.0015), parallel_foreach = TRUE
	}
	print("treedater output after outlier removal:")
	print(td.1$coef_of_variation)
	ot.1 <- outlierTips( td.1 )
	print("outliers after <.05 removal?")
	print(head( ot.1))
	plot( td.1, cex=0.7 )
	png(glue('results/18_dating_ml/aft_outl_rm_rtt_PT_{voi_id}.png'))
	rootToTipRegressionPlot(td.1)
	dev.off()
	pb <- treedater::parboot( td.1, ncpu = NCPU , nrep = NREPS, quiet=FALSE) -> pb_out #500
	print("parboot output:")
	print(pb_out)
	
	md <- md[md$Accession %in% tr.1$tip.label,]
	
	a <- getMRCA( phy=tr.1, tip=md$Accession[md$phylotype == glue("PT.B.{voi_id}.UK")] ) 
	
	tmrca <- sapply( pb$trees, function(x) node.depth.edgelength( x )[a] + x$timeOfMRCA )
	print("TMRCAs:")
	print(tmrca)
	png(glue('results/18_dating_ml/td_hist_tmrcas_{voi_id}.png'))
	hist(tmrca)
	dev.off()
	
	return(list(tr=tr.1, td=td.1, pb=pb, tmrca=tmrca))
}

combined_seqs_md_adj <- readRDS(glue("{RDS_PATH}/combined_seqs_md_adj.rds"))

NREPS <- 1000

res_ml_dating <- list()
for(i in 1:length(subtype_b_vois)) {
	print(glue("=====PT[i]===="))
	res_ml_dating[[i]] <- ml_dating(subtype_b_vois_ids_only[i], combined_seqs_md_adj[[i]])
}

saveRDS(res_ml_dating, glue("{RDS_PATH}/res_ml_dating.rds"))

hist(res_ml_dating[[1]]$tmrca, xlim = c(1980,2020), breaks=1000)
hist(res_ml_dating[[2]]$tmrca, xlim = c(1980,2020), breaks=1000)
hist(res_ml_dating[[3]]$tmrca, xlim = c(1980,2020), breaks=1000)

# join tmrcas together
list_tmrcas <- list()
for(i in 1:length(subtype_b_vois)) {
	list_tmrcas[[i]] <- data.frame(repl=seq(1,length(res_ml_dating[[i]]$tmrca)), tmrca=res_ml_dating[[i]]$tmrca, phylotype=subtype_b_vois[i], method="Treedater parboot")
}
df_tmrcas <- rbindlist(list_tmrcas)
df_tmrcas$phylotype <- factor(df_tmrcas$phylotype, levels = subtype_b_vois)
pt_pal <- c("PT.B.40.UK"="#56B4E9", "PT.B.69.UK"="#D55E00", "PT.B.133.UK"="#009E73")

saveRDS(list_tmrcas, glue("{RDS_PATH}/list_tmrcas.rds"))
saveRDS(df_tmrcas, glue("{RDS_PATH}/df_tmrcas.rds"))