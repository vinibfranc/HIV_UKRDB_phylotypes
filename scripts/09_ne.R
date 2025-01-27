libs_load <- c("dplyr", "glue","ape","treedater","mlesky","lubridate","treestructure","ggtree","ggpubr","ggplot2","viridis","data.table")
invisible( lapply(libs_load, library, character.only=TRUE) )

NCPU <- 4
extracted_clusters <- readRDS("rds/extracted_clusters.rds")
size_label_adjust <- theme(axis.text=element_text(size=7), axis.title=element_text(size=8), plot.title=element_text(size=10))
### ESTIMATE EFFECTIVE POPULATION SIZE (MLESKY) FOR EACH CLUSTER USING SKYKAPPA MODEL ###

# Estimate Ne for all phylotypes across minCladeSizes and subtypes with res=50 and tau chosen by cross-val
fit_mlesky_indiv_clusters <- function(extr_clusters, excl_paraphyletic, out_folder) {
	system(glue("mkdir -p results/09_mlesky_skygrid/{out_folder}/"))
	fits_mlesky_m2 <- pboot_mlesky_m2 <- list()
	print(glue("length clusters: {length(extr_clusters[[1]])}"))
	for(l in 1:length(extr_clusters[[1]])) { 
		if(!(l %in% excl_paraphyletic)) {
			print(glue("[l={l}]"))
			print(length(extr_clusters[[1]][[l]]$tip.label))
			# Skygrid
			class(extr_clusters[[1]][[l]]) <- "phylo"
			extr_clusters[[1]][[l]]$node.label <- NULL
			fits_mlesky_m2[[l]] <- mlskygrid(extr_clusters[[1]][[l]], sampleTimes=extr_clusters[[2]][[l]], res= 50, tau=NULL, tau_lower=.01, tau_upper=100, model=2, ncpu=NCPU)
			#print(fits_mlesky_m2[[l]] )
			if(length( unique( fits_mlesky_m2[[l]]$time)) <= 2) {
				next
			}
			pboot_mlesky_m2[[l]] <- mlesky::parboot(fits_mlesky_m2[[l]], nrep=100, ncpu=NCPU, dd=FALSE)
			#print(pboot_mlesky_m2[[l]])
			p1_ne <- plot(pboot_mlesky_m2[[l]], ggplot=TRUE, logy=FALSE) + ggplot2::ylab("Ne") + ggplot2::xlab("Time") + ggplot2::ggtitle(paste0("Phylotype ",l)) + theme_classic() + size_label_adjust #ggplot2::coord_cartesian(ylim=c(0, 75000))
			pboot_mlesky_m2[[l]]$growth <- pboot_mlesky_m2[[l]]$growthrate
			p1_gr <- plot(pboot_mlesky_m2[[l]], growth=TRUE, ggplot=TRUE, logy=FALSE)+ ggplot2::ylab("Growth rate") + ggplot2::xlab("Time") + theme_classic() + size_label_adjust
			p1 <- ggarrange(p1_ne, p1_gr, ncol=2, nrow=1)
			
			suppressMessages( ggsave(file=glue("results/09_mlesky_skygrid/{out_folder}/PT_{l}.png"),  plot=p1, dpi=600, width=7, height=5) )
		} #else {pboot_mlesky_m2[[l]] <- NULL}
	}
	list(pboot_mlesky_m2)
}

# All minCladeSize (i) and subtypes (j)
fits_mlesky <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
#fits_mlesky <- readRDS("rds/fits_mlesky.rds")
for(i in 1:length(min_cl_size_choices)) { 
	print(glue("min_cl_size_choices={min_cl_size_choices[i]}"))
	for(j in 1:length(tree_names)) {
		print(glue("tree_names={tree_names[j]}"))
		fits_mlesky[[i,j]] <- fit_mlesky_indiv_clusters(extr_clusters=extracted_clusters[[i,j]], excl_paraphyletic=c(rm_paraphyletic_pt_alns[[i,j]]), out_folder=glue("mcs{min_cl_size_choices[i]}/{tree_names[j]}"))
	}
}

saveRDS(fits_mlesky, "rds/fits_mlesky.rds")

for(i in subtype_b_vois_ids_only) {
	system(glue("mkdir -p data/09_ne_vois/"))
	write.tree(extracted_clusters[[1,4]][[1]][[i]], glue("data/09_ne_vois/phylotype_{i}_timetree.nwk"))
}


##### NOTE: do only for VOIs to see impact of choice of tau on Ne
fit_mlesky_selected_clusters_fixed_tau <- function(list_phylotypes, extr_clusters, tau_crossval=FALSE, tau, dd_fast, fprefix) {
	fits_mlesky_m2 <- pboot_mlesky_m2 <- list()
	cnt <- 1
	for(l in list_phylotypes) {
		if(!is.null(extr_clusters[[1]][[l]])) {
			print(glue("[pt={l}; tau={tau}]"))
			class(extr_clusters[[1]][[l]]) <- "phylo"
			if(tau_crossval)
				fits_mlesky_m2[[cnt]] <- mlskygrid(extr_clusters[[1]][[l]], sampleTimes=extr_clusters[[2]][[l]], res= 50, tau=NULL, tau_lower=1, tau_upper=100, model=2, ncpu=NCPU)
			else 
				fits_mlesky_m2[[cnt]] <- mlskygrid(extr_clusters[[1]][[l]], sampleTimes=extr_clusters[[2]][[l]], res= 50, tau=tau, model=2, ncpu=NCPU)
			if(length( unique( fits_mlesky_m2[[cnt]]$time)) <= 2) {
				next
			}
			pboot_mlesky_m2[[cnt]] <- mlesky::parboot(fits_mlesky_m2[[cnt]], nrep=100, ncpu=NCPU, dd=dd_fast)
			#print(pboot_mlesky_m2[[l]])
			p1_ne <- plot(pboot_mlesky_m2[[cnt]], ggplot=TRUE, logy=FALSE) + ggplot2::ylab("Ne") + ggplot2::xlab("Time") + ggplot2::ggtitle(paste0("Phylotype ",l)) + theme_classic() + size_label_adjust #ggplot2::coord_cartesian(ylim=c(0, 75000))
			pboot_mlesky_m2[[cnt]]$growth <- pboot_mlesky_m2[[cnt]]$growthrate
			p1_gr <- plot(pboot_mlesky_m2[[cnt]], growth=TRUE, ggplot=TRUE, logy=FALSE)+ ggplot2::ylab("Growth rate") + ggplot2::xlab("Time") + theme_classic() + size_label_adjust
			p1 <- ggarrange(p1_ne, p1_gr, ncol=2, nrow=1)
			
			system(glue("mkdir -p results/09_mlesky_tau_30B/{fprefix}/"))
			suppressMessages( ggsave(file=glue("results/09_mlesky_tau_30B/{fprefix}/phylotype_{l}.png"),  plot=p1, dpi=600, width=7, height=5) )
			
			cnt <- cnt+1
			print(cnt)
		}
	}
	list(pboot_mlesky_m2)
}
# add tau=50 to main figure (3A) and other tau choices to SI

# this is done for mcs=30 and subtype=B only
fits_mlesky_b_incr_res_fix_tau100 <- fit_mlesky_selected_clusters_fixed_tau(subtype_b_vois_ids_only, extracted_clusters[[1,4]], tau_crossval = FALSE, tau=100, dd_fast=FALSE, glue('tau100'))
fits_mlesky_b_incr_res_fix_tau75 <- fit_mlesky_selected_clusters_fixed_tau(subtype_b_vois_ids_only, extracted_clusters[[1,4]], tau_crossval = FALSE, tau=75, dd_fast=FALSE, glue('tau75'))
fits_mlesky_b_incr_res_fix_tau50 <- fit_mlesky_selected_clusters_fixed_tau(subtype_b_vois_ids_only, extracted_clusters[[1,4]], tau_crossval = FALSE, tau=50, dd_fast=FALSE, glue('tau50'))
fits_mlesky_b_incr_res_fix_tau25 <- fit_mlesky_selected_clusters_fixed_tau(subtype_b_vois_ids_only, extracted_clusters[[1,4]], tau_crossval = FALSE, tau=25, dd_fast=FALSE, glue('tau25'))
fits_mlesky_b_incr_res_fix_tau10 <- fit_mlesky_selected_clusters_fixed_tau(subtype_b_vois_ids_only, extracted_clusters[[1,4]], tau_crossval = FALSE, tau=10, dd_fast=FALSE, glue('tau10'))
fits_mlesky_b_incr_res_fix_tau5 <- fit_mlesky_selected_clusters_fixed_tau(subtype_b_vois_ids_only, extracted_clusters[[1,4]], tau_crossval = FALSE, tau=5, dd_fast=FALSE, glue('tau5'))
fits_mlesky_b_incr_res_fix_tau2.5 <- fit_mlesky_selected_clusters_fixed_tau(subtype_b_vois_ids_only, extracted_clusters[[1,4]], tau=2.5, dd_fast=FALSE, glue('tau2.5'))
# tau chosen by cross-validation
fits_mlesky_b_incr_res_fix_cv <- fit_mlesky_selected_clusters_fixed_tau(subtype_b_vois_ids_only, extracted_clusters[[1,4]], tau_crossval = TRUE, tau=NULL, dd_fast=FALSE, glue('tau_cv'))
# tau by cross-val, respectively: 1.28, 1.60, 5.7

saveRDS(fits_mlesky_b_incr_res_fix_tau50, glue("{RDS_PATH}/fits_mlesky_b_incr_res_fix_tau50.rds"))
#fits_all_fix_tau <- list(fits_mlesky_b_incr_res_fix_tau5, fits_mlesky_b_incr_res_fix_tau10, fits_mlesky_b_incr_res_fix_tau25, fits_mlesky_b_incr_res_fix_tau50, fits_mlesky_b_incr_res_fix_tau75, fits_mlesky_b_incr_res_fix_tau100)

### BEGIN FIGURE 3A AND S14 ###
fits_mlesky_b_incr_res_fix_tau50 <- unlist(fits_mlesky_b_incr_res_fix_tau50, recursive = F) # drop one level
f1b <- fits_1b <- list()
for(c in 1:length(subtype_b_vois_ids_only)) {
	fits_1b[[c]] <- fits_mlesky_b_incr_res_fix_tau50[[c]] #fits_mlesky_b_incr_res[[1,1]]
	f1b[[c]] <- fits_1b[[c]]$ne_ci
	f1b[[c]] <- as.data.frame(f1b[[c]])
	f1b[[c]]$time <- fits_mlesky_b_incr_res_fix_tau50[[c]]$time
	f1b[[c]]$phylotype <- subtype_b_vois_ids_only[c]
	#fit_mlesky_figs(extracted_clusters[[1,4]][[1]][[k]], extracted_clusters[[1,4]][[2]][[k]])
}

# Get Ne(t) for non-VOIs (for paraphyletic phylotypes getting merged version)
non_vois <- c(1:5,7:39,41:68,70:132,134:152,154) # already removed PT6 (empty) and backbone (153), removing 20 as well because of error in mlesky
#merge_para <- readRDS("rds/merge_para.rds")
merge_para_timetr <- readRDS(glue("{RDS_PATH}/merge_para_timetr.rds"))
rm_paraphyletic_pt_alns <- readRDS(glue("{RDS_PATH}/rm_paraphyletic_pt_alns.rds"))

subtype_b_sampleTimes <- readRDS("rds/subtype_b_sampleTimes.rds")
fit_mlesky_selected_clusters_fixed_tau_non_vois_merge_para <- function(list_phylotypes, extr_clusters, merged_para_, tau, dd_fast, fprefix) {
	fits_mlesky_m2 <- pboot_mlesky_m2 <- list()
	cnt <- 1
	for(l in list_phylotypes) { #list_phylotypes
		print(glue("[pt={l}; tau={tau}]"))
		if(l %in% rm_paraphyletic_pt_alns[[1,4]]) {
			print("PARA")
			class(merged_para_[[l]]) <- "phylo"
			fits_mlesky_m2[[cnt]] <- mlskygrid(merged_para_[[l]], sampleTimes=subtype_b_sampleTimes, res= 50, tau=tau, model=2, ncpu=NCPU)
			if(length( unique( fits_mlesky_m2[[cnt]]$time)) <= 2) {
				next
			}
		} else {
			print("MONO")
			class(extr_clusters[[1]][[l]]) <- "phylo"
			fits_mlesky_m2[[cnt]] <- mlskygrid(extr_clusters[[1]][[l]], sampleTimes=extr_clusters[[2]][[l]], res= 50, tau=tau, model=2, ncpu=NCPU)
			if(length( unique( fits_mlesky_m2[[cnt]]$time)) <= 2) {
				next
			}
		}
	
		pboot_mlesky_m2[[cnt]] <- mlesky::parboot(fits_mlesky_m2[[cnt]], nrep=100, ncpu=NCPU, dd=dd_fast)
		#print(pboot_mlesky_m2[[l]])
		p1_ne <- plot(pboot_mlesky_m2[[cnt]], ggplot=TRUE, logy=FALSE) + ggplot2::ylab("Ne") + ggplot2::xlab("Time") + ggplot2::ggtitle(paste0("Phylotype ",l)) + theme_classic() + size_label_adjust #ggplot2::coord_cartesian(ylim=c(0, 75000))
		pboot_mlesky_m2[[cnt]]$growth <- pboot_mlesky_m2[[cnt]]$growthrate
		p1_gr <- plot(pboot_mlesky_m2[[cnt]], growth=TRUE, ggplot=TRUE, logy=FALSE)+ ggplot2::ylab("Growth rate") + ggplot2::xlab("Time") + theme_classic() + size_label_adjust
		p1 <- ggarrange(p1_ne, p1_gr, ncol=2, nrow=1)
		
		system(glue("mkdir -p results/09_mlesky_tau_30B/{fprefix}/"))
		suppressMessages( ggsave(file=glue("results/09_mlesky_tau_30B/{fprefix}/phylotype_{l}.png"),  plot=p1, dpi=600, width=7, height=5) )
		
		cnt <- cnt+1
		print(cnt)
	}
	list(pboot_mlesky_m2)
}

fits_tau50_non_vois <- fit_mlesky_selected_clusters_fixed_tau_non_vois_merge_para(non_vois, extracted_clusters[[1,4]], merge_para_timetr[[1,4]], tau=50, dd_fast=FALSE, glue('tau50_non_vois'))
fits_tau50_non_vois <- unlist(fits_tau50_non_vois, recursive = FALSE)

saveRDS(fits_tau50_non_vois, glue("{RDS_PATH}/fits_tau50_non_vois.rds"))
#fits_tau50_non_vois <- readRDS(glue("{RDS_PATH}/fits_tau50_non_vois.rds"))

f1b_mean_ne <- fits_1b_mean <- list()
for(k in 1:length(non_vois)) {
	fits_1b_mean[[k]] <- fits_tau50_non_vois[[k]]
	f1b_mean_ne[[k]] <- fits_1b_mean[[k]]$ne_ci
	f1b_mean_ne[[k]] <- as.data.frame(f1b_mean_ne[[k]])
	f1b_mean_ne[[k]]$time <- fits_tau50_non_vois[[k]]$time
	f1b_mean_ne[[k]]$phylotype <- non_vois[k]
	#fit_mlesky_figs(extracted_clusters[[1,4]][[1]][[k]], extracted_clusters[[1,4]][[2]][[k]])
}

# Join all together to define common time axis
f1b_mean_ne_all <- rbindlist(f1b_mean_ne) #min(f1b_mean_ne_all$time)
common_time_ax = approx(f1b_mean_ne_all$time, f1b_mean_ne_all$ne, xout=seq(1980, max(f1b_mean_ne_all$time), length.out=nrow(f1b_mean_ne_all)/length(non_vois)), rule=2)$x

# Interpolate to get Ne estimates for the common time axis
f <- 1
ne_est_adj <- nelb_adj <- neub_adj <- ne_bnds <- list()
for(j in non_vois) {
	ne_est_adj[[f]] = approx(f1b_mean_ne[[f]]$time, f1b_mean_ne[[f]]$ne, xout=common_time_ax, rule=2)$y; ne_est_adj[[f]] = as.data.frame(ne_est_adj[[f]])
	nelb_adj[[f]] = approx(f1b_mean_ne[[f]]$time, f1b_mean_ne[[f]]$nelb, xout=common_time_ax, rule=2)$y; nelb_adj[[f]] = as.data.frame(nelb_adj[[f]])
	neub_adj[[f]] = approx(f1b_mean_ne[[f]]$time, f1b_mean_ne[[f]]$neub, xout=common_time_ax, rule=2)$y; neub_adj[[f]] = as.data.frame(neub_adj[[f]])
	ne_bnds[[f]] <- cbind(nelb_adj[[f]], ne_est_adj[[f]], neub_adj[[f]])
	ne_bnds[[f]]$time <- common_time_ax
	colnames(ne_bnds[[f]]) <- c("nelb","ne","neub","time")
	f <- f+1
}

# Average Ne estimates on the common time axis
ne_bnds_all <- rbindlist(ne_bnds)
ne_bnds_all_mean <- ne_bnds_all %>% group_by(time) %>% summarise(nelb=mean(nelb), ne=mean(ne), neub=mean(neub))
ne_bnds_all_median <- ne_bnds_all %>% group_by(time) %>% summarise(nelb=median(nelb), ne=median(ne), neub=median(neub)) # less prone to outliers
ne_bnds_all_median$phylotype <- "Median non-VOIs"
ne_bnds_all_median <- ne_bnds_all_median %>% dplyr::select( nelb,ne,neub,time,phylotype)

f1b_all <- rbindlist(f1b)
f1b_all$phylotype <- paste0("PT.B.",f1b_all$phylotype,".UK")
#f1b_all$phylotype <- factor(f1b_all$phylotype, levels = c("PT.B.9.UK","PT.B.21.UK","PT.B.23.UK","PT.B.39.UK","PT.B.49.UK","PT.B.50.UK","PT.B.68.UK"))

f1b_all_meds <- rbind(f1b_all, ne_bnds_all_median)
f1b_all_meds$phylotype <- factor(f1b_all_meds$phylotype, levels = c("Median non-VOIs", "PT.B.40.UK","PT.B.69.UK","PT.B.133.UK"),
																																	labels=c("Median non-\nVOIs", "PT.B.40.UK","PT.B.69.UK","PT.B.133.UK"))


#ne_pt_pal <- c("PT.B.5.UK"="#0099B4","PT.B.40.UK"="#088c06","PT.B.69.UK"="#ED0000","PT.B.133.UK"="#00468B","PT.B.137.UK"="#580618", "Median non-VOIs"="#ADB6B6")
ne_pt_pal <- c("PT.B.40.UK"="#56B4E9", "PT.B.69.UK"="#D55E00", "PT.B.133.UK"="#009E73", "Median non-\nVOIs"="#999999")

#f1b_all_meds
upper_y_3a <- 100
f1b_pl_t50_ok <- ggplot(data=f1b_all_meds, aes(x=time, y=ne, color=phylotype)) + geom_line(linewidth=1) +
	geom_ribbon(aes(ymin = nelb, ymax = neub, fill=phylotype), color=NA, alpha=0.10) +
	scale_color_manual(values=ne_pt_pal, name="Phylotype") + scale_fill_manual(values=ne_pt_pal, name="Phylotype") +
	#scale_alpha_manual(values=c(0.3, 0.3, 0.2, 0.25, 0.3, 0.2, 0.3, 0.2)) +
	coord_cartesian(ylim=c(0,upper_y_3a)) + labs(x="Year",y="Effective population size (Ne)") + theme_classic() + leg +
	scale_y_continuous(breaks=seq(0,upper_y_3a,by=25)) + #limits=c(0,900), 
	scale_x_continuous(limits=c(1995,2020)) +
	theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm")) #,legend.box.spacing = unit(1, "cm")
### END FIGURE 3A ###
f1b_pl_t50_ok

# Figure 3A
saveRDS(f1b_all_meds, glue("{RDS_PATH}/f1b_all_meds.rds"))
saveRDS(f1b_pl_t50_ok, glue("{RDS_PATH}/f1b_pl_t50_ok.rds"))

subtype_b_vois <- c("PT.B.40.UK","PT.B.69.UK","PT.B.133.UK")
diff_taus_plot <- function(fits_tau) {
	
	fits_tau <- unlist(fits_tau, recursive = F) # drop one level
	f1b <- fits_1b <- list()
	for(c in 1:length(subtype_b_vois_ids_only)) {
		fits_1b[[c]] <- fits_tau[[c]]
		f1b[[c]] <- fits_1b[[c]]$ne_ci
		f1b[[c]] <- as.data.frame(f1b[[c]])
		f1b[[c]]$time <- fits_tau[[c]]$time
		f1b[[c]]$phylotype <- subtype_b_vois_ids_only[c]
		#fit_mlesky_figs(extracted_clusters[[1,4]][[1]][[k]], extracted_clusters[[1,4]][[2]][[k]])
	}
	
	f1b_all <- rbindlist(f1b)
	f1b_all$phylotype <- paste0("PT.B.",f1b_all$phylotype,".UK")
	
	# f1b_all_meds <- rbind(f1b_all, ne_bnds_all_median)
	# f1b_all_meds$phylotype <- factor(f1b_all_meds$phylotype, levels = c("Median non-VOIs", "PT.B.5.UK","PT.B.40.UK","PT.B.69.UK","PT.B.133.UK","PT.B.137.UK"))
	f1b_all$phylotype <- factor(f1b_all$phylotype, levels = subtype_b_vois)
	
	upper_lim <- 400
	
	f1b_all_plot <- ggplot(data=f1b_all, aes(x=time, y=ne, color=phylotype)) + geom_line(linewidth=1) +
		geom_ribbon(aes(ymin = nelb, ymax = neub, fill=phylotype), color=NA, alpha=0.10) +
		scale_color_manual(values=ne_pt_pal, name="Phylotype") + scale_fill_manual(values=ne_pt_pal, name="Phylotype") +
		#scale_alpha_manual(values=c(0.3, 0.3, 0.2, 0.25, 0.3, 0.2, 0.3, 0.2)) +
		coord_cartesian(ylim=c(0,upper_lim)) + labs(x="Year",y="Effective population size (Ne)") + theme_classic() + leg +
		scale_y_continuous(breaks=seq(0,upper_lim,by=50)) +
		theme(plot.margin = margin(0.8, 0.8, 0.5, 0.5, "cm"))
	
	f1b_all_plot
}

f1b_pl_t100 <- diff_taus_plot(fits_mlesky_b_incr_res_fix_tau100)
f1b_pl_t75 <- diff_taus_plot(fits_mlesky_b_incr_res_fix_tau75)
#f1b_pl_t50 <- diff_taus_plot(fits_mlesky_b_incr_res_fix_tau50)
f1b_pl_t25 <- diff_taus_plot(fits_mlesky_b_incr_res_fix_tau25)
f1b_pl_t10 <- diff_taus_plot(fits_mlesky_b_incr_res_fix_tau10)
f1b_pl_t5 <- diff_taus_plot(fits_mlesky_b_incr_res_fix_tau5)
f1b_pl_t2.5 <- diff_taus_plot(fits_mlesky_b_incr_res_fix_tau2.5)

# Figure S15
mg <- theme(plot.margin = margin(0.3,0.5,0.5,0.5, "cm"))
s11 <- ggarrange(f1b_pl_t100 + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) + mg, 
																f1b_pl_t75 + theme(axis.text.x=element_blank(),axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank()) + mg, 
																f1b_pl_t25 + theme(axis.text.x=element_blank(),axis.text.y=element_blank(), axis.title.x = element_blank(), axis.title.y=element_blank()) + mg,
																f1b_pl_t10 + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + mg, 
																f1b_pl_t5 + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank()) + mg, 
																f1b_pl_t2.5 + theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank()) + mg, 
																nrow=2, ncol=3, common.legend = T, labels=c("A","B","C","D","E","F"), font.label=list(family="Helvetica", color="black",size=10)) + leg
annotate_figure(s11, left = text_grob("Effective population size (Ne)", rot = 90, vjust = 1, size=10, family = helv, face="bold"), bottom = text_grob("Year", size=10, family = helv, face="bold"))
ggsave(file=glue("{RESULTS_PATH}/figs/figS15.eps"), device=cairo_ps, dpi=600, width=8, height=6, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/figS15.jpg"), dpi=600, width=8, height=6, bg="white")

