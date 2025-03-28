libs_load <- c("glue","dplyr","forcats", "ape", "lubridate", "ggplot2",  "phytools", "data.table", "POUMM") #"motmot", "geiger", "phylolm"
invisible( lapply(libs_load, library, character.only=TRUE) )

### 1. Internode lengths, used as estimates of maximum transmission intervals, based on: ###
# https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1000590 and
# https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050050

RES_PATH <- "results"
RDS_PATH <- "rds"
# use merge_para_timetr from 04 that resolve paraphyletic PTs by tracing the MRCA of all sequences in the phylotype until a monophyletic group was formed

merge_para_timetr <- readRDS(glue("{RDS_PATH}/merge_para_timetr.rds"))
options(scipen=999)

# for mcs=30 and subtype B only
subtype_b_vois_ids_only <- c(40,69,133)
timetr_pt <- timetr_pt_only_internal <- list() #node_heights <- internode_intervals <- internode_intervals_df <- list()
internode_lengths <- internode_intervals <- internode_intervals_df <- list()
for(i in 1:length(merge_para_timetr[[1,4]])) {
	print(i)
	if(!is.null(merge_para_timetr[[1,4]][[i]])) {
		timetr_pt[[i]] <- merge_para_timetr[[1,4]][[i]]
		# get only node to node edges (ignore node to tip edges)
		timetr_pt_only_internal[[i]] <- timetr_pt[[i]]$edge[ timetr_pt[[i]]$edge[,2] >= (Ntip(timetr_pt[[i]])+1), ]
		internode_lengths[[i]] <- timetr_pt[[i]]$edge.length[ timetr_pt[[i]]$edge[,2] >= (Ntip(timetr_pt[[i]])+1) ]
		
		internode_intervals_df[[i]] <- data.frame(phylotype=i, internode_intervals=internode_lengths[[i]], status=ifelse(i %in% subtype_b_vois_ids_only, yes="VOI",no="Non-VOI"))
	}
}

# tree with tips
plot(merge_para_timetr[[1,4]][[69]]); axisPhylo(root.time=max( node.depth.edgelength(merge_para_timetr[[1,4]][[69]]) ), backward = T)

ii_all_df <- rbindlist(internode_intervals_df)

# extract mean, median, and max internode intervals for each phylotype & VOI status
ii_summary_pt <- ii_all_df %>% group_by(phylotype) %>% summarise(mean=mean(internode_intervals), median=median(internode_intervals), max=max(internode_intervals), status=unique(status))
# convert to months instead of years
ii_all_df$internode_intervals_months <- ii_all_df$internode_intervals * 12
ii_summary_pt2 <- ii_all_df %>% group_by(phylotype) %>% summarise(mean=mean(internode_intervals_months), median=median(internode_intervals_months), max=max(internode_intervals_months), status=unique(status))

res_ii_path <- glue("{RESULTS_PATH}/08_internode_intervals")
system(glue("mkdir -p {res_ii_path}"))

# plot distribution (of means, medians, and max) of the above for VOI vs non-VOI phylotypes
plot_ii_summary <- function(df, var, unit_ii) {
	pl <- ggplot(df, aes(x = !!sym(var), fill = status)) +
		geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
		labs(x = glue("Internode intervals ({unit_ii})"), y = "Count") +
		theme_classic()
	ggsave(plot=pl, file=glue("{res_ii_path}/{var}_distr_{unit_ii}.jpg"), width=8, height=6.5, dpi=300)
}
# years
plot_ii_summary(ii_summary_pt, "mean", "years")
plot_ii_summary(ii_summary_pt, "median", "years")
plot_ii_summary(ii_summary_pt, "max", "years")
# months
plot_ii_summary(ii_summary_pt2, "mean", "months")
plot_ii_summary(ii_summary_pt2, "median", "months")
plot_ii_summary(ii_summary_pt2, "max", "months")

# summary based on grouping by status only across all observations inside each phylotype (months)
ii_summary_voi_st <- ii_all_df %>% group_by(status) %>% summarise(mean=mean(internode_intervals_months), median=median(internode_intervals_months), max=max(internode_intervals_months))
print(ii_summary_voi_st)
# status   mean median   max
# Non-VOI  17.8   6.90  321.
# VOI      18.7   8.28  132.
calc_iqr2 <- function(data) {
	lower_bound <- quantile(data, 0.25, na.rm=T)  # 25th percentile
	upper_bound <- quantile(data, 0.75, na.rm=T)  # 75th percentile
	print(paste0("Lower Bound (25th Percentile):", round(lower_bound,3)))
	print(paste0("Upper Bound (75th Percentile):", round(upper_bound,3)))
}

# to use alongside medians
calc_iqr2(ii_all_df$internode_intervals_months[ii_all_df$status == "Non-VOI"]) #0.002-23.778
calc_iqr2(ii_all_df$internode_intervals_months[ii_all_df$status == "VOI"]) #0.002-26.72

# to use with means
sd(ii_all_df$internode_intervals_months[ii_all_df$status == "Non-VOI"]) #29.0
sd(ii_all_df$internode_intervals_months[ii_all_df$status == "VOI"]) #26.6

# in previous work median transmission intervals were 14 mo for MSM and 27 mo for het
# "2% of transmissions in this population were estimated to occurr within 6 months or less"
# Here: means 18 (non-VOIs) and 19 (VOIs), medians 7 (non-VOIs) and 8 (VOIs), max 321 (~27 y, non-VOIs) and 132 (~11y for VOIs)
# Comparison needs to be against MSM, our means are similar to their medians, but our medians show much lower values
# Anyway values for VOIs seem a bit lower, suggesting faster transmission (but maybe should compare against non-VOIs with matched demog data instead)
# PERHAPS could plot internode intervals (y) over time coloured by VOI status OR ii (y) by phylotype

# 2. Phylogenetic Ornstein-Uhlenbeck Mixed Model (https://github.com/venelin/POUMM)
b_ml_tree <- readRDS(glue("{RDS_PATH}/tree_b_adj2.rds"))# ML tree
vl_b <- readRDS(glue("{RDS_PATH}/vl_subtypes_comb2.rds"))[[1,4]] # mean VLs for treestructure phylotypes
vl_b_vec <- vl_b$mean_log_vl_pat
names(vl_b_vec) <- vl_b$taxon
b_ml_tree_with_vl <- keep.tip(b_ml_tree, tip=names(vl_b_vec))
b_ml_tree_with_vl$edge.length[b_ml_tree_with_vl$edge.length <= 0] <- 1e-6

#all(b_ml_tree_with_vl$tip.label %in% names(vl_b_vec))
#all(names(vl_b_vec) %in% b_ml_tree_with_vl$tip.label)
vl_b_vec_ord <- vl_b_vec[b_ml_tree_with_vl$tip.label]

fit <- POUMM(vl_b_vec_ord, b_ml_tree_with_vl, spec=list(nSamplesMCMC = 5e5))

plot_list <- plot(fit, showUnivarDensityOnDiag = TRUE, doPlot = FALSE)
plot_list$traceplot
plot_list$densplot

summary(fit, mode="short") # make sure G.R. diag stat of H2tMean is close to 1.0 and ESS big
# H2tMean GR: 1.0005123; H2tMean ESS: 614.3031 
AIC(fit) #20316
BIC(fit) #20351
coef(fit)
# alpha     theta     sigma    sigmae 
# 39.996982  4.540723  2.907523  0.707682 
logLik(fit) #-10153.19 (df=5)
summary(fit)["H2tMean"==stat, PostMean]
# 0.1668516
summary(fit)["H2tMean"==stat, unlist(HPD)]
# lower     upper 
# 0.1155296 0.2164108

#fitted(fit)
#plot(resid(fit))
#abline(h=0)

# LRT PMM vs POUMM 
specPMM <- specifyPMM(vl_b_vec_ord[1:length(vl_b_vec_ord)], b_ml_tree_with_vl)
fitPMM <- POUMM(vl_b_vec_ord[1:length(vl_b_vec_ord)], b_ml_tree_with_vl, spec = specPMM, doMCMC=FALSE)
lmtest::lrtest(fitPMM, fit)

# Model 1: fitPMM
# Model 2: fit
# #Df LogLik Df  Chisq Pr(>Chisq)    
# 1   3 -10161                         
# 2   5 -10153  2 16.112  0.0003172 ***
# POUMM outperforms PMM

# When the goal is to estimate H2tMean, it is important to specify an uninformed prior
# for it is the standard uniform distribution
specH2tMean <- specifyPOUMM_ATH2tMeanSeG0(vl_b_vec_ord[1:length(vl_b_vec_ord)], b_ml_tree_with_vl, nSamplesMCMC = 5e5)
specH2tMean$parMapping
specH2tMean$parPriorMCMC
specH2tMean$parLower
specH2tMean$parUpper
fitH2tMean <- POUMM(vl_b_vec_ord[1:length(vl_b_vec_ord)], b_ml_tree_with_vl, spec = specH2tMean)
plot(fitH2tMean, stat = c("H2tMean", "H2e", "H2tInf", "sigmae"), 
					doZoomIn = TRUE, doPlot = TRUE)

# Without setting prior
summary(fit)[stat %in% c("H2tMean", "H2e", "H2tInf", "sigmae")]
# After setting prior
summary(fitH2tMean)[stat %in% c("H2tMean", "H2e", "H2tInf", "sigmae")]
# stat     N       MLE  PostMean                 HPD      ESS     G.R.
# H2tMean  8789 0.1766551 0.1751541 0.1226263,0.2250284 471.8703 1.008458

summary(fitH2tMean, mode="short")

# stat     N           MLE      PostMean                   HPD      ESS     G.R.
# <char> <int>         <num>         <num>                <list>    <num>    <num>
# 	1:        alpha  8789  4.150998e+01  3.799054e+01     11.65091,62.88779 358.6051 1.048854
# 2:        theta  8789  4.541916e+00  4.533195e+00     4.461899,4.597685 514.2835 1.103205
# 3:      H2tMean  8789  1.766551e-01  1.751541e-01   0.1226263,0.2250284 471.8703 1.008458
# 4:       sigmae  8789  7.061503e-01  7.082375e-01   0.6833217,0.7326625 425.3933 1.018878
# 5:           g0  8789  4.986695e+00  4.904086e+00     3.964044,5.920047 900.0000 1.000845
# 6:          H2e  8789  1.680528e-01  1.628470e-01   0.1044086,0.2209526 424.5134 1.019560
# 7:       H2tInf  8789  1.774287e-01  1.807033e-01   0.1364770,0.2322744 616.5435 1.039486
# 8:       H2tMax  8789  1.774281e-01  1.795980e-01   0.1336958,0.2270891 592.7275 1.005152
# 9:        sigma  8789  2.988226e+00  2.843005e+00     1.660110,4.038409 373.7537 1.042285
# 10: sigmaG2tMean  8789  1.069889e-01  1.065564e-01 0.07808788,0.13880326 491.1143 1.006938
# 11:  sigmaG2tMax  8789  1.075581e-01  1.099179e-01 0.08099766,0.13952572 576.1890 1.014037
# 12:  sigmaG2tInf  8789  1.075584e-01  1.110518e-01 0.08072732,0.14212615 632.9163 1.129484
# 13:      logpost  8789            NA -1.016447e+04   -10167.88,-10162.06 485.8072 1.043612
# 14:       loglik  8789 -1.015321e+04            NA                 NA,NA   0.0000       NA
# 15:          AIC  8789  2.031641e+04            NA                 NA,NA   0.0000       NA
# 16:         AICc  8789  2.031642e+04            NA                 NA,NA   0.0000       NA

summary(fitH2tMean)["H2tMean"==stat, PostMean]
# 0.1751541
summary(fitH2tMean)["H2tMean"==stat, unlist(HPD)]
# lower     upper 
# 0.1226263 0.2250284 

lmtest::lrtest(fitPMM, fitH2tMean)
# Model 1: fitPMM
# Model 2: fitH2tMean
# #Df LogLik Df Chisq Pr(>Chisq)    
# 1   3 -10161                        
# 2   5 -10153  2 16.08  0.0003223 ***

# plot trait intervals vs rtt distance
data <- data.table(z = vl_b_vec_ord[1:length(vl_b_vec_ord)], t = nodeTimes(b_ml_tree_with_vl, tipsOnly = TRUE))
data <- data[, group := cut(t, breaks = 5, include.lowest = TRUE)]

ggplot(data = data, aes(x = t, y = z, group = group)) + 
	geom_violin(aes(col = group)) + geom_point(aes(col = group), size=.5)