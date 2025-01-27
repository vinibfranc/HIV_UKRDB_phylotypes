libs_load <- c("mgcv","ggplot2", "glue","lubridate","data.table","readr","scales")
invisible( lapply(libs_load, library, character.only=TRUE) )

DATA_PATH="data"
RDS_PATH="rds"
RESULTS_PATH="results"
dir.create( RESULTS_PATH )

min_cl_size_choices <- c(30, 50, 100) #250,500
tree_names <- c("A_A1","CRF_02_AG","C","B")
extracted_clusters <- readRDS("rds/extracted_clusters.rds")

# same as in 07 begin
demog_md_subtype_match <- readRDS(glue("{RDS_PATH}/demog_md_subtype_match.rds"))
demog_md_subtype_match <- demog_md_subtype_match[ (demog_md_subtype_match$status=="Naïve") & (demog_md_subtype_match$rega3subtype %in% c("A (A1)","CRF 02_AG","C","B")) & (demog_md_subtype_match$exposureid != "Not known"), ]

demog_md_subtype_match$hiv_diag_decimal_date <- decimal_date(as.Date(demog_md_subtype_match$hivpos_ymd))
demog_md_subtype_match <- demog_md_subtype_match %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y))
demog_md_subtype_match <- demog_md_subtype_match %>% mutate(
	age_group = dplyr::case_when(age_diag<=29 ~ "<29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
	age_group = factor(age_group,level = c("<29","30-39","40-49","50-59","60+")))
demog_md_subtype_match$artstart_decimal_date <- decimal_date( as.Date(gsub("\\/", "15", demog_md_subtype_match$artstart_my), "%m%d%Y") )

extr_clust_match_demog <- extr_clust_match_demog_df <- matrix(list(), nrow=length(min_cl_size_choices), ncol=length(tree_names))
for(i in 1:length(min_cl_size_choices)) {
	for(j in 1:length(tree_names)) {
		for(k in 1:length(extracted_clusters[[i,j]][[1]])) {
			print(glue("{min_cl_size_choices[i]}-{tree_names[j]}-PT{k}"))
			# Extract metadata for all patients within each phylotype (backbone included)
			if(!is.null( extracted_clusters[[i,j]][[1]][[k]]$tip.label)) {
				extr_clust_match_demog[[i,j]][[k]] <- demog_md_subtype_match[paste0("t.",demog_md_subtype_match$testindex) %in%  extracted_clusters[[i,j]][[1]][[k]]$tip.label, ]
				extr_clust_match_demog[[i,j]][[k]]$phylotype <- k
			}
		}
		extr_clust_match_demog_df[[i,j]] <- rbindlist(extr_clust_match_demog[[i,j]])
	}
}
# same as in 07 end

# focus here on mcs=30 (selected threshold), subtype B (only one with VOIs)
d_demog <- extr_clust_match_demog_df[[1,4]]

d_demog$patientindex <- as.factor( d_demog$patientindex )
d_demog$phylotype <- as.character(d_demog$phylotype)
d_demog <- d_demog[!is.na(d_demog$hiv_diag_decimal_date),]
d_demog$time <- d_demog$hiv_diag_decimal_date
d_demog <- d_demog[ order( d_demog$time), ]
d_demog <- d_demog[ !duplicated(d_demog$patientindex), ]

fit_gam <- function(d, P){
  # browser()
  d$V <- d$phylotype == P
  d$time <- floor(d$time)
  d <- d[ !is.na( d$time) , ]
  d <- d[ d$time >= min(d$time[d$V]),]
  d <- d[ d$exposureid == 'Homo/bisexual', ] # risk behavior
  
  m = mgcv::gam( V ~ s(time, bs='gp',k=5) , family = binomial(link='logit') , data = d )
  #m = mgcv::gam( V ~ s(time, bs='gp',k=1) , family = binomial(link='logit') , data = d ) #time, bs='cr'

	d$estimated = predict( m )
	estdf = data.frame( time = d$time, estimated_logodds = d$estimated )
	estdf <- estdf[ !duplicated(estdf$time), ]
	estdf <- estdf[ order( estdf$time ) , ]
	estdf <- cbind( estdf, t( sapply( 1:nrow(estdf), function(k){
		.s2 <- d[ d$time == estdf$time[k] , ]
		n <- nrow( .s2 )
		lo= estdf$estimated_logodds[k] 
		f = exp( lo) / ( 1 + exp( lo ))
		ub = qbinom( .975, size = n, prob = f  )  / n
		lb = qbinom( .025, size = n, prob = f  )  /n
		n1 <- sum( .s2$V)
		n2 <- n - n1 
		ef = n1 / n 
		c(lb= log(lb/(1-lb)), ub = log(ub/(1-ub)) , n = n
		  , weights = 1 / sqrt( f*(1-f)/n )
		  , logodds = log(ef / ( 1 - ef )) )
	})) )
        .lb  <- min(estdf$lb[!is.infinite(estdf$lb)])
        .lb <- min( .lb, min(estdf$logodds[!is.infinite(estdf$logodds)]) )
        .ub <- max( max(estdf$ub[!is.infinite(estdf$ub)])
                   ,  max(estdf$logodds[!is.infinite(estdf$logodds)]))
        estdf$lb <- pmax( estdf$lb, .lb )
        estdf$lb <- pmax( estdf$lb,.lb )
        estdf$ub <- pmin( estdf$ub,.ub )
        estdf$ub <- pmin( estdf$ub,.ub )
        estdf$logodds[ is.infinite(estdf$logodds)] <- NA
        estdf$psi <- c(NA, diff( estdf$estimated_logodds ))
 print(tail(estdf))
	pl = ggplot() + 
		geom_point( aes( x = as.POSIXct( round(date_decimal( time )) ), y = logodds, size=n ),  data= estdf)  + 
		theme_minimal() +
		theme( legend.pos='' ) +
		xlab('') + 
		ggtitle( paste0('Frequency of phylotype ', P) ) + 
		geom_path( aes(x= as.POSIXct(round(date_decimal(time))), y=estimated_logodds), data = estdf, color = 'blue' , size = 1) + 
		geom_ribbon( aes(x = as.POSIXct(round(date_decimal(time))), ymin = lb, ymax = ub ), data = estdf , alpha = .25 )
	res = list( result = estdf, m = m, plot = pl	)
 res 
}

# For all phylotypes
lgr_list30 <- list()
for(i in 1:max(as.integer(d_demog$phylotype))) {
	if(i == 6) next
	print(i)
	lgr_list30[[i]] <- fit_gam(d_demog, glue("{i}"))
	system(glue("mkdir -p {RESULTS_PATH}/10_logistic_growth_rates/30-B"))
	suppressMessages( ggsave(file=glue("{RESULTS_PATH}/10_logistic_growth_rates/30-B/{i}_lgr.png"), plot=lgr_list30[[i]]$plot, dpi=600, width=8, height=6, bg="white", limitsize=FALSE) )
}

subtype_b_vois <- c(40,69,133) #5,24,40,137

### BEGIN FIGURE 3B ### 
lgr_vois <- list()
for(k in seq_along( subtype_b_vois )) { #"50.1@B","50.2@B"
	lgr_vois[[k]] <- fit_gam( d_demog, as.character( subtype_b_vois[k]) )
	lgr_vois[[k]]$result$phylotype <- glue( "PT.B.{subtype_b_vois[k]}.UK" )
	lgr_vois[[k]] <- lgr_vois[[k]]$result
}
lgr_vois_all <- rbindlist(lgr_vois)
lgr_vois_all$phylotype <- factor(lgr_vois_all$phylotype, levels = c("PT.B.40.UK","PT.B.69.UK","PT.B.133.UK"))

saveRDS(lgr_vois_all, glue("{RDS_PATH}/lgr_vois_all.rds"))

ne3_pt_pal <- c("PT.B.40.UK"="#56B4E9", "PT.B.69.UK"="#D55E00", "PT.B.133.UK"="#009E73")

leg <- theme(text=element_text(family="Helvetica"), axis.text=element_text(size=10, color="black"), axis.title=element_text(size=10, color="black", face="bold"))

f1c_pl <- ggplot() + geom_point(data= lgr_vois_all, aes( x = as.POSIXct(round(date_decimal( time )) ), y = logodds, size=n, color=phylotype ), alpha=0.5)  + 
	geom_path( aes(x= as.POSIXct(round(date_decimal(time))), y=estimated_logodds, color=phylotype), data = lgr_vois_all, size = 1) + #color = 'blue' 
	scale_size(name="Size of matched\nclusters", limits = c(1,1500), breaks=c(5,10,100,250,500,1000,1500)) + #10,50,100,250,500,1000,1500
	scale_color_manual(values=ne3_pt_pal, name="Phylotype") + scale_fill_manual(values=ne3_pt_pal, name="Phylotype") + #scale_size_manual(name="Sequences") +
	scale_x_datetime(limits=c(as.POSIXct("1995-01-01 00:00:00 CET"),as.POSIXct("2020-01-01 00:00:00 CET")), labels = date_format("%Y", tz = "CET"), 
																		breaks=seq(as.POSIXct("1995-01-01 00:00:00 CET"), as.POSIXct("2020-01-01 00:00:00 CET"), "5 years")) +#breaks = c(1980,1990,2000,2010,2020)
	#scale_y_continuous(limits=c(-10,-3), breaks = c(-8,-7,-6,-5,-4,-3)) +
	coord_cartesian(ylim=c(-8,-4)) + #(-8,-3)
	labs(x="Year",y="Log odds of sampling") + 
	theme_classic() + leg + guides(fill="none", color="none") +
	theme(plot.margin = margin(0.8, 0.8, 0.5, 0.7, "cm"), legend.key.size = unit(1.5, 'cm'), legend.key.height = unit(0.5, 'cm'), 
							legend.key.width = unit(0.5, 'cm'), legend.text=element_text(size=10))
f1c_pl
saveRDS(f1c_pl, glue("{RDS_PATH}/f1c_pl.rds"))
### END FIGURE 3B ###