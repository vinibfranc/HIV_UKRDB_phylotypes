# load and summarise WGS calls from gofasta (closest matches of UKHSA whole-genomes when compared to UKRDB pol region)
libs_load <- c("splitstackshape", "dplyr", "glue", "readr", "stringr", "xlsx", "data.table")
invisible( lapply(libs_load, library, character.only=TRUE) )

load_summarise_muts <- function(mut_csv_file) {
	df_muts <- read.csv(mut_csv_file)
	nseqs <- nrow(df_muts)
	# replace 'nuc' for 'SYNSNP'
	df_muts$mutations <- gsub('nuc:', 'SYNSNP:', df_muts$mutations)
	# remove 'aa' when amino acid mutation
	df_muts$mutations <- gsub('aa:', '', df_muts$mutations)
	# individualise mutations: one mutation per sequence per row
	df_muts <- splitstackshape::cSplit(df_muts, 'mutations', '|', 'long') 
	df_muts$protein <- sub("\\:.*", "", df_muts$mutations)
	
	df_muts$site <- sub('.*:', "", df_muts$mutations)
	df_muts$site <- readr::parse_number(df_muts$site, na="X")
	# regex to match the site numbers (tricky for ins and del mainly)
	regex_ins_del <- "(?<=:)(\\d+)(?=:)"
	df_muts$site <- ifelse(df_muts$protein == "del", stringr::str_extract(df_muts$mutations, regex_ins_del), df_muts$site)
	df_muts$site <- ifelse(df_muts$protein == "ins", stringr::str_extract(df_muts$mutations, regex_ins_del), df_muts$site)
	df_muts$site <- as.integer(df_muts$site)
	
	# will create NAs in these columns for ins and del
	df_muts$anc_site_mut <- sub('.*:', "", df_muts$mutations)
	df_muts$mutsite <- stringr::str_sub(df_muts$mutations, -1)
	df_muts$ancsite <- stringr::str_sub(df_muts$anc_site_mut, start=1, end=1)
	
	# Adjust tat2 and rev2 mutations to be called continuously after tat1 and rev1
	df_muts$mutations <- ifelse(df_muts$protein == "tat2",
																													yes=glue("{df_muts$protein}:{df_muts$ancsite}{(df_muts$site + 72)}{df_muts$mutsite}"),
																													no=df_muts$mutations)
	
	df_muts$mutations <- ifelse(df_muts$protein == "rev2",
																													yes=glue("{df_muts$protein}:{df_muts$ancsite}{(df_muts$site + 25)}{df_muts$mutsite}"),
																													no=df_muts$mutations)
	
	# order by genomic region: non-syn muts, then insertions, deletions, and syn muts
	lvls <- c("gag","pol_pr", "pol_rt", "pol_p15", "pol_int", "vif", "vpr", "tat1", "tat2", "rev1", "rev2", "env_gp120", "env_gp41", "nef", "ins", "del", "SYNSNP")
	df_muts$protein <- factor(df_muts$protein, levels=lvls)
	# order by protein and then site
	df_muts <- df_muts %>% arrange( factor(protein, levels = lvls) , site)
	
	mut_summary <- df_muts %>% group_by(mutations) %>% summarise(n = n(), protein = unique(protein), site = unique(site)) %>% mutate(percent = 100*(n / nseqs)) %>% arrange( factor(protein, levels = lvls) , site) #%>% arrange(desc(n))
	mut_summary_prot <- mut_summary %>% group_by(protein) %>% summarise(n = n())  #%>% arrange(desc(n))
	print("Number of mutations for each protein/feature for all mutations: ")
	print(mut_summary_prot)
	
	# consider defining the ones in >=80% of the sequences
	mut_summary_def <- mut_summary %>% filter(percent >= 80)
	mut_summary_prot_def <- mut_summary_def %>% group_by(protein) %>% summarise(n = n())
	print("Number of mutations for each protein/feature defining mutations: ")
	print(mut_summary_prot_def)
	
	list(indiv_muts = df_muts, mut_summary = mut_summary, mut_summary_def = mut_summary_def, mut_summary_prot_def = mut_summary_prot_def)
}

### 1. HXB2 ref calls ### 
# `mut_calls/` has csv mut calls from gofasta
# IMPORTANT: adding PT137 as well: because NO PT69 sequence and (1) PT137 consensus is closely related to PT69 and (2) both T-test on VL and CD4 analyses suggest it could be a VOI 
# PT137 n=2 could help to narrow down muts by n=1 PT133 and n=16 PT40
csv_calls <- list.files("results/20_gofasta_outputs/ref_hxb2_B", full.names = T) 
print(csv_calls)
muts_ptB <- list()
for(i in 1:length(csv_calls)) {
	print(csv_calls[i])
	muts_ptB[[i]] <- load_summarise_muts(csv_calls[i])
	print("NROW mut_summary")
	print(nrow(muts_ptB[[i]]$mut_summary))
	print("NROW mut_summary_def")
	print(nrow(muts_ptB[[i]]$mut_summary_def))
}

# order: PT133, PT137, PT40 (not close relatives for PT69 unfortunately)
# sequences: 1 (only one seq), 2, 16
# total muts: 756, 737, 3066
# defining (>80% seqs): 756 (as only one all will be defining), 525, 219

View(muts_ptB[[1]]$mut_summary)
View(muts_ptB[[2]]$mut_summary_def)
View(muts_ptB[[3]]$mut_summary_def)

print(t(muts_ptB[[1]]$mut_summary_prot_def))
print(t(muts_ptB[[2]]$mut_summary_prot_def))
print(t(muts_ptB[[3]]$mut_summary_prot_def))

# sanity check: get a few pol_pr, pol_rt and pol_int and paste here: https://hivdb.stanford.edu/hivseq/by-patterns/
print(muts_ptB[[1]]$mut_summary_def$mutations) #percents usually not above 50%: https://hivdb.stanford.edu/hivseq/by-patterns/report/?mutations=IN%3AE11D%2CIN%3AA21T%2CIN%3AD25E%2CIN%3AK111R%2CIN%3AS119P%2CIN%3AK136N%2CIN%3AV165I%2CIN%3AD167E%2CIN%3AV201I%2CIN%3AI208L%2CIN%3AT218S%2CIN%3AD229E%2CIN%3AS230H%2CRT%3AK20R%2CRT%3AV60I%2CRT%3AD123E%2CRT%3AS162C%2CRT%3AI178L%2CRT%3AV179I%2CRT%3AR211K%2CRT%3AD250E%2CRT%3AT286A%2CRT%3AI293V%2CRT%3AP294A%2CRT%3AE297R%2CRT%3AQ334L%2CRT%3AK350R%2CRT%3AM357T%2CRT%3AA376S%2CRT%3AT377S%2CRT%3AK390R%2CRT%3AE404D%2CPR%3AI15V%2CPR%3AK20R%2CPR%3AM36L%2CPR%3AI62V%2CPR%3AL63A&name=IN%3AE11D%2BIN%3AA21T%2BIN%3AD25E%2BIN%3AK111R%2BIN%3AS119P%2BIN%3AK136N%2BIN%3AV165I%2BIN%3AD167E%2BIN%3AV201I%2BIN%3AI208L%2BIN%3AT218S%2BIN%3AD229E%2BIN%3AS230H%2BRT%3AK20R%2BRT%3AV60I%2BRT%3AD123E%2BRT%3AS162C%2BRT%3AI178L%2BRT%3AV179I%2BRT%3AR211K%2BRT%3AD250E%2BRT%3AT286A%2BRT%3AI293V%2BRT%3AP294A%2BRT%3AE297R%2BRT%3AQ334L%2BRT%3AK350R%2BRT%3AM357T%2BRT%3AA376S%2BRT%3AT377S%2BRT%3AK390R%2BRT%3AE404D%2BPR%3AI15V%2BPR%3AK20R%2BPR%3AM36L%2BPR%3AI62V%2BPR%3AL63A
print(muts_ptB[[2]]$mut_summary_def$mutations) #percents usually not above 50%: https://hivdb.stanford.edu/hivseq/by-patterns/report/?mutations=IN%3AE11D%2CIN%3AM50L%2CIN%3AS119P%2CIN%3AT122I%2CIN%3AT125A%2CIN%3AG193E%2CIN%3AV201I%2CIN%3AQ216R%2CIN%3AK219Q%2CIN%3AD256E%2CIN%3AD286N%2CPR%3AG16E%2CPR%3AL33I%2CPR%3AE35D%2CPR%3AN37A%2CPR%3AP39T%2CPR%3AL63P%2CPR%3AI64V%2CPR%3AA71V%2CPR%3AI93L%2CRT%3AD177E%2CRT%3AV245I%2CRT%3AI293V%2CRT%3AI326V%2CRT%3AI329L%2CRT%3AR356K%2CRT%3AT377Q%2CRT%3AK395R%2CRT%3AT403I&name=IN%3AE11D%2BIN%3AM50L%2BIN%3AS119P%2BIN%3AT122I%2BIN%3AT125A%2BIN%3AG193E%2BIN%3AV201I%2BIN%3AQ216R%2BIN%3AK219Q%2BIN%3AD256E%2BIN%3AD286N%2BPR%3AG16E%2BPR%3AL33I%2BPR%3AE35D%2BPR%3AN37A%2BPR%3AP39T%2BPR%3AL63P%2BPR%3AI64V%2BPR%3AA71V%2BPR%3AI93L%2BRT%3AD177E%2BRT%3AV245I%2BRT%3AI293V%2BRT%3AI326V%2BRT%3AI329L%2BRT%3AR356K%2BRT%3AT377Q%2BRT%3AK395R%2BRT%3AT403I
print(muts_ptB[[3]]$mut_summary_def$mutations) #percents usually not above 50%: https://hivdb.stanford.edu/hivseq/by-patterns/report/?mutations=IN%3AE11D%2CIN%3AL101I%2CIN%3AK156N%2CIN%3AV201I%2CPR%3AE35D%2CPR%3AI62V%2CPR%3AA71T%2CPR%3AI93L%2CPR%3AK20R%2CRT%3AR211K%2CRT%3AK390R&name=IN%3AE11D%2BIN%3AL101I%2BIN%3AK156N%2BIN%3AV201I%2BPR%3AE35D%2BPR%3AI62V%2BPR%3AA71T%2BPR%3AI93L%2BPR%3AK20R%2BRT%3AR211K%2BRT%3AK390R

# Get intersection of defining muts ("VOI defining mutations")

# all 2 vois with data
common_muts_vois <- intersect( muts_ptB[[1]]$mut_summary_def$mutations, muts_ptB[[3]]$mut_summary_def$mutations )
length(common_muts_vois) # 100
common_muts_vois_ns <- intersect( muts_ptB[[1]]$mut_summary_def$mutations[muts_ptB[[1]]$mut_summary_def$protein != "SYNSNP"], muts_ptB[[3]]$mut_summary_def$mutations[muts_ptB[[3]]$mut_summary_def$protein != "SYNSNP"])
length(common_muts_vois_ns) # 57

# 2 vois + PT137
common_muts_vois_plus_b137 <- Reduce(intersect, list(muts_ptB[[1]]$mut_summary_def$mutations, muts_ptB[[2]]$mut_summary_def$mutations, muts_ptB[[3]]$mut_summary_def$mutations ))
length(common_muts_vois_plus_b137) # 71
common_muts_vois_plus_b137_ns <- Reduce(intersect, list( muts_ptB[[1]]$mut_summary_def$mutations[muts_ptB[[1]]$mut_summary_def$protein != "SYNSNP"],
																																																									muts_ptB[[2]]$mut_summary_def$mutations[muts_ptB[[1]]$mut_summary_def$protein != "SYNSNP"],
																																																									muts_ptB[[3]]$mut_summary_def$mutations[muts_ptB[[1]]$mut_summary_def$protein != "SYNSNP"]))
length(common_muts_vois_plus_b137_ns) # 41

# diff in non-syn muts when considering PT137 and not
diff_muts <- setdiff( common_muts_vois_ns, common_muts_vois_plus_b137_ns )
# [1] "gag:K30R"        "gag:I94V"        "pol_pr:I62V"     "pol_rt:K20R"     "pol_rt:R211K"    "pol_rt:K390R"    "vif:R63K"        "env_gp120:E5G"   "env_gp120:G167D" "env_gp120:S465N"
# [11] "env_gp41:Q110D"  "env_gp41:E119Q"  "env_gp41:E151A"  "env_gp41:V182I"  "nef:A54D"        "nef:R105K"  
# 2 gag, 1 pol_pr, 3 pol_rt, 1 vif, 7 env, 2 nef REMOVED by adding PT137

common_muts_vois_df <- muts_ptB[[1]]$mut_summary_def %>%
	inner_join(muts_ptB[[3]]$mut_summary_def, by = "mutations")
nrow(common_muts_vois_df)

common_muts_vois_plus_b137_df <- muts_ptB[[1]]$mut_summary_def %>% inner_join(muts_ptB[[2]]$mut_summary_def, by = "mutations") %>% inner_join(muts_ptB[[3]]$mut_summary_def, by = "mutations")
nrow(common_muts_vois_plus_b137_df)

# check if DRMs
drms <- read.delim("data/external/drms_20240326_db_version_20221025.tsv", sep="\t") #, colnames=F
drms$mutations <- paste0(drms$protein,":",drms$ref_site, drms$site, drms$alt_site)
drms$anc_site_mut <- sub('.*:', "", drms$mutations)
drms$mutations <- factor(drms$mutations, levels=unique( drms$mutations[order(drms$protein, as.numeric(drms$site) )] ))

common_muts_vois_plus_b137_df$mutations

common_muts_vois_plus_b137_format_calls <- common_muts_vois_plus_b137_df %>% select(mutations, protein.x, site.x)
colnames(common_muts_vois_plus_b137_format_calls) <- c("mutation", "protein", "site")
common_muts_vois_plus_b137_format_calls <- as.data.frame(common_muts_vois_plus_b137_format_calls)

common_muts_vois_plus_b137_prot <- common_muts_vois_plus_b137_df %>% group_by(protein.x) %>% summarise(n=n())
common_muts_vois_plus_b137_prot
# 1 gag           6
# 2 pol_pr        1
# 3 pol_rt        2
# 4 pol_p15       1
# 5 pol_int       7
# 6 vif           1
# 7 vpr           1
# 8 tat1          2
# 9 rev1          1
# 10 env_gp120     9
# 11 env_gp41      7
# 12 nef           2
# 13 del           1
# 14 SYNSNP       30

# range of muts for each region
common_muts_vois_plus_b137_df_ns <- common_muts_vois_plus_b137_df[common_muts_vois_plus_b137_df$protein.x != "SYNSNP",]
common_muts_vois_plus_b137_df %>% group_by(protein.x) %>% summarise(min_v = min(site.x), max_v = max(site.x))
# protein.x     min_v max_v
# 1 gag         124   467
# 2 pol_pr        3     3
# 3 pol_rt      122   214
# 4 pol_p15      72    72
# 5 pol_int      10   232
# 6 vif         124   124
# 7 vpr          15    15
# 8 tat1         42    59
# 9 rev1         13    13
# 10 env_gp120     8   467
# 11 env_gp41    114   324
# 12 nef          29   124
# 13 del        5771  5771
# 14 SYNSNP     1293  9428

# format table for Data S1
format_calls <- function(df) {
	df <- df %>% dplyr::select(mutations, protein, site, n, percent)
	colnames(df) <- c("mutation", "protein", "site", "n_sequences", "percent_sequences")
	df <- as.data.frame(df)
	df
}

calls_b40 <- format_calls(muts_ptB[[3]]$mut_summary_def)
calls_b133 <- format_calls(muts_ptB[[1]]$mut_summary_def)
calls_b137 <- format_calls(muts_ptB[[2]]$mut_summary_def)

RESULTS_PATH <- "results"
write.xlsx(common_muts_vois_plus_b137_format_calls, file=glue("{RESULTS_PATH}/tables/Data_S1.xlsx"), sheetName="Common muts 2 VOIs and PT.B.137.UK", row.names=FALSE)
write.xlsx(calls_b40, file=glue("{RESULTS_PATH}/tables/Data_S1.xlsx"), sheetName="PT.B.40.UK", append=TRUE, row.names=FALSE)
write.xlsx(calls_b133, file=glue("{RESULTS_PATH}/tables/Data_S1.xlsx"), sheetName="PT.B.133.UK", append=TRUE, row.names=FALSE)
write.xlsx(calls_b137, file=glue("{RESULTS_PATH}/tables/Data_S1.xlsx"), sheetName="PT.B.137.UK", append=TRUE, row.names=FALSE)
# from analysis below: edit DataS1 to mark in bold sites found in both VOIs that overlap VB & CTL escape 

# Match POL drms
match_drms <- function(muts_cohort) {
	seq_muts_drms <- dplyr::inner_join(muts_cohort, drms, by="mutations")
	seq_muts_drms
}

# common 2 VOIs
match_drms(common_muts_vois_df_ns) # none
# common 2 VOIs + PT137
match_drms(common_muts_vois_plus_b137_df_ns) # none
# B133
match_drms(muts_ptB[[1]]$mut_summary_def) # none
# B137
match_drms(muts_ptB[[2]]$mut_summary_def) # none
# B40
match_drms(muts_ptB[[3]]$mut_summary_def) # none


library(readxl)
# check if CTLs
ctls <- read_excel(glue("data/external/table_s2-hlaescapeassociations_Carlston2012.xls"), skip=15)
colnames(ctls)[4] <- "Amino_acid"
colnames(ctls)[9] <- "Phylo_OR"
colnames(ctls)[10] <- "Direct_or_indirect"
colnames(ctls)[19] <- "Escape_Position"
colnames(ctls)[22] <- "Association_Conditions"
ctls %>% group_by(Protein) %>% summarise(min_v = min(Position), max_v = max(Position))
#nef goes until 205 but HXB2 goes up to 124, tat goes until 101 but HXB2 goes up to 87, vpr goes until 93 but HXB2 goes up to 79, vpu is defective in HXB2
# need to convert pol absolute coords to the pr, rt, p15, and int coords I have

unique(ctls$Direction)
# [1] "Adapted"    "NonAdapted"
#Direction: NonAdapted means the amino acid is less likely in the presence of the HLA than in the absence of the HLA. Adapted means it is more likely in the presence of the HLA.
#Phylo_OR: The strength of selection, as measured by the natural log of the odds ratio of escape, given the imputed transmitted/founder virus.
#Direct or Indirect: If the association is significant at q<0.2 when accounting for HIV amino acid covariation, the "Direct" result is reported. Otherwise, the association is "Indirect", and is more likely to be an artifact of covariation.

ctls_filt <- ctls[ctls$QValue <= 0.05,] #1112 rows, 137 HLAs, 336 unique positions
# Keep only Direct evidence
ctls_filt <- ctls_filt[ ctls_filt$Direct_or_indirect == "Direct" ,] #838 rows, 125 HLAs, 276 unique positions
# Split muts in adaptive (the amino acid is enriched in the presence of the noted HLA type, indicated by a positive 
# log odds for escape)  vs nonadaptive (the amino acid is depleted in the presence of the noted HLA type, indicated 
# by a negative log odds for escape)
ctls_filt$adaptive <- ifelse(ctls_filt$Phylo_OR > 0, "yes", "no")
table(ctls_filt$adaptive) #434 no, 404 yes

# pol positions ranges from 10 to 947
ctls_filt_adjcoords <- ctls_filt %>% mutate(mut_adj_coord = case_when( 
	Protein == "pol" & Position <= 99 ~ glue("pol_pr:{Consensus}{Position}{Amino_acid}"),
	Protein == "pol" & Position > 100 & Position <= 539 ~  glue("pol_rt:{Consensus}{Position-99}{Amino_acid}"),
	Protein == "pol" & Position > 540 & Position <= 659 ~ glue("pol_p15:{Consensus}{Position-539}{Amino_acid}"),
	Protein == "pol" & Position > 660 & Position <= 947 ~ glue("pol_int:{Consensus}{Position-659}{Amino_acid}"),
	TRUE ~ glue("{Protein}:{Consensus}{Position}{Amino_acid}")))
ctls_filt_adjcoords <- ctls_filt_adjcoords %>% mutate(Protein = case_when( 
	Protein == "pol" & Position <= 99 ~ "pol_pr",
	Protein == "pol" & Position > 100 & Position <= 539 ~ "pol_rt",
	Protein == "pol" & Position > 540 & Position <= 659 ~ "pol_p15",
	Protein == "pol" & Position > 660 & Position <= 947 ~ "pol_int",
	TRUE ~ glue("{Protein}")))

# Add "env_" to gp120 and gp41
ctls_filt_adjcoords$Protein[ctls_filt_adjcoords$Protein == "gp120"] <- "env_gp120"
ctls_filt_adjcoords$Protein[ctls_filt_adjcoords$Protein == "gp41"] <- "env_gp41" 
ctls_filt_adjcoords <- ctls_filt_adjcoords %>% mutate(mut_adj_coord = case_when( 
	Protein == "env_gp120" ~ glue("env_gp120:{Consensus}{Position}{Amino_acid}"),
	Protein == "env_gp41" ~ glue("env_gp41:{Consensus}{Position}{Amino_acid}"),
	TRUE ~ mut_adj_coord))
ctls_filt_adjcoords$new_coord <- sub('.*:', "", ctls_filt_adjcoords$mut_adj_coord)
ctls_filt_adjcoords$new_coord <- readr::parse_number(ctls_filt_adjcoords$new_coord, na="X")
ctls_filt_adjcoords$new_coord <- as.integer(ctls_filt_adjcoords$new_coord)
ctls_filt_adjcoords %>% group_by(Protein) %>% summarise(min_v = min(new_coord), max_v = max(new_coord))

# Matches with CTL escapes
match_ctls <- function(muts_cohort) {
	seq_muts_ctls <- dplyr::inner_join(muts_cohort, ctls_filt_adjcoords, by=c("mutations"="mut_adj_coord"), multiple = "first")
	seq_muts_ctls
}

# common 2 VOIs 
ctl_match_2v <- match_ctls(common_muts_vois_df_ns) #3 matches: gag:K30R,pol_pr:I62V,pol_int:E11D
print(ctl_match_2v$Phylo_OR) #10.657  5.677 10.973
# common 2 VOIs + PT.B.137.UK
ctl_match <- match_ctls(common_muts_vois_plus_b137_df_ns) #1 match: pol_int:E11D
print(ctl_match$Phylo_OR) #10.973 (not very big when compared to other int CTL escapes)
# B133
ctl_match_b133 <- match_ctls(muts_ptB[[1]]$mut_summary_def) # 32 matches, mostly positive
hist(ctl_match_b133$Phylo_OR, breaks=32)
# B137
ctl_match_b137 <- match_ctls(muts_ptB[[2]]$mut_summary_def) # 26 matches, mostly positive
hist(ctl_match_b137$Phylo_OR, breaks=26)
# B40
ctl_match_b40 <- match_ctls(muts_ptB[[3]]$mut_summary_def) # 7 matches, mostly positive
hist(ctl_match_b40$Phylo_OR, breaks=7)

# Check if overlapping muts with VB clade
vb_muts <- read.csv("data/external/science.abk1688_data_s2_vb_muts.csv", header=T) #275 muts
vb_muts$mutations <- paste0(vb_muts$gene,":",vb_muts$aa_hxb2,vb_muts$gene_position_hxb2,vb_muts$aa_vb)
vb_muts %>% group_by(gene) %>% summarise(min_v = min(gene_position_hxb2), max_v = max(gene_position_hxb2))

# Convert env, pol, tat and rev absolute aa coordinates to gp_120, gp_41, pol_pr, pol_rt, pol_p15, and pol_int

# pol
pol_start <- 56
vb_muts_adj <- vb_muts %>% mutate(mut_adj_coord = case_when( 
	gene == "pol" & gene_position_hxb2 > pol_start & gene_position_hxb2 <= (99 + pol_start) ~ glue("pol_pr:{aa_hxb2}{gene_position_hxb2-pol_start}{aa_vb}"),
	gene == "pol" & gene_position_hxb2 > (100 + pol_start) & gene_position_hxb2 <= (539 + pol_start) ~ glue("pol_rt:{aa_hxb2}{gene_position_hxb2-99-pol_start}{aa_vb}"),
	gene == "pol" & gene_position_hxb2 > (540 + pol_start) & gene_position_hxb2 <= (659 + pol_start) ~ glue("pol_p15:{aa_hxb2}{gene_position_hxb2-539-pol_start}{aa_vb}"),
	gene == "pol" & gene_position_hxb2 > (540 + pol_start) & gene_position_hxb2 <= (947 + pol_start) ~ glue("pol_int:{aa_hxb2}{gene_position_hxb2-659-pol_start}{aa_vb}"),
	TRUE ~ glue("{gene}:{aa_hxb2}{gene_position_hxb2}{aa_vb}")))
vb_muts_adj <- vb_muts_adj %>% mutate(gene = case_when( 
	gene == "pol" & gene_position_hxb2 > pol_start & gene_position_hxb2 <= (99 + pol_start) ~ "pol_pr",
	gene == "pol" & gene_position_hxb2 > (100 + pol_start) & gene_position_hxb2 <= (539 + pol_start) ~ "pol_rt",
	gene == "pol" & gene_position_hxb2 > (540 + pol_start) & gene_position_hxb2 <= (659 + pol_start) ~ "pol_p15",
	gene == "pol" & gene_position_hxb2 > (540 + pol_start) & gene_position_hxb2 <= (947 + pol_start) ~ "pol_int",
	TRUE ~ glue("{gene}")))

vb_muts_adj$new_coord <- sub('.*:', "", vb_muts_adj$mut_adj_coord)
vb_muts_adj$new_coord <- readr::parse_number(vb_muts_adj$new_coord, na="X"); vb_muts_adj$new_coord <- as.integer(vb_muts_adj$new_coord)
vb_muts_adj <- vb_muts_adj[vb_muts_adj$new_coord >= 0,]
vb_muts_adj %>% group_by(gene) %>% summarise(min_v = min(new_coord), max_v = max(new_coord))

# env
vb_muts_adj <- vb_muts_adj %>% mutate(mut_adj_coord = case_when( 
	gene == "env" & gene_position_hxb2 <= 511 ~ glue("env_gp120:{aa_hxb2}{gene_position_hxb2}{aa_vb}"),
	gene == "env" & gene_position_hxb2 > 512 & gene_position_hxb2 <= 856 ~ glue("env_gp41:{aa_hxb2}{gene_position_hxb2-511}{aa_vb}"),
	TRUE ~ glue("{mut_adj_coord}")))
vb_muts_adj <- vb_muts_adj %>% mutate(gene = case_when( 
	gene == "env" & gene_position_hxb2 <= 511 ~ "env_gp120",
	gene == "env" & gene_position_hxb2 > 512 & gene_position_hxb2 <= 856 ~ "env_gp41",
	TRUE ~ glue("{gene}")))

vb_muts_adj$new_coord <- sub('.*:', "", vb_muts_adj$mut_adj_coord)
vb_muts_adj$new_coord <- readr::parse_number(vb_muts_adj$new_coord, na="X"); vb_muts_adj$new_coord <- as.integer(vb_muts_adj$new_coord)
vb_muts_adj <- vb_muts_adj[vb_muts_adj$new_coord >= 0,]
vb_muts_adj %>% group_by(gene) %>% summarise(min_v = min(new_coord), max_v = max(new_coord))

# tat and rev
vb_muts_adj <- vb_muts_adj %>% mutate(mut_adj_coord = case_when( 
	gene == "tat" & gene_position_hxb2 <= 72 ~ glue("tat1:{aa_hxb2}{gene_position_hxb2}{aa_vb}"),
	gene == "tat" & gene_position_hxb2 > 72 ~ glue("tat2:{aa_hxb2}{gene_position_hxb2+72}{aa_vb}"),
	gene == "rev" & gene_position_hxb2 <= 25 ~ glue("rev1:{aa_hxb2}{gene_position_hxb2}{aa_vb}"),
	gene == "rev" & gene_position_hxb2 > 25 ~ glue("rev2:{aa_hxb2}{gene_position_hxb2+25}{aa_vb}"),
	TRUE ~ glue("{mut_adj_coord}")))
vb_muts_adj <- vb_muts_adj %>% mutate(gene = case_when( 
	gene == "tat" & gene_position_hxb2 <= 72 ~ "tat1",
	gene == "tat" & gene_position_hxb2 > 72 ~ "tat2",
	gene == "rev" & gene_position_hxb2 <= 25 ~ "rev1",
	gene == "rev" & gene_position_hxb2 > 25 ~ "rev2",
	TRUE ~ glue("{gene}")))

vb_muts_adj$new_coord <- sub('.*:', "", vb_muts_adj$mut_adj_coord)
vb_muts_adj$new_coord <- readr::parse_number(vb_muts_adj$new_coord, na="X"); vb_muts_adj$new_coord <- as.integer(vb_muts_adj$new_coord)
vb_muts_adj <- vb_muts_adj[vb_muts_adj$new_coord >= 0,]
vb_muts_adj %>% group_by(gene) %>% summarise(min_v = min(new_coord), max_v = max(new_coord))

match_vb <- function(muts_cohort) {
	seq_muts_vb <- dplyr::inner_join(muts_cohort, vb_muts_adj, by=c("mutations"="mut_adj_coord"))
	seq_muts_vb
}

# 2 VOIs and PT.B.137.UK
vb_match_2v <- match_vb(common_muts_vois_df_ns) # 1 match, pol_int:V201I
vb_match_2v$mutations 
# 2 VOIs and PT.B.137.UK
vb_match <- match_vb(common_muts_vois_plus_b137_df_ns) # 1 match, pol_int:V201I
vb_match$mutations
# vpr:H15F in VB and vpr:H15Y in UKRDB VOIs
# deletion in UKRDB VOIs is in nt 5771 at codon 77 of vpr (one nt before frameshift in HXB2 relative to other HIV-1, so it does not have the frameshift) and not found in VB
# B133
vb_match_b133 <- match_vb(muts_ptB[[1]]$mut_summary_def) # 30 matches
vb_match_b133$mutations
# B137
vb_match_b137 <- match_vb(muts_ptB[[2]]$mut_summary_def) # 25 matches
vb_match_b137$mutations
# B40
vb_match_b40 <- match_vb(muts_ptB[[3]]$mut_summary_def) # 6 matches (none particurly rare)
vb_match_b40$mutations

# Find diff muts to same site between VOIs and VB
match_vb_diff_muts_same_site <- function(muts_cohort, prot_col, site_col, muts1, muts2) {
	same_site_diff <- merge(muts_cohort, vb_muts_adj, by.x = c(prot_col, site_col), by.y = c("gene", "new_coord"))
	same_site_diff <- same_site_diff[same_site_diff[[muts1]] != same_site_diff[[muts2]] & same_site_diff$aa_hxb2 != same_site_diff$aa_vb,]
	same_site_diff
}

# 2 VOIs
vois2_ssn <- match_vb_diff_muts_same_site(common_muts_vois_df_ns, "protein.x", "site.x", "mutations.x", "mut_adj_coord")
vois2_ssn$mutations.x; vois2_ssn$mut_adj_coord # "nef:R105K/Q"    "pol_rt:R211K/G" "vpr:H15Y/F"

# 2 VOIs and PT.B.137.UK
vois_ssn <- match_vb_diff_muts_same_site(common_muts_vois_plus_b137_df_ns, "protein.x", "site.x", "mutations.x", "mut_adj_coord")
vois_ssn$mutations.x; vois_ssn$mut_adj_coord # 1 diff muts at same site (vpr); "vpr:H15Y/F"

# B133
b133_ssn <- match_vb_diff_muts_same_site(muts_ptB[[1]]$mut_summary_def, "protein", "site", "mutations.x", "mut_adj_coord")
b133_ssn$mutations.x; b133_ssn$mut_adj_coord # 70 (too many because only one seq)

# B137
b137_ssn <- match_vb_diff_muts_same_site(muts_ptB[[2]]$mut_summary_def, "protein", "site", "mutations.x", "mut_adj_coord")
b137_ssn$mutations.x; b137_ssn$mut_adj_coord # 27 (too many because only 2 seqs)

# B40
b40_ssn <- match_vb_diff_muts_same_site(muts_ptB[[3]]$mut_summary_def, "protein", "site", "mutations.x", "mut_adj_coord")
b40_ssn$mutations.x; b40_ssn$mut_adj_coord # 12 (7 env_gp120, 1 gag, 1 nef, 1 pol_rt, 2 vpr)

# Find adjacent mutations in VOIs and VB
find_adjacent_mutations <- function(muts_cohort, prot_col, site_col) {
	join_dfs <- rbind(muts_cohort, vb_muts_adj_reduced)
	join_dfs <- join_dfs[join_dfs$aa_hxb2 != join_dfs$aa_vb | is.na(join_dfs$aa_hxb2),]
	adjacent <- join_dfs %>% group_by(!!sym(prot_col)) %>%
		mutate(
			before_dataset = dataset[match(site - 1, site)],
			after_dataset = dataset[match(site + 1, site)]
		) %>%
		filter(
			(site - 1 %in% site & dataset != before_dataset) |
				(site + 1 %in% site & dataset != after_dataset)
		) %>%
		select(-before_dataset, -after_dataset)
	View(adjacent)
	# print(adjacent)
	adjacent
}

# adjust dfs to merge above
vb_muts_adj_reduced <- vb_muts_adj %>% dplyr::select(mut_adj_coord, new_coord, gene, aa_hxb2, aa_vb)
vb_muts_adj_reduced$dataset <- "VB"
colnames(vb_muts_adj_reduced) <- c("mutations", "site", "protein", "aa_hxb2", "aa_vb", "dataset")

# 2 VOIs and PT.B.137.UK
common_muts_vois_ns_reduced <- common_muts_vois_plus_b137_df_ns %>% dplyr::select(mutations, site.x, protein.x)
common_muts_vois_ns_reduced$aa_hxb2 <- NA; common_muts_vois_ns_reduced$aa_vb <- NA
common_muts_vois_ns_reduced$dataset <- "2 VOIs and PT137"
colnames(common_muts_vois_ns_reduced) <- c("mutations", "site", "protein", "aa_hxb2", "aa_vb", "dataset")

vois_ad <- find_adjacent_mutations(common_muts_vois_ns_reduced, "protein", "site")
vois_ad$mutations # 1 mut in env_gp120, 5 in env_gp41, 1 in vif

# Manually get 8-mer motifs of 28 non-syn muts and search here https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html
prev_files <- list.files(path = "results/20_prev_common_muts", pattern = "*.txt", full.names = T)

df_prev <- df_summ <- list()
for(i in 1:length(prev_files)) {
	print("=========")
	print(i)
	df_prev[[i]] <- read.delim(prev_files[i], skip=2, header=T, sep="\t")
	print(nrow(df_prev[[i]]))
	df_prev[[i]] <- df_prev[[i]] %>% distinct(Name, .keep_all = T)
	df_prev[[i]] <- df_prev[[i]] %>% mutate(mutation = str_extract(prev_files[i], "(?<=_)[^/]+(?=\\.txt)"))
	print(nrow(df_prev[[i]]))
	
	df_summ[[i]] <- df_prev[[i]] %>% summarise(mutation = unique(mutation), freq_B = sum(Subtype == 'B'), freq_incl_other_subtypes = n(), percentage_B = 100 * (freq_B / freq_incl_other_subtypes))
	print(df_summ[[i]])
}

df_summ_all <- rbindlist(df_summ)
write.csv(df_summ_all, "results/20_prev_common_muts/01_aux_freqs.csv", quote=F, row.names = F)

# Add column to Data S1 to say if mut fixed (part of consensus of ukhsa seqs) or not 
# and report results accordingly 
REF_SEQS_PATH="data/refs/"
file_cons_vs_hxb2_muts <- glue("{REF_SEQS_PATH}/cons_ukhsa_muts_vs_hxb2.csv")
muts_cons_vs_hxb2 <- load_summarise_muts(file_cons_vs_hxb2_muts)
print(nrow(muts_cons_vs_hxb2$mut_summary)) #273
print(nrow(muts_cons_vs_hxb2$mut_summary_def)) #273
View(muts_cons_vs_hxb2$mut_summary_def)

# Code below helps to manually flag consensus and non-consensus muts (Data S1 columns "present in COMPARE-HIV/INITiO consensus sequence?")

#common_muts_vois_not_in_cons <- base::setdiff(common_muts_vois, muts_cons_vs_hxb2$mut_summary_def$mutations)
#common_muts_vois_not_in_cons # 13

# sheet 1: all VOIs combined
common_muts_vois_plus_b137_not_in_cons <- base::setdiff(common_muts_vois_plus_b137, muts_cons_vs_hxb2$mut_summary_def$mutations) # 3
common_muts_vois_plus_b137_not_in_cons # 3

# sheet 2: PT40
calls_b40_not_in_cons <- base::setdiff(calls_b40$mutation, muts_cons_vs_hxb2$mut_summary_def$mutations)
calls_b40_not_in_cons #95

# sheet 3: PT133
calls_b133_not_in_cons <- base::setdiff(calls_b133$mutation, muts_cons_vs_hxb2$mut_summary_def$mutations)
calls_b133_not_in_cons #607
calls_b133_IN_cons <- intersect(calls_b133$mutation, muts_cons_vs_hxb2$mut_summary_def$mutations)
calls_b133_IN_cons #149

# sheet 4: PT137
calls_b137_not_in_cons <- base::setdiff(calls_b137$mutation, muts_cons_vs_hxb2$mut_summary_def$mutations)
calls_b137_not_in_cons #366
calls_b137_IN_cons <- intersect(calls_b137$mutation, muts_cons_vs_hxb2$mut_summary_def$mutations)
calls_b137_IN_cons #159
