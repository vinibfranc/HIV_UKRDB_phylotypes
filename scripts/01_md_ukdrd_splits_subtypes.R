libs_load <- c("dplyr","openxlsx", "glue","ape","lubridate","ggtree","ggpubr","ggplot2","viridis")
invisible( lapply(libs_load, library, character.only=TRUE) )

### PROCESS UK DR DATABASE DATA ###
demog_md <- read.csv("data/demographics.csv", header=T)
cd4_md <- read.csv("data/cd4s.csv", header=T)
vl_md <- read.csv("data/vloads.csv", header=T)
subtype_md <- read.csv("data/resistance.csv", header=T)

## LOAD VARIABLE DICTIONARIES AND REPLACE CODES
path_lkp <- "data/lookup_variables_fixed.xlsx"
lookup <- openxlsx::getSheetNames(path_lkp)
lookup_dfs <- lapply(lookup, openxlsx::read.xlsx, xlsxFile=path_lkp)
names(lookup_dfs) <- lookup

demog_md$sexid <- factor(demog_md$sexid, levels=lookup_dfs$sexid$sexid, labels=lookup_dfs$sexid$label)
demog_md$ethnicityid <- factor(demog_md$ethnicityid, levels=lookup_dfs$ethnicityid$ethnicityid, labels=lookup_dfs$ethnicityid$label)
demog_md$exposureid <- factor(demog_md$exposureid, levels=lookup_dfs$exposureid$exposureid, labels=lookup_dfs$exposureid$label)
demog_md$PHE_regiondiagnosed <- factor(demog_md$PHE_regiondiagnosed, levels=lookup_dfs$PHE_regiondiagnosed$PHE_regiondiagnosed, labels=lookup_dfs$PHE_regiondiagnosed$label)
demog_md$PHE_regiondiagnosed <- as.character(demog_md$PHE_regiondiagnosed)
demog_md %>% group_by(PHE_regiondiagnosed) %>% summarise(n=n()) %>% dplyr::arrange(desc(n))
demog_md <- demog_md %>% mutate(PHE_regiondiagnosed = case_when(
	(PHE_regiondiagnosed=="Northern Ireland") | (PHE_regiondiagnosed=="Scotland") | (PHE_regiondiagnosed=="Wales") ~ "Northern Ireland, Scotland, and Wales",
	(PHE_regiondiagnosed=="North East") | (PHE_regiondiagnosed=="North West") | (PHE_regiondiagnosed=="Yorkshire and Humber") ~ "North of England",
	(PHE_regiondiagnosed=="West Midlands") | (PHE_regiondiagnosed=="East Midlands") | (PHE_regiondiagnosed=="East of England") ~ "Midlands and East of England",
	(PHE_regiondiagnosed=="South West") | (PHE_regiondiagnosed=="South Central") | (PHE_regiondiagnosed=="South East") ~ "South of England",
	PHE_regiondiagnosed=="London" ~ "London", PHE_regiondiagnosed=="Other/unknown" ~ "Not Known", TRUE ~ "Not Known"))
demog_md %>% group_by(PHE_regiondiagnosed) %>% summarise(n=n()) %>% dplyr::arrange(desc(n))
# nrow(demog_md[demog_md$PHE_regiondiagnosed == "South of England",]) "North of England", "London", "Midlands and East of England", "Northern Ireland, Scotland, and Wales"
#demog_md$PHE_regiondiagnosed <- factor(demog_md$PHE_regiondiagnosed, levels=unique(demog_md$PHE_regiondiagnosed)) #labels=unique(demog_md$PHE_regiondiagnosed)
demog_md$hivpos_ymd <- as.Date(gsub("\\/", "15", demog_md$hivpos_my), "%m%d%Y")
demog_md$hivpos_year <- year(demog_md$hivpos_ymd)
system("mkdir -p rds/")
saveRDS(demog_md, "rds/demog_md.rds")

### STEP 2: INSPECT SUBTYPES AND SPLIT DFS AND FASTAS FOR EACH SUBTYPE ###
unique(subtype_md$rega3subtype) # 396 unique subtypes; 144 subtypes LANL
unique(subtype_md$rega3simplesubtype) # 44 unique subtypes

#subtype_md_nseqs <- subtype_md %>% add_count(patientindex, name="nseqs") %>% group_by(patientindex) %>% filter(n() == 1)

# Deduplicating only for exploratory plotting purposes
subtype_md_dedup <- subtype_md[!duplicated(subtype_md$patientindex), ] # nrow = total number of patients
subtype_md_counts <- subtype_md_dedup %>% group_by(rega3subtype) %>% summarise(n = n()) %>% arrange(desc(n))
# Remove simple subtypes with < 1000 patientindex represented
subtype_md_counts <- subtype_md_counts %>% filter(n >= 1000)
# Merge counts with patient/subtype md
subtype_md_flt <- merge(subtype_md_dedup, subtype_md_counts, by="rega3subtype")

subtype_md$status <- factor(subtype_md$status, levels=lookup_dfs$status$status, labels=lookup_dfs$status$label)

# Plot for each of the 14 subtypes with >100 patients the overall and time-stratified distributions of 
# (1) sex, (2) ethnicity, (3) exposure, (4) region_diag
demog_md_subtype_match <- demog_md %>% inner_join(subtype_md_flt, by="patientindex", multiple="first")

# Table S1
# Gender
demog_md_subtype_match_sex <- demog_md_subtype_match
demog_md_subtype_match_sex$subtype2 <- ifelse(demog_md_subtype_match_sex$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_sex$rega3subtype), "Others")
View(demog_md_subtype_match_sex %>% group_by(sexid, subtype2) %>% summarise(n=n())) #%>% mutate(n_perc=glue("{n} ({n*100/})"))

# Ethnicity
demog_md_subtype_match_eth <- demog_md_subtype_match
demog_md_subtype_match_eth$subtype2 <- ifelse(demog_md_subtype_match_eth$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_eth$rega3subtype), "Others")
demog_md_subtype_match_eth <- demog_md_subtype_match_eth %>% mutate(ethnicityid2 = case_when(
	(ethnicityid=="Black-Caribbean") | (ethnicityid=="Black-African") | (ethnicityid=="Black-other/unspecified") ~ "Black-Caribbean / African / other",
	(ethnicityid=="Indian/Pakistani/Bangladeshi") ~ "Indian/Pakistani/Bangladeshi", 
	(ethnicityid=="Other Asian/Oriental") ~ "Other Asian/Oriental", 
	(ethnicityid=="White") ~ "White",
	TRUE ~ "Other/mixed/NA"))
View(demog_md_subtype_match_eth %>% group_by(ethnicityid2, subtype2) %>% summarise(n=n()))

# Exposure
demog_md_subtype_match_exposure <- demog_md_subtype_match
demog_md_subtype_match_exposure$subtype2 <- ifelse(demog_md_subtype_match_exposure$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_exposure$rega3subtype), "Others")
View(demog_md_subtype_match_exposure %>% group_by(exposureid, rega3subtype) %>% summarise(n=n()))
demog_md_subtype_match_exposure <- demog_md_subtype_match_exposure %>% mutate(exposureid2 = case_when(
	(exposureid=="IDU") | (exposureid=="Blood products") ~ "IDU and blood products",
	(exposureid=="Homo/bisexual") ~ "Homo/bisexual", (exposureid=="Heterosexual") ~ "Heterosexual", TRUE ~ "Other/NA")) #(exposureid=="Not known") ~ "Not known"
exp1 <- demog_md_subtype_match_exposure %>% group_by(exposureid2, subtype2) %>% summarise(n=n())
exp1.1 <-  demog_md_subtype_match_exposure %>% group_by(subtype2) %>% summarise(n=n())
View(exp1)

# Region
demog_md_subtype_match_region <- demog_md_subtype_match
demog_md_subtype_match_region$subtype2 <- ifelse(demog_md_subtype_match_region$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_region$rega3subtype), "Others")
View(demog_md_subtype_match_region %>% group_by(PHE_regiondiagnosed, subtype2) %>% summarise(n=n()))

# Age at diagnosis
demog_md_subtype_match_age <- demog_md_subtype_match
demog_md_subtype_match_age$subtype2 <- ifelse(demog_md_subtype_match_age$rega3subtype %in% c("B","C","A (A1)", "CRF 02_AG"), as.character(demog_md_subtype_match_age$rega3subtype), "Others")
demog_md_subtype_match_age$hiv_diag_decimal_date <- decimal_date(as.Date(demog_md_subtype_match_age$hivpos_ymd))
demog_md_subtype_match_age <- demog_md_subtype_match_age %>% mutate(age_diag=round(hiv_diag_decimal_date-dob_y))
demog_md_subtype_match_age <- demog_md_subtype_match_age %>% mutate(
	age_group = dplyr::case_when(age_diag<=29 ~ "<29", age_diag>29 & age_diag<=39 ~ "30-39", age_diag>39 & age_diag<=49 ~ "40-49", age_diag>49 & age_diag<=59 ~ "50-59", age_diag>59 ~ "60+"),
	age_group = factor(age_group,level = c("<29","30-39","40-49","50-59","60+")))
View(demog_md_subtype_match_age %>% group_by(age_group, subtype2) %>% summarise(n=n()))

# Plots
demog_md_subtype_match <- demog_md_subtype_match %>% mutate(hivpos_year_bins = case_when(
	(hivpos_year >= 1980) & (hivpos_year <= 2000) ~ "1980-2000",(hivpos_year > 2000) & (hivpos_year <= 2010) ~ "2001-2010",
	(hivpos_year > 2010) & (hivpos_year <= 2020) ~ "2011-2020"))
#demog_md_subtype_match <- demog_md_subtype_match[!is.na(demog_md_subtype_match$hivpos_year_bins)]
demog_vars_pl <- c("sexid", "ethnicityid", "exposureid", "PHE_regiondiagnosed") #"status"
#subtype_pl <- c(rep("rega3subtype", 4))
demog_md_subtype_match_vars <- list()
system("mkdir -p results/01_exploratory_plots/")
for(i in 1:length(demog_vars_pl)) {
	demog_md_subtype_match_vars[[i]] <- demog_md_subtype_match %>% group_by_at(vars(demog_vars_pl[i], rega3subtype, hivpos_year_bins)) %>% summarise(n=n()) #subtype_pl[i])
	p <- ggplot(demog_md_subtype_match_vars[[i]], aes(fill=hivpos_year_bins, x= .data[[demog_vars_pl[i]]], y=n)) + #group=!!demog_vars_pl[i]
		geom_bar(position="stack", stat="identity") + theme(axis.text.y = element_text(size=8,color="black"), axis.text.x=element_text(color="black", angle=90, vjust=0.5, hjust=1)) + #+ theme_minimal() +
		labs(y="Count", x=demog_vars_pl[i]) + facet_wrap(~rega3subtype, scales="free") # scale_fill_manual(values=pal, name="Year bins")
	ggsave(glue("results/01_exploratory_plots/{demog_vars_pl[i]}.pdf"), plot=p, width=15, height=20, dpi=600, bg="white")
}

# Get only the >1000 freq subtypes from subtype_md_flt but keeping multiple tests / samples per patient
subtype_md <- subtype_md[subtype_md$rega3subtype %in% unique(subtype_md_flt$rega3subtype),] # 9 subtypes, 4 with >5k sequences

# Add day 15 to all and convert to yyyy-mm-dd format
subtype_md$dbsample_ymd <- as.Date(gsub("\\-", "15", subtype_md$dbsample_my), "%b%d%y")
subtype_md$dbsample_date <- decimal_date(subtype_md$dbsample_ymd)
subtype_md <- subtype_md %>% dplyr::select(rega3subtype, patientindex, testindex, dbsample_date, status, genflag) #testreasonid,n
incl_subtypes <- c("A (A1)","B","C","CRF 02_AG") #names(subtype_md_msm_naive_uq_dates_list)
subtype_md <- subtype_md[subtype_md$rega3subtype %in% incl_subtypes,] #118401 rows, 80049 unique patients

hist(subtype_md$dbsample_date)

saveRDS(subtype_md, "rds/subtype_md.rds")
# Split dfs for individual subtypes (ALL DATA)
subtype_df_list <- split(subtype_md, subtype_md$rega3subtype)

subtype_md %>% group_by(rega3subtype) %>% summarise(n=n())
subtype_md %>% group_by(status) %>% summarise(n=n())

demog_md_subtype_match_mult_sameid <- demog_md %>% inner_join(subtype_md, by="patientindex", multiple="all") # join in risk group (exposureid)
View(demog_md_subtype_match_mult_sameid %>% group_by(exposureid,rega3subtype,status) %>% summarise(n=n()))

# Get all risk groups (except 'Not known') and naive (not treated)
demog_md_subtype_match <- demog_md_subtype_match_mult_sameid[demog_md_subtype_match_mult_sameid$exposureid != "Not known",] # 118.4k to 100.5k rows
saveRDS(demog_md_subtype_match, "rds/demog_md_subtype_match.rds")
demog_md_subtype_match_naive <- demog_md_subtype_match %>% filter(status=="Naïve") %>% group_by(patientindex) %>% arrange(dbsample_date) %>% filter(row_number()==1) #unique patients: 46745

demog_md_subtype_match_naive %>% distinct(patientindex, .keep_all = T) %>% group_by(rega3subtype) %>% summarise(n=n())
# 1 A (A1)        3133
# 2 B            27938
# 3 C            12547
# 4 CRF 02_AG     3127

demog_md_subtype_match_naive_uq_dates <- demog_md_subtype_match_naive %>% dplyr::select(testindex, dbsample_date, rega3subtype, patientindex, exposureid, status)
demog_md_subtype_match_naive_uq_dates_list <- split(demog_md_subtype_match_naive_uq_dates, demog_md_subtype_match_naive_uq_dates$rega3subtype)
saveRDS(demog_md_subtype_match_naive_uq_dates_list, "rds/demog_md_subtype_match_naive_uq_dates_list.rds")
# nrows: 3133, 27938, 12547, 3127

# Read big fasta with all subtypes and write file with inputted subtype sequences only
match_seqs_subtype_md <- function(fasta_path, list_interest, subtype, seqs_folder) {
	fasta_hiv = read.FASTA(fasta_path, type = "DNA")
	fasta_hiv_filt <- fasta_hiv[( names(fasta_hiv) %in% list_interest[[subtype]]$testindex )] #subtype_md_flt_msm_uq_dates_list[[subtype]]$patientindex
	system(glue("mkdir -p {seqs_folder}"))
	#system(paste0("mkdir -p ",seqs_folder))
	write.FASTA(fasta_hiv_filt, glue("{seqs_folder}/{subtype}.fasta")) #paste0(seqs_folder,"/",subtype,".fasta")
}

seqs_folder_naive <- "data/subtype_seqs_naive"

for(i in 1:length(incl_subtypes)) {
	match_seqs_subtype_md("data/seqs_20230518.fasta", demog_md_subtype_match_naive_uq_dates_list, incl_subtypes[i], seqs_folder=seqs_folder_naive)
}

# Mean sequences per patient
fasta_all <- read.FASTA("data/seqs_20230518.fasta", type = "DNA")
names(fasta_all) <- paste0("t.",names(fasta_all)) # ok, same number of observations md and fasta
subtype_md_nseqs <- subtype_md %>% add_count(patientindex, name="nseqs") %>% group_by(patientindex) %>% filter(row_number() >= (n())) #filter(n() == 1)
mean(subtype_md_nseqs$nseqs)
sd(subtype_md_nseqs$nseqs)
median(subtype_md_nseqs$nseqs)
IQR(subtype_md_nseqs$nseqs); table(subtype_md_nseqs$nseqs); hist(subtype_md_nseqs$nseqs)

incl_subtypes_adj <- c("A_A1","CRF_02_AG","C","B")

### STEP 3: ADJUST SEQUENCES AND BUILD ALIGNMENTS (MAFFT) + MERGE ###
# Make sure seqs are the same length and mask terminal regions
# Manually rename files to match incl_subtypes_adj
for(i in 1:length(incl_subtypes_adj)) {
	system(glue("python3 scripts/padding_seqs.py {seqs_folder_naive}/{incl_subtypes_adj[i]}.fasta {seqs_folder_naive}/{incl_subtypes_adj[i]}_curated.fasta"))
}
# Also replace terminal gaps into missing char (?) manually in Aliview

# Build alignment for subtype B previosly aligned UKRDB samples with reference sequences
# https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html
# Download ref alignment (https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/align.cgi) and then get subtype C oldest sequence + remove gaps

# Subtypes A1, C, CRF02AG (alignment against reference B)
incl_subtypes_adj_notB <- c("A_A1","CRF_02_AG","C")
for(i in 1:length(incl_subtypes_adj_notB)) {
	system(glue("mafft --6merpair --thread 4 --addfragments {seqs_folder_naive}/{incl_subtypes_adj_notB[i]}_curated.fasta data/refs/HIV1_REF_2020_pol_DNA_SUBTYPE_B.fasta > {seqs_folder_naive}/{incl_subtypes_adj_notB[i]}_curated_refB_aln.fasta"))
}

# Subtybe B (alignment against reference C)
system(glue("mafft --6merpair --thread 4 --addfragments {seqs_folder_naive}/B_curated.fasta data/refs/HIV1_REF_2020_pol_DNA_SUBTYPE_C.fasta > {seqs_folder_naive}/B_curated_refC_aln.fasta"))
# IMPORTANT: Manually curate ends of alns {A1,CRF,C,B}_curated_ref{B,C}_aln (trim at length 995) and save file as same name

# IMPORTANT: exclude seqs with 120 sites gap from positions ~250 to ~365
# Filter based on length: keep only >900 nt long
for(i in 1:length(incl_subtypes_adj_notB)) {
	system(glue("python3 scripts/drop_by_length.py {seqs_folder_naive}/{incl_subtypes_adj[i]}_curated_refB_aln.fasta {seqs_folder_naive}/{incl_subtypes_adj[i]}_curated_refB_aln_len_filter.fasta 900"))
}
system(glue("python3 scripts/drop_by_length.py {seqs_folder_naive}/B_curated_refC_aln.fasta {seqs_folder_naive}/B_curated_refC_aln_len_filter.fasta 900"))

# A1 from 3133 to 2852
# CRF_02_AG from 3127 to 2866
# C from 12547 to 11785
# B from 27938 to 25338

# Add ? to 12 sites in the beginning of aln + exclude seqs if >50 ? in beginning
# excluded seqs with so many ??? (>50) in aln end
# manually mask sites with big gaps (>=10) remaining + remove seqs where gaps >50

# A1 from 2852 to 2837
# CRF_02_AG from 2866 to 2845
# C from 11785 to 11711
# B from 25338 to 25282

# Subtype B
# Get random disjoint subsets of sequences (10 sets of 2528 seqs from a total of 25282) to estimate ML trees and timetrees for outlier removal
fasta_b <- read.dna(glue("{seqs_folder_naive}/B_curated_refC_aln_len_filter.fasta"), format="fasta")
n_sets <- 10
n_seqs <- nrow(fasta_b)
sample_size <- round(n_seqs/n_sets)
# starting at 2nd position to avoid including reference
rep_seqs1 <- rep(2:sample_size)
rep_seqs2 <- rep((sample_size+1):n_seqs)
rep_seqs <- c(rep_seqs1, rep_seqs2)
rep_seqs_split <- split(rep_seqs, rep_len(1:n_sets, length(rep_seqs)) )
#set.seed(9548732)
for(i in 1:n_sets){
	aln <- fasta_b[ rep_seqs_split[[i]], ]
	write.FASTA( aln, file = glue("{seqs_folder_naive}/B_sample_{i}.fasta"))
}
# NOTE: add the reference (subtype C) back manually to all subsets

# Subtype C: 5 subsets of 2342 seqs
fasta_c <- read.dna(glue("{seqs_folder_naive}/C_curated_refB_aln_len_filter.fasta"), format="fasta")
n_sets_c <- 5; n_seqs_c <- nrow(fasta_c); sample_size_c <- round(n_seqs_c/n_sets_c)
# starting at 2nd position to avoid including reference
rep_seqs1c <- rep(2:sample_size_c); rep_seqs2c <- rep((sample_size_c+1):n_seqs_c); rep_seqs_c <- c(rep_seqs1c, rep_seqs2c)
rep_seqs_split_c <- split(rep_seqs_c, rep_len(1:n_sets_c, length(rep_seqs_c)) )
#set.seed(9548732)
for(i in 1:n_sets_c){
	aln <- fasta_c[ rep_seqs_split_c[[i]], ]
	write.FASTA( aln, file = glue("{seqs_folder_naive}/C_sample_{i}.fasta"))
}
# NOTE: add the reference (subtype B) back manually to all subsets

fasta_a1 <- read.dna(glue("{seqs_folder_naive}/A_A1_curated_refB_aln_len_filter.fasta"), format="fasta")
fasta_crf02ag <- read.dna(glue("{seqs_folder_naive}/CRF_02_AG_curated_refB_aln_len_filter.fasta"), format="fasta")

### STEP 4: BUILD ML TREE (IQTREE) AND CURATE MD, ALIGNMENT AND TREES
# Load IQTREE PATH from .Renviron file placed in HOME directory

# NOTE: running all on HPC with 5000 bootstrap runs for each tree joined in the end

# Subtype B
iqtree_bin <- Sys.getenv("IQTREE") # this points to iqtree-2.2.0-Linux
for(i in 1:n_sets) {
	system(glue("{iqtree_bin} -s {seqs_folder_naive}/B_sample_{i}.fasta -m GTR+R -nt AUTO -ntmax 2 -B 1000 -nm 5000 -bcor 0.98")) #previously -m GTR+G -ninit 50 -n 50
}

# Subtype C
for(i in 1:n_sets_c) {
	system(glue("{iqtree_bin} -s {seqs_folder_naive}/C_sample_{i}.fasta -m GTR+R -nt AUTO -ntmax 2 -B 1000 -nm 5000 -bcor 0.90")) #previously -m GTR+G -ninit 50 -n 50
}

output_iqtree_folder <- "results/iqtree"
system(glue("mkdir -p {output_iqtree_folder}"))
system(glue("mv {seqs_folder_naive}/C_sample_*.bionj {seqs_folder_naive}/C_sample_*.gz {seqs_folder_naive}/C_sample_*.iqtree {seqs_folder_naive}/C_sample_*.log {seqs_folder_naive}/C_sample_*.mldist {seqs_folder_naive}/C_sample_*.treefile {output_iqtree_folder}"))

system(glue("{iqtree_bin} -s {output_iqtree_folder}/A_A1_curated_refB_aln_len_filter.fasta -m GTR+R -nt AUTO -ntmax 2 -B 1000 -nm 5000 -bcor 0.90"))
system(glue("{iqtree_bin} -s {output_iqtree_folder}/CRF_02_AG_curated_refB_aln_len_filter.fasta -m GTR+R -nt AUTO -ntmax 2 -B 1000 -nm 5000 -bcor 0.95"))

# move resulting files to results/iqtree folder
system(glue("mkdir -p {output_iqtree_folder}"))
system(glue("mv {seqs_folder_naive}/*.bionj {seqs_folder_naive}/*.gz {seqs_folder_naive}/*.iqtree {seqs_folder_naive}/*.log {seqs_folder_naive}/*.mldist {seqs_folder_naive}/*.treefile {output_iqtree_folder}"))
