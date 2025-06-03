# HIV_UKRDB_phylotypes
Code for the paper "Evidence for Circulation of High-Virulence HIV-1 Subtype B 
Variants in the United Kingdom", which uses genomic and clinical (viral loads and CD4 counts) 
data from the UK Drug Resistance Database (UKRDB or UKDRD) to investigate HIV-1 variants with increased virulence.

This repository includes ~30 modular R scripts used in the analysis. The prefix of each file denotes its order.
The general workflow is as follows:

## 1. Identification of phylotypes with significant VL and/or CD4 differences (VOIs)

### 1a. Data preparation, alignments, and maximum likelihood (ML) trees
- `01_md_ukdrd_splits_subtypes.R`: Preparation of UK HIV Drug Resistance Database (UKDRD) data, 
splitting by subtypes and performing alignments and maximum likelihood trees
- `drop_by_length.py` and `padding_seqs.py` are helper scripts

### 1b. Time-scaled trees
- `02_subtype_a1.R`: Time-scaled tree for subtype A1
- `02_subtype_b.R`: Time-scaled tree for subtype B
- `02_subtype_c.R`: Time-scaled tree for subtype C
- `02_subtype_crf_02ag.R`: Time-scaled tree for CRF_02AG subtype

### 1c. Phylotype (i.e. cluster or variant) assignment
- `03_phylotypes.R`: Identification and classification of phylotypes
- `04_extr_clusters_flag_paraphyletic.R`: Extraction of clusters and flagging of paraphyletic groups
- `04_generate_aln_all_phylotype_seqs.R`: Generation of alignments for all phylotype sequences

### 1d. Viral load models/tests
- `stan/vl_model.stan`: Model for viral load differences implemented in Stan
- `05_crossref_test_vl.R`: Cross-referencing and models/tests for viral load differences

### 1e. CD4 decline models
- `06_crossref_cd4.R`: Cross-referencing CD4 cell count data
- `06_regression_cd4-0.1.1_ml_randeff.R`: CD4 decline model using maximum likelihood with random effects
- `06_regression_cd4-0.1.2_bayes_randeff.R`: CD4 decline model using Bayesian approach with random effects
- `06_regression_cd4-0.1.3_r2d2_prior_bayes.R`: CD4 decline model using Bayesian approach with R2D2 prior
- `06_regression_cd4-0.2_bayes_suspVOIs_vs_backbone.R`: Bayesian regression comparing suspected VOIs vs backbone phhylotype
- `07_cd4_episode_distr.R`: Distribution of CD4 measurements and related general-purpose code for paper

### 1f. Heritability of viral load and internode intervals
- `08_internodeIntervals_POUMM.R`: Analysis of internode intervals and viral load heritability using Phylogenetic Ornstein-Uhlenbeck Mixed Model (POUMM)

## 2. Targeted analyses on Variants of Interest (VOIs)

### 2a. Effective population size and logistic growth
- `09_ne.R`: Effective population size estimation
- `10_logisticgrowth-0.1.R`: Logistic growth rate estimation

### 2b. Demographic comparisons
- `11_demog_comparisons.R`: Comparison of demographics (ethnicity, region, risk group, sex)

### 2c. Drug resistance mutations (DRMs)
- `12_drm_sierra.sh`: Analysis of drug resistance mutations using sierra-local

### 2d. CD4 model and exploratory analyses for post-ART data
- `13_post_art_outcomes.R`: Analysis of potential outcomes after antiretroviral therapy (ART)

### 2e. Subtype B ML whole and subsampled tree
- `14_whole_tree_brlen_test.R`: Getting complete (e.g. >20k seqs for subtype B) trees and testing branch lengths for subtype B

### 2f. Estimate partial-pol phylotype consensus sequences and plot tree/distance matrix
- `15_pt_consensus.R`: Estimation of partial-pol phylotype consensus sequences
- `16_dist_tree_consensus.R`: Plotting tree and distance matrix for consensus sequences

### 2g. BLAST VOI sequences and date trees
- `17_blastdb_lanl_tree_global.R`: BLAST VOI sequences against Los Alamos National Laboratory (LANL) database matches at >= 95% identity
- `18_dating-0.1_ml.R`: ML-based tree dating (treedater)
- `18_dating-0.2_beast_prep.R`: Preparation for BEAST tree dating (BEAST1)

### 2h. Plot VOI+LANL ML trees and timetrees with annotations
- `19_treevis.R`: Visualization of trees with sociodemographic and clinical annotations

### 2i. UKHSA whole genomes and mutation extraction
- `20_ukhsa_trees.R`: Find UKHSA whole genomes closely related to UKRDB sequences at partial-pol level
- `20_muts_wgs_ref_hxb2.R`: Extract mutations relative to HXB2 and UKRDB consensus

### 2j. Coreceptor usage
- `21_coreceptor_usage.R`: Analysis of HIV coreceptor usage

## Utilities
- `00_join_plots.R`: Utility script for joining main paper plots

## UKRDB data access

Individual-level HIV genetic, demographic, and clinical data from the UKRDB can be accessed for collaborative projects, 
subject to the approval of a research proposal by the UKRDB Steering Committee.

## What was possible to make public

Phylotype-level consensus sequences have been deposited in Zenodo (https://doi.org/10.5281/zenodo.14792962).

## Getting Started

Packages needed are at the top of each script and need to be installed if e.g. not part of base R. The version of R used during paper preparation was: R 4.1.3.

See the session info below for specific package versions used:

```r
─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.1.3 (2022-03-10)
 os       Ubuntu 20.04.6 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language en_US:en
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Isle_of_Man
 date     2025-05-30
 rstudio  2024.12.0+467 Kousa Dogwood (desktop)
 pandoc   2.5 @ /usr/bin/pandoc
 quarto   1.5.57 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto

─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 ! package           * version    date (UTC) lib source
   abind               1.4-5      2016-07-21 [1] CRAN (R 4.1.3)
   adaptMCMC           1.5        2024-01-29 [1] CRAN (R 4.1.3)
   ade4                1.7-18     2021-09-16 [1] CRAN (R 4.1.3)
   ape               * 5.8-0.4    2024-10-15 [1] Github (emmanuelparadis/ape@5c14a70)
   aplot               0.1.3      2022-04-01 [1] CRAN (R 4.1.3)
   assertthat          0.2.1      2019-03-21 [1] CRAN (R 4.1.3)
   backports           1.4.1      2021-12-13 [1] CRAN (R 4.1.3)
   base64enc           0.1-3      2015-07-28 [1] CRAN (R 4.1.3)
   bayesplot         * 1.10.0     2022-11-16 [1] CRAN (R 4.1.3)
   bayestestR          0.13.0     2022-09-18 [1] CRAN (R 4.1.3)
   betareg           * 3.2-1      2024-09-12 [1] CRAN (R 4.1.3)
   BiocGenerics      * 0.40.0     2021-10-26 [1] Bioconductor
   Biostrings        * 2.62.0     2021-10-26 [1] Bioconductor
   bitops              1.0-7      2021-04-24 [1] CRAN (R 4.1.3)
   boot              * 1.3-31     2024-08-28 [1] CRAN (R 4.1.3)
   bridgesampling      1.1-2      2021-04-16 [1] CRAN (R 4.1.3)
   brms              * 2.19.0     2023-03-14 [1] CRAN (R 4.1.3)
   Brobdingnag         1.2-9      2022-10-19 [1] CRAN (R 4.1.3)
   broom               0.7.12     2022-01-28 [1] CRAN (R 4.1.3)
   callr               3.7.0      2021-04-20 [1] CRAN (R 4.1.3)
   caper             * 1.0.1      2018-04-17 [1] CRAN (R 4.1.3)
   car               * 3.0-12     2021-11-06 [1] CRAN (R 4.1.3)
   carData           * 3.0-5      2022-01-06 [1] CRAN (R 4.1.3)
   cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.1.3)
   checkmate           2.1.0      2022-04-21 [1] CRAN (R 4.1.3)
   cli                 3.6.1      2023-03-23 [1] CRAN (R 4.1.3)
   clusterGeneration   1.3.7      2020-12-15 [1] CRAN (R 4.1.3)
   coda                0.19-4     2020-09-30 [1] CRAN (R 4.1.3)
   codetools           0.2-20     2024-03-31 [1] CRAN (R 4.1.3)
   colorspace          2.1-0      2023-01-23 [1] CRAN (R 4.1.3)
   colourpicker        1.2.0      2022-10-28 [1] CRAN (R 4.1.3)
   combinat            0.0-8      2012-10-29 [1] CRAN (R 4.1.3)
   cowplot           * 1.1.1      2020-12-30 [1] CRAN (R 4.1.3)
   crayon              1.5.2      2022-09-29 [1] CRAN (R 4.1.3)
   crosstalk           1.2.0      2021-11-04 [1] CRAN (R 4.1.3)
   curl                4.3.2      2021-06-23 [1] CRAN (R 4.1.3)
   data.table        * 1.16.2     2024-10-10 [1] CRAN (R 4.1.3)
   datawizard          0.6.3      2022-10-22 [1] CRAN (R 4.1.3)
   DBI                 1.1.3      2022-06-18 [1] CRAN (R 4.1.3)
   dbplyr              2.1.1      2021-04-06 [1] CRAN (R 4.1.3)
   digest              0.6.37     2024-08-19 [1] CRAN (R 4.1.3)
   distributional      0.3.2      2023-03-22 [1] CRAN (R 4.1.3)
   dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.1.3)
   DT                  0.21       2022-02-26 [1] CRAN (R 4.1.3)
   dygraphs            1.1.1.6    2018-07-11 [1] CRAN (R 4.1.3)
   effectsize          0.8.1      2022-10-18 [1] CRAN (R 4.1.3)
   ellipsis            0.3.2      2021-04-29 [1] CRAN (R 4.1.3)
   emmeans             1.8.4-1    2023-01-17 [1] CRAN (R 4.1.3)
   estimability        1.4.1      2022-08-05 [1] CRAN (R 4.1.3)
   expm                0.999-6    2021-01-13 [1] CRAN (R 4.1.3)
   fansi               1.0.4      2023-01-22 [1] CRAN (R 4.1.3)
   farver              2.1.1      2022-07-06 [1] CRAN (R 4.1.3)
   fastbaps          * 1.0.8      2024-11-18 [1] Github (gtonkinhill/fastbaps@1c32218)
   fastmap             1.1.1      2023-02-24 [1] CRAN (R 4.1.3)
   fastmatch           1.1-3      2021-07-23 [1] CRAN (R 4.1.3)
   flexmix             2.3-18     2022-06-07 [1] CRAN (R 4.1.3)
   forcats           * 0.5.1      2021-01-27 [1] CRAN (R 4.1.3)
   foreach             1.5.2      2022-02-02 [1] CRAN (R 4.1.3)
   Formula             1.2-4      2020-10-16 [1] CRAN (R 4.1.3)
   fs                  1.6.1      2023-02-06 [1] CRAN (R 4.1.3)
   generics            0.1.3      2022-07-05 [1] CRAN (R 4.1.3)
   GenomeInfoDb      * 1.30.1     2022-01-30 [1] Bioconductor
   GenomeInfoDbData    1.2.7      2022-04-26 [1] Bioconductor
   ggeffects           1.1.4      2022-10-23 [1] CRAN (R 4.1.3)
   ggforce           * 0.3.3      2021-03-05 [1] CRAN (R 4.1.3)
   ggfun               0.0.6      2022-04-01 [1] CRAN (R 4.1.3)
   ggnewscale        * 0.4.7      2022-03-25 [1] CRAN (R 4.1.3)
   ggplot2           * 3.5.1      2024-04-23 [1] CRAN (R 4.1.3)
   ggplotify           0.1.0      2021-09-02 [1] CRAN (R 4.1.3)
   ggpubr            * 0.6.0      2023-02-10 [1] CRAN (R 4.1.3)
   ggsci             * 3.0.0      2023-03-08 [1] CRAN (R 4.1.3)
   ggsignif            0.6.3      2021-09-09 [1] CRAN (R 4.1.3)
   ggtree            * 3.2.1      2021-11-16 [1] Bioconductor
   glue              * 1.8.0      2024-09-30 [1] CRAN (R 4.1.3)
   gridExtra           2.3        2017-09-09 [1] CRAN (R 4.1.3)
   gridGraphics        0.5-1      2020-12-13 [1] CRAN (R 4.1.3)
   gtable              0.3.1      2022-09-01 [1] CRAN (R 4.1.3)
   gtools            * 3.9.5      2023-11-20 [1] CRAN (R 4.1.3)
   haven               2.4.3      2021-08-04 [1] CRAN (R 4.1.3)
   hms                 1.1.2      2022-08-19 [1] CRAN (R 4.1.3)
   htmltools           0.5.4      2022-12-07 [1] CRAN (R 4.1.3)
   htmlwidgets         1.6.2      2023-03-17 [1] CRAN (R 4.1.3)
   httpuv              1.6.9      2023-02-14 [1] CRAN (R 4.1.3)
   httr                1.4.4      2022-08-17 [1] CRAN (R 4.1.3)
   igraph              1.4.1      2023-02-24 [1] CRAN (R 4.1.3)
   inline              0.3.19     2021-05-31 [1] CRAN (R 4.1.3)
   insight             0.18.6     2022-10-23 [1] CRAN (R 4.1.3)
   IRanges           * 2.28.0     2021-10-26 [1] Bioconductor
   iterators           1.0.14     2022-02-05 [1] CRAN (R 4.1.3)
   jsonlite            1.8.4      2022-12-06 [1] CRAN (R 4.1.3)
   knitr               1.42       2023-01-25 [1] CRAN (R 4.1.3)
   lamW                2.2.4      2024-06-17 [1] CRAN (R 4.1.3)
   later               1.3.0      2021-08-18 [1] CRAN (R 4.1.3)
   lattice             0.22-6     2024-03-20 [1] CRAN (R 4.1.3)
   lazyeval            0.2.2      2019-03-15 [1] CRAN (R 4.1.3)
   lifecycle           1.0.3      2022-10-07 [1] CRAN (R 4.1.3)
   limSolve          * 1.5.6      2019-11-14 [1] CRAN (R 4.1.3)
   lme4              * 1.1-30     2022-07-08 [1] CRAN (R 4.1.3)
   lmtest              0.9-39     2021-11-07 [1] CRAN (R 4.1.3)
   loo                 2.6.0      2023-03-31 [1] CRAN (R 4.1.3)
   lpSolve             5.6.18     2023-02-01 [1] CRAN (R 4.1.3)
   lubridate         * 1.9.2      2023-02-10 [1] CRAN (R 4.1.3)
   magrittr          * 2.0.3      2022-03-30 [1] CRAN (R 4.1.3)
   maps              * 3.4.0      2021-09-25 [1] CRAN (R 4.1.3)
   markdown            1.13       2024-06-04 [1] CRAN (R 4.1.3)
   MASS              * 7.3-60     2023-05-04 [1] CRAN (R 4.1.3)
   Matrix            * 1.3-4      2021-06-01 [1] CRAN (R 4.1.3)
   matrixStats         0.61.0     2021-09-17 [1] CRAN (R 4.1.3)
   mgcv              * 1.9-1      2023-12-21 [1] CRAN (R 4.1.3)
   mime                0.12       2021-09-28 [1] CRAN (R 4.1.3)
   miniUI              0.1.1.1    2018-05-18 [1] CRAN (R 4.1.3)
   minqa               1.2.4      2014-10-09 [1] CRAN (R 4.1.3)
   mlesky            * 0.1.7      2025-01-31 [1] Github (emvolz-phylodynamics/mlesky@28d97b8)
   mnormt              2.0.2      2020-09-01 [1] CRAN (R 4.1.3)
   modelr              0.1.8      2020-05-19 [1] CRAN (R 4.1.3)
   modeltools          0.2-23     2020-03-05 [1] CRAN (R 4.1.3)
   munsell             0.5.0      2018-06-12 [1] CRAN (R 4.1.3)
   mvtnorm           * 1.1-3      2021-10-08 [1] CRAN (R 4.1.3)
   nlme              * 3.1-166    2024-08-14 [1] CRAN (R 4.1.3)
   nloptr              2.0.0      2022-01-26 [1] CRAN (R 4.1.3)
   nnet                7.3-19     2023-05-03 [1] CRAN (R 4.1.3)
   numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 4.1.3)
   openxlsx          * 4.2.5.2    2023-02-06 [1] CRAN (R 4.1.3)
   parameters          0.19.0     2022-10-05 [1] CRAN (R 4.1.3)
   patchwork         * 1.1.1      2020-12-17 [1] CRAN (R 4.1.3)
   pbmcapply         * 1.5.1      2022-04-28 [1] CRAN (R 4.1.3)
   performance       * 0.10.0     2022-10-03 [1] CRAN (R 4.1.3)
   phangorn            2.11.1     2023-01-23 [1] CRAN (R 4.1.3)
   phytools          * 1.0-1      2022-01-03 [1] CRAN (R 4.1.3)
   pillar              1.9.0      2023-03-22 [1] CRAN (R 4.1.3)
   pkgbuild            1.3.1      2021-12-20 [1] CRAN (R 4.1.3)
   pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.1.3)
   plotrix             3.8-2      2021-09-08 [1] CRAN (R 4.1.3)
   plyr                1.8.7      2022-03-24 [1] CRAN (R 4.1.3)
   polyclip            1.10-0     2019-03-14 [1] CRAN (R 4.1.3)
   posterior         * 1.4.1      2023-03-14 [1] CRAN (R 4.1.3)
   POUMM             * 2.1.7      2020-10-27 [1] CRAN (R 4.1.3)
   prettyunits         1.1.1      2020-01-24 [1] CRAN (R 4.1.3)
   processx            3.5.2      2021-04-30 [1] CRAN (R 4.1.3)
   promises            1.2.0.1    2021-02-11 [1] CRAN (R 4.1.3)
   ps                  1.6.0      2021-02-28 [1] CRAN (R 4.1.3)
   purrr             * 1.0.1      2023-01-10 [1] CRAN (R 4.1.3)
   quadprog            1.5-8      2019-11-20 [1] CRAN (R 4.1.3)
   QuickJSR            1.2.2      2024-06-07 [1] CRAN (R 4.1.3)
   R6                  2.5.1      2021-08-19 [1] CRAN (R 4.1.3)
   rBLAST            * 0.99.2     2023-06-26 [1] https://mhahsler.r-universe.dev (R 4.1.3)
   RColorBrewer      * 1.1-3      2022-04-03 [1] CRAN (R 4.1.3)
   Rcpp              * 1.0.14     2025-01-12 [1] CRAN (R 4.1.3)
   RcppParallel        5.1.5      2022-01-05 [1] CRAN (R 4.1.3)
   RCurl               1.98-1.6   2022-02-08 [1] CRAN (R 4.1.3)
   readr             * 2.1.4      2023-02-10 [1] CRAN (R 4.1.3)
   readxl              1.4.0      2022-03-28 [1] CRAN (R 4.1.3)
   reprex              2.0.1      2021-08-05 [1] CRAN (R 4.1.3)
   reshape2          * 1.4.4      2020-04-09 [1] CRAN (R 4.1.3)
   rJava               1.0-6      2021-12-10 [1] CRAN (R 4.1.3)
   rlang               1.1.4      2024-06-04 [1] CRAN (R 4.1.3)
   rlist             * 0.4.6.2    2021-09-03 [1] CRAN (R 4.1.3)
   rstan             * 2.32.6     2024-03-05 [1] CRAN (R 4.1.3)
   rstantools          2.3.1      2023-03-30 [1] CRAN (R 4.1.3)
   rstatix             0.7.2      2023-02-01 [1] CRAN (R 4.1.3)
   rstudioapi          0.13       2020-11-12 [1] CRAN (R 4.1.3)
   rsvg              * 2.6.2      2025-03-23 [1] CRAN (R 4.1.3)
   rvest               1.0.2      2021-10-16 [1] CRAN (R 4.1.3)
   S4Vectors         * 0.32.4     2022-03-24 [1] Bioconductor
   sandwich            3.0-1      2021-05-18 [1] CRAN (R 4.1.3)
   scales            * 1.3.0      2023-11-28 [1] CRAN (R 4.1.3)
   scatterplot3d       0.3-41     2018-03-14 [1] CRAN (R 4.1.3)
   seqinr            * 4.2-8      2021-06-09 [1] CRAN (R 4.1.3)
   sessioninfo       * 1.2.3      2025-02-05 [1] CRAN (R 4.1.3)
   shiny               1.7.4      2022-12-15 [1] CRAN (R 4.1.3)
   shinyjs             2.1.0      2021-12-23 [1] CRAN (R 4.1.3)
   shinystan           2.6.0      2022-03-03 [1] CRAN (R 4.1.3)
   shinythemes         1.2.0      2021-01-25 [1] CRAN (R 4.1.3)
   sjlabelled          1.2.0      2022-04-10 [1] CRAN (R 4.1.3)
   sjmisc              2.8.9      2021-12-03 [1] CRAN (R 4.1.3)
   sjPlot            * 2.8.11     2022-08-07 [1] CRAN (R 4.1.3)
   sjstats             0.18.1     2021-01-09 [1] CRAN (R 4.1.3)
   splitstackshape   * 1.4.8      2019-04-21 [1] CRAN (R 4.1.3)
   StanHeaders       * 2.32.9     2024-05-29 [1] CRAN (R 4.1.3)
   stringi             1.7.12     2023-01-11 [1] CRAN (R 4.1.3)
   stringr           * 1.5.0      2022-12-02 [1] CRAN (R 4.1.3)
   tensorA             0.36.2     2020-11-19 [1] CRAN (R 4.1.3)
   threejs             0.3.3      2020-01-21 [1] CRAN (R 4.1.3)
   tibble            * 3.2.1      2023-03-20 [1] CRAN (R 4.1.3)
   tidyr             * 1.3.0      2023-01-24 [1] CRAN (R 4.1.3)
   tidyselect          1.2.0      2022-10-10 [1] CRAN (R 4.1.3)
   tidytree            0.3.9      2022-03-04 [1] CRAN (R 4.1.3)
   tidyverse         * 1.3.1      2021-04-15 [1] CRAN (R 4.1.3)
   timechange          0.2.0      2023-01-11 [1] CRAN (R 4.1.3)
   tmvnsim             1.0-2      2016-12-15 [1] CRAN (R 4.1.3)
   treedater         * 0.5.3      2022-03-19 [1] Github (emvolz/treedater@7b8a72a)
   treeio            * 1.18.1     2021-11-14 [1] Bioconductor
   treestructure     * 0.4.0      2024-11-09 [1] Github (emvolz-phylodynamics/treestructure@d215cd1)
   tweenr              1.0.2      2021-03-23 [1] CRAN (R 4.1.3)
   tzdb                0.3.0      2022-03-28 [1] CRAN (R 4.1.3)
   utf8                1.2.3      2023-01-31 [1] CRAN (R 4.1.3)
   V8                  4.4.2      2024-02-15 [1] CRAN (R 4.1.3)
   vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.1.3)
   viridis           * 0.6.4      2023-07-22 [1] CRAN (R 4.1.3)
   viridisLite       * 0.4.1      2022-08-22 [1] CRAN (R 4.1.3)
   withr               2.5.0      2022-03-03 [1] CRAN (R 4.1.3)
   xlsx              * 0.6.5      2020-11-10 [1] CRAN (R 4.1.3)
   xlsxjars            0.6.1      2014-08-22 [1] CRAN (R 4.1.3)
   xml2                1.3.5      2023-07-06 [1] CRAN (R 4.1.3)
   xtable              1.8-4      2019-04-21 [1] CRAN (R 4.1.3)
   xts                 0.12.1     2020-09-09 [1] CRAN (R 4.1.3)
   XVector           * 0.34.0     2021-10-26 [1] Bioconductor
   yulab.utils         0.0.4      2021-10-09 [1] CRAN (R 4.1.3)
   zip                 2.3.0      2023-04-17 [1] CRAN (R 4.1.3)
   zlibbioc            1.40.0     2021-10-26 [1] Bioconductor
   zoo                 1.8-9      2021-03-09 [1] CRAN (R 4.1.3)

 * ── Packages attached to the search path.
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

## Citation

Provisionally, you can cite this preprint:

Franceschi, Vinicius B. and Drake, Kieran O. and Bibby, David F. and Sabin, Caroline and Dunn, David T. and Mbisa, Jean L. and Volz, Erik. 
Evidence for Circulation of High-Virulence HIV-1 Subtype B Variants in the United Kingdom. 
Available at SSRN: https://ssrn.com/abstract=4929798

The analyses presented previously in the preprint were improved here, including, in particular, 
comparisons against other partitioning methods, enhanced viral load and CD4 decline models, 
and a new selection of VOIs. Despite these enhancements, the results remain consistent
with those reported in the preprint. The more robust and updated analyses reported here 
have been submitted to a journal.