library(ggpubr)
library(ggplot2)
library(cowplot)
library(glue)

helv="Helvetica"
RDS_PATH="rds"
RESULTS_PATH="results"

# Figure 1: ML tree of consensus phylotypes from partial-pol sequences with CD4 and VL annotations
f3_5 <- readRDS(glue("{RDS_PATH}/f3_5.rds"))
ggsave(file=glue("{RESULTS_PATH}/figs/fig1_to_edit.svg"), dpi=600, width=8, height=10, bg="white")
# improved aes on inkscape and saved as fig1.svg

# Figure 2
f2a_vl <- readRDS(glue("{RDS_PATH}/f2a_vl.rds"))
f2b <- readRDS(glue("{RDS_PATH}/f2b.rds"))
f2_age_bxp <- readRDS(glue("{RDS_PATH}/f2_age_bxp.rds"))
f1c <- readRDS(glue("{RDS_PATH}/f1c.rds"))
row1 <- cowplot::plot_grid(f2a_vl + theme(legend.position = "none"), f2b, ncol = 2, rel_widths = c(2, 3), labels = c("A", "B"))
row2 <- cowplot::plot_grid(f2_age_bxp + theme(legend.position = "none"), f1c, ncol = 2, rel_widths = c(2, 5), labels = c("C", "D"))
f2 <- cowplot::plot_grid(row1, row2, nrow = 2,label_fontfamily = helv,label_fontface = "plain",label_size = 10)
ggsave(plot=f2, file=glue("{RESULTS_PATH}/figs/fig2.eps"), device=cairo_ps, dpi=600, width=8, height=10, bg="white")
ggsave(plot=f2, file=glue("{RESULTS_PATH}/figs/fig2.jpg"), dpi=600, width=8, height=10, bg="white")

# Figure 3: growth (Ne and logistic growth)
f1b_pl_t50_ok <- readRDS(glue("{RDS_PATH}/f1b_pl_t50_ok.rds"))
f1c_pl <- readRDS(glue("{RDS_PATH}/f1c_pl.rds"))
ggarrange(f1b_pl_t50_ok + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), f1c_pl, nrow=2, ncol=1, labels = c("A","B"), font.label=list(family="Helvetica", color="black",size=10)) #legend="top"
ggsave(file=glue("{RESULTS_PATH}/figs/fig3.eps"), device = cairo_ps, dpi=600, width=8, height=9, bg="white")
ggsave(file=glue("{RESULTS_PATH}/figs/fig3.jpg"), dpi=600, width=8, height=9, bg="white")

# Figure 4: ML tree with global and timetre with annots
f4ac <- readRDS(glue("{RDS_PATH}/f4ac.rds"))
f4bd_pl <- readRDS(glue("{RDS_PATH}/f4bd_pl.rds"))
ggarrange(f4ac, f4bd_pl, nrow=1, ncol=2, font.label=list(family=helv, color="black",size=10), heights=c(1,1))
ggsave(file=glue("{RESULTS_PATH}/figs/fig4_to_edit.svg"), dpi=600, width=10, height=12, bg="white")
# improved aes on inkscape and saved as fig4.svg