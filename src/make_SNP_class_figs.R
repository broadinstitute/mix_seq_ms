library(tidyverse)
library(here)
library(ggrepel)
library(magrittr)

source(here::here('src', 'MixSeq_helpers.R'))
source(here::here('src', 'expt_meta_data.R'))
results_dir <- here::here('data')
fig_dir <- here::here('figures')

source(here::here('src', 'make_singlet_SNP_heatmap.R'))
make_singlet_SNP_heatmap()

source(here::here('src', 'make_SNP_GE_comparison_fig.R'))
res <- make_SNP_GE_comparison_fig('expt1')
print(res$perc_correct)
res <- make_SNP_GE_comparison_fig('expt3')
print(res$perc_correct)

source(here::here('src', 'make_oop_classification_fig.R'))
res <- make_oop_classification_fig()
print(res)

source(here::here('src', 'make_cell_dist_fig.R'))
make_cell_dist_fig()

source(here::here('src', 'make_demux_SNP_comparison_fig.R'))
res <- make_demux_SNP_comparison_fig()
print(res)

source(here::here('src', 'make_cell_quality_example_fig.R'))
make_cell_quality_example_fig(sc_expts$DMSO_6hr_expt1)

source(here::here('src', 'make_doublet_in_pool_fig.R'))
make_doublet_in_pool_fig()
