library(tidyverse)
library(ggrepel)
library(cowplot)
library(Seurat)
library(here)

source(here::here('src', 'expt_meta_data.R'))
source(here::here('src', 'global_params.R'))
source(here::here('src', 'MixSeq_helpers.R'))

results_dir <- here::here('data')
fig_dir <- here::here('figures')

all_CL_features <- read_rds(file.path(here::here('data', 'all_CL_features.rds')))

gsc_data <- read_rds('~/CPDS/data/MSigDB/gsc_data.rds')
gsc_data$combined <- GSEABase::GeneSetCollection(c(gsc_data$hallmark, gsc_data$canonical))


source(here::here('src', 'make_tsne_example_fig.R'))
make_tsne_example_fig(sc_DE_meta$idasanutlin_24hr_expt1)
make_tsne_example_fig(sc_DE_meta$trametinib_24hr_expt3)

source(here::here('src', 'make_nutlin_G1_arrest_fig.R'))
make_nutlin_G1_arrest_fig()

source(here::here('src', 'make_nutlin_volcanos_heatmap.R'))
make_nutlin_volcanos_heatmap()

source(here::here('src', 'make_example_MOA_plots.R'))
make_example_MOA_plots()

source(here::here('src', 'make_GPX4_figs.R'))
make_GPX4_figs()

source(here::here('src', 'make_hash_tag_qc_figs.R'))
make_hash_tag_qc_figs()

source(here::here('src', 'make_trametinib_tc_figs.R'))
make_trametinib_tc_figs()


source(here::here('src', 'make_viability_sig_plots.R'))
make_viability_sig_plots(sc_DE_meta$trametinib_24hr_expt3)

source(here::here('src', 'make_cell_rep_fig.R'))
make_cell_rep_fig(sc_DE_meta$trametinib_24hr_expt3)

source(here::here('src', 'make_PCA_example_figures.R'))
make_PCA_example_figures(sc_DE_meta$trametinib_24hr_expt3)
make_PCA_example_figures(sc_DE_meta$Trametinib_expt10)

source(here::here('src', 'make_subclone_response_figs.R'))
make_subclone_response_figs(sc_DE_meta$trametinib_24hr_expt1, targ_CL = 'IALM_LUNG')

source(here::here('src', 'make_bortezomib_subpop_figs.R'))
make_bortezomib_subpop_figs(c('RCC10RGB_KIDNEY', 'SNU1079_BILIARY_TRACT', 'RCM1_LARGE_INTESTINE', 'RERFLCAD1_LUNG', 'TEN_ENDOMETRIUM'))

source(here::here('src', 'make_dabrafenib_heterogeneity_figs.R'))
make_dabrafenib_heterogeneity_figs()

source(here::here('src', 'run_CL_downsampling_analysis.R'))
run_CL_downsampling_analysis(sc_DE_meta$trametinib_24hr_expt3)

source(here::here('src', 'run_cell_downsampling_analysis.R'))
run_cell_downsampling_analysis(sc_DE_meta$trametinib_24hr_expt3)


#GLOBAL ANALYSES
source(here::here('src', 'make_CCchange_AUC_plot.R'))
make_CCchange_AUC_plot()

source(here::here('src', 'make_rel_abundance_plot.R'))
make_rel_abundance_plot()

source(here::here('src', 'make_PC_AUC_compare_fig.R'))
make_PC_AUC_compare_fig()

source(here::here('src', 'make_global_dimred.R'))
make_global_dimred()

source(here::here('src', 'predictive_model_compare.R'))
run_predictive_model_compare()
make_predictive_model_figs()

source(here::here('src', 'make_viability_signature_heatmaps.R'))
make_viability_signature_heatmaps()
