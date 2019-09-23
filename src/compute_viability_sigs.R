library(magrittr)
library(tidyverse)
library(taigr)
library(here)
source(here::here('src', 'MixSeq_helpers.R'))
source(here::here('src', 'expt_meta_data.R'))
source(here::here('src', 'global_params.R'))
results_dir <- here::here('data')

all_CL_features <- read_rds(file.path(results_dir, 'all_CL_features.rds'))

#PARAMS
min_cells <- 5 #minimum cells per condition to include cell line in analysis

compute_viability_sig <- function(cur_expt) {
  print(sprintf('Processing expt %s', cur_expt$expt_name))
  out_dir <- file.path(results_dir, cur_expt$expt_name)

  CL_df <- all_CL_features[[cur_expt$expt_name]]
  
  dat <- load_collapsed_profiles(cur_expt, results_dir, type = 'sum')
  # dat <- load_collapsed_profiles(cur_expt, results_dir, type = 'avg_cpm') 
  
  #get usable cell lines
  cell_df <- load_cell_info(cur_expt) %>% 
    tidyr::separate(col = 'condition', into = c('condition', 'batch'), sep = '_')
  usable_cell_lines <- cell_df %>% 
    dplyr::group_by(singlet_ID, condition) %>% 
    dplyr::summarise(n = n()) %>% 
    reshape2::dcast(singlet_ID ~ condition, value.var = 'n') %>% 
    dplyr::filter(control >= min_cells, treat >= min_cells) %>% 
    .[['singlet_ID']]
  print(sprintf('Using %d/%d cell lines with sufficient cells', length(usable_cell_lines), length(unique(cell_df$singlet_ID))))
  
  #restrict to usable cell lines
  usamps <- dat$sample_info %>% dplyr::filter(CCLE_ID %in% usable_cell_lines) %>% .[['sample_name']]
  dat$profile_mat <- dat$profile_mat[, usamps]
  dat$sample_info %<>% dplyr::filter(sample_name %in% usamps)
  
  limma_res <- fit_viability_models(dat$profile_mat, dat$sample_info, CL_df, prior_cnt = globals$prior_cnt, use_voom = globals$use_voom)
  write_rds(limma_res, file.path(out_dir, 'limma_res.rds'))
}
for (ii in sc_DE_meta) {
  compute_viability_sig(ii)
}
