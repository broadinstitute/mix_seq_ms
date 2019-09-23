library(magrittr)
library(tidyverse)
library(here)

source(here::here('src', 'MixSeq_helpers.R'))
source(here::here('src', 'expt_meta_data.R'))
results_dir <- here::here('data')

run_sum_collapse <- function(cur_expt) {
  sprintf('Processing expt %s', cur_expt$expt_name)
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  seuObj <- load_sc_data(cur_expt)
  summed_counts <- get_summed_counts(seuObj)
  readr::write_rds(summed_counts, file.path(out_dir, 'summed_counts.rds'))
  
  #compute per condition average expression
  seuObj <- Seurat::NormalizeData(object = seuObj, 
                                  normalization.method = "LogNormalize", 
                                  scale.factor = 1e6) #CPM normalization
  
  cond_avgs <- Seurat::AverageExpression(seuObj, use.scale = FALSE, return.seurat = FALSE)$RNA %>% t() 
  readr::write_rds(cond_avgs, file.path(out_dir, 'avg_profiles_cpm.rds'))
}
for (cur_expt in sc_expts) {
  run_sum_collapse(cur_expt)
}
