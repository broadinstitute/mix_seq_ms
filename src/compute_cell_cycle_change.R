library(magrittr)
library(tidyverse)
library(Seurat)
library(here)

source(here::here('src', 'MixSeq_helpers.R'))
source(here::here('src', 'expt_meta_data.R'))
source(here::here('src', 'global_params.R'))

results_dir <- here::here('data')

#compute estimates of the change in fraction of cells in each cell cycle phase, as well as change in median cell cycle phase scores (comparing treat vs control groups)
compute_cell_cycle_change <- function(cur_expt) {
  sprintf('Processing expt %s', cur_expt$expt_name)
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  seuObj <- load_sc_data(cur_expt, sc_expts)

  seuObj <- NormalizeData(object = seuObj, 
                          normalization.method = "LogNormalize", 
                          scale.factor = globals$scale_fac)
  
  seuObj <- CellCycleScoring(object = seuObj,
                             s.features = Seurat::cc.genes$s.genes,
                             g2m.features = Seurat::cc.genes$g2m.genes,
                             set.ident = FALSE)
  seuObj <- relabel_cell_cycle_phase(seuObj)
  
  cc_df <- seuObj@meta.data %>%
    as.data.frame() %>%
    tidyr::separate(col = 'condition', into = c('condition', 'sample'), sep = '_') %>%
    dplyr::mutate(condition = factor(condition, levels = c('control', 'treat')))
  cc_phases <- unique(cc_df$Phase)
  
  cc_changes <- cc_df %>% 
    ddply(.(singlet_ID), function(cur_df) {
      ldply(cc_phases, function(cur_phase) {
        cur_tab <- table(cur_df$Phase == cur_phase, cur_df$condition)
        if (nrow(cur_tab) == 1 | any(colSums(cur_tab) == 0)) {
          return(NULL)
        }
        prop.test(cur_tab[2,], colSums(cur_tab),
                  conf.level = 0.95) %>% 
          broom::tidy() %>% 
          dplyr::mutate(delta_frac = estimate2 - estimate1,
                 delta_low = -conf.low, 
                 delta_high = -conf.high) %>% 
          dplyr::select(delta_frac, control_prop = estimate1, treat_prop = estimate2, p.value, delta_low, delta_high) %>% 
          mutate(Phase = cur_phase)
      })
    }) %>% 
    dplyr::rename(CCLE_ID = singlet_ID)
  
   median_scores <- full_join(  
     cc_df %>% dplyr::filter(condition == 'control') %>% 
      dplyr::group_by(singlet_ID) %>% 
      dplyr::summarise(control_S = median(S.Score),
                       control_G2M = median(G2M.Score)),
     cc_df %>% dplyr::filter(condition == 'treat') %>% 
      dplyr::group_by(singlet_ID) %>% 
      dplyr::summarise(treat_S = median(S.Score),
                       treat_G2M = median(G2M.Score)),
    by = 'singlet_ID'
  ) %>% mutate(delta_S = treat_S - control_S,
               delta_G2M = treat_G2M - control_G2M) %>% 
    dplyr::rename(CCLE_ID = singlet_ID)
    
  module_scores_df <- full_join(median_scores, cc_changes, by = 'CCLE_ID')
  write_csv(module_scores_df, file.path(out_dir, 'module_scores.csv'))
}
for (ee in sc_DE_meta) {
  print(sprintf('Processing experiment %s', ee$expt_name))
  compute_cell_cycle_change(ee)
}
