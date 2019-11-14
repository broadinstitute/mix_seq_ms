library(magrittr)
library(tidyverse)
library(Seurat)
library(here)

source(here::here('src', 'expt_meta_data.R'))
source(here::here('src', 'global_params.R'))
source(here::here('src', 'MixSeq_helpers.R'))

results_dir <- here::here('data')

#PARAMS
pca_prior_counts <- 10 
min_cells <- 5 #min cells per condition to include cell line in PCA calc

compute_lfc_pca <- function(cur_expt) {
  sprintf('Processing expt %s', cur_expt$expt_name)
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  cell_df <- load_cell_info(cur_expt) %>% 
    tidyr::separate(col = 'condition', into = c('condition', 'batch'), sep = '_')
  usable_cell_lines <- cell_df %>% 
    dplyr::group_by(singlet_ID, condition) %>% 
    dplyr::summarise(n = n()) %>% 
    reshape2::dcast(singlet_ID ~ condition, value.var = 'n') %>% 
    dplyr::filter(control >= min_cells, treat >= min_cells) %>% 
    .[['singlet_ID']]
  
  dat <- load_compute_CPM_stats(cur_expt, results_dir, prior_counts = globals$pca_prior_cnt, type = 'avg_cpm', use_CLs = usable_cell_lines)
  
  #get top most variable genes (by LFC)
  gene_SDs <- dat$LFC_mat %>% 
    apply(1, sd, na.rm=T)
  use_genes <- gene_SDs %>% 
    sort(decreasing = T) %>% 
    head(globals$n_highvar_genes) %>% 
    names() %>% 
    as.character()
  
  ucols <- which(apply(dat$LFC_mat[use_genes,], 2, sd) > 0)
  pca_umat <- scale(t(dat$LFC_mat[use_genes, ucols]), center = T, scale = F) %>% t()
  pc_res <- prcomp(pca_umat, center = TRUE, scale. = FALSE)
  
  write_rds(pc_res, file.path(out_dir, 'LFC_PCA_avgcpm.rds'))
}
for (cur_expt in sc_DE_meta) {
  print(sprintf('Processing experiment %s', cur_expt$expt_name))
  compute_lfc_pca(cur_expt)
}

