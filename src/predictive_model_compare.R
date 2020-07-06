
run_predictive_model_compare <- function() {
  #PARAMS
  min_cells <- 5 #minimum cells per condition to include cell line in analysis
  min_R2 <- -1
  n_top_feat <- 1000
  random_seed <- 1
  kfold <- 10
  max_AUC <- 1.5
  
  feat_dat <- list(
    GE = list(
      data.name="depmap-rnaseq-expression-data-ccd0",
      data.version = 14,
      data.file = 'CCLE_depMap_19Q3_TPM',
      transpose = F),
    # CN = list(data.name = 'depmap-wes-cn-data-97cc', 
    #           data.version = 18, 
    #           data.file = 'public_19Q3_gene_cn',
    #           transpose = F),
    MUT_DAM = list(
      data.name = 'depmap-mutation-calls-9a1a',
      data.version = 12,
      data.file = 'damaging_mutation',
      transpose = F),
    MUT_HOT = list(
      data.name = 'depmap-mutation-calls-9a1a',
      data.version = 12,
      data.file = 'hotspot_mutation',
      transpose = F)) %>% 
    taigr::load.all.from.taiga() 
  
  doMC::registerDoMC(cores = 3)
  
  targ_dsets <- sc_DE_meta[str_match(names(sc_DE_meta), '(expt[0-9]+)$')[,2] %in% c('expt1', 'expt3', 'expt10')]
  all_pmod_results <- ldply(targ_dsets, function(cur_expt) {
    print(sprintf('Processing expt %s', cur_expt$expt_name))
    out_dir <- file.path(results_dir, cur_expt$expt_name)
  
    cell_df <- load_cell_info(cur_expt, QC_filter = TRUE) %>% 
      tidyr::separate(col = 'condition', into = c('condition', 'batch'), sep = '_')
    cell_counts <- cell_df %>% 
      dplyr::group_by(singlet_ID, condition) %>% 
      dplyr::summarise(n = n()) %>% 
      reshape2::dcast(singlet_ID ~ condition, value.var = 'n')
    usable_CLs <- cell_counts %>% 
      dplyr::filter(control >= min_cells, treat >= min_cells) %>% 
      .[['singlet_ID']] %>% 
      celllinemapr::ccle.to.arxspan()
    
    CL_df <- all_CL_features[[cur_expt$expt_name]]

    dat <- load_compute_CPM_stats(cur_expt, results_dir, prior_counts = globals$pca_prior_cnt, type = 'avg_cpm')
    colnames(dat$LFC_mat) <- celllinemapr::ccle.to.arxspan(colnames(dat$LFC_mat))
    
    y <- CL_df$AUC_avg %>% set_names(CL_df$DEPMAP_ID)
    y <- y[!is.na(y)]
    y <- pmin(y, max_AUC)
    
    #MODEL BASED PREDICTION USING LFC  
    X <- list(
      LFC = list(data = t(dat$LFC_mat))
    )
    input_data <- proc_input_data(X, y)
    fids <- caret::createFolds(input_data$y, k = kfold, list = FALSE) %>% 
      set_names(names(input_data$y))
    n_LFC_CLs <- length(input_data$y)
    LFC_CLs <- names(input_data$y)
    LFC_results <- get_rf_xv_results(input_data, kfold, fids = fids, n_top_feat = n_top_feat)
    
    #NOW TRY MODELS USING BASELINE FEATURES
    X <- list(
      GE = list(data = feat_dat$GE),
      MUT_DAM = list(data = feat_dat$MUT_DAM),
      MUT_HOT = list(data = feat_dat$MUT_HOT)
      # CN = list(data = feat_dat$CN)
    )
    input_data <- proc_input_data(X, y[usable_CLs])
    baseline_matched_CLs <- names(input_data$y)
    baseline_matched_results <- get_rf_xv_results(input_data, kfold, fids = fids[names(input_data$y)], n_top_feat = n_top_feat)
      
    #USING ALL AVAILABLE CELL LINES
   input_data <- proc_input_data(X, y)
   fids <- caret::createFolds(input_data$y, k = kfold, list = FALSE) %>% 
     set_names(names(input_data$y))
   baseline_all_CLs <- names(input_data$y)
   baseline_all_results <- get_rf_xv_results(input_data, kfold, fids = fids, n_top_feat = n_top_feat)
      
      pred_mod_errs <- rbind(
        data_frame(type = 'matched_sample_pred',
                   n = length(baseline_matched_CLs),
                   corr_R2 = caret::R2(baseline_matched_results$prediction, baseline_matched_results$obs, form = 'corr'),
                   trad_R2 = caret::R2(baseline_matched_results$prediction, baseline_matched_results$obs, form = 'traditional')),
        data_frame(type = 'all_sample_pred',
                   n = length(baseline_all_CLs),
                   corr_R2 = caret::R2(baseline_all_results$prediction, baseline_all_results$obs, form = 'corr'),
                   trad_R2 = caret::R2(baseline_all_results$prediction, baseline_all_results$obs, form = 'traditional')),
        data_frame(type = 'LFC_pred',
                   n = length(LFC_CLs),
                   corr_R2 = caret::R2(LFC_results$prediction, LFC_results$obs, form = 'corr'),
                   trad_R2 = caret::R2(LFC_results$prediction, LFC_results$obs, form = 'traditional'))
      ) %>% dplyr::mutate(expt = cur_expt$expt_name)
      return(pred_mod_errs)
  })
  
  write_rds(all_pmod_results, file.path(results_dir, 'all_pmod_results.rds'))
} 



make_predictive_model_figs <- function() {
  all_pmod_results <- read_rds(file.path(results_dir, 'all_pmod_results.rds'))
  all_pmod_results %<>% dplyr::mutate(
    expt = str_replace_all(expt, '_expt10', '_24hr_expt10')
  )
  used_expts <- c('Trametinib_6hr_expt1',
                  'Bortezomib_6hr_expt1',
                  'Bortezomib_24hr_expt1',
                  'Idasanutlin_6hr_expt1',
                  'Idasanutlin_24hr_expt1',
                  'Trametinib_24hr_expt3', 
                  'Dabrafenib_24hr_expt3',
                  'BRD3379_6hr_expt3',
                  'BRD3379_24hr_expt3',
                  'Afatinib_24hr_expt10',
                  'Navitoclax_24hr_expt3',
                  'AZD5591_24hr_expt10',
                  'Everolimus_24hr_expt10',
                  'Gemcitabine_24hr_expt10',
                  'Prexasertib_24hr_expt10',
                  'Taselisib_24hr_expt10',
                  'JQ1_24hr_expt10')
  
  pmod_results_wide <- 
    all_pmod_results %>% 
      dplyr::select(expt, type, trad_R2) %>% 
      tidyr::spread(type, trad_R2) %>% 
    dplyr::filter(expt %in% used_expts) %>% 
    dplyr::mutate(time = str_match(expt, '^[:alnum:]+_([:alnum:]+)')[,2],
                  drug = str_match(expt, '^([:alnum:]+)_')[,2]) 
  
  pmod_results_wide %>% 
    ggplot(aes(all_sample_pred, LFC_pred)) + 
    geom_point(aes(fill = time), pch = 21, size = 3, color = 'white', stroke = 0.2) +
    geom_abline() +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggrepel::geom_label_repel(aes(label = drug, color = time), size = 3, label.padding = 0.1) +
    # coord_cartesian(xlim = c(-0.1, 0.6), ylim = c(-0.1, 0.6)) +
    xlab('Baseline features\n(R2)') +
    ylab('Transcriptional responses\n(R2)') +
    guides(color = F) + 
    theme_Publication() +
    scale_color_Publication() +
    scale_fill_Publication()
  ggsave(file.path(fig_dir, 'pred_R2_scatter_all.png'), width = 4, height = 4)
  ggsave(file.path(fig_dir, 'pred_R2_scatter_all.pdf'), width = 4, height = 4)
  
  with(pmod_results_wide, {
    print(sum(!is.na(matched_sample_pred) & !is.na(LFC_pred)))
    print(wilcox.test(matched_sample_pred, LFC_pred, paired = TRUE))
  })
  pmod_results_wide %>% 
    ggplot(aes(matched_sample_pred, LFC_pred)) + 
    geom_point(aes(fill = time), pch = 21, size = 3, color = 'white', stroke = 0.2) +
    geom_abline() +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggrepel::geom_label_repel(aes(label = drug, color = time), size = 3, label.padding = 0.1) +
    # coord_cartesian(xlim = c(-0.1, 0.6), ylim = c(-0.1, 0.6)) +
    xlab('Baseline features\n(R2)') +
    ylab('Transcriptional responses\n(R2)') +
    guides(color = F) + 
    theme_Publication() +
    scale_color_Publication() +
    scale_fill_Publication()
  ggsave(file.path(fig_dir, 'pred_R2_scatter_matched.png'), width = 4, height = 4)
  ggsave(file.path(fig_dir, 'pred_R2_scatter_matched.pdf'), width = 4, height = 4)
  
  n_results_wide <- 
    all_pmod_results %>% 
    dplyr::select(expt, type, n) %>% 
    tidyr::spread(type, n) %>% 
    dplyr::filter(expt %in% used_expts) %>% 
    left_join(pmod_results_wide, by = 'expt', suffix = c('_n', '_R2')) %>% 
    dplyr::mutate(LFC_baseline_diff = LFC_pred_R2 - all_sample_pred_R2)

  ggplot(n_results_wide %>% dplyr::filter(time == '24hr'), aes(LFC_pred_n, LFC_baseline_diff)) + 
    geom_point(fill = 'black', pch = 21, size = 3, color = 'white', stroke = 0.2) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggrepel::geom_label_repel(aes(label = drug), size = 3, label.padding = 0.1) +
    xlab('LFC sample size') +
    ylab('Prediction accuracy diff.\n(LFC_R2 - Baseline_R2)') +
    guides(color = F) + 
    theme_Publication() 
    # scale_color_Publication() +
    # scale_fill_Publication()
  ggsave(file.path(fig_dir, 'pred_R2_diff_vs_n.png'), width = 4, height = 4)
  
}