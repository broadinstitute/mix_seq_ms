make_nutlin_volcanos_heatmap <- function() {
  
  #PARAMS
  n_top_genes <- 40
  sat_quantile <- 0.99
  LFC_thresh <- 1
  n_label <- 15
  

  cur_expt <- sc_DE_meta$idasanutlin_24hr_expt1
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  
  dat <- load_collapsed_profiles(cur_expt, results_dir, type = 'sum') 
  
  #compute DE for TP53 WT cell lines
  cur_samps <- dat$sample_info %>% 
    dplyr::filter((CCLE_ID %in% globals$TP53_WT_cls_pool22)) %>% 
    .[['sample_name']]
  WT_DE <- fit_viability_models(dat$profile_mat[, cur_samps], 
                                dat$sample_info[match(cur_samps, dat$sample_info$sample_name),], 
                                sensitivity_df = NULL, prior_cnt = globals$prior_cnt)$res_avg 
  
  #compute DE for TP53 null cell lines
  cur_samps <- dat$sample_info %>% 
    dplyr::filter(!(CCLE_ID %in% globals$TP53_WT_cls_pool22)) %>% 
    .[['sample_name']]
  null_DE <- fit_viability_models(dat$profile_mat[, cur_samps], 
                                dat$sample_info[match(cur_samps, dat$sample_info$sample_name),], 
                                sensitivity_df = NULL, prior_cnt = globals$prior_cnt)$res_avg 
  
  
  #make volcano plots of WT and Null DE analysis
  max_val <- max(max(WT_DE$logFC), max(null_DE$logFC))
  min_val <- min(min(WT_DE$logFC), min(null_DE$logFC))
  g_WT <- WT_DE %>% 
    dplyr::mutate(is_sig = abs(logFC) > LFC_thresh, adj.P.Val < globals$q_thresh) %>% 
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(data = . %>% filter(!is_sig), fill = 'black', pch = 21, size = 1.5, color = 'gray', stroke = 0.1, alpha = 0.5) +
    geom_point(data = . %>% filter(is_sig), fill = 'darkred', pch = 21, size = 1.5, alpha = 0.8, color = 'gray', stroke = 0.1) +
    geom_label_repel(data = . %>% 
                       dplyr::filter(is_sig) %>% 
                       dplyr::arrange(dplyr::desc(abs(logFC))) %>% head(n_label),
                     aes(label = Gene),
                     size = 2.5, label.padding = 0.2) +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_vline(xintercept = c(-LFC_thresh, LFC_thresh), linetype = 'dashed', size = 0.25) + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    theme_Publication() +
    xlim(min_val, max_val) +
    ggtitle('TP53 WT')
  
  
  g_Null <- null_DE %>% 
    dplyr::mutate(is_sig = abs(logFC) > LFC_thresh, adj.P.Val < globals$q_thresh) %>% 
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(data = . %>% filter(!is_sig), fill = 'black', pch = 21, size = 1.5, color = 'gray', stroke = 0.1, alpha = 0.5) +
    geom_point(data = . %>% filter(is_sig), fill = 'darkred', pch = 21, size = 1.5, alpha = 0.8, color = 'gray', stroke = 0.1) +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_vline(xintercept = c(-LFC_thresh, LFC_thresh), linetype = 'dashed', size = 0.25) + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    theme_Publication() +
    xlim(min_val, max_val) +
    ggtitle('TP53 Mutant')
  ggsave(file.path(fig_dir, 'nutlin_TP53WT_avg_volcano.png'), width = 4, height = 3.5, plot = g_WT)
  ggsave(file.path(fig_dir, 'nutlin_TP53Null_avg_volcano.png'), width = 4, height = 3.5, plot = g_Null)
  
  
  #now create heatmap plot
  CL_df <- all_CL_features[[cur_expt$expt_name]]
  CL_df %<>% 
    dplyr::mutate(TP53_status = 'WT',
                      TP53_status = ifelse(TP53_DAM_HOT > 0 | TP53_DAM_MUT > 0, 'Mutant', TP53_status),
                      TP53_status = factor(TP53_status, levels = c('WT', 'Mutant'))) %>% 
      dplyr::select(CCLE_ID, AUC = AUC_avg, TP53_status)
  
  sc_dat <- load_compute_CPM_stats(cur_expt, results_dir, type = 'sum', prior_counts = globals$prior_cnt)
  top_genes <- WT_DE %>% 
    dplyr::filter(adj.P.Val < globals$q_thresh) %>% 
    dplyr::arrange(dplyr::desc(abs(logFC))) %>% 
    head(n_top_genes) %>% 
    .[['Gene']]
  
  
  CL_ord <- CL_df %>% 
    dplyr::filter(CCLE_ID %in% colnames(sc_dat$LFC_mat)) %>% 
    dplyr::arrange(TP53_status) %>% 
    .[['CCLE_ID']]
  
  CL_df %<>% 
    dplyr::mutate(short_name = str_match(CCLE_ID, '^([:alnum:]+)_')[,2]) 
  name_remap <- with(CL_df %>% dplyr::filter(CCLE_ID %in% colnames(sc_dat$LFC_mat)),
                     short_name %>% set_names(CCLE_ID))
  colnames(sc_dat$LFC_mat) <- plyr::revalue(colnames(sc_dat$LFC_mat), replace = name_remap)
  cmax <- max(abs(quantile(sc_dat$LFC_mat[top_genes,], c(1-sat_quantile, sat_quantile), na.rm=T)))
  make_LFC_heatmap(sc_dat$LFC_mat, top_genes, 
                   CL_ann_df = CL_df %>% 
                     dplyr::filter(short_name %in% colnames(sc_dat$LFC_mat)) %>% 
                     dplyr::mutate(sens = 1 - AUC) %>% 
                     dplyr::select(`nutlin sens.` = sens, `TP53 status` = TP53_status, short_name) %>% 
                     column_to_rownames(var = 'short_name') %>% 
                     .[colnames(sc_dat$LFC_mat),,drop=FALSE],
                   CL_list = CL_ord %>% plyr::revalue(name_remap),
                   color_lims = c(-cmax, cmax),
                   cluster_rows = TRUE,
                   transpose = TRUE,
                   fontsize_col = 6, 
                   fontsize_row = 6, 
                   fontsize = 8,
                   treeheight_col = 0, 
                   treeheight_row = 10,
                   filename = file.path(fig_dir, sprintf('%s_LFC_heatmap.png', cur_expt$expt_name)),
                   show_rownames = TRUE,
                   width = 5.5, height = 3,
                   annotation_legend = TRUE)
  
  gene_stat <- WT_DE$logFC %>% set_names(WT_DE$Gene)
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = 50, n_lab_per = 5, lab_size = 3)
  ggsave(file.path(fig_dir, 'nutlin_P53WT_GSEA_stem.png'), width = 5.5, height = 2.5)

}