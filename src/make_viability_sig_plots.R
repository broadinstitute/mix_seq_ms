
make_viability_sig_plots <- function(cur_expt) {
  library(latex2exp)
  
  #PARAMS
  n_genes_heatmap <- 30
  gsea_dir <- 'both'
  volc_point_size <- 2
  n_label <- 10
  volc_text_size <- 2.5
  sat_quantile <- 0.98 #what quantile to saturate colormap
  trametinib_example_genes <- c('EGR1', 'RRM2')
  
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  limma_res <-  read_rds(file.path(out_dir, 'limma_res.rds'))
  
  CL_df <- all_CL_features[[cur_expt$expt_name]]
  
  limma_res$res_int %>% 
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(fill = 'black', pch = 21, size = volc_point_size, alpha = 0.8, color = 'white', stroke = 0.1) +
    geom_label_repel(data = . %>% arrange(dplyr::desc(abs(logFC))) %>% head(n_label),
                     aes(label = Gene),
                     size = volc_text_size, label.padding = 0.1) +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    xlab(TeX('$\\beta_0')) +
    cdsr::theme_Publication() 
  ggsave(file.path(fig_dir, sprintf('%s_int_volcano.png', cur_expt$expt_name)),
         width = 4, height = 3.5)
  
  limma_res$res_slope %>% 
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(fill = 'black', pch = 21, size = volc_point_size, alpha = 0.8, color = 'white', stroke = 0.1) +
    geom_label_repel(data = . %>% arrange(dplyr::desc(abs(logFC))) %>% head(n_label),
                     aes(label = Gene),
                     size = volc_text_size, label.padding = 0.1) +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    xlab(TeX('$\\beta_1')) +
    cdsr::theme_Publication() 
  ggsave(file.path(fig_dir, sprintf('%s_slope_volcano.png', cur_expt$expt_name)),
         width = 4, height = 3.5)

  limma_res$res_avg %>% 
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(fill = 'black', pch = 21, size = volc_point_size, alpha = 0.8, color = 'white', stroke = 0.1) +
    geom_label_repel(data = . %>% arrange(dplyr::desc(abs(logFC))) %>% head(n_label),
                     aes(label = Gene),
                     size = volc_text_size, label.padding = 0.1) +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    xlab('Avg logFC') +
    cdsr::theme_Publication() 
  ggsave(file.path(fig_dir, sprintf('%s_avg_volcano.png', cur_expt$expt_name)),
         width = 4, height = 3.5)
    
  #make GSEA plots
  gene_stat <- limma_res$res_int$logFC %>% set_names(limma_res$res_int$Gene)
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = globals$gsea_top_n, dir = gsea_dir, n_lab_per = 5, lab_size = 3, max_chars = 30)
  ggsave(file.path(fig_dir, sprintf('%s_int_GSEA.png', cur_expt$expt_name)),
         width = 5.5, height = 2.5)
  
  gene_stat <- limma_res$res_slope$logFC %>% set_names(limma_res$res_slope$Gene)
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = globals$gsea_top_n, dir = gsea_dir, n_lab_per = 5, lab_size = 3, max_chars = 30)
  ggsave(file.path(fig_dir, sprintf('%s_slope_GSEA.png', cur_expt$expt_name)),
         width = 5.5, height = 2.5)
 
  gene_stat <- limma_res$res_avg$logFC %>% set_names(limma_res$res_avg$Gene)
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = globals$gsea_top_n, dir = gsea_dir, n_lab_per = 5, lab_size = 3, max_chars = 30)
  ggsave(file.path(fig_dir, sprintf('%s_avg_GSEA.png', cur_expt$expt_name)),
         width = 5.5, height = 2.5)
  
  #load LFC matrix
  dat <- load_compute_CPM_stats(cur_expt, results_dir, type = 'sum', prior_counts = globals$prior_cnt)
  LFC_mat <- t(dat$LFC_mat)
  
  #make scatterplots of individual gene LFC vs drug sensitivity
  make_gene_LFC_sens_scatter <- function(cur_gene, CL_df, LFC_mat) {
    cur_df <- CL_df %>% inner_join(data.frame(CCLE_ID = rownames(LFC_mat), LFC = LFC_mat[, cur_gene]))
    cur_slope <- limma_res$res_slope %>% filter(Gene == cur_gene) %>% .[['logFC']]
    cur_int <- limma_res$res_int %>% filter(Gene == cur_gene) %>% .[['logFC']]
    ggplot(cur_df, aes(sens, LFC)) + 
      geom_point(pch = 21, size = 2, alpha = 0.75, fill = 'black', color = 'white', stroke = 0.1) + 
      geom_smooth(method = 'lm') +
      geom_hline(yintercept = 0, linetype = 'dashed') + 
      geom_vline(xintercept = 0, linetype = 'dashed') +
      xlab(sprintf('%s sensitivity', cur_expt$drug_name)) + 
      ylab(sprintf('%s LFC', cur_gene)) +
      cdsr::theme_Publication()
    ggsave(file.path(fig_dir, sprintf('%s_%s_sens_LFC_scatter.png', cur_expt$drug_name, cur_gene)), width = 3.5, height = 3)
  }
  if (cur_expt$expt_name == 'Trametinib_24hr_expt3') {
    l_ply(trametinib_example_genes, make_gene_LFC_sens_scatter, CL_df, LFC_mat)
  }
  
}