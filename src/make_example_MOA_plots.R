make_example_MOA_plots <- function() {
  
  #PARAMS
  n_label <- 15
  min_LFC <- 0
  n_lab_GSEA <- 6
  top_genes_GSEA <- 50
  
  ## EVEROLIMUS 24hr
  cur_expt <- sc_DE_meta$Everolimus_expt10
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  limma_res <-  read_rds(file.path(out_dir, 'limma_res.rds'))
  limma_res$res_avg %>%
    dplyr::mutate(is_sig = abs(logFC) > min_LFC, adj.P.Val < globals$q_thresh) %>%
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(pch = 21, fill = 'black', color = 'white', stroke = 0.1, size = 1.5, alpha = 1) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_label_repel(data = . %>% filter(is_sig) %>% arrange(dplyr::desc(abs(logFC))) %>% head(n_label),
                     aes(label = Gene),
                     size = 2.5) +
    xlab('Avg logFC') +
    cdsr::theme_Publication()
  ggsave(file.path(fig_dir, sprintf('%s_avg_volcano.png', cur_expt$expt_name)),
         width = 4, height = 3.5)
  
  gene_stat <- limma_res$res_avg$logFC %>% set_names(limma_res$res_avg$Gene)
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = top_genes_GSEA, n_lab_per = n_lab_GSEA, lab_size = 3)
  ggsave(file.path(fig_dir, sprintf('%s_avg_GSEA.png', cur_expt$expt_name)),
         width = 6.5, height = 3)
  
  
  ## GEMCITABINE 24hr
  cur_expt <- sc_DE_meta$Gemcitabine_expt10
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  limma_res <-  read_rds(file.path(out_dir, 'limma_res.rds'))
  limma_res$res_avg %>%
    dplyr::mutate(is_sig = abs(logFC) > min_LFC, adj.P.Val < globals$q_thresh) %>%
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(pch = 21, fill = 'black', color = 'white', stroke = 0.1, size = 1.5, alpha = 1) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_label_repel(data = . %>% filter(is_sig) %>% arrange(dplyr::desc(abs(logFC))) %>% head(n_label),
                     aes(label = Gene),
                     size = 2.5) +
    xlab('Avg logFC') +
    cdsr::theme_Publication()
  ggsave(file.path(fig_dir, sprintf('%s_avg_volcano.png', cur_expt$expt_name)),
         width = 4, height = 3.5)
  
  gene_stat <- limma_res$res_avg$logFC %>% set_names(limma_res$res_avg$Gene)
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = top_genes_GSEA, n_lab_per = n_lab_GSEA, lab_size = 3)
  ggsave(file.path(fig_dir, sprintf('%s_avg_GSEA.png', cur_expt$expt_name)),
         width = 6.5, height = 3)
  
  
  ## BORTEZOMIB 6hr
  cur_expt <- sc_DE_meta$bortezomib_6hr_expt1
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  limma_res <-  read_rds(file.path(out_dir, 'limma_res.rds'))
  limma_res$res_avg %>% 
    dplyr::mutate(is_sig = abs(logFC) > min_LFC, adj.P.Val < globals$q_thresh) %>% 
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(pch = 21, fill = 'black', color = 'white', stroke = 0.1, size = 1.5, alpha = 1) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_label_repel(data = . %>% filter(is_sig) %>% arrange(dplyr::desc(abs(logFC))) %>% head(n_label),
                     aes(label = Gene),
                     size = 2.5) +
    xlab('Avg logFC') +
    cdsr::theme_Publication() 
  ggsave(file.path(fig_dir, sprintf('%s_avg_volcano.png', cur_expt$expt_name)),
         width = 4, height = 3.5)
  
  gene_stat <- limma_res$res_avg$logFC %>% set_names(limma_res$res_avg$Gene)
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = top_genes_GSEA, n_lab_per = n_lab_GSEA, lab_size = 3)
  ggsave(file.path(fig_dir, sprintf('%s_avg_GSEA.png', cur_expt$expt_name)),
         width = 6.5, height = 3)
}
