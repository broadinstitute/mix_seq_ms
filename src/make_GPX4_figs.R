
make_GPX4_figs <- function() {
  #PARAMS
  n_label <- 5
  
  cur_expt <- sc_DE_meta$GPX4_expt2
  sprintf('Processing expt %s', cur_expt$expt_name)
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  
  seuObj <- load_sc_data(cur_expt, sc_expts = sc_expts)
  
  seuObj <- NormalizeData(object = seuObj, 
                          normalization.method = "LogNormalize", 
                          scale.factor = globals$scale_fac)
  
  #get dataframe of single-cell GPX4 expression levels
  df <- data.frame(Seurat::FetchData(object = seuObj, vars = 'GPX4'), check.names = F) %>% 
    cbind(seuObj@meta.data) %>% 
    tidyr::separate('condition', into = c('treat_cond', 'batch'), sep = '_', remove = FALSE) %>% 
    dplyr::mutate(condition = revalue(condition, 
                               c(control_1 = 'sgOR2J2',
                                 control_2 = 'sgLACZ',
                                 treat_1 = 'sgGPX4_1',
                                 treat_2 = 'sgGPX4_2')),
           condition = factor(condition, levels = c('sgLACZ', 'sgOR2J2', 'sgGPX4_1', 'sgGPX4_2')))
  
  #plot GP4 expression per guide for example cell line
  targ_CL <- 'YD38_UPPER_AERODIGESTIVE_TRACT'
  df %>% 
    dplyr::filter(singlet_ID == targ_CL) %>% 
    ggplot(aes(condition, GPX4)) + 
    geom_violin(aes(fill = condition), alpha = 0.5) +
    ggbeeswarm::geom_beeswarm(size = 0.1, alpha = 0.6) +
    guides(fill = FALSE) + 
    cdsr::theme_Publication() +
    ylab('GPX4 expression (log2CPM)') +
    theme(axis.text.x = element_text(angle = 80, hjust = 1),
          axis.title.x = element_blank(),
          title = element_text(size = 8),
          axis.title.y = element_text(size = 10)) +
    ggtitle(targ_CL)
  ggsave(file.path(fig_dir, sprintf('GPX4_expression_%s.png', targ_CL)), width = 3.5, height = 3)
  
  #scatterplot of average expression in control vs treat across cell lines
  df %>% 
    dplyr::group_by(singlet_ID, treat_cond) %>% 
    dplyr::summarise(GPX4 = mean(GPX4)) %>% 
    reshape2::dcast(singlet_ID ~ treat_cond, value.var = 'GPX4') %>% 
    ggplot(aes(control, treat)) + 
    geom_point(pch = 21, alpha = 0.75, fill = 'black', color = 'white', size = 4) + 
    geom_point(data = . %>% filter(singlet_ID == targ_CL), 
               pch = 21, alpha = 1, fill = 'red', color = 'white', size = 4) + 
    geom_abline() +
    xlim(2, 5.75) + ylim(2, 5.75) + 
    xlab('Control GPX4 expression\n(logCPM)') + 
    ylab('Post-KO GPX4 expression\n(logCPM)') + 
    cdsr::theme_Publication()
  ggsave(file.path(fig_dir, 'GPX4_expression_scatter.png'), width = 3.5, height = 3)
  
  #barplot of GPX4 LFC across cell lines
  aa <- df %>% 
    dplyr::group_by(singlet_ID, treat_cond) %>% 
    dplyr::summarise(GPX4 = mean(GPX4)) %>% 
    reshape2::dcast(singlet_ID ~ treat_cond, value.var = 'GPX4') %>% 
    dplyr::mutate(LFC = treat - control,
           singlet_ID = str_match(singlet_ID, '([:alnum:]+)_')[,2])
  aa %<>% dplyr::mutate(singlet_ID = factor(singlet_ID, levels = aa %>% arrange((LFC)) %>% .[['singlet_ID']])) 
  ggplot(aa, aes(singlet_ID, LFC)) + 
    geom_bar(stat = 'identity', width = 0.7, fill = 'darkred') +
    cdsr::theme_Publication() +
    geom_hline(yintercept = 0, lwd = 1) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab('Cell Line') + 
    ylab("GPX4 logFC") +
    ylim(-1.75, 1.75)
  ggsave(file.path(fig_dir, 'GPX4_LFC_bar.png'), width = 3.5, height = 3.)
  
  
  #average response
  limma_res <-  read_rds(file.path(out_dir, 'limma_res.rds'))
  
  limma_res$res_avg %>% 
    dplyr::mutate(is_sig = adj.P.Val < globals$q_thresh) %>% 
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(fill = 'black', color = 'white', stroke = 0.1, pch = 21, size = 2, alpha = 0.75) +
    geom_label_repel(data = . %>% 
                       dplyr::arrange(dplyr::desc(abs(logFC))) %>% head(10), 
                     aes(label = Gene),
                     size = 2.5,
                     label.padding = 0.15) +
    geom_point(data = . %>% dplyr::filter(Gene == 'GPX4'), color = 'red', size = 2, alpha = 1) +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    cdsr::theme_Publication() 
  ggsave(file.path(fig_dir, 'GPX4_avg_volcano.png'),
         width = 3.5, height = 3.5)
  
  #Now compare response in dependent vs non-dependent lines
  #load GPX4 dependency data
  gene.dep <- load.from.taiga(data.name='avana-public-19q3-0900', data.version=5, data.file='gene_dependency') %>% 
    cdsr::extract_hugo_symbol_colnames() %>% 
    cdsr::map_arxspan_to_ccle()
  GPX4_dep <- gene.dep[, 'GPX4', drop=F] %>%
    as.data.frame() %>%
    rownames_to_column(var = 'CCLE_ID')
  
  dat <- load_collapsed_profiles(cur_expt, results_dir, type = 'sum')
  all_CLs <- unique(dat$sample_info$CCLE_ID)
  
  #test association between GPX4 dependency status and transcriptional response
  dep_status <- GPX4_dep %>% 
    dplyr::filter(CCLE_ID %in% all_CLs) 
  dep_status <- (dep_status$GPX4 > 0.5) %>% set_names(dep_status$CCLE_ID)
  res <- run_DE_association_test(cur_expt, results_dir, dep_status, prior_cnt = globals$prior_cnt)
  res %>% 
    dplyr::mutate(is_sig = adj.P.Val < globals$q_thresh) %>% 
    ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(color = 'black', size = 1, alpha = 0.5) +
    geom_point(data = . %>% filter(Gene == 'GPX4'), color = 'red', size = 2, alpha = 1) +
    geom_label_repel(data = . %>% 
                       dplyr::arrange(dplyr::desc(abs(logFC))) %>% 
                       head(n_label),
                     aes(label = Gene),
                     size = 2.5) +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    cdsr::theme_Publication() 
  ggsave(file.path(fig_dir, 'GPX4_dep_diff_volcano.png'),
         width = 3.5, height = 3)
  
  
  #make GSEA plot of avg DE genes
  gene_stat <- limma_res$res_avg$logFC %>% set_names(limma_res$res_avg$Gene)
  gene_stat['GPX4'] <- NA
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = 50, n_lab_per = 5, lab_size = 2.5)
  ggsave(file.path(fig_dir, 'GPX4_avg_GSEA.png'),
         width = 5.5, height = 2.5)
  
  gene_stat <- res$logFC %>% set_names(res$Gene)
  gene_stat['GPX4'] <- NA
  make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = 50, n_lab_per = 25)
  ggsave(file.path(fig_dir, 'GPX4_slope_GSEA.png'),
         width = 4, height = 3)
}