make_PC_AUC_compare_fig <- function() {
    n_pcs_analyze <- 3
    used_expts <- c('Trametinib_6hr_expt1',
                    'Bortezomib_6hr_expt1',
                    'Bortezomib_24hr_expt1',
                    'Idasanutlin_6hr_expt1',
                    'Idasanutlin_24hr_expt1',
                    'Trametinib_24hr_expt3', 
                    'Dabrafenib_24hr_expt3',
                    'Navitoclax_24hr_expt3', 
                    'BRD3379_6hr_expt3',
                    'BRD3379_24hr_expt3',
                    'Afatinib_expt10',
                    'AZD5591_expt10',
                    'Everolimus_expt10',
                    'Gemcitabine_expt10',
                    'Prexasertib_expt10',
                    'Taselisib_expt10',
                    'JQ1_expt10')
    
    print(sprintf('Not using the following expts:'))
    print(setdiff((names(sc_DE_meta)), (used_expts)))
    
    elist <- sc_DE_meta[sapply(sc_DE_meta, function(x) {x$expt_name %in% used_expts})]
    all_CL_PCs <- llply(elist, function(cur_expt) {
      if (!(cur_expt$expt_name %in% used_expts)) {return(NULL)}
      print(sprintf('Processing expt %s', cur_expt$expt_name))
      out_dir <- file.path(results_dir, cur_expt$expt_name)
      #load PCA res
      pca_res <- read_rds(file.path(out_dir, 'LFC_PCA_avgcpm.rds'))
      
      CL_df <- pca_res$rotation[, seq(n_pcs_analyze), drop = F] %>% 
        as.data.frame() %>% 
        rownames_to_column(var = 'CCLE_ID') %>% 
        mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(CCLE_ID)) %>% 
        left_join(all_CL_features[[cur_expt$expt_name]] %>% dplyr::select(-CCLE_ID), by = 'DEPMAP_ID')
       return(CL_df)
    }) 

  PC_PRISM_cors <- all_CL_PCs %>% 
    ldply(function(df) {
      aa <- psych::corr.test(df[, grepl('^PC', colnames(df)), drop = F], df[, 'sens', drop=F], adjust = 'none') 
      full_join(aa$r %>% as.data.frame() %>% rownames_to_column(var = 'PC'), 
                aa$p %>% as.data.frame() %>% rownames_to_column(var = 'PC'),
                by = 'PC', suffix = c('_r', '_p')) %>% 
        cbind(aa$se %>% set_colnames(paste0(colnames(.), '_se'))) 
    }, .id = 'expt') %>% 
    mutate(expt = str_replace_all(expt, '_expt10', '_24hr_expt10')) %>% 
    mutate(sens_rmag = abs(sens_r),
           drug = str_match(expt, '^([:alnum:]+)_')[,2],
           treat = str_replace_all(as.character(expt), '_expt[0-9]$', ''),
           expt_name = str_replace(expt, '_expt[0-9]+$', ''))

  PC_PRISM_cors %>% 
    dplyr::filter(PC == 'PC1') %>% 
    dplyr::mutate(q = p.adjust(sens_p, method = 'BH')) %>% 
    dplyr::mutate(expt_name = factor(expt_name, levels = PC_PRISM_cors %>% 
                                       dplyr::distinct(expt_name, .keep_all=T) %>% 
                                       dplyr::arrange(dplyr::desc(sens_rmag)) %>%
                                       .[['expt_name']])) %>% 
    ggplot(aes(expt_name, sens_rmag)) + 
    geom_bar(aes(fill = q < 0.1), stat = 'identity') +
    cdsr::theme_Publication() +
    ylab('PC1-sensitivity\ncorrelation magnitude') +
    scale_fill_manual(values = c(`FALSE` = 'darkgray', `TRUE` = 'darkred')) +
    guides(fill = guide_legend(title = 'FDR < 0.1')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank()) +
    ylim(0, 0.8)
  ggsave(file.path(fig_dir, 'all_PC1_sens_corrs.png'), width = 4, height = 3.5)

  PC_PRISM_cors %>% 
    dplyr::filter(PC == 'PC2') %>% 
    dplyr::mutate(q = p.adjust(sens_p, method = 'BH')) %>% 
    dplyr::mutate(expt_name = factor(expt_name, levels = PC_PRISM_cors %>% 
                                       dplyr::distinct(expt_name, .keep_all=T) %>% 
                                       dplyr::arrange(dplyr::desc(sens_rmag)) %>%
                                       .[['expt_name']])) %>% 
    ggplot(aes(expt_name, sens_rmag)) + 
    geom_bar(aes(fill = q < 0.1), stat = 'identity') +
    cdsr::theme_Publication() +
    ylab('PC2-sensitivity\ncorrelation magnitude') +
    scale_fill_manual(values = c(`FALSE` = 'darkgray', `TRUE` = 'darkred')) +
    guides(fill = guide_legend(title = 'FDR < 0.1')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank())+
    ylim(0, 0.8)
  ggsave(file.path(fig_dir, 'all_PC2_sens_corrs.png'), width = 4, height = 3.5)
 
  PC_PRISM_cors %>% 
    dplyr::filter(PC == 'PC3') %>% 
    dplyr::mutate(q = p.adjust(sens_p, method = 'BH')) %>% 
    dplyr::mutate(expt_name = factor(expt_name, levels = PC_PRISM_cors %>% 
                                       dplyr::distinct(expt_name, .keep_all=T) %>% 
                                       dplyr::arrange(dplyr::desc(sens_rmag)) %>%
                                       .[['expt_name']])) %>% 
    ggplot(aes(expt_name, sens_rmag)) + 
    geom_bar(aes(fill = q < 0.1), stat = 'identity') +
    cdsr::theme_Publication() +
    ylab('PC3-sensitivity\ncorrelation magnitude') +
    scale_fill_manual(values = c(`FALSE` = 'darkgray', `TRUE` = 'darkred')) +
    guides(fill = guide_legend(title = 'FDR < 0.1')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank())+
    ylim(0, 0.8)
  ggsave(file.path(fig_dir, 'all_PC3_sens_corrs.png'), width = 4, height = 3.5)
  
}
