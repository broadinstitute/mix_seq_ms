make_rel_abundance_plot <- function() {
  
  used_expts <- c(
    'Idasanutlin_24hr_expt1',
    'Idasanutlin_6hr_expt1',
    'Trametinib_6hr_expt1',
    'Trametinib_24hr_expt3', 
    'Bortezomib_24hr_expt1',
    'Bortezomib_6hr_expt1',
    'Navitoclax_24hr_expt3',
    'Dabrafenib_24hr_expt3',
    'BRD3379_24hr_expt3',
    'BRD3379_6hr_expt3',
    'Afatinib_expt10',
    'AZD5591_expt10',
    'Everolimus_expt10',
    'Taselisib_expt10',
    'Gemcitabine_expt10',
    'Prexasertib_expt10',
    'JQ1_expt10')
  
  pseudo_cnt <- 1 #add this many 'fake' cells to stabilize lfc estimates
  
  expt_list <- sc_DE_meta
  names(expt_list) <- lapply(expt_list, function(x) x$expt_name)
  res <- llply(expt_list, function(cur_expt) {
    print(sprintf('Processing expt %s', cur_expt$expt_name))
    out_dir <- file.path(results_dir, cur_expt$expt_name)
    
    #load cell cycle
    cell_info <- load_cell_info(cur_expt, QC_filter = TRUE)
    
    counts_tab <- cell_info %>% 
      dplyr::group_by(singlet_ID, condition) %>% 
      dplyr::summarise(n = n()) %>% 
      reshape2::acast(singlet_ID ~ condition, value.var = 'n') %>% 
      magrittr::add(pseudo_cnt)
    rel_probs <- counts_tab %>% 
      scale(center = F, scale = colSums(., na.rm=T)) %>% 
      as.data.frame() 
    LFC_df <- data_frame(CCLE_ID = rownames(rel_probs),
                      LFC = log2(rowMeans(rel_probs[, grepl('treat', colnames(rel_probs)), drop = FALSE])) - log2(rowMeans(rel_probs[, grepl('control', colnames(rel_probs)), drop = FALSE])))
    
    CL_df <- LFC_df %>%
      dplyr::mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(CCLE_ID)) %>% 
      dplyr::select(DEPMAP_ID, LFC) %>% 
      left_join(all_CL_features[[cur_expt$expt_name]], by = 'DEPMAP_ID') 
    
    return(CL_df)
  })
  
  LFC_cor <- ldply(res, function(df) {
    cor.test(df$sens, df$LFC, use = 'pairwise.complete.obs') %>% 
      broom::tidy()
  }, .id = 'expt') %>% 
    dplyr::filter(expt %in% used_expts) %>% 
    dplyr::mutate(expt = str_replace_all(expt, 'expt10', '24hr_expt10'),
                  time = str_match(expt, '_([0-9]+hr)_')[,2],
                  expt = str_replace_all(expt, '_expt[0-9]+$', ''))  
  ggplot(LFC_cor, aes(estimate, -log10(p.value))) + 
    geom_point(aes(fill = time), pch = 21, size = 3, color = 'white', stroke = 0.1) + 
    geom_text_repel(aes(label = expt), size = 2.5) +
    cdsr::scale_fill_Publication() +
    xlab('Sensitivity-Cell count\ncorrelation') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    cdsr::theme_Publication()
  ggsave(file.path(fig_dir, 'rel_abundance_corr_volcano.png'), width = 4, height = 4)

  
  cc_arrest_stats <- read_rds(file.path(results_dir, 'CC_arrest_stats.rds'))
  comb <- left_join(cc_arrest_stats, LFC_cor)  %>% 
    dplyr::mutate(expt = str_replace(expt, '_[0-9]+hr', ''))
  ggplot(comb, aes(estimate, G1_cor)) + 
    geom_point(pch = 21, size = 3, fill = 'black', color = 'white', stroke = 0.2) +
    geom_errorbar(aes(ymax = G1_uci, ymin = G1_lci), lwd = 0.2) +
    geom_errorbarh(aes(xmax = conf.high, xmin = conf.low), lwd = 0.2) +
    geom_text_repel(aes(label = expt), size = 2.5) +
    xlab('Sensitivity/-Cell count\ncorrelation') +
    ylab('Sensitivity-G0/G1-arrest\ncorrelation') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    cdsr::theme_Publication()
    ggsave(file.path(fig_dir, 'rel_abundance_G1arrest_scatter.png'), width = 4.5, height = 4)
}