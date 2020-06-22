make_cell_rep_fig <- function(cur_expt) {
  cell_info <- load_cell_info(cur_expt)
  
  counts_tab <- cell_info %>% 
    dplyr::mutate(condition = str_match(condition, '([:alnum:]+)_[12]')[,2]) %>% 
    dplyr::group_by(singlet_ID, condition) %>% 
    dplyr::summarise(n = n()) %>% 
    reshape2::acast(singlet_ID ~ condition, value.var = 'n')
  
  cord <- rowSums(counts_tab, na.rm=T) %>% 
    sort(decreasing = TRUE) %>% 
    names(.)
  
  counts_tab %>% 
    reshape2::melt() %>% 
    set_colnames(c('CCLE_ID', 'condition', 'cells')) %>% 
    dplyr::mutate(CCLE_ID = factor(CCLE_ID, levels = cord)) %>% 
    ggplot(aes(CCLE_ID, cells, fill = condition)) + 
    geom_bar(stat = 'identity') +
    theme_Publication() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab('Cell Line') +
    # scale_fill_Publication() 
    scale_fill_manual(values = c('darkgray', 'darkred')) 
    
  ggsave(file.path(fig_dir, sprintf('%s_cell_rep.png', cur_expt$expt_name)),
         width = 5, height = 3.5)
}
