#make figure showing distribution of cell classification types across experiments
make_cell_dist_fig <- function() {
  used_expts <- names(sc_expts)[grepl('expt[1|3|2]$', names(sc_expts))]
  used_expts %<>% union(names(sc_expts)[grepl('expt10$', names(sc_expts))])
  res <- ldply(sc_expts[used_expts], function(cur_expt) {
    ip <- load_cell_info(cur_expt, QC_filter = FALSE) %>% 
      mutate(expt = cur_expt$expt_name)
  })
  res %<>% mutate(batch = str_match(expt, '_(expt[0-9]$)')[,2],
                  expt = factor(expt, levels = res %>% dplyr::group_by(expt) %>% dplyr::summarise(n = n()) %>% arrange(n) %>% .[['expt']]))

  wierd_cell_stats <- res %>% 
    dplyr::filter(!(cell_quality %in% c('normal', 'doublet'))) %>% 
    dplyr::group_by(cell_quality) %>% 
    dplyr::summarise(n = n())
  
  res %<>% mutate(cell_quality = plyr::revalue(cell_quality, replace = c(normal = 'singlet', low_confidence = 'low_quality', empty_droplet = 'low_quality')))
  
  #distribution of cell classification
  ggplot(res, aes(expt, fill = cell_quality)) + 
    geom_histogram(stat = 'count', position = 'stack') +
    cdsr::theme_Publication() +
    ylab('Number of cells') +
    cdsr::scale_fill_Publication() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    guides(fill = guide_legend(title = element_blank(), nrow = 2))
  ggsave(file.path(fig_dir, 'cell_quality_dists.png'),
         width = 4, height = 3)
  
}
