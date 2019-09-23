make_cell_quality_example_fig <- function(cur_expt) {
  cell_df <- load_cell_info(cur_expt, QC_filter = FALSE)
  cell_df %<>% dplyr::mutate(mod_fit = singlet_dev + doublet_dev_imp)
  cell_df %<>% dplyr::mutate(cell_quality = plyr::revalue(cell_quality, c(normal = 'singlet', empty_droplet = 'low_quality')))
  ggplot(cell_df, aes(mod_fit, doublet_dev_imp)) + 
    geom_point(aes(fill = cell_quality), pch = 21, size = 2, color = 'white', stroke = 0.2) +
    cdsr::theme_Publication() +
    guides(fill = guide_legend(title = element_blank(), override.aes = list(size = 3))) +
    ylab('Doublet model improvement') +
    xlab('Best model fit')
  ggsave(file.path(fig_dir, sprintf('%s_cell_quality_scatter.png', cur_expt$expt_name)), width = 4, height = 3.75)
}