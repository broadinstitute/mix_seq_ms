#make figure showing SNP classification LLs for all cells in an example experiment
make_singlet_SNP_heatmap <- function() {
  cur_expt <- sc_expts$Untreated_6hr_expt1
  
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  all_stats <- read_csv('/Volumes/GoogleDrive/My Drive/Project Apollo/Apollo drugs /scRNAseq Pilot 001/Untreated_6hr/classification_stats.csv')
  cl_df <- load_cell_info(cur_expt)
  
  df <- all_stats %>% 
    dplyr::select(ref_ID = singlet_ID, coef, dev_ratio, dev_ratio_z, barcode) %>% 
    filter(barcode %in% cl_df$barcode) %>% 
    left_join(cl_df, by = 'barcode') %>% 
    arrange(singlet_ID) 
  
  df %<>% mutate(barcode = factor(barcode, levels = df %>% .[['barcode']] %>% unique()),
                 name_short = str_match(ref_ID, '(^[:alnum:]+)_')[,2])
  
  ggplot(df, aes(barcode, name_short)) + 
    geom_tile(aes(fill = dev_ratio)) + 
    scale_fill_gradient(low = "white", high = 'darkblue', name = 'Normalized LL\n(Deviance ratio)') +
    cdsr::theme_Publication(legend_bottom = FALSE) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6)) +
    xlab(sprintf('Cells (n = %d)', length(unique(df$barcode)))) +
    ylab('Reference cell line')  +
    guides(fill = guide_colorbar(barwidth = 0.75, barheight = 6, 
                                 title.theme = element_text(size = 11), 
                                 label.theme = element_text(size = 9))) +
    coord_cartesian(clip = 'off')
  ggsave(file.path(fig_dir, 'SNP_heatmap_example.png'), width = 5, height = 3.)
}
