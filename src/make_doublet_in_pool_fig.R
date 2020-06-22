#make figure showing example of doublet out-ot-pool classification
make_doublet_in_pool_fig <- function() {
  cur_expt <- sc_expts$Untreated_6hr_expt1
  df <- read_csv('/Volumes/GoogleDrive/My Drive/Project Apollo/Apollo drugs /scRNAseq Pilot 001/Untreated_6hr/doublet_out_of_pool_class.csv')
  act_df <- load_cell_info(cur_expt, QC_filter = F)
  stopifnot(nrow(df) == nrow(act_df))
  in_pool <- unique(act_df %>% dplyr::filter(cell_quality == 'normal') %>% .[['singlet_ID']])
  df$cell_quality <- act_df$cell_quality
  df$doublet_dev_imp <- act_df$doublet_dev_imp
  df %<>% dplyr::filter(cell_quality %in% c('doublet', 'normal'))
  df %<>% dplyr::mutate(both_in_pool = (doublet_CL1 %in% in_pool) & (doublet_CL2 %in% in_pool))
  df %<>% dplyr::mutate(cell_quality = plyr::revalue(cell_quality, c(normal = 'singlet')))
  
  cnts <- df %>% 
    dplyr::group_by(cell_quality, both_in_pool) %>% 
    dplyr::summarise(n = n())
  fracs <- df %>% 
    dplyr::group_by(cell_quality) %>% 
    dplyr::summarise(pp = mean(both_in_pool)*100,
              tot = n())
  cnts %<>% full_join(fracs)
  print(fracs)
  ggplot(cnts, aes(cell_quality)) + 
    geom_bar(aes(y = n, group = both_in_pool, fill = both_in_pool), stat = 'identity') +
    geom_text(data = fracs, size=3.5, aes(y = tot, label = paste0(round(pp, 1), ' %')), vjust = 0, size = 3.5) +
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank()) +
    ylab('Number of cells') +
    guides(fill = guide_legend(label.theme = element_text(size = 10))) +
    coord_cartesian(ylim = c(0, max(fracs$tot)), clip = 'off')
  ggsave(file.path(fig_dir, sprintf('%s_doublet_in_pool_bar.png', cur_expt$expt_name)), width = 3, height = 2.75)
  
}

