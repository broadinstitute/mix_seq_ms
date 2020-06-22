make_demux_SNP_comparison_fig <- function() {
  #PARAMS
  doublet_prob_thresh <- 0.5
  
  #linear fit of doublet fraction to loading density
  #taken from supplementary figure 1a (https://media.nature.com/original/nature-assets/ncomms/2017/170116/ncomms14049/extref/ncomms14049-s1.pdf)
  rep_slope = 0.0067/1000
  rep_offset = 0.015
  
  used_expts <- names(sc_expts)[grepl('expt[1|3]$', names(sc_expts))]
  all_res <- ldply(sc_expts[used_expts], function(cur_expt) {
    our <- load_cell_info(cur_expt, QC_filter = FALSE) %>% 
            dplyr::mutate(our_call = ifelse(cell_quality == 'normal', 'singlet', NA),
                          our_call = ifelse(cell_quality == 'doublet', 'doublet', our_call))
    expt_short <- str_replace(cur_expt$expt_name, '_expt[1-3]', '')
    dem_path <- file.path('/Volumes/GoogleDrive/My Drive/demuxlet_tests', paste0(expt_short, '.best'))
    if (file.exists(dem_path) & cur_expt$expt_batch == 'expt1') {
      dem <- read_tsv(dem_path) %>% 
            dplyr::mutate(dem_call = ifelse(grepl('^SNG', BEST), 'singlet', NA),
                          dem_call = ifelse(grepl('^DBL', BEST), 'doublet', dem_call))
      
      #combine into one dataset
      comb <- full_join(our, dem, by = c('barcode' = 'BARCODE'))
    } else {
      comb <- our
      print(sprintf('Missing data for %s', cur_expt$expt_name))
    }
    return(comb %>% dplyr::mutate(expt = cur_expt$expt_name))
  })
  all_res %<>% dplyr::mutate(
       FN = our_call == 'singlet' & dem_call == 'doublet',
       FP = our_call == 'doublet' & dem_call == 'singlet',
       doublet_agree = our_call == dem_call,
       singlet_agree = SNG.1ST == singlet_ID,
       singlet_agree = ifelse(our_call == 'singlet' & dem_call == 'singlet', singlet_agree, NA))
  
  expt_stats <- all_res %>% 
    dplyr::group_by(expt) %>% 
    dplyr::summarise(n_cells = n(),
              tot_doublets_ours = sum(cell_quality == 'doublet', na.rm=T),
              tot_singlets_ours = sum(cell_quality == 'normal', na.rm=T),
              tot_doublets_dem = sum(dem_call == 'doublet', na.rm=T),
              tot_doublets_dem_noLQ = sum(dem_call == 'doublet' & !is.na(our_call), na.rm=T),
              tot_singlets_dem = sum(dem_call == 'singlet', na.rm=T),
              frac_singlet_agree = mean(singlet_agree, na.rm=TRUE),
              frac_doublet_agree = mean(doublet_agree, na.rm=TRUE),
              tot_FPs = sum(FP, na.rm=T),
              tot_FNs = sum(FN, na.rm=T),
              tot_doublets = sum(dem_call == 'doublet', na.rm=T),
              tot_singlets = sum(dem_call == 'singlet', na.rm=T)) %>% 
    dplyr::mutate(FP_rate = 100 * tot_FPs / tot_singlets,
           FN_rate = 100 * tot_FNs / tot_doublets,
           frac_doublet_ours = tot_doublets_ours / n_cells,
           frac_doublets_dem = tot_doublets_dem / n_cells,
           frac_doublets_dem_noLQ = tot_doublets_dem_noLQ / n_cells,
           frac_doublets_dem = ifelse(tot_doublets_dem == 0, NA, frac_doublets_dem),
           frac_doublets_dem_noLQ = ifelse(frac_doublets_dem_noLQ == 0, NA, frac_doublets_dem_noLQ))
  
  #doublet FP vs FN 
  ggplot(expt_stats, aes(FP_rate, FN_rate)) + 
    geom_point(size = 2.5) +
    xlab('False-positive\ndoublets (%)') +
    ylab('False-negative\ndoublets (%)') + 
    theme_Publication()
  ggsave(file.path(fig_dir, 'demux_SNP_doublet_agreement.png'),
         width = 3, height = 2.5)
  
  
  rep_dr <- data_frame(n_cells = seq(2500, 12500, by = 500)) %>% 
    mutate(frac_doublet = n_cells*rep_slope + rep_offset)
  
  ggplot(expt_stats, aes(n_cells, frac_doublet_ours)) + 
    geom_smooth(method = 'lm', color = 'black') +
    geom_point(pch = 21, size = 2.5, fill = 'black', color = 'white', stroke = 0.2) + 
    geom_point(aes(y = frac_doublets_dem), pch = 21, size = 2.5, fill = 'blue', color = 'white', stroke = 0.2) + 
    geom_point(aes(y = frac_doublets_dem_noLQ), pch = 21, size = 2.5, fill = 'green', color = 'white', stroke = 0.2) + 
    theme_Publication() +
    xlab('Number of cells recovered') +
    ylab('Fraction of doublets') +
    geom_line(data = rep_dr, aes(y = frac_doublet), color = 'red') 
  ggsave(file.path(fig_dir, 'demux_SNP_comparison_doubletrate.png'),
         width = 3.5, height = 3)
  
  results <- list(
    expt_stats = expt_stats
  )
  return(results)
}