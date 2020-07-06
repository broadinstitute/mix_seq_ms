make_nutlin_G1_arrest_fig <- function() {
  cur_expt <- sc_DE_meta$idasanutlin_24hr_expt1
  
  CL_df <- all_CL_features[[cur_expt$expt_name]]
  
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  
  mod_scores <- read_csv(file.path(out_dir, 'module_scores.csv'))
  mod_scores %<>% full_join(CL_df, by = "CCLE_ID")
  mod_scores %<>% mutate(TP53_status = ifelse(CCLE_ID %in% globals$TP53_WT_cls_pool22, 'TP53 WT', 'TP53 Mut'))
  
  ggplot(mod_scores %>% filter(Phase == 'G0/G1'), 
         aes(AUC_avg, delta_frac, fill = TP53_status)) +
    geom_errorbar(aes(ymax = delta_high, ymin = delta_low), width = 0.01) +
    geom_point(alpha = 0.75, size = 4, pch = 21, color = 'white', stroke = 0.4) +
    ylab('Delta G0/G1 fraction') +
    xlab('Nutlin sensitivity (AUC)') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    theme_Publication() +
    scale_fill_manual(values = c(`TP53 Mut` = 'black', `TP53 WT` = 'red')) +
    guides(fill = guide_legend(title = element_blank()))
  ggsave(file.path(fig_dir, paste0(cur_expt$expt_name, '_', 'delta_G1_scatter.png')), width = 4, height = 3.5)
  ggsave(file.path(fig_dir, paste0(cur_expt$expt_name, '_', 'delta_G1_scatter.pdf')), width = 4, height = 3.5)
}
