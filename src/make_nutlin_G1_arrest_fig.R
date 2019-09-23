make_nutlin_G1_arrest_fig <- function() {
  cur_expt <- sc_DE_meta$idasanutlin_24hr_expt1
  
  CL_df <- all_CL_features[[cur_expt$expt_name]]
  
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  
  mod_scores <- read_csv(file.path(out_dir, 'module_scores.csv'))
  mod_scores %<>% full_join(CL_df, by = "CCLE_ID")
  mod_scores %<>% mutate(TP53_status = ifelse(CCLE_ID %in% globals$TP53_WT_cls_pool22, 'TP53 WT', 'TP53 Mut'))
  
  ggplot(mod_scores, aes(AUC_avg, delta_G1_frac, fill = TP53_status)) +
    geom_errorbar(aes(ymax = delta_G1_high, ymin = delta_G1_low), width = 0.01) +
    geom_point(alpha = 0.75, size = 4, pch = 21, color = 'white', stroke = 0.4) +
    ylab('Delta G1 fraction') +
    xlab('Nutlin sensitivity (AUC)') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    cdsr::theme_Publication() +
    guides(fill = guide_legend(title = element_blank()))
  ggsave(file.path(fig_dir, paste0(cur_expt$expt_name, '_', 'delta_G1_scatter.png')), width = 4, height = 3.5)
}
