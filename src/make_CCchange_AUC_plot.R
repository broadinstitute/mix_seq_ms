make_CCchange_AUC_plot <- function() {
  sens_max <- 1; sens_min <- 0 #min and max values of sensitivity for estimating weighted averages
  
  used_expts <- c(
    'Idasanutlin_24hr_expt1',
    # 'Idasanutlin_6hr_expt1',
    # 'Trametinib_6hr_expt1',
    'Trametinib_24hr_expt3', 
    'Bortezomib_24hr_expt1',
    'Navitoclax_24hr_expt3',
    'Dabrafenib_24hr_expt3',
    'BRD3379_24hr_expt3',
    # 'BRD3379_6hr_expt3',
    'Afatinib_expt10',
    'AZD5591_expt10',
    'Everolimus_expt10',
    'Taselisib_expt10',
    'Gemcitabine_expt10',
    'Prexasertib_expt10',
    'JQ1_expt10')
  
  expt_list <- sc_DE_meta
  names(expt_list) <- lapply(expt_list, function(x) x$expt_name)
  all_cc_features <- llply(expt_list, function(cur_expt) {
    print(sprintf('Processing expt %s', cur_expt$expt_name))
    out_dir <- file.path(results_dir, cur_expt$expt_name)
  
    #load cell cycle
    mod_stats <- read_csv(file.path(out_dir, 'module_scores.csv'))
    
    CL_df <- mod_stats %>%
      dplyr::mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(CCLE_ID)) %>% 
      dplyr::select(DEPMAP_ID, delta_frac, Phase) %>% 
      left_join(all_CL_features[[cur_expt$expt_name]], by = 'DEPMAP_ID') 
    
     return(CL_df)
  })
  
  aa <- all_cc_features %>% 
    ldply(function(df) {
      #use bounded sensitivity for weighted averages
      df %<>% dplyr::mutate(weights = pmin(sens, sens_max),
                     weights = pmax(weights, sens_min)) %>% 
        tidyr::spread(Phase, delta_frac)
      G2M_c <- cor.test(df$`G2/M`, df$sens)
      G1_c <- cor.test(df$`G0/G1`, df$sens)
      S_c <- cor.test(df$S, df$sens)
      wtd_t_G1 <- weights::wtd.t.test(df$`G0/G1`, weight = df$weights)
      wtd_t_G2M <- weights::wtd.t.test(df$`G2/M`, weight = df$weights)
      wtd_t_S <- weights::wtd.t.test(df$`S`, weight = df$weights)
      data_frame(wmean_G2M = weighted.mean(df$`G2/M`[!is.na(df$weights)], w = df$weights[!is.na(df$weights)], na.rm=T),
                 wmean_S = weighted.mean(df$S[!is.na(df$weights)], w = df$weights[!is.na(df$weights)], na.rm=T),
                 wmean_G1 = weighted.mean(df$`G0/G1`[!is.na(df$weights)], w = df$weights[!is.na(df$weights)], na.rm=T),
                 G1_t_pval = wtd_t_G1$coefficients['p.value'],
                 G2M_t_pval = wtd_t_G2M$coefficients['p.value'],
                 S_t_pval = wtd_t_S$coefficients['p.value'],
                 G2M_cor = G2M_c$estimate,
                 G2M_lci = G2M_c$conf.int[1],
                 G2M_uci = G2M_c$conf.int[2],
                 G1_cor = G1_c$estimate,
                 G1_lci = G1_c$conf.int[1],
                 G1_uci = G1_c$conf.int[2],
                 S_cor = S_c$estimate,
                 S_lci = S_c$conf.int[1],
                 S_uci = S_c$conf.int[2])
    }, .id = 'expt') %>% 
    dplyr::mutate(drug = str_match(expt, '(^[:alnum:]+)_')[,2])
  
  #PLOT WEIGHTED AVERAGES OF CC PHASE CHANGES
  ord <- aa %>% 
    dplyr::filter(expt %in% used_expts) %>% 
    dplyr::mutate(tot = abs(wmean_S) + abs(wmean_G1) + abs(wmean_G2M)) %>% 
    dplyr::arrange(tot) %>% 
    .[['drug']]
  aa %>% 
    dplyr::filter(expt %in% used_expts) %>% 
    dplyr::mutate(drug = factor(drug, levels = ord)) %>% 
    dplyr::select(drug, wmean_S, wmean_G1, wmean_G2M) %>%
    reshape2::melt(id.vars = 'drug') %>% 
    dplyr::mutate(variable = str_replace(variable, 'wmean_', ''),
                  variable = plyr::revalue(variable, replace = c(G1 = 'G0/G1', G2M = 'G2/M')),
           variable = factor(variable, levels = c('G0/G1', 'S', 'G2/M'))) %>% 
    ggplot(aes(drug, value, fill = variable, group = variable)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    theme_Publication() +
    ylab('Delta cell fraction') +
    guides(fill = guide_legend(title = 'Phase')) +
    scale_fill_Publication() +
    theme(axis.text = element_text(angle = 90, hjust = 1)) 
  ggsave(file.path(fig_dir, 'CC_wphase_all_compounds.png'), width = 4.5, height = 3.5)
  
 
  #MAKE SCATTERPLOT OF DELTA-PHASE SENS CORRELATIONS
  aa %<>% 
    dplyr::filter(expt %in% used_expts) %>% 
    dplyr::mutate(expt = str_replace_all(expt, 'expt10', '24hr_expt10')) %>% 
    dplyr::filter(!(expt %in% c('Trametinib_24hr_expt1'))) %>% #pick the trametinib 24hr expt with more cell lines
    dplyr::mutate(expt = str_replace_all(expt, '_expt[0-9]+', ''),
           drug = str_match(expt, '^([:alnum:]+)_')[,2])
  
  ggplot(aa, aes(S_cor, G2M_cor)) + 
    ggrepel::geom_label_repel(aes(label = drug, color = drug), size = 2.5) +
    geom_errorbar(aes(ymax = G2M_uci, ymin = G2M_lci), alpha = 0.5) +
    geom_errorbarh(aes(xmax = S_uci, xmin = S_lci), alpha = 0.5) +
    geom_point(aes(fill = drug), pch = 21, size = 3, color = 'gray', stroke = 0.2) +
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    geom_vline(xintercept = 0, linetype = 'dashed') +
    guides(color = F, fill = F) + 
    theme_Publication() +
    xlab('S-phase sensitivity corr.') +
    ylab('G2/M-phase sensitivity corr.')
  ggsave(file.path(fig_dir, 'G2M_S_sens_corrs.png'), width = 4, height = 3.5)

  write_rds(aa, file.path(results_dir, 'CC_arrest_stats.rds'))

  
  library(pals)
  library(plot3D)
  aa %<>% mutate(drug = factor(drug))
  color = pals::alphabet(nlevels(aa$drug)) %>% set_names(levels(aa$drug))
  aa %<>% mutate(col = color[drug])
  x <- aa$S_cor
  y <- aa$G2M_cor
  z <- aa$G1_cor
  # Compute the linear regression (z = ax + by + d)
  fit <- lm(z ~ x + y)
  # predict values on regular xy grid
  grid.lines = 8
  x.pred <- seq(min(x), max(x), length.out = grid.lines)
  y.pred <- seq(min(y), max(y), length.out = grid.lines)
  xy <- expand.grid( x = x.pred, y = y.pred)
  z.pred <- matrix(predict(fit, newdata = xy), 
                   nrow = grid.lines, ncol = grid.lines)
  # fitted points for droplines to surface
  fitpoints <- predict(fit)
  # dev.off()
  # png(file.path(fig_dir, 'G2M_S_G1_sens_all_corrs.png'), width = 4.5, height = 4.5, units = 'in', res = 300)
  pdf(file.path(fig_dir, 'G2M_S_G1_sens_all_corrs.pdf'), width = 4.5, height = 4.5)
  # scatter plot with regression plane
  scatter3D(x, y, z, pch = 16, cex = 1.5, col = aa$col, colvar = NULL,
            theta = 110, phi = 20, ticktype = "detailed",
            xlab = "S cor.", ylab = "G2/M cor.", zlab = "G0/G1 cor.",
            surf = list(x = x.pred, y = y.pred, z = z.pred,
                        facets = NA, fit = fitpoints, lwd = 0.5, col = 'darkgray'),
            bty = 'b2')
  text3D(x+0.1,y+0.075,z, pch=10,labels=aa$drug, col = aa$col, add=TRUE, cex=0.75)
  dev.off()
  

}
