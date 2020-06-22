make_oop_classification_fig <- function() {
  all_ref_lines <- read_csv(here::here('data', 'bulk_reference_CLs.csv'))
  total_panel_size <- nrow(all_ref_lines)
  
  used_expts <- names(sc_expts)[grepl('expt[1|3]$', names(sc_expts))]
  res <- ldply(sc_expts[used_expts], function(cur_expt) {
    print(cur_expt$expt_name)
    base_df <- load_cell_info(cur_expt, QC_filter = FALSE)
    oop_df <- load_out_of_pool_class(cur_expt)
    stopifnot(nrow(base_df) == nrow(oop_df))
    base_df$oop_class <- oop_df$singlet_ID
    base_df %<>% mutate(class_in_pool = singlet_ID == oop_class) %>% 
      filter(cell_quality == 'normal') %>% 
      mutate(expt = cur_expt$expt_name)
  })
  
  # #remove errors due to isogenic pairs
  # isogenic_pairs <- list(c('SW480_LARGE_INTESTINE', 'SW620_LARGE_INTESTINE'),
  #                        c('FTC133_THYROID', 'FTC238_THYROID'),
  #                        # c('HEC1A_ENDOMETRIUM', 'HEC1B_ENDOMETRIUM'),
  #                        c('PATU8988S_PANCREAS', 'PATU8988T_PANCREAS'))
  # to_remove <- c()
  # for (cur_pair in isogenic_pairs) {
  #   cur_iso_inds <- c(which(res$singlet_ID == cur_pair[[1]] & res$oop_class == cur_pair[[2]]),
  #                     which(res$singlet_ID == cur_pair[[2]] & res$oop_class == cur_pair[[1]]))
  #   to_remove %<>% c(cur_iso_inds)  
  # }
  # sprintf('%d/%d isogenic pair errors removed', length(to_remove), nrow(res))
  # res <- res[-to_remove,]
  
  #estimate error rates, apply correction
  st <- ddply(res, .(expt), function(df) {
    cur_e <- prop.test(sum(!df$class_in_pool), nrow(df), conf.level = 0.95)
    data.frame(n_samps = nrow(df),
               err_est = cur_e$estimate,
               err_lb = cur_e$conf.int[[1]],
               err_ub = cur_e$conf.int[[2]],
               pool_size = length(unique(df$singlet_ID)))
  }) %>% mutate(err_est = err_est / (1 - (pool_size-1) / total_panel_size),
                err_lb = err_lb / (1 - (pool_size-1) / total_panel_size),
                err_ub = err_ub / (1 - (pool_size-1) / total_panel_size),
                expt_batch = str_match(expt, '_(expt[0-9])')[,2])
  st %<>% mutate(expt = factor(expt, levels = st %>% 
                                 arrange(expt_batch, expt) %>% 
                                 .[['expt']]))

  st %<>% mutate(pool = ifelse(expt_batch == 'expt3', '99', '24'))
  ggplot(st, aes(expt, err_est, fill = pool)) + 
    geom_bar(stat = 'identity', color = 'black', width = 0.8, lwd = 1) +
    geom_errorbar(aes(color = pool, ymax = err_ub, ymin = err_lb), width = 0.4, lwd = 0.5) +
    # ylim(-0.01, 0.01) +
    theme_Publication() + 
    scale_fill_Publication() + 
    scale_color_Publication() +
    ylab('Error rate') +
    xlab('Experiment') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    guides(color = F,
           fill = guide_legend(title = 'pool size')) 
  ggsave(file.path(fig_dir, 'out_of_pool_err_rate.png'), width = 4, height = 3)

  
  results = list(ip_tab = table(res$class_in_pool),
                 mean_err = mean(st$err_est),
                 sem_err = sd(st$err_est)/sqrt(nrow(st))
  )
  return(results)
}