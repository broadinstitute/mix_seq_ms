
run_CL_downsampling_analysis <- function(cur_expt, recompute = FALSE) {
  
  #PARAMS
  min_LFC <- 1
  min_counts_per_gene <- 5
  min_det_samples <- 0.05
  n_reps <- 100
  min_size <- 5
  size_step <- 5
  examp_n <- c(5, 10, 20, 40)
  
  res_file <- file.path(results_dir, sprintf('CL_downsampling_results_%s.rds', cur_expt$expt_name))
  
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  
  CL_df <- all_CL_features[[cur_expt$expt_name]]
  limma_res <- read_rds(file.path(out_dir, 'limma_res.rds'))
  
  dat <- load_collapsed_profiles(cur_expt, results_dir, type = 'sum')
  all_CLs <- unique(dat$sample_info$CCLE_ID)
  used_genes <- which(rowSums(dat$profile_mat > min_counts_per_gene) > ncol(dat$profile_mat)*min_det_samples)
  
  poss_sizes <- seq(from = min_size, to = length(all_CLs), by = size_step)
  poss_sizes <- poss_sizes[poss_sizes != length(all_CLs)]
 
  get_subsample_results <- function(n_CLs, used_genes = NULL) {
    cur_CLs <- sample_n(data_frame(CL = all_CLs), size = n_CLs, replace = FALSE)$CL
    usamps <- which(dat$sample_info$CCLE_ID %in% cur_CLs)
    cur_limma_res <- fit_viability_models(dat$profile_mat[, usamps,drop=FALSE], dat$sample_info[usamps,], CL_df, used_genes = used_genes)
    return(cur_limma_res)
  }
  
  if (!file.exists(res_file) | recompute) {
    doMC::registerDoMC(cores = 2)
    res <- ldply(poss_sizes, function(cur_size) {
      print(cur_size)
      ldply(seq(n_reps), function(ii) {
        cur_limma_res <- get_subsample_results(cur_size, used_genes = used_genes)
        avg_cor <- cor(cur_limma_res$res_avg$logFC, limma_res$res_avg[match(cur_limma_res$res_avg$Gene, limma_res$res_avg$Gene), 'logFC'], use = 'pairwise.complete.obs')[1]
        int_cor <- cor(cur_limma_res$res_int$logFC, limma_res$res_int[match(cur_limma_res$res_int$Gene, limma_res$res_int$Gene), 'logFC'], use = 'pairwise.complete.obs')[1]
        slope_cor <- cor(cur_limma_res$res_slope$logFC, limma_res$res_slope[match(cur_limma_res$res_slope$Gene, limma_res$res_slope$Gene), 'logFC'], use = 'pairwise.complete.obs')[1]
        
        data.frame(size = cur_size, 
                   it = ii,
                   n_sig_avg = sum(cur_limma_res$res_avg$adj.P.Val < globals$q_thresh & abs(cur_limma_res$res_avg$logFC) > min_LFC),
                   n_sig_int = sum(cur_limma_res$res_int$adj.P.Val < globals$q_thresh & abs(cur_limma_res$res_int$logFC) > min_LFC),
                   n_sig_slope = sum(cur_limma_res$res_slope$adj.P.Val < globals$q_thresh & abs(cur_limma_res$res_slope$logFC) > min_LFC),
                   avg_cor = avg_cor,
                   int_cor = int_cor,
                   slope_cor = slope_cor)
      }, .parallel = TRUE)
    })
      
    #add data point with all 
    res %<>% rbind(
      data.frame(size = length(all_CLs), 
                 it = 1,
                 n_sig_avg = sum(limma_res$res_avg$adj.P.Val < globals$q_thresh & abs(limma_res$res_avg$logFC) > min_LFC),
                 n_sig_int = sum(limma_res$res_int$adj.P.Val < globals$q_thresh & abs(limma_res$res_int$logFC) > min_LFC),
                 n_sig_slope = sum(limma_res$res_slope$adj.P.Val < globals$q_thresh & abs(limma_res$res_slope$logFC) > min_LFC),
                 avg_cor = 1,
                 int_cor = 1,
                 slope_cor = 1))
    write_rds(res, res_file)
  } else {
    res <- read_rds(res_file)
  }
  
  #aggregate stats across reps
  res_avg <- res %>% 
    dplyr::group_by(size) %>% 
    dplyr::summarise(
      avg_avg = mean(n_sig_avg),
      avg_int = mean(n_sig_int),
      avg_slope = mean(n_sig_slope),
      se_avg = sd(n_sig_avg)/sqrt(n_reps),
      se_int = sd(n_sig_int)/sqrt(n_reps),
      se_slope = sd(n_sig_slope)/sqrt(n_reps),
      avg_avgC = mean(avg_cor),
      avg_intC = mean(int_cor),
      avg_slopeC = mean(slope_cor),
      se_avgC = sd(avg_cor)/sqrt(n_reps),
      se_intC = sd(int_cor)/sqrt(n_reps),
      se_slopeC = sd(slope_cor)/sqrt(n_reps)
    )
  
  ggplot(res_avg, aes(size, avg_avgC)) + 
    geom_point() + 
    geom_line() +
    geom_errorbar(data = res_avg %>% filter(!is.na(se_avgC)),
                  aes(ymax = avg_avgC + se_avgC, ymin = avg_avgC - se_avgC)) +
    ggtitle('Avg. response') +
    ylab('Correlation with full profile') +
    xlab('Num. cell lines') +
    ylim(0, 1) +
    cdsr::scale_color_Publication() + 
    cdsr::theme_Publication() +
    geom_vline(xintercept = examp_n, linetype = 'dashed', color = 'red')
  ggsave(file.path(fig_dir, 'trametinib_downsampling_avgC.png'), width = 3.5, height = 3.5)
  
  ggplot(res_avg, aes(size, avg_intC)) + 
    geom_point() + 
    geom_line() +
    geom_errorbar(data = res_avg %>% filter(!is.na(se_intC)),
                  aes(ymax = avg_intC + se_intC, ymin = avg_intC - se_intC)) +
    ggtitle('Viability-independent') +
    ylab('Correlation with full profile') +
    xlab('Num. cell lines') +
    ylim(0, 1) +
    cdsr::scale_color_Publication() + 
    cdsr::theme_Publication() +
    geom_vline(xintercept = examp_n, linetype = 'dashed', color = 'red')
  ggsave(file.path(fig_dir, 'trametinib_downsampling_intC.png'), width = 3.5, height = 3.5)
  
  
  ggplot(res_avg, aes(size, avg_slopeC)) + 
    geom_point() + 
    geom_line() +
    geom_errorbar(data = res_avg %>% filter(!is.na(se_slopeC)),
                  aes(ymax = avg_slopeC + se_slopeC, ymin = avg_slopeC - se_slopeC)) +
    ggtitle('Viability-related') +
    ylab('Correlation with full profile') +
    xlab('Num. cell lines')+
    ylim(0, 1) +
    cdsr::scale_color_Publication() + 
    cdsr::theme_Publication() +
    geom_vline(xintercept = examp_n, linetype = 'dashed', color = 'red')
  ggsave(file.path(fig_dir, 'trametinib_downsampling_slopeC.png'), width = 3.5, height = 3.5)
  

  make_subsample_scatter <- function(cur_N) {
    cur_limma_res <- get_subsample_results(cur_N, used_genes = used_genes)
    comb <- full_join(cur_limma_res$res_avg, limma_res$res_avg, by = 'Gene', suffix = c('_sub', '_full'))
    g1 <- ggplot(comb, aes(logFC_full, logFC_sub)) + 
      geom_point(pch = 21, size = 1.5, fill = 'black', color = 'white', stroke = 0.1) + 
      geom_abline() +
      xlab('logFC (all cell lines)') +
      ylab('logFC (subset)') +
      geom_hline(yintercept = 0, linetype = 'dashed') + 
      geom_vline(xintercept = 0, linetype = 'dashed') +
      cdsr::theme_Publication() +
      ggtitle('Average')
    comb <- full_join(cur_limma_res$res_slope, limma_res$res_slope, by = 'Gene', suffix = c('_sub', '_full'))
    g2 <- ggplot(comb, aes(logFC_full, logFC_sub)) + 
      geom_point(pch = 21, size = 1.5, fill = 'black', color = 'white', stroke = 0.1) + 
      geom_abline() +
      xlab('logFC (all cell lines)') +
      ylab('logFC (subset)') +
      geom_hline(yintercept = 0, linetype = 'dashed') + 
      geom_vline(xintercept = 0, linetype = 'dashed') +
      cdsr::theme_Publication() +
      ggtitle('Viability-related')
    cowplot::plot_grid(g1, g2)
  }
  
  set.seed(1)
  for (n in examp_n) {
    make_subsample_scatter(n)
    ggsave(file.path(fig_dir, sprintf('trametinib_CL_subsample_%d.png', n)), width = 6, height = 3)
  }
}


