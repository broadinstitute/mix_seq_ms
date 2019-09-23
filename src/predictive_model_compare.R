
# Takes a list of predictor sets and cbinds to a single predictor matrix
proc_list_data <- function(X_list, used_CLs) {
  X <- lapply(names(X_list), function(x_name) {
    x <- X_list[[x_name]]
    colnames(x$data) <- paste0(colnames(x$data), '_', x_name)
    return(x$data[used_CLs, ])
  }) %>% do.call(cbind, .)
  return(X)
}

#find CL list intersecting all rownames of X and target y
get_used_CLs <- function(X, y) {
  extract_CLs <- function(d) {
    if (class(d) == 'list') {
      if (length(d) > 1) {
        used_CLs <- Reduce(intersect, lapply(d, function(x) {rownames(x$data %>% na.omit)}))
      } else {
        used_CLs <- rownames(d[[1]]$data %>% na.omit)
      }
    } else if (class(d) == 'matrix') {
      used_CLs <- rownames(d)
    } else {
      used_CLs <- names(d)
    }
    return(used_CLs)
  }
  data_sets <- list(X = X, y = y)
  return(Reduce(intersect, lapply(data_sets, extract_CLs)))
}


#' Take a set of inputs which may be in list or matrix form, and process into model-function inputs
#'
#' @param X
#' @param y
#' @param Z
#' @param family
#'
#' @return
#' @export
#'
#' @examples
proc_input_data <- function(X, y) {
  used_CLs <- get_used_CLs(X, y[!is.na(y)])
  if (class(X) == 'list') {
    X <- proc_list_data(X, used_CLs)
  } else if (class(X) == 'matrix') {
    X <- X[used_CLs, , drop=FALSE]
  }
  y <- y[used_CLs]

  X %<>% set_colnames(make.names(colnames(X)))
  return(list(X = X, y = y))
}

#get top marginal correlations between X and y in input data
get_top_feat <- function(input_data, n_top_feat, sample_set = NULL) {
  if (!is.null(sample_set)) {
    X = input_data$X[sample_set,]
    y = input_data$y[sample_set]
  } else {
    X <- input_data$X
    y <- input_data$y
  }
  cor(X, as.matrix(y, ncol = 1), use = 'pairwise.complete.obs') %>% 
    set_colnames('cor') %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'Gene') %>% 
    dplyr::filter(!is.na(cor)) %>% 
    dplyr::arrange(dplyr::desc(abs(cor))) %>% 
    head(n_top_feat) %>% 
    .[['Gene']]
}

# Get random forest out-of-sample predictions using k-fold CV
get_rf_xv_results <- function(input_data, kfold, fids, n_top_feat) {
  LFC_results <- ldply(seq(kfold), function(xv) {
    cur_train <- which(fids != xv)
    cur_test <- which(fids == xv)
    used_feat <- get_top_feat(input_data, n_top_feat = n_top_feat, sample_set = cur_train)
    rf_df <- cbind(input_data$X[, used_feat], data_frame(response = input_data$y))
    mod_fit <- ranger::ranger('response ~ .', data = rf_df[cur_train,])
    test_preds <- predict(mod_fit, data = rf_df[cur_test,])$predictions
    data_frame(
      xv_id = xv,
      DEPMAP_ID = names(cur_test),
      prediction = test_preds,
      obs = rf_df[cur_test, 'response']
    )
  }, .parallel = TRUE)
}


run_predictive_model_compare <- function() {
  #PARAMS
  min_cells <- 5 #minimum cells per condition to include cell line in analysis
  min_R2 <- -1
  n_top_feat <- 1000
  random_seed <- 1
  kfold <- 10
  model_type <- 'regression'
  models <- c('rf')
 
  feat_dat <- list(
    GE = list(
      data.name="depmap-rnaseq-expression-data-ccd0",
      data.version = 14,
      data.file = 'CCLE_depMap_19Q3_TPM',
      transpose = F),
    # CN = list(data.name = 'depmap-wes-cn-data-97cc', 
    #           data.version = 18, 
    #           data.file = 'public_19Q3_gene_cn',
    #           transpose = F),
    MUT_DAM = list(
      data.name = 'depmap-mutation-calls-9a1a',
      data.version = 12,
      data.file = 'damaging_mutation',
      transpose = F),
    MUT_HOT = list(
      data.name = 'depmap-mutation-calls-9a1a',
      data.version = 12,
      data.file = 'hotspot_mutation',
      transpose = F)) %>% 
    taigr::load.all.from.taiga() 
  
  doMC::registerDoMC(cores = 3)
  
  targ_dsets <- sc_DE_meta[str_match(names(sc_DE_meta), '(expt[0-9]+)$')[,2] %in% c('expt1', 'expt3', 'expt10')]
  all_pmod_results <- ldply(targ_dsets, function(cur_expt) {
    print(sprintf('Processing expt %s', cur_expt$expt_name))
    out_dir <- file.path(results_dir, cur_expt$expt_name)
  
    cell_df <- load_cell_info(cur_expt, QC_filter = TRUE) %>% 
      tidyr::separate(col = 'condition', into = c('condition', 'batch'), sep = '_')
    cell_counts <- cell_df %>% 
      dplyr::group_by(singlet_ID, condition) %>% 
      dplyr::summarise(n = n()) %>% 
      reshape2::dcast(singlet_ID ~ condition, value.var = 'n')
    usable_CLs <- cell_counts %>% 
      dplyr::filter(control >= min_cells, treat >= min_cells) %>% 
      .[['singlet_ID']] %>% 
      celllinemapr::ccle.to.arxspan()
    
    CL_df <- all_CL_features[[cur_expt$expt_name]]

    dat <- load_compute_CPM_stats(cur_expt, results_dir, prior_counts = globals$pca_prior_cnt, type = 'avg_cpm')
    colnames(dat$LFC_mat) <- celllinemapr::ccle.to.arxspan(colnames(dat$LFC_mat))
    
    y <- CL_df$AUC_avg %>% set_names(CL_df$DEPMAP_ID)
    y <- y[!is.na(y)]

    #MODEL BASED PREDICTION USING LFC  
    X <- list(
      LFC = list(data = t(dat$LFC_mat))
    )
    input_data <- proc_input_data(X, y)
    fids <- caret::createFolds(input_data$y, k = kfold, list = FALSE) %>% 
      set_names(names(input_data$y))
    n_LFC_CLs <- length(input_data$y)
    LFC_CLs <- names(input_data$y)
    LFC_results <- get_rf_xv_results(input_data, kfold, fids = fids, n_top_feat = n_top_feat)

    
    #NOW TRY MODELS USING BASELINE FEATURES
    X <- list(
      GE = list(data = feat_dat$GE),
      MUT_DAM = list(data = feat_dat$MUT_DAM),
      MUT_HOT = list(data = feat_dat$MUT_HOT)
      # CN = list(data = feat_dat$CN)
    )
    input_data <- proc_input_data(X, y[usable_CLs])
    baseline_matched_CLs <- names(input_data$y)
    baseline_matched_results <- get_rf_xv_results(input_data, kfold, fids = fids[names(input_data$y)], n_top_feat = n_top_feat)
      
    #USING ALL AVAILABLE CELL LINES
   input_data <- proc_input_data(X, y)
   fids <- caret::createFolds(input_data$y, k = kfold, list = FALSE) %>% 
     set_names(names(input_data$y))
   baseline_all_CLs <- names(input_data$y)
   baseline_all_results <- get_rf_xv_results(input_data, kfold, fids = fids, n_top_feat = n_top_feat)
      
      pred_mod_errs <- rbind(
        data_frame(type = 'matched_sample_pred',
                   n = length(baseline_matched_CLs),
                   pear_cor = cor.test(baseline_matched_results$prediction, baseline_matched_results$obs)$estimate,
                   pear_cor_uci = cor.test(baseline_matched_results$prediction, baseline_matched_results$obs)$conf.int[2],
                   pear_cor_lci = cor.test(baseline_matched_results$prediction, baseline_matched_results$obs)$conf.int[1]),
        data_frame(type = 'all_sample_pred',
                   n = length(baseline_all_CLs),
                   pear_cor = cor.test(baseline_all_results$prediction, baseline_all_results$obs)$estimate,
                   pear_cor_uci = cor.test(baseline_all_results$prediction, baseline_all_results$obs)$conf.int[2],
                   pear_cor_lci = cor.test(baseline_all_results$prediction, baseline_all_results$obs)$conf.int[1]),
        data_frame(type = 'LFC_pred',
                   n = length(LFC_CLs),
                   pear_cor = cor.test(LFC_results$prediction, LFC_results$obs)$estimate,
                   pear_cor_uci = cor.test(LFC_results$prediction, LFC_results$obs)$conf.int[2],
                   pear_cor_lci = cor.test(LFC_results$prediction, LFC_results$obs)$conf.int[1])
      ) %>% mutate(expt = cur_expt$expt_name)
      return(pred_mod_errs)
  })
  
  write_rds(all_pmod_results, file.path(results_dir, 'all_pmod_results.rds'))
} 

make_predictive_model_figs <- function() {
  all_pmod_results <- read_rds(file.path(results_dir, 'all_pmod_results.rds'))
  all_pmod_results %<>% mutate(
    expt = str_replace_all(expt, '_expt10', '_24hr_expt10')
  )
  used_expts <- c('Trametinib_6hr_expt1',
                  'Bortezomib_6hr_expt1',
                  'Bortezomib_24hr_expt1',
                  'Idasanutlin_6hr_expt1',
                  'Idasanutlin_24hr_expt1',
                  'Trametinib_24hr_expt3', 
                  'Dabrafenib_24hr_expt3',
                  # 'CGS15943_6hr_expt3',
                  'BRD3379_6hr_expt3',
                  'BRD3379_24hr_expt3',
                  'Afatinib_24hr_expt10',
                  'AZD5591_24hr_expt10',
                  'Everolimus_24hr_expt10',
                  'Gemcitabine_24hr_expt10',
                  'Prexasertib_24hr_expt10',
                  'Taselisib_24hr_expt10',
                  'JQ1_24hr_expt10')
  
  
  ov_R2_vals <- all_pmod_results %>% 
    filter(expt %in% used_expts) %>% 
    distinct(expt, type, .keep_all=TRUE) %>% 
    dplyr::select(type, expt, ov_R2) %>% 
    reshape2::dcast(expt ~ type, value.var = 'ov_R2') %>% 
    mutate(time = str_match(expt, '^[:alnum:]+_([:alnum:]+)')[,2],
           drug = str_match(expt, '^([:alnum:]+)_')[,2]) 
  
  #scatterplot both times, all samples
  ov_R2_vals %>% 
    ggplot(aes(all_sample_pred, LFC_pred)) + 
    geom_point(aes(fill = time), pch = 21, size = 3, color = 'white', stroke = 0.2) +
    geom_abline() +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggrepel::geom_label_repel(aes(label = drug, color = time), size = 3, label.padding = 0.1) +
    coord_cartesian(xlim = c(-0.1, 0.6), ylim = c(-0.1, 0.6)) +
    xlab('Baseline features\n(R2)') +
    ylab('Transcriptional responses\n(R2)') +
    guides(color = F) + 
    cdsr::theme_Publication() +
    cdsr::scale_color_Publication() +
    cdsr::scale_fill_Publication()
  ggsave(file.path(fig_dir, 'pred_R2_scatter_all.png'), width = 4, height = 4)
  
  
  n_vals <- all_pmod_results %>% 
    filter(expt %in% used_expts) %>% 
    distinct(expt, type, .keep_all=T) %>% 
    dplyr::select(type, expt, n) %>% 
    reshape2::dcast(expt ~ type, value.var = 'n')%>% 
    mutate(time = str_match(expt, '^[:alnum:]+_([:alnum:]+)')[,2],
           drug = str_match(expt, '^([:alnum:]+)_')[,2]) %>% 
    left_join(ov_R2_vals %>% 
                mutate(LFC_baseline_diff = LFC_pred - all_sample_pred) %>% 
                dplyr::select(expt, LFC_baseline_diff), 
              by = 'expt')
  ggplot(n_vals, aes(LFC_pred, LFC_baseline_diff)) + 
    geom_point(aes(fill = time), pch = 21, size = 3, color = 'white', stroke = 0.2) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggrepel::geom_label_repel(aes(label = drug, color = time), size = 3, label.padding = 0.1) +
    xlab('LFC sample size') +
    ylab('Prediction accuracy diff.\n(LFC_R2 - Baseline_R2)') +
    guides(color = F) + 
    cdsr::theme_Publication() +
    cdsr::scale_color_Publication() +
    cdsr::scale_fill_Publication()
  ggsave(file.path(fig_dir, 'pred_R2_diff_vs_n.png'), width = 4, height = 4)
  
  
  #scatterplot boht times, matched samples
  ov_R2_vals %>% 
    ggplot(aes(matched_sample_pred, LFC_pred)) + 
    geom_point(aes(fill = time), pch = 21, size = 3, color = 'white', stroke = 0.2) +
    geom_abline() +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggrepel::geom_label_repel(data = ov_R2_vals,
                              aes(label = drug, color = time), size = 3, label.padding = 0.1) +
    coord_cartesian(xlim = c(-0.1, 0.6), ylim = c(-0.1, 0.6)) +
    xlab('Baseline features\n(R2)') +
    ylab('Transcriptional responses\n(R2)') +
    guides(color = F) +
    cdsr::theme_Publication() +
    cdsr::scale_color_Publication() +
    cdsr::scale_fill_Publication()
  ggsave(file.path(fig_dir, 'pred_R2_scatter_matched.png'), width = 4, height = 4)
  
  #matched samples, 24hr only
  ov_R2_vals %>% 
    filter(time == '24hr') %>% 
    ggplot(aes(matched_sample_pred, LFC_pred)) + 
    geom_point(aes(fill = time), pch = 21, size = 3, color = 'white', stroke = 0.2) +
    geom_abline() +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ggrepel::geom_label_repel(aes(label = drug, color = time), size = 2.5, label.padding = 0.1) +
    coord_cartesian(xlim = c(-0.1, 0.6), ylim = c(-0.1, 0.6)) +
    xlab('Baseline features\n(R2)') +
    ylab('Transcriptional responses\n(R2)') +
    cdsr::theme_Publication() +
    guides(color = F) + 
    cdsr::scale_color_Publication() +
    cdsr::scale_fill_Publication() 
  ggsave(file.path(fig_dir, 'pred_R2_scatter_matched_24hr.png'), width = 4, height = 4)
  
  
  ov_R2_vals %>% 
    dplyr::filter(grepl('6hr', expt)) %>% 
    ggplot(aes(all_sample_pred, LFC_pred)) + 
    geom_point() +
    geom_abline() +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed')

}