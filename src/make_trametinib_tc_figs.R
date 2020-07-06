make_trametinib_tc_figs <- function() {
  library(scales)
  
  #PARAMS
  n_top_sig_genes <- 20
  n_PCs <- 100
  min_cells_per_cond <- 5
  min_counts_per_gene <- 1
  min_det_samps_per_gene <- 5
  n_top_genes = 50 #number of top genes per category to avg for making avg relative time course plot
  
  #umap params
  umap_n_neighbors <- 15
  umap_min_dist <- 0.6
  
  #for heatmap plots
  n_int_genes <- 10 
  n_slope_genes <- 10
  
  #for time-course example plots
  gene_1 <- 'EGR1'
  gene_2 <- 'MCM7'
  CL1 <- 'RCM1_LARGE_INTESTINE'
  CL2 <- 'TEN_ENDOMETRIUM'
  
  tp_levs <- c('Tram_0', 'Tram_3', 'Tram_6', 'Tram_12', 'Tram_24', 'Tram_48')
  
  cur_expt <- sc_expts$tramet_tc_expt5
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  
  
  ##MAKE EXAMPLE TIME COURSE PLOTS
  make_gene_vlnplot <- function(object, cur_CL, cur_gene) {
    cond_col_pal <- RColorBrewer::brewer.pal(name = 'Blues', n = 6) %>% set_names(c('0hr', '3hr', '6hr', '12hr', '24hr', '48hr')) %>% c(DMSO = 'pink')
    ids <- Seurat::FetchData(object, 'singlet_ID')
    object <- object[, ids$singlet_ID == cur_CL]
    df <- Seurat::FetchData(object, cur_gene) %>% 
      cbind(object@meta.data[, setdiff(colnames(object@meta.data), cur_gene)]) %>% 
      dplyr::mutate(hash_tag = as.character(hash_tag),
             hash_tag = ifelse(grepl('DMSO', hash_tag), 'DMSO', hash_tag),
             hash_tag = plyr::revalue(hash_tag, c(`Untreated_48hr` = 'Tram_0hr')),
             hash_tag = str_replace_all(hash_tag, 'Tram_', ''),
             hash_tag = factor(hash_tag, levels = c('DMSO', '0hr', '3hr', '6hr', '12hr', '24hr', '48hr'))) %>% 
      dplyr::filter(!is.na(hash_tag))  
    
    #add slight jitter to expression levels (taken from Seurat function)
    df[, cur_gene] <- df[, cur_gene] + rnorm(n = nrow(df))/1e+05
    
    df$GE = df[[cur_gene]]
    #get average and SEM per group
    avgs <- df %>% 
      dplyr::group_by(hash_tag) %>% 
      dplyr::summarise(SD = stats::sd(GE),
                       GE = mean(GE, na.rm=T),
                       n = n()) %>% 
      dplyr::mutate(SEM = SD/sqrt(n))
    ggplot(df, aes(hash_tag, GE)) +
      geom_violin(aes(fill = hash_tag), alpha = 0.5, scale = 'width') +
      geom_jitter(height = 0, size = 0.75, width = 0.2, alpha = 0.75) +
      ylab(paste0(cur_gene, ' (logCPM)')) +
      geom_point(data = avgs, size = 3, alpha = 0.75, color = 'red') + 
      geom_errorbar(data = avgs, aes(ymax = GE + SEM, ymin = GE - SEM), alpha = 0.75, color = 'red', width = 0.2) + 
      geom_line(data = avgs, aes(group = 1)) +
      theme_Publication() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
            axis.title.x = element_blank()) +
      scale_fill_manual(values = cond_col_pal) +
      guides(fill = F) 
  }
  
  
  CL_df <- all_CL_features[[sc_DE_meta$trametinib_24hr_expt1$expt_name]] #CL features are same for expt1 trametinib data
  
  seuObj <- load_sc_data(cur_expt, sc_expts = sc_expts)
  
  #get rid of any unlabeled cells or hash tag multiplets
  hash_assign <- FetchData(seuObj, 'hash_tag')
  seuObj <- seuObj[, !(hash_assign$hash_tag %in% c('unknown', 'multiplet'))]
  
  #set DMSO and untreated to control condition for normalization purposes
  df <- seuObj@meta.data %>% 
    dplyr::mutate(condition = ifelse(grepl("DMSO", hash_tag), 'control', hash_tag),
           condition = ifelse(hash_tag == 'Untreated_48hr', 'control', condition))
  seuObj[['condition']] <- df$condition
  
  seuObj <- NormalizeData(object = seuObj, 
                          normalization.method = "LogNormalize", 
                          scale.factor = globals$scale_fac)
  
  seuObj <- CellCycleScoring(object = seuObj,
                             s.features = Seurat::cc.genes$s.genes,
                             g2m.features = Seurat::cc.genes$g2m.genes,
                             set.ident = FALSE)
  
  n_tot_cells <- nrow(seuObj@meta.data)
  print(sprintf('%s total cells', n_tot_cells))
  
  ## MAKE embedding FIGURE
  seuObj <- ScaleData(object = seuObj)
  seuObj <- FindVariableFeatures(object = seuObj,
                                 nfeatures = globals$n_highvar_genes,
                                 selection.method = 'vst')
  seuObj <- RunPCA(object = seuObj,
                   features = VariableFeatures(seuObj),
                   seed.use = 1,
                   npcs = n_PCs,
                   verbose = FALSE)
  
  #make TSNE figure
  # seuObj <- RunTSNE(object = seuObj, 
  #                   dims = 1:n_PCs, 
  #                   check_duplicates = FALSE, 
  #                   seed.use = 1, 
  #                   perplexity = globals$tsne_perplexity)
  seuObj <- RunUMAP(object = seuObj, 
                    dims = 1:n_PCs, 
                    seed.use = 1, 
                    n.neighbors = umap_n_neighbors, 
                    min.dist = umap_min_dist)
  
  df <- Embeddings(object = seuObj, reduction = 'umap') %>% 
    as.data.frame() %>% 
    set_colnames(c('t1', 't2')) %>% 
    rownames_to_column(var = 'barcode') %>% 
    cbind(seuObj@meta.data) %>% 
    dplyr::mutate(hash_tag = as.character(hash_tag),
                  hash_tag = ifelse(grepl("DMSO", hash_tag), 'control', hash_tag),
           hash_tag = ifelse(hash_tag == 'Untreated_48hr', 'control', hash_tag),
           time_point = str_match(hash_tag, 'Tram_([0-9]+)hr')[,2],
           time_point = factor(time_point, levels = c(3, 6, 12, 24, 48)))
  
  #get avgs per cell line
  avgs <-  df %>% 
    dplyr::group_by(singlet_ID) %>% 
    dplyr::summarise(ct1 = median(t1, na.rm=T),
                     ct2 = median(t2, na.rm=T))
  
  y_offset = 0.2 #add slight shift to cell line name y-locs
  avgs %<>% dplyr::mutate(short_name = str_match(singlet_ID, '^([:alnum:]+)_')[,2],
                   ct2 = ct2 + y_offset)
  ggplot(df, aes(t1, t2)) + 
    geom_point(data = df %>% filter(hash_tag == 'control'), fill = 'pink', pch = 21, size = 0.5, alpha = 0.75, stroke = 0.2, color = 'gray') +
    geom_point(data = df %>% filter(hash_tag != 'unknown', !is.na(time_point)), aes(fill = time_point), pch = 21, size = 0.75, color = 'gray', alpha = 0.8, stroke = 0.1) +
    scale_fill_brewer(palette = 1, type = 'seq') +
    geom_text(data = avgs, mapping = aes(x = ct1, y = ct2, label = short_name), fontface = 'bold', color = 'black', size = 2.) + 
    theme_Publication() + 
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    guides(fill = guide_legend(override.aes = list(size = 3)))
  ggsave(file.path(fig_dir, 'trametinib_tc_fulltsne.png'),
         width = 4.5, height = 4)
  ggsave(file.path(fig_dir, 'trametinib_tc_fulltsne.pdf'),
         width = 4.5, height = 4)  
  
  ## MAKE EXAMPLE TIME COURSE PLOTS
  #set order of hash_tag factor levels
  seuObj@meta.data$hash_tag %<>% factor(
    levels = c('Untreated_48hr',
               'DMSO_3hr', 
               'DMSO_6hr',
               'DMSO_12hr',
               'DMSO_24hr',
               'DMSO_48hr',
               'Tram_3hr',
               'Tram_6hr',
               'Tram_12hr',
               'Tram_24hr',
               'Tram_48hr')
  )
  g1 <- make_gene_vlnplot(seuObj, CL1, gene_1)
  g2 <- make_gene_vlnplot(seuObj, CL1, gene_2)
  ggsave(file = file.path(fig_dir, sprintf('tc_violin_%s_%s.png', CL1, gene_1)),
         plot = g1, width = 4., height = 3.5)
  ggsave(file = file.path(fig_dir, sprintf('tc_violin_%s_%s.pdf', CL1, gene_1)),
         plot = g1, width = 4., height = 3.5)
  ggsave(file = file.path(fig_dir, sprintf('tc_violin_%s_%s.png', CL1, gene_2)),
         plot = g2, width = 4., height = 3.5)
  ggsave(file = file.path(fig_dir, sprintf('tc_violin_%s_%s.pdf', CL1, gene_2)),
         plot = g2, width = 4., height = 3.5)  
  
  g1 <- make_gene_vlnplot(seuObj, CL2, gene_1)
  g2 <- make_gene_vlnplot(seuObj, CL2, gene_2)
  ggsave(file = file.path(fig_dir, sprintf('tc_violin_%s_%s.png', CL2, gene_1)),
         plot = g1, width = 4., height = 3.5)
  ggsave(file = file.path(fig_dir, sprintf('tc_violin_%s_%s.pdf', CL2, gene_1)),
         plot = g1, width = 4., height = 3.5)
  ggsave(file = file.path(fig_dir, sprintf('tc_violin_%s_%s.png', CL2, gene_2)),
         plot = g2, width = 4., height = 3.5)
  ggsave(file = file.path(fig_dir, sprintf('tc_violin_%s_%s.pdf', CL2, gene_2)),
         plot = g2, width = 4., height = 3.5)  
  
  ## GET SUM-COLLAPSED PROFILES
  #group by cell line/condition
  cell_df <- seuObj@meta.data %>% 
    rownames_to_column(var = 'cell_ID') %>% 
    dplyr::filter(!is.na(hash_tag)) %>% 
    dplyr::mutate(ID = paste0(singlet_ID, '.', hash_tag))
  
  #identify usable cell-line/conditions
  all_groups <- cell_df %>% 
    dplyr::group_by(ID) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::filter(n >= min_cells_per_cond) %>% 
    .[['ID']]
  
  summed_counts <- plyr::laply(all_groups, function(cur_group) {
      cur_cells <- cell_df %>% 
        dplyr::filter(ID == cur_group) %>% 
        .[['cell_ID']]
      Matrix::rowSums(GetAssayData(seuObj, 'counts')[, cur_cells, drop = FALSE])
    }) %>% set_rownames(all_groups)
  summed_counts %<>% t()
    
  seuObj_anorm <- Seurat::NormalizeData(object = seuObj, 
                                  normalization.method = "LogNormalize", 
                                  scale.factor = 1e6) #CPM normalization
  cond_avgs <- plyr::laply(all_groups, function(cur_group) {
    cur_cells <- cell_df %>% 
      dplyr::filter(ID == cur_group) %>% 
      .[['cell_ID']]
    log2(globals$prior_cnt + rowMeans(expm1(GetAssayData(seuObj_anorm, slot = 'data')[,cur_cells])))
  }) %>% set_rownames(all_groups)
  cond_avgs %<>% t()
  
  
  sample_info <- data_frame(sample_name = colnames(summed_counts)) %>% 
    tidyr::separate(sample_name, into = c('CCLE_ID', 'condition'), sep = '\\.', remove = F) %>% 
    dplyr::mutate(condition = plyr::revalue(condition, replace = c(`Untreated_48hr` = 'Tram_0hr'))) %>% 
    tidyr::separate(condition, into = c('treat_cond', 'time_point'), sep = '_', remove = F) %>% 
    dplyr::mutate(condition = ifelse(treat_cond %in% c('DMSO'), 'DMSO', condition))
  sample_info %<>% left_join(CL_df, by = 'CCLE_ID')
  
  used_genes <- which(rowSums(summed_counts >= min_counts_per_gene) >= min_det_samps_per_gene)
  print(sprintf('%d used genes', length(used_genes)))
  dge <- edgeR::DGEList(summed_counts[used_genes,])
  dge <- edgeR::calcNormFactors(dge, method = 'TMMwzp')
  dd <- edgeR::cpm(dge, prior.count = globals$prior_cnt, log = TRUE)
  
  #FIT AVG MODEL
  design <- model.matrix(~condition + CCLE_ID, data = sample_info)
  fit <- limma::lmFit(dd, design)
  fit2 <- limma::eBayes(fit, trend = TRUE)
  ctests <- c('conditionTram_0hr', 
              'conditionTram_3hr', 
              'conditionTram_6hr', 
              'conditionTram_12hr',
              'conditionTram_24hr', 
              'conditionTram_48hr')
  all_res_avg <- ldply(ctests, function(cur_coef) {
    limma::topTable(fit2, number = Inf, coef = cur_coef, confint = TRUE) %>%
      rownames_to_column(var = 'Gene') %>% 
      mutate(condition = cur_coef)
  }) %>% dplyr::mutate(condition = str_match(condition, '(Tram_[0-9]+)hr')[,2],
                condition = factor(condition, levels = tp_levs))
  
  #repeat using mean-collapsed (rather than sum-collapsed) data
  fit_avg_collapse <- limma::lmFit(cond_avgs[used_genes,], design) %>% 
    limma::eBayes(trend = TRUE)
  all_res_avg_avg_collapse <- ldply(ctests, function(cur_coef) {
    limma::topTable(fit_avg_collapse, number = Inf, coef = cur_coef, confint = TRUE) %>%
      rownames_to_column(var = 'Gene') %>% 
      mutate(condition = cur_coef)
  }) %>% dplyr::mutate(condition = str_match(condition, '(Tram_[0-9]+)hr')[,2],
                       condition = factor(condition, levels = tp_levs))
  
  #FIT VIABILITY SIG MODEL
  design <- model.matrix(~condition + CCLE_ID + condition:sens, data = sample_info)
  colnames(design) %<>% make.names
  design <- design[, !grepl('conditionDMSO.sens', colnames(design))] #don't include DMSO-sensitivity interaction term 
  
  fit <- limma::lmFit(dd, design)
  fit2 <- limma::eBayes(fit, trend = TRUE)
  ctests_sens <- c('conditionTram_0hr.sens', 
                   'conditionTram_3hr.sens', 
                   'conditionTram_6hr.sens', 
                   'conditionTram_12hr.sens', 
                   'conditionTram_24hr.sens', 
                   'conditionTram_48hr.sens')
  all_res_slope <- ldply(ctests_sens, function(cur_coef) {
    limma::topTable(fit2, number = Inf, coef = cur_coef, confint = TRUE) %>%
      rownames_to_column(var = 'Gene') %>% 
      mutate(condition = cur_coef)
  }) %>% dplyr::mutate(condition = str_match(condition, '(Tram_[0-9]+)hr.sens')[,2],
                condition = factor(condition, levels = tp_levs))
  
  #fit model with mean-collapsed data
  fit_avg_collapse <- limma::lmFit(cond_avgs[used_genes,], design) %>% 
    limma::eBayes(trend = TRUE)
  all_res_slope_avg_collapse <- ldply(ctests_sens, function(cur_coef) {
    limma::topTable(fit_avg_collapse, number = Inf, coef = cur_coef, confint = TRUE) %>%
      rownames_to_column(var = 'Gene') %>% 
      mutate(condition = cur_coef)
  }) %>% dplyr::mutate(condition = str_match(condition, '(Tram_[0-9]+)hr.sens')[,2],
                       condition = factor(condition, levels = tp_levs))
  
  ctests_int <- c('conditionTram_0hr', 
                  'conditionTram_3hr', 
                  'conditionTram_6hr', 
                  'conditionTram_12hr', 
                  'conditionTram_24hr', 
                  'conditionTram_48hr')
  all_res_int <- ldply(ctests_int, function(cur_coef) {
    limma::topTable(fit2, number = Inf, coef = cur_coef, confint = TRUE) %>%
      rownames_to_column(var = 'Gene') %>% 
      mutate(condition = cur_coef)
  }) %>% dplyr::mutate(condition = str_match(condition, '(Tram_[0-9]+)hr')[,2],
                condition = factor(condition, levels = tp_levs))
  all_res_int_avg_collapse <- ldply(ctests_int, function(cur_coef) {
    limma::topTable(fit_avg_collapse, number = Inf, coef = cur_coef, confint = TRUE) %>%
      rownames_to_column(var = 'Gene') %>% 
      mutate(condition = cur_coef)
  }) %>% dplyr::mutate(condition = str_match(condition, '(Tram_[0-9]+)hr')[,2],
                       condition = factor(condition, levels = tp_levs))
  
  #compare sum and avg collapsed results
  comb_slope_results <- full_join(all_res_slope, all_res_slope_avg_collapse, by = c('Gene', 'condition'), suffix = c('_sum', '_mean')) %>% 
    dplyr::mutate(condition = paste0(str_match(condition, 'Tram_([0-9]+)')[,2], 'hr'),
                  condition = factor(condition, levels = c('0hr', '3hr', '6hr', '12hr', '24hr', '48hr')))
  comb_slope_results %>% 
    ggplot(aes(logFC_sum, logFC_mean)) + 
    geom_point(pch = 21, size = 2, color = 'white', fill = 'black', stroke = 0.1) + 
    geom_abline() + 
    labs(x = 'Effect Size (sum-collapse)',
         y = 'Effect Size (mean-collapse)') +
    ggpubr::stat_cor() + 
    theme_Publication() +
    facet_wrap(~condition) +
    ggtitle('Viability related')
  ggsave(file.path(fig_dir, 'tc_slope_sum_mean_compare.png'), width = 6, height = 5)
  ggsave(file.path(fig_dir, 'tc_slope_sum_mean_compare.pdf'), width = 6, height = 5)
  
  comb_int_results <- full_join(all_res_int, all_res_int_avg_collapse, by = c('Gene', 'condition'), suffix = c('_sum', '_mean')) %>% 
    dplyr::mutate(condition = paste0(str_match(condition, 'Tram_([0-9]+)')[,2], 'hr'),
                  condition = factor(condition, levels = c('0hr', '3hr', '6hr', '12hr', '24hr', '48hr')))
  comb_int_results %>% 
    ggplot(aes(logFC_sum, logFC_mean)) + 
    geom_point(pch = 21, size = 2, color = 'white', fill = 'black', stroke = 0.1) + 
    geom_abline() + 
    labs(x = 'Effect Size (sum-collapse)',
         y = 'Effect Size (mean-collapse)') +
    ggpubr::stat_cor() + 
    theme_Publication() +
    facet_wrap(~condition) +
    ggtitle('Viability independent')
  ggsave(file.path(fig_dir, 'tc_int_sum_mean_compare.png'), width = 6, height = 5)
  ggsave(file.path(fig_dir, 'tc_int_sum_mean_compare.pdf'), width = 6, height = 5)
  
  
  #HEATMAP PLOTS
  #get genes which are significantly down-regulated, sort by effect size
  int_genes <- all_res_int %>% 
    dplyr::filter(adj.P.Val < 0.1, logFC < 0) %>% 
    dplyr::arrange(logFC) %>% 
    .[['Gene']] %>% 
    unique() %>% 
    head(n_int_genes)
  
  slope_genes <- all_res_slope %>% 
    dplyr::filter(adj.P.Val < 0.1, logFC < 0) %>% 
    dplyr::arrange(logFC) %>%
    # arrange(desc(abs(logFC))) %>% 
    .[['Gene']] %>% 
    unique() %>% 
    head(n_slope_genes)
  
  all_res_int %>% 
    dplyr::mutate(condition = as.character(condition),
           condition = paste0(str_replace_all(condition, 'Tram_', ''), 'hr'),
           condition = factor(condition, levels = c('0hr', '3hr', '6hr', '12hr', '24hr', '48hr'))) %>% 
    dplyr::filter(Gene %in% int_genes) %>% 
    dplyr::mutate(Gene = factor(Gene, levels = int_genes)) %>% 
    ggplot(aes(condition, Gene, fill = logFC)) + 
    geom_tile() +
    xlab('Time post-treatment') +
    labs(fill = 'Effect Size') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    scale_fill_gradient2(limits = c(-3,3), oob = squish, low = muted("blue"), high = muted("red"))
  ggsave(file.path(fig_dir, 'tc_int_hmap.png'), width = 3.5, height = 3.5)
  ggsave(file.path(fig_dir, 'tc_int_hmap.pdf'), width = 3.5, height = 3.5) 
  
  all_res_slope %>% 
    dplyr::mutate(condition = as.character(condition),
           condition = paste0(str_replace_all(condition, 'Tram_', ''), 'hr'),
           condition = factor(condition, levels = c('0hr', '3hr', '6hr', '12hr', '24hr', '48hr'))) %>% 
    dplyr::filter(Gene %in% slope_genes) %>% 
    dplyr::mutate(Gene = factor(Gene, levels = slope_genes)) %>% 
    ggplot(aes(condition, Gene, fill = logFC)) + 
    geom_tile() +
    xlab('Time post-treatment') +
    guides(fill = guide_colorbar(title = 'Effect Size')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    scale_fill_gradient2(limits = c(-4,4), oob = squish, low = muted("blue"), high = muted("red"))
  ggsave(file.path(fig_dir, 'tc_slope_hmap.png'), width = 3.5, height = 3.5)
  ggsave(file.path(fig_dir, 'tc_slope_hmap.pdf'), width = 3.5, height = 3.5)
  

  #repeat above using mean-collapsed model estimates
  all_res_int_avg_collapse %>% 
    dplyr::mutate(condition = as.character(condition),
                  condition = paste0(str_replace_all(condition, 'Tram_', ''), 'hr'),
                  condition = factor(condition, levels = c('0hr', '3hr', '6hr', '12hr', '24hr', '48hr'))) %>% 
    dplyr::filter(Gene %in% int_genes) %>% 
    dplyr::mutate(Gene = factor(Gene, levels = int_genes)) %>% 
    ggplot(aes(condition, Gene, fill = logFC)) + 
    geom_tile() +
    xlab('Time post-treatment') +
    labs(fill = 'Effect Size') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    scale_fill_gradient2(limits = c(-3,3), oob = squish, low = muted("blue"), high = muted("red"))
  ggsave(file.path(fig_dir, 'tc_int_hmap_avg_collapsed.png'), width = 3.5, height = 3.5)
  ggsave(file.path(fig_dir, 'tc_int_hmap_avg_collapsed.pdf'), width = 3.5, height = 3.5)  
  
  all_res_slope_avg_collapse %>% 
    dplyr::mutate(condition = as.character(condition),
                  condition = paste0(str_replace_all(condition, 'Tram_', ''), 'hr'),
                  condition = factor(condition, levels = c('0hr', '3hr', '6hr', '12hr', '24hr', '48hr'))) %>% 
    dplyr::filter(Gene %in% slope_genes) %>% 
    dplyr::mutate(Gene = factor(Gene, levels = slope_genes)) %>% 
    ggplot(aes(condition, Gene, fill = logFC)) + 
    geom_tile() +
    xlab('Time post-treatment') +
    guides(fill = guide_colorbar(title = 'Effect Size')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    scale_fill_gradient2(limits = c(-4,4), oob = squish, low = muted("blue"), high = muted("red"))
  ggsave(file.path(fig_dir, 'tc_slope_hmap_avg_collapsed.png'), width = 3.5, height = 3.5)
  ggsave(file.path(fig_dir, 'tc_slope_hmap_avg_collapsed.pdf'), width = 3.5, height = 3.5)  
    
  #GSEA on these components
  get_condition_GSEA <- function(df, gsc, n_top_genes) {
    all_genes <- unique(df$Gene)
    ldply(unique(df$condition), function(cur_cond) {
      cur_hits_down <- df %>% 
        dplyr::filter(condition == cur_cond) %>% 
        dplyr::arrange(logFC) %>%
        head(n_top_genes) %>% 
        .[['Gene']] %>% 
        as.character()
      cur_hits_up <- df %>% 
        dplyr::filter(condition == cur_cond) %>% 
        dplyr::arrange(desc(logFC)) %>%
        head(n_top_genes) %>% 
        .[['Gene']] %>% 
        as.character()
      rbind(
        run_GSAhyper(cur_hits_down, all_genes, gsc) %>% 
          dplyr::mutate(condition = cur_cond,
                 direction = 'down'),
        run_GSAhyper(cur_hits_up, all_genes, gsc) %>% 
          dplyr::mutate(condition = cur_cond,
                 direction = 'up')
      )
    }) %>% 
      dplyr::mutate(gset_dir = paste0(gene_set, '.', direction),
             enrich = -log10(`p-value`),
             condition = as.character(condition),
             condition = paste0(str_replace_all(condition, 'Tram_', ''), 'hr'),
             condition = factor(condition, levels = c('0hr', '3hr', '6hr', '12hr', '24hr', '48hr')))
  }
  gsea_int <- get_condition_GSEA(all_res_int, gsc_data$hallmark, n_top_genes)
  gsea_slope <- get_condition_GSEA(all_res_slope, gsc_data$hallmark, n_top_genes)
  gsea_avg <- get_condition_GSEA(all_res_avg, gsc_data$hallmark, n_top_genes)
  
  gsea_int %>% 
    dplyr::filter(gset_dir %in% c('HALLMARK_KRAS_SIGNALING_UP.down')) %>% 
    dplyr::mutate(gene_set = str_replace(gene_set, 'HALLMARK_', '')) %>% 
    ggplot(aes(condition, enrich, group = gset_dir)) +
    geom_line() +
    geom_point(size = 3) + 
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.text = element_text(size = 9)) +
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    # guides(color = FALSE) +
    ylab('Gene set enrichment\n-log10(p-value)') +
    xlab('Time post-treatment')
  ggsave(file.path(fig_dir, 'trametinib_tc_int_gsea_new.png'),
         width = 3, height = 2.25)
  ggsave(file.path(fig_dir, 'trametinib_tc_int_gsea_new.pdf'),
         width = 3, height = 2.25)
  
  gsea_slope %>% 
    dplyr::filter(gset_dir %in% c('HALLMARK_G2M_CHECKPOINT.down')) %>% 
    dplyr::mutate(gene_set = str_replace(gene_set, 'HALLMARK_', '')) %>% 
    ggplot(aes(condition, enrich, group = gset_dir)) +
    geom_line() +
    geom_point(size = 3) + 
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.text = element_text(size = 9)) +
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    guides(color = guide_legend(nrow = 3, title = element_blank())) +
    # guides(color = FALSE) +
    ylab('Gene set enrichment\n-log10(p-value)') +
    xlab('Time post-treatment') 
  ggsave(file.path(fig_dir, 'trametinib_tc_slope_gsea_new.png'),
         width = 3, height = 2.25)
  ggsave(file.path(fig_dir, 'trametinib_tc_slope_gsea_new.pdf'),
         width = 3, height = 2.25)
  
  #repeat with avg-collapsing  
  gsea_int_avg_collapse <- get_condition_GSEA(all_res_int_avg_collapse, gsc_data$hallmark, n_top_genes)
  gsea_slope_avg_collapse <- get_condition_GSEA(all_res_slope_avg_collapse, gsc_data$hallmark, n_top_genes)

  gsea_int_avg_collapse %>% 
    dplyr::filter(gset_dir %in% c('HALLMARK_KRAS_SIGNALING_UP.down')) %>% 
    dplyr::mutate(gene_set = str_replace(gene_set, 'HALLMARK_', '')) %>% 
    ggplot(aes(condition, enrich, group = gset_dir)) +
    geom_line() +
    geom_point(size = 3) + 
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.text = element_text(size = 9)) +
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    # guides(color = FALSE) +
    ylab('Gene set enrichment\n-log10(p-value)') +
    xlab('Time post-treatment')
  ggsave(file.path(fig_dir, 'trametinib_tc_int_gsea_new_avg_collapse.png'),
         width = 3, height = 2.25)
  ggsave(file.path(fig_dir, 'trametinib_tc_int_gsea_new_avg_collapse.pdf'),
         width = 3, height = 2.25)
  
  gsea_slope_avg_collapse %>% 
    dplyr::filter(gset_dir %in% c('HALLMARK_G2M_CHECKPOINT.down')) %>% 
    dplyr::mutate(gene_set = str_replace(gene_set, 'HALLMARK_', '')) %>% 
    ggplot(aes(condition, enrich, group = gset_dir)) +
    geom_line() +
    geom_point(size = 3) + 
    theme_Publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.text = element_text(size = 9)) +
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    guides(color = guide_legend(nrow = 3, title = element_blank())) +
    # guides(color = FALSE) +
    ylab('Gene set enrichment\n-log10(p-value)') +
    xlab('Time post-treatment') 
  ggsave(file.path(fig_dir, 'trametinib_tc_slope_gsea_new_avg_collapse.png'),
         width = 3, height = 2.25)
  ggsave(file.path(fig_dir, 'trametinib_tc_slope_gsea_new_avg_collapse.pdf'),
         width = 3, height = 2.25)
  
  
  # CC ARREST
  cell_df %<>% dplyr::mutate(condition = plyr::revalue(hash_tag, replace = c(`Untreated_48hr` = 'Tram_0hr')))
  df <- cell_df %>% 
    dplyr::mutate(condition = as.character(condition)) %>% 
    dplyr::mutate(condition = ifelse(grepl('DMSO', condition), 'DMSO', condition)) %>% 
    dplyr::mutate(is_G1 = Phase == "G1") %>% 
    left_join(CL_df, by = c('singlet_ID' = 'CCLE_ID'))
  
  #compute avg G1-fraction per cell line and time point
  avgs <- df %>% 
    dplyr::filter(!is.na(condition)) %>% 
    dplyr::group_by(singlet_ID, condition) %>% 
    dplyr::summarise(G1_prop = mean(is_G1, na.rm=T),
              n = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(n >= min_cells_per_cond)
  
  #compute and subtract baseline G1-proportion per cell line
  avgs %<>% left_join(avgs %>% 
                        dplyr::filter(condition == 'DMSO') %>% 
                        dplyr::select(singlet_ID, baseline = G1_prop),
                      by = 'singlet_ID') %>% 
    dplyr::mutate(G1_prop = G1_prop - baseline) %>% 
    dplyr::group_by(condition) %>% 
    dplyr::summarise(G1_prop_avg = mean(G1_prop, na.rm=T),
              sd = sd(G1_prop, na.rm=T),
              n = sum(!is.na(G1_prop)),
              se = sd/sqrt(n)) %>% 
    dplyr::mutate(condition = str_replace_all(condition, 'Tram_', ''),
           condition = factor(condition, levels = c('DMSO', '0hr', '3hr', '6hr', '12hr', '24hr', '48hr')))
  
  ggplot(avgs %>% dplyr::filter(condition != 'DMSO'), 
         aes(condition, G1_prop_avg)) + 
    geom_point() + 
    geom_line(aes(group = 1)) +
    geom_errorbar(aes(ymax = G1_prop_avg + se, 
                      ymin = G1_prop_avg - se),
                  width = 0.25) +
    ylab('Delta G0/G1 fraction') +
    theme_Publication() + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          axis.title.x = element_blank())
  ggsave(file.path(fig_dir, 'tc_G1_frac.png'), width = 3.5, height = 3)
  ggsave(file.path(fig_dir, 'tc_G1_frac.pdf'), width = 3.5, height = 3)
  
  
  #FIGURE WITH ALL VOLCANOS
  sample_info <- data_frame(sample_name = colnames(summed_counts)) %>%
    tidyr::separate(sample_name, into = c('CCLE_ID', 'condition'), sep = '\\.', remove = F) %>%
    tidyr::separate(condition, into = c('treat_cond', 'time_point'), sep = '_', remove = F)
  sample_info %<>% left_join(CL_df)
  
  design <- model.matrix(~0 + condition + CCLE_ID, data = sample_info)
  fit <- limma::lmFit(dd, design)
  ctests <- c('conditionTram_3hr',
              'conditionTram_6hr',
              'conditionTram_12hr',
              'conditionTram_24hr',
              'conditionTram_48hr',
              'conditionDMSO_3hr',
              'conditionDMSO_6hr',
              'conditionDMSO_12hr',
              'conditionDMSO_24hr',
              'conditionDMSO_48hr')
  cm <- laply(ctests, function(cur_coef) {
    contr_vec <- (colnames(design) == cur_coef) - (colnames(design) == 'conditionUntreated_48hr')
  }) %>%
    set_rownames(ctests) %>%
    set_colnames(colnames(design)) %>%
    t()
  fit <- limma::contrasts.fit(fit, contrasts = cm)
  fit2 <- limma::eBayes(fit, trend = TRUE)
  
  
  all_res_avg <- ldply(ctests, function(cur_coef) {
    limma::topTable(fit2, number = Inf, coef = cur_coef) %>%
      rownames_to_column(var = 'Gene') %>%
      dplyr::mutate(condition = cur_coef)
  }) %>%
    dplyr::mutate(condition = str_replace_all(condition, 'condition', ''),
           condition = factor(condition, levels = c('DMSO_3hr', 'DMSO_6hr', 'DMSO_12hr', 'DMSO_24hr', 'DMSO_48hr', 'Tram_3hr', 'Tram_6hr', 'Tram_12hr', 'Tram_24hr', 'Tram_48hr'))) %>% 
    dplyr::mutate(is_sig = adj.P.Val < 0.1)
  
  ggplot(all_res_avg, aes(logFC, -log10(P.Value))) +
    geom_point(data = all_res_avg %>% filter(is_sig), pch = 21, fill = 'darkred', color = 'white', stroke = 0.1, size = 2) +
    geom_point(data = all_res_avg %>% filter(!is_sig), pch = 21, fill = 'black', color = 'white', stroke = 0.1, size = 2) +
    geom_text_repel(data = all_res_avg %>%
                      filter(abs(logFC) > 1, adj.P.Val < 0.1) %>%
                      group_by(condition) %>%
                      top_n(n = 5, wt = abs(logFC)),
                    aes(label = Gene), size = 2.5) +
    theme_Publication() +
    facet_wrap(~condition, nrow = 2)
  ggsave(file.path(fig_dir, 'tramet_tc_all_volcanos.png'),
         width = 9, height = 4)
  ggsave(file.path(fig_dir, 'tramet_tc_all_volcanos.pdf'),
         width = 9, height = 4)
  }
