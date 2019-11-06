make_viability_signature_heatmaps <- function() {
  
  min_cells <- 5
  n_genes <- 1000
  
  all_res <- ldply(sc_DE_meta[setdiff(names(sc_DE_meta), 'GPX4_expt2')], function(cur_expt) {
    
    sprintf('Processing expt %s', cur_expt$expt_name)
    out_dir <- file.path(results_dir, cur_expt$expt_name)
    CL_df <- all_CL_features[[cur_expt$expt_name]]
   
    cell_df <- load_cell_info(cur_expt, QC_filter = TRUE) %>% 
      tidyr::separate(col = 'condition', into = c('condition', 'batch'), sep = '_')
    usable_cell_lines <- cell_df %>% 
      dplyr::group_by(singlet_ID, condition) %>% 
      dplyr::summarise(n = n()) %>% 
      reshape2::dcast(singlet_ID ~ condition, value.var = 'n') %>% 
      dplyr::filter(control >= min_cells, treat >= min_cells) %>% 
      .[['singlet_ID']] %>% 
      intersect(unique(CL_df$CCLE_ID))
    
    dat <- load_compute_CPM_stats(cur_expt, results_dir, prior_counts = globals$prior_cnt, type = 'avg_cpm', use_CLs = usable_cell_lines)
    LFC_mat <- dat$LFC_mat %>% t()
    
    #get viability-LFC correlations
    res <- biomarkVS::corr.test(LFC_mat, CL_df[match(rownames(LFC_mat), CL_df$CCLE_ID), 'sens'], adjust = 'none') 
    res_df <- data_frame(Gene = names(res$r), r = res$r, p = res$p) %>% 
      mutate(q = p.adjust(p, method = 'BH'), 
             expt = cur_expt$expt_name)
    
    #also join in average and slope of lm fits
    limma_res <- read_rds(file.path(out_dir, 'limma_res.rds'))
    res_df %<>% left_join(limma_res$res_slope, by = 'Gene') %>% 
      left_join(limma_res$res_avg, by = 'Gene', suffix = c('', '_avg')) 
    return(res_df)
  })
  
  #all 24hr post-pert
  used_expts_selective <- c(
    'Idasanutlin_24hr_expt1',
    'Trametinib_24hr_expt3', 
    'Dabrafenib_24hr_expt3',
    'BRD3379_24hr_expt3',
    'Afatinib_expt10',
    'Navitoclax_24hr_expt3',
    'AZD5591_expt10',
    'Everolimus_expt10',
    'Taselisib_expt10'
    )
  used_chemos <- c(
    'Gemcitabine_expt10',
    'Prexasertib_expt10',
    'Bortezomib_24hr_expt1',
    'JQ1_expt10')
  
  #get per-gene stats 
  gavgs <- all_res %>% 
    dplyr::filter(expt %in% c(used_expts_selective, used_chemos)) %>% 
    dplyr::group_by(Gene) %>% 
    dplyr::summarise(mean_ar = mean(abs(r), na.rm=T),
              mean_r = mean(r, na.rm=T),
              num_NA = sum(is.na(logFC_avg)),
              n = sum(!is.na(t))) %>% 
    dplyr::filter(num_NA == 0)
  
  used_genes <- all_res %>% 
    dplyr::filter(Gene %in% unique(gavgs$Gene)) %>% 
    dplyr::filter(expt %in% used_expts_selective, q < globals$q_thresh) %>%
    dplyr::arrange(dplyr::desc(abs(r))) %>%
    dplyr::distinct(Gene) %>%
    head(n_genes) %>%
    .[['Gene']]
  
  #create matrix of correlation values for selective compounds
  mat_sel <- all_res %>% 
    dplyr::filter(Gene %in% used_genes, expt %in% used_expts_selective) %>% 
    dplyr::mutate(expt = str_match(expt, '(^[:alnum:]+)_')[,2]) %>% 
    reshape2::acast(expt ~ Gene, value.var = 'r') %>% 
    t() %>% 
    scale(center = F, scale = F) %>% 
    t()
  
  #create matrix of avg LFCs for chemo compounds
  mat_chem <- all_res %>%
    filter(Gene %in% used_genes, expt %in% used_chemos) %>%
    mutate(expt = str_match(expt, '(^[:alnum:]+)_')[,2]) %>% 
    reshape2::acast(expt ~ Gene, value.var = 'logFC_avg') %>%
    t() %>%
    scale(center = F, scale = F) %>%
    t()
  
  
  library(circlize)
  library(ComplexHeatmap)
  library(GSEABase)
  height_per <- 0.6
  hmap_width <- 10
  col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
  hmap_sel <- Heatmap(mat_sel, col = col_fun, show_column_names = F, 
          width = unit(hmap_width, 'cm'),
          height = unit(height_per*length(used_expts_selective), 'cm'),
          column_split = NULL,
          border = TRUE,
          row_title = 'Selective',
          column_title = 'Genes',
          heatmap_legend_param = list(title = 'Sensitivity-LFC\ncorrelation', 
                                      legend_height = unit(3, 'cm'),
                                      legend_width = unit(0.5, 'cm')))
  
  
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  hmap_chem <- Heatmap(mat_chem, 
                      col = col_fun, 
                      column_order = unlist(column_order(hmap_sel)),
                      show_column_names = F, 
                      border = TRUE,
                      row_title = 'pan-toxic',
                      column_title = 'Genes',
                      width = unit(hmap_width, 'cm'),
                      height = unit(height_per*length(used_chemos), 'cm'),
                      row_gap = unit(10, "mm"),
                      heatmap_legend_param = list(title = 'avg logFC', 
                                                  legend_height = unit(3, 'cm'),
                                                  legend_width = unit(0.5, 'cm')))
  
  
  #add gene set annotations
  cc_genes <- colnames(mat_sel) %in% geneIds(gsc_data$hallmark[['HALLMARK_G2M_CHECKPOINT']])
  translation_genes <- colnames(mat_sel) %in% geneIds(gsc_data$canonical[['REACTOME_TRANSLATION']])
  p53_genes <- colnames(mat_sel) %in% geneIds(gsc_data$hallmark[['HALLMARK_P53_PATHWAY']]) 
  gs_mat <- (0 + rbind(cc_genes, translation_genes, p53_genes)) %>% 
    set_rownames(c('cell cycle', 'translation', 'P53 signaling'))
  hmap_gs <- Heatmap(gs_mat, name = "gene sets", 
          col = c("0" = "white", "1" = "black"), 
          show_heatmap_legend = FALSE, 
          height = unit(2, "cm"), 
          row_title = 'gene sets',
          cluster_rows = F, 
          border = T,
          cluster_columns = F) 
  hmap <- hmap_sel %v% hmap_chem %v% hmap_gs
  pdf(file.path(fig_dir, 'viability_sig_heatmap.pdf'),
      width = 8, height = 5)
  draw(hmap)
  dev.off()
  
   

  #make correlation matrix plot  
  library(Hmisc)
  comb_mat <- rbind(mat_sel[row_order(hmap_sel),], mat_chem[row_order(hmap_chem),])
  col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
  cmat <- cor(t(comb_mat))
  pvals <- rcorr(t(comb_mat))$P
  diag(pvals) <- 1
  cor_hmap <- Heatmap(cmat, col = col_fun, cluster_rows = F, cluster_columns = F,
          heatmap_legend_param = list(title = 'Correlation', 
                                      legend_height = unit(4, 'cm'),
                                      legend_width = unit(0.5, 'cm')),
          rect_gp = gpar(col = "white", lwd = 1)
          )
  pdf(file.path(fig_dir, 'viability_sig_corr_heatmap.pdf'),
      width = 5, height = 3.7)
  draw(cor_hmap)
  dev.off()

}
