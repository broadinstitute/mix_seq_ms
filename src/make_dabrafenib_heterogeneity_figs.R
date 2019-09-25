make_dabrafenib_heterogeneity_figs <- function() {
  cur_expt <- sc_DE_meta$dabrafenib_24hr_expt3
  
  #example cell lines
  CL1 <- 'A375_SKIN'
  CL2 <- 'IGR1_SKIN'
  
  seuObj <- load_sc_data(cur_expt, sc_expts = sc_expts)
  all_CLs <- unique(seuObj@meta.data$singlet_ID)
  
  CL_annotations <- all_CL_features[[cur_expt$expt_name]]
  #merge in primary site annotations
  CL_meta = load.from.taiga(
    data.name='master-cell-line-export-0306',
    data.version=435, 
    data.file='masterfile_2019-09-23') %>% 
    dplyr::select(DEPMAP_ID = DepMap_ID, Disease, Subtype = 'Disease Subtype') %>% 
    dplyr::distinct(DEPMAP_ID, .keep_all=T)
  CL_annotations %<>% left_join(CL_meta, by = 'DEPMAP_ID')
  CL_annotations %<>% 
    dplyr::mutate(BRAF_mel = ifelse((BRAF_MUT > 0) & (Subtype == 'melanoma'), TRUE, FALSE))
  BRAF_mut_mel <- CL_annotations %>% dplyr::filter(BRAF_mel, CCLE_ID %in% all_CLs) %>% .[['CCLE_ID']]
  

  #now compare the differential response between two cell lines
  tt_comb <- compare_CL_responses(seuObj, CL1, CL2)
  ggplot(tt_comb %>% mutate(FDR = FDR_diff), 
         aes(logFC_CL1, logFC_CL2, color = -log10(FDR))) + 
    xlab(sprintf('logFC (%s)', CL1)) + 
    ylab(sprintf('logFC (%s)', CL2)) +
    geom_point(alpha = 0.8, size = 1) + 
    geom_abline(linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    geom_vline(xintercept = 0, linetype = 'dashed') +
    ggrepel::geom_label_repel(data = . %>% 
                                dplyr::filter(FDR_diff < 0.1) %>% 
                                arrange(dplyr::desc(abs(logFC_diff))) %>% 
                                head(15),
                              aes(label = Gene), 
                              size = 2) +
    scale_color_gradient(limits = c(0, 8) ,oob = scales::squish, breaks = c(0, 5), low = 'darkgrey', high = 'red') +
    cdsr::theme_Publication() 
  ggsave(file.path(fig_dir, sprintf('%s_%s_dabrafenib_compare.png', CL1, CL2)),
         width = 4, height = 4)

  #make heatmap 
  dat <- load_compute_CPM_stats(cur_expt, results_dir, type = 'avg_cpm', use_CLs = BRAF_mut_mel, prior_counts = globals$pca_prior_cnt)
  top_genes <- apply(dat$LFC_mat, 1, sd, na.rm=T) %>% 
    sort(decreasing = T) %>% 
    names() %>% 
    head(50)
  make_LFC_heatmap(LFC_mat = dat$LFC_mat, gene_list = top_genes,
                   transpose = TRUE,
                   cluster_cols = TRUE,
                   cluster_rows = TRUE,
                   filename = file.path(fig_dir, 'BRAF_melanoma_dabrafenib_heatmap.png'),
                   width = 7, height = 3,
                   fontsize = 8, treeheight_row = 15, treeheight_col = 15)

}
