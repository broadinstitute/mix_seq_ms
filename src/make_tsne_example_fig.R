make_tsne_example_fig <- function(expt_params) {

  #PARAMS
  show_labels <- FALSE
  vtr <- c()
  fwidth <- 5
  fheight <- 4.5

  
  expt_name <- expt_params$expt_name
  seuObj <- load_sc_data(expt_params, sc_expts)
  CL_annotations <- all_CL_features[[expt_name]]
  
  seuObj <- NormalizeData(object = seuObj, 
                          normalization.method = "LogNormalize", 
                          scale.factor = globals$scale_fac)
  
  seuObj <- CellCycleScoring(object = seuObj,
                             s.features = Seurat::cc.genes$s.genes,
                             g2m.features = Seurat::cc.genes$g2m.genes,
                             set.ident = FALSE)
  seuObj <- relabel_cell_cycle_phase(seuObj)
  
  cell_df <- seuObj@meta.data %>% 
    rownames_to_column(var = 'barcode') %>% 
    tidyr::separate(col = 'condition', into = c('treat_cond', 'batch'), sep = '_')
  seuObj[['treat_cond']] <- cell_df$treat_cond
  
  seuObj <- ScaleData(object = seuObj, vars.to.regress = vtr)
  
  seuObj <- FindVariableFeatures(object = seuObj,
                              nfeatures = globals$n_highvar_genes,
                              do.plot = FALSE,
                              selection.method = 'vst')
  
  n_CLs <- length(levels(seuObj))
  n_PCs <- n_CLs * globals$n_pcs_per_CL
  seuObj <- RunPCA(object = seuObj,
                   features = VariableFeatures(seuObj),
                   seed.use = 1,
                   npcs = n_PCs,
                   verbose = FALSE)
  
  seuObj <- RunTSNE(object = seuObj, dims = 1:n_PCs, check_duplicates = FALSE, 
                    seed.use = 0, perplexity = globals$tsne_perplexity)
  # seuObj <- RunUMAP(object = seuObj, dims = 1:n_PCs, min_dist = 0.7, n_neighbors = 100)
  
  dred <- 'tsne'
  df <- Embeddings(object = seuObj, reduction = dred) %>% 
    as.data.frame() %>% 
    set_colnames(c('t1', 't2')) %>% 
    rownames_to_column(var = 'barcode') %>% 
    cbind(seuObj@meta.data) %>% 
    tidyr::separate(col = 'condition', into = c('condition', 'batch'), sep = '_') %>% 
    dplyr::mutate(condition = factor(condition)) 
  
  #compute avg coordinates per cell line and condition
  avgs <- full_join(
    df %>% 
      dplyr::filter(condition == 'control') %>% 
      dplyr::group_by(singlet_ID) %>% 
      dplyr::summarise(ct1 = median(t1, na.rm=T),
                ct2 = median(t2, na.rm=T)),
    df %>% 
      dplyr::filter(condition == 'treat') %>% 
      dplyr::group_by(singlet_ID) %>% 
      dplyr::summarise(tt1 = median(t1, na.rm=T),
                tt2 = median(t2, na.rm=T))
  ) %>% dplyr::mutate(CL_short = str_match(singlet_ID, '^([:alnum:]+)_')[,2])
  
  #add in TP53 status
  df %<>% dplyr::mutate(TP53_WT = singlet_ID %in% globals$TP53_WT_cls_pool22)
  avgs %<>% dplyr::mutate(TP53_WT = singlet_ID %in% globals$TP53_WT_cls_pool22)
  
  # #make control-only tsne
  # g <- ggplot(df %>% filter(condition == 'control'), aes(t1, t2)) +
  #   geom_point(aes(color = singlet_ID), alpha = 0.75, size = 0.6) +
  #   guides(color = F) +
  #   xlab('tSNE 1') + ylab('tSNE 2') 
  # if (show_labels) {
  #   g <- g + geom_point(data = avgs, mapping = aes(x = ct1, y = ct2), size = 0, alpha = 0)+
  #     geom_text(data = avgs, mapping = aes(x = ct1, y = ct2, label = singlet_ID), size = 1.5)
  # }
  # g
  # ggsave(file.path(fig_dir, paste0(expt_name, '_', 'control_tsne.png')), width = fwidth, height = fheight, dpi = 400)
  
  # #make combined tsne colored by cell line
  # g <- ggplot(df %>% filter(cell_quality == 'normal'), aes(t1, t2)) +
  #   geom_point(aes(alpha = condition, color = singlet_ID, size = condition)) +
  #   scale_alpha_manual(values = c(control = 0.25, treat = 0.8)) +
  #   scale_size_manual(values = c(control = 0.1, treat = 0.3)) +
  #   guides(color = F, alpha = F, size = F) +
  #   xlab('tSNE 1') + ylab('tSNE 2') +
  #   geom_segment(data = avgs, aes(x = ct1, y = ct2, xend = tt1, yend = tt2),
  #                arrow = arrow(length = unit(0.2,"cm")), size = 1) 
  # if (show_labels) {
  #   g <- g + geom_point(data = avgs, mapping = aes(x = ct1, y = ct2), size = 0, alpha = 0) + 
  #     geom_text(data = avgs, mapping = aes(x = ct1, y = ct2, label = singlet_ID), size = 1.5)
  # }
  # g
  # ggsave(file.path(fig_dir, paste0(expt_name, '_', 'combined_tsne.png')), width = fwidth, height = fheight, dpi = 400)
  
  #make combined tsne colored by treatment condition
  df %<>% dplyr::mutate(tcond = plyr::revalue(condition, c(`control` = 'DMSO', `treat` = expt_params$drug_name)))
  g <- ggplot(df %>% filter(cell_quality == 'normal'), aes(t1, t2)) +
    geom_point(  aes(fill = tcond), pch = 21, size = 0.75, color = 'white', alpha = 0.8, stroke = 0.1) +
    guides(alpha = F, size = F, fill = guide_legend(title = element_blank(), override.aes = list(size = 3))) +
    scale_fill_manual(values = c('darkgray', 'darkred')) + 
    xlab('tSNE 1') + ylab('tSNE 2') +
    geom_segment(data = avgs, aes(x = ct1, y = ct2, xend = tt1, yend = tt2),
                 arrow = arrow(length = unit(0.2,"cm")), size = 0.75) +
    cdsr::theme_Publication()
  if (expt_params$drug_name == 'idasanutlin') {
    g <- g + geom_point(data = avgs, mapping = aes(x = ct1, y = ct2), size = 0, alpha = 0) + 
      geom_text(data = avgs %>% filter(!TP53_WT), mapping = aes(x = ct1, y = ct2, label = CL_short), color = 'black', size = 2.5) +
      geom_text(data = avgs %>% filter(TP53_WT), mapping = aes(x = ct1, y = ct2, label = CL_short), color = 'darkblue', size = 2.5)
    }
  g
  ggsave(file.path(fig_dir, paste0(expt_name, '_', 'combined_tsne_condcol.png')), width = fwidth, height = fheight, dpi = 400)
  
  #cell cycle phase example (specifically for nutlin)
  df %<>% dplyr::mutate(Phase = factor(Phase, levels = c('G0/G1', 'S', 'G2/M')))
  g <- ggplot(df %>% 
                # dplyr::filter(cell_quality == 'normal', condition == 'treat'),
              dplyr::filter(cell_quality == 'normal'),
              aes(t1, t2)) +
    geom_point(aes(fill = Phase, size = condition, alpha = condition), pch = 21, color = 'white', alpha = 0.7, stroke = 0.1) +
    scale_alpha_manual(values = c(`control` = 0.4, `treat` = 0.8)) +
    # scale_size_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.3)) +
    scale_size_manual(values = c(`control` = 0.4, `treat` = 0.8)) +
    guides(alpha = F, size = F, fill = guide_legend(override.aes = list(size = 3), nrow = 3),
           color = guide_legend(nrow = 2)) +
    xlab('tSNE 1') + ylab('tSNE 2') +
    cdsr::theme_Publication() + 
    cdsr::scale_fill_Publication()+
    theme(legend.text=element_text(size=8)) +
    geom_segment(data = avgs, aes(x = ct1, y = ct2, xend = tt1, yend = tt2, color = TP53_WT),
                  arrow = arrow(length = unit(0.2,"cm")), size = 0.7, alpha = 0.75) 
  g
  ggsave(file.path(fig_dir, paste0(expt_name, '_', 'combined_tsne_CCP.png')), width = 4.2, height = 4.25)

}

