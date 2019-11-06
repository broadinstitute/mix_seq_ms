make_bortezomib_subpop_figs <- function(targ_CLs) {
  
  #PARAMS
  n_PCs <- 10
  use_SCTransform <- TRUE
  prior_cnt <- 0.125
  clust_k <- 20
  clust_res <- 0.25
  
  seuObj <- load_sc_data(sc_DE_meta$bortezomib_24hr_expt1, sc_expts = sc_expts)
  
  seuObj <- NormalizeData(object = seuObj,
                          normalization.method = "LogNormalize",
                          scale.factor = globals$scale_fac)
  seuObj <- CellCycleScoring(object = seuObj,
                             s.features = Seurat::cc.genes$s.genes,
                             g2m.features = Seurat::cc.genes$g2m.genes,
                             set.ident = FALSE)
  seuObj <- relabel_cell_cycle_phase(seuObj)
  
  make_subpop_plots <- function(cur_CL) {
    si <- FetchData(seuObj, vars = 'singlet_ID')
    seuSub <- seuObj[, which(si == cur_CL)]
    if (use_SCTransform) {
      seuSub <- SCTransform(seuSub, variable.features.n = globals$n_highvar_genes)
    } else {
      seuSub <- ScaleData(object = seuSub)
      seuSub <- FindVariableGenes(object = seuSub,
                                  top.genes = globals$n_highvar_genes,
                                  do.plot = FALSE,
                                  selection.method = 'vst')
    }
    
    seuSub <- RunPCA(object = seuSub,
                     features = VariableFeatures(seuSub),
                     npcs = n_PCs,
                     do.print = FALSE, 
                     verbose = F)
    
    seuSub <- RunTSNE(object = seuSub, dims = 1:n_PCs, perplexity = globals$tsne_perplexity, seed.use = 1)
    
    controlSub <- subset(seuSub, subset = condition %in% c('control_1', 'control_2'))
    treatSub <- subset(seuSub, subset = condition == 'treat_1')
    
    #apply clustering within treated cells
    treatSub <- FindNeighbors(treatSub, dims = 1:n_PCs, k.param = clust_k, verbose = FALSE)  
    treatSub <- Seurat::FindClusters(treatSub, resolution = clust_res, verbose = FALSE)
    
    d_use <- 'tsne'
    df <- Embeddings(object = seuSub, reduction = d_use) %>% 
      as.data.frame() %>% 
      .[, c(1,2)] %>% 
      set_colnames(c('t1', 't2')) %>% 
      rownames_to_column(var = 'barcode') %>% 
      cbind(seuSub@meta.data) %>% 
      dplyr::mutate(condition = ifelse(condition %in% c('control_1', 'control_2'), 'DMSO', condition),
             condition = plyr::revalue(condition, c(treat_1 = 'Bort. 24hr', treat_2 = 'Bort. 6hr')),
             condition = factor(condition, levels = c('DMSO', 'Bort. 6hr', 'Bort. 24hr'))) %>% 
      left_join(data_frame(barcode = names(Idents(treatSub)),
                           cluster = Idents(treatSub)), by = 'barcode') %>% 
      dplyr::mutate(cluster = plyr::revalue(cluster, replace = c(`0` = 1, `1` = 2)),
             class = ifelse(is.na(cluster), 'DMSO', paste0('bortezomib_', cluster)),
             class = factor(class, levels = c('DMSO', 'bortezomib_1', 'bortezomib_2')),
             Phase = factor(Phase, levels = c('G0/G1', 'S', 'G2/M')))
    
    stopifnot(nlevels(treatSub) == 2) #assumes that there are two clusters in the treated population
    
    #set cluster naming according to G1 proportion
    clust_phase <- df %>% 
      dplyr::filter(!is.na(cluster)) %>% 
      dplyr::group_by(cluster) %>% 
      dplyr::mutate(is_G1 = Phase == 'G0/G1') %>% 
      dplyr::summarise(avg_G1 = mean(is_G1))
    if (clust_phase %>% dplyr::arrange(dplyr::desc(avg_G1)) %>% head(1) %>% .[['cluster']] == 2) {
      df %<>% dplyr::mutate(cluster = plyr::revalue(cluster, replace = c(`1` = 2, `2` = 1)))
    }
    df %<>% dplyr::mutate(cluster = factor(cluster, levels = c(1, 2)))
    
    xr <- with(df %>% filter(cell_quality == 'normal'), range(df$t1))
    yr <- with(df %>% filter(cell_quality == 'normal'), range(df$t2))
    avgs <- df %>% 
      dplyr::filter(condition == 'Bort. 24hr') %>% 
      dplyr::group_by(cluster) %>% 
      dplyr::summarise(t1 = median(t1), t2 = median(t2)) %>% 
      dplyr::mutate(cluster_label = paste0('cluster ', cluster))
    g <- ggplot(df %>% dplyr::filter(cell_quality == 'normal'),
                     aes(t1, t2)) +
      geom_point(aes(fill = condition), pch = 21, color = 'white', stroke = 0.2, size = 2.5) +
      coord_cartesian(xlim = xr, ylim = yr, clip = 'off') +
      guides(alpha = F, size = F) +
      xlab(paste0(d_use, ' 1')) + ylab(paste0(d_use, ' 2')) +
      scale_fill_manual(values = c(`DMSO` = 'darkgray',
                                    `Bort. 24hr` = 'indianred4')) +
      geom_text(data = avgs, aes(x = t1-2, y = t2 + 5, label = cluster_label), size = 8) +
      cdsr::theme_Publication()
    
    g_phase <- ggplot(df %>% filter(cell_quality == 'normal'),
                      aes(t1, t2)) +
      geom_point(aes(stroke = condition, fill = Phase, color = condition), pch = 21, size = 2.5, alpha = 0.75) +
      guides(alpha = F, stroke = F,
             fill = guide_legend(nrow = 3, override.aes = list(stroke = 0, size =3)), 
             color = guide_legend(nrow = 2, override.aes = list(stroke = 1, size = 3),
                                  title = element_blank())) +
      scale_color_manual(values = c('darkgray', 'indianred4')) +
      scale_discrete_manual(aesthetics = "stroke", values = c(0.5, 1)) +
      xlab(paste0(d_use, ' 1')) + ylab(paste0(d_use, ' 2')) +
      geom_text(data = avgs, aes(x = t1-2, y = t2 + 5, label = cluster_label), size = 8) +
      coord_cartesian(xlim = xr, ylim = yr, clip = 'off') +
      cdsr::scale_fill_Publication() +
      cdsr::theme_Publication()
    

    ggsave(file.path(fig_dir, sprintf('%s_bortezomib_subpop_phase.png', cur_CL)), width = 3.5, height = 4)
    ggsave(file.path(fig_dir, sprintf('%s_bortezomib_subpop.png', cur_CL)), plot = g, width = 4, height = 4)
    
    
    #Run DE analysis comparing the two treatment clusters
    used_genes <- which(Matrix::rowSums(GetAssayData(treatSub, 'counts')[, rownames(treatSub@meta.data)] > 0) > nrow(treatSub@meta.data)*globals$min_frac_cells_det)
    dge <- DGEList(GetAssayData(treatSub, 'counts')[used_genes,rownames(treatSub@meta.data)])
    dge <- calcNormFactors(dge, method = 'TMMwzp')
    design <- model.matrix(~cluster + cell_det_rate , data = df)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design, prior.count = prior_cnt)
    qlf <- glmQLFTest(fit, coef = 'cluster2')
    tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
    tt_res %>% 
      ggplot(aes(logFC, -log10(PValue))) +
      geom_point(aes(fill = FDR < 0.1), pch = 21, stroke = 0.1, color = 'white', size = 2) +
      geom_label_repel(data = . %>% arrange(dplyr::desc(abs(logFC))) %>% head(10),
                       aes(label = Gene),
                       size = 2.5) +
      geom_vline(xintercept = 0, linetype = 'dashed') + 
      geom_hline(yintercept = 0, linetype = 'dashed') +
      cdsr::theme_Publication() +
      scale_fill_manual(values = c(`FALSE` = 'black', `TRUE` = 'darkred'))
    ggsave(file.path(fig_dir, sprintf('%s_bortezomib_treat_cluster_volcano.png', cur_CL)),
           width = 3.5, height = 3)
    
    gene_stat <- tt_res$logFC %>% set_names(tt_res$Gene)
    gsea_plot <- make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = 50, n_lab_per = 5, lab_size = 2.5, max_chars = 30)  
    ggsave(file.path(fig_dir, sprintf('%s_bortezomib_treat_cluster_GSEA.png', cur_CL)),
           width = 5.5, height = 2.5)
  }
  
  for (cur_CL in targ_CLs) {
    make_subpop_plots(cur_CL)
  }
}