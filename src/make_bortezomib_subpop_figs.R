make_bortezomib_subpop_figs <- function(targ_CLs, dred = 'umap') {
  library(magrittr)
  library(edgeR)
  
  #PARAMS
  n_PCs <- 10
  use_SCTransform <- FALSE
  clust_k <- 20
  clust_res <- 0.25
  #umap params
  umap_n_neighbors <- 15
  umap_min_dist <- 1
  
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
      seuSub <- FindVariableFeatures(object = seuSub,
                                  nfeatures = globals$n_highvar_genes,
                                  do.plot = FALSE,
                                  selection.method = 'vst')
    }
    
    seuSub <- RunPCA(object = seuSub,
                     features = VariableFeatures(seuSub),
                     npcs = n_PCs,
                     do.print = FALSE, 
                     verbose = F)
    
    if (dred == 'tsne') {
      seuSub <- RunTSNE(object = seuSub, dims = 1:n_PCs, perplexity = globals$tsne_perplexity, seed.use = 1)
    } else if (dred == 'umap') {
      seuSub <- RunUMAP(object = seuSub, dims = 1:n_PCs, min.dist = umap_min_dist, n.neighbors = umap_n_neighbors, seed.use = 1)
    } else {
      stop('dred value not accepted')
    }
        
    controlSub <- subset(seuSub, subset = condition %in% c('control_1', 'control_2'))
    treatSub <- subset(seuSub, subset = condition == 'treat_1')
    
    #apply clustering within treated cells
    treatSub <- FindNeighbors(treatSub, dims = 1:n_PCs, k.param = clust_k, verbose = FALSE)  
    treatSub <- Seurat::FindClusters(treatSub, resolution = clust_res, verbose = FALSE)
    
    df <- Embeddings(object = seuSub, reduction = dred) %>% 
      as.data.frame() %>% 
      .[, c(1,2)] %>% 
      set_colnames(c('t1', 't2')) %>% 
      rownames_to_column(var = 'barcode') %>% 
      cbind(seuSub@meta.data) %>% 
      dplyr::mutate(condition = ifelse(condition %in% c('control_1', 'control_2'), 'DMSO', condition),
             condition = plyr::revalue(condition, c(treat_1 = 'bortezomib')),
             condition = factor(condition, levels = c('DMSO', 'bortezomib'))) %>% 
      left_join(data_frame(barcode = names(Idents(treatSub)),
                           cluster = Idents(treatSub)), by = 'barcode') %>% 
      dplyr::mutate(cluster = plyr::revalue(cluster, replace = c(`0` = 1, `1` = 2)),
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
    df %<>% dplyr::mutate(cluster = factor(cluster, levels = c(1, 2)),
                        cluster_cond = ifelse(is.na(cluster), 'DMSO', paste0('bortezomib_', cluster)),
                    cluster_cond = factor(cluster_cond, levels = c('DMSO', 'bortezomib_1', 'bortezomib_2')))
    
    xr <- with(df, range(df$t1))
    yr <- with(df, range(df$t2))
    avgs <- df %>% 
      # dplyr::filter(condition == 'Bort. 24hr') %>% 
      dplyr::group_by(cluster_cond) %>% 
      dplyr::summarise(t1 = median(t1), t2 = median(t2)) 
    
    seg_avgs <- rbind(
      cbind(avgs %>% dplyr::filter(cluster_cond == 'DMSO') %>% dplyr::rename(c1 = t1, c2 = t2), avgs %>% filter(cluster_cond == 'bortezomib_1') %>% dplyr::select(-cluster_cond)),
      cbind(avgs %>% dplyr::filter(cluster_cond == 'DMSO') %>% dplyr::rename(c1 = t1, c2 = t2), avgs %>% filter(cluster_cond == 'bortezomib_2') %>% dplyr::select(-cluster_cond)
    ))
      
    
    g_phase <- ggplot(df, aes(t1, t2)) +
      geom_point(aes(fill = Phase, size = condition), pch = 21, alpha = 0.8, stroke = 0.25) +
      guides(
             size = guide_legend(nrow = 3, title = element_blank()), 
             fill = guide_legend(nrow = 3, override.aes = list(stroke = 0, size =3))) +
      scale_size_manual(values = c(1.25, 2.75)) +
      xlab(paste0(dred, ' 1')) + ylab(paste0(dred, ' 2')) +
      geom_text(data = avgs, aes(x = t1-2, y = t2 + 5, label = cluster_cond), size = 6) +
      geom_segment(data = seg_avgs, aes(x = c1, y = c2, xend = t1, yend = t2),
                   arrow = arrow(length = unit(0.8,"cm")), size = 1.5) +
      coord_cartesian(xlim = xr, ylim = yr, clip = 'off') +
      scale_fill_Publication() +
      theme_Publication()

    ggsave(file.path(fig_dir, sprintf('%s_%s_bortezomib_subpop_phase.png', cur_CL, dred)), width = 3.5, height = 4, plot = g_phase)
    ggsave(file.path(fig_dir, sprintf('%s_%s_bortezomib_subpop_phase.pdf', cur_CL, dred)), width = 3.5, height = 4, plot = g_phase)

    
    #Run DE analysis comparing the two treatment clusters
    used_genes <- which(Matrix::rowSums(GetAssayData(treatSub, 'counts')[, rownames(treatSub@meta.data)] > 0) > nrow(treatSub@meta.data)*globals$min_frac_cells_det)
    dge <- DGEList(GetAssayData(treatSub, 'counts')[used_genes,rownames(treatSub@meta.data)])
    dge <- calcNormFactors(dge, method = 'TMMwzp')
    design <- model.matrix(~cluster + cell_det_rate , data = df)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
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
      theme_Publication() +
      scale_fill_manual(values = c(`FALSE` = 'black', `TRUE` = 'darkred'))
    ggsave(file.path(fig_dir, sprintf('%s_bortezomib_treat_cluster_volcano.png', cur_CL)),
           width = 3.5, height = 3)
    
    gene_stat <- tt_res$logFC %>% set_names(tt_res$Gene)
    gsea_plot <- make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = 50, n_lab_per = 5, lab_size = 2.5, max_chars = 30)  
    ggsave(file.path(fig_dir, sprintf('%s_bortezomib_treat_cluster_GSEA.png', cur_CL)),
           width = 5.5, height = 2.5)
    
    #Run DE analysis comparing the treatment responses of each cluster to control
    comb_data <- merge(controlSub, treatSub)
    comb_data[['tcond']] <- ifelse(grepl('treat', comb_data[['condition']]$condition), 'treat', 'control')
    comb_data[['cond_cluster']] <- paste0(comb_data[['tcond']]$tcond, '_', comb_data[['seurat_clusters']]$seurat_clusters)
    
    used_genes <- which(Matrix::rowSums(GetAssayData(comb_data, 'counts')[, rownames(comb_data@meta.data)] > 0) > nrow(comb_data@meta.data)*globals$min_frac_cells_det)
    dge <- DGEList(GetAssayData(comb_data, 'counts')[used_genes,rownames(comb_data@meta.data)])
    dge <- calcNormFactors(dge, method = 'TMMwzp')
    design <- model.matrix(~0 + cond_cluster + cell_det_rate , data = comb_data@meta.data)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    
    treat_group1 <- (grepl('treat_0', colnames(design))) %>% magrittr::divide_by(., sum(.))
    treat_group2 <- (grepl('treat_1', colnames(design))) %>% magrittr::divide_by(., sum(.))  
    control_group <- (grepl('control', colnames(design))) %>% magrittr::divide_by(., sum(.))   
    
    #get differential treatment effect
    contr_vec <- treat_group2-control_group - (treat_group1 -control_group)
    qlf <- glmQLFTest(fit, contrast = contr_vec)
    tt_diff <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
    
    #cluster1 treatment effect
    contr_vec <- treat_group1-control_group 
    qlf <- glmQLFTest(fit, contrast = contr_vec)
    tt_CL1 <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
    
    #cluster2 treatment effect
    contr_vec <- treat_group2-control_group
    qlf <- glmQLFTest(fit, contrast = contr_vec)
    tt_CL2 <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
    
    tt_comb <- full_join(tt_CL1, tt_CL2, by = 'Gene', suffix = c('_CL1', '_CL2')) %>% 
      left_join(tt_diff %>% 
                  dplyr::select(Gene, PValue_diff = PValue, logFC_diff = logFC, FDR_diff = FDR), by = "Gene")
    
    df <- tt_comb %>% dplyr::mutate(FDR = FDR_diff)
    ggplot(df, 
           aes(logFC_CL1, logFC_CL2)) + 
      xlab('logFC cluster 1') +
      ylab('logFC cluster 2') +
      geom_point(aes(fill = -log10(FDR)), pch = 21, alpha = 0.8, size = 2, stroke = 0.1, color = 'white') + 
      geom_abline(linetype = 'dashed') + 
      geom_hline(yintercept = 0, linetype = 'dashed') + 
      geom_vline(xintercept = 0, linetype = 'dashed') +
      ggrepel::geom_label_repel(data = df %>% 
                                  dplyr::filter(FDR_diff < 0.1) %>% 
                                  arrange(dplyr::desc(abs(logFC_diff))) %>% 
                                  head(15) %>% 
                                  rbind(
                                    df %>% 
                                      dplyr::arrange(dplyr::desc(abs(logFC_CL1 + logFC_CL2))) %>% 
                                      head(5)) %>% 
                                  dplyr::distinct(Gene, .keep_all=T),
                                aes(label = Gene, color = -log10(FDR)), 
                                size = 2.5, label.padding = 0.1) +
      scale_fill_gradient(limits = c(0, 3), oob = scales::squish, breaks = c(0, 3), labels = c('0', '>3'), low = 'darkgrey', high = 'red') +
      scale_color_gradient(limits = c(0, 3), oob = scales::squish, breaks = c(0, 3), labels = c('0', '>3'), low = 'darkgrey', high = 'red') +
      guides(color = FALSE) +
      theme_Publication() 
    ggsave(file.path(fig_dir, sprintf('%s_bortezomib_treat_cluster_LFC_scatter.png', cur_CL)),
           width = 3.5, height = 3.5)
    
  }
  
  for (cur_CL in targ_CLs) {
    make_subpop_plots(cur_CL)
  }
}
