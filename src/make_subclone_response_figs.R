make_subclone_response_figs <- function(expt, targ_CL) {
  library(edgeR)
  
  #PARAMS
  n_pcs_ind <- 10
  clust_k <- 20
  clust_res <- 0.25
  use_SCTransform <- TRUE
  
  seuObj <- load_sc_data(expt, sc_expts = sc_expts)
  cell_df <- seuObj@meta.data %>% 
    rownames_to_column(var = 'barcode') %>% 
    tidyr::separate(col = 'condition', into = c('treat_cond', 'batch'), sep = '_')
  seuObj[['treat_cond']] <- cell_df$treat_cond
  
  #process data to align for sub-pop clustering
  seuObj <- NormalizeData(object = seuObj, 
                          normalization.method = "LogNormalize", 
                          scale.factor = globals$scale_fac)
  seuObj <- CellCycleScoring(object = seuObj,
                             s.features = Seurat::cc.genes$s.genes,
                             g2m.features = Seurat::cc.genes$g2m.genes,
                             set.ident = FALSE)
  seuObj <- relabel_cell_cycle_phase(seuObj)
  
  si <- FetchData(seuObj, 'singlet_ID')
  seuSub <- seuObj[, si$singlet_ID == targ_CL]
 
    if (use_SCTransform) {
      seuSub <- SCTransform(seuSub, variable.features.n = globals$n_highvar_genes)
      seuSub_al <- SCTransform(seuSub, vars.to.regress = c('condition'), variable.features.n = globals$n_highvar_genes)
    } else {
      seuSub <- ScaleData(object = seuSub)
      seuSub <- FindVariableFeatures(object = seuSub,
                                     nfeatures = globals$n_highvar_genes,
                                     do.plot = FALSE,
                                     selection.method = 'vst')
      seuSub_al <- ScaleData(object = seuSub, vars.to.regress = c('condition'))
      seuSub_al <- FindVariableFeatures(object = seuSub_al,
                                     nfeatures = globals$n_highvar_genes,
                                     do.plot = FALSE,
                                     selection.method = 'vst')
    }
    seuSub <- RunPCA(object = seuSub,
                     npcs = n_pcs_ind,
                     features = VariableFeatures(seuSub),
                     do.print = FALSE,
                     verbose = FALSE,
                     seed.use = 42)
    seuSub_al <- RunPCA(object = seuSub_al,
                        features = VariableFeatures(seuSub_al),
                     npcs = n_pcs_ind,
                     do.print = FALSE,
                     verbose = FALSE,
                     seed.use = 42)
    seuSub <- RunTSNE(object = seuSub, dims = 1:n_pcs_ind, 
                      check_duplicates = FALSE,
                      perplexity = globals$tsne_perplexity,
                      seed.use = 42)
    seuSub_al <- Seurat::FindNeighbors(seuSub_al, 
                                   dims = 1:n_pcs_ind,
                                   k.param = clust_k,
                                   verbose = FALSE)
    seuSub_al <- Seurat::FindClusters(seuSub_al, 
                                      resolution = clust_res,
                                      random.seed = 1,
                                      verbose = FALSE)
    
    seuSub[['cluster']] <- Idents(seuSub_al)
    stopifnot(nlevels(seuSub_al) == 2) #code assumes we find two clusters
    
    df <- Embeddings(seuSub, reduction = 'tsne') %>% 
      cbind(seuSub@meta.data)
    df %<>% dplyr::mutate(tcond = plyr::revalue(treat_cond, c(`control` = 'DMSO', `treat` = expt$drug_name)),
                          cluster = plyr::revalue(cluster, replace = c(`0` = 1, `1` = 2)),
                          Phase = factor(Phase, levels = c('G0/G1', 'S', 'G2/M')))
    
    avgs <- df %>% 
      dplyr::filter(tcond == 'DMSO') %>% 
      dplyr::group_by(cluster) %>% 
      dplyr::summarise(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) %>% 
      dplyr::mutate(cluster_label = paste0('cluster ', cluster))
    tavgs <- df %>% 
      dplyr::filter(tcond != 'DMSO') %>% 
      dplyr::group_by(cluster) %>% 
      dplyr::summarise(tSNE_t1 = median(tSNE_1), tSNE_t2 = median(tSNE_2)) %>% 
      dplyr::mutate(cluster_label = paste0('cluster ', cluster))
    seg_avgs <- full_join(avgs, tavgs, by = 'cluster')
    
    ggplot(df, aes(tSNE_1, tSNE_2)) +
      geom_point(aes(fill = tcond), size = 2, pch = 21, color = 'white', stroke = 0.1, alpha = 0.75) +
      xlab('tSNE 1') + ylab('tSNE 2') +
      scale_fill_manual(values = c('darkgray', 'indianred4')) +
      cdsr::theme_Publication() +
      geom_text(data = avgs, aes(x = tSNE_1, y = tSNE_2, label = cluster_label), size = 7) +
      geom_segment(data = seg_avgs, aes(x = tSNE_1, y = tSNE_2, xend = tSNE_t1, yend = tSNE_t2),
                   arrow = arrow(length = unit(0.2,"cm")), size = 0.75) +
      coord_cartesian(clip = 'off') +
      guides(fill = guide_legend(title = element_blank(), override.aes = list(size = 3)))
    ggsave(file.path(fig_dir, sprintf('%s_%s_cluster.png', targ_CL, expt$expt_name)), width = 4, height = 4)

    ggplot(df, aes(tSNE_1, tSNE_2)) +
      geom_point(aes(color = tcond, fill = Phase, stroke = tcond), pch = 21, size = 2.25, alpha = 0.75) +
      xlab('tSNE 1') + ylab('tSNE 2') +
      cdsr::scale_fill_Publication() +
      scale_color_manual(values = c('darkgray', 'indianred4')) +
      scale_discrete_manual(aesthetics = "stroke", values = c(0.5, 1)) +
      cdsr::theme_Publication() +
      geom_text(data = avgs, aes(x = tSNE_1, y = tSNE_2, label = cluster_label), size = 7) +
      geom_segment(data = seg_avgs, aes(x = tSNE_1, y = tSNE_2, xend = tSNE_t1, yend = tSNE_t2),
                   arrow = arrow(length = unit(0.2,"cm")), size = 0.75) +
      coord_cartesian(clip = 'off') +
      guides(stroke = FALSE,
             fill = guide_legend(nrow = 3, override.aes = list(stroke = 0, size = 3)), 
             color = guide_legend(title = element_blank(), override.aes = list(size = 3), nrow = 2))
    ggsave(file.path(fig_dir, sprintf('%s_%s_cluster_phase.png', targ_CL, expt$expt_name)), width = 3.5, height = 4)

    
    #RUN DE analysis to test cluster-treatment interaction
    seuSub[['cluster_condition']] <- paste0(seuSub[['condition']][,1], '.', seuSub[['cluster']][,1])
    
    used_genes <- which(Matrix::rowSums(GetAssayData(seuSub, 'counts')[, rownames(seuSub@meta.data)] > 0) > nrow(seuSub@meta.data)*globals$min_frac_cells_det)
    dge <- DGEList(GetAssayData(seuSub, 'counts')[used_genes,rownames(seuSub@meta.data)])
    dge <- calcNormFactors(dge, method = 'TMMwzp')
    
    #build design matrix for treat-vs-control comparison (avg across cell lines)
    design <- model.matrix(~0 + cell_det_rate + cluster_condition, data = seuSub@meta.data)
    colnames(design) <- str_replace(colnames(design), 'cluster_condition', '')
    
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    
    treat_group1 <- as.numeric(grepl('treat', colnames(design)) & grepl('0', colnames(design))) %>% magrittr::divide_by(., sum(.))
    control_group1 <- as.numeric(grepl('control', colnames(design)) & grepl('0', colnames(design))) %>% magrittr::divide_by(., sum(.))
    treat_group2 <- as.numeric(grepl('treat', colnames(design)) & grepl('1', colnames(design))) %>% magrittr::divide_by(., sum(.))
    control_group2 <- as.numeric(grepl('control', colnames(design)) & grepl('1', colnames(design))) %>% magrittr::divide_by(., sum(.))
    
    #get differential treatment effect
    contr_vec <- treat_group2-control_group2 - (treat_group1 -control_group1)
    qlf <- glmQLFTest(fit, contrast = contr_vec)
    tt_diff <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
    
    #cluster1 treatment effect
    contr_vec <- treat_group1-control_group1 
    qlf <- glmQLFTest(fit, contrast = contr_vec)
    tt_CL1 <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
    
    #cluster2 treatment effect
    contr_vec <- treat_group2-control_group2
    qlf <- glmQLFTest(fit, contrast = contr_vec)
    tt_CL2 <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
    
    #cluster2 vs cluster1 baseline
    contr_vec <- control_group2-control_group1
    qlf <- glmQLFTest(fit, contrast = contr_vec)
    tt_CL2_CL1 <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
      
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
      guides(color = F) +
      cdsr::theme_Publication() 
    ggsave(file.path(fig_dir, sprintf('%s_%s_cluster_diffLFC.png', targ_CL, expt$expt_name)), width = 3, height = 3)
    
    tt_CL2_CL1 %>% 
      ggplot(aes(logFC, -log10(PValue))) +
      geom_point(aes(color = FDR < 0.1), size = 0.75, alpha = 0.5) +
      geom_label_repel(data = . %>% dplyr::arrange(dplyr::desc(abs(logFC))) %>% head(10),
                       aes(label = Gene),
                       size = 2.5) +
      geom_vline(xintercept = 0, linetype = 'dashed') + 
      geom_hline(yintercept = 0, linetype = 'dashed') +
      cdsr::theme_Publication() +
      scale_color_manual(values = c(`FALSE` = 'black', `TRUE` = 'darkred'))
    ggsave(file.path(fig_dir, sprintf('%s_%s_baseline_cluster_volcano.png', targ_CL, expt$expt_name)),
           width = 4, height = 3.5)
    
    gene_stat <- tt_CL2_CL1$logFC %>% set_names(tt_CL2_CL1$Gene)
    make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = 50, n_lab_per = 5, lab_size = 3, max_chars = 30)  
    ggsave(file.path(fig_dir, sprintf('%s_%s_baseline_cluster_GSEA.png', targ_CL, expt$expt_name)),
           width = 5.5, height = 2.5)
    
    tt_diff %>% 
      ggplot(aes(logFC, -log10(PValue))) +
      geom_point(aes(color = FDR < 0.1), size = 0.75, alpha = 0.5) +
      geom_label_repel(data = . %>% dplyr::arrange(dplyr::desc(abs(logFC))) %>% head(10),
                       aes(label = Gene),
                       size = 2.5) +
      geom_vline(xintercept = 0, linetype = 'dashed') + 
      geom_hline(yintercept = 0, linetype = 'dashed') +
      cdsr::theme_Publication() +
      scale_color_manual(values = c(`FALSE` = 'black', `TRUE` = 'darkred'))
    ggsave(file.path(fig_dir, sprintf('%s_%s_diff_cluster_volcano.png', targ_CL, expt$expt_name)),
           width = 4, height = 3.5)
    gene_stat <- tt_diff$logFC %>% set_names(tt_diff$Gene)
    make_hyper_gsa_plot(gene_stat, gsc_data$combined, top_n = 50, n_lab_per = 5, lab_size = 3, max_chars = 30)  
    ggsave(file.path(fig_dir, sprintf('%s_%s_diff_cluster_GSEA.png', targ_CL, expt$expt_name)),
           width = 5.5, height = 2.5)
}