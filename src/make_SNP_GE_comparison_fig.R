#make tSNE plot comparing SNP and GE-based classification for a given experiment batch (neg-control data)
make_SNP_GE_comparison_fig <- function(expt_batch, dred) {
  
  all_ref_lines <- read_csv(here::here('data', 'bulk_reference_CLs.csv'))
  source(here::here('src', 'global_params.R'))
  
  #PARAMS
  n_highvar_genes_compare <- 5000 #number of genes to use for comparing CCLE and cluster avgs
  vtr <- c()
  #umap params
  umap_n_neighbors <- 15
  umap_min_dist <- 1
  
  if (expt_batch == 'expt1') {
    cur_expt <- list(
      data_sets = list(Untreated = "Untreated_6hr_expt1",
                       DMSO_6hr = 'DMSO_6hr_expt1',
                       DMSO_24hr = 'DMSO_24hr_expt1'),
      expt_batch = 'expt1',
      drug_name = NULL,
      load_from_taiga = TRUE
    )
    clust_res <- 1
    dot_size <- 1
  } else if (expt_batch == 'expt3') {
    cur_expt <- list(
      data_sets = list(
        DMSO_6hr = 'DMSO_6hr_expt3',
        DMSO_24hr = 'DMSO_24hr_expt3'),
      expt_batch = 'expt3',
      drug_name = NULL,
      load_from_taiga = TRUE
    )
    clust_res <- 4
    dot_size <- 0.5
  }
  
  CCLE_GE <- load.from.taiga(data.name="depmap-rnaseq-expression-data-ccd0",
                             data.version = 14,
                             data.file = 'CCLE_depMap_19Q3_TPM') 
  
  seuObj <- load_sc_data(cur_expt, sc_expts = sc_expts, QC_filter = FALSE)
  seuObj[['DEPMAP_ID']] <- celllinemapr::ccle.to.arxspan(Seurat::FetchData(seuObj, 'singlet_ID')$singlet_ID)
  cq <- FetchData(seuObj, vars = c('cell_quality'))
  seuObj <- seuObj[, which(cq$cell_quality %in% c('normal', 'doublet'))]
  
  n_cls <- nlevels(seuObj)
  n_pcs <- globals$n_pcs_per_CL*n_cls
  
  seuObj <- NormalizeData(object = seuObj, 
                          normalization.method = "LogNormalize", 
                          scale.factor = globals$scale_fac)
  
  seuObj <- ScaleData(object = seuObj, vars.to.regress = vtr)
  
  seuObj <- FindVariableFeatures(object = seuObj,
                                 nfeatures = globals$n_highvar_genes,
                                 selection.method = 'vst')
  
  seuObj <- RunPCA(object = seuObj,
                   features = VariableFeatures(seuObj),
                   seed.use = 1,
                   npcs = n_pcs,
                   verbose = FALSE)
  
  if (dred == 'tsne') {
    seuObj <- RunTSNE(object = seuObj, dims = 1:n_pcs, check_duplicates = FALSE, seed.use = 1)
  } else if (dred == 'umap') {
    seuObj <- RunUMAP(object = seuObj, dims = 1:n_pcs, min.dist = umap_min_dist, n.neighbors = umap_n_neighbors)
  } else {
    stop('dred value not accepted')
  }
  
  full_df <- Embeddings(object = seuObj, reduction = dred) %>% 
    as.data.frame() %>% 
    set_colnames(c('t1', 't2')) %>% 
    rownames_to_column(var = 'barcode') %>% 
    cbind(seuObj@meta.data) %>% 
    mutate(condition = factor(condition))
  
  #now reprocess using just good singlets
  cq <- FetchData(seuObj, vars = c('cell_quality'))
  seuObj <- seuObj[, which(cq$cell_quality == 'normal')]
  seuObj <- ScaleData(object = seuObj, vars.to.regress = vtr)
  seuObj <- FindVariableFeatures(object = seuObj,
                                 nfeatures = globals$n_highvar_genes,
                                 selection.method = 'vst')
  seuObj <- RunPCA(object = seuObj,
                   features = VariableFeatures(seuObj),
                   seed.use = 1,
                   npcs = n_pcs,
                   verbose = FALSE)
  
  #cluster cells
  seuObj <- Seurat::FindNeighbors(seuObj, reduction = 'pca',
                                  dims = 1:n_pcs,
                                  k.param = globals$nn_cluster_k, 
                                  force.recalc = TRUE,
                                  verbose = FALSE)
  seuObj <- Seurat::FindClusters(seuObj, resolution = clust_res, verbose = FALSE)
  
  #compute summed counts across cells for each cluster, then CPM transform
  clust_cpm <- laply(levels(seuObj), function(clust) {
    Matrix::rowSums(GetAssayData(seuObj, slot = 'counts')[, names(Idents(seuObj)[Idents(seuObj) == clust])])
  }) %>%
    t() %>%
    edgeR::cpm(prior.count = globals$prior_cnt, log = TRUE)

  in_pool_lines <- unique(seuObj@meta.data$singlet_ID)
  in_pool_lines_depmap <- celllinemapr::ccle.to.arxspan(in_pool_lines)
  
  colnames(CCLE_GE) <- str_match(colnames(CCLE_GE), '\\((ENSG[0-9]+)\\)')[,2]
  sc_ensemble <- seuObj@misc$Ensembl_ID
  common_genes <- intersect(colnames(CCLE_GE), sc_ensemble)
  
  high_var_genes <- apply(CCLE_GE[in_pool_lines_depmap, match(common_genes, colnames(CCLE_GE))], 2, sd) %>% 
    sort(decreasing = T) %>% 
    head(globals$n_highvar_genes) %>% 
    names()  
  
  #mean-center both CCLE and cluster expression profiles
  CCLE_GE_norm <- CCLE_GE[in_pool_lines_depmap, high_var_genes] %>% scale(center = T, scale = F) %>% t() 
  clust_cpm_norm <- clust_cpm[match(high_var_genes, seuObj@misc$Ensembl_ID),] %>% t() %>% scale(center = T, scale = F) %>% t()
  
  cc <- cor(CCLE_GE_norm, clust_cpm_norm)
  best_matches <- in_pool_lines_depmap[apply(cc, 2, which.max)]
  GE_class <- best_matches[as.numeric(Idents(seuObj))]
  
  perc_correct <- mean(GE_class == seuObj@meta.data$DEPMAP_ID)
  
  full_df %<>% left_join(data_frame(barcode = rownames(seuObj@meta.data),
                                    GE_class = GE_class), 
                         by = 'barcode')
  full_df %<>% mutate(agree = GE_class == DEPMAP_ID)
  
  ggplot(full_df, aes(t1, t2)) +
    geom_point(data = full_df %>% filter(agree),
               aes(color = DEPMAP_ID), alpha = 0.75, size = 0.15) +
    geom_point(data = full_df %>% filter(cell_quality == 'doublet'), color = 'black', alpha = 0.75, size = dot_size) +
    geom_point(data = full_df %>% filter(cell_quality == 'low_quality'), color = 'gray', alpha = 0.75, size = dot_size) +
    geom_point(data = full_df %>% filter(!agree, cell_quality == 'normal'), color = 'red', alpha = 0.75, size = dot_size) +
    guides(color = F) +
    xlab(paste0(dred, ' 1')) + ylab(paste0(dred, ' 2')) +
    # ggtitle(sprintf('%.1f%% agreement', perc_correct*100)) +
    cdsr::theme_Publication()
  ggsave(file.path(fig_dir, sprintf('GE_SNP_agree_%s_%s.png', dred, cur_expt$expt_batch)), width = 6, height = 5.5)
  
  results = list(perc_correct = perc_correct)
  return(results)
}
