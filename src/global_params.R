# Params
globals <- list(
  scale_fac = 1e5,
  n_highvar_genes = 5000, #number of high-var genes to use for PCA/clustering
  nn_cluster_k = 10, #nearest-neighbor param for clustering
  vtr = c(), #vars to regress
  n_pcs_per_CL = 2, #num PCs per cell line
  tsne_perplexity = 25, 
  use_voom = FALSE,
  prior_cnt = 1,
  pca_prior_cnt = 10,
  q_thresh = 0.1,
  gsea_top_n = 50,
  min_frac_cells_det = 0.05,
  TP53_WT_cls_pool22 = c('LNCAPCLONEFGC_PROSTATE',
                   'DKMG_CENTRAL_NERVOUS_SYSTEM',
                   'NCIH226_LUNG',
                   'RCC10RGB_KIDNEY',
                   'SNU1079_BILIARY_TRACT',
                   'CCFSTTG1_CENTRAL_NERVOUS_SYSTEM',
                   'COV434_OVARY')
)
