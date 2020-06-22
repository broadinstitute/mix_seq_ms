make_bortezomib_heterogeneity_figs <- function() {
  library(magrittr)
  library(reshape2)
  library(ComplexHeatmap)
  library(Seurat)
  library(tidyverse)

  #PARAMS
  treatment <- "bortezomib_24hr_expt1"
  control <- "dmso_24hr_expt1"
  min_cells_per_cond <- 10
  cluster_res <- .25

  #LOAD DATA
  sc_data <- load_sc_data(sc_DE_meta$bortezomib_24hr_expt1, sc_expts = sc_expts)
  #sc_data <- cdsc::load_and_process_sc_data_by_expt(c(treatment,control), quality = c("normal"))

  n_cls <- sc_data$singlet_ID %>% unique() %>% length()
  n_condition <- sc_data@meta.data %>% count(singlet_ID,condition) %>% acast(singlet_ID ~ condition, value.var = 'n')
  usable_cls <- rownames(n_condition)[rowSums(n_condition >= min_cells_per_cond) == 2]
  sc_data <- sc_data[,sc_data$singlet_ID %in% usable_cls]
  n_PCs <- length(usable_cls) * globals$n_pcs_per_CL

  #PROCESS DATA
  sc_data  <- NormalizeData(object = sc_data,
                          normalization.method = "LogNormalize",
                          scale.factor = globals$scale_fac)
  sc_data <- CellCycleScoring(object = sc_data,
                             s.features = Seurat::cc.genes$s.genes,
                             g2m.features = Seurat::cc.genes$g2m.genes,
                             set.ident = FALSE)
  sc_data <- relabel_cell_cycle_phase(sc_data)
  sc_data <- ScaleData(object = sc_data)
  sc_data <- FindVariableFeatures(object = sc_data,
                              top.genes = globals$n_highvar_genes,
                              do.plot = FALSE,
                              selection.method = 'vst')
  sc_data <- RunPCA(object = sc_data,
                   features = VariableFeatures(sc_data),
                   npcs = n_PCs,
                   seed.use = 1,
                   do.print = FALSE,
                   verbose = F)

  #IDENTIFY CLUSTERS
  sc_cluster <- NULL
  for (cl in unique(sc_data$singlet_ID)){
    sc_cl <- sc_data[,sc_data$singlet_ID == cl & sc_data$condition == treatment]
    sc_cl %<>% FindNeighbors(dims = 1:n_PCs,k.param = 28)
    sc_cl %<>% FindClusters(resolution = .25)
    if (is.null(sc_cluster)){
      sc_cluster <- sc_cl
    }
    else {
      sc_cluster %<>% merge(sc_cl)
    }
  }
  sc_data<- sc_data[,sc_data$condition == control] %>% merge(sc_cluster)

  #ORDER CLUSTERS BY CELL CYCLE
  hetro_cls <- sc_data[,sc_data$condition == treatment]@meta.data %>% group_by(singlet_ID) %>%
    summarise(n_cluster = n_distinct(seurat_clusters)) %>% filter(n_cluster > 1) %>% .[["singlet_ID"]]
  sc_hetro <- sc_data[,sc_data$condition == treatment & sc_data$singlet_ID %in% hetro_cls]
  ordered_clusters <- sc_hetro@meta.data %>% count(singlet_ID,Phase,seurat_clusters) %>% spread(value = "n",key = "Phase",fill = 0) %>%
    mutate(frac_S = S/(S+`G0/G1`+`G2/M`)) %>% group_by(singlet_ID) %>%
    mutate(flip = ifelse(first(frac_S)>last(frac_S),TRUE,FALSE)) %>% ungroup() %>%
    mutate(new_cluster = ifelse(flip,ifelse(seurat_clusters == 0,1,0),ifelse(seurat_clusters == 0,0,1)))
  md <-  sc_hetro@meta.data %>% left_join(ordered_clusters)
  sc_hetro %<>% AddMetaData(as.character(md$new_cluster),col.name = "cluster")


  #DIFFERENTIAL EXPRESSION
  md <- sc_hetro@meta.data %>% rownames_to_column(var = 'cell_ID') %>% group_by(singlet_ID,cluster)
  group_sum <- function(group) {
    Matrix::rowSums(GetAssayData(sc_hetro, slot = 'counts')[, group$cell_ID, drop = FALSE])
  }
  summed_counts <- md %>% group_split() %>% plyr::laply(.fun = group_sum)
  sample_info <- md %>% dplyr::summarise() %>% transmute(CCLE_ID = singlet_ID,expt = ifelse(cluster == "1","treat","control"),
                                          condition = expt) %>% as.data.frame()
  rownames(summed_counts) <- md %>% dplyr::summarise() %>% tidyr::unite(ID,singlet_ID,cluster,sep = ',') %>% .[["ID"]]
  limma_res <- fit_viability_models(summed_counts %>% t(),sample_info = sample_info)
  limma_avg <- limma_res$res_avg
  hits <- limma_avg %>% filter(P.Value < .05) %>% arrange(abs(logFC)) %>% tail(n_genes_to_show)


  #VOLCANO PLOT
  volcano <- limma_avg %>% ggplot(aes(logFC, -log10(P.Value))) +
    geom_point(data = limma_avg %>% filter(!Gene %in% Seurat::cc.genes$s.genes), fill = 'black', pch = 21, size = 1.5, alpha = 0.8, color = 'gray', stroke = 0.1) +
    geom_point(data = limma_avg %>% filter(Gene %in% Seurat::cc.genes$s.genes), fill = 'darkred', pch = 21, size = 3, alpha = 0.8, color = 'gray', stroke = 0.1) +
    geom_label_repel(data = . %>% dplyr::arrange(dplyr::desc(logFC)) %>% head(10),
                     aes(label = Gene),
                     size = 2.5, label.padding = 0.2) +
    geom_label_repel(data = . %>% dplyr::arrange(logFC) %>% head(10),
                     aes(label = Gene),
                     size = 2.5, label.padding = 0.2) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    cdsr::theme_Publication()

  ggsave(file.path(fig_dir,'bortezomib_volcano.png'), width = 4, height = 4, plot = volcano)

  #HEATMAP
  heatmap_genes <- limma_avg %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% head(25) %>%  .[["Gene"]]
  heatmap_genes <- c(heatmap_genes,limma_avg %>% filter(logFC < 0) %>% arrange(adj.P.Val) %>% head(25) %>%  .[["Gene"]])

  mean_control <- get_summed_counts(sc_data[,sc_data$condition == control]) %>% t() %>% cpm(log = TRUE,prior.count = 1)
  treatment_counts <- GetAssayData(sc_data[,sc_data$condition == treatment],slot = "counts") %>%
    cpm(log = TRUE,prior.count = .05)
  logFC <- NULL
  for (cl in unique(sc_data$singlet_ID)){
    cells = sc_data@meta.data %>% rownames_to_column(var = "barcode") %>%
      filter(condition == treatment,singlet_ID == cl) %>% .[["barcode"]]
    if(is.null(logFC)){
      logFC <- treatment_counts[heatmap_genes,cells] -  mean_control[heatmap_genes,cl]
    }
    else{
      logFC %<>% cbind(treatment_counts[heatmap_genes,cells] -  mean_control[heatmap_genes,cl])
    }
  }

  sc_hetro %<>% AddMetaData(str_c(sc_hetro$singlet_ID,sc_hetro$cluster),col.name = "cl_cluster")

  md <- sc_hetro[,sc_hetro$condition == treatment]@meta.data
  md %<>% rownames_to_column(var = "barcode")
  md <- sc_hetro@meta.data %>% rownames_to_column(var = "barcode") %>% select(barcode,cluster) %>%
    mutate(fixed_cluster = cluster) %>% select(-cluster) %>% right_join(md,by = "barcode") %>%
    mutate(replace_na(fixed_cluster,0))
  df_ana <- md %>% mutate(`Cell line` = str_split_fixed(singlet_ID,pattern = "_",n = 2)[,1],
                          Cluster = ifelse(fixed_cluster == 0,1,2)) %>%  select(`Cell line`,Cluster,Phase)
  row_split <- factor(md$cl_cluster,levels= c("RCC10RGB_KIDNEY0","RCC10RGB_KIDNEY1","SNU1079_BILIARY_TRACT0",
                                              "SNU1079_BILIARY_TRACT1","TEN_ENDOMETRIUM0","TEN_ENDOMETRIUM1",
                                              "RERFLCAD1_LUNG0","RERFLCAD1_LUNG1","NCIH226_LUNG0","NCIH226_LUNG1",
                                              "SKMEL3_SKIN0","SKMEL3_SKIN1","COLO680N_OESOPHAGUS0","COLO680N_OESOPHAGUS1",
                                              "IALM_LUNG0","IALM_LUNG1","RCM1_LARGE_INTESTINE0","RCM1_LARGE_INTESTINE1"))
  col_ano <- rowAnnotation(df = df_ana,col = list(Phase = c("G0/G1" = "#386cb0", "G2/M" = "#7fc97f", "S" = "#fdb462"),
                                                  Cluster = c(`1` = "gray", `2` = "brown")))
  pdf(file.path(fig_dir,'bortezomib_heatmap.pdf'), width = 8, height = 7)
  Heatmap(t(logFC[,md$barcode]),name = "logFC",cluster_rows = F,cluster_row_slices = T,show_row_names = F,
          left_annotation = col_ano,split = row_split,column_names_gp = gpar(fontsize = 10), row_title = " ",
          show_column_dend = F,row_gap = unit(c(1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3), "mm"),border = TRUE)
  dev.off()
  return()
}
