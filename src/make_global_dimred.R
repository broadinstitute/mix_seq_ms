
make_global_dimred <- function() {
  #PARAMS
  n_pca_genes <- 5000
  min_cells_per_cond <- 5
  min_tot_cells <- 40
  min_per_cond <- 10
  npcs <- 25
  sens_clims <- c(0, 0.5)
  metric <- 'cosine'
  umap_nneighbors <- 15
  umap_mindist <- 0.7
  
  used_expts <- setdiff(names(sc_DE_meta), c('GPX4_expt2', 
                                             'navitoclax_24hr_expt3',
                                             'AZD5591_expt10',
                                             'everolimus_6hr_expt3', 
                                             'Trametinib_expt10'))
  all_res <- llply(sc_DE_meta[used_expts], function(cur_expt) {
    print(sprintf('Processing expt %s', cur_expt$expt_name))
    out_dir <- file.path(results_dir, cur_expt$expt_name)
    
    drug_sens <- all_CL_features[[cur_expt$expt_name]] %>% 
      dplyr::select(DEPMAP_ID, sens)
    
    dat <- load_compute_CPM_stats(cur_expt, results_dir, prior_counts = globals$pca_prior_cnt, type = 'avg_cpm')
    cell_df <- load_cell_info(cur_expt)
    cell_counts <- cell_df %>% 
      dplyr::mutate(condition = ifelse(grepl('treat', condition), 'treat', 'control')) %>% 
      dplyr::group_by(singlet_ID, condition) %>% 
      dplyr::summarise(n = n()) %>% 
      reshape2::dcast(singlet_ID ~ condition, value.var = 'n') %>% 
      mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(singlet_ID))
    u_cell_lines <- cell_counts %>% 
      dplyr::filter(singlet_ID %in% colnames(dat$LFC_mat)) %>% 
      dplyr::filter(treat >= min_cells_per_cond, control >= min_cells_per_cond) %>% 
      .[['singlet_ID']]
    return(list(
      LFCs = dat$LFC_mat[, u_cell_lines],
      sample_df = data_frame(
        drug = cur_expt$drug_name,
        expt = cur_expt$expt_name,
        CCLE_ID = colnames(dat$LFC_mat[, u_cell_lines])
      ) %>% 
        dplyr::mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(CCLE_ID)) %>% 
        left_join(drug_sens, by = 'DEPMAP_ID') %>% 
        left_join(cell_counts, by = c('DEPMAP_ID'))
    ))
  })
  
  #make combined matrices 
  all_LFCs <- do.call(cbind, lapply(all_res, function(x) x$LFCs))
  all_info <- do.call(rbind, lapply(all_res, function(x) x$sample_df)) %>% 
    dplyr::mutate(tot_cells = treat + control)
  
  usamps <- which(all_info[,'tot_cells'] >= min_tot_cells & all_info[, 'control'] >= min_per_cond & all_info[, 'treat'] >= min_per_cond)
  all_info <- all_info[usamps,]
  all_LFCs <- all_LFCs[, usamps]
  
  #get genes for PCA
  comp_genes <- all_LFCs %>% 
    apply(1, sd, na.rm=T) %>% 
    sort(decreasing = T) %>% 
    head(n_pca_genes) %>% 
    names() %>% 
    as.character()
  
  #bundle object
  xx <- all_LFCs[comp_genes,]
  colnames(xx) <- paste0(all_info$expt, '_', all_info$CCLE_ID)
  rownames(all_info) <- colnames(xx)
  new_obj <- Seurat::CreateSeuratObject(xx,
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = all_info)
  
  new_obj <- Seurat::ScaleData(new_obj, features = rownames(GetAssayData(new_obj)), do.scale = F, do.center = T, verbose = F)
  
  new_obj %<>% Seurat::RunPCA(assay='RNA',
                              features = rownames(GetAssayData(new_obj)),
                              npcs = npcs, verbose = F)
  new_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:npcs,
                               reduction = 'pca',
                               n.neighbors = umap_nneighbors,
                               min.dist =  umap_mindist, 
                               metric = metric,
                               seed.use = 3,
                               verbose=F)
  
  df <- Embeddings(new_obj, reduction = 'umap') %>% 
    cbind(all_info) %>% 
    dplyr::mutate(expt = str_replace_all(expt, 'expt10', '24hr_expt10'),
           time = str_match(expt, '^[:alnum:]+_([:alnum:]+)')[,2],
           drug = str_match(expt, '^([:alnum:]+)_')[,2]) %>% 
    dplyr::mutate(TP53_status = ifelse(CCLE_ID %in% globals$TP53_WT_cls_pool22, 'TP53WT', 'TP53null'),
             drug = str_replace_all(drug, 'Idasanutlin', 'Nutlin'),
             drug_status = ifelse(drug == 'Nutlin', paste0(drug, '_', TP53_status), drug),
             drug_time = paste0(drug, '_', time),
             drug_status_time = paste0(drug_status, '_', time))
  
  avgs <- df %>% 
    dplyr::group_by(drug_status_time, drug_time) %>% 
    dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm=T),
              UMAP_2 = median(UMAP_2, na.rm=T))
  
  cols <- c(Afatinib_24hr = '#C77CFF',#C77CFF
            AZD5591_24hr = '#E68613',
            Bortezomib_24hr = '#CD9600',
            Bortezomib_6hr = '#ABA300',
            BRD3379_24hr = '#7CAE00',
            BRD3379_6hr = '#0CB702',
            Dabrafenib_24hr = '#00BE67',
            Everolimus_24hr = '#00C19A',
            Gemcitabine_24hr = '#FF68A1',
            Navitoclax_24hr = 'darkgreen',
            JQ1_24hr = '#00B8E7',
            Nutlin_24hr = '#00A9FF',
            Nutlin_6hr = '#8494FF',
            Prexasertib_24hr = '#F8766D', 
            Taselisib_24hr = '#ED68ED',
            Trametinib_24hr = '#FF61CC',
            Trametinib_6hr = '#00BFC4')
  stopifnot(all(df$drug_time %in% names(cols)))
  
  zoom_x <- c(-7.75, -1.5); zoom_y <- c(-7, 0)
  
  #Make overall plot colored by treatment
  g <- df %>% 
    ggplot(aes(UMAP_1, UMAP_2)) + 
    geom_point(aes(fill = drug_time, text = sprintf('%s\n%s\n%.3f', drug_time, CCLE_ID, sens)), color = 'white', stroke = 0.1, pch = 21, size = 1.8) +
    guides(fill = F, color = F) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) + 
    geom_label_repel(data = avgs, aes(label = drug_status_time, color = drug_time), size = 2.5, label.padding = 0.1) +
    theme_Publication() +
    geom_rect(xmin = zoom_x[1], xmax = zoom_x[2], ymin = zoom_y[1], ymax = zoom_y[2], 
              fill = NA, color = 'black', lwd = 0.2)
  ggsave(file.path(fig_dir, 'full_LFC_umap_treat.png'), plot = g, width = 5, height = 4)
  # plotly::ggplotly(g)
  
  #now make zoomed plot
  omics_BRAF <- all_CL_features$Dabrafenib_24hr_expt3
  avgs_z <- df %>% 
    dplyr::filter(UMAP_1 > zoom_x[1], UMAP_1 < zoom_x[2],
           UMAP_2 > zoom_y[1], UMAP_2 < zoom_y[2],
           drug_time %in% c('Trametinib_6hr', 'Trametinib_24hr', 'Afatinib_24hr')) %>% 
    dplyr::group_by(drug_time) %>% 
    dplyr::summarise(UMAP_1 = median(UMAP_1),
              UMAP_2 = median(UMAP_2)) 
  g <- df %>% 
    left_join(omics_BRAF, by = 'CCLE_ID') %>% 
    mutate(sens_BRAF = drug == 'Dabrafenib' & BRAF_MUT > 0) %>% 
    ggplot(aes(UMAP_1, UMAP_2)) + 
    geom_point(aes(fill = drug_time), color = 'white', stroke = 0.1, pch = 21, size = 2) +
    geom_point(data = . %>% filter(sens_BRAF),
               aes(fill = drug_time), color = 'white', stroke = 0.1, pch = 21, size = 4) +
    guides(fill = F, color = F) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    # geom_text(data = avgs_z, aes(label = drug_time, color = drug_time), size = 4) +
    theme_Publication()  +
    theme_void() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank()) +
    # geom_label_repel(data = avgs, aes(label = drug_status_time, color = drug_time), size = 2.25, label.padding = 0.1) +
    coord_cartesian(xlim = zoom_x, ylim = zoom_y)
  ggsave(file.path(fig_dir, 'full_LFC_umap_treat_zoom.png'), plot = g, width = 2.5, height = 2)
  
  #now color by sensitivity
  g <- df %>% 
    ggplot(aes(UMAP_1, UMAP_2)) + 
    geom_point(aes(fill = sens, text = sprintf('%s\n%s\n%.3f', drug_time, CCLE_ID, sens)), color = 'white', stroke = 0.1, pch = 21, size = 1.8) +
    guides(color = F, fill = guide_colorbar(title = 'sensitivity')) +
    # geom_label_repel(data = avgs, aes(label = drug_status_time, color = drug_time), size = 2.25, label.padding = 0.1) +
    theme_Publication()  +
    theme(legend.key.width = unit(1.25, 'cm')) +
    scale_color_manual(values = cols) +
    scale_fill_gradient(limits = sens_clims, oob = scales::squish) 
  ggsave(file.path(fig_dir, 'full_LFC_umap_sens.png'), plot = g, width = 4, height = 4)

}
