
make_hash_tag_qc_figs <- function() {
  base_data_dir <- '/Volumes/GoogleDrive/My Drive/Project Apollo/Apollo drugs /scRNAseq_005'
  
  #load cell info
  cur_expt <- sc_expts$tramet_tc_expt5
  cell_df <- load_cell_info(cur_expt, QC_filter = FALSE)
  cell_df %<>% dplyr::mutate(hash_quality = ifelse(!(hash_tag %in% c('multiplet', 'unknown')), 'singlet', hash_tag),
                             cell_quality = ifelse(!(cell_quality %in% c('normal', 'doublet')), 'low_quality', cell_quality))
  
  #load hash count data
  hash_countsA <- data.table::fread(file.path(base_data_dir, 'channelA/channelA_combined_hashing.csv'))
  hashA_df <- data.frame(t(hash_countsA[, 2:ncol(hash_countsA)]))
  colnames(hashA_df) <- hash_countsA[[1]]
  hashA_df %<>% rownames_to_column(var = 'barcode')
  
  hash_countsB <- data.table::fread(file.path(base_data_dir, 'channelB/channelB_combined_hashing.csv'))
  hashB_df <- data.frame(t(hash_countsB[, 2:ncol(hash_countsB)]))
  colnames(hashB_df) <- hash_countsB[[1]]
  hashB_df %<>% rownames_to_column(var = 'barcode')
  
  #label cell barcodes with 10x channel and combine
  hashA_df %<>% mutate(barcode = paste0(barcode, '_A'))
  hashB_df %<>% mutate(barcode = paste0(barcode, '_B'))
  hash_df <- rbind(hashA_df, hashB_df)
  
  
  #plot hash tag distributions per cell line
  #exclude low quality cells
  cell_df %<>% dplyr::mutate(
    barcode_short = str_sub(barcode, 1, 16),
    barcode_short = paste0(barcode_short, '_', channel),
    short_name = str_match(singlet_ID, '^([:alnum:]+)_')[,2],
    hash_tag = factor(hash_tag, levels = c(
      'DMSO_3hr', 'DMSO_6hr', 'DMSO_12hr', 'DMSO_24hr', 'DMSO_48hr',
      'Tram_3hr', 'Tram_6hr', 'Tram_12hr', 'Tram_24hr', 'Tram_48hr',
      'Untreated_48hr', 'multiplet', 'unknown')))
  CL_ord <- cell_df %>% 
    dplyr::filter(cell_quality != 'low_quality') %>% 
    dplyr::group_by(short_name) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::arrange(dplyr::desc(n)) %>% 
    .[['short_name']]
  
  ggplot(cell_df %>% 
           dplyr::filter(cell_quality != 'low_quality') %>% 
           dplyr::mutate(short_name = factor(short_name, levels = CL_ord)), aes(short_name, fill = hash_tag, color = hash_tag == 'unknown')) + 
    geom_histogram(stat = 'count', lwd = 0.25) +
    guides(color = F) + 
    scale_color_manual(values = c(`FALSE` = 'black', `TRUE` = 'green')) +
    ylab('Cell Count') + 
    cdsr::theme_Publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 9)) + 
    guides(fill = guide_legend(nrow = 5)) 
  ggsave(file.path(fig_dir, 'hash_tag_dist_by_CL.png'),
         width = 4.5, height = 4)
  
  
  ## NOW MAKE TSNE PLOT OF HASH TAG COUNTS DATA
  #apply CLR transform to hash tag data
  CLR <- function(x) {
    x <- x + 1
    log(x / (exp(sum(log(x[x > 0]), na.rm=TRUE) / length(x))))
  }
  hash_norm <- apply(hash_df %>% column_to_rownames(var = 'barcode'),
                     1, CLR) %>% t()
  used_bcs <- intersect(cell_df %>% 
                    dplyr::filter(hash_quality != 'unknown', cell_quality != 'low_quality') %>%
                    .[['barcode_short']], 
                  rownames(hash_norm))
  hto <- CreateSeuratObject(t(hash_norm[used_bcs, ]))
  hto.dist.mtx <- as.matrix(dist(hash_norm[used_bcs,]))
  
  aa <- with(cell_df[match(used_bcs, cell_df$barcode_short),], as.character(hash_tag) %>% set_names(barcode_short)) 
  hto %<>% Seurat::AddMetaData(aa, 'hash_tag')
  hto <- RunTSNE(hto, 
                 genes.use = rownames(hto@raw.data), 
                 distance.matrix = hto.dist.mtx,
                 perplexity = globals$tsne_perplexity)
  
  df <- Embeddings(object = hto, reduction = 'tsne') %>% 
    as.data.frame() %>% 
    set_colnames(c('t1', 't2')) %>% 
    rownames_to_column(var = 'barcode_short') %>% 
    left_join(cell_df)
  
  ggplot(df %>% dplyr::mutate(hash_tag = plyr::revalue(hash_tag, c(`Untreated_48hr` = 'Untreated'))), 
         aes(t1, t2)) + 
    geom_point(data = . %>% dplyr::filter(hash_tag != 'multiplet'), aes(fill = hash_tag), size = 0.5, alpha = 0.75, size = 2, pch = 21, stroke = 0.1, color = 'white') + 
    geom_point(data = . %>% dplyr::filter(hash_tag == 'multiplet'), pch = 21, size = 1, fill = 'darkgray', color = 'white', stroke = 0.1) +
    cdsr::theme_Publication() +
    xlim('tSNE 1') + ylab('tSNE 2') +
    guides(fill = guide_legend(override.aes = list(size = 2), nrow = 4)) +
    theme(legend.text = element_text(size = 9),
          legend.title = element_blank())
  
  ggsave(file.path(fig_dir, 'hash_tag_tSNE.png'), width = 4., height = 4.25)
}