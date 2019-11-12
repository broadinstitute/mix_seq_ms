
run_cell_downsampling_analysis <- function(cur_expt) {

n_reps <- 5
n_cells <- c(1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

sprintf('Processing expt %s', cur_expt$expt_name)
out_dir <- file.path(results_dir, cur_expt$expt_name)
res_file <- file.path(results_dir, sprintf('cell_downsampling_results_%s.rds', cur_expt$expt_name))

limma_res <- read_rds(file.path(out_dir, 'limma_res.rds'))

CL_df <- all_CL_features[[cur_expt$expt_name]]


dat <- load_collapsed_profiles(cur_expt, results_dir, type = 'sum')
used_genes <- which(rowSums(dat$profile_mat > min_counts_per_gene) > ncol(dat$profile_mat)*min_det_samples)

seuObj <- load_sc_data(cur_expt, sc_expts = sc_expts)

seuObj <- Seurat::NormalizeData(object = seuObj, 
                                normalization.method = "LogNormalize", 
                                scale.factor = 1e6) #CPM normalization

cell_df <- seuObj@meta.data %>% 
  rownames_to_column(var = 'barcode') %>% 
  dplyr::mutate(condition = ifelse(grepl('control', condition), 'control', 'treat')) 

#print avg cells per condition
tot_cells_df <- cell_df %>% 
  dplyr::group_by(singlet_ID, condition) %>% 
  dplyr::summarise(n = n()) %>% 
  spread(key = 'condition', value = 'n') %>% 
  dplyr::mutate(min_n = pmin(control, treat)) %>% 
  dplyr::mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(singlet_ID)) %>% 
  left_join(CL_df)

avg_cells <- mean(c(mean(tot_cells_df$control), mean(tot_cells_df$treat)))
print(avg_cells)

#calculate avg LFC
calc_LFC <- function(cur_Seu, ugenes) {
  cond_avgs <- Seurat::AverageExpression(cur_Seu, features = ugenes, 
                                         use.scale = FALSE, return.seurat = FALSE, add.ident = 'condition', verbose = FALSE)$RNA %>% t() 
  cond_lcpm <- log2(globals$prior_cnt + cond_avgs)
  avg_control <- colMeans(cond_lcpm[grepl('control', rownames(cond_lcpm)), ,drop=F], na.rm=T)
  avg_treat <- colMeans(cond_lcpm[!grepl('control', rownames(cond_lcpm)), ,drop=F], na.rm=T)
  avg_LFC <- avg_treat - avg_control
  return(avg_LFC)
}

targ_CLs <- tot_cells_df %>% 
  dplyr::filter(min_n >= 100) %>% 
  dplyr::pull(singlet_ID)

all_avg_cors <- ldply(targ_CLs, function(cur_CL) {
  all_cells = cell_df %>% 
    dplyr::filter(singlet_ID == cur_CL) %>% 
    group_by(condition) %>% 
    sample_n(size = 100) %>% 
    .[['barcode']]
  seuSub <- subset(seuObj, cells = all_cells)
  seuSub <- FindVariableFeatures(object = seuSub,
                                 nfeatures = globals$n_highvar_genes,
                                 do.plot = FALSE,
                                 selection.method = 'vst')
  ugenes <- VariableFeatures(seuSub)
  avg_LFC <- calc_LFC(seuSub, ugenes = ugenes)
  
  calc_sub_cor <- function(n_per_group, n_reps) {
    map_dbl(seq(n_reps), function(ii) {
      set.seed(ii)
      cur_cells <- cell_df %>% 
        dplyr::filter(barcode %in% all_cells) %>% 
        dplyr::group_by(condition) %>% 
        dplyr::sample_n(size = n_per_group) %>% 
        .[['barcode']]
      seuSubs <- seuSub[, cur_cells]
      sub_avg_LFC <- calc_LFC(seuSubs, ugenes = ugenes)
      return(cor(sub_avg_LFC, avg_LFC, use = 'pairwise.complete.obs'))
    }) %>% mean()
  }
  avg_sub_cor <- map_dbl(n_cells, calc_sub_cor, n_reps = n_reps)
  return(tibble(cor = avg_sub_cor, n_cells = n_cells, CCLE_ID = cur_CL))
})
all_avg_cors %<>% left_join(CL_df, by = "CCLE_ID")

write_rds(all_avg_cors, res_file)

ggplot(all_avg_cors, aes(n_cells, cor, color = sens, group = CCLE_ID)) + 
  geom_point() + 
  geom_line() + 
  guides(color = guide_colorbar(title = 'sensitivity', barwidth = 8)) +
  ylim(0, 1) +
  labs(x = 'N. cells per group', y = 'Correlation with\nall-cell profile') +
  cdsr::theme_Publication()
ggsave(file.path(fig_dir, 'trametinib_cell_downsampling.png'), width = 3.5, height = 3.5)

}