
run_predictive_model_compare_pooled <- function() {
  
#PARAMS
min_cells <- 5 #minimum cells per condition to include cell line in analysis
n_top_feat <- 1000
random_seed <- 1
kfold <- 10
max_AUC <- 1.5


doMC::registerDoMC(cores = 3)

used_expts <- c('bortezomib_24hr_expt1',
                'idasanutlin_24hr_expt1',
                'trametinib_24hr_expt3', 
                'dabrafenib_24hr_expt3',
                'BRD3379_24hr_expt3',
                'Afatinib_expt10',
                'navitoclax_24hr_expt3',
                'AZD5591_expt10',
                'Everolimus_expt10',
                'Gemcitabine_expt10',
                'Prexasertib_expt10',
                'Taselisib_expt10',
                'JQ1_expt10')

targ_dsets <- sc_DE_meta[used_expts]

#gather LFC profiles and AUC across target experiments
all_res <- llply(targ_dsets, function(cur_expt) {
  print(sprintf('Processing expt %s', cur_expt$expt_name))
  out_dir <- file.path(results_dir, cur_expt$expt_name)
  
  cell_df <- load_cell_info(cur_expt, QC_filter = TRUE) %>% 
    tidyr::separate(col = 'condition', into = c('condition', 'batch'), sep = '_')
  cell_counts <- cell_df %>% 
    dplyr::group_by(singlet_ID, condition) %>% 
    dplyr::summarise(n = n()) %>% 
    reshape2::dcast(singlet_ID ~ condition, value.var = 'n') %>% 
    dplyr::mutate(tot = control + treat,
                  depMapID = celllinemapr::ccle.to.arxspan(singlet_ID))
  usable_CLs <- cell_counts %>% 
    dplyr::filter(control >= min_cells, treat >= min_cells) %>% 
    .[['singlet_ID']] %>% 
    celllinemapr::ccle.to.arxspan()
  
  CL_df <- all_CL_features[[cur_expt$expt_name]]
  
  dat <- load_compute_CPM_stats(cur_expt, results_dir, prior_counts = globals$pca_prior_cnt, type = 'avg_cpm')
  colnames(dat$LFC_mat) <- celllinemapr::ccle.to.arxspan(colnames(dat$LFC_mat))
  
  y <- CL_df$AUC_avg %>% set_names(CL_df$DEPMAP_ID)
  y <- y[!is.na(y)]
  common_CLs <- intersect(names(y), colnames(dat$LFC_mat))
  res <- list(X = t(dat$LFC_mat[, common_CLs]), 
              y = y[common_CLs], 
              name = cur_expt$expt_name,
              num_cells = cell_counts[match(common_CLs, cell_counts$depMapID), 'tot'])
  
return(res)})

#aggregate X and y data
all_y <- do.call(c, lapply(all_res, function(x) x$y))  
all_y <- pmin(all_y, max_AUC)
comb_X <- do.call(rbind, llply(all_res, function(x) x$X)) %>% 
  set_rownames(names(all_y))

#build and fit global model
all_X <- list(LFC = list(data = comb_X))
input_data <- proc_input_data(all_X, all_y)
fids <- caret::createFolds(input_data$y, k = kfold, list = FALSE) %>% 
  set_names(names(input_data$y))

LFC_results <- get_rf_xv_results(input_data, kfold, fids = fids, n_top_feat = n_top_feat)
tot_corr_R2 <- caret::R2(pred = LFC_results$prediction, obs = LFC_results$obs, formula = 'corr')
tot_trad_R2 <- caret::R2(pred = LFC_results$prediction, obs = LFC_results$obs, formula = 'traditional')

LFC_results %<>% separate(DEPMAP_ID, into = c('treat', 'cl'), sep = '\\.')



#fit single model to all data to assess feature importance
used_feat <- which(apply(input_data$X, 2, sd) > 0) #only use genes where there is non-zero variance
rf_df <- cbind(input_data$X[, used_feat], data_frame(response = input_data$y))
mod_fit <- ranger::ranger(data = rf_df, importance = 'impurity', dependent.variable.name = 'response')

cc <- cor(input_data$X[, used_feat], as.matrix(input_data$y, ncol = 1), use = 'pairwise.complete.obs')

gene_stats <- tibble(imp = mod_fit$variable.importance, 
                    Gene = names(mod_fit$variable.importance)) %>% 
                      dplyr::mutate(cor = -cc) %>% #note minus to account for sensitivity = 1-AUC 
 dplyr::mutate(Gene = str_replace(Gene, '_LFC', ''))

imp_threshold <- 0.05
ggplot(gene_stats, aes(cor, imp)) + 
  geom_point(fill = 'black', pch = 21, size = 2, alpha = 0.8, color = 'white', stroke = 0.1) +
  geom_label_repel(data = gene_stats %>% arrange(dplyr::desc(imp)) %>% head(15), aes(label = Gene), size = 2.5, label.padding = 0.1) +
  geom_hline(yintercept = imp_threshold, linetype = 'dashed') +
  labs(x = 'LFC-sensitivity\ncorrelation', y = 'Feature importance') +
  theme_Publication()
ggsave(file.path(fig_dir, 'global_model_feature_import.png'),
       width = 4, height = 3.5)

up_genes <- gene_stats %>% 
 filter(cor > 0, imp > imp_threshold) %>% 
 arrange(dplyr::desc(imp)) %>% 
 .[['Gene']]
down_genes <- gene_stats %>% 
 filter(cor < 0, imp > imp_threshold) %>% 
 arrange(dplyr::desc(imp)) %>% 
 .[['Gene']]
all_genes <- str_replace_all(colnames(input_data$X), '_LFC', '')
res_down <- run_GSAhyper(down_genes, all_genes, gsc_data$combined)
res_up <- run_GSAhyper(up_genes, all_genes, gsc_data$combined)

make_stem_plot_precom(res_up, res_down, n_lab_per = 5)
ggsave(file.path(fig_dir, 'global_feat_imp_gsea.png'),
       width = 5.5, height = 2.5)


## COMPARE TO L1k analysis
# andy_res <- read_csv('~/Downloads/global_viability_sig_l1k.csv')
# gene.info <- load.from.taiga(data.name='l1000-metadata-8489', data.version=1, data.file='gene_info')
# andy_res %<>% left_join(gene.info %>% dplyr::select(Gene = pr_gene_symbol, pr_is_lm), by = 'Gene')
# 
# comb <- inner_join(gene_stats, andy_res, by = 'Gene', suffix = c('_MIXseq', 'L1k'))
# ggplot(comb, aes(cor, logFC)) + 
#   geom_point() + 
#   geom_abline()
# res <- lin_associations(input_data$X[, used_feat], input_data$y, W = NULL, robust.se = F, scale.A = F, scale.y = F, is.dependent = F, var.th = 0, coef.var = F)
# res %<>% mutate(Gene = str_replace(feature, '_LFC', ''),
#                 betahat = -betahat)
# comb %<>% inner_join(res, by = 'Gene')
# df <- comb %>% filter(pr_is_lm == 1)
# top_genes <- 10
# lab_genes <- c(
#   df %>%
#     dplyr::arrange(desc(abs(logFC))) %>%
#     head(top_genes) %>%
#     .[['Gene']],
#   df %>%
#     dplyr::arrange(desc(abs(betahat))) %>%
#     head(top_genes) %>%
#     .[['Gene']]
# )
# ggplot(comb, aes(logFC, betahat)) + 
#   geom_point(pch = 21, fill = 'lightgray', color = 'white', stroke = 0.1, size = 1.5, alpha = 0.8) +
#   geom_point(data = df %>% filter(pr_is_lm == 1), pch = 21, fill = 'black', color = 'white', stroke = 0.1, size = 2) +
#   # geom_abline() +
#   ggplot2::geom_vline(xintercept = 0, linetype = 'dashed') +
#   ggplot2::geom_hline(yintercept = 0, linetype = 'dashed') +
#   theme_Publication() +
#   xlab('L1k') + ylab('MIX-Seq') +
#   geom_label_repel(data = df %>% filter(Gene %in% lab_genes), aes(label = Gene), size = 2)
# ggsave(file.path(fig_dir, 'L1k_vs_MIXseq_viability_sig.png'),
#        width = 4, height = 3.5)
}
