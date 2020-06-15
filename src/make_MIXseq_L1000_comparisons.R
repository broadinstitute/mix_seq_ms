

# read in L1000 data, just for overlapping drugs using the cmapR package
get_L1000_subset <- function(sig_file, LINCS_level5_file, hgnc_data,
                             drugs = c('bortezomib', 'navitoclax', 'gemcitabine',
                                        'trametinib', 'everolimus', 'JQ1')) {
  library(cmapR)
  library(dplyr)
  sig_info <- data.table::fread(sig_file) %>% as.data.frame()
  sample_names <- c()
  for(drug in drugs) {
    sample_names <- c(sample_names, sig_info$sig_id[grep(drug, sig_info$pert_iname, ignore.case = T)])
  }
  drug_mat_sample_names <- cmapR::read_gctx_ids(LINCS_level5_file, dim='column')
  drug_mat_column_inds <- which(drug_mat_sample_names %in% sample_names)
  
  L1000_drug_mat <- cmapR::parse_gctx(LINCS_level5_file, cid = drug_mat_column_inds)
  common_genes <- which(rownames(L1000_drug_mat@mat) %in% hgnc_data$entrez_id)
  L1000_drug_mat <- L1000_drug_mat@mat[common_genes,]
  hgnc_data <- filter(hgnc_data, entrez_id %in% rownames(L1000_drug_mat))
  rownames(hgnc_data) <- hgnc_data$entrez_id
  rownames(L1000_drug_mat) <- hgnc_data[rownames(L1000_drug_mat),'symbol']
  return(L1000_drug_mat)
  
}
# make plots comparing MIXseq and L1000 LINCS Phase 1 data downloaded from clue.io
# input: MIXseq avg logFC matrices per experiment, L1000 signatures downloaded from clue.io, L1000 gene info. 
# marking genes as inferred or landmark, current MIXseq experiment, current L1000 drug, plot title and parameters.
L1000_Phase1_avg_response_comparison <- function(MIXseq_avg_response, l1000_drug_mat, l1000_sig_file, l1000_gene_info, cur_experiment,
                                          cur_drug, title="", text_xlim=c(NA,NA), text_ylim = c(NA,NA)) {
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  
  l1000_sig_info <- data.table::fread(l1000_sig_file) %>% as.data.frame()
  # get MIXseq avg logFC values for the current experiment
  MIXseq_exp_avg_response <- MIXseq_avg_response[[cur_experiment]]$logFC %>% set_names(MIXseq_avg_response[[cur_experiment]]$Gene)
  
  cur_l1000_samples <- filter(l1000_sig_info, pert_iname == cur_drug)$sig_id
  l1000_drug_avg_response <- rowMeans(L1000_drug_mat[,cur_l1000_samples]) %>% set_names(rownames(L1000_drug_mat))

  common_genes <- intersect(names(l1000_drug_avg_response), names(MIXseq_exp_avg_response))
  
  cur_avg_response <- data.frame(Gene = common_genes, MIXseq = MIXseq_exp_avg_response[common_genes], L1000 = l1000_drug_avg_response[common_genes],
                        gene = ifelse(common_genes %in% filter(l1000_gene_info, pr_is_lm == 1)$pr_gene_symbol, 'landmark', 'inferred'))
  
  cur_avg_response$max.logFC <- apply(cur_avg_response, 1, function(x) max(abs(as.numeric(x[2])), abs(as.numeric(x[3]))))
  cur_avg_response$gene <- factor(cur_avg_response$gene, levels = c('inferred', 'landmark'))
  cur_avg_response$gene_label <- ""
  cur_avg_response <- cur_avg_response %>% dplyr::arrange(dplyr::desc(max.logFC))
  cur_avg_response[1:20,'gene_label'] <- cur_avg_response[1:20,'Gene'] %>% as.character()
  
  LFC_comparison_plot <- ggplot(cur_avg_response,  aes(MIXseq, L1000, color=gene)) + 
    geom_point(alpha=0.6, size=.75) +
    geom_point(data=filter(cur_avg_response, gene=='landmark'),  aes(MIXseq, L1000, color=gene), size=.75, alpha=0.5) + 
    ggpubr::stat_cor(label.y.npc = 'top', size=2, show.legend = F) +
    geom_smooth(method = 'lm') +
    geom_vline(xintercept=0, linetype="dashed", color = "gray80") +
    geom_hline(yintercept=0, linetype="dashed", color = "gray80") +
    xlab("MIX-Seq avg LFC") +
    ylab("L1000 avg mod-z") +
    #ggrepel::geom_text_repel(data=cur_avg_response,
    #                         aes(label=gene_label), size=2, show.legend = FALSE, xlim=text_xlim, ylim=text_ylim) +
    ggtitle(title) +
    publication_theme() + 
    theme(text=element_text(size=6), 
          axis.text = element_text(size=6), 
          axis.title = element_text(size=7))
  
  return(LFC_comparison_plot)
  
}

# make plots comparing MIXseq and L1000 LINCS Phase 2 data using data downloaded from from http://amp.pharm.mssm.edu/Slicr
# input: MIXseq avg logFC matrices per experiment, L1000 signatures downloaded from LINCS, L1000 meta data matrix downloaded from LINCS,
# L1000 gene info. marking genes as inferred or landmark,current MIXseq experiment, current L1000 drug, plot title and parameters.
L1000_Phase2_avg_response_comparison <- function(MIXseq_avg_response, l1000_LINCS_drug_mat, l1000_LINCS_meta, l1000_gene_info, cur_experiment,
                                          cur_drug, title="", text_xlim=c(NA,NA), text_ylim = c(NA,NA)) {
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(ggpubr)
  
  # get MIXseq avg logFC values for the current experiment
  MIXseq_exp_avg_response <- MIXseq_avg_response[[cur_experiment]]$logFC %>% set_names(MIXseq_avg_response[[cur_experiment]]$Gene)
  
  # parse L1000 matrices from LINCS
  l1000_LINCS_drug_mat <- l1000_LINCS_drug_mat %>% dplyr::select(-`signatures in which the gene is significant`)
  l1000_LINCS_drug_mat <- filter(l1000_LINCS_drug_mat, `gene symbol` != "-666") # remove genes with label -666 (NAs in the data)
  gene_list <- l1000_LINCS_drug_mat$`gene symbol`
  l1000_LINCS_drug_mat <- l1000_LINCS_drug_mat %>%
    dplyr::select(-`gene symbol`)
  l1000_LINCS_drug_mat <- aggregate(l1000_LINCS_drug_mat, by=list(gene_list), FUN=mean) # average values across duplicated genes
  l1000_LINCS_drug_mat <- l1000_LINCS_drug_mat %>% column_to_rownames("Group.1")
  
  # get L1000 average values for the current drug
  l1000_samples <- l1000_LINCS_meta$sig_id[grep(cur_drug, l1000_LINCS_meta$pert_desc)]
  l1000_drug_response <- l1000_LINCS_drug_mat[,l1000_samples] %>% as.matrix()
  l1000_drug_avg_response <- rowMeans(l1000_drug_response) %>% set_names(rownames(l1000_drug_response))
  
  common_genes <- intersect(names(l1000_drug_avg_response), names(MIXseq_exp_avg_response))
  
  cur_avg_response <- data.frame(Gene = common_genes, MIXseq = MIXseq_exp_avg_response[common_genes], L1000 = l1000_drug_avg_response[common_genes],
                                 gene = ifelse(common_genes %in% filter(l1000_gene_info, pr_is_lm == 1)$pr_gene_symbol, 'landmark', 'inferred'))
  
  cur_avg_response$max.logFC <- apply(cur_avg_response, 1, function(x) max(abs(as.numeric(x[2])), abs(as.numeric(x[3]))))
  cur_avg_response$gene <- factor(cur_avg_response$gene, levels = c('inferred', 'landmark'))
  cur_avg_response$gene_label <- ""
  cur_avg_response <- cur_avg_response %>% arrange(desc(max.logFC))
  cur_avg_response[1:20,'gene_label'] <- cur_avg_response[1:20,'Gene'] %>% as.character()
  
  LFC_comparison_plot <- ggplot(cur_avg_response, aes(MIXseq, L1000, color=gene)) + 
    geom_point(alpha=0.6, size=.75) +
    geom_point(data=filter(cur_avg_response, gene=='landmark'), aes(MIXseq, L1000, color=gene), size=.75, alpha=0.5) + 
    ggpubr::stat_cor(label.y.npc = 'top', size=2, show.legend = F) +
    ggplot2::geom_smooth(method = 'lm') +
    geom_vline(xintercept=0, linetype="dashed", color = "gray50") +
    geom_hline(yintercept=0, linetype="dashed", color = "gray50") +
    xlab("MIX-Seq avg LFC") +
    ylab("L1000 avg mod-z") +
    #ggrepel::geom_text_repel(data=cur_avg_response,
    #                         aes(label=gene_label), size=2, show.legend = FALSE, xlim=text_xlim, ylim=text_ylim) +
    ggtitle(title) +
    publication_theme() + 
    theme(text=element_text(size=6), 
                                      axis.text = element_text(size=6), 
                                      axis.title = element_text(size=7))
  
  return(LFC_comparison_plot)
  
}

publication_theme <- function (base_size = 12, base_family = "Helvetica") 
{
  library(grid)
  library(ggthemes)
  library(ggplot2)
  (theme_foundation(base_size = base_size, base_family = base_family) + 
      theme(plot.title = element_text(face = "bold", size = rel(1.2), 
                                      hjust = 0.5), text = element_text(), panel.background = element_rect(colour = NA), 
            plot.background = element_rect(colour = NA), panel.border = element_rect(colour = NA), 
            axis.title = element_text(face = "bold", size = rel(1)), 
            axis.title.y = element_text(angle = 90, vjust = 2), 
            axis.title.x = element_text(vjust = -0.2), axis.text = element_text(), 
            axis.line = element_line(colour = "black"), axis.ticks = element_line(), 
            panel.grid.major = element_line(colour = "#f0f0f0"), 
            panel.grid.minor = element_blank(), legend.key = element_rect(colour = NA), 
            legend.position = "bottom", legend.text = element_text(size = rel(1.2)), 
            legend.direction = "horizontal", legend.key.size = unit(0.3, 
                                                                    "cm"), legend.margin = unit(0, "cm"), legend.title = element_text(face = "italic"), 
            plot.margin = unit(c(5, 5, 5, 5), "mm"), strip.background = element_rect(colour = "#f0f0f0", 
                                                                                     fill = "#f0f0f0"), strip.text = element_text(face = "bold")))
}

# input files required:
# MIXSEQ_all_avg_responses.rds - calculation of average, per experiment log-fold change values for the MIXseq data
# hgnc.complete.set.txt - gene symbol matrix, downloaded from ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt 9/18/2019
# GSE92742_Broad_LINCS_sig_info.txt - matrix of L1000 phase 1 sample ids and drug info, downloaded from clue.io
# GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx - matrix of L1000 level 5 LINCS Phase 1 data, downloaded from clue.io
# GSE92742_Broad_LINCS_gene_info.txt - matrix of gene info, downloaded from clue.io
# level5_CD_matrix.csv - matrix of L1000 level 5 LINCS Phase 2, downloaded from http://amp.pharm.mssm.edu/Slicr
# level5_CD_meta.csv - matrix of L1000 level 5 LINCS Phase 2 metadata, downloaded from http://amp.pharm.mssm.edu/Slicr
# inputs:
# data_folder: folder to read required data files from and write plot to
create_L1000_MIXseq_comparison_fig <- function(data_folder) {
  
  MIXseq_avg_response <- readRDS(file.path(data_folder, "MIXSEQ_all_avg_responses.rds"))
  hgnc_data <- data.table::fread(file.path(data_folder, 'hgnc.complete.set.txt')) %>% as.data.frame()
  
  l1000_drug_mat <- get_L1000_subset(file.path(data_folder, "GSE92742_Broad_LINCS_sig_info.txt"), 
                                     file.path(data_folder, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"),
                                     hgnc_data)
  
  l1000_gene_info <- data.table::fread(file.path(data_folder, "GSE92742_Broad_LINCS_gene_info.txt")) %>% as.data.frame()
  
  sig_file <- file.path(data_folder, "GSE92742_Broad_LINCS_sig_info.txt")
  
  l1000_LINCS_drug_mat <- read_csv(file.path(data_folder, "level5_CD_matrix.csv"))
  l1000_LINCS_meta <- read_csv(file.path(data_folder, "level5_CD_meta.csv"))
  
  # Bortezomib plot using LINCS Phase 1 data
  g1 <- L1000_Phase1_avg_response_comparison(MIXseq_avg_response, l1000_drug_mat,
                                      sig_file, l1000_gene_info,
                                      'bortezomib_24hr_expt1', 'bortezomib', 'Bortezomib',
                                      text_xlim=c(2,11))
  # Navitoclax plot using LINCS Phase 1 data
  g2 <- L1000_Phase1_avg_response_comparison(MIXseq_avg_response, l1000_drug_mat, 
                                      sig_file, l1000_gene_info, 
                                      'navitoclax_24hr_expt3', 'navitoclax', 'Navitoclax',
                                      text_ylim = c(-1.25,1.1))
  # Gemcitabine plot using LINCS Phase 1 data
  g3 <- L1000_Phase1_avg_response_comparison(MIXseq_avg_response, l1000_drug_mat, 
                                      sig_file, l1000_gene_info, 
                                      'Gemcitabine_expt10', 'gemcitabine', 'Gemcitabine')
  # Trametinib plot using LINCS Phase 2 data
  g4 <- L1000_Phase2_avg_response_comparison(MIXseq_avg_response, l1000_LINCS_drug_mat, 
                                            l1000_LINCS_meta, l1000_gene_info, 
                                            'trametinib_24hr_expt3', 'trametinib', 'Trametinib')
  # Everolimus plot using LINCS Phase 2 data
  g5 <- L1000_Phase2_avg_response_comparison(MIXseq_avg_response, l1000_LINCS_drug_mat, 
                                            l1000_LINCS_meta, l1000_gene_info, 
                                            'Everolimus_expt10', 'everolimus', 'Everolimus',
                                            text_xlim = c(0,1.5))
  # JQ1 plot using LINCS Phase 2 data
  g6 <- L1000_Phase2_avg_response_comparison(MIXseq_avg_response, l1000_LINCS_drug_mat, l1000_LINCS_meta, l1000_gene_info, 
                                            'JQ1_expt10', '(+)JQ1', 'JQ1',
                                            text_xlim = c(0.1,4))
  
  l1000_comparison_fig <- grid.arrange(g1, g2, g3, g4, g5, g6)
  
  ggsave(file.path(data_folder, "L1000_MIXseq_comparison.png"), l1000_comparison_fig, width=5.6, height=7.7)
  
}
