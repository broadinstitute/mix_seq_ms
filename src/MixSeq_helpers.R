library(magrittr)
library(plyr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(taigr)


###------------------------- DATA IO -----------------------------------

#' load specified single-cell dataset into seurat object
#'
#' @param expt_params: list of parameters defining the desired dataset  
#' @param min_cells_per_gene: optional minimum number of detected genes per cell 
#' @param min_genes_per_cell: optional minimum number of cells where a gene is detected 
#' @param QC_filter: use only cells called 'normal' (not doublet or low quality) 
#' @param local: force loading data from local data directory (default FALSE)
#'
#' @return: Seurat object with dataset
#' @export
#'
#' @examples
load_sc_data <- function(expt_params, sc_expts, min_cells_per_gene = 0, min_genes_per_cell = 0, QC_filter = TRUE, local = globals$local_data) {
  if ('data_sets' %in% names(expt_params)) {
    #if loading a list of datasets
    if (expt_params$load_from_taiga & !local) {
      data_sets <- lapply(expt_params$data_sets, function(dset) {
        cdsc::load_sc_data_with_expt_string(sc_expts[[dset]]$taiga_name, data.version = sc_expts[[dset]]$taiga_version)
      })
    } else {
      data_dir <- here::here('data')
      data_sets <- lapply(expt_params$data_sets, function(dset) {
        load_expt_data(file.path(data_dir, dset))
      })
    }
    dat <- merge_sc_data(data_sets)
  } else {
    if (!is.null(expt_params$taiga_name) & !local) {
      dat <- cdsc::load_sc_data_with_expt_string(expt_params$taiga_name, data.version = expt_params$taiga_version)
    } else {
      data_dir <- here::here('data')
      dat <- load_expt_data(file.path(data_dir, expt_params$expt_name))
    }
    dat$cell_info %<>% tibble::column_to_rownames(var = 'barcode')
  }
  dat %<>% convert_to_hugo_symbols() #convert genes to hugo symbols (and only keep unique hugo symbol genes)
  
  #make into Seurat object
  seuObj <- CreateSeuratObject(dat$counts, 
                               min.cells = min_cells_per_gene,
                               min.features = min_genes_per_cell,
                               meta.data = dat$cell_info)
  #make singlet classifications the cell identifiers
  seuObj <- SetIdent(seuObj, value = 'singlet_ID')
  
  if (QC_filter) {
    #get rid of doublets and low-quality cells
    cq <- Seurat::FetchData(seuObj, vars = 'cell_quality')
    seuObj <- seuObj[, which(cq$cell_quality == 'normal')]
  }
  
  #store gene info here
  seuObj@misc <- dat$gene_info
  
  #add mitochondrial gene fraction
  seuObj[["percent.mito"]] <- PercentageFeatureSet(seuObj, pattern = "^MT-")
  
  #add cellular detection rate (fraction of genes detected in a given cell)
  seuObj[['cell_det_rate']] <- seuObj$nFeature_RNA/nrow(GetAssayData(seuObj))
  return(seuObj)
}


#' Load counts data and associated cell stats for an experiment
#'
#' @param expt_dir: Full path to the experiment directory 
#'
#' @return: List containing counts, cell_info (with SNP classifications) and gene info
#' @export
#'
#' @examples
load_expt_data <- function(expt_dir) {
  library(readr)
  #Load single-cell counts matrix
  rMat <- Matrix::readMM(file.path(expt_dir, 'matrix.mtx'))
  genes <- read_tsv(file.path(expt_dir, 'genes.tsv'), col_names = F, col_types = cols()) %>% 
    set_colnames(c('Ensembl_ID', 'Gene_Symbol'))
  barcodes <- read_tsv(file.path(expt_dir, 'barcodes.tsv'), col_names = F, col_types = cols()) 
  classifications <- read_csv(file.path(expt_dir, 'classifications.csv'), col_types = cols()) %>% 
    as.data.frame() 
  stopifnot('barcode' %in% colnames(classifications))
  stopifnot(nrow(classifications) == nrow(barcodes))
  colnames(rMat) <- barcodes$X1
  rownames(rMat) <- genes$Ensembl_ID
  num_cells <- ncol(rMat)
  num_genes <- nrow(rMat)  
  cat(sprintf('Loaded matrix with %d cells and %d genes\n', num_cells, num_genes))
  return(list(counts = rMat, gene_info = genes, cell_info = classifications))
}


#' Load table of cell info from a specified set of experiments
#'
#' @param expt_params: parameters defining desired experiment 
#' @param QC_filter: whether or not to use only cells called 'normal' 
#' @param local: force loading data from local data directory (default = FALSE)
#'
#' @return
#' @export
#'
#' @examples
load_cell_info <- function(expt_params, QC_filter = TRUE, local = globals$local_data) {
  if ('data_sets' %in% names(expt_params)) {
    if (expt_params$load_from_taiga & !local) {
      df <- ldply(expt_params$data_sets, function(dset) {
        taiga_name <- cdsc:::get_permaname_from_expt_string(sc_expts[[dset]]$taiga_name) 
        taigr::load.from.taiga(data.name=taiga_name, data.file = 'classifications', data.version = sc_expts[[dset]]$taiga_version)
      }, .id = 'condition')
    } else {
      data_dir <- here::here('data')
      df <- ldply(expt_params$data_sets, function(dset) {
        read_csv(file.path(data_dir, dset, 'classifications.csv'), col_types = cols()) %>% 
          as.data.frame() 
      }, .id = 'condition')
    }
  } else {
    if (!is.null(expt_params$taiga_name) & !local) {
      taiga_name <- cdsc:::get_permaname_from_expt_string(expt_params$taiga_name) 
    df <- taigr::load.from.taiga(data.name=taiga_name, data.file = 'classifications', data.version = expt_params$taiga_version)
    } else {
      data_dir <- here::here('data')
      df <- read_csv(file.path(data_dir, expt_params$expt_name, 'classifications.csv'), col_types = cols()) %>% 
        as.data.frame()  
    }
  }
  
  stopifnot('barcode' %in% colnames(df))
  if (QC_filter) {
    df %<>% dplyr::filter(cell_quality == 'normal')
  }
  return(df)
}


# helper function that gets appropriate path for data when stored on Google Drive
get_Gdrive_path <- function(expt_params) {
  if (dir.exists("~/google_drive/Project Apollo/")) {
    gdrive_path <- "~/google_drive/Project Apollo/"
  } else if (dir.exists('/Volumes/GoogleDrive/My Drive/Project Apollo')) {
    gdrive_path <- '/Volumes/GoogleDrive/My Drive/Project Apollo'
  } else {
    stop('Need to add Gdrive path')
  }
  if (expt_params$expt_batch == 'expt3') {
    data_dir <- paste0(gdrive_path, '/Apollo drugs /scRNAseq_003_004')
  } else if (expt_params$expt_batch == 'expt2') {
    data_dir <- paste0(gdrive_path, '/Apollo genetics/Apollo_scRNAseq_002')
  } else if (expt_params$expt_batch == 'expt1') {
    data_dir <- paste0(gdrive_path, '/Apollo drugs /scRNAseq Pilot 001')
  } else if (expt_params$expt_batch == 'expt5') {
    data_dir <-paste0(gdrive_path, '/Apollo drugs /scRNAseq_005')
  } else if(expt_params$expt_batch == 'expt6') {
    data_dir <- paste0(gdrive_path, '/Apollo genetics/Apollo_scRNAseq_006')
  } else if(expt_params$expt_batch == 'expt9') {
    data_dir <- paste0(gdrive_path, '/Prensner')
  } else if (expt_params$expt_batch == 'expt10') {
    data_dir <- paste0(gdrive_path,  '/Apollo drugs /scRNAseq_010')
  } else if(expt_params$expt_batch == 'expt11') {
    data_dir <- paste0(gdrive_path,  '/Apollo genetics/Apollo_scRNAseq_011')
  }
  return(data_dir)
}


#' Load out-of-pool classifications file
#'
#' @param expt_params 
#'
#' @return
#' @export
#'
#' @examples
load_out_of_pool_class <- function(expt_params) {
  if (!is.null(expt_params$taiga_name)) {
    warning('out-of-pool classifications not yet on Taiga, pulling from Google Drive')
  }
  data_dir <- get_Gdrive_path(expt_params)
  if ('data_sets' %in% names(expt_params)) {
      df <- ldply(expt_params$data_sets, function(dset) {
      read_csv(file.path(data_dir, dset, 'out_of_pool_classifications.csv')) %>% 
        as.data.frame() 
    }, .id = 'condition')
  } else{
    short_dname <- str_replace(expt_params$expt_name, '_expt[0-9]+$', '')
    df <- read_csv(file.path(data_dir, short_dname, 'out_of_pool_classifications.csv')) %>% 
      as.data.frame() 
  }
  return(df)
}


#' Load sum-collapsed counts data from a specified expt
#'
#' @param expt_data 
#' @param results_dir 
#' @param type either 'sum' or 'avg' to specify which collapsed profiles
#'
#' @return list with combined counts_mat and sample info dataframe
#' @export
#'
#' @examples
load_collapsed_profiles <- function(expt_data, results_dir, type) {
  stopifnot(type %in% c('sum', 'avg_cpm'))  
  if (type %in% c('sum')) {
    dfile <- 'summed_counts.rds'
   } else if (type == 'avg_cpm') {
    dfile <- 'avg_profiles_cpm.rds'
  }
  all_mats <- llply(names(expt_data$data_sets) %>% set_names(names(expt_data$data_sets)), function(cur_dset) {
    cur_file <- file.path(results_dir,
                          expt_data$data_sets[[cur_dset]],
                          dfile) %>% stringr::str_replace_all('_combined', '') #ignore "combined" suffix
    cur_mat <- read_rds(cur_file)
  })

  common_genes <- Reduce(intersect, lapply(all_mats, colnames))
  all_profiles <- do.call(rbind, lapply(all_mats, function(mat) {mat[, common_genes]})) %>% t()
  sample_info <- ldply(names(all_mats), function(cur_name) {
    data.frame(CCLE_ID = rownames(all_mats[[cur_name]]), stringsAsFactors = F, check.names = F) %>%
      mutate(expt = cur_name)
  }) %>% mutate(sample_name = paste0(CCLE_ID, '.', expt)) %>%
    tidyr::separate(col = 'expt', into = c('condition', 'batch'), sep = '_', remove = FALSE)
  colnames(all_profiles) <- sample_info$sample_name
  
  return(list(profile_mat = all_profiles, sample_info = sample_info))
}


#' Load a set of collapsed expression profiles across a set of treat and control experiments. Either sum- or avg-collapsed across single cells
#'
#' @param expt_info list specifying expt details
#' @param results_dir directory with saved summed_counts.rds results
#' @param prior_counts prior_counts (default 1)
#' @param use_CLs subset to list of used cell lines (optional)
#' @param type Whether to use summed or averaged data per cell line (must be either 'sum' or 'avg')
#'
#' @return
#' @export
#'
#' @examples
load_compute_CPM_stats <- function(expt_info, results_dir, type = 'sum', prior_counts = 1, use_CLs = NULL) {
  
  dat <- load_collapsed_profiles(expt_info, results_dir, type)
  
  if (!is.null(use_CLs)) {
    #restrict to usable cell lines
    usamps <- dat$sample_info %>% dplyr::filter(CCLE_ID %in% use_CLs) %>% .[['sample_name']]
    dat$profile_mat <- dat$profile_mat[, usamps]  
    dat$sample_info %<>% dplyr::filter(sample_name %in% usamps)  
  }
  
  all_CLs <- unique(dat$sample_info$CCLE_ID)
  if (type == 'sum') {
    dat$profile_mat <- edgeR::DGEList(dat$profile_mat)
    dat$profile_mat <- edgeR::calcNormFactors(dat$profile_mat, method = 'TMMwzp')
    dat$profile_mat %<>% edgeR::cpm(prior.count = prior_counts, log = TRUE)
  }
  if (type == 'avg_cpm') {
    dat$profile_mat <- log2(prior_counts + dat$profile_mat)
  }
  
  #average CPM per sample across set of control and treatment conditions
  control_CPM <- ldply(all_CLs, function(cur_CL) {
    rowMeans(dat$profile_mat[, dat$sample_info %>% dplyr::filter(condition == 'control', CCLE_ID == cur_CL) %>% .[['sample_name']],drop=F])
  }) %>% t() %>% set_colnames(all_CLs)
  treat_CPM <- ldply(all_CLs, function(cur_CL) {
    rowMeans(dat$profile_mat[, dat$sample_info %>% dplyr::filter(condition == 'treat', CCLE_ID == cur_CL) %>% .[['sample_name']],drop=F])
  }) %>% t() %>% set_colnames(all_CLs)
  
  #compute LFC per sample
  stopifnot(all.equal(rownames(treat_CPM), rownames(control_CPM)))
  stopifnot(all.equal(colnames(treat_CPM), colnames(control_CPM)))
  LFC_mat <- treat_CPM - control_CPM
  return(list(LFC_mat = LFC_mat, control_CPM = control_CPM, treat_CPM = treat_CPM))
}


###------------------------- DATA MANIPULATION -----------------------------------

#' Take a list of data objects and merge them into a single object
#'
#' @param data_list: Named list of sc data sets (either preloaded, or list of data paths)
#' @param target_CLs: List of CLs to use (based on SNP singlet classifications)
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#'
#' @return: Data list with counts, gene_info and cell_info
#' @export
merge_sc_data <- function(data_list, target_CLs = NULL) {
  
  tot_cells <- sum(sapply(data_list, function(dset) {ncol(dset$counts)}))
  
  all_sub_data <- plyr::llply(seq(length(data_list)), function(ii) {
    cur_data <- data_list[[ii]]
    if (is.null(target_CLs)) {
      target_CLs <- unique(cur_data$cell_info$singlet_ID)
      keep_cells <- with(cur_data$cell_info, which(singlet_ID %in% target_CLs))
      cur_data$counts <- cur_data$counts[, keep_cells]
      cur_data$cell_info <- cur_data$cell_info[keep_cells,]
    }
    #append batch number to cell barcode to ensure there are no conflicts after merge
    colnames(cur_data$counts) %<>% paste0(., '_', ii)
    cur_data$cell_info %<>% dplyr::mutate(barcode = paste0(barcode, '_', ii))
    return(cur_data)
  }) %>% set_names(names(data_list))
  
  #make sure gene info is aligned across all matrices being merged
  if (length(data_list) > 1) {
    for (ii in seq(2, length(data_list)))
    stopifnot(all.equal(data_list[[1]]$gene_info, data_list[[ii]]$gene_info))
  }
  all_gene_info <- data_list[[1]]$gene_info
  
  #merge counts
  all_counts <- Reduce(cbind, plyr::llply(all_sub_data, function(x) x$counts))
  
  #rbind all cell infos
  all_cell_info <- plyr::ldply(names(all_sub_data), function(name) {
    all_sub_data[[name]]$cell_info %>% mutate(condition = name)
  }) %>% tibble::column_to_rownames(var = 'barcode')
  
  print(sprintf('Loaded %d of %d total cells', ncol(all_counts), tot_cells))
  return(list(counts = all_counts, gene_info = all_gene_info, cell_info = all_cell_info))
}


#' Convert gene labels to hgnc symbols (rather than ensemble ids)
#'
#' @param dat Seurat object
#'
#' @return
#' @export
#'
#' @examples
convert_to_hugo_symbols <- function(dat) {
  dat$gene_info %<>% mutate(Gene_Symbol = make.unique(Gene_Symbol))
  ensemble_to_hugo <- dat$gene_info$Gene_Symbol %>% set_names(dat$gene_info$Ensembl_ID)
  rownames(dat$counts) %<>% ensemble_to_hugo[.] %>% set_names(NULL)
  stopifnot(!any(is.na(rownames(dat$counts))))
  return(dat)
}


#' Calculate total reads (summed across cells) for each group of cells
#'
#' @param seuDat Seurat object
#' @param group_var grouping variable
#'
#' @return
#' @export
#'
#' @examples
get_summed_counts <- function(seuDat, group_var = 'singlet_ID') {
  md <- seuDat@meta.data %>% 
    rownames_to_column(var = 'cell_ID')
  stopifnot(group_var %in% colnames(md))
  md$gvar <- md[[group_var]]
  all_groups <- unique(md$gvar) 
  summed_counts <- plyr::laply(all_groups, function(cur_group) {
    cur_cells <- md %>% 
      dplyr::filter(gvar == cur_group) %>% 
      .[['cell_ID']]
    Matrix::rowSums(GetAssayData(seuDat, slot = 'counts')[, cur_cells, drop = FALSE])
  }) %>% set_rownames(all_groups)
  return(summed_counts)
}


#' Apply log counts-per million transformation, with optional library size adjustment
#'
#' @param counts 
#' @param prior_cnt 
#' @param nf: Optional vector of normalization factors 
#'
#' @return
#' @export
#'
#' @examples
lcpm_trans <- function(counts, prior_cnt, nf = 1) {
  cpm <- counts + prior_cnt
  lcpm <- scale(cpm, center = F, scale = colSums(cpm) * nf / 1e6) %>%
    log2()
  return(lcpm)
}


#' Take Seurat object with cell cycle phase classified, and relabel G1 to G0/G1 and G2M to G2/M
#'
#' @param seuDat 
#'
#' @return seuRat object with Phase metadata relabeled 
#' @export
#'
#' @examples
relabel_cell_cycle_phase <- function(seuDat) {
  stopifnot('Phase' %in% colnames(seuDat@meta.data))
  seuDat@meta.data$Phase <- plyr::revalue(seuDat@meta.data$Phase, replace = c(G1 = 'G0/G1', G2M = 'G2/M'))
  return(seuDat)
}


###------------------------- STATS -----------------------------------

#' Run differential expression analysis to compute viability-related and -independent components of response
#'
#' @param all_counts matrix of sum-collapsed gene counts across samples
#' @param sample_info #sample info table
#' @param sensitivity_df #dataframe of perturbation sensitivity values (with CCLE_ID cell line names)
#' @param prior_cnt 
#' @param use_voom #default FALSE for limma-trend
#' @param norm_meth 
#' @param min_counts_per_gene 
#' @param min_det_samples 
#'
#' @return
#' @export
#'
#' @examples
  fit_viability_models <- function(all_counts, sample_info, sensitivity_df = NULL, prior_cnt = 1, used_genes = NULL,
                                 use_voom = FALSE, use_edgeR = FALSE, norm_meth = 'none', min_counts_per_gene = 5, min_det_samples = 0.05) {
  
  if (is.null(used_genes)) {
    used_genes <- which(rowSums(all_counts > min_counts_per_gene) > ncol(all_counts)*min_det_samples)
  }
  print(sprintf('Using %d/%d genes', length(used_genes), nrow(all_counts)))
  dge <- edgeR::DGEList(all_counts[used_genes,])
  dge <- edgeR::calcNormFactors(dge, method = 'TMMwzp')

  design <- model.matrix(~0 + expt + CCLE_ID, data = sample_info)

  treat_group <- as.numeric(grepl('expttreat', colnames(design))) %>% magrittr::divide_by(., sum(.))
  control_group <- as.numeric(grepl('exptcontrol', colnames(design))) %>% magrittr::divide_by(., sum(.))
  avg_contrast <- treat_group-control_group
  if (use_edgeR) {
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design, prior.count = prior.count)
    qlf <- glmQLFTest(fit, contrast = avg_contrast)
    res_avg <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
  } else {
    if (use_voom) {
      dd <- limma::voom(dge, design, plot=FALSE, normalize.method = norm_meth)
      utrend = FALSE
    } else {
      # dd <- lcpm_trans(dge$counts, prior_cnt = prior_cnt, nf = dge$samples$norm.factors)
      dd <- edgeR::cpm(dge, prior.count = prior_cnt, log = TRUE)
      utrend = TRUE
    }
    fit <- limma::lmFit(dd, design)
    fit <- limma::contrasts.fit(fit, avg_contrast)
    fit2 <- limma::eBayes(fit, trend = utrend)
    res_avg <- limma::topTable(fit2, number = Inf) %>%
      rownames_to_column(var = 'Gene')
  }
  if (!is.null(sensitivity_df)) {
    sample_info %<>% left_join(sensitivity_df, by = "CCLE_ID")
    usamps <- which(!is.na(sample_info$sens))
    sample_info %<>% mutate(is_treat = ifelse(condition == 'treat',1, 0))
    design <- model.matrix(~0 + expt + CCLE_ID + I(is_treat*sens), data = sample_info[usamps,])

    colnames(design)[ncol(design)] <- 'treat_sens'
    colnames(design) %<>% make.names()
    treat_group <- as.numeric(grepl('expttreat', colnames(design))) %>% magrittr::divide_by(., sum(.))
    control_group <- as.numeric(grepl('exptcontrol', colnames(design))) %>% magrittr::divide_by(., sum(.))
    int_contrast <- treat_group-control_group
    slope_contrast <- as.numeric(colnames(design) == 'treat_sens')
    contr_mat <- rbind(int_contrast, slope_contrast) %>% set_rownames(c('intercept', 'slope')) %>% t()

    if (use_edgeR) {
      dge <- estimateDisp(dge[, usamps], design = design)
      fit <- glmQLFit(dge, design = design, prior.count = prior.count)
      qlf_int <- glmQLFTest(fit, coef = 'intercept')
      res_int <- topTags(qlf_int, n = Inf)$table %>% rownames_to_column(var = 'Gene')
      qlf_slope <- glmQLFTest(fit, coef = 'slope')
      res_slope <- topTags(qlf_slope, n = Inf)$table %>% rownames_to_column(var = 'Gene')
    } else {
      dd <- dd[, usamps] #restrict to samples with sensitivity data
    
      #run limma analysis
      fit <- limma::lmFit(dd, design)
      fit <- limma::contrasts.fit(fit, contrasts = contr_mat)
      fit2 <- limma::eBayes(fit, trend = utrend)
      
      res_int <- limma::topTable(fit2, number = Inf, coef = 'intercept') %>%
        rownames_to_column(var = 'Gene')
      res_slope <- limma::topTable(fit2, number = Inf, coef = 'slope') %>%
        rownames_to_column(var = 'Gene')
    }
  } else {
    res_int <- NULL
    res_slope <- NULL
  }
  return(list(res_avg = res_avg, res_int = res_int, res_slope = res_slope))
}


#' Title: Run test of associations between treatment DE response and a specified cell line feature
#'
#' @param cur_expt 
#' @param results_dir 
#' @param CL_feature named vector of cell line feature values
#' @param prior_cnt 
#'
#' @return
#' @export
#'
#' @examples
run_DE_association_test <- function(cur_expt, results_dir, CL_feature, prior_cnt = 1, covar = NULL) {
  dat_sum <- load_collapsed_profiles(cur_expt, results_dir, type = 'sum')
  dat_sum$sample_info %<>% filter(CCLE_ID %in% names(CL_feature[!is.na(CL_feature)]))
  dat_sum$sample_info %<>% left_join(data_frame(CCLE_ID = names(CL_feature), CL_feature = CL_feature), by = 'CCLE_ID')
  dat_sum$sample_info %<>% mutate(CL_feature = ifelse(condition == 'treat', CL_feature, 0)) #make interaction term with treatment
  if (!is.null(covar)) {
    dat_sum$sample_info %<>% left_join(data_frame(CCLE_ID = names(covar), covar = covar), by = 'CCLE_ID')
    dat_sum$sample_info %<>% mutate(covar = ifelse(condition == 'treat', covar, 0)) #make interaction term with treatment
  }
  dge <- edgeR::DGEList(dat_sum$profile_mat)
  dge <- edgeR::calcNormFactors(dge, method = 'TMMwzp')
  # dd <- lcpm_trans(dge$counts, prior_cnt = prior_cnt, nf = dge$samples$norm.factors)
  dd <- edgeR::cpm(dge, prior.count = prior_cnt, log = TRUE)
  
  dat_sum$sample_info %<>% filter(!is.na(CL_feature))
  frm <- '~condition + CCLE_ID + CL_feature'
  if (!is.null(covar)) {
    frm <- '~condition + CCLE_ID + CL_feature + covar'
  }
  design <- model.matrix(as.formula(frm), data = dat_sum$sample_info)
  fit <- limma::lmFit(dd[, dat_sum$sample_info$sample_name], design)
  fit2 <- limma::eBayes(fit, trend = TRUE)
  res <- limma::topTable(fit2, number = Inf, coef = 'CL_feature') %>%
    rownames_to_column(var = 'Gene')
  return(res)
}


#' Title: Compare the treat-and-control responses of one cell line
#'
#' @param seuDat 
#' @param min_frac_cells_det 
#' @param prior.count 
#' @param use_norm 
#' @param CL_name 
#'
#' @return
#' @export
#'
#' @examples
estimate_CL_responses <- function(seuDat, CL_name, min_frac_cells_det = 0.05, prior.count = 0.125, 
                                  use_norm = TRUE, method = 'edgeRQLF', cond2 = 'treat', cond1 = 'control') {
  library(edgeR)
  library(Seurat)
  cn <- FetchData(seuDat, 'singlet_ID')
  seuSub <- seuDat[, which(cn == CL_name)]
  used_genes <- which(Matrix::rowSums(GetAssayData(seuSub, 'counts')[, rownames(seuSub@meta.data)] > 0) > nrow(seuSub@meta.data)*min_frac_cells_det)
  dge <- DGEList(GetAssayData(seuSub, 'counts')[used_genes,rownames(seuSub@meta.data)])
  if (use_norm) {
    dge <- calcNormFactors(dge, method = 'TMMwzp')
  }
  #build design matrix for treat-vs-control comparison (avg across cell lines)
  df <- seuSub@meta.data 
  design <- model.matrix(~0 + condition + cell_det_rate, data = df)
  treat_group <- as.numeric(grepl(cond2, colnames(design))) %>% magrittr::divide_by(., sum(.))
  control_group <- as.numeric(grepl(cond1, colnames(design))) %>% magrittr::divide_by(., sum(.))
  #get differential treatment effect
  contr_vec <- treat_group-control_group

  if (method == 'edgeRQLF') {
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design, prior.count = prior.count)
    qlf <- glmQLFTest(fit, contrast = contr_vec)
    tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
  } else if (method == 'limma-voom') {
    dd <- limma::voom(dge, design, plot=FALSE)
    fit <- limma::lmFit(dd, design)
    fit <- limma::contrasts.fit(fit, contr_vec)
    fit2 <- limma::eBayes(fit)
    tt_res <- limma::topTable(fit2, number = Inf) %>%
      rownames_to_column(var = 'Gene')
  } else {
    stop('Unsupported method')
  }
  return(tt_res)
}


#' Title: Compare the treat-vs-control responses between two cell lines using edgeR
#'
#' @param seuDat 
#' @param CL1 
#' @param CL2 
#' @param min_frac_cells_det 
#' @param prior.count 
#' @param use_norm 
#'
#' @return
#' @export
#'
#' @examples
compare_CL_responses <- function(seuDat, CL1, CL2, min_frac_cells_det = 0.05, prior.count = 0.125, use_norm = TRUE) {
  library(edgeR)
  library(Seurat)
  cid <- FetchData(seuDat, 'singlet_ID')
  seuSub <- seuDat[, cid$singlet_ID %in% c(CL1, CL2)]
  used_genes <- which(Matrix::rowSums(GetAssayData(seuSub, 'counts')[, rownames(seuSub@meta.data)] > 0) > nrow(seuSub@meta.data)*min_frac_cells_det)
  dge <- DGEList(GetAssayData(seuSub, 'counts')[used_genes,rownames(seuSub@meta.data)])
  if (use_norm) {
    dge <- calcNormFactors(dge, method = 'TMMwzp')
  }
  #build design matrix for treat-vs-control comparison (avg across cell lines)
  df <- seuSub@meta.data %>% 
    mutate(CL_cond = paste0(condition, '_', singlet_ID))
  design <- model.matrix(~0 + CL_cond + cell_det_rate, data = df)
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design, prior.count = prior.count)
  
  treat_group1 <- as.numeric(grepl('treat', colnames(design)) & grepl(CL1, colnames(design))) %>% magrittr::divide_by(., sum(.))
  control_group1 <- as.numeric(grepl('control', colnames(design)) & grepl(CL1, colnames(design))) %>% magrittr::divide_by(., sum(.))
  treat_group2 <- as.numeric(grepl('treat', colnames(design)) & grepl(CL2, colnames(design))) %>% magrittr::divide_by(., sum(.))
  control_group2 <- as.numeric(grepl('control', colnames(design)) & grepl(CL2, colnames(design))) %>% magrittr::divide_by(., sum(.))
  
  #get differential treatment effect
  contr_vec <- treat_group1-control_group1 - (treat_group2 -control_group2)
  qlf <- glmQLFTest(fit, contrast = contr_vec)
  tt_diff <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
  
  #CL1 treatment effect
  contr_vec <- treat_group1-control_group1 
  qlf <- glmQLFTest(fit, contrast = contr_vec)
  tt_CL1 <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
  
  #CL2 treatment effect
  contr_vec <- treat_group2-control_group2
  qlf <- glmQLFTest(fit, contrast = contr_vec)
  tt_CL2 <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
  
  #combine results
  tt_comb <- full_join(tt_CL1, tt_CL2, by = 'Gene', suffix = c('_CL1', '_CL2'))
  tt_comb %<>% left_join(tt_diff %>% 
                           dplyr::select(Gene, PValue_diff = PValue, logFC_diff = logFC, FDR_diff = FDR))
  return(tt_comb)
}




###------------------------- PLOTTING -----------------------------------

#' Title
#'
#' @param gene_stat: named vector of gene-level stats
#' @param gsc: Gene-set collection object 
#' @param top_n Number of top up- and down- genes to use for hyper-geom enrichment test
#' @param gsSizeLim Two element vector with min and max gsc sizes
#' @param n_lab_per Number of top gs to label
#'
#' @return
#' @export
#'
#' @examples
make_hyper_gsa_plot <- function(gene_stat, gsc, top_n = 50, lfc_thresh = NULL, gsSizeLim = c(1, Inf), n_lab_per = 15, dir = 'both', universe = NULL, lab_size = 2.5, return_stats = FALSE, max_chars = 35) {
  mod_hgnc_symbols <- grepl('\\.', names(gene_stat))
  print(sprintf('Dropping %s/%s genes with modified hgnc symbols', sum(mod_hgnc_symbols), length(gene_stat)))
  gene_stat <- gene_stat[!mod_hgnc_symbols]
  if (is.null(lfc_thresh)) {
    top_down_hits <- gene_stat %>% 
      sort() %>% 
      head(top_n) %>% 
      names()
    top_up_hits <- gene_stat %>% 
      sort(decreasing = TRUE) %>% 
      head(top_n) %>% 
      names()
  } else {
    top_down_hits <- names(gene_stat[gene_stat <= -lfc_thresh])
    top_up_hits <- names(gene_stat[gene_stat >= lfc_thresh])
  }
  if (is.null(universe)) {
    universe <- names(gene_stat)
  }
  if (dir == 'up' | dir == 'both') {
    cur_up_GSEA <- cdsr::run_GSAhyper(top_up_hits, universe, gsc, gsSizeLim = gsSizeLim)
  }
  if (dir == 'down' | dir == 'both') {
    cur_down_GSEA <- cdsr::run_GSAhyper(top_down_hits, universe, gsc, gsSizeLim = gsSizeLim) 
  }
  if (return_stats) {
    stats_return <- list(up_GSEA = cur_up_GSEA, down_GSEA = cur_down_GSEA)
    return(stats_return)
  }
    if (dir == 'both' | dir == 'up') {
      cur_up_GSEA %<>%
        arrange(`p-value`) %>% 
        head(n_lab_per) %>% 
        mutate(gs_id = paste0(gene_set, '.up'),
               direction = 'up')
      cur_up_GSEA %<>% 
        map_df(rev) %>% 
        mutate(gs_id = factor(gs_id, levels = .[['gs_id']]))
    }
    if (dir == 'both' | dir == 'down') {
      cur_down_GSEA %<>%
        arrange(`p-value`) %>% 
        head(n_lab_per) %>% 
        mutate(gs_id = paste0(gene_set, '.down'),
               direction = 'down')
      cur_down_GSEA %<>% 
        map_df(rev) %>% 
        mutate(gs_id = factor(gs_id, levels = .[['gs_id']]))
    }
    if (dir == 'both') {
      
      pmax <- max(-log10(rbind(cur_up_GSEA, cur_down_GSEA)[['p-value']]))
      
      #add in points at 0 for stems
      cur_up_GSEA <- rbind(cur_up_GSEA %>% mutate(is_buff = FALSE), cur_up_GSEA %>% mutate(`p-value` = 1, is_buff = TRUE))
      cur_up_GSEA %<>% mutate(gene_set = str_replace_all(gene_set, '^REACTOME_|^PID_|^KEGG_|^HALLMARK_|^GO_', ''))
      cur_up_GSEA %<>% mutate(gene_set = str_trunc(gene_set, max_chars))
      g_up <- ggplot(cur_up_GSEA, 
                     aes(gs_id, -log10(`p-value`))) + 
        geom_point(data = cur_up_GSEA %>% filter(!is_buff), aes(color = direction), size = 2, alpha = 0.75) + 
        geom_line(aes(group = gs_id, color = direction)) +
        geom_text(data = cur_up_GSEA %>% filter(is_buff),
                  aes(label = gene_set), 
                  angle = 0, vjust = -0.5, hjust = 0, nudge_y = -0., size = lab_size) +
        # geom_vline(xintercept = n_gene_sets + 0.5, linetype = 'dashed') +
        cdsr::theme_Publication() + 
        ylab('-log10(p-value)') +
        ylim(0, pmax) +
        # xlab('gene set') +
        scale_color_manual(values = c(down = 'darkblue', up = 'darkred')) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank()) +
        coord_flip(expand = TRUE, clip = 'off') +
        guides(color = FALSE) +
        ggtitle('Up-regulated')
      
      col <- 'blue'
      cur_down_GSEA <- rbind(cur_down_GSEA %>% mutate(is_buff = FALSE), cur_down_GSEA %>% mutate(`p-value` = 1, is_buff = TRUE))
      cur_down_GSEA %<>% mutate(gene_set = str_replace_all(gene_set, '^REACTOME_|^PID_|^KEGG_|^HALLMARK_|^GO_', ''))
      cur_down_GSEA %<>% mutate(gene_set = str_replace_all(gene_set, '^REACTOME_|^PID_|^KEGG_|^HALLMARK_|^GO_', ''))
      cur_down_GSEA %<>% mutate(gene_set = str_trunc(gene_set, max_chars))
      g_down <- ggplot(cur_down_GSEA, 
                     aes(gs_id, -log10(`p-value`))) + 
        geom_point(data = cur_down_GSEA %>% filter(!is_buff), aes(color = direction), size = 2, alpha = 0.75) + 
        geom_line(aes(group = gs_id, color = direction)) +
        geom_text(data = cur_down_GSEA %>% filter(is_buff),
                  aes(label = gene_set), 
                  angle = 0, vjust = -0.5, hjust = 1, nudge_y = -0., size = lab_size) +
        # geom_vline(xintercept = n_gene_sets + 0.5, linetype = 'dashed') +
        cdsr::theme_Publication() + 
        ylab('-log10(p-value)') +
        xlab('gene set') +
        scale_color_manual(values = c(down = 'darkblue', up = 'darkred')) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        scale_y_reverse(limits = c(pmax, 0)) +
        coord_flip(expand = TRUE, clip = 'off') +
        guides(color = FALSE) +
        ggtitle('Down-regulated') 
        
      g <- cowplot::plot_grid(g_down, g_up)
      
    }
  return(g)
}


#make stem plot from precomputed hyper gsea scores
make_stem_plot_precom <- function(cur_up_GSEA, cur_down_GSEA, n_lab_per = 10, lab_size = 2.5) {
  cur_up_GSEA %<>%
    arrange(`p-value`) %>% 
    head(n_lab_per) %>% 
    mutate(gs_id = paste0(gene_set, '.up'),
           direction = 'up')
  cur_up_GSEA %<>% 
    map_df(rev) %>% 
    mutate(gs_id = factor(gs_id, levels = .[['gs_id']]))
  cur_down_GSEA %<>%
    arrange(`p-value`) %>% 
    head(n_lab_per) %>% 
    mutate(gs_id = paste0(gene_set, '.down'),
           direction = 'down')
  cur_down_GSEA %<>% 
    map_df(rev) %>% 
    mutate(gs_id = factor(gs_id, levels = .[['gs_id']]))
  pmax <- max(-log10(rbind(cur_up_GSEA, cur_down_GSEA)[['p-value']]))
  
  #add in points at 0 for stems
  cur_up_GSEA <- rbind(cur_up_GSEA %>% mutate(is_buff = FALSE), cur_up_GSEA %>% mutate(`p-value` = 1, is_buff = TRUE))
  cur_up_GSEA %<>% mutate(gene_set = str_replace_all(gene_set, '^REACTOME_|^PID_|^KEGG_|^HALLMARK_|^GO_', ''))
  g_up <- ggplot(cur_up_GSEA, 
                 aes(gs_id, -log10(`p-value`))) + 
    geom_point(data = cur_up_GSEA %>% filter(!is_buff), aes(color = direction), size = 2, alpha = 0.75) + 
    geom_line(aes(group = gs_id, color = direction)) +
    geom_text(data = cur_up_GSEA %>% filter(is_buff),
              aes(label = gene_set), 
              angle = 0, vjust = -0.5, hjust = 0, nudge_y = -0., size = lab_size) +
    # geom_vline(xintercept = n_gene_sets + 0.5, linetype = 'dashed') +
    cdsr::theme_Publication() + 
    ylab('-log10(p-value)') +
    ylim(0, pmax) +
    # xlab('gene set') +
    scale_color_manual(values = c(down = 'darkblue', up = 'darkred')) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    coord_flip(expand = TRUE, clip = 'off') +
    guides(color = FALSE) +
    ggtitle('Up-regulated')
  
  col <- 'blue'
  cur_down_GSEA <- rbind(cur_down_GSEA %>% mutate(is_buff = FALSE), cur_down_GSEA %>% mutate(`p-value` = 1, is_buff = TRUE))
  cur_down_GSEA %<>% mutate(gene_set = str_replace_all(gene_set, '^REACTOME_|^PID_|^KEGG_|^HALLMARK_|^GO_', ''))
  g_down <- ggplot(cur_down_GSEA, 
                   aes(gs_id, -log10(`p-value`))) + 
    geom_point(data = cur_down_GSEA %>% filter(!is_buff), aes(color = direction), size = 2, alpha = 0.75) + 
    geom_line(aes(group = gs_id, color = direction)) +
    geom_text(data = cur_down_GSEA %>% filter(is_buff),
              aes(label = gene_set), 
              angle = 0, vjust = -0.5, hjust = 1, nudge_y = -0., size = lab_size) +
    # geom_vline(xintercept = n_gene_sets + 0.5, linetype = 'dashed') +
    cdsr::theme_Publication() + 
    ylab('-log10(p-value)') +
    xlab('gene set') +
    ylim(0, pmax) +
    scale_color_manual(values = c(down = 'darkblue', up = 'darkred')) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_y_reverse() +
    coord_flip(expand = TRUE, clip = 'off') +
    guides(color = FALSE) +
    ggtitle('Down-regulated')
  
  g <- cowplot::plot_grid(g_down, g_up)
  return(g)
}


#' Title Wrapper around pheatmap function
#'
#' @param LFC_mat 
#' @param gene_list 
#' @param cluster_genes 
#' @param CL_ann_df 
#' @param color_lims 
#' @param CL_list 
#' @param cluster_rows 
#' @param cluster_cols 
#' @param transpose 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
make_LFC_heatmap <- function(LFC_mat,
         gene_list,
         cluster_genes = NULL,
         CL_ann_df = NULL,
         color_lims = NULL,
         CL_list = NULL,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         dist_meth = 'euclidean',
         transpose = FALSE, ...) {
  library(pheatmap)
  library(grDevices)
  library(RColorBrewer)
  colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                       "RdYlBu")))(100)
  if (is.null(cluster_genes)) {
    cluster_genes <- gene_list
  }
  if (is.null(CL_list) & cluster_cols) {
    cclust <- hclust(dist(t(LFC_mat[cluster_genes,])))
  } else {
    cclust <- cluster_cols
  }
  if (cluster_rows) {
    cluster_rows = hclust(dist(LFC_mat[gene_list,]))
  }
  if (is.null(color_lims)) {
    breaksList <- NA
  } else {
    breaksList = seq(color_lims[1], color_lims[2], length.out = 100)
    LFC_mat[LFC_mat > color_lims[2]] <- color_lims[2]
    LFC_mat[LFC_mat < color_lims[1]] <- color_lims[1]
  }
  if (is.null(CL_list)) {
    if (transpose) {
      pheatmap::pheatmap(t(LFC_mat[gene_list,]),
                         annotation_row = CL_ann_df,
                         cluster_cols = cluster_rows,
                         cluster_rows = cclust,
                         color = colors,
                         breaks = breaksList,
                         ...)
    } else {
      pheatmap::pheatmap(LFC_mat[gene_list,],
                         annotation_col = CL_ann_df,
                         cluster_rows = cluster_rows,
                         cluster_cols = cclust,
                         color = colors,
                         breaks = breaksList,
                         ...)
    }
  } else {
    if (transpose) {
      pheatmap::pheatmap(t(LFC_mat[gene_list,CL_list]),
                         annotation_row = CL_ann_df,
                         cluster_cols = cluster_rows,
                         cluster_rows = F,
                         color = colors,
                         breaks = breaksList,
                         ...)
    } else {
      pheatmap::pheatmap(LFC_mat[gene_list,CL_list],
                         annotation_col = CL_ann_df,
                         cluster_rows = cluster_rows, cluster_cols = F,
                         color = colors,
                         breaks = breaksList,
                         ...)
    }
  }
}




