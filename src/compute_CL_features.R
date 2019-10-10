library(magrittr)
library(tidyverse)
library(taigr)
library(here)

source(here::here('src', 'MixSeq_helpers.R'))
source(here::here('src', 'expt_meta_data.R'))
results_dir <- here::here('data')

#OMICS data
omics_dsets <- list(
  GE = list(data.name="depmap-rnaseq-expression-data-ccd0",
            data.version = 14,
            data.file = 'CCLE_depMap_19Q3_TPM_ProteinCoding'),
  MUT_HOT = list(data.name='depmap-mutation-calls-9a1a', data.version=12, data.file='hotspot_mutation'),
  MUT_DAM = list(data.name='depmap-mutation-calls-9a1a', data.version=12, data.file='damaging_mutation')
)
omics_dat <- taigr::load.all.from.taiga(omics_dsets) %>%    
  cdsr::extract_hugo_symbol_colnames() 

#PRISM Repurposing data
all_PRISM_AUC <- load.from.taiga(data.name='secondary-screen-0854', data.version=9, data.file='secondary_dose_response_curve_parameters') %>% 
  dplyr::rename(DEPMAP_ID = depmap_id) %>% 
  filter(!is.na(DEPMAP_ID))

#GDSC data
all_GDSC_AUC <- load.from.taiga(data.name = 'sanger-drug-sensitivities-d265',
                            data.version = 5,
                            data.file = 'v17.3_fitted_dose_response',
                            transpose = F)
arxspan <- load.from.taiga(data.name='arxspan-cell-line-export-f808', data.version=322, data.file='ACH') %>% 
  dplyr::select(COSMIC_ID = COSMIC.ID, DEPMAP_ID = arxspan_id)
all_GDSC_AUC %<>% left_join(arxspan, by = 'COSMIC_ID')

#CRISPR data
CRISPR <- load.from.taiga(data.name='avana-public-19q3-0900', data.version=5, data.file='gene_effect') %>% 
  cdsr::extract_hugo_symbol_colnames()


#' Title
#'
#' @param drug_name 
#'
#' @return
#' @export
#'
#' @examples
get_drug_sensitivity_data <- function(drug_name) {
  drug_name <- tolower(drug_name)
  if (is.null(drug_name)) {
    AUC_data <- data.frame(DEPMAP_ID = character(), stringsAsFactors = F)
    return(AUC_data)
  }
  
  #PRISM data
  if (drug_name == 'brd3379') {
    PRISM_AUC <- read_csv('~/CPDS/Apollo/data/assaydata_BRD-K01353379-034-08-4.csv') %>%
      dplyr::rename(CCLE_ID = Cell_line, PRISM_AUC = `10uM`) %>%
      dplyr::mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(CCLE_ID, ignore.problems = T), 
             PRISM_AUC = log2(PRISM_AUC / 100),
             PRISM_AUC = (pmax(PRISM_AUC, -3) + 3)/3) %>% #cap at -3 and map this to 0-1 range
      dplyr::select(DEPMAP_ID, PRISM_AUC) %>%
      dplyr::filter(!is.na(DEPMAP_ID)) %>% 
      as.data.frame()
  } else if (drug_name == 'azd5591') {
    PRISM_AUC <- read_csv(here::here('data', 'MCL1_Inhibitors_MTS011/BRD-K00005264-001-01-9/DRC_table.csv')) %>% 
      dplyr::select(CCLE_ID = ccle_name,
                    PRISM_AUC = auc) %>% 
      dplyr::mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(CCLE_ID, ignore.problems = T)) %>% 
      dplyr::filter(!is.na(DEPMAP_ID)) %>% 
      dplyr::select(DEPMAP_ID, PRISM_AUC) %>% 
      dplyr::group_by(DEPMAP_ID) %>% 
      dplyr::summarise(PRISM_AUC = mean(PRISM_AUC, na.rm=T)) %>% 
      dplyr::ungroup()
     
    # PRISM_AUC <- read_csv('~/CPDS/jmm-scratch/MCL1/results/all_comb_curve_fits_D3.csv') %>% 
    #   filter(brd_id == 'Servier-Nat') %>% 
    #   dplyr::select(CCLE_ID = cell_line_id,
    #                 PRISM_AUC = auc) %>% 
    #   dplyr::mutate(DEPMAP_ID = celllinemapr::ccle.to.arxspan(CCLE_ID, ignore.problems = F)) %>% 
    #   dplyr::select(DEPMAP_ID, PRISM_AUC)
  } else {
    cur_broad_id <- all_PRISM_AUC %>% 
      distinct(broad_id, name) %>% 
      filter(grepl(drug_name, name, ignore.case = TRUE)) %>% 
      .[['broad_id']]
    stopifnot(length(cur_broad_id) < 2)

    if (length(cur_broad_id) == 1) {
      PRISM_AUC <- all_PRISM_AUC %>% 
        dplyr::filter(broad_id == cur_broad_id) %>%
        dplyr::group_by(DEPMAP_ID) %>% 
        dplyr::summarise(PRISM_AUC = mean(auc, na.rm=T)) %>% 
        dplyr::ungroup()
    } else {
      print('Couldnt find matching compound in PRISM data')
      PRISM_AUC <- data.frame(DEPMAP_ID = character(),
                              PRISM_AUC = numeric())
    }
  } 
  
  if (drug_name == 'idasanutlin') {
    cur_gdsc_drug <- 'nutlin'
  } else {
    cur_gdsc_drug <- drug_name
  }
  gdsc_match_compounds <- grep(cur_gdsc_drug, unique(all_GDSC_AUC$DRUG_NAME), ignore.case = TRUE, value = T)
  if (drug_name == 'jq1') {gdsc_match_compounds <-  'JQ1'}
  if (length(gdsc_match_compounds) > 1) {
    print('warning, more than one matching GDSC compounds detected')
    stop()
    gdsc_match_compounds <- gdsc_match_compounds[[1]]
  }
  if (length(gdsc_match_compounds) > 0) {
    GDSC_AUC <- all_GDSC_AUC %>% 
      dplyr::filter(DRUG_NAME == gdsc_match_compounds) %>% 
      dplyr::group_by(DEPMAP_ID) %>% 
      dplyr::summarise(GDSC_AUC = mean(AUC, na.rm=T)) %>%  #in some cases there are multiple profiles for each cell line for a given drug
      dplyr::ungroup()
  } else {
    print('Couldnt find matching compound in GDSC data')
    GDSC_AUC <- data.frame(DEPMAP_ID = character(), GDSC_AUC = numeric())
  }
  
  AUC_data <- full_join(PRISM_AUC, GDSC_AUC, by = 'DEPMAP_ID')
  
  #normalize then average together the GDSC and PRISM AUCs
  AUC_mat <- as.matrix(AUC_data[, c('PRISM_AUC', 'GDSC_AUC')])
  if (max(rowSums(!is.na(AUC_mat))) > 1) {
    AUC_data$AUC_avg <- rowMeans(limma::normalizeQuantiles(AUC_mat), na.rm=T)
  } else {
    AUC_data$AUC_avg <- rowMeans(AUC_mat, na.rm=T)
  }
  return(AUC_data)
} 


#' Title
#'
#' @param gene_deps
#'
#' @return
#' @export
#'
#' @examples
get_gene_dep_data <- function(gene_deps) {
  if (is.null(gene_deps)) {
    return(data.frame(CCLE_ID = character()))
  }
  df <- CRISPR[, gene_deps, drop=F] %>%
    as.data.frame() %>%
    rownames_to_column(var = 'DEPMAP_ID')
  colnames(df)[colnames(df) != 'DEPMAP_ID'] %<>% paste0('CRISPR_', .)
  return(df)
}


#' Title
#'
#' @param omics_features: list of lists with dataset and gene name specified for each feature
#'
#' @return
#' @export
#'
#' @examples
get_omics_features <- function(omics_features) {
  if (!is.list(omics_features)) {
    return(data.frame(DEPMAP_ID = character()))
  }
  
  omics_df <- llply(names(omics_features), function(cur_feat) {
    df <- data.frame(DEPMAP_ID = rownames(omics_dat[[omics_features[[cur_feat]]$dset]]))
    df[[cur_feat]] <- omics_dat[[omics_features[[cur_feat]]$dset]][, omics_features[[cur_feat]]$gene]
    return(df)
  }) %>% join_all(by = 'DEPMAP_ID', type = 'full') 
  return(omics_df)
}



all_CL_features <- llply(sc_DE_meta, function(cur_expt) {
  print(cur_expt$expt_name)
  sc_df <- load_cell_info(cur_expt)
  unique_CLs <- unique(sc_df$singlet_ID)
  in_pool_df <- data_frame(CCLE_ID = unique_CLs,
                           in_pool = TRUE,
                      DEPMAP_ID = celllinemapr::ccle.to.arxspan(unique_CLs))
  if (!is.null(cur_expt$drug_name)) {
    CL_df <- get_drug_sensitivity_data(cur_expt$drug_name) %>%
        dplyr::mutate(sens = 1 - AUC_avg) %>% 
        dplyr::filter(!is.na(DEPMAP_ID))
  } else {
    CL_df <- get_gene_dep_data(cur_expt$gene_name) 
    CL_df$sens <- -CL_df[[paste0('CRISPR_', cur_expt$gene_name)]]
  }
  CL_df %<>% full_join(get_omics_features(cur_expt$annotate_omics), by = "DEPMAP_ID")
  CL_df %<>% full_join(in_pool_df, by = 'DEPMAP_ID') 
  CL_df %<>% dplyr::mutate(in_pool = ifelse(is.na(in_pool), FALSE, in_pool))
  return(CL_df)
}) %>% set_names(sapply(sc_DE_meta, function(x) x$expt_name))
write_rds(all_CL_features, file.path(results_dir, 'all_CL_features.rds'))
